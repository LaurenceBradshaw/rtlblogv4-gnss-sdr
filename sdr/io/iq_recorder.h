#pragma once
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include "iq_format.h" // Iq_sample_format + encode_iq (source-agnostic codec)
#include "types.h"

// Asynchronous IQ recorder: a PASSIVE tap that writes the pipeline's sample stream to a file without
// affecting the live pipeline. record() encodes a copy to the chosen format and returns immediately; a
// writer thread drains the queue to disk. The queue is byte-bounded - if the disk cannot keep up it DROPS
// blocks (and warns at the end) rather than block the caller, so recording degrades, never the receiver.
//
// It taps the stream that feeds the channels, i.e. POST-decimation: with --decimate it records the
// decimated stream (a one-shot decimation pre-process - replay the smaller/slower file directly), and with
// no decimation it records the raw source stream (e.g. a long live RTL-SDR capture). The samples are encoded
// to `format` (encode_iq, the inverse of the file reader's decode), headerless, so a recording replays as a
// normal capture with `--format <that format> --sample-rate <processing rate>`.
class Iq_recorder
{
public:
    // Opens the file + starts the writer thread; throws on open failure. format = the on-disk sample format.
    Iq_recorder( const std::string& path, Iq_sample_format format );
    ~Iq_recorder(); // drains the queue, stops the thread, closes the file

    Iq_recorder( const Iq_recorder& )            = delete;
    Iq_recorder& operator=( const Iq_recorder& ) = delete;

    void record( const Complex_buf& samples ); // encode + enqueue (non-blocking); drops + counts on overflow

private:
    void writer_loop();

    std::ofstream                     out_;
    Iq_sample_format                  format_;
    std::thread                       writer_;
    std::mutex                        mtx_;
    std::condition_variable           cv_;
    std::deque<std::vector<uint8_t>>  queue_; // encoded byte blocks awaiting the writer
    size_t                            queued_bytes_   = 0;
    bool                              stop_           = false;
    uint64_t                          dropped_blocks_ = 0;

    static constexpr size_t MAX_QUEUE_BYTES = 256u * 1024 * 1024; // bound so a slow disk never grows RAM unbounded
};

// Ensure an IQ-recording path ends in a known extension - .f32 or .iq (case-insensitive) - defaulting to .f32
// if it has neither. Empty stays empty. Shared by the CLI and the GUI so a recording always has a clear,
// replayable extension. (The actual sample format is chosen separately - see Iq_recorder / record_format.)
std::string with_iq_extension( const std::string& path );

#pragma once

#include <atomic>
#include <cstdint>
#include <thread>
#include <vector>
#include "stream_device.h"

// rtl-sdr.h is only needed in the .cpp; forward-declare the opaque handle here so including this
// header does not require librtlsdr's headers.
struct rtlsdr_dev;

// Live RTL-SDR front-end. Tune (centre freq / sample rate / gain / AGC) BEFORE start_streaming();
// the librtlsdr async read runs on an internal worker thread, delivering converted samples through
// the Stream_device callback, so start_streaming() returns immediately (same contract as the file
// source). Configure for GNSS L1 (1575.42 MHz) before streaming.
class Rtlsdr_device : public Stream_device
{
public:
    explicit Rtlsdr_device( int device_index = 0 );

    // The device handle is a unique resource - not copyable.
    Rtlsdr_device( const Rtlsdr_device& )            = delete;
    Rtlsdr_device& operator=( const Rtlsdr_device& ) = delete;

    ~Rtlsdr_device() override;

    void     start_streaming( Sample_callback callback ) override;
    void     stop_streaming() override;
    bool     is_streaming() const override;
    uint32_t sample_rate_hz() const override;
    uint64_t samples_consumed() const override;

    void set_centre_freq_hz( uint32_t freq_hz ) override;
    void set_sample_rate_hz( uint32_t rate_hz ) override;
    void set_gain_tenths_db( int gain_tenths_db ) override;
    void set_agc( bool enable ) override;
    void set_bias_tee( bool enable ) override; // power an active antenna's LNA (rtlsdr_set_bias_tee)

    // Read back what the hardware ACTUALLY settled on (the device rounds rate/freq/gain to discrete
    // values) - used to log/verify the requested config really reached the device.
    uint32_t centre_freq_hz() const;
    int      tuner_gain_tenths_db() const;

private:
    // Runs librtlsdr's blocking async read; cancelled by stop_streaming().
    void        streaming_thread();
    // librtlsdr C callback: raw interleaved uint8 I/Q -> Complex_buf -> user callback.
    static void raw_callback( unsigned char* buf, uint32_t len, void* ctx );

    rtlsdr_dev*           dev_;
    Sample_callback       callback_;
    std::thread           worker_;
    std::atomic<bool>     streaming_ { false };
    std::atomic<uint64_t> samples_consumed_ { 0 };
};

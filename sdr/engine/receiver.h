#pragma once
#include <atomic>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <set>
#include <string>
#include <vector>
#include "acquisition_aiding.h"
#include "channel.h"
#include "fir_decimator.h"
#include "iq_file_device.h" // Iq_sample_format
#include "observation.h"
#include "position.h"
#include "signal_selection.h" // Signal_id

class Stream_device;
class Sample_buffer;
class Signal;
class Thread_pool;
class Scheduler;

// How to bring the receiver up: which source, and the front-end settings. Mirrors the CLI flags,
// but is decoupled from cxxopts so the receiver can be constructed from a GUI, a test, etc.
struct Receiver_config
{
    bool             use_rtlsdr     = false;
    std::string      file_path;                              // file source
    uint32_t         sample_rate_hz = 2048000;               // SOURCE (native) rate
    Iq_sample_format format         = Iq_sample_format::INT16; // file source
    int              device_index   = 0;                       // RTL-SDR
    double           gain_db        = -1.0;                    // RTL-SDR; <0 => hardware AGC
    // FIR-decimate the source by this integer factor before processing (1 = none). The whole pipeline
    // then runs at sample_rate_hz / decimation - cheaper, and lets a high-rate capture (e.g. 25 MHz)
    // run near real time. The fractional rate that an indivisible factor leaves is harmless.
    uint32_t         decimation     = 1;

    // Which signals to search, each with its own PRN allowlist (empty set = all PRNs in that signal's
    // range). Empty list -> default_signal_selection() (GPS L1 C/A + Galileo E1-B, all PRNs).
    // (Not named `signals` - that is a Qt macro in the GUI build.)
    std::vector<Signal_selection> selected_signals;
};

// The source / run parameters the GUI's Source tab can change between runs (the Receiver_config source
// fields, minus the signal selection which has its own setter). Applied via set_source_params(); takes
// effect on the next run() like the signal selection.
struct Source_params
{
    bool             use_rtlsdr     = false;
    std::string      file_path;
    Iq_sample_format format         = Iq_sample_format::INT16;
    uint32_t         sample_rate_hz = 2048000; // SOURCE (native) rate
    int              device_index   = 0;       // RTL-SDR
    double           gain_db        = -1.0;    // RTL-SDR; <0 => hardware AGC
    uint32_t         decimation     = 1;       // FIR-decimate the source by this factor (1 = none)
};

// Timing/progress of the processing loop, published for the GUI status bar. exec = wall-clock since
// streaming started; stream = seconds of IQ pushed into the ring buffer (can outrun exec - disk beats
// real time); channel = the slowest channel's processed position (what is actually correlated). When
// channel lags exec the receiver is behind real time.
struct Receiver_status
{
    double exec_s    = 0.0;
    double stream_s  = 0.0;
    double channel_s = 0.0;
};

// A (constellation, PRN) the receiver is configured to search. Known from construction (before Start),
// so the GUI can list every satellite - idle until it actually acquires.
struct Configured_satellite
{
    Constellation constellation = Constellation::Unknown;
    int           prn           = 0;
    Code          code          = Code::CA; // signal component - so the GUI idle row matches its code
};

// Owns the whole receiver pipeline (source -> channels -> observables -> PVT) and drives its
// processing loop. run() blocks until the stream ends or stop() is called, so it is meant to be
// run on its own thread when a GUI is attached, or called directly for headless operation.
//
// The read accessors (latest_position / channel_snapshots) are thread-safe: a GUI polling on
// another thread reads a coherent copy of the published state, never the worker-mutated internals.
class Receiver
{
public:
    explicit Receiver( Receiver_config config );
    ~Receiver();

    Receiver( const Receiver& )            = delete;
    Receiver& operator=( const Receiver& ) = delete;

    // Build the pipeline and run the processing loop. Blocks. Throws on setup failure
    // (e.g. missing file / no device) - the caller wraps it.
    void run();
    // Ask run() to exit at the next tick. Thread-safe and idempotent.
    void stop();

    // ---- thread-safe read-only views (for a GUI on another thread) ----
    // The most recent PVT solution (nullopt until the filter is seeded).
    std::optional<Position_solution> latest_position() const;
    // A coherent snapshot of every channel's observable state, one entry per channel.
    std::vector<Channel_snapshot> channel_snapshots() const;
    // The receiver-wide accumulated almanac (sv_key(con,prn) -> coarse orbit), aggregated across all channels and
    // persisting across the run. For acquisition aiding (Phase 3) and the GUI. Thread-safe copy.
    std::map<int, Almanac> almanac() const
    {
        std::lock_guard<std::mutex> lock( state_mutex_ );
        return almanac_;
    }
    // Loop timing/progress for the status bar.
    Receiver_status status() const;
    // Every (constellation, PRN) the receiver searches - available before Start (does not need the
    // pipeline up), so the GUI can list all satellites even while idle.
    const std::vector<Configured_satellite>& configured_satellites() const
    {
        return configured_sats_;
    }

    // The signals + PRNs the next run will search. The GUI selection tab reads this to initialise and
    // writes it back via set_signal_selection() before a (re)Start; takes effect on the next setup().
    const std::vector<Signal_selection>& signal_selection() const
    {
        return signal_selection_;
    }
    // Replace the signal selection (set exactly; empty = search nothing) and re-enumerate
    // configured_satellites(). Only meaningful while stopped (the next run() rebuilds channels from it,
    // and the GUI reads configured_satellites() for its live pre-Start list); ignored while running.
    void set_signal_selection( std::vector<Signal_selection> selection );

    // Current config - the GUI Source tab reads this to initialise its fields (so they reflect the CLI
    // args / defaults). Available before Start (does not need the pipeline up).
    const Receiver_config& config() const
    {
        return config_;
    }
    // Replace the source/run params (file vs RTL-SDR, file path, format, native rate, device, gain,
    // decimation) for the NEXT run. Only meaningful while stopped; ignored mid-run, takes effect on the
    // next setup(). The signal selection is set separately via set_signal_selection().
    void set_source_params( Source_params params );

    // ---- Per-SV graph history (subscription model) ----
    // The GUI subscribes the SVs whose graph windows are open; the run loop publishes only those each
    // tick. So the GUI never touches a live channel - it reads published copies, like the snapshots. One
    // subscription per SV publishes the whole history (I/Q + Doppler), so every graph type shares it.
    // const + mutable state (guarded by state_mutex_) so a const Receiver& suffices for the GUI.
    // Keyed by (constellation, PRN, CODE) so one SV tracked on two components (e.g. GPS L1CA + L1C) has
    // an independent history per code - the GUI shows a separate row + graphs for each.
    void subscribe_history( Constellation constellation, int prn, Code code, bool on ) const;
    std::optional<Tracking_history::Snapshot> published_history( Constellation constellation, int prn, Code code ) const;

private:
    void setup();               // construct device/buffer/signals/channels/pool/scheduler
    void teardown();            // tear them down in dependency order
    void enumerate_configured_sats(); // fill configured_sats_ from signal_selection_ (+ PRN filters)
    // From a coarse fix + the decoded almanac, predict each configured SV's visibility + LOS Doppler and
    // push them to aiding_ (so acquisition can skip below-horizon SVs and tightly window the rest).
    void update_acquisition_predictions( const Position_solution& fix );
    void publish_histories();   // run-loop helper: copy subscribed channels' history to the published map
    void freeze_histories();    // at EOF: snapshot ALL data-having channels so graphs keep their last state
    void clear_published_state(); // reset the GUI-visible published state (on Stop / fresh Start)

    static int history_key( Constellation constellation, int prn, Code code )
    {
        // sv_key * 10 + code: Code has < 10 values, so this stays collision-free and keeps the two codes
        // of one SV (e.g. GPS L1CA / L1C) on distinct keys.
        return sv_key( constellation, prn ) * 10 + static_cast<int>( code );
    }

    Receiver_config                   config_;
    std::vector<Signal_selection>     signal_selection_; // resolved from config_ (or the default set)
    std::vector<Configured_satellite> configured_sats_;  // fixed at construction
    std::atomic<bool>                 running_ { false };

    std::unique_ptr<Sample_buffer>        sample_buffer_;
    std::unique_ptr<Stream_device>        device_;
    std::unique_ptr<Fir_decimator>        decimator_; // optional layer: device -> decimator -> buffer
    std::vector<std::unique_ptr<Signal>>  signals_;
    Acquisition_aiding                    aiding_;
    std::vector<std::unique_ptr<Channel>> channels_;
    std::vector<Channel*>                 channel_ptrs_;
    std::unique_ptr<Thread_pool>          pool_;
    std::unique_ptr<Scheduler>            scheduler_;
    Observation_engine                    obs_engine_;
    Position_solver                       pos_solver_;

    // Published state for GUI readers: the run thread refreshes these under state_mutex_ every
    // tick; the accessors return copies. Readers thus never touch the live channels (which the run
    // thread destroys at teardown), so there is no data race / use-after-free across threads.
    mutable std::mutex               state_mutex_;
    std::optional<Position_solution> latest_position_;
    std::vector<Channel_snapshot>    latest_snapshots_;
    Receiver_status                  latest_status_;
    std::map<int, Almanac>       almanac_; // receiver-wide accumulated almanac (sv_key(con,prn) -> coarse orbit)

    // History subscriptions + published copies (guarded by state_mutex_; mutable for the const API).
    // Subscriptions are REFERENCE-COUNTED per SV key: several graph windows can watch the same SV, so a
    // key stays published until the LAST of them unsubscribes.
    mutable std::map<int, int>                        history_subscriptions_; // key -> refcount
    mutable std::map<int, Tracking_history::Snapshot> published_histories_;
};

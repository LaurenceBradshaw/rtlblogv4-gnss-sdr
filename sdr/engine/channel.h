#pragma once
#include <atomic>
#include <chrono>
#include <cmath>
#include <memory>
#include <mutex>
#include "acquisition.h"
#include "signal_aiding.h"
#include "constants.h"
#include "ephemeris.h"
#include "gps_l1ca_navigation.h"
#include "sample_buffer.h"
#include "signal.h"
#include "tracker.h"
#include "tracking_history.h"
#include "types.h"

enum class Channel_state
{
    IDLE,
    ACQUIRING,
    TRACKING
};

// A coherent, read-only snapshot of a channel's observable state, published by the owning worker
// at the end of each process() quantum and read by other threads (the Observation_engine on the
// main thread; later the GUI). Reading the live channel state across threads is a data race - the
// worker mutates the ephemeris / tracking / sample counters while a reader is part-way through -
// so EVERYTHING the readers need is captured here, under the channel's snapshot mutex, in one go.
struct Channel_snapshot
{
    Satellite_id  satellite_id   = 0;
    Constellation constellation  = Constellation::Unknown; // set from the channel on publish
    Code          code           = Code::CA;               // signal component (e.g. GPS L1CA vs L1Cd)
    Band          band           = Band::L1;               // RF band - so (con,sv) dedup can pick one component
    bool          has_observable = false;                  // valid ephemeris + TOW anchored at snapshot time

    Channel_state state     = Channel_state::IDLE; // acquisition / tracking lifecycle stage
    bool          has_lock  = false;               // carrier (PLL) lock; only meaningful while TRACKING
    double        cn0_db_hz = 0.0;                 // C/N0 measured at acquisition (dB-Hz); 0 until acquired

    Sample_index next_sample          = 0;   // channel's tracking position when snapshotted
    double       transmission_time_s  = 0.0; // SV transmit time (s) at next_sample
    double       carrier_doppler_hz   = 0.0; // tracked carrier Doppler at next_sample
    double       carrier_acceleration = 0.0; // d(Doppler)/dt, for projecting to the common sample
    double       wavelength_m         = 0.0; // c / carrier frequency (constant; cached for convenience)
    double       carrier_phase_cycles = 0.0; // accumulated carrier phase (cycles) at next_sample (Hatch smoothing)
    uint32_t     lock_session         = 0;   // ++ each (re-)acquisition; Hatch resets its arc when this changes
    int          carrier_lock_breaks  = 0;   // carrier-domain cycle-slip counter this arc (loss-of-lock events)

    Ephemeris eph;  // broadcast orbit/clock (base fields are all the orbit model needs)
    Iono      iono; // broadcast Klobuchar iono (iono.valid false unless this SV decoded SF4 p18)

    // Project the snapshotted SV transmit time forward/back to the common observation sample. The transmit
    // clock does NOT advance 1:1 with receiver time: from pr = (t_rx - t_tx)*c with d(pr)/dt_rx = range_rate,
    // dt_tx/dt_rx = 1 - range_rate/c (a receding SV stretches the signal -> transmit time advances slower).
    // local_s is RECEIVER seconds (rx_sample is the common observation sample = the slowest channel, so
    // faster channels project BACK; an L1Cd 10 ms-epoch channel can be ~10 ms stale), so scale it by that
    // ratio. Small (a few m at the staleness extremes) but correct - and it keeps L1CA/L1Cd consistent when
    // their next_sample grids differ. The slowest channel (local_s = 0) is unaffected.
    double transmit_time_at( Sample_index rx_sample, double sample_rate_hz ) const
    {
        const double local_s = ( static_cast<int64_t>( rx_sample ) - static_cast<int64_t>( next_sample ) ) / sample_rate_hz;
        const double range_rate_m_s = -wavelength_m * carrier_doppler_hz;
        return transmission_time_s + local_s * ( 1.0 - range_rate_m_s / constants::SPEED_OF_LIGHT_M_S );
    }

    // Pseudorange rate = -lambda * Doppler (textbook; the conj-wipe carrier_freq_ is already the
    // negated-physical-Doppler that makes this come out as the static receiver's ~0 m/s velocity).
    // The Doppler is taken AT next_sample and NOT projected to rx_sample: over the snapshot's
    // staleness the true range rate barely moves (<~0.1 m/s), but extrapolating by carrier_acceleration
    // amplifies a marginally-locked SV's noisy acceleration into a large prr error that can wreck the
    // velocity seed. carrier_acceleration is kept on the snapshot for display, not used here.
    double range_rate_at( Sample_index /*rx_sample*/, double /*sample_rate_hz*/ ) const
    {
        return -wavelength_m * carrier_doppler_hz;
    }

    // Accumulated carrier phase as a RANGE in metres at the common observation sample, sign-aligned to the
    // pseudorange (so a phase delta matches a pseudorange delta): = -lambda * (cycles projected to rx_sample
    // via the Doppler, like transmit_time_at). Only DELTAS are used (Hatch carrier-smoothing); the absolute
    // value carries the unknown integer-cycle ambiguity.
    double carrier_phase_range_m( Sample_index rx_sample, double sample_rate_hz ) const
    {
        const double local_s = ( static_cast<int64_t>( rx_sample ) - static_cast<int64_t>( next_sample ) ) / sample_rate_hz;
        return -wavelength_m * ( carrier_phase_cycles + carrier_doppler_hz * local_s );
    }
};

// A single satellite-tracking channel. Passive: it owns no thread. The Scheduler
// hands it to a Thread_pool worker as a work-item when it has data to process.
class Channel
{
public:
    // signal: the constellation/band this channel tracks (shared across channels of the
    //         same signal). Supplies the code replicas, nav decoder, and signal physics.
    // aiding: receiver-wide Doppler-recenter estimate, shared by all channels. This channel
    //         reads it to aim its acquisition search and reports its own Doppler back on lock.
    Channel(
        const Signal&       signal,
        Satellite_id        satellite_id,
        uint32_t            sample_rate_hz,
        Sample_buffer&      sample_buffer,
        Signal_aiding& aiding
    );

    // Not copyable or movable (holds FFTW plans, a Sample_buffer reference, etc.)
    Channel( const Channel& )            = delete;
    Channel& operator=( const Channel& ) = delete;

    // Read-only, called by the dispatcher thread
    // True if there is work this channel can do right now (data available, and not
    // in acquisition back-off). Only ever consulted while the channel is NOT
    // scheduled, so next_sample_/state_/retry_after_ are settled (published by the
    // owning worker's release of scheduled_).
    bool has_pending_work() const;
    bool is_acquiring() const
    {
        return state_ == Channel_state::ACQUIRING;
    }
    Sample_index next_sample() const
    {
        return next_sample_.load( std::memory_order_relaxed );
    }

    // Oldest sample this channel still needs - used to advance the buffer's oldest
    // index. A tracking channel needs from next_sample_ onward; an acquiring channel
    // only ever reads the newest block (it jumps there), so it never pins old data
    // (critically, a backed-off channel must NOT freeze the buffer at a stale index).
    Sample_index buffer_floor() const;

    // True once the stream has ended at `eof` and this channel has consumed as far as it can:
    // it cannot form another full block (one acq_block(), which also bounds a tracking epoch)
    // within acq_block() of the end. Uses the channel's OWN block size, so signals with
    // different code periods (GPS 1 ms vs Galileo 4 ms) each drain to their own limit - a
    // single shared threshold would strand the longer-block channels and stall shutdown.
    bool is_drained( Sample_index eof ) const
    {
        const Sample_index need = static_cast<Sample_index>( acq_block() );
        return buffer_floor() >= ( eof > need ? eof - need : Sample_index { 0 } );
    }

    // Single-owner flag: the dispatcher sets it before submitting a task; the worker
    // clears it when process() returns. Guarantees <=1 worker per channel at a time.
    std::atomic<bool>& scheduled()
    {
        return scheduled_;
    }

    // Called by the dispatcher (sole accessor - channel not yet scheduled) right before
    // an acquisition task is enqueued. Publishes the block the worker will jump to, so
    // that while the task sits queued/running and the dispatcher pins next_sample(), the
    // floor reflects where it WILL read (~ newest), not a value left stale by backoff.
    // process_acquisition() re-applies the same forward jump, keeping next_sample_
    // monotonic; without this the status briefly reads the stale (backed-off) index.
    void refresh_acquire_floor();

    // Called by a pool worker (the channel's sole owner while running)
    // Process one bounded quantum: one acquisition integrate (which may transition
    // straight into tracking), or a capped batch of tracking epochs.
    void process();

    bool has_observable() const; // valid eph + TOW anchored
    // Thread-safe, coherent copy of this channel's observable state (see Channel_snapshot). The
    // Observation_engine / GUI read THIS, never the live members the owning worker mutates.
    Channel_snapshot snapshot() const;
    // Thread-safe copy of this channel's graph history (I/Q heatmap + Doppler). The owning worker writes
    // history_ each epoch; this read is guarded by the history's own mutex.
    Tracking_history::Snapshot history_snapshot() const
    {
        return history_.snapshot();
    }
    // Thread-safe copy of the almanac entries this channel's decoder has gathered (PRN -> coarse orbit).
    // Published under snapshot_mutex_ by the owning worker (in publish_snapshot); the Receiver aggregates
    // these across channels for acquisition aiding. Constellation-wide (any tracked SV fills it in).
    std::map<int, Almanac> almanac_snapshot() const
    {
        std::lock_guard<std::mutex> lk( snapshot_mutex_ );
        return published_almanac_;
    }
    Satellite_id satellite_id() const
    {
        return satellite_id_;
    }
    const Ephemeris& ephemeris() const
    {
        return navigation_->ephemeris();
    }
    const Signal& signal() const
    {
        return signal_;
    }

private:
    void process_acquisition();
    void process_tracking();

    // Capture the current observable state into snapshot_ under snapshot_mutex_. Called by the
    // owning worker at the end of process() (it is the sole mutator of the channel's state).
    void publish_snapshot();

    const double current_transmission_time_s() const;
    // Acquisition block = one FFT window (acq_fft_factor * one code period); must match
    // Acquisition_engine::m_ so integrate() has exactly the samples its FFT consumes. ROUND (not
    // truncate) the samples-per-period: the code replica is round(rate*period) samples, so at a
    // fractional sample rate (e.g. 4.166 MHz from --decimate 6) truncating here left the block one row
    // short and integrate() read past it -> segfault.
    size_t acq_block() const
    {
        return static_cast<size_t>( signal_.params().acq_fft_factor ) * static_cast<size_t>( std::llround( EPOCH_SAMPLES ) );
    }

    const Signal&                signal_;
    Signal_aiding&          aiding_;
    Sample_buffer&               sample_buffer_;
    Satellite_id                 satellite_id_;
    Channel_state                state_;
    std::atomic<Sample_index>    next_sample_; // written by owning worker, read by dispatcher
    Acquisition_engine           acquisition_;
    std::unique_ptr<Tracker>     tracking_;
    std::unique_ptr<Nav_decoder> navigation_;
    uint32_t                     sample_rate_hz_;
    const double                 EPOCH_SAMPLES; // samples per 1 ms code period
    bool                         lock_reported_;
    int                          acq_attempts_       = 0;   // full acquisition attempts (DUMP_ACQ telemetry)
    double                       acq_cn0_db_hz_      = 0.0; // C/N0 from the acquisition that locked this channel
    Sample_index                 track_start_sample_ = 0;   // next_sample_ at tracking handoff (frame-sync timeout)
    // Cross-code/-frequency Doppler aiding throttles (see process_tracking): when this channel last published
    // its measured Doppler (donor) and last re-centered its NCO from a sibling (recipient).
    Sample_index     last_aid_report_sample_ = 0;
    Sample_index     last_nudge_sample_      = 0;
    Sample_index     last_clock_ff_sample_   = 0; // throttle of the common-clock-drift feedforward (Fix 3)
    uint32_t         track_session_          = 0; // ++ on each (re-)acquisition; published so Hatch resets its arc
    Tracking_history history_; // per-epoch graph history (prompt I/Q, ...)

    // Acquisition back-off (#1): after a full attempt finds nothing, don't retry
    // until retry_after_, with exponentially growing cooldown. Written by the owning
    // worker; read by has_pending_work() - synchronised by the scheduled_ flag.
    std::chrono::steady_clock::time_point retry_after_ {};
    std::chrono::milliseconds             backoff_ { 0 };

    std::atomic<bool> scheduled_ { false };

    // Coherent observable snapshot: written by the owning worker (publish_snapshot), read by other
    // threads (snapshot()). A mutex - not a lock-free seqlock - because publish/read are infrequent
    // (per quantum / per tick) and brief (a small struct copy), so contention is negligible.
    Channel_snapshot       snapshot_;
    std::map<int, Almanac> published_almanac_; // decoded almanac, published under snapshot_mutex_
    mutable std::mutex     snapshot_mutex_;
};

#include "channel.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "constants.h"
#include "logging.h"
#include "timing.h"

namespace
{
// Cap on tracking epochs processed in one quantum, so a lagging channel can't hog a
// worker. The real-time gate normally bounds a quantum to ~the dispatcher interval
// (a handful of epochs); this is just a safety ceiling.
constexpr int MAX_EPOCHS_PER_QUANTUM = 64;

// Acquisition back-off schedule: after a failed full attempt, wait this long before
// retrying, doubling each time up to the cap. Satellites that aren't in view go quiet.
constexpr std::chrono::milliseconds ACQ_BACKOFF_START { 500 };
constexpr std::chrono::milliseconds ACQ_BACKOFF_MAX { 5000 };

// Cross-code / cross-frequency Doppler aiding (see process_tracking).
constexpr double AID_REPORT_PERIOD_S = 0.02; // throttle a donor's sibling-Doppler reports (Doppler drifts slowly)
constexpr double NUDGE_MIN_TRACK_S   = 0.5;  // don't re-center a freshly handed-off channel still pulling in normally
constexpr double NUDGE_DWELL_S       = 1.0;  // min gap between re-centers, so the loop gets time to pull in
constexpr double NUDGE_MIN_DELTA_HZ  = 100.0; // only re-center if the aided Doppler differs enough to be worth it
constexpr double CLOCK_FF_PERIOD_S   = 0.1;   // throttle the common-clock-drift feedforward (PVT-rate state)
} // namespace

Channel::Channel(
    const Signal&       signal,
    Satellite_id        satellite_id,
    uint32_t            sample_rate_hz,
    Sample_buffer&      sample_buffer,
    Signal_aiding& aiding
)
    : signal_( signal ),
      aiding_( aiding ),
      sample_buffer_( sample_buffer ),
      satellite_id_( satellite_id ),
      state_( Channel_state::ACQUIRING ),
      next_sample_( 0 ),
      acquisition_(
          signal.code_samples( satellite_id, sample_rate_hz ), satellite_id, static_cast<double>( sample_rate_hz ),
          signal.params(), aiding
      ),
      tracking_( signal.make_tracker( satellite_id, static_cast<double>( sample_rate_hz ) ) ),
      navigation_( signal.make_nav_decoder( satellite_id ) ),
      sample_rate_hz_( sample_rate_hz ),
      EPOCH_SAMPLES( static_cast<double>( sample_rate_hz ) * signal.params().code_period_s ),
      lock_reported_( false )
{
}

// buffer_floor
Sample_index Channel::buffer_floor() const
{
    if( state_ == Channel_state::TRACKING )
    {
        return next_sample_.load( std::memory_order_relaxed );
    }

    // Acquiring: needs only the newest block (process_acquisition jumps there), so it
    // pins nothing older than newest+1-ACQ_BLOCK, regardless of its frozen next_sample_.
    const Sample_index newest = sample_buffer_.newest_valid_index();
    const Sample_index blk    = static_cast<Sample_index>( acq_block() );
    return ( newest + 1 >= blk ) ? newest + 1 - blk : 0;
}

// refresh_acquire_floor
void Channel::refresh_acquire_floor()
{
    // buffer_floor() while ACQUIRING returns newest+1-ACQ_BLOCK (the jump target).
    const Sample_index floor = buffer_floor();
    if( floor > next_sample_.load( std::memory_order_relaxed ) )
    {
        next_sample_.store( floor, std::memory_order_relaxed );
    }
}

// has_pending_work
bool Channel::has_pending_work() const
{
    if( state_ == Channel_state::ACQUIRING )
    {
        if( std::chrono::steady_clock::now() < retry_after_ )
        {
            return false; // backing off - not in view, don't burn a worker on it
        }
        // Almanac horizon skip: if the broadcast almanac + the current fix place this SV below the horizon,
        // don't search it at all - just don't schedule the channel. Re-evaluated each poll, so it resumes
        // the moment the prediction says it has risen. Unknown SVs (no prediction) stay searchable.
        if( !aiding_.searchable( signal_.params().constellation, static_cast<int>( satellite_id_ ) ) )
        {
            return false;
        }
        const Sample_index newest = sample_buffer_.newest_valid_index();
        return newest + 1 >= static_cast<Sample_index>( acq_block() );
    }

    // TRACKING - is the next epoch's full block real-time-available?
    const Sample_index ns   = next_sample_.load( std::memory_order_relaxed );
    const int          need = tracking_->compute_samples_needed();
    return sample_buffer_.is_available( ns + static_cast<Sample_index>( need ) - 1 );
}

// process
// One bounded quantum. Sole owner of the channel while this runs (guaranteed by the
// Scheduler's scheduled_ flag).
void Channel::process()
{
    if( state_ == Channel_state::ACQUIRING )
    {
        process_acquisition(); // one integrate; may transition to TRACKING
    }

    // Drain the tracking epochs that are available this round (covers both normal
    // tracking and continuing straight after a successful acquire in the same quantum).
    if( state_ == Channel_state::TRACKING )
    {
        for( int i = 0; i < MAX_EPOCHS_PER_QUANTUM; i++ )
        {
            const Sample_index before = next_sample_.load( std::memory_order_relaxed );
            process_tracking();
            if( next_sample_.load( std::memory_order_relaxed ) == before )
            {
                break; // no progress -> no more data available this round
            }
        }
    }

    // Publish a coherent observable snapshot for cross-thread readers (Observation_engine / GUI).
    publish_snapshot();
}

bool Channel::has_observable() const
{
    return navigation_->ephemeris().valid && navigation_->tow_anchored();
}

// Capture this channel's observable state into snapshot_. Runs on the owning worker (the sole
// mutator of the channel), so reading the live members here is race-free; only the final copy into
// snapshot_ is locked, against concurrent readers. The transmit-time / range-rate PROJECTION to a
// common observation sample now lives in Channel_snapshot (computed by the reader from these
// captured values), so it is coherent with the rest of the snapshot.
void Channel::publish_snapshot()
{
    Channel_snapshot s;
    s.satellite_id         = satellite_id_;
    s.constellation        = navigation_->ephemeris().constellation;
    s.code                 = signal_.params().code;
    s.band                 = signal_.params().band;
    s.has_observable       = has_observable();
    s.next_sample          = next_sample_.load( std::memory_order_relaxed );
    s.transmission_time_s  = current_transmission_time_s();
    s.carrier_doppler_hz   = tracking_->get_carrier_doppler_hz();
    s.carrier_acceleration = tracking_->get_carrier_acceleration();
    s.carrier_phase_cycles = tracking_->get_carrier_phase_cycles();
    s.lock_session         = track_session_;
    s.carrier_lock_breaks  = tracking_->get_carrier_lock_breaks();
    s.state                = state_;
    s.has_lock             = ( state_ == Channel_state::TRACKING ) && tracking_->has_lock();
    // Live tracking C/N0 (M2M4) once it has a window; fall back to the acquisition C/N0 until then.
    s.cn0_db_hz = ( state_ == Channel_state::TRACKING && tracking_->get_cn0_db_hz() > 0.0 ) ? tracking_->get_cn0_db_hz()
                                                                                            : acq_cn0_db_hz_;
    s.wavelength_m         = constants::SPEED_OF_LIGHT_M_S / signal_.params().carrier_freq_hz;
    s.eph                  = navigation_->ephemeris(); // copies the base Ephemeris (all the orbit model needs)
    s.iono                 = *navigation_->iono();     // broadcast iono (iono.valid false unless this SV decoded it)

    std::lock_guard<std::mutex> lk( snapshot_mutex_ );
    snapshot_          = s;
    published_almanac_ = navigation_->almanac(); // publish under the lock; the worker owns navigation_
}

Channel_snapshot Channel::snapshot() const
{
    std::lock_guard<std::mutex> lk( snapshot_mutex_ );
    return snapshot_;
}

const double Channel::current_transmission_time_s() const
{
    // Get the coarse 6-second subframe master anchor
    const double t_tow = navigation_->ephemeris().tow;

    // Whole code periods (epochs) since the eph.tow anchor. ms_since_tow_update() is an epoch_count_ delta
    // (1 ms/epoch L1CA, 10 ms/epoch L1C), so it is RESET-ROBUST (it does not depend on the per-bit/per-ms nav
    // counters that get reset on a TOW update - the reason the old t_bits/t_ms terms were abandoned). It is
    // whole-epoch only; the sub-epoch part is recovered to ~+/-1 sample by the snapshot's whole-sample
    // next_sample projection (transmit_time_at). Adding tracking_->get_fractional_chip_time() here does NOT
    // improve this and breaks it: that getter returns elapsed-WITHIN-period (~= a FULL period right after an
    // epoch completes), which ms_since_tow_update already counts -> double-counts one period (verified: +9 ms
    // L1CA-vs-L1C, position destroyed); its complement (the sub-chip residual) only WORSENED the fix. L1CA
    // PVT is already accurate (~tens of m) without any t_frac, so the dropped fractional code phase is a
    // minor effect, not the dual-code-seed cause. See the backlog "Sub-ms transmit-time precision" note.
    const double t_since_tow_update = navigation_->ms_since_tow_update() * 0.001;

    // Sub-sample piece 1 (jitter): add the live code-NCO sub-sample offset of next_sample from the code
    // boundary, so t_tx isn't quantised to the nearest whole sample. Validated via the TRUTH_LLH harness
    // (jitter target 42 m -> ~m). A per-channel acquisition-rounding BIAS remains (the acq seed, piece 2).
    return t_tow + t_since_tow_update + tracking_->code_phase_offset_s();
}

// process_acquisition
// mirrors GNSS-SDRLIB sdraccuisition() - see acquisition.cpp comments for detail.
void Channel::process_acquisition()
{
    const size_t ACQ_BLOCK = acq_block();
    Sample_index ns        = next_sample_.load( std::memory_order_relaxed );

    // Acquisition doesn't need continuous data: jump to the most recent block so the
    // buffer's oldest index can advance and free space.
    const Sample_index newest = sample_buffer_.newest_valid_index();
    if( newest + 1 >= static_cast<Sample_index>( ACQ_BLOCK ) )
        ns = std::max( ns, newest + 1 - static_cast<Sample_index>( ACQ_BLOCK ) );

    // Check the LAST sample in the block (next_sample_+ACQ_BLOCK-1), not one past end.
    if( !sample_buffer_.is_available( ns + static_cast<Sample_index>( ACQ_BLOCK ) - 1 ) )
    {
        next_sample_.store( ns, std::memory_order_relaxed ); // publish the jump
        return;
    }

    const Sample_block block = sample_buffer_.at( ns, ACQ_BLOCK );

    bool acquired;
    {
        timing::Scoped t( timing::Bucket::Acquisition );
        acquired = acquisition_.integrate( block );
    }

    if( acquired )
    {
        const Acquisition_result& r = acquisition_.result();
        ++acq_attempts_;

        if( std::getenv( "DUMP_ACQ" ) != nullptr ) // per-SV acq telemetry (sensitivity investigation)
        {
            logging::log(
                logging::Level::Info,
                fmt::format(
                    "ACQDIAG {} PRN {:2d} attempt {:3d} PASS ratio {:.2f} (thr {:.1f}) C/N0 {:.1f} dopp {:+.0f} code {:.0f}",
                    signal_.params().name, satellite_id_, acq_attempts_, r.metric, Acquisition_engine::ACQTH,
                    r.cn0_db_hz, r.doppler_hz, r.code_phase
                )
            );
        }

        logging::log(
            logging::Level::Info,
            fmt::format(
                "Acquisition - {} PRN {:2d}  Doppler {:+.0f} Hz  Code phase {:.0f} samples  Ratio {:.1f}  C/N0 {:.1f} dB-Hz",
                signal_.params().name,
                satellite_id_,
                r.doppler_hz,
                r.code_phase,
                r.metric,
                r.cn0_db_hz
            )
        );

        // Feed this satellite's Doppler into the receiver-wide estimate so later acquisitions
        // recenter their search on the common clock offset (pre-PVT bridge).
        aiding_.report( r.doppler_hz, signal_.params().carrier_freq_hz );

        acq_cn0_db_hz_ = r.cn0_db_hz;
        tracking_->initialise( r );
        ++track_session_; // new lock arc - the carrier phase accumulator just reset, so bump the id Hatch keys on
        history_.reset(); // fresh graph history for this lock
        state_ = Channel_state::TRACKING;

        // Align next_sample_ to the detected code boundary (buffloc = acq_start + acqcodei).
        ns += static_cast<Sample_index>( std::round( r.code_phase ) );
        next_sample_.store( ns, std::memory_order_relaxed );
        track_start_sample_ = ns; // for the frame-sync timeout (see process_tracking)

        acquisition_.reset();
        return;
    }

    ns += static_cast<Sample_index>( EPOCH_SAMPLES );
    next_sample_.store( ns, std::memory_order_relaxed );

    if( acquisition_.intg_count() >= acquisition_.target_integrations() )
    {
        ++acq_attempts_;
        if( std::getenv( "DUMP_ACQ" ) != nullptr ) // per-SV acq telemetry: log FAILED attempts (the misses)
        {
            const Acquisition_result& r = acquisition_.result(); // populated by the last check_acquisition()
            logging::log(
                logging::Level::Info,
                fmt::format(
                    "ACQDIAG {} PRN {:2d} attempt {:3d} FAIL ratio {:.2f} (thr {:.1f}) C/N0 {:.1f} dopp {:+.0f} code {:.0f}",
                    signal_.params().name, satellite_id_, acq_attempts_, r.metric, Acquisition_engine::ACQTH,
                    r.cn0_db_hz, r.doppler_hz, r.code_phase
                )
            );
        }
        // Full attempt finished with no hit -> back off before retrying so satellites
        // that aren't in view stop consuming workers (#1).
        acquisition_.reset();
        backoff_     = ( backoff_.count() == 0 ) ? ACQ_BACKOFF_START : std::min( backoff_ * 2, ACQ_BACKOFF_MAX );
        retry_after_ = std::chrono::steady_clock::now() + backoff_;
    }
}

// process_tracking
// mirrors GNSS-SDRLIB sdrthread() tracking section:
//   bufflocnow = sdrtracking(sdr, buffloc, cnt)
//   if flagtrk: cumsumcorr, pll, dll, clearcumsumcorr, buffloc += samples_consumed
void Channel::process_tracking()
{
    const Sample_index ns         = next_sample_.load( std::memory_order_relaxed );
    const int          curr_nsamp = tracking_->compute_samples_needed();

    // Loss-of-lock / divergence guard. A healthy epoch is ~one code period
    // (EPOCH_SAMPLES); if the code NCO has run away, curr_nsamp goes absurd. Drop the
    // channel back to acquisition so it stops tracking garbage (and, importantly, stops
    // pinning the sample buffer at a stale index, which stalls the whole pipeline).
    const int nominal = static_cast<int>( EPOCH_SAMPLES );
    if( curr_nsamp < nominal / 2 || curr_nsamp > nominal * 2 )
    {
        logging::log(
            logging::Level::Info,
            fmt::format( "Tracking diverged (nsamp={}), re-acquiring PRN {:2d}", curr_nsamp, satellite_id_ )
        );
        state_         = Channel_state::ACQUIRING;
        lock_reported_ = false;
        acquisition_.reset();
        navigation_ = signal_.make_nav_decoder( satellite_id_ ); // reset bit-sync / frame state
        return;                                                  // next_sample_ unchanged; process() will re-acquire
    }

    // Frame-sync timeout. A channel that locked onto a false / cross-correlation peak gets bit sync on
    // the garbage but never frame sync (subframe preamble / I/NAV page). If we have tracked past the
    // time a clean lock decodes its first frame and still have no frame sync, the lock is bad: drop
    // back to acquire. Re-acquisition jumps to fresh data and almost always lands on the true peak,
    // so the SV is recovered instead of tracking noise for the whole run.
    const double track_s = static_cast<double>( ns - track_start_sample_ ) / static_cast<double>( sample_rate_hz_ );
    if( track_s > signal_.params().frame_sync_timeout_s && !navigation_->frame_synced() )
    {
        logging::log(
            logging::Level::Info, fmt::format( "No frame sync after {:.0f}s, re-acquiring PRN {:2d}", track_s, satellite_id_ )
        );
        state_         = Channel_state::ACQUIRING;
        lock_reported_ = false;
        acquisition_.reset();
        navigation_ = signal_.make_nav_decoder( satellite_id_ );
        return;
    }

    // Check the last sample (n-1), not one past end.
    if( !sample_buffer_.is_available( ns + static_cast<Sample_index>( curr_nsamp ) - 1 ) )
    {
        return;
    }

    const Sample_block block = sample_buffer_.at( ns, static_cast<size_t>( curr_nsamp ) );

    // Per-epoch order mirrors GNSS-SDRLIB: correlate -> navigation -> loop filters.
    // Navigation MUST run before run_loops() so that bit_sync / sw_loop reflect THIS
    // epoch (else the 10 ms coherent integration straddles data-bit boundaries).
    Tracking_output out;
    {
        timing::Scoped t( timing::Bucket::TrackCorrelate );
        out = tracking_->correlate_epoch( block );
    }

    if( !out.valid )
    {
        return;
    }

    navigation_->process( out.prompt_i, out.prompt_i_prev );

    {
        timing::Scoped t( timing::Bucket::TrackLoops );
        tracking_->run_loops( navigation_->bit_sync_found(), navigation_->sw_loop(), satellite_id_ );
    }

    if( tracking_->has_lock() && !lock_reported_ )
    {
        lock_reported_ = true;
        logging::log(
            logging::Level::Info, fmt::format( "Tracking - {} PRN {:2d}  LOCK", signal_.params().name, satellite_id_ )
        );
    }
    if( !tracking_->has_lock() && lock_reported_ )
    {
        lock_reported_ = false;
        logging::log(
            logging::Level::Info, fmt::format( "Tracking - {} PRN {:2d}  lock lost", signal_.params().name, satellite_id_ )
        );
    }

    // ---- Coupled-clock feedforward (Fix 3) ----
    // Feed the receiver-wide common clock drift (PVT) forward into our NCOs so the loops don't chase the
    // coupled-oscillator drift ramp reactively (helps a weak channel hold lock under thermal drift). No-op
    // pre-PVT and when unchanged; throttled (the common drift updates only at the PVT rate).
    if( ns - last_clock_ff_sample_ >= static_cast<Sample_index>( CLOCK_FF_PERIOD_S * sample_rate_hz_ ) )
    {
        tracking_->apply_common_clock_drift( aiding_.clock_fraction() );
        last_clock_ff_sample_ = ns;
    }

    // ---- Cross-code / cross-frequency per-SV Doppler aiding ----
    // Keyed by (constellation,PRN) only, so an SV's signals on different codes/bands (e.g. GPS L1CA and
    // L1Cd, or L1 and L5) share one entry: a clean lock on one aids a struggling sibling on another.
    const Constellation con = signal_.params().constellation;
    const int           prn = static_cast<int>( satellite_id_ );

    // DONOR: once frame-synced (a genuinely correct lock - safe against the cross-correlation false locks
    // that the frame-sync timeout evicts, which have high C/N0 + carrier lock but never frame-sync), publish
    // this SV's measured Doppler as df/f. Throttled. get_carrier_doppler_hz() is the tracker's negated
    // carrier_freq_ (conj-wipe), so negate it back to the PHYSICAL Doppler the aiding stores.
    if( navigation_->frame_synced()
        && ns - last_aid_report_sample_ >= static_cast<Sample_index>( AID_REPORT_PERIOD_S * sample_rate_hz_ ) )
    {
        aiding_.report_sv_doppler( con, prn, -tracking_->get_carrier_doppler_hz(), signal_.params().carrier_freq_hz );
        last_aid_report_sample_ = ns;
    }

    // RECIPIENT: a channel that is TRACKING but has neither carrier lock nor frame sync is struggling. If a
    // sibling has measured this SV's Doppler, re-center our carrier NCO on it (scaled to our carrier; df/f is
    // band-agnostic) and then leave the loop alone for NUDGE_DWELL_S to pull in (sticky). Gated so we never
    // disturb a freshly handed-off channel still pulling in normally, and only when the aided Doppler differs
    // enough to be worth the reset.
    double       sib_fraction = 0.0;
    const double track_s2     = static_cast<double>( ns - track_start_sample_ ) / sample_rate_hz_;
    if( !tracking_->has_lock() && !navigation_->frame_synced() && track_s2 > NUDGE_MIN_TRACK_S
        && ns - last_nudge_sample_ >= static_cast<Sample_index>( NUDGE_DWELL_S * sample_rate_hz_ )
        && aiding_.sv_doppler_fraction( con, prn, sib_fraction ) )
    {
        const double target_hz  = sib_fraction * signal_.params().carrier_freq_hz; // physical Doppler
        const double current_hz = -tracking_->get_carrier_doppler_hz();            // physical Doppler
        if( std::abs( target_hz - current_hz ) > NUDGE_MIN_DELTA_HZ )
        {
            tracking_->steer_carrier_doppler( target_hz ); // physical Doppler in; tracker negates internally
            last_nudge_sample_ = ns;
            logging::log(
                logging::Level::Info,
                fmt::format(
                    "Tracking aid - {} PRN {:2d} carrier re-centered {:+.0f} -> {:+.0f} Hz from sibling lock",
                    signal_.params().name, satellite_id_, current_hz, target_hz
                )
            );
        }
    }

    next_sample_.store( ns + static_cast<Sample_index>( out.samples_consumed ), std::memory_order_relaxed );

    // Per-epoch history for the GUI graphs (prompt I/Q constellation; Doppler later). Recorded always
    // so a graph opened later already has data; cheap (one cell + occasional fade, under its own mutex).
    history_.record(
        static_cast<double>( ns ) / sample_rate_hz_,
        tracking_->get_carrier_doppler_hz(),
        tracking_->get_prompt_i(),
        tracking_->get_prompt_q()
    );
}

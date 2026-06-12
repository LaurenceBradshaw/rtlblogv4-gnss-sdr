#include "channel.h"
#include <algorithm>
#include <cmath>
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
} // namespace

Channel::Channel(
    const Signal&       signal,
    Satellite_id        satellite_id,
    uint32_t            sample_rate_hz,
    Sample_buffer&      sample_buffer,
    Acquisition_aiding& aiding
)
    : signal_( signal ),
      aiding_( aiding ),
      sample_buffer_( sample_buffer ),
      satellite_id_( satellite_id ),
      state_( Channel_state::ACQUIRING ),
      next_sample_( 0 ),
      acquisition_(
          signal.code_samples( satellite_id, sample_rate_hz ), static_cast<double>( sample_rate_hz ), signal.params(), aiding
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
    s.has_observable       = has_observable();
    s.next_sample          = next_sample_.load( std::memory_order_relaxed );
    s.transmission_time_s  = current_transmission_time_s();
    s.carrier_doppler_hz   = tracking_->get_carrier_doppler_hz();
    s.carrier_acceleration = tracking_->get_carrier_acceleration();
    s.state                = state_;
    s.has_lock             = ( state_ == Channel_state::TRACKING ) && tracking_->has_lock();
    // Live tracking C/N0 (M2M4) once it has a window; fall back to the acquisition C/N0 until then.
    s.cn0_db_hz = ( state_ == Channel_state::TRACKING && tracking_->get_cn0_db_hz() > 0.0 ) ? tracking_->get_cn0_db_hz()
                                                                                            : acq_cn0_db_hz_;
    s.wavelength_m         = 299792458.0 / signal_.params().carrier_freq_hz;
    s.eph                  = navigation_->ephemeris(); // copies the base Ephemeris (all the orbit model needs)
    s.iono                 = *navigation_->iono();     // broadcast iono (iono.valid false unless this SV decoded it)

    std::lock_guard<std::mutex> lk( snapshot_mutex_ );
    snapshot_ = s;
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

    const double t_since_tow_update =
        navigation_->ms_since_tow_update() * 0.001; // TODO: Figure out if this includes t_bits, t_ms and t_frac or not
    // Note: This is used so that in the event that a satellites ephemeris is not updated while t_bits, t_ms and t_frac have all
    // had their values reset. There could be a better way to do this

    // // Add the elapsed whole navigation bits (scaled to seconds)
    // const double t_bits = static_cast<double>( navigation_->get_current_bit_index() ) *
    //                       ( static_cast<double>( signal_.params().nav_bit_ms ) / 1000.0 );

    // // Add the elapsed milliseconds inside the current bit (explicitly 0.001s per tick)
    // const double t_ms = static_cast<double>( navigation_->get_current_ms_tick() ) * 0.001;

    // // Add the precise sub-millisecond chip fraction from the NCO tracker
    // const double t_frac = tracking_->get_fractional_chip_time();

    // Combine all components into the ultra-precise time of transmission
    return t_tow + t_since_tow_update; // + t_bits + t_ms + t_frac;
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

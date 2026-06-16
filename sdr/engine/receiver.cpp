#include "receiver.h"
#include <chrono>
#include <cmath>
#include <fmt/format.h>
#include <thread>
#include "geodesy.h"
#include "logging.h"
#include "orbit.h"
#include "rtlsdr_device.h"
#include "sample_buffer.h"
#include "scheduler.h"
#include "signal.h"
#include "thread_pool.h"
#include "types.h"

namespace
{
// Round n up to the next power of two.
size_t next_power_of_two( size_t n )
{
    size_t p = 1;
    while( p < n )
    {
        p <<= 1;
    }
    return p;
}

constexpr uint32_t GNSS_L1_HZ = 1575420000U; // GPS L1 C/A + Galileo E1 band centre

// Build the Signal object for one selection entry (make_signal per (constellation, band, code)).
// GPS anchors the fix + receiver clock; Galileo adds SVs via the EKF inter-system-bias state.
std::unique_ptr<Signal> make_signal_for( const Signal_selection& sel )
{
    return make_signal( sel.id.constellation, sel.id.band, sel.id.code );
}
} // namespace

Receiver::Receiver( Receiver_config config )
    : config_( std::move( config ) )
{
    // Resolve the signal selection once (empty config -> the default GPS + Galileo set), so the
    // constructor's satellite enumeration and setup()'s channel creation always agree.
    signal_selection_ = config_.selected_signals.empty() ? default_signal_selection() : config_.selected_signals;

    enumerate_configured_sats();
}

void Receiver::enumerate_configured_sats()
{
    // Enumerate the configured satellites, so the GUI can list them all before Start. Honour each
    // signal's PRN allowlist so the list matches the channels setup() will actually create.
    configured_sats_.clear();
    for( const Signal_selection& sel : signal_selection_ )
    {
        const auto sig            = make_signal_for( sel );
        const auto [sv_lo, sv_hi] = sig->sv_range();
        for( int sv = sv_lo; sv <= sv_hi; ++sv )
        {
            if( prn_selected( sel.prns, sv ) )
            {
                configured_sats_.push_back( { sig->params().constellation, sv, sig->params().code } );
            }
        }
    }
}

void Receiver::update_acquisition_predictions( const Position_solution& fix )
{
    // GPS almanac only for now. The receiver-wide almanac_ was just refreshed on this same (run) thread,
    // so reading it here without the lock is safe (the lock only guards GUI readers).
    constexpr double C        = 299792458.0;
    constexpr double EL_MASK  = -2.0 * 3.14159265358979 / 180.0; // skip only SVs clearly below the horizon
    const Ecef       user { fix.ecef_x_m, fix.ecef_y_m, fix.ecef_z_m };
    const double     t = obs_engine_.reception_time_s(); // current GPS time of week (s)

    for( const Configured_satellite& sat : configured_sats_ )
    {
        if( sat.constellation != Constellation::Gps )
        {
            continue;
        }
        const auto it = almanac_.find( sat.prn );
        if( it == almanac_.end() || !it->second.valid )
        {
            continue; // no almanac for this SV yet -> leave it searchable / common-mode aided
        }
        const Ephemeris eph = it->second.as_ephemeris();
        const Ecef      sp  = orbit::satellite_ecef_pos( eph, t, Constellation::Gps );
        const Ecef      sv  = orbit::satellite_ecef_vel( eph, t, Constellation::Gps );
        double          el = 0.0, az = 0.0;
        look_angles( user, sp, el, az );
        // Predicted line-of-sight Doppler as a fraction of carrier: range-rate = sat_vel . LOS (static rx),
        // physical Doppler fraction = -range_rate/c (matches Channel_snapshot::range_rate_at's convention).
        const double dx = sp.x - user.x, dy = sp.y - user.y, dz = sp.z - user.z;
        const double r  = std::sqrt( dx * dx + dy * dy + dz * dz );
        const double rr = ( r > 0.0 ) ? ( sv.x * dx + sv.y * dy + sv.z * dz ) / r : 0.0;
        aiding_.set_prediction( Constellation::Gps, sat.prn, el > EL_MASK, -rr / C );
    }
}

void Receiver::set_signal_selection( std::vector<Signal_selection> selection )
{
    if( running_.load() )
    {
        return; // takes effect only on the next run; ignore a change mid-run
    }
    // Set exactly what was asked (empty = search nothing) so the GUI's live preview is faithful when all
    // signals are unchecked - the default set is only the CONSTRUCTOR's fallback for an unspecified config.
    signal_selection_ = std::move( selection );
    enumerate_configured_sats();
}

void Receiver::set_source_params( Source_params params )
{
    if( running_.load() )
    {
        return; // takes effect only on the next run; ignore a change mid-run
    }
    config_.use_rtlsdr     = params.use_rtlsdr;
    config_.file_path      = std::move( params.file_path );
    config_.format         = params.format;
    config_.sample_rate_hz = params.sample_rate_hz;
    config_.device_index   = params.device_index;
    config_.gain_db        = params.gain_db;
    config_.decimation     = std::max( 1u, params.decimation );
}

Receiver::~Receiver() = default;

void Receiver::setup()
{
    // Start from a clean slate so a restart (GUI Start after Stop) runs afresh: drop anything a
    // previous run left, and reset the stateful pieces that teardown() does not own (the EKF, the
    // observation clock anchor, the acquisition aiding). An IQ file replays from the start because a
    // brand-new Iq_file_device is created below.
    signals_.clear();
    channels_.clear();
    channel_ptrs_.clear();
    aiding_.reset();
    obs_engine_ = Observation_engine {};
    pos_solver_ = Position_solver {};
    clear_published_state(); // a fresh run starts from a clean GUI state (drops any frozen EOF state)

    // The SOURCE runs at config_.sample_rate_hz; the rest of the pipeline runs at the (optionally
    // decimated) PROCESSING rate. The Fir_decimator sits in the streaming callback (device -> decimator
    // -> buffer), so the buffer, channels, scheduler and timing all use the processing rate.
    const uint32_t decim         = std::max( 1u, config_.decimation );
    const uint32_t native_rate   = config_.sample_rate_hz;
    const uint32_t sample_rate_hz = native_rate / decim; // processing rate
    decimator_                    = decim > 1 ? std::make_unique<Fir_decimator>( static_cast<int>( decim ) ) : nullptr;

    // Buffer must hold at least 2 s of samples so that push()'s 1 s sleep
    // (see sample_buffer.cpp) always wakes up to find meaningful space freed.
    const size_t buffer_capacity = next_power_of_two( 2 * static_cast<size_t>( sample_rate_hz ) );
    sample_buffer_ = std::make_unique<Sample_buffer>( buffer_capacity, static_cast<double>( sample_rate_hz ) );

    // Source: live RTL-SDR or recorded file, behind the Stream_device interface.
    std::string source_desc;
    if( config_.use_rtlsdr )
    {
        auto rtl = std::make_unique<Rtlsdr_device>( config_.device_index );
        rtl->set_sample_rate_hz( native_rate );
        rtl->set_centre_freq_hz( GNSS_L1_HZ );
        if( config_.gain_db >= 0.0 )
        {
            rtl->set_gain_tenths_db( static_cast<int>( config_.gain_db * 10.0 ) );
        }
        else
        {
            rtl->set_agc( true );
        }
        source_desc = fmt::format( "RTL-SDR device {} (centre {} Hz)", config_.device_index, GNSS_L1_HZ );
        device_     = std::move( rtl );
    }
    else
    {
        device_     = std::make_unique<Iq_file_device>( config_.file_path, native_rate, config_.format );
        source_desc = config_.file_path;
    }

    // Signals to search for. The objects must outlive the channels (channels hold a const reference).
    // signals_[i] pairs with signal_selection_[i] (same order), so each uses its own PRN allowlist.
    for( const Signal_selection& sel : signal_selection_ )
    {
        signals_.push_back( make_signal_for( sel ) );
    }

    // Passive channels (no per-channel thread): one per (signal, selected SV). A Thread_pool of
    // N workers does the correlation work; the Scheduler hands ready channels to it.
    for( size_t i = 0; i < signals_.size(); ++i )
    {
        const Signal&             sig         = *signals_[i];
        const std::set<int>&      prn_filter  = signal_selection_[i].prns;
        const auto [sv_lo, sv_hi]             = sig.sv_range();
        for( int sv = sv_lo; sv <= sv_hi; sv++ )
        {
            if( !prn_selected( prn_filter, sv ) )
            {
                continue;
            }
            if( !sig.broadcasts( static_cast<Satellite_id>( sv ) ) )
            {
                continue; // SV doesn't transmit this signal - don't waste a channel that can never lock
            }
            channels_.push_back( std::make_unique<Channel>(
                sig, static_cast<Satellite_id>( sv ), sample_rate_hz, *sample_buffer_, aiding_
            ) );
        }
    }

    std::string sig_list;
    for( const Signal_selection& s : signal_selection_ )
    {
        const std::string prns = s.prns.empty() ? "all" : fmt::format( "{}", s.prns.size() );
        sig_list += fmt::format(
            "{}{}:{}({} PRNs)", sig_list.empty() ? "" : ", ", signal_token( s.id ), signal_name( s.id ), prns
        );
    }
    const std::string rate_desc = decim > 1
                                      ? fmt::format( "{} Hz (decimated /{} from {} Hz)", sample_rate_hz, decim, native_rate )
                                      : fmt::format( "{} Hz", sample_rate_hz );
    logging::log(
        logging::Level::Info,
        fmt::format(
            "Streaming {} at {} - signals [{}] - launching {} channels", source_desc, rate_desc, sig_list, channels_.size()
        )
    );

    const unsigned num_workers = std::max( 1u, std::thread::hardware_concurrency() );
    pool_                      = std::make_unique<Thread_pool>( num_workers );

    channel_ptrs_.reserve( channels_.size() );
    for( auto& ch : channels_ )
    {
        channel_ptrs_.push_back( ch.get() );
    }

    scheduler_ = std::make_unique<Scheduler>( channel_ptrs_, *pool_, *sample_buffer_ );

    logging::log( logging::Level::Info, fmt::format( "Thread pool: {} workers", num_workers ) );
}

void Receiver::teardown()
{
    // Teardown order matters: stop producing, stop dispatching (no new tasks),
    // drain the pool (in-flight tasks still reference channels), then channels.
    if( device_ )
    {
        device_->stop_streaming();
    }
    decimator_.reset(); // safe now: the streaming thread (its only user) has stopped
    scheduler_.reset();
    pool_.reset();
    channels_.clear();
    channel_ptrs_.clear();
}

void Receiver::run()
{
    // running_ set BEFORE setup() so a stop() during setup is not lost (setup can't be interrupted
    // mid-way, but we bail before streaming if a stop arrived).
    running_ = true;
    setup();
    if( !running_ )
    {
        teardown();
        return;
    }

    // The buffer + channels run at the PROCESSING rate (after any decimation); the device still streams
    // at the native rate. So sample-index -> time uses the processing rate, but the device's consumed
    // count (native samples) divides by the native rate to give wall seconds.
    const double native_rate_hz   = static_cast<double>( config_.sample_rate_hz );
    const double sample_rate_hz   = native_rate_hz / std::max( 1u, config_.decimation );

    constexpr auto TICK_INTERVAL = std::chrono::milliseconds( 100 );
    constexpr int  LOG_EVERY_N   = 10;

    // Capture wall_start BEFORE start_streaming() so it is always <=
    // sample_buffer's internal stream_start_ (set on first push).
    // This guarantees: elapsed_samples <= wall_s * sample_rate_hz at every
    // point, so channel_s (derived from elapsed_samples) never exceeds wall_s.
    const auto wall_start = std::chrono::steady_clock::now();

    // Decimator layer (if enabled) sits here, between the source and the buffer: device -> decimator -> buffer.
    device_->start_streaming( [this]( const Complex_buf& samples ) {
        if( decimator_ )
        {
            sample_buffer_->push( decimator_->process( samples ) );
        }
        else
        {
            sample_buffer_->push( samples );
        }
    } );

    std::optional<Position_solution> current_position = std::nullopt;
    int                              tick             = 0;
    while( running_ )
    {
        std::this_thread::sleep_for( TICK_INTERVAL );

        // Feed the previous fix in so iono/tropo corrections can use it (off until one exists).
        const Ecef user_ecef =
            current_position ? Ecef { current_position->ecef_x_m, current_position->ecef_y_m, current_position->ecef_z_m }
                             : Ecef {};
        obs_engine_.generate(
            channel_ptrs_, scheduler_->min_next_sample(), sample_rate_hz, current_position ? &user_ecef : nullptr
        );
        current_position = pos_solver_.compute_solution( obs_engine_.get_measurements(), obs_engine_.reception_time_s() );

        // Loop timing/progress (computed every tick so the GUI status bar stays responsive).
        // exec: wall-clock since streaming started. stream: seconds of IQ pushed into the ring
        // buffer (can exceed exec - disk I/O outruns real-time). channel: the slowest channel's
        // processed position (what is actually correlated).
        const Receiver_status status {
            std::chrono::duration<double>( std::chrono::steady_clock::now() - wall_start ).count(),
            static_cast<double>( device_->samples_consumed() ) / native_rate_hz, // source (native) seconds streamed
            static_cast<double>( scheduler_->min_next_sample() ) / sample_rate_hz, // channel position (processing rate)
        };

        // Publish a coherent copy of the GUI-visible state (PVT + per-channel snapshots + status) under
        // the lock. Snapshotting here, on the run thread that owns the channels, keeps GUI readers off
        // the live channels entirely.
        {
            std::vector<Channel_snapshot> snaps;
            snaps.reserve( channel_ptrs_.size() );
            std::map<int, Gps_almanac> alm_updates; // gathered outside state_mutex_ (uses snapshot_mutex_)
            for( Channel* ch : channel_ptrs_ )
            {
                snaps.push_back( ch->snapshot() );
                for( const auto& [prn, a] : ch->almanac_snapshot() )
                {
                    alm_updates[prn] = a; // any channel's copy of an SV's almanac is the same broadcast data
                }
            }
            std::lock_guard<std::mutex> lock( state_mutex_ );
            latest_position_  = current_position;
            latest_snapshots_ = std::move( snaps );
            latest_status_    = status;
            // Receiver-wide almanac: accumulate (never drop), so an SV's coarse orbit persists for
            // acquisition aiding even after the channel that decoded it stops tracking.
            for( const auto& [prn, a] : alm_updates )
            {
                almanac_[prn] = a;
            }
        }

        publish_histories();

        // No clock steering here: the EKF carries the receiver clock bias + drift as states
        // and tracks them itself, so the master time anchor stays at its coarse value and we
        // do not feed the solved bias back (that would fight the filter's clock model).

        // Feed the solved clock DRIFT back into acquisition aiding as a rigorous common-mode
        // recenter. The drift is a range-rate (m/s); the equivalent fractional carrier offset
        // (df/f) the aiding wants is cd/c. Refreshed each fix so the recenter tracks the
        // (slowly varying) LO drift; takes precedence over the pre-PVT satellite-mean bridge.
        if( current_position && current_position->valid )
        {
            aiding_.set_clock_fraction( current_position->clock_drift_m_s / 299792458.0 );
            update_acquisition_predictions( *current_position );
        }

        if( ++tick % LOG_EVERY_N == 0 )
        {
            if( current_position && current_position->valid )
            {
                // Reference constellation always; an ISB only for non-reference constellations present.
                std::string ref_isb = fmt::format( "  ref={}", constellation_name( current_position->reference ) );
                for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
                {
                    if( current_position->isb_present[c] )
                    {
                        ref_isb += fmt::format(
                            "  isb[{}]={:.1f} m", constellation_name( static_cast<Constellation>( c ) ), current_position->isb_m[c]
                        );
                    }
                }
                logging::log(
                    logging::Level::Info,
                    fmt::format(
                        "Position - ECEF x={:.1f} m  y={:.1f} m  z={:.1f} m  clock bias={:.6f} m{}",
                        current_position->ecef_x_m,
                        current_position->ecef_y_m,
                        current_position->ecef_z_m,
                        current_position->clock_bias_m,
                        ref_isb
                    )
                );
                logging::log(
                    logging::Level::Info,
                    fmt::format(
                        "Velocity - ECEF x={:.1f} m/s  y={:.1f} m/s  z={:.1f} m/s  clock drift={:.6f} m/s",
                        current_position->ecef_x_m_s,
                        current_position->ecef_y_m_s,
                        current_position->ecef_z_m_s,
                        current_position->clock_drift_m_s
                    )
                );
            }

            logging::log(
                logging::Level::Info,
                fmt::format(
                    "exec={:.1f}s  streamed={:.1f}s  channels={:.1f}s", status.exec_s, status.stream_s, status.channel_s
                )
            );

            if( status.channel_s < status.exec_s - 0.2 )
            {
                logging::log(
                    logging::Level::Warning,
                    fmt::format( "Channels are {:.1f}s behind real time", status.exec_s - status.channel_s )
                );
            }
        }

        // Exit once the file is fully streamed AND every channel has consumed as far as
        // its own block size allows (acq_block differs per signal - GPS 1 ms vs Galileo
        // 4 ms - so each channel reports drained against its own limit).
        if( !device_->is_streaming() )
        {
            const Sample_index eof      = sample_buffer_->write_index();
            bool               all_done = true;
            for( const auto& ch : channels_ )
            {
                if( !ch->is_drained( eof ) )
                {
                    all_done = false;
                    break;
                }
            }
            if( all_done )
            {
                break;
            }
        }
    }

    // Distinguish a natural END OF FILE (loop broke while still running) from a user STOP (stop() cleared
    // running_). On EOF we FREEZE: keep the last published snapshots and capture every channel's graph
    // history before teardown destroys the channels, so panels/graphs opened afterwards show the last
    // state. On Stop we RESET the published state.
    const bool eof = running_.load( std::memory_order_relaxed );
    if( eof )
    {
        freeze_histories();
    }
    running_ = false;
    teardown();
    if( !eof )
    {
        clear_published_state();
    }
    logging::log( logging::Level::Info, eof ? "Reached end of file." : "Stopped." );
}

void Receiver::freeze_histories()
{
    // Snapshot every channel that has graph data (was tracking), keyed by SV, into the published map -
    // overriding the subscription filter so any SV's graph opened after EOF finds its last history. Runs
    // on the run thread before teardown (channels still alive).
    std::map<int, Tracking_history::Snapshot> frozen;
    for( Channel* ch : channel_ptrs_ )
    {
        Tracking_history::Snapshot snap = ch->history_snapshot();
        if( snap.iq.has_data || !snap.doppler.fast.hz.empty() )
        {
            frozen[history_key(
                ch->signal().params().constellation, static_cast<int>( ch->satellite_id() ), ch->signal().params().code
            )] = std::move( snap );
        }
    }
    std::lock_guard<std::mutex> lock( state_mutex_ );
    published_histories_ = std::move( frozen );
}

void Receiver::clear_published_state()
{
    std::lock_guard<std::mutex> lock( state_mutex_ );
    latest_position_ = std::nullopt;
    latest_snapshots_.clear();
    latest_status_ = {};
    published_histories_.clear();
    almanac_.clear(); // fresh run starts with an empty almanac (warm-start persistence is a future nicety)
    // history_subscriptions_ is left intact: graph windows that are still open want data again on a restart.
}

void Receiver::stop()
{
    running_ = false;
}

std::optional<Position_solution> Receiver::latest_position() const
{
    std::lock_guard<std::mutex> lock( state_mutex_ );
    return latest_position_;
}

std::vector<Channel_snapshot> Receiver::channel_snapshots() const
{
    std::lock_guard<std::mutex> lock( state_mutex_ );
    return latest_snapshots_;
}

Receiver_status Receiver::status() const
{
    std::lock_guard<std::mutex> lock( state_mutex_ );
    return latest_status_;
}

void Receiver::subscribe_history( Constellation constellation, int prn, Code code, bool on ) const
{
    const int                   key = history_key( constellation, prn, code );
    std::lock_guard<std::mutex> lock( state_mutex_ );
    if( on )
    {
        ++history_subscriptions_[key];
    }
    else if( auto it = history_subscriptions_.find( key ); it != history_subscriptions_.end() && --it->second <= 0 )
    {
        // last watcher for this SV closed -> drop the subscription + its published copy
        history_subscriptions_.erase( it );
        published_histories_.erase( key );
    }
}

std::optional<Tracking_history::Snapshot>
Receiver::published_history( Constellation constellation, int prn, Code code ) const
{
    std::lock_guard<std::mutex> lock( state_mutex_ );
    if( auto it = published_histories_.find( history_key( constellation, prn, code ) ); it != published_histories_.end() )
    {
        return it->second;
    }
    return std::nullopt;
}

void Receiver::publish_histories()
{
    // Snapshot the subscription set, copy each subscribed channel's history WITHOUT holding state_mutex_
    // (only the per-channel history mutexes), then swap the published map in under the lock. Runs on the
    // run thread, which owns the channels - the GUI only ever reads the published copies.
    std::set<int> subs;
    {
        std::lock_guard<std::mutex> lock( state_mutex_ );
        for( const auto& [key, count] : history_subscriptions_ )
        {
            subs.insert( key );
        }
    }
    if( subs.empty() )
    {
        return;
    }

    std::map<int, Tracking_history::Snapshot> fresh;
    for( Channel* ch : channel_ptrs_ )
    {
        const int key = history_key(
            ch->signal().params().constellation, static_cast<int>( ch->satellite_id() ), ch->signal().params().code
        );
        if( subs.count( key ) )
        {
            fresh[key] = ch->history_snapshot();
        }
    }

    std::lock_guard<std::mutex> lock( state_mutex_ );
    published_histories_ = std::move( fresh );
}

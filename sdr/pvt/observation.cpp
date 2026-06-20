#include "observation.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <utility>
#include "atmosphere.h"
#include "constants.h"
#include "geodesy.h"
#include "logging.h"

void Observation_engine::generate(
    const std::vector<Channel*>& channels, Sample_index rx_sample, double sample_rate_hz, const Ecef* user_ecef
)
{
    measurements_.clear();

    // Take ONE coherent snapshot per channel up front (the owning workers publish these); every read
    // below uses the snapshots, never the live channel state the workers concurrently mutate. This is
    // the fix for the async data race - the observables are now formed from race-free, self-consistent
    // per-channel state.
    std::vector<Channel_snapshot> snaps;
    snaps.reserve( channels.size() );
    for( const Channel* ch : channels )
    {
        snaps.push_back( ch->snapshot() );
    }


    // Cache the broadcast Klobuchar iono once any channel has decoded SF4 page 18. It is
    // GPS-system-wide, so the first SV's copy serves the whole fix and persists thereafter.
    if( !iono_.valid )
    {
        for( const Channel_snapshot& s : snaps )
        {
            if( s.iono.valid )
            {
                iono_ = s.iono;
                break;
            }
        }
    }

    // Receiver geodetic position from the previous fix - drives the look-angle geometry.
    // Stays off until the first fix gives us a user_ecef.
    const Geodetic user_geo = ( user_ecef != nullptr ) ? ecef_to_geodetic( *user_ecef ) : Geodetic { 0.0, 0.0, 0.0 };

    double t_rx_common_s = static_cast<double>( rx_sample ) / sample_rate_hz;

    if( !has_clock_offset_ )
    {
        // Coarse time base: anchor on the FIRST satellite that has an observable, ANY constellation.
        // This only needs to put t_rx_gps_tow in the right ballpark (~70 ms nominal transit); the
        // solver's clock-bias + per-constellation inter-system-bias states absorb the exact offsets.
        // So the anchor does not privilege a constellation or depend on acquisition order - a
        // different anchor constellation just shifts the (free) clock-bias state.
        for( const Channel_snapshot& s : snaps )
        {
            if( s.has_observable && s.has_lock )
            {
                // Raw pseudorange is roughly 0.07 s of signal travel time from space.
                master_clock_offset_s_ = s.transmit_time_at( rx_sample, sample_rate_hz ) + 0.070 - t_rx_common_s;
                has_clock_offset_      = true;
                break;
            }
        }
    }

    // 3. Skip generation if we haven't found our first coarse time anchor yet
    if( !has_clock_offset_ )
    {
        return;
    }

    double t_rx_gps_tow_s  = t_rx_common_s + master_clock_offset_s_;
    last_reception_time_s_ = t_rx_gps_tow_s;

    // TRUTH-RESIDUAL HARNESS (env TRUTH_LLH="lat_deg,lon_deg,alt_m"): validate per-SV observables against a
    // KNOWN static position (no EKF, clock-independent). Per snapshot (PRE-dedup, so per-code) it logs
    // resid = pseudorange - true_geometric_range = c*rx_clock (common) + per-SV ranging error; the offline
    // analyzer (tools/truth_residuals.py) removes the per-epoch common term (median) -> per-SV bias+jitter.
    // This is the metric for any observable work (sub-sample DLL, code biases) - NOT the +/-150 m EKF position.
    bool have_truth = false;
    Ecef truth_rx {};
    if( const char* truth_env = std::getenv( "TRUTH_LLH" ) )
    {
        double lat_deg = 0.0, lon_deg = 0.0, h_m = 0.0;
        if( std::sscanf( truth_env, "%lf,%lf,%lf", &lat_deg, &lon_deg, &h_m ) == 3 )
        {
            const double lat = lat_deg * M_PI / 180.0, lon = lon_deg * M_PI / 180.0;
            const double f = 1.0 / 298.257223563, e2 = f * ( 2.0 - f ), a = 6378137.0; // WGS-84
            const double N = a / std::sqrt( 1.0 - e2 * std::sin( lat ) * std::sin( lat ) );
            truth_rx   = { ( N + h_m ) * std::cos( lat ) * std::cos( lon ),
                         ( N + h_m ) * std::cos( lat ) * std::sin( lon ),
                         ( N * ( 1.0 - e2 ) + h_m ) * std::sin( lat ) };
            have_truth = true;
        }
    }
    if( have_truth )
    {
        for( const Channel_snapshot& s : snaps )
        {
            if( !s.has_observable || !s.has_lock )
            {
                continue;
            }
            const Constellation con     = s.constellation;
            const double        t_tx    = s.transmit_time_at( rx_sample, sample_rate_hz );
            const double        transit = t_rx_gps_tow_s - t_tx;
            const Ecef          sv      = orbit::satellite_ecef_pos( s.eph, t_tx, con );
            const double        cs = std::cos( constants::EARTH_ROTATION_RATE_RAD_S * transit ), sn = std::sin( constants::EARTH_ROTATION_RATE_RAD_S * transit );
            const Ecef          svr { sv.x * cs + sv.y * sn, -sv.x * sn + sv.y * cs, sv.z }; // Sagnac to rx epoch
            const double range = std::sqrt( ( svr.x - truth_rx.x ) * ( svr.x - truth_rx.x )
                                            + ( svr.y - truth_rx.y ) * ( svr.y - truth_rx.y )
                                            + ( svr.z - truth_rx.z ) * ( svr.z - truth_rx.z ) );
            const double sv_clk = orbit::satellite_clock_offset( s.eph, t_tx, con );
            const double resid  = transit * constants::SPEED_OF_LIGHT_M_S + sv_clk * constants::SPEED_OF_LIGHT_M_S - range;
            logging::log(
                logging::Level::Info,
                fmt::format(
                    "TRUTHDIAG rx={} con={} prn={:2d} code={} resid={:.3f} range={:.1f}",
                    static_cast<uint64_t>( rx_sample ), static_cast<int>( con ), static_cast<int>( s.satellite_id ),
                    static_cast<int>( s.code ), resid, range
                )
            );
        }
    }

    // De-duplicate: a satellite tracked on two components (e.g. GPS L1CA + L1Cd) must contribute ONE row
    // to PVT, not two. Group the eligible snapshots by (constellation, satellite_id) and keep one per SV.
    // A snapshot is eligible only if it can produce an observable (valid ephemeris + anchored TOW) AND is
    // carrier-locked - gating on lock (the standard gnss-sdr policy) keeps a cycle-slipping / fading
    // channel's corrupted range + range-rate out of the EKF; the lock detector's hysteresis means a brief
    // dip won't drop a healthy SV. (Cross-correlation false tracks keep a real carrier, so they still pass
    // here - the frame-sync timeout is what evicts those.)
    //
    // Selection is by a FIXED component priority, NOT C/N0 - deliberately, even though C/N0 is now comparable
    // across components (the SNV estimator replaced M2M4, which used to saturate ~38 dB-Hz on the long-epoch
    // components and made cross-component C/N0 meaningless). Two reasons to keep fixed priority for SAME-BAND
    // (L1CA vs L1Cd): (1) the L1C PILOT carries only a fraction of the signal power, so its C/N0 is
    // SYSTEMATICALLY ~2-3 dB below L1CA's (measured) - a C/N0 race would just always pick L1CA anyway; (2)
    // fixed priority is DETERMINISTIC, whereas a C/N0 race could flip component mid-run on a momentary
    // fluctuation, and each flip jumps the pseudorange by ~1 sample and jostles the EKF. So pick by code_rank
    // and let the has_observable && has_lock gate handle failover (preferred component when healthy, else the
    // other carries the SV). The comparable C/N0 instead pays off CROSS-BAND (future L1+L5), where the
    // components have genuinely different quality worth racing on.
    auto code_rank = []( Code c )
    {
        // higher = preferred; only same-constellation codes ever compete here
        switch( c )
        {
        case Code::CA:
            return 1; // GPS L1 C/A - proven + carries the broadcast iono (primary)
        case Code::Cd:
            return 2; // GPS L1Cd - BOC pilot, used when L1CA is unavailable for the SV
        default:
            return 0; // Galileo E1-B / BeiDou B1I - no sibling code competes yet
        }
    };

    std::map<std::pair<Constellation, Satellite_id>, const Channel_snapshot*> best;
    for( const Channel_snapshot& s : snaps )
    {
        if( !s.has_observable || !s.has_lock )
        {
            continue;
        }
        const auto key = std::make_pair( s.constellation, s.satellite_id );
        const auto it  = best.find( key );
        if( it == best.end() || code_rank( s.code ) > code_rank( it->second->code ) )
        {
            best[key] = &s;
        }
    }

    for( const auto& [key, winner] : best )
    {
        const Channel_snapshot& s = *winner;

        Satellite_measurement m;
        m.satellite_id  = s.satellite_id;
        m.constellation = s.constellation;

        const double t_tx      = s.transmit_time_at( rx_sample, sample_rate_hz );
        m.transmit_time_s      = t_tx;
        m.pseudorange_m        = ( t_rx_gps_tow_s - t_tx ) * constants::SPEED_OF_LIGHT_M_S;
        m.pseudorange_rate_m_s = s.range_rate_at( rx_sample, sample_rate_hz );

        // Satellite ECEF state + clock from the broadcast ephemeris, dispatched per constellation
        // (GPS GPST / Galileo GST; same Keplerian model, different mu).
        const Constellation con = s.constellation;

        Ecef sat_pos                = orbit::satellite_ecef_pos( s.eph, t_tx, con );
        m.satellite_pos_x           = sat_pos.x;
        m.satellite_pos_y           = sat_pos.y;
        m.satellite_pos_z           = sat_pos.z;
        m.satellite_clock_offset_s  = orbit::satellite_clock_offset( s.eph, t_tx, con );
        Ecef sat_vel                = orbit::satellite_ecef_vel( s.eph, t_tx, con );
        m.satellite_vel_x           = sat_vel.x;
        m.satellite_vel_y           = sat_vel.y;
        m.satellite_vel_z           = sat_vel.z;
        m.satellite_clock_drift_s_s = orbit::satellite_clock_drift( s.eph, t_tx, con );

        // SV clock correction. IS-GPS-200: GPS_time = SV_time - dt_sv, so the true GPS transmit time
        // is t_tx - dt_sv, giving pr_corrected = pr_raw + dt_sv*c (and likewise the rate). (Sign was
        // inverted before, which spread the per-SV residuals by ~2*dt_sv*c and threw the fix ~1000 km off.)
        m.pseudorange_m += m.satellite_clock_offset_s * constants::SPEED_OF_LIGHT_M_S;
        m.pseudorange_rate_m_s += m.satellite_clock_drift_s_s * constants::SPEED_OF_LIGHT_M_S;

        double transit_time_s = t_rx_gps_tow_s - t_tx;
        double earth_spin_rad = constants::EARTH_ROTATION_RATE_RAD_S * transit_time_s;

        double cos_spin = std::cos( earth_spin_rad );
        double sin_spin = std::sin( earth_spin_rad );

        // Rotate Position
        m.satellite_pos_x = sat_pos.x * cos_spin + sat_pos.y * sin_spin;
        m.satellite_pos_y = -sat_pos.x * sin_spin + sat_pos.y * cos_spin;
        m.satellite_pos_z = sat_pos.z;

        // Rotate Velocity
        m.satellite_vel_x = sat_vel.x * cos_spin + sat_vel.y * sin_spin;
        m.satellite_vel_y = -sat_vel.x * sin_spin + sat_vel.y * cos_spin;
        m.satellite_vel_z = sat_vel.z;

        // Look angles need only a receiver position (the previous fix), so compute them whenever we
        // have one - the solver weights each satellite by elevation. The atmospheric delays both
        // lengthen the measured range, so subtract them: tropo needs only the geometry; the Klobuchar
        // iono also needs the broadcast SF4-p18 coefficients (so it waits for those, tropo does not).
        if( user_ecef != nullptr )
        {
            const Ecef sat_pos_ecef { m.satellite_pos_x, m.satellite_pos_y, m.satellite_pos_z };
            double     elevation_rad = 0.0;
            double     azimuth_rad   = 0.0;
            look_angles( *user_ecef, sat_pos_ecef, elevation_rad, azimuth_rad );
            m.elevation_rad = elevation_rad;

            double correction_m = tropo_delay_m( elevation_rad, user_geo.alt_m );
            if( iono_.valid && iono_.model == Iono::Model::Klobuchar )
            {
                correction_m += klobuchar_iono_delay_m(
                    iono_.alpha, iono_.beta, user_geo.lat_rad, user_geo.lon_rad, elevation_rad, azimuth_rad, t_rx_gps_tow_s
                );
            }
            m.pseudorange_m -= correction_m;
        }

        // Hatch carrier-smoothing (when enabled): within one continuous lock arc, blend the noisy code
        // pseudorange with the precise carrier-phase delta (integrated Doppler). N ramps to HATCH_WINDOW (the
        // iono-divergence bound). Start/restart the arc on the first epoch, a re-acquisition (lock_session
        // change), or a dedup component flip (code change); lock loss is already gated out (only locked SVs
        // reach here). Disabled -> the raw code pseudorange passes through untouched.
        if( hatch_enabled_ )
        {
            const int    sv      = sv_key( s.constellation, static_cast<int>( s.satellite_id ) );
            const double phase_m = s.carrier_phase_range_m( rx_sample, sample_rate_hz );
            Hatch_state& h       = hatch_[sv];
            if( h.n == 0 || h.session != s.lock_session || h.code != s.code )
            {
                h.pr_smooth = m.pseudorange_m; // (re)start the arc on the raw code pr
                h.n         = 1;
            }
            else
            {
                h.n = std::min( h.n + 1, HATCH_WINDOW );
                const double predicted = h.pr_smooth + ( phase_m - h.phase_prev_m ); // carrier-propagated from k-1
                h.pr_smooth            = m.pseudorange_m / h.n + ( 1.0 - 1.0 / h.n ) * predicted;
            }
            h.phase_prev_m  = phase_m;
            h.session       = s.lock_session;
            h.code          = s.code;
            m.pseudorange_m = h.pr_smooth;
        }

        measurements_.push_back( m );
    }
}

const std::vector<Satellite_measurement>& Observation_engine::get_measurements() const
{
    return measurements_;
}

void Observation_engine::adjust_master_clock_offset( double correction_s )
{
    master_clock_offset_s_ -= correction_s;
}

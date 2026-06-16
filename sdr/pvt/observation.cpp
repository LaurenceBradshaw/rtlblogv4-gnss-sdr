#include "observation.h"
#include <map>
#include <utility>
#include "atmosphere.h"
#include "geodesy.h"
#include "logging.h"

void Observation_engine::generate(
    const std::vector<Channel*>& channels, Sample_index rx_sample, double sample_rate_hz, const Ecef* user_ecef
)
{
    const double c                    = 299792458.0; // speed of light in m/s
    const double EARTH_ROTATION_SPEED = 7.292115e-5; // rad/s
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

    // De-duplicate: a satellite tracked on two components (e.g. GPS L1CA + L1Cd) must contribute ONE row
    // to PVT, not two. Group the eligible snapshots by (constellation, satellite_id) and keep one per SV.
    // A snapshot is eligible only if it can produce an observable (valid ephemeris + anchored TOW) AND is
    // carrier-locked - gating on lock (the standard gnss-sdr policy) keeps a cycle-slipping / fading
    // channel's corrupted range + range-rate out of the EKF; the lock detector's hysteresis means a brief
    // dip won't drop a healthy SV. (Cross-correlation false tracks keep a real carrier, so they still pass
    // here - the frame-sync timeout is what evicts those.)
    //
    // Selection is by a FIXED component priority, NOT C/N0. C/N0 is not comparable across components: the
    // M2M4 estimator saturates at high per-epoch SNR, and a component's per-epoch SNR scales with its
    // coherent integration, so GPS L1Cd (10 ms code) reads a flat ~38 dB-Hz regardless of true strength
    // while L1CA (1 ms) reads its true 41-46 (verified same-SV: CA-Cd gap 3-8 dB, Cd flat across SVs).
    // Comparing those magnitudes would unfairly favour L1CA. So pick deterministically by code_rank and let
    // the has_observable && has_lock gate handle failover: the preferred component is used whenever it is
    // healthy, else the other carries the SV. The M2 fix made the per-component transmit times agree, so
    // the choice - and failover - is measurement-consistent. (A high-SNR-robust C/N0 estimator could later
    // restore quality-weighted selection; see the multi-code-fusion backlog.)
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
        m.pseudorange_m        = ( t_rx_gps_tow_s - t_tx ) * c;
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
        m.pseudorange_m += m.satellite_clock_offset_s * c;
        m.pseudorange_rate_m_s += m.satellite_clock_drift_s_s * c;

        double transit_time_s = t_rx_gps_tow_s - t_tx;
        double earth_spin_rad = EARTH_ROTATION_SPEED * transit_time_s;

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

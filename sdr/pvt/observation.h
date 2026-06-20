#pragma once
#include <cstdint>
#include <map>
#include <vector>
#include "channel.h"
#include "constellations.h"
#include "orbit.h"

struct Satellite_measurement
{
    Satellite_id  satellite_id;
    Constellation constellation;
    double        satellite_pos_x;
    double        satellite_pos_y;
    double        satellite_pos_z;
    double        satellite_clock_offset_s;
    double        satellite_vel_x;
    double        satellite_vel_y;
    double        satellite_vel_z;
    double        satellite_clock_drift_s_s;
    double        transmit_time_s;
    double        pseudorange_m;
    // double       carrier_phase_cycles; // Later
    // double doppler_hz; // Later
    double pseudorange_rate_m_s; // Calculated from Doppler, but converted to m/s for position solver convenience
    double elevation_rad = 0.0;  // satellite elevation from the last user fix (0 until one exists); for weighting
};

// Iono (Klobuchar) + tropo corrections are applied here, per satellite, in generate(): the
// models live in atmosphere.h and the geometry in geodesy.h. They need a receiver position,
// so the previous fix (user_ecef) is fed in - corrections stay off until the first fix exists.

class Observation_engine
{
public:
    // common rx_sample (e.g. scheduler->min_next_sample()); raw PR = (t_rx_common - t_tx)*c.
    // user_ecef: the previous position fix, used for the iono/tropo look-angles + lat/lon;
    //            nullptr before the first fix (atmospheric corrections are then skipped).
    void generate(
        const std::vector<Channel*>& channels,
        Sample_index                 rx_sample,
        double                       sample_rate_hz,
        const Ecef*                  user_ecef
    );

    const std::vector<Satellite_measurement>& get_measurements() const;

    // Common GPS reception time (s) of the most recent generate() - the EKF uses it for dt.
    double reception_time_s() const
    {
        return last_reception_time_s_;
    }

    void adjust_master_clock_offset( double correction_s );

private:
    std::vector<Satellite_measurement> measurements_;

    bool   has_clock_offset_      = false;
    double master_clock_offset_s_ = 0.0; // TODO: Currently doesn't account for week rollover.
    // TODO: Is also currently GPS specific.

    // Broadcast Klobuchar iono, cached from the first channel that decodes SF4 page 18 (it is
    // GPS-system-wide, so any one SV's copy serves every satellite in the fix).
    Iono iono_;

    double last_reception_time_s_ = 0.0; // t_rx_gps_tow_s of the last generate()

    // Hatch carrier-smoothing of the code pseudorange (per selected SV). Within one continuous carrier-phase
    // arc, blend the noisy code pr with the precise carrier-phase delta. HATCH_WINDOW bounds the arc so the
    // single-frequency code-carrier IONO divergence (~2*delta_iono*N) stays small (~100 s on L1); harmless at
    // any N on an iono-free capture. N=1 disables smoothing. Re-tune / go divergence-free with dual frequency.
    static constexpr int HATCH_WINDOW = 100; // max samples averaged (PVT epochs; ~100 s at a 1 Hz fix rate)
    struct Hatch_state
    {
        double   pr_smooth    = 0.0;       // last smoothed pseudorange (m)
        double   phase_prev_m = 0.0;       // carrier-phase range (m) at the previous epoch
        int      n            = 0;         // window count so far (ramps to HATCH_WINDOW; 0 = no arc yet)
        uint32_t session      = 0;         // lock_session of the arc this state belongs to (reset on change)
        Code     code         = Code::CA;  // selected component of the arc (reset on a dedup failover flip)
    };
    std::map<int, Hatch_state> hatch_; // keyed by sv_key(constellation, prn)
};

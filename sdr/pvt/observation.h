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
    // hatch_enabled: carrier-smooth the code pseudorange (HATCH_WINDOW). Off -> raw code pseudorange.
    explicit Observation_engine( bool hatch_enabled = true ) : hatch_enabled_( hatch_enabled ) {}

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
    bool                 hatch_enabled_ = true; // off -> raw code pseudorange (smoothing block skipped)
    struct Hatch_state
    {
        double   pr_smooth    = 0.0;       // last smoothed pseudorange (m)
        double   phase_prev_m = 0.0;       // carrier-phase range (m) at the previous epoch
        int      n            = 0;         // window count so far (ramps to HATCH_WINDOW; 0 = no arc yet)
        uint32_t session      = 0;         // lock_session of the arc this state belongs to (reset on change)
        Code     code         = Code::CA;  // selected component of the arc (reset on a dedup failover flip)
    };
    std::map<int, Hatch_state> hatch_; // keyed by sv_key(constellation, prn)

    // Cycle-slip / observable-discontinuity detection (code-minus-carrier step). CMC = code_pr -
    // carrier_phase_range is ~constant within one continuous carrier arc (only slow iono moves it; the integer
    // ambiguity is a per-arc constant that cancels in the step), so an abrupt CMC step between epochs flags a
    // discontinuity. Runs per-SV REGARDLESS of Hatch (reusable). Reseeds on lock_session change so a re-acquire
    // is not mistaken for a slip. Single-frequency limit: catches the IMPACTFUL events (multi-cycle / loss-of-
    // lock, exceeding code noise), not sub-metre ones - which are also low-impact.
    //
    // The threshold ADAPTS per SV to that signal's own ΔCMC noise. Measured RMS(ΔCMC) is ~1.7 m on GPS L1CA
    // (BPSK) but ~3.7 m on Galileo E1 (BOC + double-estimator, ~2x noisier per epoch), so one fixed threshold
    // cannot fit both - tight enough for GPS slips trips on BOC noise tails (~15 m peaks). So flag when
    // |ΔCMC| > SLIP_SIGMA_K * running-σ(ΔCMC), floored. K=6 sits well above the ~3σ noise peaks of both signals
    // (zero false triggers on the clean sim) while still catching the slips that matter, and tracks up
    // automatically in a high-multipath (real-world) environment. NOTE: on BOC this is really an OBSERVABLE-
    // discontinuity detector (the DE can step the code independently of the carrier) - conservative-safe for
    // Hatch, but a pure carrier-cycle-slip signal (a tracker PLL loss-of-lock indicator) is the separate next
    // layer, needed for PPP carrier ambiguities.
    static constexpr double SLIP_SIGMA_K = 6.0;  // flag a step beyond this many σ of the running ΔCMC noise
    static constexpr double SLIP_FLOOR_M = 8.0;  // ...but never below this (small slips are low-impact)
    static constexpr double SLIP_EMA     = 0.05; // EMA weight for the running ΔCMC variance (~20-epoch window)
    static constexpr int    SLIP_WARMUP  = 20;   // epochs to let σ settle before detecting (arc start = pull-in)
    struct Cmc_state
    {
        double   prev_cmc = 0.0;
        double   var_ema  = 0.0; // EMA of ΔCMC² -> running noise variance backing the adaptive threshold
        uint32_t session  = 0;   // arc id (lock_session); a change reseeds (a re-acquire is not a slip)
        int      n        = 0;   // ΔCMC samples seen this arc (0 = unprimed; also gates the warm-up)
    };
    std::map<int, Cmc_state> cmc_;                                               // keyed by sv_key
    bool cmc_slip( int sv, double code_pr_m, double phase_m, uint32_t session ); // layer-1 (code-domain) detector

    // Layer 2: carrier-domain cycle-slip / loss-of-lock indicator. The tracker counts carrier phase-lock breaks
    // per arc (cos(2φ) dropping out of lock - independent of the code / double-estimator, so it catches the
    // Costas full-cycle slips the CMC test relies on too); a rise in that count between epochs is a slip. This
    // is the pure carrier signal PPP's ambiguities key on; Hatch ORs it with the CMC test. Reseeds on arc change.
    struct Lockbreak_state
    {
        int      breaks  = 0;     // last-seen carrier_lock_breaks for this SV (this arc)
        uint32_t session = 0;     // arc id (lock_session)
        bool     primed  = false; // false until the first sample of an arc (which only seeds)
    };
    std::map<int, Lockbreak_state> lockbreak_;                           // keyed by sv_key
    bool carrier_lock_slip( int sv, int lock_breaks, uint32_t session ); // updates lockbreak_
};

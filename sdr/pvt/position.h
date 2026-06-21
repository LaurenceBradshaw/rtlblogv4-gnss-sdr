#pragma once
#include <Eigen/Dense>
#include <map>
#include <optional>
#include <tuple>
#include <vector>
#include "observation.h"

struct Position_solution
{
    bool valid;

    double ecef_x_m;
    double ecef_y_m;
    double ecef_z_m;
    double clock_bias_m;
    double ecef_x_m_s;
    double ecef_y_m_s;
    double ecef_z_m_s;
    double clock_drift_m_s;

    // The reference constellation (the receiver clock bias is relative to ITS system time); always set.
    Constellation reference = Constellation::Unknown;
    // Per-constellation inter-system bias (m). isb_present[c] is true ONLY for a non-reference
    // constellation that is actually contributing - so a single-constellation fix has none to show.
    double isb_m[NUM_CONSTELLATIONS]      = {};
    bool   isb_present[NUM_CONSTELLATIONS] = {};
};

// Extended Kalman filter PVT estimator (multi-constellation).
//
// State (N = 8 + NUM_CONSTELLATIONS): ECEF position (m), ECEF velocity (m/s), receiver clock bias
// (m), clock drift (m/s), then one inter-system bias (m) PER constellation:
//     x = [px py pz  vx vy vz  cb cd  isb[Gps] isb[Galileo] isb[Beidou]]
// One constellation is the REFERENCE (highest-priority one present, see constellations.h): cb is
// the receiver clock relative to its system time, and its own isb slot is held at 0 (folded into
// cb - you can only observe inter-constellation DIFFERENCES, so N constellations need N-1 biases).
// Every other constellation's pseudorange predicts range + cb + isb[c]; that isb absorbs the
// system-time offset + receiver inter-system hardware delay (near-constant -> a slow random walk).
// The filter is SEEDED from a MIXED least-squares fix (position + clock + one bias per non-reference
// constellation present), needing 4 + (#non-reference) sats and passing a PDOP geometry gate - so it
// is not blocked on any one constellation or on acquisition order, only on having a well-conditioned
// fix. Biases solved in the seed start moderately known; a constellation that appears LATER bootstraps
// from a large initial covariance so the innovation gate accepts its first (possibly ms-level) update.
//
// Process model (predict): constant-velocity for position (p += v*dt) with a slowly varying
// velocity, and the standard coupled receiver-clock model (cb += cd*dt). The clock states are
// coupled because this SDR's single reference oscillator drives both the RF LO and the ADC
// sample clock, so its offset shows up consistently in the carrier Doppler AND the pseudorange
// (the clock bias is seen drifting at ~-1346 m/s on this capture, matching cd). cd thus absorbs
// the common-mode Doppler, leaving receiver velocity ~0 for the static CTTC capture.
//
// Measurement (update): pseudoranges (-> position + clock bias) and pseudorange rates
// (-> velocity + clock drift), each weighted by satellite elevation (low elevation = noisier,
// the larger atmospheric path + multipath). The very first fix is seeded by a plain
// least-squares solve, because the filter needs an initial state + covariance to start from;
// every fix after that is a predict + update cycle.
class Position_solver
{
public:
    // measurements: this epoch's per-satellite observables.
    // rx_time_s:    common GPS reception time of this set (s) - used for the predict-step dt.
    // Returns nullopt until the filter has been seeded (enough sats + good enough geometry).
    std::optional<Position_solution>
    compute_solution( const std::vector<Satellite_measurement>& measurements, double rx_time_s );

    // Enable tightly-coupled carrier-phase fusion (float per-arc ambiguities). OFF by default: with BROADCAST
    // products on SHORT arcs the ambiguities can't converge (little geometry change to separate bias from
    // position) and it degrades the code+TDCP fix - it only pays off with precise products + long arcs (PPP).
    // The machinery is in place; flip this on when those land.
    void set_carrier_phase_enabled( bool on )
    {
        carrier_phase_enabled_ = on;
    }

private:
    static constexpr int ISB_BASE = 8;                          // first inter-system-bias slot
    static constexpr int N        = ISB_BASE + NUM_CONSTELLATIONS; // state dimension (one isb per constellation)
    // Seed geometry gate: reject a seed whose position dilution of precision exceeds this, so the
    // filter is not locked onto a sparse, badly-conditioned fix (wait for more/better sats instead).
    static constexpr double MAX_SEED_PDOP = 6.0;
    // Seed pseudorange RAIM: drop an SV from the seed least-squares if its post-fit range residual
    // exceeds this (m). A healthy seed residual is metres-to-tens-of-metres (measurement noise +
    // unmodelled iono/tropo); a cross-correlation / false-lock track is orders of magnitude larger and
    // would otherwise drag the seed position hundreds of km off (PDOP doesn't catch it - geometry is
    // fine, the measurement is not). Generous so a marginal-but-real SV is never dropped.
    static constexpr double SEED_PR_RAIM_RESIDUAL_M = 10000.0;

    // --- tuning (process noise spectral densities + measurement noise) -------------------
    // Spatial: continuous white-noise-acceleration model. ACCEL_PSD is how hard the platform
    // can accelerate (m^2/s^3); small = trust the constant-velocity model (good for a static
    // or slow receiver), large = let position/velocity move freely.
    static constexpr double ACCEL_PSD = 1.0;
    // Clock random-walk PSDs (m^2/s for the bias, m^2/s^3 for the drift). Both small: the
    // bias and drift are near-constant on this capture.
    static constexpr double CLOCK_BIAS_PSD  = 1.0;
    static constexpr double CLOCK_DRIFT_PSD = 1.0;
    // Inter-system bias random-walk PSD (m^2/s). Tiny - the GPS/Galileo offset is near-constant.
    static constexpr double ISB_PSD = 0.1;
    // Base measurement standard deviations (at zenith); scaled by 1/sin(elevation).
    static constexpr double PR_STD_M    = 5.0; // pseudorange (m)
    static constexpr double PRR_STD_M_S = 0.5; // pseudorange rate from instantaneous Doppler (m/s)
    // Pseudorange rate from TDCP (time-differenced carrier phase) - much lower noise than Doppler, so the
    // solver trusts it far more. Conservative fixed value; the true σ scales as ~σ_phase/Δt (a refinement
    // parked with the acceleration-model backlog item). Set per measurement via Satellite_measurement::prr_from_tdcp.
    static constexpr double PRR_TDCP_STD_M_S = 0.05;
    // Elevation floor for the 1/sin(el) weighting, so low/unknown-elevation sats stay finite.
    static constexpr double MIN_ELEVATION_RAD = 5.0 * 3.14159265358979323846 / 180.0;
    // Innovation gate: reject a measurement if normalised innovation^2 exceeds this (~6 sigma).
    // Keeps end-of-file garbage / gross outliers from dragging the filter off.
    static constexpr double INNOV_GATE2 = 36.0;
    // Velocity RAIM: a pseudorange-rate whose robust-LS residual exceeds this (m/s) is a gross
    // outlier (a marginally-locked SV with a biased carrier Doppler) and is excluded from the
    // velocity/clock-drift solve. Good SVs sit well under 1 m/s; a biased one is tens of m/s.
    static constexpr double VELOCITY_RAIM_RESIDUAL_M_S = 5.0;

    // --- carrier-phase float ambiguities (tightly-coupled code+carrier) ------------------
    // One bias state (m) per SV-arc, appended to x_ after the base N states, keyed by (constellation, prn,
    // code). The carrier-phase observable is range + clock + isb + bias; the bias absorbs the unknown
    // integer-cycle ambiguity (+ the arc-start range/clock, since carrier_phase_range starts at 0 each arc).
    // A cycle slip or a lock_session change (re-acquire) re-initialises that bias. Slots are not freed
    // (bounded by the channel count); pruning long-dead arcs is a future refinement.
    static constexpr double PHASE_STD_M     = 0.01;   // carrier-phase measurement noise (m, at zenith)
    static constexpr double AMB_INIT_VAR_M2 = 1.0e6;  // initial bias variance (covers the clock-bias unknown)
    static constexpr double AMB_PSD         = 1.0e-6; // bias random-walk PSD (m^2/s) - essentially constant

    bool          carrier_phase_enabled_ = false; // tightly-coupled carrier-phase fusion (see setter)
    bool          initialised_ = false;
    double        last_time_s_ = 0.0;
    Constellation reference_   = Constellation::Gps; // set at init; its isb slot stays 0

    // Which non-reference constellations have contributed since the seed (latched, so a momentary
    // dropout doesn't blink the ISB display off). Reset on (re)initialise.
    bool isb_present_[NUM_CONSTELLATIONS] = {};

    // State + covariance are DYNAMIC-size: the base N states today, growing by one float carrier-phase
    // ambiguity per active SV-arc once carrier-phase fusion lands. The base layout (PX/VX/CB/CD/ISB) keeps
    // its fixed indices; ambiguities occupy [N..]. Sized in initialise().
    Eigen::VectorXd x_; // state estimate
    Eigen::MatrixXd P_; // state covariance (sized to match x_)

    // Active carrier-phase ambiguities: key (constellation, prn, code) -> its bias state index in x_ (>= N)
    // and the lock_session of its current arc (a change re-initialises the bias). One entry per tracked SV-arc.
    struct Amb_key
    {
        Constellation con;
        Satellite_id  prn;
        Code          code;
        bool          operator<( const Amb_key& o ) const
        {
            return std::tie( con, prn, code ) < std::tie( o.con, o.prn, o.code );
        }
    };
    struct Amb_state
    {
        int      index;
        uint32_t session;
    };
    std::map<Amb_key, Amb_state> amb_;

    // Seed the filter with a least-squares fix on the reference constellation. Returns false if
    // no constellation yet has MIN_REFERENCE_SATS satellites.
    bool initialise( const std::vector<Satellite_measurement>& measurements );
    // Propagate state + covariance forward by dt under the process model.
    void predict( double dt );
    // Fuse this epoch's pseudoranges + pseudorange rates into the state.
    void update( const std::vector<Satellite_measurement>& measurements );
    // Mark each non-reference constellation present in `measurements` as active (latched).
    void latch_present( const std::vector<Satellite_measurement>& measurements );

    // Carrier-phase ambiguity management (Phase 2). allocate_ambiguity appends a bias state (x_/P_ grow by
    // one); reset_ambiguity re-initialises an existing one (cycle slip / re-acquire); ensure_ambiguity does
    // whichever is needed for this measurement's arc and returns the bias state index.
    int  allocate_ambiguity( double init_bias );
    void reset_ambiguity( int index, double init_bias );
    int  ensure_ambiguity( const Satellite_measurement& sm );

    Position_solution as_solution() const;
};

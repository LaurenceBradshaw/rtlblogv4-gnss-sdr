#pragma once
#include "tracking_core.h"

// Tracker for pilot+data BOC signals (Galileo E1, GPS L1C, BeiDou B1C). Tracks the data-
// free PILOT component: Costas + FLL pull-in while syncing the secondary (e.g. CS25) code,
// then wipes the secondary off and runs a pure 4-quadrant PLL (no FLL) on the now-coherent
// pilot for a clean lock. (VEML 5-tap DLL and extended coherent integration are Phases D/E;
// the carrier/code NCO + correlator come from Tracking_core.)
class Pilot_tracker : public Tracking_core
{
public:
    using Tracking_core::Tracking_core; // inherit the (code, sample_rate, sig, secondary) ctor

    void run_loops( bool bit_sync, bool sw_loop, Satellite_id prn ) override;

protected:
    void   configure_taps( double ci ) override; // DE 5-tap: P, code E/L, subcarrier E/L
    double code_error() const override;          // envelope code DLL from the code E/L gates (subcarrier-indep.)

private:
    // Double-estimator SUBCARRIER loop (SLL): drives subcarrier_offset_ (Tracking_core) from the subcarrier
    // E/L gates so the subcarrier tracks the received signal independently of the code phase. It may lock onto
    // any subcarrier lobe; the integer ambiguity is resolved in code_phase_offset_s() by rounding against the
    // code phase (NOT here - so do NOT clamp/wrap subcarrier_offset_). Run only once carrier+secondary synced.
    void                    subcarrier_update();
    static constexpr double SUBC_GAIN = 0.2; // 1st-order SLL gain (code elements per unit discriminator)
};

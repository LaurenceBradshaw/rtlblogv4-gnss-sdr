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
    void   configure_taps( double ci ) override; // 5-tap VE/E/P/L/VL (BOC false-lock handling)
    double code_error() const override;          // VEMLP: (sqrt(|VE|^2+|E|^2)-sqrt(|L|^2+|VL|^2))/sum
};

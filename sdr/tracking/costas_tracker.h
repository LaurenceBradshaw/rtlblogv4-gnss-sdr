#pragma once
#include "tracking_core.h"

// Tracker for data-bearing BPSK signals (GPS L1 C/A, L2C, BeiDou B1I, Galileo E1-B data):
// 3-tap E/P/L correlator + Costas PLL + (cross/dot or atan) FLL. The carrier/code NCO and
// the correlator live in Tracking_core; this class only adds the per-epoch loop strategy.
class Costas_tracker : public Tracking_core
{
public:
    using Tracking_core::Tracking_core; // inherit the (code, sample_rate, sig, secondary) ctor

    void run_loops( bool bit_sync, bool sw_loop, Satellite_id prn ) override;

protected:
    void   correlate( const Sample_block& block, int n ) override; // BPSK replica (reads the code array)
    void   configure_taps( double ci ) override;                   // 3-tap P/E/L at +/- corr_spacing_ samples
    double code_error() const override;                            // normalised (|E|-|L|)/(|E|+|L|)
};

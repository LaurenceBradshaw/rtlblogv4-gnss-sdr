#pragma once
#include "tracking_core.h"

// Tracker for pilot+data BOC signals (Galileo E1, GPS L1C, BeiDou B1C). Tracks the data-free PILOT component:
// Costas + FLL pull-in while syncing the secondary (e.g. CS25) code, then wipes the secondary off and runs a
// pure 4-quadrant PLL on the now-coherent pilot. Ranging uses the Hodgart/Blunt DOUBLE ESTIMATOR: a 7-tap
// correlator (Boc_de_replica) tracks the code (envelope) and subcarrier as two independent delays; the precise
// subcarrier estimate resolves the BOC ambiguity against the code phase (see code_phase_offset_s). The
// carrier/code NCO + the shared correlate_impl come from Tracking_core.
class Pilot_tracker : public Tracking_core
{
public:
    using Tracking_core::Tracking_core; // inherit the (code, sample_rate, sig, secondary) ctor

    void   run_loops( bool bit_sync, bool sw_loop, Satellite_id prn ) override;
    double code_phase_offset_s() const override; // code phase + the DE subcarrier refinement (Hodgart Eq.4)

protected:
    void correlate( const Sample_block& block, int n ) override; // double-estimator (factored BOC) replica
    void configure_taps( double ci ) override;                   // DE 7-tap: P, code E/L (x in-phase/quadrature subcarrier),
                                                                 // subcarrier E/L
    double code_error() const override;                          // subcarrier-phase-independent envelope code DLL

private:
    // Double-estimator SUBCARRIER loop (SLL): drives subcarrier_offset_ from the subcarrier E/L gates so the
    // subcarrier tracks the signal independently of the code phase. It may lock onto any subcarrier lobe; the
    // integer ambiguity is resolved in code_phase_offset_s() by rounding against the code phase (NOT here - so
    // do NOT clamp/wrap subcarrier_offset_). Run only once carrier+secondary synced.
    void                    subcarrier_update();
    static constexpr double SUBC_GAIN = 0.2; // 1st-order SLL gain (code elements per unit discriminator)

    // Double-estimator state (BOC-only, so it lives here, not in the shared core). subcarrier_offset_ is the
    // subcarrier-minus-code delay (SLL state); tap_sc_offset_ the per-tap subcarrier offsets (set by
    // configure_taps); ext_count_ counts epochs in the current extended-integration window.
    double                       subcarrier_offset_ = 0.0;
    std::array<double, MAX_TAPS> tap_sc_offset_ {};
    int                          ext_count_ = 0;
};

#pragma once
#include "acquisition.h" // Acquisition_result
#include "types.h"       // Sample_block, Satellite_id

// Per-epoch output of a tracker's correlate step (consumed by the Channel + nav decoder).
struct Tracking_output
{
    bool   valid;
    int    samples_consumed;
    double prompt_i;      // current epoch raw prompt (data channel)
    double prompt_i_prev; // previous epoch raw prompt (for nav bit sync)
};

// Abstract carrier+code tracker for one channel. The Channel drives it each epoch:
//   compute_samples_needed() -> correlate_epoch(block) -> [navigation] -> run_loops().
// Concrete strategies share a common NCO/correlation core (Tracking_core):
//   - Costas_tracker : data-bearing BPSK signals (GPS L1 C/A, L2C, BeiDou B1I, ...).
//   - Pilot_tracker  : pilot+data BOC signals (Galileo E1, GPS L1C, BeiDou B1C, ...).
// Constructed per (signal, SV) via Signal::make_tracker().
class Tracker
{
public:
    virtual ~Tracker() = default;

    // Initialise carrier/code NCOs from an acquisition result.
    virtual void initialise( const Acquisition_result& acq ) = 0;

    // Samples the next epoch needs (one code period at the current code frequency).
    virtual int compute_samples_needed() const = 0;

    // Correlate one epoch; block must hold >= compute_samples_needed() samples.
    virtual Tracking_output correlate_epoch( const Sample_block& block ) = 0;

    // Run the loop filters for this epoch (bit_sync / sw_loop come from the nav decoder).
    virtual void run_loops( bool bit_sync, bool sw_loop, Satellite_id prn ) = 0;

    virtual bool   has_lock() const     = 0;
    virtual double get_cn0_db_hz() const = 0; // live tracking C/N0 (M2M4); 0 until the first window

    virtual double get_fractional_chip_time() const = 0;

    virtual double get_code_freq() const            = 0;
    virtual double get_carrier_doppler_hz() const   = 0;
    virtual double get_carrier_acceleration() const = 0;

    // Raw prompt correlator output this epoch (the I/Q constellation point). For history/graphs.
    virtual double get_prompt_i() const = 0;
    virtual double get_prompt_q() const = 0;
};

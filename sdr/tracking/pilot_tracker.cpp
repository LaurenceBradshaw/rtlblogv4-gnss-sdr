#include "pilot_tracker.h"
#include <cmath>
#include "logging.h"

// VEML 5-tap layout (code-element offsets from the prompt): P[0], E[1]/L[2] at +/- the
// resolvable Costas spacing (~1 sample), and the wide VE[3]/VL[4] at +/- 0.6 ranging chip.
// The wide VE/VL span the BOC autocorrelation so the DLL is pulled to the main peak rather
// than a side peak. (0.6 chip * code_elements_per_chip_ -> code elements; ~2-3 samples here.)
void Pilot_tracker::configure_taps( double ci )
{
    const double s       = corr_spacing_ * ci;            // E/L: ~1 sample (resolvable at 4 MHz)
    const double vs      = 0.6 * code_elements_per_chip_; // VE/VL: 0.6 ranging chip
    n_taps_              = 5;
    tap_offset_chips_[0] = 0.0; // P
    tap_offset_chips_[1] = -s;  // E
    tap_offset_chips_[2] = s;   // L
    tap_offset_chips_[3] = -vs; // VE
    tap_offset_chips_[4] = vs;  // VL
}

// VEMLP discriminator (GNSS-SDR dll_nc_vemlp_normalized): combine (VE,E) into the early
// power and (L,VL) into the late power, then normalised early-minus-late.
double Pilot_tracker::code_error() const
{
    const double early = std::sqrt(
        sum_I_[3] * sum_I_[3] + sum_Q_[3] * sum_Q_[3] + // VE
        sum_I_[1] * sum_I_[1] + sum_Q_[1] * sum_Q_[1]
    ); // E
    const double late = std::sqrt(
        sum_I_[2] * sum_I_[2] + sum_Q_[2] * sum_Q_[2] + // L
        sum_I_[4] * sum_I_[4] + sum_Q_[4] * sum_Q_[4]
    ); // VL
    return ( early - late ) / ( early + late + 1e-10 );
}

// run_loops
// Pilot strategy (mirrors GNSS-SDR Galileo_E1 DLL_PLL VEML carrier loop, track_pilot=true):
// pull in with Costas + atan FLL while searching the secondary-code phase; once the
// secondary is synced, wipe it off (per-epoch polarity) and run a pure 4-quadrant PLL with
// no FLL on the now data-free, phase-coherent pilot.
void Pilot_tracker::run_loops( bool /*bit_sync*/, bool /*sw_loop*/, Satellite_id /*prn*/ )
{
    int                      polarity = 1;
    const Tracking_loop_prm* prm      = nullptr;
    bool                     use_fll  = true;
    bool                     pure_pll = false;

    if( secondary_sync_ )
    {
        const float chip = secondary_[secondary_index_];
        polarity         = ( chip >= 0.0f ? 1 : -1 ) * secondary_polarity_;
        secondary_index_ = ( secondary_index_ + 1 ) % static_cast<int>( secondary_.size() );
        prm              = &prm2_; // narrow, pure PLL
        // Keep the FLL assisting the pure PLL: on the now data-free, secondary-wiped pilot the
        // atan(I/Q) FLL gives a clean frequency error, and it holds lock through the residual
        // frequency that a phase-only PLL can lose (the cause of the intermittent unlock seen
        // on weaker SVs after secondary sync). FLL-assisted PLL is the standard robust choice.
        use_fll  = true;
        pure_pll = true;
    }
    else
    {
        advance_secondary_sync(); // uses this epoch's raw prompt (before cumsum below)
        // Grab frequency with the strong FLL (prm1) for ~0.8 s, then HOLD with the tighter
        // prm2 (Costas+FLL) so the prompt signs are clean enough for the secondary to sync.
        prm = ( epoch_count_ < 200 ) ? &prm1_ : &prm2_;
    }

    cumsum_corr( polarity );

    pll_update( *prm, epoch_period_, use_fll, pure_pll );
    dll_update( *prm, epoch_period_ );
    clear_cumsum();

    ++epoch_count_;
}

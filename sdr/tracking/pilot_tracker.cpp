#include "pilot_tracker.h"
#include <cmath>
#include "logging.h"

// DOUBLE-ESTIMATOR (Hodgart/Blunt) 7-tap layout. The BOC replica factors as primary[chip] * subcarrier(elem),
// so we track the code and subcarrier as two independent delays (de_mode_ in Tracking_core does the factoring).
// The CODE (envelope) gates are correlated against BOTH the in-phase subcarrier AND a QUADRATURE subcarrier
// (shifted T_s/2 = 0.5 element = a quarter subcarrier period, expressed here as a +0.5 element subcarrier
// offset). Combining the two as sqrt(I^2+Q^2) makes the envelope SUBCARRIER-PHASE-INDEPENDENT, so the code DLL
// tracks the true primary triangle no matter where the subcarrier sits - this is what stops the code loop
// settling into the false 2-D equilibrium (code at the half-lobe) that biased L1C by ~0.5 element.
//   P[0]            - code@prompt, subcarrier@prompt: full BOC prompt for carrier/data.
//   cEi[1]/cEq[2]   - CODE early, subcarrier in-phase / quadrature.
//   cLi[3]/cLq[4]   - CODE late,  subcarrier in-phase / quadrature.
//   sE[5]/sL[6]     - SUBCARRIER early/late (code@prompt): the precise SLL; locks to any lobe, the integer
//                     ambiguity is resolved against the code phase in code_phase_offset_s().
// Double-estimator correlator: the BOC replica is factored into primary[chip] x subcarrier(elem) with the
// code and subcarrier at independent delays (subcarrier_offset_ is the subcarrier-minus-code delay).
void Pilot_tracker::correlate( const Sample_block& block, int n )
{
    correlate_impl(
        block,
        n,
        Boc_de_replica { code_.data(), code_len_, tap_offset_chips_.data(), tap_sc_offset_.data(), subcarrier_offset_ }
    );
}

void Pilot_tracker::configure_taps( double ci )
{
    const double s = corr_spacing_ * ci; // code E/L spacing (on the wide ~1-chip primary triangle)
    const double q = 0.5;                // quadrature subcarrier shift (T_s/2, quarter subcarrier period)
    // Subcarrier E/L spacing must be well inside the subcarrier correlation half-width (T_s/2 = 0.5 element:
    // |R_sc| is a triangle with its apex every 1 element and zeros every 0.5 element). 0.3 element keeps both
    // gates on the apex slope for a clean, high-gain discriminator.
    const double s_sc    = 0.3;
    n_taps_              = 7;
    tap_offset_chips_[0] = 0.0;
    tap_sc_offset_[0]    = 0.0; // P
    tap_offset_chips_[1] = -s;
    tap_sc_offset_[1]    = 0.0; // cE in-phase
    tap_offset_chips_[2] = -s;
    tap_sc_offset_[2]    = q; // cE quadrature
    tap_offset_chips_[3] = s;
    tap_sc_offset_[3]    = 0.0; // cL in-phase
    tap_offset_chips_[4] = s;
    tap_sc_offset_[4]    = q; // cL quadrature
    tap_offset_chips_[5] = 0.0;
    tap_sc_offset_[5]    = -s_sc; // sE (subcarrier early)
    tap_offset_chips_[6] = 0.0;
    tap_sc_offset_[6]    = s_sc; // sL (subcarrier late)
}

// CODE (envelope) discriminator = normalised early-minus-late on the SUBCARRIER-PHASE-INDEPENDENT envelopes:
// each code gate's power sums in-phase + quadrature subcarrier (and both carrier phases), so sqrt(.) is the
// primary triangle regardless of subcarrier alignment -> a single peak, no false code lock.
double Pilot_tracker::code_error() const
{
    const double early = std::sqrt(
        sum_I_[1] * sum_I_[1] + sum_Q_[1] * sum_Q_[1] + sum_I_[2] * sum_I_[2] + sum_Q_[2] * sum_Q_[2]
    ); // |cE| over in-phase+quadrature subcarrier
    const double late =
        std::sqrt( sum_I_[3] * sum_I_[3] + sum_Q_[3] * sum_Q_[3] + sum_I_[4] * sum_I_[4] + sum_Q_[4] * sum_Q_[4] ); // |cL|
    return ( early - late ) / ( early + late + 1e-10 );
}

// SUBCARRIER discriminator + SLL = normalised early-minus-late on the SUBCARRIER gates (sE/sL). Non-coherent
// magnitudes are fine: the loop only has to lock onto the NEAREST subcarrier lobe (its integer-T_s ambiguity is
// resolved against the code phase in code_phase_offset_s()), so subcarrier_offset_ is deliberately NOT
// clamped/wrapped - constraining it would fight the natural lobe lock and corrupt the precise estimate.
void Pilot_tracker::subcarrier_update()
{
    const double early = std::hypot( sum_I_[5], sum_Q_[5] ); // sE
    const double late  = std::hypot( sum_I_[6], sum_Q_[6] ); // sL
    const double d     = ( early - late ) / ( early + late + 1e-10 );
    subcarrier_offset_ -= SUBC_GAIN * d; // early>late (d>0): replica subcarrier is late -> advance it (earlier)
}

// Double-estimator combine (Hodgart Eq.4). The SUBCARRIER delay is the precise estimate but is ambiguous mod
// T_s (= 1 code element here: 2 elements/chip, T_s = T_c/2). subcarrier_offset_ is the subcarrier-minus-code
// delay and may have locked onto any lobe (= n*T_s + epsilon). The code phase is unambiguous, so round the
// integer-T_s part away and keep only the precise sub-element refinement eps = subcarrier_offset_ -
// round(subcarrier_offset_) (T_s = 1 element), added to the code phase. Correct as long as |code error| <
// T_s/2 (Hodgart Eq.5), which the code DLL holds.
double Pilot_tracker::code_phase_offset_s() const
{
    const double eps = subcarrier_offset_ - std::round( subcarrier_offset_ );
    return ( remaining_code_ + eps ) / code_rate_;
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

    // Verify a completed sync: a CORRECT secondary phase lets the pure PLL clean the pilot, so cos(2*phi)
    // climbs within a few seconds. A FALSE sync (picked at medium C/N0) leaves it stuck low forever. If it
    // hasn't cleaned up after the grace window, drop the sync and re-search - this both kills the stuck
    // "energy on Q" state and gives the SV more chances to sync on the right phase.
    // Slow check: a dirty pre-sync lock that still hasn't cleaned up past the grace window was a false sync.
    const bool slow_false =
        epoch_count_ - epoch_at_sync_ > SECONDARY_SYNC_VERIFY_EPOCHS && carrier_lock_test_ < SECONDARY_SYNC_MIN_LOCK;
    // Fast check: the lock was clean BEFORE this sync but the sync COLLAPSED it - a false sync on a long/weak
    // overlay (GPS L1C's 1800-chip overlay). Catch it quickly so the clean Costas lock isn't degraded for long.
    const bool fast_broke = cos2phi_at_sync_ > SECONDARY_SYNC_WAS_CLEAN &&
                            epoch_count_ - epoch_at_sync_ > SECONDARY_SYNC_FAST_EPOCHS &&
                            carrier_lock_test_ < cos2phi_at_sync_ - SECONDARY_SYNC_DEGRADE;
    if( secondary_sync_ && ( slow_false || fast_broke ) )
    {
        secondary_sync_ = false;
        sec_i_hist_.clear();
        sec_q_hist_.clear();
        ext_count_ = 0;
    }

    if( secondary_sync_ )
    {
        const float chip = secondary_[secondary_index_];
        polarity         = ( chip >= 0.0f ? 1 : -1 ) * secondary_polarity_;
        secondary_index_ = ( secondary_index_ + 1 ) % static_cast<int>( secondary_.size() );
        cumsum_corr( polarity ); // accumulate the secondary-wiped correlators

        if( carrier_lock_test_ > SECONDARY_SYNC_MIN_LOCK )
        {
            // EXTENDED coherent integration: the pilot is cleanly locked, so keep accumulating over
            // EXTEND_SYMBOLS epochs and run a NARROW-bandwidth FLL-assisted pure PLL once per window.
            // The coherent SNR boost steadies the medium-C/N0 lock that otherwise churns; the FLL nulls
            // the small residual frequency that would otherwise rotate the held-NCO prompt within a window.
            if( ++ext_count_ >= EXTEND_SYMBOLS )
            {
                const double dt = EXTEND_SYMBOLS * epoch_period_;
                pll_update( prm2_, dt, /*use_fll=*/true, /*pure_pll=*/true );
                dll_update( prm2_, dt );
                subcarrier_update();
                clear_cumsum();
                ext_count_ = 0;
            }
            // else: still filling the window - hold the NCO, update nothing this epoch.
        }
        else
        {
            // Just synced / still cleaning up: per-epoch FLL-assisted pure PLL pulls the pilot clean
            // (FLL holds the residual frequency a phase-only PLL could lose), before extending.
            ext_count_ = 0;
            pll_update( prm2_, epoch_period_, /*use_fll=*/true, /*pure_pll=*/true );
            dll_update( prm2_, epoch_period_ );
            subcarrier_update();
            clear_cumsum();
        }
    }
    else
    {
        advance_secondary_sync(); // uses this epoch's raw prompt (before cumsum below)
        if( secondary_sync_ )
        {
            epoch_at_sync_   = epoch_count_;       // record when (re-)sync happened, for the verify above
            cos2phi_at_sync_ = carrier_lock_test_; // baseline lock for the degradation check above
            ext_count_       = 0;
        }
        // Grab frequency with the strong FLL (prm1) for ~0.8 s, then HOLD with the tighter
        // prm2 (Costas+FLL) so the prompt signs are clean enough for the secondary to sync.
        prm = ( epoch_count_ < 200 ) ? &prm1_ : &prm2_;
        cumsum_corr( polarity );
        pll_update( *prm, epoch_period_, use_fll, pure_pll );
        dll_update( *prm, epoch_period_ );
        clear_cumsum();
    }

    ++epoch_count_;
}

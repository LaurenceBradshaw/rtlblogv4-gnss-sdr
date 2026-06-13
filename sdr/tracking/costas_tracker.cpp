#include "costas_tracker.h"
#include <cmath>

namespace
{
// Normalised early-minus-late envelope discriminator (GNSS-SDRLIB dll()). pE, pL are the early/late
// correlator magnitudes. Result in (-1, 1): >0 when the code replica is late (early stronger), 0 at
// the correlation peak (pE == pL). The 1e-10 guards the all-zero (no-signal) case.
double normalised_early_late( double pE, double pL )
{
    return ( pE - pL ) / ( pE + pL + 1e-10 );
}
} // namespace

// 3-tap correlator: Prompt[0], Early[1] at -s, Late[2] at +s (s = corr_spacing_ samples).
void Costas_tracker::configure_taps( double ci )
{
    const double s       = corr_spacing_ * ci; // E/L spacing in code elements
    n_taps_              = 3;
    tap_offset_chips_[0] = 0.0;
    tap_offset_chips_[1] = -s;
    tap_offset_chips_[2] = s;
}

// Normalised early-minus-late envelope discriminator. mirrors GNSS-SDRLIB dll().
double Costas_tracker::code_error() const
{
    const double pE = std::sqrt( sum_I_[1] * sum_I_[1] + sum_Q_[1] * sum_Q_[1] );
    const double pL = std::sqrt( sum_I_[2] * sum_I_[2] + sum_Q_[2] * sum_Q_[2] );
    return normalised_early_late( pE, pL );
}

// run_loops
// Costas/FLL strategy for data-bearing signals. mirrors GNSS-SDRLIB sdrthread():
//   cumsumcorr(polarity=1) -> before bit sync run prm1 every epoch; after bit sync run
//   prm2 every loopms (sw_loop); clear the accumulators when a loop ran.
//   GPS (fll_active_): cross/dot FLL, bit_sync/sw_loop driven.
//   Galileo E1-B data (1 symbol/epoch, !fll_active_): atan FLL, time-based prm1 -> prm2.
void Costas_tracker::run_loops( bool bit_sync, bool sw_loop, Satellite_id /*prn*/ )
{
    const Tracking_loop_prm* prm = nullptr;
    double                   dt  = epoch_period_;

    if( !fll_active_ )
    {
        // E1-B data channel (no pilot): Costas + atan FLL, time-based prm1 -> prm2 (~2 s).
        constexpr int E1_PULLIN_EPOCHS = 500;
        prm                            = ( epoch_count_ < E1_PULLIN_EPOCHS ) ? &prm1_ : &prm2_;
    }
    // GPS path (matches GNSS-SDRLIB sdrthread()): before bit sync run every epoch with prm1;
    // after bit sync run every loopms=10 ms (sw_loop) with prm2.
    else if( !bit_sync )
    {
        prm = &prm1_;
    }
    else if( sw_loop )
    {
        prm = &prm2_;
        dt  = 1e-2;
    }

    cumsum_corr( 1 ); // ocode is all 1s for these signals

    if( prm != nullptr )
    {
        pll_update( *prm, dt, /*use_fll=*/true, /*pure_pll=*/false );
        dll_update( *prm, dt );
        clear_cumsum();
    }

    ++epoch_count_;
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "dll_normalised_early_late_discriminator", "[tracking][dll]" )
{
    // Peak: equal early/late envelopes -> zero error.
    REQUIRE( normalised_early_late( 5.0, 5.0 ) == Catch::Approx( 0.0 ) );
    // Early stronger (replica late) -> positive; late stronger -> negative; antisymmetric.
    REQUIRE( normalised_early_late( 3.0, 1.0 ) == Catch::Approx( 0.5 ) );
    REQUIRE( normalised_early_late( 1.0, 3.0 ) == Catch::Approx( -0.5 ) );
    // Bounded in [-1, 1] at the extremes.
    REQUIRE( normalised_early_late( 1.0, 0.0 ) == Catch::Approx( 1.0 ).margin( 1e-9 ) );
    REQUIRE( normalised_early_late( 0.0, 1.0 ) == Catch::Approx( -1.0 ).margin( 1e-9 ) );
    // No signal -> guarded to 0 (no divide-by-zero).
    REQUIRE( normalised_early_late( 0.0, 0.0 ) == Catch::Approx( 0.0 ) );
}
#endif

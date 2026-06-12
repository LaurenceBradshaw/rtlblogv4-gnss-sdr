#include "beidou_b1i_code.h"

#include <cmath>
#include <stdexcept>

namespace beidou
{

namespace
{
// G2 phase assignment: each PRN selects two G2 register taps (1-indexed); the G2
// output is their product (== XOR in +/-1 arithmetic). BeiDou SIS-ICD Table; identical
// to GNSS-SDRLIB gencode_B1IB2I phase[37][2].
constexpr int PHASE[B1i_code::MAX_PRN][2] = {
    { 1, 3 },  { 1, 4 },  { 1, 5 },  { 1, 6 },  { 1, 8 },  { 1, 9 },  { 1, 10 },  { 1, 11 }, { 2, 7 }, { 3, 4 },
    { 3, 5 },  { 3, 6 },  { 3, 8 },  { 3, 9 },  { 3, 10 }, { 3, 11 }, { 4, 5 },   { 4, 6 },  { 4, 8 }, { 4, 9 },
    { 4, 10 }, { 4, 11 }, { 5, 6 },  { 5, 8 },  { 5, 9 },  { 5, 10 }, { 5, 11 },  { 6, 8 },  { 6, 9 }, { 6, 10 },
    { 6, 11 }, { 8, 9 },  { 8, 10 }, { 8, 11 }, { 9, 10 }, { 9, 11 }, { 10, 11 },
};

// Generate the 2046-chip B1I primary code as +/-1 floats.
// Mirrors GNSS-SDRLIB gencode_B1IB2I:
//   G1 feedback poly:  1 + x + x^7 + x^8 + x^9 + x^10 + x^11   -> taps R[0,6,7,8,9,10]
//   G2 feedback poly:  1 + x + x^2 + x^3 + x^4 + x^5 + x^8 + x^9 + x^11 -> R[0,1,2,3,4,7,8,10]
//   init phase: 01010101010 (stored as +/-1 with bit 1->-1, bit 0->+1)
//   chip = -G1 * G2     (registers reset each period -> truncated to 2046 chips)
std::vector<float> generate_b1i_chips( Satellite_id prn_id )
{
    if( prn_id < 1 || prn_id > static_cast<Satellite_id>( B1i_code::MAX_PRN ) )
    {
        throw std::out_of_range( "BeiDou B1I PRN out of range (1.." + std::to_string( B1i_code::MAX_PRN ) + ")" );
    }
    const int t0 = PHASE[prn_id - 1][0];
    const int t1 = PHASE[prn_id - 1][1];

    // +/-1 registers; init 01010101010 -> {+1,-1,+1,-1,...} (bit 0->+1, bit 1->-1).
    int R1[11] = { 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1 };
    int R2[11] = { 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1 };

    std::vector<float> code( B1i_code::PERIOD_CHIPS );
    for( int i = 0; i < B1i_code::PERIOD_CHIPS; i++ )
    {
        const int g1 = R1[10];
        const int g2 = R2[t0 - 1] * R2[t1 - 1];
        code[i]      = static_cast<float>( -g1 * g2 );

        const int c1 = R1[0] * R1[6] * R1[7] * R1[8] * R1[9] * R1[10];
        const int c2 = R2[0] * R2[1] * R2[2] * R2[3] * R2[4] * R2[7] * R2[8] * R2[10];
        for( int j = 10; j > 0; j-- )
        {
            R1[j] = R1[j - 1];
            R2[j] = R2[j - 1];
        }
        R1[0] = c1;
        R2[0] = c2;
    }
    return code;
}
} // namespace

Complex_buf B1i_code::generate( Satellite_id prn_id, uint32_t sampling_rate )
{
    const int  samples_per_period = static_cast<int>( std::round( PERIOD_TIME * sampling_rate ) );
    const auto chip_levels        = generate_b1i_chips( prn_id );

    Complex_buf out( samples_per_period );

    // Resample chips -> samples by accumulating chip phase, flooring to a chip index.
    // Same convention as gps::L1ca_code::generate (rescode in GNSS-SDRLIB).
    double       coff = 0.0;
    const double ci   = static_cast<double>( PERIOD_CHIPS ) / samples_per_period;
    for( int n = 0; n < samples_per_period; n++, coff += ci )
    {
        if( coff >= PERIOD_CHIPS )
        {
            coff -= PERIOD_CHIPS;
        }
        out[n] = Complex_sample( chip_levels[static_cast<int>( coff )], 0.0f );
    }
    return out;
}

std::vector<float> B1i_code::chips( Satellite_id prn_id )
{
    return generate_b1i_chips( prn_id );
}

} // namespace beidou

// Unit tests
#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <numeric>

// These tests validate the generator via the *structural* properties every valid Gold
// ranging code must have - a sharp periodic autocorrelation peak with bounded sidelobes,
// bounded cross-correlation between SVs, near-balance, and SV distinctness. They are
// deliberately NOT compared against chips emitted by the same generator (that would be
// circular). A wrong feedback polynomial, tap pair, or init state destroys these
// properties immediately. The algorithm mirrors the BeiDou SIS-ICD (and GNSS-SDRLIB
// gencode_B1IB2I) - see generate_b1i_chips above.
//
// Bound: for Gold codes from 11-stage registers the cross/auto sidelobes are three-
// valued with magnitude <= t(11) = 2^6 + 1 = 65. B1I truncates the underlying sequences
// from 2047->2046 chips, which perturbs this (measured worst ~140). We require every
// off-peak correlation and the DC level to stay under 10% of the 2046 main peak - a
// usable-spreading-code threshold the truncated code clears comfortably (~7%) while a
// wrong polynomial / tap / init collapses the autocorrelation far past it.

namespace
{
constexpr int B1I_N            = beidou::B1i_code::PERIOD_CHIPS; // 2046
constexpr int B1I_SIDELOBE_MAX = B1I_N / 10;                     // 204 - 10% of the main peak

// Periodic correlation of two +/-1 chip sequences at a given circular shift.
int b1i_periodic_corr( const std::vector<float>& a, const std::vector<float>& b, int tau )
{
    int acc = 0;
    for( int i = 0; i < B1I_N; i++ )
    {
        acc += static_cast<int>( a[i] ) * static_cast<int>( b[( i + tau ) % B1I_N] );
    }
    return acc;
}
} // namespace

TEST_CASE( "B1i_code_length_and_values", "[beidou][code]" )
{
    const auto c = beidou::B1i_code::chips( 6 );
    REQUIRE( c.size() == static_cast<size_t>( B1I_N ) );
    for( const float v : c )
    {
        REQUIRE( ( v == 1.0f || v == -1.0f ) );
    }
}

TEST_CASE( "B1i_prn_range_bounds", "[beidou][code]" )
{
    REQUIRE_NOTHROW( beidou::B1i_code::chips( 1 ) );
    REQUIRE_NOTHROW( beidou::B1i_code::chips( beidou::B1i_code::MAX_PRN ) );
    REQUIRE_THROWS( beidou::B1i_code::chips( 0 ) );
    REQUIRE_THROWS( beidou::B1i_code::chips( beidou::B1i_code::MAX_PRN + 1 ) );
}

TEST_CASE( "B1i_codes_near_balanced", "[beidou][code]" )
{
    // A balanced spreading code has near-equal +/-1 counts; |sum| must be tiny vs N=2046.
    for( int prn = 1; prn <= beidou::B1i_code::MAX_PRN; prn++ )
    {
        const auto c   = beidou::B1i_code::chips( prn );
        const int  sum = static_cast<int>( std::accumulate( c.begin(), c.end(), 0.0f ) );
        INFO( "PRN " << prn << " chip sum = " << sum );
        REQUIRE( std::abs( sum ) <= B1I_SIDELOBE_MAX );
    }
}

TEST_CASE( "B1i_autocorrelation_peak_and_sidelobes", "[beidou][code]" )
{
    for( const int prn : { 1, 6, 19, 37 } )
    {
        const auto c = beidou::B1i_code::chips( prn );
        REQUIRE( b1i_periodic_corr( c, c, 0 ) == B1I_N ); // perfect peak at zero lag

        int max_sidelobe = 0;
        for( int tau = 1; tau < B1I_N; tau++ )
        {
            max_sidelobe = std::max( max_sidelobe, std::abs( b1i_periodic_corr( c, c, tau ) ) );
        }
        INFO( "PRN " << prn << " max abs autocorr sidelobe = " << max_sidelobe );
        REQUIRE( max_sidelobe <= B1I_SIDELOBE_MAX );
    }
}

TEST_CASE( "B1i_cross_correlation_bounded", "[beidou][code]" )
{
    const auto a = beidou::B1i_code::chips( 1 );
    for( const int prn : { 2, 6, 19, 37 } )
    {
        const auto b = beidou::B1i_code::chips( prn );

        int max_cross = 0;
        for( int tau = 0; tau < B1I_N; tau++ )
        {
            max_cross = std::max( max_cross, std::abs( b1i_periodic_corr( a, b, tau ) ) );
        }
        INFO( "PRN 1 x PRN " << prn << " max abs cross-corr = " << max_cross );
        REQUIRE( max_cross <= B1I_SIDELOBE_MAX );
    }
}

TEST_CASE( "B1i_codes_distinct", "[beidou][code]" )
{
    for( int prn = 2; prn <= beidou::B1i_code::MAX_PRN; prn++ )
    {
        REQUIRE( beidou::B1i_code::chips( prn ) != beidou::B1i_code::chips( 1 ) );
    }
}

TEST_CASE( "B1i_resampling_one_period", "[beidou][code]" )
{
    constexpr uint32_t fs       = 4'000'000;
    const auto         samples  = beidou::B1i_code::generate( 6, fs );
    const auto         expected = static_cast<size_t>( std::lround( beidou::B1i_code::PERIOD_TIME * fs ) );
    REQUIRE( samples.size() == expected ); // 4000 samples for 1 ms @ 4 MS/s
    for( const auto& s : samples )
    {
        REQUIRE( s.imag() == 0.0f );
        REQUIRE( ( s.real() == 1.0f || s.real() == -1.0f ) );
    }
}
#endif // ENABLE_UNIT_TESTS

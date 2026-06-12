#include "gps_l1ca_code.h"

#include <vector>

namespace gps
{
static std::vector<int> generate_ca_bits( Satellite_id prn_id )
{
    auto shift = []( std::vector<int>& reg, const std::vector<int>& feedback, const std::vector<int>& output ) -> int
    {
        int out = 0;
        for( int o : output )
        {
            out ^= reg[o - 1];
        }
        int fb = 0;
        for( int f : feedback )
        {
            fb ^= reg[f - 1];
        }
        for( int i = (int) reg.size() - 1; i >= 1; --i )
        {
            reg[i] = reg[i - 1];
        }
        reg[0] = fb;
        return out;
    };

    std::vector<std::vector<int>> tap_list = { { 2, 6 }, { 3, 7 },  { 4, 8 }, { 5, 9 }, { 1, 9 },  { 2, 10 }, { 1, 8 },
                                               { 2, 9 }, { 3, 10 }, { 2, 3 }, { 3, 4 }, { 5, 6 },  { 6, 7 },  { 7, 8 },
                                               { 8, 9 }, { 9, 10 }, { 1, 4 }, { 2, 5 }, { 3, 6 },  { 4, 7 },  { 5, 8 },
                                               { 6, 9 }, { 1, 3 },  { 4, 6 }, { 5, 7 }, { 6, 8 },  { 7, 9 },  { 8, 10 },
                                               { 1, 6 }, { 2, 7 },  { 3, 8 }, { 4, 9 }, { 5, 10 }, { 4, 10 }, { 1, 7 },
                                               { 2, 8 }, { 4, 10 } };

    std::vector<int> G1( 10, 1 ), G2( 10, 1 );
    std::vector<int> ca;
    auto             tap = tap_list[prn_id - 1];
    for( int i = 0; i < 1023; i++ )
    {
        int g1 = shift( G1, { 3, 10 }, { 10 } );
        int g2 = shift( G2, { 2, 3, 6, 8, 9, 10 }, tap );
        ca.push_back( ( g1 + g2 ) % 2 );
    }
    return ca;
}

Complex_buf L1ca_code::generate( Satellite_id prn_id, uint32_t sampling_rate )
{
    const int  samples_per_period = static_cast<int>( std::round( PERIOD_TIME * sampling_rate ) );
    const auto ca_bits            = generate_ca_bits( prn_id );

    std::vector<int> ca_levels( PERIOD_CHIPS );
    for( int i = 0; i < PERIOD_CHIPS; i++ )
    {
        ca_levels[i] = ca_bits[i] ? 1 : -1;
    }

    Complex_buf out( samples_per_period );

    // mirrors GNSS-SDRLIB rescode(code, 1023, coff=0, smax=0, ci, nsamp, rcode):
    //   for each sample: *p = code[(int)coff]; coff += ci;
    // Uses floor (integer cast), not round, and accumulates coff to match chip boundaries.
    double       coff = 0.0;
    const double ci   = static_cast<double>( PERIOD_CHIPS ) / samples_per_period;
    for( int n = 0; n < samples_per_period; n++, coff += ci )
    {
        if( coff >= PERIOD_CHIPS )
            coff -= PERIOD_CHIPS;
        out[n] = Complex_sample( static_cast<float>( ca_levels[static_cast<int>( coff )] ), 0.0f );
    }

    return out;
}

std::vector<float> L1ca_code::chips( Satellite_id prn_id )
{
    const auto         bits = generate_ca_bits( prn_id );
    std::vector<float> out( PERIOD_CHIPS );
    for( int i = 0; i < PERIOD_CHIPS; i++ )
    {
        out[i] = bits[i] ? 1.0f : -1.0f;
    }
    return out;
}

} // namespace gps

// Unit tests
#ifdef ENABLE_UNIT_TESTS
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <numeric>

// GPS L1 C/A is a full (untruncated) 1023-chip Gold code, so the structural bounds are
// exact: it is balanced (|sum| == 1), and the periodic auto/cross-correlation is three-
// valued {-65, -1, 63} -> magnitude <= 65. On top of that we pin the generator against
// an INDEPENDENT published reference: the IS-GPS-200 Table 3-Ia "first 10 chips" octal
// values. Matching those is a non-circular check that PRN -> code is exactly correct.

namespace
{
constexpr int GPS_N        = gps::L1ca_code::PERIOD_CHIPS; // 1023
constexpr int GPS_MAX_PRN  = 32;
constexpr int GPS_CORR_MAX = 65; // max |off-peak| for an n=10 Gold code: 2^6 + 1

// Periodic correlation of two +/-1 chip sequences at a given circular shift.
int gps_periodic_corr( const std::vector<float>& a, const std::vector<float>& b, int tau )
{
    int acc = 0;
    for( int i = 0; i < GPS_N; i++ )
    {
        acc += static_cast<int>( a[i] ) * static_cast<int>( b[( i + tau ) % GPS_N] );
    }
    return acc;
}
} // namespace

TEST_CASE( "gps_l1ca_code_length_and_values", "[gps][code]" )
{
    const auto c = gps::L1ca_code::chips( 7 );
    REQUIRE( c.size() == static_cast<size_t>( GPS_N ) );
    for( const float v : c )
    {
        REQUIRE( ( v == 1.0f || v == -1.0f ) );
    }
}

TEST_CASE( "gps_l1ca_first_ten_chips_match_icd", "[gps][code]" )
{
    // IS-GPS-200 Table 3-Ia, "First 10 Chips" (octal) - independent published reference.
    constexpr int ref_octal[10] = { 01440, 01620, 01710, 01744, 01133, 01455, 01131, 01454, 01626, 01504 };
    for( int prn = 1; prn <= 10; prn++ )
    {
        const auto c = gps::L1ca_code::chips( prn );
        int        v = 0;
        for( int i = 0; i < 10; i++ )
        {
            v = ( v << 1 ) | ( ( c[i] > 0.0f ) ? 1 : 0 ); // chip[0] is the MSB
        }
        INFO( "PRN " << prn );
        REQUIRE( v == ref_octal[prn - 1] );
    }
}

TEST_CASE( "gps_l1ca_codes_balanced", "[gps][code]" )
{
    // A full-length Gold code has exactly 512 vs 511 chips -> |sum of +/-1| == 1.
    for( int prn = 1; prn <= GPS_MAX_PRN; prn++ )
    {
        const auto c   = gps::L1ca_code::chips( prn );
        const int  sum = static_cast<int>( std::accumulate( c.begin(), c.end(), 0.0f ) );
        INFO( "PRN " << prn << " chip sum = " << sum );
        REQUIRE( std::abs( sum ) == 1 );
    }
}

TEST_CASE( "gps_l1ca_autocorrelation_peak_and_sidelobes", "[gps][code]" )
{
    for( const int prn : { 1, 7, 19, 32 } )
    {
        const auto c = gps::L1ca_code::chips( prn );
        REQUIRE( gps_periodic_corr( c, c, 0 ) == GPS_N ); // perfect peak at zero lag

        int max_sidelobe = 0;
        for( int tau = 1; tau < GPS_N; tau++ )
        {
            max_sidelobe = std::max( max_sidelobe, std::abs( gps_periodic_corr( c, c, tau ) ) );
        }
        INFO( "PRN " << prn << " max abs autocorr sidelobe = " << max_sidelobe );
        REQUIRE( max_sidelobe <= GPS_CORR_MAX );
    }
}

TEST_CASE( "gps_l1ca_cross_correlation_bounded", "[gps][code]" )
{
    const auto a = gps::L1ca_code::chips( 1 );
    for( const int prn : { 2, 7, 19, 32 } )
    {
        const auto b = gps::L1ca_code::chips( prn );

        int max_cross = 0;
        for( int tau = 0; tau < GPS_N; tau++ )
        {
            max_cross = std::max( max_cross, std::abs( gps_periodic_corr( a, b, tau ) ) );
        }
        INFO( "PRN 1 x PRN " << prn << " max abs cross-corr = " << max_cross );
        REQUIRE( max_cross <= GPS_CORR_MAX );
    }
}

TEST_CASE( "gps_l1ca_codes_distinct", "[gps][code]" )
{
    for( int prn = 2; prn <= GPS_MAX_PRN; prn++ )
    {
        REQUIRE( gps::L1ca_code::chips( prn ) != gps::L1ca_code::chips( 1 ) );
    }
}

TEST_CASE( "gps_l1ca_resampling_one_period", "[gps][code]" )
{
    constexpr uint32_t fs       = 4'000'000;
    const auto         samples  = gps::L1ca_code::generate( 7, fs );
    const auto         expected = static_cast<size_t>( std::lround( gps::L1ca_code::PERIOD_TIME * fs ) );
    REQUIRE( samples.size() == expected ); // 4000 samples for 1 ms @ 4 MS/s
    for( const auto& s : samples )
    {
        REQUIRE( s.imag() == 0.0f );
        REQUIRE( ( s.real() == 1.0f || s.real() == -1.0f ) );
    }
}
#endif // ENABLE_UNIT_TESTS

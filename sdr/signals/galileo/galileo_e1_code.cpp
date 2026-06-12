#include "galileo_e1_code.h"

#include <cmath>
#include <stdexcept>
#include "galileo_e1b_hex.h"
#include "galileo_e1c_hex.h"

namespace galileo
{

namespace
{
// Hex character -> 0..15.
int hex_digit( char c )
{
    if( c >= '0' && c <= '9' )
    {
        return c - '0';
    }
    if( c >= 'A' && c <= 'F' )
    {
        return c - 'A' + 10;
    }
    if( c >= 'a' && c <= 'f' )
    {
        return c - 'a' + 10;
    }
    throw std::invalid_argument( "Galileo E1 hex: bad digit" );
}

// Decode a 1023-hex-digit memory code to 4092 +/-1 primary chips. Convention mirrors
// GNSS-SDRLIB hex2bin + the gencode_E1B/E1C inversion: bits are MSB-first within each
// hex digit, and chip = (bit == 1) ? +1 : -1.
std::vector<float> decode_primary_hex( const char* hex )
{
    std::vector<float> code( 4092 );
    for( int i = 0; i < 4092 / 4; i++ )
    {
        const int d = hex_digit( hex[i] );
        for( int k = 0; k < 4; k++ )
        {
            const int bit   = ( d >> ( 3 - k ) ) & 1; // MSB first
            code[4 * i + k] = bit ? 1.0f : -1.0f;
        }
    }
    return code;
}

// BOC(1,1): chip c -> half-chip pair (-c, +c). Mirrors GNSS-SDRLIB boc(code,len,crate,1,1)
// (replicate each chip into 2 half-chips, then negate the even half-chip).
std::vector<float> apply_boc11( const std::vector<float>& primary )
{
    std::vector<float> boc( primary.size() * 2 );
    for( size_t i = 0; i < primary.size(); i++ )
    {
        boc[2 * i]     = -primary[i];
        boc[2 * i + 1] = primary[i];
    }
    return boc;
}

// Resample a BOC half-chip sequence to one period of samples at sample_rate_hz (floor of
// accumulated chip phase). Shared by E1-B / E1-C generate().
Complex_buf resample_boc( const std::vector<float>& boc, double period_s, uint32_t sample_rate_hz )
{
    const int    n   = static_cast<int>( std::round( period_s * sample_rate_hz ) );
    const int    len = static_cast<int>( boc.size() );
    Complex_buf  out( n );
    double       coff = 0.0;
    const double ci   = static_cast<double>( len ) / n;
    for( int i = 0; i < n; i++, coff += ci )
    {
        if( coff >= len )
        {
            coff -= len;
        }
        out[i] = Complex_sample( boc[static_cast<int>( coff )], 0.0f );
    }
    return out;
}
} // namespace

std::vector<float> E1b_code::primary_chips( Satellite_id prn_id )
{
    if( prn_id < 1 || prn_id > static_cast<Satellite_id>( MAX_PRN ) )
    {
        throw std::out_of_range( "Galileo E1B PRN out of range (1.." + std::to_string( MAX_PRN ) + ")" );
    }
    return decode_primary_hex( E1B_HEX[prn_id - 1] );
}

std::vector<float> E1b_code::chips( Satellite_id prn_id )
{
    return apply_boc11( primary_chips( prn_id ) );
}

Complex_buf E1b_code::generate( Satellite_id prn_id, uint32_t sampling_rate )
{
    return resample_boc( chips( prn_id ), PERIOD_TIME, sampling_rate );
}

// E1-C (pilot)
std::vector<float> E1c_code::primary_chips( Satellite_id prn_id )
{
    if( prn_id < 1 || prn_id > static_cast<Satellite_id>( MAX_PRN ) )
    {
        throw std::out_of_range( "Galileo E1C PRN out of range (1.." + std::to_string( MAX_PRN ) + ")" );
    }
    return decode_primary_hex( E1C_HEX[prn_id - 1] );
}

std::vector<float> E1c_code::chips( Satellite_id prn_id )
{
    return apply_boc11( primary_chips( prn_id ) );
}

Complex_buf E1c_code::generate( Satellite_id prn_id, uint32_t sampling_rate )
{
    return resample_boc( chips( prn_id ), PERIOD_TIME, sampling_rate );
}

std::vector<float> E1c_code::secondary_chips()
{
    // CS25 secondary: hex 380AD90 -> first 25 of 28 bits (skip last 3), chip =
    // (bit == 1) ? -1 : +1 (NOT inverted - matches GNSS-SDRLIB gencode_E1CO).
    static constexpr char hex[] = "380AD90";
    std::vector<float>    sec;
    sec.reserve( SECONDARY_CHIPS );
    for( int i = 0; i < 7 && static_cast<int>( sec.size() ) < SECONDARY_CHIPS; i++ )
    {
        const int d = hex_digit( hex[i] );
        for( int k = 0; k < 4 && static_cast<int>( sec.size() ) < SECONDARY_CHIPS; k++ )
        {
            const int bit = ( d >> ( 3 - k ) ) & 1;
            sec.push_back( bit ? -1.0f : 1.0f );
        }
    }
    return sec;
}

} // namespace galileo

// Unit tests
#ifdef ENABLE_UNIT_TESTS
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>

// Galileo E1-B primary codes are MEMORY codes lifted verbatim from the Galileo SIS-ICD
// (via GNSS-SDRLIB), so there is no LFSR/Gold structure to re-derive independently.
// These tests validate (a) the hex decode (length, strictly +/-1), (b) the BOC(1,1)
// transform exactly (each chip -> opposite-sign half-chip pair), (c) that the codes
// still have the sharp autocorrelation / low cross-correlation a usable ranging code
// needs, and (d) SV distinctness + resampling. Memory codes do not satisfy the 3-valued
// Gold bound, so the correlation threshold is empirical (well under the main peak).

using galileo::E1b_code;

namespace
{
constexpr int GAL_PRIMARY_N = E1b_code::PRIMARY_CHIPS; // 4092
constexpr int GAL_CORR_MAX  = 320;                     // memory-code sidelobe ceiling, < 8% of the
                                                       // 4092 peak (measured worst ~208 for these SVs)

int gal_periodic_corr( const std::vector<float>& a, const std::vector<float>& b, int tau )
{
    int acc = 0;
    for( int i = 0; i < GAL_PRIMARY_N; i++ )
    {
        acc += static_cast<int>( a[i] ) * static_cast<int>( b[( i + tau ) % GAL_PRIMARY_N] );
    }
    return acc;
}
} // namespace

TEST_CASE( "galileo_e1b_primary_length_and_values", "[galileo][code]" )
{
    const auto c = E1b_code::primary_chips( 11 );
    REQUIRE( c.size() == static_cast<size_t>( GAL_PRIMARY_N ) );
    for( const float v : c )
    {
        REQUIRE( ( v == 1.0f || v == -1.0f ) );
    }
}

TEST_CASE( "galileo_e1b_boc11_structure", "[galileo][code]" )
{
    // BOC(1,1): each primary chip c becomes the half-chip pair (-c, +c). So the BOC code
    // is twice as long, its odd half-chips reproduce the primary, and every adjacent
    // half-chip pair is opposite-signed.
    const auto prim = E1b_code::primary_chips( 11 );
    const auto boc  = E1b_code::chips( 11 );
    REQUIRE( boc.size() == static_cast<size_t>( E1b_code::BOC_CHIPS ) );
    for( size_t i = 0; i < prim.size(); i++ )
    {
        REQUIRE( boc[2 * i] == -prim[i] );
        REQUIRE( boc[2 * i + 1] == prim[i] );
        REQUIRE( boc[2 * i] == -boc[2 * i + 1] );
    }
}

TEST_CASE( "galileo_e1b_prn_range_bounds", "[galileo][code]" )
{
    REQUIRE_NOTHROW( E1b_code::primary_chips( 1 ) );
    REQUIRE_NOTHROW( E1b_code::primary_chips( E1b_code::MAX_PRN ) );
    REQUIRE_THROWS( E1b_code::primary_chips( 0 ) );
    REQUIRE_THROWS( E1b_code::primary_chips( E1b_code::MAX_PRN + 1 ) );
}

TEST_CASE( "galileo_e1b_autocorrelation_peak_and_sidelobes", "[galileo][code]" )
{
    for( const int prn : { 1, 11, 19 } )
    {
        const auto c = E1b_code::primary_chips( prn );
        REQUIRE( gal_periodic_corr( c, c, 0 ) == GAL_PRIMARY_N ); // perfect peak at zero lag

        int max_sidelobe = 0;
        for( int tau = 1; tau < GAL_PRIMARY_N; tau++ )
        {
            max_sidelobe = std::max( max_sidelobe, std::abs( gal_periodic_corr( c, c, tau ) ) );
        }
        INFO( "PRN " << prn << " max abs autocorr sidelobe = " << max_sidelobe );
        REQUIRE( max_sidelobe <= GAL_CORR_MAX );
    }
}

TEST_CASE( "galileo_e1b_cross_correlation_bounded", "[galileo][code]" )
{
    const auto a = E1b_code::primary_chips( 1 );
    for( const int prn : { 11, 19 } )
    {
        const auto b = E1b_code::primary_chips( prn );

        int max_cross = 0;
        for( int tau = 0; tau < GAL_PRIMARY_N; tau++ )
        {
            max_cross = std::max( max_cross, std::abs( gal_periodic_corr( a, b, tau ) ) );
        }
        INFO( "PRN 1 x PRN " << prn << " max abs cross-corr = " << max_cross );
        REQUIRE( max_cross <= GAL_CORR_MAX );
    }
}

TEST_CASE( "galileo_e1b_codes_distinct", "[galileo][code]" )
{
    for( int prn = 2; prn <= E1b_code::MAX_PRN; prn++ )
    {
        REQUIRE( E1b_code::primary_chips( prn ) != E1b_code::primary_chips( 1 ) );
    }
}

TEST_CASE( "galileo_e1b_resampling_one_period", "[galileo][code]" )
{
    constexpr uint32_t fs       = 4'000'000;
    const auto         samples  = E1b_code::generate( 11, fs );
    const auto         expected = static_cast<size_t>( std::lround( E1b_code::PERIOD_TIME * fs ) );
    REQUIRE( samples.size() == expected ); // 16000 samples for 4 ms @ 4 MS/s
    for( const auto& s : samples )
    {
        REQUIRE( s.imag() == 0.0f );
        REQUIRE( ( s.real() == 1.0f || s.real() == -1.0f ) );
    }
}

// E1-C (pilot)
using galileo::E1c_code;

TEST_CASE( "galileo_e1c_primary_shape_and_boc", "[galileo][code]" )
{
    const auto prim = E1c_code::primary_chips( 7 );
    const auto boc  = E1c_code::chips( 7 );
    REQUIRE( prim.size() == static_cast<size_t>( E1c_code::PRIMARY_CHIPS ) );
    REQUIRE( boc.size() == static_cast<size_t>( E1c_code::BOC_CHIPS ) );
    for( size_t i = 0; i < prim.size(); i++ )
    {
        REQUIRE( ( prim[i] == 1.0f || prim[i] == -1.0f ) );
        REQUIRE( boc[2 * i] == -prim[i] ); // BOC(1,1): chip -> (-c, +c)
        REQUIRE( boc[2 * i + 1] == prim[i] );
    }
}

TEST_CASE( "galileo_e1c_distinct_from_e1b", "[galileo][code]" )
{
    // The pilot and data channels use different memory codes.
    for( const int prn : { 1, 11, 19, 50 } )
    {
        REQUIRE( E1c_code::primary_chips( prn ) != E1b_code::primary_chips( prn ) );
    }
}

TEST_CASE( "galileo_e1c_cs25_secondary_matches_icd", "[galileo][code]" )
{
    // CS25 is the same 25-chip sequence for all SVs. Bit pattern 0011100000001010110110010
    // (Galileo OS SIS-ICD), mapped bit 1 -> -1, bit 0 -> +1. Independent published reference.
    const char* cs25_bits = "0011100000001010110110010";
    const auto  sec       = E1c_code::secondary_chips();
    REQUIRE( sec.size() == static_cast<size_t>( E1c_code::SECONDARY_CHIPS ) );
    for( int i = 0; i < E1c_code::SECONDARY_CHIPS; i++ )
    {
        const float expected = ( cs25_bits[i] == '1' ) ? -1.0f : 1.0f;
        INFO( "CS25 chip " << i );
        REQUIRE( sec[i] == expected );
    }
}
#endif // ENABLE_UNIT_TESTS

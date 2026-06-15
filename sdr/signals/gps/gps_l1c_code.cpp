#include "gps_l1c_code.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include "gps_l1c_tables.h"

namespace gps
{
namespace
{
using namespace gps_l1c_tables;

// Legendre sequence of length 10223 (+/-1): l[(i*i) % 10223] = -1, l[0] = 1, rest +1. Built once.
const std::vector<int>& legendre()
{
    static const std::vector<int> L = [] {
        std::vector<int> l( 10223, 1 );
        for( int i = 0; i < 10224; ++i )
        {
            l[( i * i ) % 10223] = -1;
        }
        l[0] = 1;
        return l;
    }();
    return L;
}

// One 10230-chip Weil primary (+/-1): weilcode[i] = -L[i]*L[(i+w) mod 10223], with a fixed 7-chip
// sequence inserted at index (insert-1). Shared by L1Cp/L1Cd (different weil/insert tables).
std::vector<float> weil_primary( int prn, const int* weil_tbl, const int* insert_tbl )
{
    const int          w   = weil_tbl[prn - 1];
    const int          p   = insert_tbl[prn - 1] - 1;
    const std::vector<int>& leg = legendre();

    std::vector<float> weilcode( 10223 );
    for( int i = 0; i < 10223; ++i )
    {
        weilcode[i] = -static_cast<float>( leg[i] ) * static_cast<float>( leg[( i + w ) % 10223] );
    }

    static const int   insertbit[7] = { -1, 1, 1, -1, 1, -1, -1 };
    std::vector<float> code( 10230 );
    int                i = 0;
    for( ; i < p; ++i )
    {
        code[i] = weilcode[i];
    }
    for( int j = 0; j < 7; ++j )
    {
        code[i++] = static_cast<float>( insertbit[j] );
    }
    for( ; i < 10230; ++i )
    {
        code[i] = weilcode[i - 7];
    }
    return code;
}

// BOC(1,1): chip c -> half-chip pair (-c, +c). Same convention as galileo_e1_code.
std::vector<float> apply_boc11( const std::vector<float>& primary )
{
    std::vector<float> boc( primary.size() * 2 );
    for( size_t i = 0; i < primary.size(); ++i )
    {
        boc[2 * i]     = -primary[i];
        boc[2 * i + 1] = primary[i];
    }
    return boc;
}

// Resample a BOC half-chip sequence to one period at sample_rate_hz (floor of accumulated chip phase).
Complex_buf resample_boc( const std::vector<float>& boc, double period_s, uint32_t sample_rate_hz )
{
    const int    n   = static_cast<int>( std::round( period_s * sample_rate_hz ) );
    const int    len = static_cast<int>( boc.size() );
    Complex_buf  out( n );
    double       coff = 0.0;
    const double ci   = static_cast<double>( len ) / n;
    for( int i = 0; i < n; ++i, coff += ci )
    {
        if( coff >= len )
        {
            coff -= len;
        }
        out[i] = Complex_sample( boc[static_cast<int>( coff )], 0.0f );
    }
    return out;
}

// Octal string -> 11-element +/-1 array (GNSS-SDRLIB oct2bin): each octal digit -> 3 bits MSB-first
// with logical 0 -> +1, 1 -> -1; `skiplast` drops the spare bits at the end vs the start; `flip`
// reverses. n=4 octal digits -> 12 bits, nbit=11 so one bit is skipped.
std::array<int, 11> oct2bin( const char* oct, bool skiplast, bool flip )
{
    static const int octlist[8][3] = {
        { 1, 1, 1 }, { 1, 1, -1 }, { 1, -1, 1 }, { 1, -1, -1 }, { -1, 1, 1 }, { -1, 1, -1 }, { -1, -1, 1 }, { -1, -1, -1 }
    };
    constexpr int       n = 4, nbit = 11, skip = 3 * n - nbit; // skip = 1
    std::array<int, 11> bin {};
    int                 j = 0;
    for( int i = 0; i < n; ++i )
    {
        for( int k = 0; k < 3; ++k )
        {
            if( !skiplast && i == 0 && k < skip )
            {
                continue;
            }
            if( skiplast && i == n - 1 && k >= 3 - skip )
            {
                continue;
            }
            bin[j++] = octlist[oct[i] - '0'][k];
        }
    }
    if( flip )
    {
        std::reverse( bin.begin(), bin.end() );
    }
    return bin;
}

// L1Co overlay (1800 chips, +/-1) for a GPS PRN (<64): the S1 Fibonacci LFSR (multiplicative +/-1 form).
std::vector<float> overlay_l1co( int prn )
{
    std::array<int, 11> r1 = oct2bin( L1CO_S1INIT[prn - 1], /*skiplast=*/false, /*flip=*/true );
    std::array<int, 11> t1 = oct2bin( L1CO_S1POLY[prn - 1], /*skiplast=*/true, /*flip=*/true );
    t1[10]                 = -1; // last tap always set

    std::vector<float> code( 1800 );
    for( int i = 0; i < 1800; ++i )
    {
        const int s1 = r1[10];
        int       c1 = 1;
        for( int j = 0; j < 11; ++j )
        {
            if( t1[j] == -1 )
            {
                c1 *= r1[j];
            }
        }
        for( int j = 10; j > 0; --j )
        {
            r1[j] = r1[j - 1];
        }
        r1[0]   = c1;
        code[i] = static_cast<float>( -s1 );
    }
    return code;
}

void check_prn( Satellite_id prn, int max_prn, const char* which )
{
    if( prn < 1 || prn > static_cast<Satellite_id>( max_prn ) )
    {
        throw std::out_of_range( std::string( "GPS " ) + which + " PRN out of range (1.." + std::to_string( max_prn ) + ")" );
    }
}
} // namespace

std::vector<float> L1cp_code::primary_chips( Satellite_id prn_id )
{
    check_prn( prn_id, MAX_PRN, "L1Cp" );
    return weil_primary( static_cast<int>( prn_id ), L1CP_WEIL, L1CP_INSERT );
}
std::vector<float> L1cp_code::chips( Satellite_id prn_id )
{
    return apply_boc11( primary_chips( prn_id ) );
}
std::vector<float> L1cp_code::overlay_chips( Satellite_id prn_id )
{
    check_prn( prn_id, MAX_PRN, "L1Cp" );
    return overlay_l1co( static_cast<int>( prn_id ) );
}
Complex_buf L1cp_code::generate( Satellite_id prn_id, uint32_t sampling_rate )
{
    return resample_boc( chips( prn_id ), PERIOD_TIME, sampling_rate );
}

std::vector<float> L1cd_code::primary_chips( Satellite_id prn_id )
{
    check_prn( prn_id, MAX_PRN, "L1Cd" );
    return weil_primary( static_cast<int>( prn_id ), L1CD_WEIL, L1CD_INSERT );
}
std::vector<float> L1cd_code::chips( Satellite_id prn_id )
{
    return apply_boc11( primary_chips( prn_id ) );
}
Complex_buf L1cd_code::generate( Satellite_id prn_id, uint32_t sampling_rate )
{
    return resample_boc( chips( prn_id ), PERIOD_TIME, sampling_rate );
}

} // namespace gps

#ifdef ENABLE_UNIT_TESTS
#include <numeric>
#include <catch2/catch_test_macros.hpp>

namespace
{
// Reference first-16 chips, computed by an INDEPENDENT Python reimplementation of the IS-GPS-800
// Legendre/Weil/insert (and the S1 LFSR) reading the SAME GNSS-SDRLIB tables. Matching these pins the
// hand-ported algorithm (not the verbatim-extracted tables).
constexpr float L1CP_PRN1[16] = { -1, -1, -1, 1, -1, 1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1 };
constexpr float L1CD_PRN1[16] = { 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1 };
constexpr float L1CO_PRN1[16] = { 1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, -1, -1, -1, -1 };
constexpr float L1CP_PRN4[16] = { 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1 };
constexpr float L1CD_PRN4[16] = { -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

void require_pm1( const std::vector<float>& c )
{
    for( float v : c )
    {
        REQUIRE( ( v == 1.0f || v == -1.0f ) );
    }
}
} // namespace

TEST_CASE( "gps_l1c_primary_codes_match_reference", "[gps][l1c][code]" )
{
    const auto cp1 = gps::L1cp_code::primary_chips( 1 );
    const auto cd1 = gps::L1cd_code::primary_chips( 1 );
    const auto cp4 = gps::L1cp_code::primary_chips( 4 );
    const auto cd4 = gps::L1cd_code::primary_chips( 4 );
    REQUIRE( cp1.size() == 10230 );
    require_pm1( cp1 );
    require_pm1( cd1 );
    for( int i = 0; i < 16; ++i )
    {
        INFO( "chip " << i );
        REQUIRE( cp1[i] == L1CP_PRN1[i] );
        REQUIRE( cd1[i] == L1CD_PRN1[i] );
        REQUIRE( cp4[i] == L1CP_PRN4[i] );
        REQUIRE( cd4[i] == L1CD_PRN4[i] );
    }
}

TEST_CASE( "gps_l1c_primary_codes_balanced", "[gps][l1c][code]" )
{
    // The 10230-chip L1C primaries are exactly balanced (5115 of each sign -> sum 0).
    for( Satellite_id prn : { 1, 4, 11, 20, 31 } )
    {
        const auto cp = gps::L1cp_code::primary_chips( prn );
        const auto cd = gps::L1cd_code::primary_chips( prn );
        REQUIRE( std::accumulate( cp.begin(), cp.end(), 0.0f ) == 0.0f );
        REQUIRE( std::accumulate( cd.begin(), cd.end(), 0.0f ) == 0.0f );
    }
}

TEST_CASE( "gps_l1c_overlay_code", "[gps][l1c][code]" )
{
    const auto ov1 = gps::L1cp_code::overlay_chips( 1 );
    REQUIRE( ov1.size() == 1800 );
    require_pm1( ov1 );
    REQUIRE( std::accumulate( ov1.begin(), ov1.end(), 0.0f ) == 0.0f ); // balanced (900/900)
    for( int i = 0; i < 16; ++i )
    {
        INFO( "overlay chip " << i );
        REQUIRE( ov1[i] == L1CO_PRN1[i] );
    }
}

TEST_CASE( "gps_l1c_boc11_structure", "[gps][l1c][code]" )
{
    const auto prim = gps::L1cp_code::primary_chips( 4 );
    const auto boc  = gps::L1cp_code::chips( 4 );
    REQUIRE( boc.size() == 2 * prim.size() ); // 20460 half-chips
    for( size_t i = 0; i < prim.size(); ++i ) // each chip c -> (-c, +c)
    {
        REQUIRE( boc[2 * i] == -prim[i] );
        REQUIRE( boc[2 * i + 1] == prim[i] );
    }
    // generate() samples one 10 ms period at the requested rate.
    REQUIRE( gps::L1cp_code::generate( 4, 2048000 ).size() == 20480 );
}

TEST_CASE( "gps_l1c_prn_range_checked", "[gps][l1c][code]" )
{
    REQUIRE_THROWS_AS( gps::L1cp_code::primary_chips( 0 ), std::out_of_range );
    REQUIRE_THROWS_AS( gps::L1cp_code::primary_chips( 64 ), std::out_of_range );
    REQUIRE_NOTHROW( gps::L1cp_code::primary_chips( 63 ) );
}
#endif


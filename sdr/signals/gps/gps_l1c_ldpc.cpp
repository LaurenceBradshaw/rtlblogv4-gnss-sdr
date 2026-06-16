#include "gps_l1c_ldpc.h"
#include <cmath>
#include <vector>
#include "gps_l1c_ldpc_tables.h"

namespace gps::l1c::ldpc
{
namespace
{
constexpr int    MAX_ITER = 250;
constexpr double ERR_PROB = 1e-5;                            // assumed hard-symbol error probability
const double     LCH      = std::log( ( 1.0 - ERR_PROB ) / ERR_PROB ); // |channel LLR| for a hard symbol

// Sparse parity-check matrix as bipartite adjacency: which variable nodes each check touches, and which
// checks each variable touches, plus a flat edge list. Built once per code from the IS-GPS-800E tables.
struct Ldpc_code
{
    int                           m = 0, n = 0; // checks, codeword length
    std::vector<std::vector<int>> chk_edges;    // edge ids per check
    std::vector<std::vector<int>> var_edges;    // edge ids per variable
    std::vector<int>              e_var;         // variable node of each edge
    std::vector<int>              e_chk;         // check node of each edge
    bool                          built = false;
};

// Add the 1-positions of a sub-table at (row_off + r-1, col_off + c-1). Mirrors gen_B_LDPC_H.
template <size_t N>
void add_block(
    std::vector<std::pair<int, int>>& ones, const uint16_t ( &tbl )[N][2], int row_off, int col_off
)
{
    for( size_t i = 0; i < N; i++ )
    {
        ones.emplace_back( row_off + tbl[i][0] - 1, col_off + tbl[i][1] - 1 );
    }
}

void build( Ldpc_code& c, int m, int n, int g, const std::vector<std::pair<int, int>>& ones )
{
    c.m = m;
    c.n = n;
    c.chk_edges.assign( m, {} );
    c.var_edges.assign( n, {} );
    for( const auto& [r, col] : ones )
    {
        const int e = static_cast<int>( c.e_chk.size() );
        c.e_chk.push_back( r );
        c.e_var.push_back( col );
        c.chk_edges[r].push_back( e );
        c.var_edges[col].push_back( e );
    }
    c.built = true;
    (void) g;
}

// IS-GPS-800E LDPC structure (gen_B_LDPC_H): A at (0,0); B at (0,m); C at (m-g,0); D at (m-g,m);
// E at (m-g,m+g); T at (0,m+g). g is the gap (1 for CNAV-2).
template <size_t NA, size_t NB, size_t NC, size_t ND, size_t NE, size_t NT>
Ldpc_code make_code(
    int m, int n, int g, const uint16_t ( &A )[NA][2], const uint16_t ( &B )[NB][2],
    const uint16_t ( &C )[NC][2], const uint16_t ( &D )[ND][2], const uint16_t ( &E )[NE][2],
    const uint16_t ( &T )[NT][2]
)
{
    std::vector<std::pair<int, int>> ones;
    ones.reserve( NA + NB + NC + ND + NE + NT );
    add_block( ones, A, 0, 0 );
    add_block( ones, B, 0, m );
    add_block( ones, C, m - g, 0 );
    add_block( ones, D, m - g, m );
    add_block( ones, E, m - g, m + g );
    add_block( ones, T, 0, m + g );
    Ldpc_code c;
    build( c, m, n, g, ones );
    return c;
}

const Ldpc_code& sf2_code()
{
    static const Ldpc_code c = make_code(
        600, 1200, 1, H_CNV2_SF2_A, H_CNV2_SF2_B, H_CNV2_SF2_C, H_CNV2_SF2_D, H_CNV2_SF2_E, H_CNV2_SF2_T
    );
    return c;
}
const Ldpc_code& sf3_code()
{
    static const Ldpc_code c = make_code(
        274, 548, 1, H_CNV2_SF3_A, H_CNV2_SF3_B, H_CNV2_SF3_C, H_CNV2_SF3_D, H_CNV2_SF3_E, H_CNV2_SF3_T
    );
    return c;
}

// Log-domain sum-product (belief-propagation) decode. LLR sign convention: positive => bit 0. Returns
// true iff all parity checks pass; fills `msg` with the first m systematic bits.
bool decode( const Ldpc_code& c, const uint8_t* syms, uint8_t* msg )
{
    const int           E = static_cast<int>( c.e_chk.size() );
    std::vector<double> Lvc( E ), Lcv( E, 0.0 ); // variable->check and check->variable messages

    // Resolve the global polarity ambiguity (BPSK 180-deg lock + the demod's sign convention): the hard
    // symbols are either the codeword or its bitwise complement. The complement is generally NOT a
    // codeword (the code has odd-degree checks), so pick the polarity whose hard syndrome is smaller -
    // that is the true codeword. The decoded message then comes out in the correct polarity directly.
    auto hard_syndrome = [&]( int inv )
    {
        int u = 0;
        for( int ch = 0; ch < c.m; ch++ )
        {
            int p = 0;
            for( int e : c.chk_edges[ch] )
            {
                p ^= ( syms[c.e_var[e]] ^ inv ) & 1;
            }
            u += p;
        }
        return u;
    };
    const int inv = ( hard_syndrome( 1 ) < hard_syndrome( 0 ) ) ? 1 : 0;

    std::vector<double> Lch( c.n );
    for( int v = 0; v < c.n; v++ )
    {
        Lch[v] = ( ( syms[v] ^ inv ) == 0 ) ? LCH : -LCH;
    }
    for( int e = 0; e < E; e++ )
    {
        Lvc[e] = Lch[c.e_var[e]]; // init: variable sends its channel LLR
    }

    std::vector<uint8_t> bit( c.n, 0 );
    for( int iter = 0; iter < MAX_ITER; iter++ )
    {
        // Check-node update (tanh rule), numerically stable via prefix/suffix products of tanh(L/2).
        for( int ch = 0; ch < c.m; ch++ )
        {
            const auto& es = c.chk_edges[ch];
            const int   k  = static_cast<int>( es.size() );
            if( k == 0 )
            {
                continue;
            }
            std::vector<double> t( k );
            for( int j = 0; j < k; j++ )
            {
                double x = std::tanh( 0.5 * Lvc[es[j]] );
                t[j]     = std::min( std::max( x, -0.999999999 ), 0.999999999 ); // clamp for atanh
            }
            std::vector<double> pre( k + 1, 1.0 ), suf( k + 1, 1.0 );
            for( int j = 0; j < k; j++ )
            {
                pre[j + 1] = pre[j] * t[j];
            }
            for( int j = k - 1; j >= 0; j-- )
            {
                suf[j] = suf[j + 1] * t[j];
            }
            for( int j = 0; j < k; j++ )
            {
                const double prod_others = pre[j] * suf[j + 1]; // product over all edges except j
                Lcv[es[j]]               = 2.0 * std::atanh( std::min( std::max( prod_others, -0.999999999 ), 0.999999999 ) );
            }
        }

        // Variable-node update + hard decision.
        for( int v = 0; v < c.n; v++ )
        {
            double sum = Lch[v];
            for( int e : c.var_edges[v] )
            {
                sum += Lcv[e];
            }
            for( int e : c.var_edges[v] )
            {
                Lvc[e] = sum - Lcv[e];
            }
            bit[v] = ( sum < 0.0 ) ? 1 : 0;
        }

        // Parity check: stop as soon as every check is satisfied.
        bool ok = true;
        for( int ch = 0; ch < c.m && ok; ch++ )
        {
            int parity = 0;
            for( int e : c.chk_edges[ch] )
            {
                parity ^= bit[c.e_var[e]];
            }
            if( parity )
            {
                ok = false;
            }
        }
        if( ok )
        {
            for( int i = 0; i < c.m; i++ )
            {
                msg[i] = bit[i];
            }
            return true;
        }
    }

    for( int i = 0; i < c.m; i++ ) // not converged: still return the best hard decision
    {
        msg[i] = bit[i];
    }
    return false;
}
} // namespace

bool decode_sf2( const uint8_t* syms, uint8_t* msg )
{
    return decode( sf2_code(), syms, msg );
}
bool decode_sf3( const uint8_t* syms, uint8_t* msg )
{
    return decode( sf3_code(), syms, msg );
}

} // namespace gps::l1c::ldpc

#ifdef ENABLE_UNIT_TESTS
#include <array>
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "gps_l1c_ldpc_all_zero_codeword", "[gps][l1c][ldpc]" )
{
    // The all-zero word is a valid codeword of any linear code (H*0 = 0): must decode + report valid.
    std::array<uint8_t, 1200> z2 {};
    std::array<uint8_t, 600>  m2 {};
    REQUIRE( gps::l1c::ldpc::decode_sf2( z2.data(), m2.data() ) );
    std::array<uint8_t, 548> z3 {};
    std::array<uint8_t, 274> m3 {};
    REQUIRE( gps::l1c::ldpc::decode_sf3( z3.data(), m3.data() ) );
}

TEST_CASE( "gps_l1c_ldpc_corrects_errors_sf2", "[gps][l1c][ldpc]" )
{
    // All-zero codeword with N flipped symbols -> belief propagation must correct back to all-zero.
    for( int nerr : { 1, 5, 10, 20 } )
    {
        std::array<uint8_t, 1200> s {};
        for( int i = 0; i < nerr; i++ )
        {
            s[( i * 97 + 13 ) % 1200] = 1; // spread-out deterministic error positions
        }
        std::array<uint8_t, 600> m {};
        INFO( "SF2 nerr=" << nerr );
        REQUIRE( gps::l1c::ldpc::decode_sf2( s.data(), m.data() ) );
        for( auto b : m )
        {
            REQUIRE( b == 0 );
        }
    }
}

TEST_CASE( "gps_l1c_ldpc_corrects_errors_sf3", "[gps][l1c][ldpc]" )
{
    for( int nerr : { 1, 5, 10 } )
    {
        std::array<uint8_t, 548> s {};
        for( int i = 0; i < nerr; i++ )
        {
            s[( i * 53 + 7 ) % 548] = 1;
        }
        std::array<uint8_t, 274> m {};
        INFO( "SF3 nerr=" << nerr );
        REQUIRE( gps::l1c::ldpc::decode_sf3( s.data(), m.data() ) );
        for( auto b : m )
        {
            REQUIRE( b == 0 );
        }
    }
}
#endif

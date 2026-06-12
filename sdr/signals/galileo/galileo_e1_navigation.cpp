#include "galileo_e1_navigation.h"

#include <algorithm>
#include <cstdlib>
#include "logging.h"

namespace
{
// Soft-decision Viterbi, rate 1/2, K=7
// Ports GNSS-SDR Viterbi_Decoder for the Galileo FEC: G1=0171, G2=0133 (octal), the
// non-systematic convolutional encoder, terminated with MM=6 tail bits. Branch metric is
// the soft correlation (Gamma): an output bit of 1 adds the matching received soft value.
constexpr int   VIT_K      = 7;
constexpr int   VIT_MM     = 6;
constexpr int   VIT_STATES = 1 << VIT_MM; // 64
constexpr int   VIT_G1     = 0171;        // octal 171
constexpr int   VIT_G2     = 0133;        // octal 133
constexpr float VIT_MAXLOG = 1e7F;

int parity( int x )
{
    int p = 0;
    while( x )
    {
        p ^= ( x & 1 );
        x >>= 1;
    }
    return p;
}

// Precomputed encoder transition tables (one set, shared by all decoders).
struct Viterbi_tables
{
    int out0[VIT_STATES], out1[VIT_STATES], state0[VIT_STATES], state1[VIT_STATES];
    Viterbi_tables()
    {
        for( int s = 0; s < VIT_STATES; s++ )
        {
            out0[s] = enc( 0, s, &state0[s] );
            out1[s] = enc( 1, s, &state1[s] );
        }
    }
    static int enc( int input, int state_in, int* next )
    {
        const int reg = ( input << ( VIT_K - 1 ) ) | state_in; // 7-bit register
        const int o   = ( parity( reg & VIT_G1 ) << 1 ) | parity( reg & VIT_G2 );
        *next         = reg >> 1;
        return o;
    }
};

// Decode 2*(LL+MM) soft symbols (G1,G2 interleaved, G2 already un-inverted) -> LL bits.
void viterbi_decode( const std::vector<float>& in, std::vector<int>& out )
{
    static const Viterbi_tables T;
    const int                   LL = static_cast<int>( out.size() );
    const int                   N  = LL + VIT_MM;

    std::vector<float> prev( VIT_STATES, -VIT_MAXLOG ), next( VIT_STATES, -VIT_MAXLOG );
    std::vector<int>   pbit( static_cast<size_t>( N ) * VIT_STATES );
    std::vector<int>   pstate( static_cast<size_t>( N ) * VIT_STATES );
    prev[0] = 0.0F; // start in the all-zeros state

    for( int t = 0; t < N; t++ )
    {
        const float r0 = in[2 * t];     // G1 received soft
        const float r1 = in[2 * t + 1]; // G2 received soft
        // Branch metric for each candidate output symbol (bit1=G1, bit0=G2).
        const float mc[4] = { 0.0F, r1, r0, r0 + r1 };

        for( int s = 0; s < VIT_STATES; s++ )
        {
            const float m0 = prev[s] + mc[T.out0[s]];
            if( m0 > next[T.state0[s]] )
            {
                next[T.state0[s]]                                           = m0;
                pstate[static_cast<size_t>( t ) * VIT_STATES + T.state0[s]] = s;
                pbit[static_cast<size_t>( t ) * VIT_STATES + T.state0[s]]   = 0;
            }
            const float m1 = prev[s] + mc[T.out1[s]];
            if( m1 > next[T.state1[s]] )
            {
                next[T.state1[s]]                                           = m1;
                pstate[static_cast<size_t>( t ) * VIT_STATES + T.state1[s]] = s;
                pbit[static_cast<size_t>( t ) * VIT_STATES + T.state1[s]]   = 1;
            }
        }

        const float mx = *std::max_element( next.begin(), next.end() );
        for( int s = 0; s < VIT_STATES; s++ )
        {
            prev[s] = next[s] - mx;
            next[s] = -VIT_MAXLOG;
        }
    }

    // Trace back: tail steps (no output), then the LL data bits.
    int s = 0;
    for( int t = N - 1; t >= LL; t-- )
    {
        s = pstate[static_cast<size_t>( t ) * VIT_STATES + s];
    }
    for( int t = LL - 1; t >= 0; t-- )
    {
        out[t] = pbit[static_cast<size_t>( t ) * VIT_STATES + s];
        s      = pstate[static_cast<size_t>( t ) * VIT_STATES + s];
    }
}

// CRC-24Q (Galileo I/NAV): poly 0x1864CFB, init 0, MSB-first, no reflection/final-xor. Run
// over n message bits (leading zeros pad to a byte boundary but don't affect a zero-init CRC).
uint32_t crc24q( const int* bits, int n )
{
    uint32_t crc = 0;
    for( int i = 0; i < n; i++ )
    {
        const uint32_t top = ( crc >> 23 ) & 1u;
        crc                = ( crc << 1 ) & 0xFFFFFFu;
        if( top ^ static_cast<uint32_t>( bits[i] & 1 ) )
        {
            crc ^= 0x864CFBu;
        }
    }
    return crc;
}

uint32_t bits_to_uint( const int* bits, int start, int len )
{
    uint32_t v = 0;
    for( int i = 0; i < len; i++ )
    {
        v = ( v << 1 ) | static_cast<uint32_t>( bits[start + i] & 1 );
    }
    return v;
}

// I/NAV field readers. The Galileo SIS-ICD bit tables are 1-based from the start of the
// 128-bit Data_jk word; jk[] here is 0-based MSB-first, so position p maps to index p-1.
uint32_t read_u( const int* jk, int start, int len )
{
    return bits_to_uint( jk, start - 1, len );
}
int64_t read_s( const int* jk, int start, int len )
{
    const uint32_t u = read_u( jk, start, len );
    if( ( u >> ( len - 1 ) ) & 1u ) // sign bit set -> two's complement
    {
        return static_cast<int64_t>( u ) - ( static_cast<int64_t>( 1 ) << len );
    }
    return static_cast<int64_t>( u );
}

// LSB scale factors (Galileo OS SIS-ICD), from gnss-sdr MATH_CONSTANTS.h.
constexpr double PI_TWO_N31 = 1.462918079267160e-009; // pi * 2^-31
constexpr double PI_TWO_N43 = 3.571577341960839e-013; // pi * 2^-43
constexpr double TWO_N5     = 0.03125;                // 2^-5
constexpr double TWO_N19    = 1.907348632812500e-006; // 2^-19
constexpr double TWO_N29    = 1.862645149230957e-009; // 2^-29
constexpr double TWO_N33    = 1.164153218269348e-010; // 2^-33
constexpr double TWO_N34    = 5.82076609134674e-011;  // 2^-34
constexpr double TWO_N46    = 1.4210854715202e-014;   // 2^-46
constexpr double TWO_N59    = 1.73472347597681e-018;  // 2^-59
} // namespace

Galileo_e1b_decoder::Galileo_e1b_decoder( Satellite_id prn, const Signal_params& params )
    : satellite_id_( prn ),
      params_( params )
{
    // Stamp constellation + PRN so the observables/PVT pipeline dispatches the Galileo orbit/clock.
    eph_.constellation = Constellation::Galileo;
    eph_.prn           = prn;
}

void Galileo_e1b_decoder::process( double prompt_i, double /*prompt_i_prev*/ )
{
    symbol_count_++;
    const float sym  = static_cast<float>( prompt_i );
    const int   sign = ( sym >= 0.0 ) ? 1 : -1;

    // Slide the sign into the preamble window.
    for( int i = 0; i < PREAMBLE_LEN - 1; i++ )
    {
        pre_win_[i] = pre_win_[i + 1];
    }
    pre_win_[PREAMBLE_LEN - 1] = sign;
    if( pre_fill_ < PREAMBLE_LEN )
    {
        pre_fill_++;
    }

    if( !synced_ )
    {
        if( pre_fill_ < PREAMBLE_LEN )
        {
            return;
        }
        int corr = 0;
        for( int i = 0; i < PREAMBLE_LEN; i++ )
        {
            corr += pre_win_[i] * PREAMBLE[i];
        }
        if( std::abs( corr ) >= PREAMBLE_LEN - 1 ) // allow 1 symbol error
        {
            // Vote by page-part phase: the true preamble recurs at the SAME phase (mod 250)
            // every page part, so its phase accumulates votes while random false matches
            // (~1% of windows) scatter across phases. Sync when one phase reaches the vote
            // threshold. Robust where the old "exactly 250 from the last detection" broke - a
            // single false match between two true preambles would reset that and block sync.
            const int phase = static_cast<int>( symbol_count_ % PAGE_PART_SYMS );
            pre_pol_[phase] += ( corr >= 0 ) ? 1 : -1;
            if( ++pre_votes_[phase] >= PREAMBLE_VOTE_TH )
            {
                synced_         = true;
                frame_polarity_ = ( pre_pol_[phase] >= 0 ) ? 1 : -1;
                frame_pos_      = PREAMBLE_LEN - 1; // just consumed the 10th preamble symbol
                data_buf_.clear();
                have_even_ = false;
                logging::log(
                    logging::Level::Info,
                    fmt::format(
                        "Navigation - {} PRN {:2d}  preamble found  polarity={:+d}",
                        params_.name,
                        satellite_id_,
                        frame_polarity_
                    )
                );
            }
        }
        return;
    }

    // Synced: collect the 240 data symbols of each page part (skip the 10 preamble symbols).
    frame_pos_ = ( frame_pos_ + 1 ) % PAGE_PART_SYMS;
    if( frame_pos_ >= PREAMBLE_LEN )
    {
        data_buf_.push_back( frame_polarity_ * sym );
        if( static_cast<int>( data_buf_.size() ) == DATA_SYMS )
        {
            decode_page_part();
            data_buf_.clear();
        }
    }
}

void Galileo_e1b_decoder::decode_page_part()
{
    // 1. De-interleave (8 rows x 30 cols): out[c*8+r] = in[r*30+c].
    std::vector<float> deint( DATA_SYMS );
    for( int r = 0; r < 8; r++ )
    {
        for( int c = 0; c < 30; c++ )
        {
            deint[c * 8 + r] = data_buf_[r * 30 + c];
        }
    }
    // 2. Un-invert the G2 output (NOT gate in the encoder): negate every odd-index symbol.
    for( int i = 0; i < DATA_SYMS; i++ )
    {
        if( i & 1 )
        {
            deint[i] = -deint[i];
        }
    }
    // 3. Viterbi -> 114 bits.
    std::vector<int> bits( PAGE_PART_BITS );
    viterbi_decode( deint, bits );

    // 4. Even/odd page assembly (bit 0 is the even/odd flag: 0=even, 1=odd).
    if( bits[0] == 0 )
    {
        page_even_ = bits; // store the even page part
        have_even_ = true;
        return;
    }
    if( !have_even_ )
    {
        return; // odd without a preceding even
    }

    // Assemble even(114) + odd(114) = 228 bits; CRC-24Q over bits [0,196), checksum at [196,220).
    std::vector<int> page( 2 * PAGE_PART_BITS );
    std::copy( page_even_.begin(), page_even_.end(), page.begin() );
    std::copy( bits.begin(), bits.end(), page.begin() + PAGE_PART_BITS );
    have_even_ = false;

    const uint32_t computed = crc24q( page.data(), 196 );
    const uint32_t checksum = bits_to_uint( page.data(), 196, 24 );
    if( computed == checksum )
    {
        crc_ok_++;

        // Rebuild the 128-bit Data_jk word from the CRC-valid page and extract its fields.
        // Layout (gnss-sdr split_page): Data_k = even bits [2,114), Data_j = odd bits [116,132).
        int jk[128];
        for( int i = 0; i < 112; i++ )
        {
            jk[i] = page[2 + i];
        }
        for( int i = 0; i < 16; i++ )
        {
            jk[112 + i] = page[PAGE_PART_BITS + 2 + i];
        }
        extract_word( jk );
    }
    else
    {
        crc_fail_++;
    }
}

// Extract one I/NAV word (128-bit Data_jk) into eph_. Word types 1-4 carry the orbit/clock
// parameters (gated on a common IODnav); word 5 carries the GST week/TOW. Mirrors gnss-sdr
// galileo_inav_message.cc read_page_1..5; bit positions/scales from Galileo_INAV.h.
void Galileo_e1b_decoder::extract_word( const int* jk )
{
    const int word_type = static_cast<int>( read_u( jk, 1, 6 ) );

    // Words 1-4 share an IODnav; a change starts a fresh ephemeris set (word 5 is independent).
    auto check_iod = [&]( int iod )
    {
        if( iod != eph_iodnav_ )
        {
            eph_iodnav_ = iod;
            eph_words_ &= ~0xF; // drop words 1-4, keep GST (bit 4)
            eph_.valid = false;
        }
        eph_.iod_nav = iod;
    };

    switch( word_type )
    {
    case 1: // Ephemeris (1/4)
        check_iod( static_cast<int>( read_u( jk, 7, 10 ) ) );
        eph_.toe    = static_cast<double>( read_u( jk, 17, 14 ) ) * 60.0;
        eph_.m0     = static_cast<double>( read_s( jk, 31, 32 ) ) * PI_TWO_N31;
        eph_.e      = static_cast<double>( read_u( jk, 63, 32 ) ) * TWO_N33;
        eph_.sqrt_a = static_cast<double>( read_u( jk, 95, 32 ) ) * TWO_N19;
        eph_words_ |= 0x1;
        break;
    case 2: // Ephemeris (2/4)
        check_iod( static_cast<int>( read_u( jk, 7, 10 ) ) );
        eph_.omega0 = static_cast<double>( read_s( jk, 17, 32 ) ) * PI_TWO_N31;
        eph_.i0     = static_cast<double>( read_s( jk, 49, 32 ) ) * PI_TWO_N31;
        eph_.omega  = static_cast<double>( read_s( jk, 81, 32 ) ) * PI_TWO_N31;
        eph_.idot   = static_cast<double>( read_s( jk, 113, 14 ) ) * PI_TWO_N43;
        eph_words_ |= 0x2;
        break;
    case 3: // Ephemeris (3/4) + SISA
        check_iod( static_cast<int>( read_u( jk, 7, 10 ) ) );
        eph_.omegadot = static_cast<double>( read_s( jk, 17, 24 ) ) * PI_TWO_N43;
        eph_.delta_n  = static_cast<double>( read_s( jk, 41, 16 ) ) * PI_TWO_N43;
        eph_.cuc      = static_cast<double>( read_s( jk, 57, 16 ) ) * TWO_N29;
        eph_.cus      = static_cast<double>( read_s( jk, 73, 16 ) ) * TWO_N29;
        eph_.crc      = static_cast<double>( read_s( jk, 89, 16 ) ) * TWO_N5;
        eph_.crs      = static_cast<double>( read_s( jk, 105, 16 ) ) * TWO_N5;
        eph_.sisa     = static_cast<int>( read_u( jk, 121, 8 ) );
        eph_words_ |= 0x4;
        break;
    case 4: // Ephemeris (4/4) + clock correction
        check_iod( static_cast<int>( read_u( jk, 7, 10 ) ) );
        eph_.cic = static_cast<double>( read_s( jk, 23, 16 ) ) * TWO_N29;
        eph_.cis = static_cast<double>( read_s( jk, 39, 16 ) ) * TWO_N29;
        eph_.toc = static_cast<double>( read_u( jk, 55, 14 ) ) * 60.0;
        eph_.af0 = static_cast<double>( read_s( jk, 69, 31 ) ) * TWO_N34;
        eph_.af1 = static_cast<double>( read_s( jk, 100, 21 ) ) * TWO_N46;
        eph_.af2 = static_cast<double>( read_s( jk, 121, 6 ) ) * TWO_N59;
        eph_words_ |= 0x8;
        break;
    case 5: // GST week + time of week (also iono/BGD/health, not extracted)
        eph_.week = static_cast<int>( read_u( jk, 74, 12 ) );
        eph_.tow  = static_cast<double>( read_u( jk, 86, 20 ) );
        eph_words_ |= 0x10;
        // Anchor the SV transmit-time clock. TOW5 is the GST at the even-page-part preamble; we
        // are now at the end of the odd part, (2*PAGE_PART_SYMS - 1) = 499 symbols later (see
        // ms_since_tow_update). The page CRC-24 already validated this, so trust it directly -
        // no GPS-style TOW-continuity gate needed.
        tow_ref_symbol_ = symbol_count_ - ( 2 * PAGE_PART_SYMS - 1 );
        tow_anchored_   = true;
        tow_confirmed_  = true;
        break;
    default:
        break;
    }

    // Per-page log, mirroring the GPS per-subframe line: one line for every CRC-valid I/NAV page,
    // showing the word type, the GST tow/week (known once word 5 is in), and the ephemeris-word
    // collection mask (0x1F = complete) so a stalled SV's missing word type is visible.
    const std::string gst = ( eph_words_ & 0x10 )
                                ? fmt::format( "tow={:.1f}s  week={:d}", eph_.tow, eph_.week )
                                : std::string( "tow=?  week=?" );
    logging::log(
        logging::Level::Info,
        fmt::format(
            "Navigation - {} PRN {:2d}  word {:2d}  {}  eph=0x{:02X}", params_.name, satellite_id_, word_type, gst, eph_words_
        )
    );

    // Words 1-4 (0xF) + GST (0x10) all collected -> a complete ephemeris.
    if( !eph_.valid && ( eph_words_ & 0x1F ) == 0x1F )
    {
        eph_.valid = true;
        logging::log(
            logging::Level::Info,
            fmt::format(
                "Navigation - {} PRN {:2d}  ephemeris  tow={:.1f}s  week={:d}", params_.name, satellite_id_, eph_.tow, eph_.week
            )
        );
    }
}

#ifdef ENABLE_UNIT_TESTS
#include <numeric>
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "galileo_crc24q_zero_message", "[galileo][nav][crc]" )
{
    // Zero-init, no final XOR: an all-zeros message has CRC 0.
    std::vector<int> zeros( 200, 0 );
    REQUIRE( crc24q( zeros.data(), 200 ) == 0u );
}

TEST_CASE( "galileo_crc24q_remainder_property", "[galileo][nav][crc]" )
{
    // The defining property of a Q-CRC (zero init, MSB-first, no final XOR): appending the
    // 24-bit remainder to the message makes the CRC of the whole thing zero. This is exactly
    // how the receiver validates a page, and it is independent of the internal shift logic.
    std::vector<int> msg;
    uint32_t         lfsr = 0xACE1u;
    for( int i = 0; i < 196; ++i )
    {
        lfsr = ( lfsr >> 1 ) ^ ( ( lfsr & 1u ) ? 0xB400u : 0u ); // arbitrary repeatable bits
        msg.push_back( static_cast<int>( lfsr & 1u ) );
    }

    const uint32_t crc = crc24q( msg.data(), static_cast<int>( msg.size() ) );

    std::vector<int> with_crc = msg;
    for( int i = 23; i >= 0; --i )
    {
        with_crc.push_back( static_cast<int>( ( crc >> i ) & 1u ) );
    }
    REQUIRE( crc24q( with_crc.data(), static_cast<int>( with_crc.size() ) ) == 0u );

    // A single flipped message bit must break that (error detection).
    with_crc[50] ^= 1;
    REQUIRE( crc24q( with_crc.data(), static_cast<int>( with_crc.size() ) ) != 0u );
}

TEST_CASE( "galileo_viterbi_round_trip", "[galileo][nav][fec]" )
{
    // Encode known bits with the SAME G1=0171/G2=0133 convolutional encoder the decoder models,
    // terminate with MM zero tail bits, map each output bit to a +/-1 soft value, and confirm
    // the soft-decision Viterbi recovers the original bits exactly.
    const std::vector<int> data = { 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1 };

    std::vector<float> soft;
    int                state = 0;
    auto               emit  = [&]( int input ) {
        int       next = 0;
        const int o    = Viterbi_tables::enc( input, state, &next );
        state          = next;
        soft.push_back( ( ( o >> 1 ) & 1 ) ? 1.0F : -1.0F ); // G1
        soft.push_back( ( o & 1 ) ? 1.0F : -1.0F );          // G2
    };
    for( int b : data )
    {
        emit( b );
    }
    for( int i = 0; i < 6; ++i ) // MM=6 tail bits flush the register
    {
        emit( 0 );
    }

    std::vector<int> out( data.size() );
    viterbi_decode( soft, out );
    REQUIRE( out == data );
}

TEST_CASE( "galileo_parity_is_popcount_mod_2", "[galileo][nav][fec]" )
{
    REQUIRE( parity( 0 ) == 0 );
    REQUIRE( parity( 0b1 ) == 1 );
    REQUIRE( parity( 0b11 ) == 0 );    // two set bits
    REQUIRE( parity( 0b1011 ) == 1 );  // three set bits
    REQUIRE( parity( 0xFF ) == 0 );    // eight set bits
    REQUIRE( parity( 0xFE ) == 1 );    // seven set bits
}
#endif

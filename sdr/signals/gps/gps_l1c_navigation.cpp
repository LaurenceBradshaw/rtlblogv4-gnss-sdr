#include "gps_l1c_navigation.h"
#include <cmath>
#include <cstdint>
#include "gps_l1c_ldpc.h"
#include "logging.h"

namespace
{
// Bit-array readers: each element of `b` is a single 0/1 bit, MSB-first. (Our LDPC output is unpacked.)
uint64_t bu( const uint8_t* b, int p, int n )
{
    uint64_t v = 0;
    for( int i = 0; i < n; i++ )
    {
        v = ( v << 1 ) | ( b[p + i] & 1u );
    }
    return v;
}
int64_t bs( const uint8_t* b, int p, int n ) // two's-complement sign-extend
{
    uint64_t v = bu( b, p, n );
    if( n < 64 && ( v >> ( n - 1 ) ) )
    {
        v |= ~0ull << n;
    }
    return static_cast<int64_t>( v );
}

// CRC-24Q (poly 0x1864CFB, init 0, MSB-first) over n message bits. A frame [data | 24-bit CRC] is valid
// iff the CRC over ALL its bits is 0 (the remainder property), so no right-pad/offset logic is needed.
uint32_t crc24q( const uint8_t* bits, int n )
{
    uint32_t crc = 0;
    for( int i = 0; i < n; i++ )
    {
        const uint32_t top = ( crc >> 23 ) & 1u;
        crc                = ( crc << 1 ) & 0xFFFFFFu;
        if( top ^ ( bits[i] & 1u ) )
        {
            crc ^= 0x864CFBu;
        }
    }
    return crc;
}
} // namespace

// Subframe-1 template generator (PocketSDR sdr_code.LFSR + rev_reg + sync_CNV2_frame): the 9-bit TOI is
// coded into 52 symbols - the MSB (bit 9) is symbol 0, and the low 8 bits seed an 8-bit LFSR (polynomial
// 0b10011111) that emits 51 symbols, each XORed with the MSB. The whole 400-entry table is built once.
namespace
{
uint32_t rev_reg( uint32_t R, int N )
{
    uint32_t RR = 0;
    for( int i = 0; i < N; i++ )
    {
        RR = ( RR << 1 ) | ( ( R >> i ) & 1u );
    }
    return RR;
}
} // namespace

void Gps_l1c_decoder::build_sf1_templates()
{
    constexpr uint32_t TAP = 0b10011111u; // 8-bit LFSR polynomial (PocketSDR)
    for( int toi = 0; toi < N_TOI; toi++ )
    {
        const uint8_t bit9 = ( toi >> 8 ) & 1u;
        uint32_t      R    = rev_reg( static_cast<uint32_t>( toi ) & 0xFFu, 8 );
        sf1_[toi][0]       = bit9;
        for( int i = 0; i < SF1_SYMS - 1; i++ ) // 51 LFSR symbols
        {
            // PocketSDR maps CHIP=(1,-1): bit = (CHIP[R&1]+1)/2 = !(R&1). Then XOR the TOI MSB.
            const uint8_t b = static_cast<uint8_t>( ( R & 1u ) ^ 1u );
            sf1_[toi][1 + i] = b ^ bit9;
            R = ( static_cast<uint32_t>( __builtin_parity( R & TAP ) ) << 7 ) | ( R >> 1 );
        }
    }
}

Gps_l1c_decoder::Gps_l1c_decoder( Satellite_id prn )
    : satellite_id_( prn )
{
    // Stamp the constellation + PRN on eph_ up front (eph_ stays INVALID until a frame decodes - this is
    // only identity, not an observable) so the channel snapshot reports GPS/PRN from the first epoch. Else
    // an undecoded L1C channel publishes constellation=Unknown -> the GUI shows "?04", a grey C/N0 bar, and
    // no graph data (the per-SV history subscription keys on the constellation). Mirrors the L1CA decoder.
    eph_.constellation = Constellation::Gps;
    eph_.prn           = prn;
    build_sf1_templates();
    syms_.reserve( 4 * WINDOW );
}

// match a WINDOW-symbol view against subframe-1 of `toi` (at the start) and toi+1 (at +1800 = the next
// frame's subframe-1). +1 normal, 0 reversed (180-deg polarity), -1 no match. Tolerant (<= MATCH_TOL
// symbol mismatches per 52-block) so an occasional bit error doesn't block sync, while a random
// alignment (floor ~10-13 mismatches) is rejected.
int Gps_l1c_decoder::match_frame( const uint8_t* w, int toi ) const
{
    const auto& SF1 = sf1_[toi];
    const auto& SFn = sf1_[( toi + 1 ) % N_TOI];
    int         miss_first = 0, miss_last = 0; // mismatches assuming NORMAL polarity
    for( int i = 0; i < SF1_SYMS; i++ )
    {
        miss_first += ( w[i] != SF1[i] );
        miss_last += ( w[FRAME_SYMS + i] != SFn[i] );
    }
    if( miss_first <= MATCH_TOL && miss_last <= MATCH_TOL )
    {
        return 1; // normal polarity (almost all symbols agree)
    }
    if( miss_first >= SF1_SYMS - MATCH_TOL && miss_last >= SF1_SYMS - MATCH_TOL )
    {
        return 0; // reversed polarity (almost all symbols disagree -> whole frame inverted)
    }
    return -1;
}

// De-interleave subframes 2+3 (38x46 block) of the 1800-symbol frame at `w`, polarity-correct (XOR rev_),
// LDPC-decode SF2 (1200->600) + SF3 (548->274), and CRC-check both. (PocketSDR decode_CNV2.)
bool Gps_l1c_decoder::decode_frame( const uint8_t* w )
{
    // PocketSDR: syms_d = (frame ^ rev)[52:].reshape(46,38).T.ravel(). XOR commutes with the permute.
    std::array<uint8_t, FRAME_SYMS - SF1_SYMS> deint; // 1748 = 38*46
    for( int row = 0; row < 46; row++ )
    {
        for( int col = 0; col < 38; col++ )
        {
            deint[col * 46 + row] = w[SF1_SYMS + row * 38 + col] ^ static_cast<uint8_t>( rev_ );
        }
    }
    std::array<uint8_t, 600> sf2;
    std::array<uint8_t, 274> sf3;
    const bool ldpc2 = gps::l1c::ldpc::decode_sf2( deint.data(), sf2.data() );        // deint[0..1199]
    const bool ldpc3 = gps::l1c::ldpc::decode_sf3( deint.data() + 1200, sf3.data() ); // deint[1200..1747]
    const bool ok    = ldpc2 && ldpc3 && crc24q( sf2.data(), 600 ) == 0 && crc24q( sf3.data(), 274 ) == 0;
    logging::log(
        logging::Level::Info,
        fmt::format(
            "Navigation - GPS L1C PRN {:2d}  CNAV-2 frame TOI={:d} {} (ldpc {}/{})", satellite_id_, toi_,
            ok ? "DECODED OK" : "CRC FAIL", ldpc2 ? 1 : 0, ldpc3 ? 1 : 0
        )
    );
    if( ok )
    {
        parse_sf2( sf2.data() ); // CRC valid -> trustworthy; fill eph_ (Phase 4c)
    }
    return ok;
}

// Parse CNAV-2 subframe-2 (IS-GPS-800E Table 3.5-2) into eph_. Bit offsets are 0-based from SF2 bit 0;
// scales are the CNAV set (gnss-sdr GPS_CNAV.h): A = AREF + dA, n = n0 + dn0, OmegaDot = ref + dOmegaDot,
// angles in semicircles (x pi -> rad). The Adot / dn0dot / dn0dot rate terms are dropped for now (their
// effect over the ~2 h fit is small); they're a later refinement of the orbit model.
void Gps_l1c_decoder::parse_sf2( const uint8_t* s )
{
    constexpr double PI      = 3.14159265358979323846;
    constexpr double AREF    = 26559710.0;        // semi-major axis reference (m)
    constexpr double ODOTREF = -2.6e-9 * PI;      // OmegaDot reference (rad/s)
    auto             p2      = []( int n ) { return std::ldexp( 1.0, -n ); };

    eph_.constellation = Constellation::Gps;
    eph_.prn           = satellite_id_;
    eph_.week          = static_cast<int>( bu( s, 0, 13 ) );          // WN
    eph_.toe           = static_cast<double>( bu( s, 38, 11 ) ) * 300.0; // tOE
    eph_.toc           = eph_.toe;                                    // CNAV-2 clock ref = toe
    const double A     = AREF + static_cast<double>( bs( s, 49, 26 ) ) * p2( 9 ); // dA
    eph_.sqrt_a        = std::sqrt( A );
    eph_.delta_n       = static_cast<double>( bs( s, 100, 17 ) ) * p2( 44 ) * PI;
    eph_.m0            = static_cast<double>( bs( s, 140, 33 ) ) * p2( 32 ) * PI;
    eph_.e             = static_cast<double>( bu( s, 173, 33 ) ) * p2( 34 );
    eph_.omega         = static_cast<double>( bs( s, 206, 33 ) ) * p2( 32 ) * PI;
    eph_.omega0        = static_cast<double>( bs( s, 239, 33 ) ) * p2( 32 ) * PI;
    eph_.i0            = static_cast<double>( bs( s, 272, 33 ) ) * p2( 32 ) * PI;
    eph_.omegadot      = ODOTREF + static_cast<double>( bs( s, 305, 17 ) ) * p2( 44 ) * PI;
    eph_.idot          = static_cast<double>( bs( s, 322, 15 ) ) * p2( 44 ) * PI;
    eph_.cis           = static_cast<double>( bs( s, 337, 16 ) ) * p2( 30 );
    eph_.cic           = static_cast<double>( bs( s, 353, 16 ) ) * p2( 30 );
    eph_.crs           = static_cast<double>( bs( s, 369, 24 ) ) * p2( 8 );
    eph_.crc           = static_cast<double>( bs( s, 393, 24 ) ) * p2( 8 );
    eph_.cus           = static_cast<double>( bs( s, 417, 21 ) ) * p2( 30 );
    eph_.cuc           = static_cast<double>( bs( s, 438, 21 ) ) * p2( 30 );
    eph_.af0           = static_cast<double>( bs( s, 470, 26 ) ) * p2( 35 );
    eph_.af1           = static_cast<double>( bs( s, 496, 20 ) ) * p2( 48 );
    eph_.af2           = static_cast<double>( bs( s, 516, 10 ) ) * p2( 60 );
    eph_.group_delay   = static_cast<double>( bs( s, 526, 13 ) ) * p2( 35 );
    eph_.valid         = true;

    // TOW / transmit-time anchor. Absolute TOW = ITOW (2-hour interval of week) * 7200 + TOI * 18. Per
    // IS-GPS-800 the TOI broadcast in a frame is the count at the start of the *next* frame, so this TOW is
    // the leading edge of the NEXT frame - whose subframe 1 sits at w[1800..1851], i.e. SF1_SYMS (52) epochs
    // before the current symbol. (Anchoring it to w[0] instead, 1851 epochs back, put every SV's transmit
    // time a full frame (~18 s) late; that cancels in the pseudorange but evaluates each SV's orbit ~18 s
    // ahead - ~18000 km of arc - which left the fix ~30 km off. Verified against L1CA: the per-SV t_tx then
    // agree to a common ~tens of metres.)
    //
    // Resolve the subframe-1 TOI alias: the 52-symbol SF1 template for toi^256 (bit 9 flipped) is the exact
    // bitwise COMPLEMENT of toi's, so frame sync accepts both toi and toi^256 (at opposite carrier polarity)
    // and keeps whichever is smaller - wrong whenever the true TOI is >=256. LDPC+CRC pass either way (they
    // don't depend on the TOI), and the symbols alone can't tell the pair apart, so it's silent until the
    // TOW lands 256*18=4608 s off and the SV position (computed at t_tx) is ~18000 km out. Disambiguate with
    // toe (this same CRC-valid message): the transmit time lies within the ephemeris curve-fit interval of
    // toe, while the 4608 s alias falls far outside - so pick the candidate whose TOW is closest to toe.
    const int    itow  = static_cast<int>( bu( s, 13, 8 ) );
    const int    alias = toi_ ^ 0x100;
    const double base  = static_cast<double>( itow ) * 7200.0;
    const double tow0  = base + static_cast<double>( toi_ ) * 18.0;
    if( alias < N_TOI )
    {
        const double tow1 = base + static_cast<double>( alias ) * 18.0;
        if( std::fabs( tow1 - eph_.toe ) < std::fabs( tow0 - eph_.toe ) )
        {
            toi_ = alias; // the aliased frame number is the true one
        }
    }
    eph_.tow                     = base + static_cast<double>( toi_ ) * 18.0;
    tow_anchored_                = true;
    tow_confirmed_               = true;
    active_eph_tow_anchor_epoch_ = epoch_count_ - SF1_SYMS; // next-frame edge: 52 epochs before current
}

void Gps_l1c_decoder::process( double prompt_i, double /*prompt_i_prev*/ )
{
    // One L1Cd data symbol per epoch (10 ms code period). Hard decision; the residual polarity is
    // resolved by the frame-sync `rev` match below (PocketSDR decode_L1CD: 1 if P.real>=0 else 0).
    syms_.push_back( prompt_i >= 0.0 ? 1 : 0 );
    if( syms_.size() > 4 * static_cast<size_t>( WINDOW ) ) // amortised trim; keep well over WINDOW
    {
        syms_.erase( syms_.begin(), syms_.begin() + 2 * WINDOW );
    }
    ++epoch_count_;

    if( static_cast<int>( syms_.size() ) < WINDOW )
    {
        return;
    }
    const uint8_t* w = &syms_[syms_.size() - WINDOW];

    if( frame_synced_ )
    {
        // Re-confirm at each frame boundary (PocketSDR: re-sync at lock == fsync + 1800).
        if( epoch_count_ - epoch_at_sync_ >= static_cast<uint64_t>( FRAME_SYMS ) )
        {
            const int toi_next = ( toi_ + 1 ) % N_TOI;
            const int r        = match_frame( w, toi_next );
            if( r >= 0 )
            {
                toi_           = toi_next;
                rev_           = r;
                epoch_at_sync_ = epoch_count_;
                decode_frame( w ); // LDPC+CRC the frame (the proper per-frame confirmation, PocketSDR)
            }
            else
            {
                frame_synced_ = false; // failed re-confirm -> drop (the search re-establishes sync)
                logging::log(
                    logging::Level::Info,
                    fmt::format( "Navigation - GPS L1C PRN {:2d}  CNAV-2 frame sync LOST", satellite_id_ )
                );
            }
        }
        return;
    }

    // Not synced: search all 400 TOIs for the frame boundary in the current window.
    for( int toi = 0; toi < N_TOI; toi++ )
    {
        const int r = match_frame( w, toi );
        if( r >= 0 )
        {
            frame_synced_  = true;
            toi_           = toi;
            rev_           = r;
            epoch_at_sync_ = epoch_count_;
            logging::log(
                logging::Level::Info,
                fmt::format(
                    "Navigation - GPS L1C PRN {:2d}  CNAV-2 frame sync TOI={:d} ({})", satellite_id_, toi,
                    r == 1 ? "normal" : "reversed"
                )
            );
            decode_frame( w ); // LDPC+CRC the just-synced frame
            break;
        }
    }
}

bool Gps_l1c_decoder::bit_sync_found() const
{
    return epoch_count_ > 0; // one symbol per epoch; no sub-symbol bit sync needed for L1C
}

bool Gps_l1c_decoder::sw_loop() const
{
    return true; // unused by Pilot_tracker (it runs its own loop schedule)
}

bool Gps_l1c_decoder::frame_synced() const
{
    return frame_synced_;
}

const Ephemeris& Gps_l1c_decoder::ephemeris() const
{
    return eph_; // invalid until Phase 3 (LDPC decode + ephemeris parse)
}

const Iono* Gps_l1c_decoder::iono() const
{
    static const Iono dummy_iono {};
    return &dummy_iono; // no L1C iono decode yet
}

int Gps_l1c_decoder::get_current_bit_index() const
{
    return frame_synced_ ? static_cast<int>( ( epoch_count_ - epoch_at_sync_ ) % FRAME_SYMS ) : 0;
}

int Gps_l1c_decoder::get_current_ms_tick() const
{
    return 0;
}

const double Gps_l1c_decoder::ms_since_tow_update() const
{
    // ms elapsed since the eph.tow anchor epoch. One L1C epoch = one L1Cd code period = 10 ms (vs 1 ms
    // for L1CA), so scale the epoch delta by 10. The channel adds eph.tow to this (*1e-3) and projects to
    // the exact sample, so this only needs whole-epoch (10 ms) resolution.
    return static_cast<double>( epoch_count_ - active_eph_tow_anchor_epoch_ ) * 10.0;
}

#include "gps_l1ca_navigation.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include "almanac.h"
#include "logging.h"

// File-scope helpers (mirrors sdrnav.c / sdrnav_gps.c utilities)

// Unsigned bit extraction from byte buffer (MSB-first).
// mirrors RTKLIB getbitu() used throughout sdrnav_gps.c
static uint32_t getbitu( const uint8_t* buf, int pos, int len )
{
    uint32_t bits = 0;
    for( int i = pos; i < pos + len; i++ )
    {
        bits = ( bits << 1 ) | ( ( buf[i / 8] >> ( 7 - i % 8 ) ) & 1u );
    }

    return bits;
}

// Signed bit extraction (sign-extends result).
// mirrors RTKLIB getbits()
static int32_t getbits( const uint8_t* buf, int pos, int len )
{
    const uint32_t u = getbitu( buf, pos, len );
    if( len > 0 && ( u >> ( len - 1 ) ) )
    {
        return static_cast<int32_t>( u | ( ~0u << len ) );
    }

    return static_cast<int32_t>( u );
}

// Two-field unsigned concatenation.
// mirrors getbitu2() in sdrnav.c
static uint32_t getbitu2( const uint8_t* b, int p1, int l1, int p2, int l2 )
{
    return ( getbitu( b, p1, l1 ) << l2 ) | getbitu( b, p2, l2 );
}

// Two-field signed concatenation.
// mirrors getbits2() in sdrnav.c
static int32_t getbits2( const uint8_t* b, int p1, int l1, int p2, int l2 )
{
    if( getbitu( b, p1, 1 ) )
    {
        return static_cast<int32_t>( ( getbits( b, p1, l1 ) << l2 ) | getbitu( b, p2, l2 ) );
    }

    return static_cast<int32_t>( getbitu2( b, p1, l1, p2, l2 ) );
}

// Pack +/-1 bit array into bytes (MSB-first, -1->bit 1, +1->bit 0).
// mirrors bits2byte(bits, nbits, nbytes, right=0, out) in sdrnav.c
static void bits_to_bytes( const int* bits, int nbits, uint8_t* out, int nbytes )
{
    std::memset( out, 0, static_cast<size_t>( nbytes ) );
    for( int i = 0; i < nbits && i < nbytes * 8; i++ )
    {
        if( bits[i] < 0 )
        {
            out[i / 8] |= static_cast<uint8_t>( 0x80u >> ( i % 8 ) );
        }
    }
}

// GPS L1CA 30-bit word parity check (Hamming).
// Input: word[0..31] in +/-1, where word[0]=D29*, word[1]=D30* (carry from prev word),
//        word[2..25]=24 data bits, word[26..31]=6 parity bits.
// mirrors paritycheck_l1ca() in sdrnav_gps.c
static bool parity_check_word( const int* d )
{
    int p[6];
    // clang-format off
    p[0] = d[0] * d[2] * d[3] * d[4] * d[6] * d[7]  * d[11] * d[12] * d[13] * d[14] * d[15] * d[18] * d[19] * d[21] * d[24];
    p[1] = d[1] * d[3] * d[4] * d[5] * d[7] * d[8]  * d[12] * d[13] * d[14] * d[15] * d[16] * d[19] * d[20] * d[22] * d[25];
    p[2] = d[0] * d[2] * d[4] * d[5] * d[6] * d[8]  * d[9]  * d[13] * d[14] * d[15] * d[16] * d[17] * d[20] * d[21] * d[23];
    p[3] = d[1] * d[3] * d[5] * d[6] * d[7] * d[9]  * d[10] * d[14] * d[15] * d[16] * d[17] * d[18] * d[21] * d[22] * d[24];
    p[4] = d[1] * d[2] * d[4] * d[6] * d[7] * d[8]  * d[10] * d[11] * d[15] * d[16] * d[17] * d[18] * d[19] * d[22] * d[23] * d[25];
    p[5] = d[0] * d[4] * d[6] * d[7] * d[9] * d[10] * d[11] * d[12] * d[14] * d[16] * d[20] * d[23] * d[24] * d[25];
    // clang-format on

    int stat = 0;
    for( int i = 0; i < 6; i++ )
    {
        stat += ( p[i] - d[26 + i] );
    }

    return stat == 0;
}

// GPS ICD scale factors used in subframe decode
static constexpr double P2_5   = 1.0 / ( 1LL << 5 );
static constexpr double P2_11  = 1.0 / ( 1LL << 11 );
static constexpr double P2_17  = 1.0 / ( 1LL << 17 );
static constexpr double P2_19  = 1.0 / ( 1LL << 19 );
static constexpr double P2_20  = 1.0 / ( 1LL << 20 );
static constexpr double P2_21  = 1.0 / ( 1LL << 21 );
static constexpr double P2_23  = 1.0 / ( 1LL << 23 );
static constexpr double P2_38  = 1.0 / ( 1LL << 38 );
static constexpr double P2_24  = 1.0 / ( 1LL << 24 );
static constexpr double P2_27  = 1.0 / ( 1LL << 27 );
static constexpr double P2_29  = 1.0 / ( 1LL << 29 );
static constexpr double P2_30  = 1.0 / ( 1LL << 30 );
static constexpr double P2_31  = 1.0 / ( 1LL << 31 );
static constexpr double P2_33  = 1.0 / ( 1LL << 33 );
static constexpr double P2_43  = 1.0 / ( 1LL << 43 );
static constexpr double P2_55  = 1.0 / ( 1LL << 55 );
static constexpr double SC2RAD = 3.1415926535897932; // pi (semi-circles to radians)

// 10-bit GPS week rollover adjustment.
// mirrors adjgpsweek() in RTKLIB - uses system clock to pick the right epoch.
// mirrors the same logic called in decode_subfrm1()
static int adjust_gps_week( int transmitted_week ) // TODO: Using system clock doesn't work for IQ file...
{
    // GPS epoch: Jan 6, 1980 = Unix time 315964800 s
    const auto now              = static_cast<int64_t>( std::time( nullptr ) );
    const int  current_gps_week = static_cast<int>( ( now - 315964800LL ) / ( 7 * 24 * 3600 ) );
    const int  w                = std::max( current_gps_week, 1560 );

    return transmitted_week + ( w - transmitted_week % 1024 + 512 ) / 1024 * 1024;
}

// GPS L1CA subframe decoders.  All bit positions are 0-based into the 38-byte packed
// subframe buffer (300 bits, left-aligned).  Matches decode_subfrm1/2/3 in sdrnav_gps.c.

static void decode_sf1( const uint8_t* buf, Gps_ephemeris& e )
{
    e.tow         = getbitu( buf, 30, 17 ) * 6.0;
    e.week        = adjust_gps_week( static_cast<int>( getbitu( buf, 60, 10 ) ) );
    e.iodc        = static_cast<int>( getbitu2( buf, 82, 2, 210, 8 ) );
    e.group_delay = getbits( buf, 196, 8 ) * P2_31;
    e.toc         = getbitu( buf, 218, 16 ) * 16.0;
    e.af2         = getbits( buf, 240, 8 ) * P2_55;
    e.af1         = getbits( buf, 248, 16 ) * P2_43;
    e.af0         = getbits( buf, 270, 22 ) * P2_31;
}

static void decode_sf2( const uint8_t* buf, Gps_ephemeris& e )
{
    e.tow     = getbitu( buf, 30, 17 ) * 6.0;
    e.iode    = static_cast<int>( getbitu( buf, 60, 8 ) );
    e.crs     = getbits( buf, 68, 16 ) * P2_5;
    e.delta_n = getbits( buf, 90, 16 ) * P2_43 * SC2RAD;
    e.m0      = getbits2( buf, 106, 8, 120, 24 ) * P2_31 * SC2RAD;
    e.cuc     = getbits( buf, 150, 16 ) * P2_29;
    e.e       = getbitu2( buf, 166, 8, 180, 24 ) * P2_33;
    e.cus     = getbits( buf, 210, 16 ) * P2_29;
    e.sqrt_a  = getbitu2( buf, 226, 8, 240, 24 ) * P2_19;
    e.toe     = getbitu( buf, 270, 16 ) * 16.0;
}

static void decode_sf3( const uint8_t* buf, Gps_ephemeris& e )
{
    e.tow      = getbitu( buf, 30, 17 ) * 6.0;
    e.cic      = getbits( buf, 60, 16 ) * P2_29;
    e.omega0   = getbits2( buf, 76, 8, 90, 24 ) * P2_31 * SC2RAD;
    e.cis      = getbits( buf, 120, 16 ) * P2_29;
    e.i0       = getbits2( buf, 136, 8, 150, 24 ) * P2_31 * SC2RAD;
    e.crc      = getbits( buf, 180, 16 ) * P2_5;
    e.omega    = getbits2( buf, 196, 8, 210, 24 ) * P2_31 * SC2RAD;
    e.omegadot = getbits( buf, 240, 24 ) * P2_43 * SC2RAD;
    e.iode     = static_cast<int>( getbitu( buf, 270, 8 ) );
    e.idot     = getbits( buf, 278, 14 ) * P2_43 * SC2RAD;
}

// Subframe 4 is subcommutated over 25 pages; only page 18 (identified by SV ID 56 in word 3)
// carries the Klobuchar ionosphere model + leap seconds. Every other SF4/SF5 page is almanac
// (coarse orbits for all SVs), which PVT does not need yet. Bit positions per IS-GPS-200
// 20.3.3.5.1.7 / 20.3.3.5.2.4, scale factors mirror RTKLIB decode_subfrm4().
static void decode_sf4( const uint8_t* buf, Iono& iono )
{
    const int sv_id = static_cast<int>( getbitu( buf, 62, 6 ) ); // page id (data ID is bits 60-61)
    if( sv_id != 56 )
    {
        return; // almanac / other page
    }

    iono.alpha[0]     = getbits( buf, 68, 8 ) * P2_30;
    iono.alpha[1]     = getbits( buf, 76, 8 ) * P2_27;
    iono.alpha[2]     = getbits( buf, 90, 8 ) * P2_24;
    iono.alpha[3]     = getbits( buf, 98, 8 ) * P2_24;
    iono.beta[0]      = getbits( buf, 106, 8 ) * 2048.0;  // 2^11
    iono.beta[1]      = getbits( buf, 120, 8 ) * 16384.0; // 2^14
    iono.beta[2]      = getbits( buf, 128, 8 ) * 65536.0; // 2^16
    iono.beta[3]      = getbits( buf, 136, 8 ) * 65536.0; // 2^16
    iono.leap_seconds = getbits( buf, 240, 8 );           // dt_LS (word 9)
    iono.model        = Iono::Model::Klobuchar;
    iono.valid        = true;
}

// Decode a GPS ALMANAC page (the coarse constellation-wide orbits subcommutated over SF4/SF5). Per
// IS-GPS-200 20.3.3.5: the per-SV almanac is SF5 pages 1-24 (SVs 1-24) + SF4 pages 2-10 (SVs 25-32),
// each identified by the SV ID in word 3 (bits 62-67). Word-aligned bit positions (this codebase keeps the
// 6 parity bits per 30-bit word, so words start every 30 bits); field order + scales mirror RTKLIB
// decode_alm_sat(). `week` is the current GPS week (from the decoded ephemeris) - the almanac page carries
// only toa-within-week. Returns the PRN (and fills alm[prn]) if this was an almanac page, else 0.
static int decode_almanac( const uint8_t* buf, int sfn, int week, std::map<int, Almanac>& alm )
{
    const int  svid        = static_cast<int>( getbitu( buf, 62, 6 ) ); // page SV ID (== the PRN for almanac pages)
    const bool is_alm_page = ( sfn == 5 && svid >= 1 && svid <= 24 ) || ( sfn == 4 && svid >= 25 && svid <= 32 );
    if( !is_alm_page )
    {
        return 0; // iono (p18), health/config (p25), or a reserved page
    }

    Almanac a;
    a.constellation = Constellation::Gps;
    a.prn           = static_cast<Satellite_id>( svid );
    a.e             = getbitu( buf, 68, 16 ) * P2_21;
    a.toa           = getbitu( buf, 90, 8 ) * 4096.0;                    // 2^12
    a.i0            = ( 0.3 + getbits( buf, 98, 16 ) * P2_19 ) * SC2RAD; // 0.3 semicircles reference + delta_i
    a.omegadot      = getbits( buf, 120, 16 ) * P2_38 * SC2RAD;
    a.health        = static_cast<int>( getbitu( buf, 136, 8 ) );
    a.sqrt_a        = getbitu( buf, 150, 24 ) * P2_11;
    a.omega0        = getbits( buf, 180, 24 ) * P2_23 * SC2RAD;
    a.omega         = getbits( buf, 210, 24 ) * P2_23 * SC2RAD;
    a.m0            = getbits( buf, 240, 24 ) * P2_23 * SC2RAD;
    // af0 (11-bit, signed) is split: 8 MSBs at 270, 3 LSBs at 289 (af1 sits between). Mirror RTKLIB exactly.
    const int af0_8 = static_cast<int>( getbits( buf, 270, 8 ) );
    a.af1           = getbits( buf, 278, 11 ) * P2_38;
    a.af0           = getbitu( buf, 289, 3 ) * P2_17 + af0_8 * P2_20;
    a.week          = week;
    a.valid         = true;
    alm[svid]       = a;
    return svid;
}

// Dispatch to the appropriate subframe decoder based on the HOW subframe ID.
// mirrors decode_frame_l1ca() in sdrnav_gps.c
// Returns the subframe number (1-5), or 0 on error.
static int decode_gps_frame( const uint8_t* buf, Gps_ephemeris& eph, Iono& iono )
{
    const int id = static_cast<int>( getbitu( buf, 49, 3 ) );
    switch( id )
    {
    case 1:
        decode_sf1( buf, eph );
        break;
    case 2:
        decode_sf2( buf, eph );
        break;
    case 3:
        decode_sf3( buf, eph );
        break;
    case 4:
        decode_sf4( buf, iono ); // page 18 only (iono + leap seconds); rest is almanac
        break;
    default:
        break; // SF5 - almanac, not needed for pseudorange yet
    }
    return id;
}

// Gps_l1ca_decoder implementation

Gps_l1ca_decoder::Gps_l1ca_decoder( Satellite_id prn, const Signal_params& params )
    : satellite_id_( prn ),
      params_( params )
{
    // Stamp constellation + PRN on both ephemeris copies so the observables/PVT pipeline can
    // dispatch the orbit/clock per constellation (and not rely on uninitialised fields).
    eph_.constellation = eph_current_.constellation = Constellation::Gps;
    eph_.prn = eph_current_.prn = prn;
}

// process
// Top-level entry point called every tracking epoch.
// mirrors sdrnavigation() in sdrnav.c
void Gps_l1ca_decoder::process( double prompt_i, double prompt_i_prev )
{
    // bit_count_ = current position within the 20-epoch bit period
    // mirrors: sdr->nav.biti = cnt % sdr->nav.rate
    bit_count_ = static_cast<int>( epoch_count_ % NAV_RATE );

    // Bit synchronisation: wait 2 seconds before starting (mirrors: cnt>2000/(ctime*1000))
    if( !bit_sync_found_ && epoch_count_ > 2000 )
    {
        check_sync( prompt_i, prompt_i_prev );
    }

    if( bit_sync_found_ )
    {
        check_bit( prompt_i );

        // When a bit was just decided, try to find the preamble and decode.
        // For L1CA, predecodefec() is a no-op (fbitsdec = fbits), so we go
        // straight to findpreamble().
        if( sw_sync_ )
        {
            if( preamble_found_ )
            {
                nav_bit_count_++;

                if( nav_bit_count_ >= 300 )
                {
                    nav_bit_count_ = 0;
                }
            }

            if( find_preamble() )
            {
                decode_subframe();
            }
        }
    }

    // NOTE: epoch_count_ is incremented at the END of process() here (L1C does it at the TOP). This
    // asymmetry is load-bearing: ms_since_tow_update() is read by the Channel AFTER process() returns, so
    // the eph TOW anchor uses epoch_count_+1 (see parse in decode_subframe) to match. Keep them in sync.
    ++epoch_count_;
}

int Gps_l1ca_decoder::get_current_bit_index() const
{
    return nav_bit_count_;
}

int Gps_l1ca_decoder::get_current_ms_tick() const
{
    if( !bit_sync_found_ )
    {
        return 0;
    }

    return ( bit_count_ - sync_index_ + 20 ) % NAV_RATE;
}

const double Gps_l1ca_decoder::ms_since_tow_update() const
{
    return ( epoch_count_ - active_eph_tow_anchor_epoch_ );
}

// check_sync
// Detect nav bit boundaries by voting on sign-change positions.
// mirrors checksync() in sdrnav.c (non-BeiDou path)
void Gps_l1ca_decoder::check_sync( double prompt_i, double prompt_i_prev )
{
    if( prompt_i_prev * prompt_i < 0.0 )
    {
        bit_sync_votes_[bit_count_] += 1;

        int max_votes = 0;
        int max_idx   = 0;
        for( int i = 0; i < NAV_RATE; i++ )
        {
            if( bit_sync_votes_[i] > max_votes )
            {
                max_votes = bit_sync_votes_[i];
                max_idx   = i;
            }
        }

        if( max_votes > NAV_SYNCTH )
        {
            // mirrors: nav->synci = biti - 1 (wrap)
            sync_index_ = max_idx - 1;
            if( sync_index_ < 0 )
            {
                sync_index_ += NAV_RATE;
            }

            bit_sync_found_ = true;
        }
    }
}

// check_bit
// Accumulate prompt I over a 20-epoch bit period and decide each bit.
// mirrors checkbit() in sdrnav.c (loopms = NAV_LOOP_MS = 10)
void Gps_l1ca_decoder::check_bit( double prompt_i )
{
    sw_sync_ = false;
    sw_loop_ = false;

    const int diffi = bit_count_ - sync_index_;

    // diffi == 1 (normal): start of a new bit
    // diffi == -(NAV_RATE-1) (wrapped): same condition with bit_count_=0, sync_index_=19
    if( diffi == 1 || diffi == -( NAV_RATE - 1 ) )
    {
        bit_accum_ = prompt_i;
        bit_cnt_   = 1;
    }
    else
    {
        bit_accum_ += prompt_i;
    }

    // Loop filter trigger - mirrors: if (nav->cnt % loopms == 0) nav->swloop = ON
    sw_loop_ = ( bit_cnt_ % NAV_LOOP_MS == 0 );

    // diffi == 0: end of bit - decide the bit value and store it
    if( diffi == 0 )
    {
        // Store the RAW correlator-sign bit, NOT polarity-corrected. The 180-deg carrier ambiguity is
        // resolved uniformly at decode time (find_preamble / decode_subframe multiply the whole frame
        // by the preamble-matched polarity, and GPS parity is invariant to a UNIFORM inversion via the
        // D30* rule). Baking the per-subframe polarity_ in here instead made the buffer non-uniform
        // across subframe boundaries - the 2 carry bits (D29*/D30*) belong to the PREVIOUS subframe and
        // were stored under its polarity, so word 0 (TLM) parity failed every time the data inversion
        // flipped between subframes -> ~half the subframes silently dropped -> eph delayed a full 30 s
        // cycle. Mirrors GNSS-SDRLIB checkbit (flagpol is SBAS-only; L1CA stores the raw bit).
        const int bit = ( bit_accum_ < 0.0 ) ? -1 : 1;

        for( size_t i = 0; i < frame_bits_.size() - 1; ++i )
        {
            frame_bits_[i] = frame_bits_[i + 1];
        }
        frame_bits_.back() = bit;

        sw_sync_ = true;
    }

    ++bit_cnt_;
}

// find_preamble
// Search for the GPS L1CA TLM preamble and verify parity on all 10 words.
// mirrors findpreamble() + paritycheck() for L1CA in sdrnav.c
bool Gps_l1ca_decoder::find_preamble()
{
    // Correlate frame_bits_[NAV_ADDFLEN .. NAV_ADDFLEN+NAV_PRELEN-1] with preamble.
    // mirrors: for (i=0; i<prelen; i++) corr += fbitsdec[addflen+i] * prebits[i]
    int corr = 0;
    for( int i = 0; i < NAV_PRELEN; i++ )
    {
        corr += frame_bits_[NAV_ADDFLEN + i] * PREAMBLE[i];
    }

    if( std::abs( corr ) != NAV_PRELEN )
    {
        return false;
    }

    const int candidate_polarity = ( corr > 0 ) ? 1 : -1;

    // Parity check all 10 words with this polarity candidate.
    // mirrors: paritycheck(nav) for CTYPE_L1CA
    // Each 32-bit window frame_bits_[w*30 .. w*30+31] provides:
    //   [0..1] = D29*/D30* carry from previous word, [2..31] = 30-bit word
    int ok = 0;
    for( int w = 0; w < 10; w++ )
    {
        int word[32];
        for( int j = 0; j < 32; j++ )
        {
            word[j] = candidate_polarity * frame_bits_[w * 30 + j];
        }

        // D30* correction: if carry bit 1 is -1, invert the 24 data bits
        // mirrors: if (bits[i*30+1]==-1) for (j=2;j<26;j++) bits[i*30+j]*=-1
        if( word[1] == -1 )
        {
            for( int j = 2; j < 26; j++ )
            {
                word[j] = -word[j];
            }
        }

        if( parity_check_word( word ) )
        {
            ok++;
        }
    }

    if( ok != 10 )
    {
        return false;
    }

    // Cadence gate: once locked onto the preamble, a genuine TLM preamble recurs EXACTLY every 300 bits
    // (one subframe), i.e. at nav_bit_count_ == 0. Reject a parity-passing match that lands OFF that
    // boundary - it is a false preamble (the 8-bit pattern + 10 words happening to pass parity on
    // misaligned bits). Such matches are rare on a clean capture but rampant on some sims (Skydel), and
    // each decodes a garbage HOW TOW that re-bootstraps the TOW anchor, so tow_confirmed_ never latches
    // and no ephemeris is ever published. The first preamble (preamble_found_ == false) bootstraps the
    // cadence and is accepted at any phase; a missed (bit-errored) subframe is harmless - nav_bit_count_
    // free-runs mod 300, so the next real preamble still lands at 0.
    if( preamble_found_ && nav_bit_count_ != 0 )
    {
        return false;
    }

    polarity_ = candidate_polarity;

    if( !preamble_found_ )
    {
        preamble_found_ = true;
        logging::log(
            logging::Level::Info,
            fmt::format( "Navigation - {} PRN {:2d}  preamble found  polarity={:+d}", params_.name, satellite_id_, polarity_ )
        );
    }

    return true;
}

// decode_subframe
// Apply D* bit inversion, pack frame buffer to bytes, dispatch to subframe decoder.
// mirrors decode_l1ca() + decode_frame_l1ca() in sdrnav_gps.c
void Gps_l1ca_decoder::decode_subframe()
{
    // Working copy with polarity and per-word D* correction applied.
    // mirrors: for (i=0;i<10;i++) { if (fbitsdec[i*30+1]==-1) invert fbitsdec[i*30+2..25] }
    std::array<int, NAV_FLEN + NAV_ADDFLEN> bits;
    for( int i = 0; i < NAV_FLEN + NAV_ADDFLEN; i++ )
    {
        bits[i] = polarity_ * frame_bits_[i];
    }

    for( int w = 0; w < 10; w++ )
    {
        if( bits[w * 30 + 1] == -1 )
        {
            for( int j = 2; j < 26; j++ )
            {
                bits[w * 30 + j] = -bits[w * 30 + j];
            }
        }
    }

    // Pack to bytes: skip NAV_ADDFLEN carry bits, pack the 300 data bits.
    // mirrors: bits2byte(&fbitsdec[addflen], flen, 38, 0, bin)
    uint8_t bin[38] = {};
    bits_to_bytes( &bits[NAV_ADDFLEN], NAV_FLEN, bin, 38 );

    // The frame already passed per-word parity in find_preamble(), so decode it into the WORKING
    // ephemeris (eph_current_) immediately - including the very first subframe after lock, so a clean
    // SF1 is captured on the first pass instead of being discarded (which cost a whole ~30 s frame
    // before the next SF1, dominating TTFF). The TOW check now only DROPS the corrupted-HOW glitch
    // (valid parity, garbage TOW); publishing eph_current_ -> eph_ is what waits for a trusted set.
    // The TOW is the HOW count (bits 30..46) x 6 s. A valid TOW-count is 0..100799 (one week); an
    // out-of-range count is an impossible time, so the subframe is corrupted / a false-preamble match -
    // drop it WITHOUT touching the TOW anchor (belt-and-braces with the cadence gate in find_preamble).
    const uint32_t tow_count = static_cast<uint32_t>( getbitu( bin, 30, 17 ) );
    if( tow_count > 100799 )
    {
        return;
    }
    const double     tow        = tow_count * 6.0;
    const Tow_status tow_status = classify_tow( tow );
    if( tow_status == Tow_status::CorruptRejected )
    {
        return;
    }

    const bool had_iono = iono_.valid;
    const int  sfn      = decode_gps_frame( bin, eph_current_, iono_ );

    // Almanac (SF4/SF5 pages) - coarse constellation-wide orbits for acquisition aiding (not PVT). One
    // tracked SV slowly fills the whole set; aggregated receiver-wide by the engine via almanac().
    if( sfn == 4 || sfn == 5 )
    {
        const size_t known_before = almanac_.size();
        if( int aprn = decode_almanac( bin, sfn, eph_current_.week, almanac_ ) )
        {
            if( almanac_.size() > known_before ) // log once per newly-known SV (pages re-decode every cycle)
            {
                const Almanac& a = almanac_[aprn];
                logging::log(
                    logging::Level::Info,
                    fmt::format(
                        "Navigation - GPS PRN {:2d} decoded ALMANAC for PRN {:2d} (toa={:.0f} sqrtA={:.2f} "
                        "e={:.3e} health={}); {} SVs known",
                        satellite_id_,
                        a.prn,
                        a.toa,
                        a.sqrt_a,
                        a.e,
                        a.health,
                        almanac_.size()
                    )
                );
            }
        }
    }

    if( sfn >= 1 && sfn <= 3 )
    {
        sf_decoded_ |= ( 1 << ( sfn - 1 ) );
        if( sf_decoded_ == 0x7 )
        {
            eph_current_.valid = true;
        }
        // Publish the working copy to the PVT-visible eph_ only once the set is COMPLETE (SF1+2+3),
        // the TOW lock is CONFIRMED (two consecutive consistent subframes - so a lone unconfirmed lock
        // or a glitch never reaches PVT), and the set is self-consistent: IODE (8-bit, SF2/SF3) ==
        // low 8 bits of IODC (10-bit, SF1) per IS-GPS-200 20.3.4.4. That consistency check is also the
        // safety net that makes decoding the first subframe optimistically safe: a glitched SF1 (valid
        // parity, bad data) yields a mismatched IODC and is simply not published, then re-decoded next
        // pass. (Comparing the full iodc==iode only matched when iodc < 256, so the commit almost never
        // fired - which gated has_observable() and the TOW-anchor refresh.)
        if( eph_current_.valid && tow_confirmed_ && ( eph_current_.iodc & 0xFF ) == eph_current_.iode )
        {
            eph_ = eph_current_;
            // epoch_count_ is incremented at the END of process() (the ++ below the decode call), but the
            // channel reads ms_since_tow_update() = epoch_count_ - this_anchor AFTER process() returns, i.e.
            // with the post-increment value. Anchor to epoch_count_ + 1 so it matches that convention; else
            // every L1CA transmit time is one epoch (1 ms) high. Invisible to L1CA-only PVT (a common bias
            // absorbed into the clock state) but a 300 km inconsistency when combined with L1C for the same
            // SV. (L1C's decoder ++s epoch_count_ at the TOP of its process(), so it has no such skew.)
            active_eph_tow_anchor_epoch_ = epoch_count_ + 1;
        }
    }

    // Log every subframe (1-5) so SF4/SF5 reception is visible, even though only 1-3 feed the
    // ephemeris. Week comes only from subframe 1; show it as unknown until SF1 has been decoded.
    if( sfn >= 1 && sfn <= 5 )
    {
        const std::string week_str =
            ( sf_decoded_ & 0x1 ) ? fmt::format( "week={:d}", eph_current_.week ) : std::string( "week=?" );
        logging::log(
            logging::Level::Info,
            fmt::format(
                "Navigation - {} PRN {:2d}  subframe {:d}  tow={:.1f}s  {}",
                params_.name,
                satellite_id_,
                sfn,
                tow, // HOW TOW (set for every subframe, including SF4/5 which don't fill eph_)
                week_str
            )
        );
    }

    // Klobuchar iono is broadcast on SF4 page 18 only (every ~12.5 min); announce it once.
    if( !had_iono && iono_.valid )
    {
        logging::log(
            logging::Level::Info,
            fmt::format(
                "Navigation - {} PRN {:2d}  iono alpha=[{:.2e} {:.2e} {:.2e} {:.2e}] beta=[{:.2e} {:.2e} {:.2e} "
                "{:.2e}]  leap={:d}s",
                params_.name,
                satellite_id_,
                iono_.alpha[0],
                iono_.alpha[1],
                iono_.alpha[2],
                iono_.alpha[3],
                iono_.beta[0],
                iono_.beta[1],
                iono_.beta[2],
                iono_.beta[3],
                iono_.leap_seconds
            )
        );
    }
}

// Classify a subframe's HOW TOW against the continuity anchor. The tracked epoch count is locked to
// the signal (one epoch = 1 ms of transmit time), so a valid HOW TOW must equal the last good TOW
// plus the epochs elapsed since. Logic:
//   - First TOW: a tentative anchor (Provisional) - decode it (parity already validated it), but it
//     is not yet trusted for PVT; the published-ephemeris gate also requires tow_confirmed_.
//   - A later TOW that AGREES with the anchor: Continuous - confirm + refresh the anchor.
//   - A later TOW that DISAGREES while the anchor is CONFIRMED: it is the corrupted one (the clock is
//     trustworthy) - CorruptRejected; the caller DROPS it WITHOUT touching the anchor, so the next
//     good subframe is still validated. This is the fix for a glitch poisoning the anchor.
//   - A TOW that disagrees while the anchor is only tentative: can't tell which is bad, so re-bootstrap
//     on the newer one (Provisional).
// Net: parity-valid subframes are decoded into the working ephemeris immediately (a clean SF1 is not
// thrown away), but only a confirmed, consistent set is published to PVT, and a corrupt TOW costs
// only itself.
Gps_l1ca_decoder::Tow_status Gps_l1ca_decoder::classify_tow( double tow )
{
    const double     epoch_period_s = params_.code_period_s; // one tracking epoch (1 ms for L1 C/A)
    constexpr double TOW_TOL_S      = 1.0;                   // << 6 s subframe spacing, >> tracked-time jitter

    if( !tow_anchored_ )
    {
        tow_anchor_value_ = tow;
        tow_anchor_epoch_ = epoch_count_;
        tow_anchored_     = true;
        return Tow_status::Provisional; // tentative anchor, not yet trusted
    }

    const double elapsed_s = static_cast<double>( epoch_count_ - tow_anchor_epoch_ ) * epoch_period_s;
    const double predicted = tow_anchor_value_ + elapsed_s;
    if( std::fabs( tow - predicted ) < TOW_TOL_S )
    {
        tow_anchor_value_ = tow; // agrees -> trusted; confirm + refresh
        tow_anchor_epoch_ = epoch_count_;
        tow_confirmed_    = true;
        return Tow_status::Continuous;
    }

    if( tow_confirmed_ )
    {
        // Trustworthy anchor -> this TOW is the corrupted one. Reject, keep the anchor intact.
        logging::log(
            logging::Level::Warning,
            fmt::format(
                "Navigation - {} PRN {:2d}  rejected subframe: TOW {:.1f}s != expected {:.1f}s (corrupted HOW)",
                params_.name,
                satellite_id_,
                tow,
                predicted
            )
        );
        return Tow_status::CorruptRejected;
    }

    // Tentative anchor disagreed -> re-bootstrap on the newer observation.
    tow_anchor_value_ = tow;
    tow_anchor_epoch_ = epoch_count_;
    return Tow_status::Provisional;
}

#ifdef ENABLE_UNIT_TESTS
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include "orbit.h"

TEST_CASE( "gps_getbitu_extracts_msb_first", "[gps][nav][bits]" )
{
    // 0xB4 = 1011 0100, 0xF0 = 1111 0000.
    const uint8_t buf[2] = { 0xB4, 0xF0 };
    REQUIRE( getbitu( buf, 0, 8 ) == 0xB4u );
    REQUIRE( getbitu( buf, 0, 4 ) == 0b1011u );
    REQUIRE( getbitu( buf, 4, 4 ) == 0b0100u );
    REQUIRE( getbitu( buf, 1, 3 ) == 0b011u );
    REQUIRE( getbitu( buf, 6, 4 ) == 0b0011u ); // crosses the byte boundary (..00|11..)
    REQUIRE( getbitu( buf, 0, 16 ) == 0xB4F0u );
}

TEST_CASE( "gps_getbits_sign_extends", "[gps][nav][bits]" )
{
    const uint8_t ff[1]  = { 0xFF }; // 8-bit -1
    const uint8_t x80[1] = { 0x80 }; // 8-bit -128
    const uint8_t x40[1] = { 0x40 }; // 8-bit +64
    REQUIRE( getbits( ff, 0, 8 ) == -1 );
    REQUIRE( getbits( x80, 0, 8 ) == -128 );
    REQUIRE( getbits( x40, 0, 8 ) == 64 );
    REQUIRE( getbits( ff, 0, 4 ) == -1 );    // 0b1111 sign-extended
    REQUIRE( getbits( x40, 1, 2 ) == -2 );   // bits 1..2 of 0100 0000 = 0b10 = -2
    REQUIRE( getbitu( x80, 0, 8 ) == 128u ); // unsigned counterpart differs
}

TEST_CASE( "gps_parity_check_word_detects_single_bit_errors", "[gps][nav][parity]" )
{
    // The all-+1 word is self-consistent: every parity equation is a product of +1 -> +1, so
    // the six parity bits must be +1 for parity_check_word to pass. This is hand-verifiable.
    int word[32];
    for( int i = 0; i < 32; ++i )
    {
        word[i] = 1;
    }
    REQUIRE( parity_check_word( word ) );

    // Corrupting any single bit (data, carry, or parity) must be detected.
    for( int i = 0; i < 32; ++i )
    {
        int w[32];
        std::copy( std::begin( word ), std::end( word ), std::begin( w ) );
        w[i] = -1;
        INFO( "flipped bit " << i );
        REQUIRE_FALSE( parity_check_word( w ) );
    }
}

TEST_CASE( "gps_parity_round_trip_over_data", "[gps][nav][parity]" )
{
    // Build a valid word for arbitrary data + carries by evaluating the ICD parity equations
    // (IS-GPS-200 Table 20-XIV), then confirm the decoder accepts it and rejects a data flip.
    auto encode = []( int* d )
    {
        // clang-format off
        d[26] = d[0]*d[2]*d[3]*d[4]*d[6]*d[7]*d[11]*d[12]*d[13]*d[14]*d[15]*d[18]*d[19]*d[21]*d[24];
        d[27] = d[1]*d[3]*d[4]*d[5]*d[7]*d[8]*d[12]*d[13]*d[14]*d[15]*d[16]*d[19]*d[20]*d[22]*d[25];
        d[28] = d[0]*d[2]*d[4]*d[5]*d[6]*d[8]*d[9]*d[13]*d[14]*d[15]*d[16]*d[17]*d[20]*d[21]*d[23];
        d[29] = d[1]*d[3]*d[5]*d[6]*d[7]*d[9]*d[10]*d[14]*d[15]*d[16]*d[17]*d[18]*d[21]*d[22]*d[24];
        d[30] = d[1]*d[2]*d[4]*d[6]*d[7]*d[8]*d[10]*d[11]*d[15]*d[16]*d[17]*d[18]*d[19]*d[22]*d[23]*d[25];
        d[31] = d[0]*d[4]*d[6]*d[7]*d[9]*d[10]*d[11]*d[12]*d[14]*d[16]*d[20]*d[23]*d[24]*d[25];
        // clang-format on
    };

    // A pseudo-random but fixed bit pattern (both carry polarities exercised across the cases).
    const uint32_t patterns[2] = { 0xA5C3F0u, 0x5A3C0Fu };
    for( int carry = 0; carry < 2; ++carry )
    {
        for( uint32_t pat : patterns )
        {
            int d[32];
            d[0] = ( carry & 1 ) ? -1 : 1; // D29*
            d[1] = ( carry & 1 ) ? 1 : -1; // D30*
            for( int i = 0; i < 24; ++i )
            {
                d[2 + i] = ( ( pat >> ( 23 - i ) ) & 1u ) ? -1 : 1;
            }
            encode( d );
            REQUIRE( parity_check_word( d ) );

            int flipped[32];
            std::copy( std::begin( d ), std::end( d ), std::begin( flipped ) );
            flipped[10] = -flipped[10]; // flip one data bit
            REQUIRE_FALSE( parity_check_word( flipped ) );
        }
    }
}

namespace
{
// Parse a 76-char hex string into the 38-byte packed subframe buffer.
std::array<uint8_t, 38> sf_from_hex( const char* hex )
{
    auto nib = []( char ch ) -> int
    {
        if( ch >= '0' && ch <= '9' )
            return ch - '0';
        return ( ch | 0x20 ) - 'a' + 10; // lower-case a-f
    };
    std::array<uint8_t, 38> buf {};
    for( int i = 0; i < 38; ++i )
    {
        buf[i] = static_cast<uint8_t>( ( nib( hex[2 * i] ) << 4 ) | nib( hex[2 * i + 1] ) );
    }
    return buf;
}

// Real GPS L1 C/A subframes 1/2/3 for PRN 1, captured from the CTTC capture (one consistent
// ephemeris set: IODC 11 == IODE 11). Used as a decode regression + physical-plausibility fixture.
constexpr const char* CTTC_PRN1_SF1 = "8B10B435DFF89A3B19000009C420F3E2494088818F0567D8A120C2D6DA0D00001E0806E136B0";
constexpr const char* CTTC_PRN1_SF2 = "8B10B435DFFAACB0B017558CD3DD81791C4F140558021DF36C20855528460D9823456DA156F0";
constexpr const char* CTTC_PRN1_SF3 = "8B10B435DFFCBBBFFECED7081E8E430010278C872DEB516CC0C77370D3F1FFA8A8082C3A6E70";
} // namespace

TEST_CASE( "gps_decode_real_subframes_prn1", "[gps][nav][ephemeris][cttc]" )
{
    const auto sf1 = sf_from_hex( CTTC_PRN1_SF1 );
    const auto sf2 = sf_from_hex( CTTC_PRN1_SF2 );
    const auto sf3 = sf_from_hex( CTTC_PRN1_SF3 );

    Gps_ephemeris e;
    Iono          iono;
    REQUIRE( decode_gps_frame( sf1.data(), e, iono ) == 1 ); // dispatch by the HOW subframe id
    REQUIRE( decode_gps_frame( sf2.data(), e, iono ) == 2 );
    REQUIRE( decode_gps_frame( sf3.data(), e, iono ) == 3 );

    // Issue-of-data consistency (the gate the decoder uses to publish an ephemeris).
    REQUIRE( e.iodc == 11 );
    REQUIRE( e.iode == 11 );
    REQUIRE( ( e.iodc & 0xFF ) == e.iode );

    // Regression: pin the decoded fields (scale factors + bit positions) against the captured values.
    // The 10-bit broadcast week is deterministic; the absolute week from adjust_gps_week() is NOT (it
    // resolves the rollover against the system clock), so check the broadcast value modulo 1024.
    REQUIRE( ( e.week & 1023 ) == 710 );
    REQUIRE( e.toc == Catch::Approx( 374400.0 ) );
    REQUIRE( e.toe == Catch::Approx( 374400.0 ) );
    REQUIRE( e.sqrt_a == Catch::Approx( 5153.699286 ).epsilon( 1e-9 ) );
    REQUIRE( e.e == Catch::Approx( 1.702986890e-3 ).epsilon( 1e-7 ) );
    REQUIRE( e.m0 == Catch::Approx( 2.907767059 ).epsilon( 1e-7 ) );
    REQUIRE( e.delta_n == Catch::Approx( 4.691266839e-9 ).epsilon( 1e-6 ) );
    REQUIRE( e.omega0 == Catch::Approx( -0.4632164247 ).epsilon( 1e-7 ) );
    REQUIRE( e.omega == Catch::Approx( 0.3142515846 ).epsilon( 1e-7 ) );
    REQUIRE( e.omegadot == Catch::Approx( -7.986046937e-9 ).epsilon( 1e-6 ) );
    REQUIRE( e.i0 == Catch::Approx( 0.9604440504 ).epsilon( 1e-7 ) );
    REQUIRE( e.idot == Catch::Approx( 3.335853237e-10 ).epsilon( 1e-6 ) );
    REQUIRE( e.crs == Catch::Approx( 11.656250 ) );
    REQUIRE( e.crc == Catch::Approx( 182.375000 ) );
    REQUIRE( e.af0 == Catch::Approx( 1.312186942e-5 ).epsilon( 1e-7 ) );
    REQUIRE( e.group_delay == Catch::Approx( 8.381903172e-9 ).epsilon( 1e-6 ) );

    // Physical plausibility (non-circular): the decoded Keplerian elements must describe a real GPS
    // orbit. If any scale factor / bit field were wrong, these would not hold.
    REQUIRE( e.e < 0.03 );                                                      // GPS eccentricity is small
    REQUIRE( e.sqrt_a * e.sqrt_a == Catch::Approx( 26.56e6 ).epsilon( 0.01 ) ); // semi-major ~26 560 km
    REQUIRE( e.i0 > 0.9 );                                                      // inclination ~55 deg
    REQUIRE( e.i0 < 1.05 );
    e.constellation = Constellation::Gps;
    const Ecef   p  = orbit::satellite_ecef_pos( e, e.toe, Constellation::Gps );
    const double r  = std::sqrt( p.x * p.x + p.y * p.y + p.z * p.z );
    REQUIRE( r == Catch::Approx( 26.56e6 ).epsilon( 0.02 ) ); // a real orbit radius from decoded elements
}

TEST_CASE( "gps_how_tow_count_range_guard", "[gps][nav][tow]" )
{
    // decode_subframe() reads the HOW TOW-count from bits 30..46 (17 bits) of the packed subframe and
    // drops the subframe if the count exceeds 100799 (an impossible time) - this is what stops a
    // false-preamble's garbage HOW from poisoning the TOW anchor (the Skydel GPS-PVT fix). This test
    // pins the field position, the physical cadence, and the threshold boundary - all on real data.
    const auto sf1 = sf_from_hex( CTTC_PRN1_SF1 );
    const auto sf2 = sf_from_hex( CTTC_PRN1_SF2 );
    const auto sf3 = sf_from_hex( CTTC_PRN1_SF3 );

    // (1) Real subframes carry a VALID TOW-count, so the guard never drops genuine data.
    const uint32_t c1 = getbitu( sf1.data(), 30, 17 );
    const uint32_t c2 = getbitu( sf2.data(), 30, 17 );
    const uint32_t c3 = getbitu( sf3.data(), 30, 17 );
    REQUIRE( c1 <= 100799u );
    REQUIRE( c2 <= 100799u );
    REQUIRE( c3 <= 100799u );

    // (2) Physical truth (non-circular): the HOW TOW-count is in 6 s units and names the NEXT subframe,
    // so three consecutive subframes increment by exactly 1. If the bit position were wrong this fails.
    REQUIRE( c2 == c1 + 1 );
    REQUIRE( c3 == c2 + 1 );

    // (3) Threshold boundary. Set bits 30..46 (RTKLIB MSB-first, matching getbitu) to chosen counts and
    // confirm the > 100799 guard fires exactly at 100800 (one full week of 6 s counts), not at 100799.
    auto set_bits = []( uint8_t* buf, int pos, int len, uint32_t val )
    {
        for( int i = 0; i < len; ++i )
        {
            const int      bitpos = pos + i;
            const uint32_t bit    = ( val >> ( len - 1 - i ) ) & 1u;
            const uint8_t  mask   = static_cast<uint8_t>( 1u << ( 7 - bitpos % 8 ) );
            if( bit )
                buf[bitpos / 8] |= mask;
            else
                buf[bitpos / 8] = static_cast<uint8_t>( buf[bitpos / 8] & ~mask );
        }
    };
    auto bad = sf_from_hex( CTTC_PRN1_SF1 );
    set_bits( bad.data(), 30, 17, 100799u ); // max valid count -> accepted (NOT > 100799)
    REQUIRE( getbitu( bad.data(), 30, 17 ) == 100799u );
    REQUIRE_FALSE( getbitu( bad.data(), 30, 17 ) > 100799u );
    set_bits( bad.data(), 30, 17, 100800u ); // first impossible count -> dropped
    REQUIRE( getbitu( bad.data(), 30, 17 ) > 100799u );
    set_bits( bad.data(), 30, 17, 0x1FFFFu ); // all-ones 17-bit field (a typical garbage decode)
    REQUIRE( getbitu( bad.data(), 30, 17 ) > 100799u );
}
#endif

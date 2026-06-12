#pragma once
#include <array>
#include <cstdint>
#include "ephemeris.h"
#include "ionospheric.h"
#include "signal.h"
#include "types.h"

// GPS L1CA navigation constants - mirrors sdr.h / sdrinit.c initnavstruct() values
static constexpr int NAV_RATE    = 20;  // NAVRATE_L1CA:  epochs (ms) per nav bit
static constexpr int NAV_FLEN    = 300; // NAVFLEN_L1CA:  bits per subframe
static constexpr int NAV_ADDFLEN = 2;   // NAVADDFLEN_L1CA: D29*/D30* carry bits
static constexpr int NAV_PRELEN  = 8;   // NAVPRELEN_L1CA: preamble length (bits)
static constexpr int NAV_SYNCTH  = 50;  // NAVSYNCTH:     votes to declare bit sync
static constexpr int NAV_LOOP_MS = 10;  // LOOP_L1CA: loop filter interval (ms)

// Decoded GPS L1CA ephemeris (subframes 1, 2, 3).
// Field names and formulas follow GPS ICD IS-GPS-200 / GNSS-SDRLIB sdrnav_gps.c.
struct Gps_ephemeris : public Ephemeris
{
    int iodc = -1; // issue of data clock
    int iode = -1; // issue of data ephemeris (should match iodc)
    // group delay Tgd lives in the base Ephemeris::group_delay now.
}; // TODO: Figure out what the whole curve fit thing in the GPS ICD means.

// GPS broadcast ionosphere (Klobuchar) + leap-second offset, from subframe 4 page 18
// (SV ID 56). These are GPS-system-wide (every SV broadcasts the same set), so one decoded
// page serves the whole position fix. alpha/beta drive the Klobuchar slant-delay model;
// leap_seconds (dt_LS) is only needed to print UTC, not for the GPS-time PVT itself.
// struct Iono
// {
//     bool   valid        = false;
//     double alpha[4]     = {}; // amplitude coeffs (s, s/sc, s/sc^2, s/sc^3)
//     double beta[4]      = {}; // period coeffs    (s, s/sc, s/sc^2, s/sc^3)
//     int    leap_seconds = 0;  // dt_LS: GPS-UTC offset (s)
// };
// TODO: Above should follow the same format as Ephemeris and be a base class for the virtual iono() method in Nav_decoder, but
// this is simpler for now since only GPS has iono.

// GPS L1CA navigation bit synchronisation + subframe decoder.
// mirrors GNSS-SDRLIB sdrnavigation() + checksync() + checkbit() + findpreamble() +
// decode_l1ca() + decode_frame_l1ca() in sdrnav.c / sdrnav_gps.c
class Gps_l1ca_decoder : public Nav_decoder
{
public:
    Gps_l1ca_decoder( Satellite_id prn, const Signal_params& params );

    // Called every tracking epoch immediately after correlation.
    // prompt_i      = II_[0]:   current epoch raw prompt I
    // prompt_i_prev = old_I_[0]: previous epoch raw prompt I
    // mirrors: sdrnavigation(sdr, buffloc, cnt)
    void process( double prompt_i, double prompt_i_prev ) override;

    int          get_current_bit_index() const override;
    int          get_current_ms_tick() const override;
    const double ms_since_tow_update() const override;

    bool bit_sync_found() const override
    {
        return bit_sync_found_;
    }
    bool preamble_found() const
    {
        return preamble_found_;
    }
    bool frame_synced() const override
    {
        return preamble_found_;
    }
    bool ephemeris_valid() const
    {
        return eph_.valid;
    }

    const Ephemeris& ephemeris() const override
    {
        return eph_;
    }

    // Broadcast Klobuchar iono + leap seconds (subframe 4 page 18). System-wide, so the
    // Observation_engine can take it from whichever channel decoded the page first.
    bool iono_valid() const
    {
        return iono_.valid;
    }
    const Iono* iono() const override
    {
        return &iono_;
    }

    // True in the epoch where the loop filter should fire (prm2 mode post-nav-sync).
    // mirrors nav->swloop
    bool sw_loop() const override
    {
        return sw_loop_;
    }

private:
    // mirrors checksync() - detects bit boundaries via sign changes
    void check_sync( double prompt_i, double prompt_i_prev );
    // mirrors checkbit() - accumulates epochs into bit decisions
    void check_bit( double prompt_i );
    // mirrors findpreamble() for L1CA - checks frame buffer for preamble + parity
    bool find_preamble();
    // mirrors decode_l1ca() + decode_frame_l1ca() + decode_subfrm1/2/3()
    void decode_subframe();

    // Classify a subframe's HOW TOW against the continuity anchor (see classify_tow in the .cpp).
    // Decides whether a (parity-valid) subframe is trustworthy, the first/tentative one, or the
    // corrupted-HOW glitch that must be dropped - so a clean first subframe is still decoded.
    enum class Tow_status
    {
        Continuous,      // agrees with a confirmed anchor -> trustworthy
        Provisional,     // first anchor, or a re-bootstrap -> decode, but not yet confirmed
        CorruptRejected, // confirmed anchor but TOW disagrees -> the glitch; drop this subframe
    };
    Tow_status classify_tow( double tow );

    // GPS L1CA preamble bits in +/-1 notation - mirrors pre_l1ca in initnavstruct()
    static constexpr std::array<int, NAV_PRELEN> PREAMBLE = { 1, -1, -1, -1, 1, -1, 1, 1 };

    Satellite_id         satellite_id_;
    const Signal_params& params_;          // this signal's physics: name (logs), code_period_s (epoch), ...
    uint64_t             epoch_count_ = 0; // total epochs since construction (== cnt in GNSS-SDRLIB)

    // Bit synchronisation state - mirrors sdrnav_t flagsync/synci/bitsync/biti
    int                       bit_count_      = 0;     // biti = epoch_count_ % NAV_RATE
    int                       sync_index_     = 0;     // synci: which bit_count_ value marks end of a bit
    bool                      bit_sync_found_ = false; // flagsync
    std::array<int, NAV_RATE> bit_sync_votes_ {};      // bitsync[rate]: vote counts

    // Bit accumulation state - mirrors sdrnav_t bitIP/cnt/swreset/swsync/swloop
    double bit_accum_ = 0.0;   // bitIP: accumulated prompt I over current bit
    int    bit_cnt_   = 0;     // nav->cnt: epoch counter for loop filter timing
    bool   sw_sync_   = false; // swsync: true when a bit was decided this epoch
    bool   sw_loop_   = false; // swloop: true when loop filter should fire

    // Frame buffer: last NAV_FLEN + NAV_ADDFLEN bits as +/-1 (fbits/fbitsdec for L1CA)
    std::array<int, NAV_FLEN + NAV_ADDFLEN> frame_bits_ {};
    int                                     nav_bit_count_ = 0; // Tracks 0 to 299 bits in the current 6s subframe

    // Polarity: +1 or -1, set when preamble is found  (nav->polarity)
    int  polarity_       = 1;
    bool preamble_found_ = false; // flagtow
    int  sf_decoded_     = 0;     // bitmask: bits 0/1/2 = SF1/SF2/SF3 decoded

    Gps_ephemeris eph_current_; // most recently decoded subframe (may be invalid)
    Gps_ephemeris eph_;         // current ephemeris (valid only if eph_.valid == true and iodc/iode match)
    Iono      iono_;        // broadcast Klobuchar iono + leap seconds (SF4 page 18)
};

#pragma once
#include <array>
#include <cstdint>
#include <vector>
#include "ephemeris.h"
#include "ionospheric.h"
#include "signal.h"

// GPS L1C CNAV-2 navigation decoder. CNAV-2 is a 1800-symbol frame (18 s @ 100 sym/s, on the L1Cd data
// component) = subframe 1 (9-bit TOI, 52 symbols, LFSR/BCH-coded) + subframe 2 (600 bits, LDPC(1200,600)
// + CRC) + subframe 3 (274 bits, LDPC(548,274) + CRC), with subframes 2+3 block-interleaved (38x46).
// Ported from PocketSDR (BSD-2; python/sdr_nav.py decode_L1CD/sync_CNV2_frame/decode_CNV2) - see
// [[reference-repos]], [[gps-l1c-plan]].
//
// PHASE 1 (this file): frame SYNC only. Frame timing comes from correlating the received 52-symbol
// subframe 1 against the LFSR-coded TOI patterns for all 400 TOIs (NOT from the L1Co overlay), so a clean
// data lock is enough. Once synced, frame_synced() is true and the TOI counts up every 18 s. The LDPC
// decode of subframes 2/3 + the ephemeris parse are Phase 2/3 (eph_ stays invalid until then).
class Gps_l1c_decoder : public Nav_decoder
{
public:
    explicit Gps_l1c_decoder( Satellite_id prn );

    void             process( double prompt_i, double prompt_i_prev ) override;
    bool             bit_sync_found() const override;
    bool             sw_loop() const override;
    bool             frame_synced() const override;
    const Ephemeris& ephemeris() const override;
    const Iono*      iono() const override;

    int          get_current_bit_index() const override;
    int          get_current_ms_tick() const override;
    const double ms_since_tow_update() const override;

    // CNAV-2 frame layout (symbols).
    static constexpr int SF1_SYMS   = 52;                    // subframe 1 (9-bit TOI, LFSR/BCH-coded)
    static constexpr int FRAME_SYMS = 1800;                  // one CNAV-2 frame (18 s)
    static constexpr int WINDOW     = FRAME_SYMS + SF1_SYMS; // 1852: a frame + the next frame's SF1
    static constexpr int N_TOI      = 400;                   // TOI range 0..399 (frame number within the 2-hour cycle)
    // Frame-sync match tolerance: max symbol mismatches allowed per 52-symbol subframe-1 block. The
    // true alignment gives 0-2 mismatches (25%-power data, occasional bit error); a random alignment
    // floors at ~10-13, so 5 cleanly separates. Requiring BOTH ends (this frame + next) within tolerance
    // makes a false sync astronomically unlikely. (PocketSDR matches exactly because its full LDPC+CRC
    // decode confirms the sync; our Phase 1 has no LDPC yet, so the tolerant match stands on its own.)
    static constexpr int MATCH_TOL = 5;

private:
    // sync_CNV2_frame: does the WINDOW-symbol view match the subframe-1 patterns for `toi` (start) and
    // toi+1 (end)? Returns +1 normal, 0 reversed (180-deg polarity), -1 no match. Mirrors PocketSDR.
    int  match_frame( const uint8_t* window, int toi ) const;
    void build_sf1_templates(); // fill sf1_[toi][0..51] once
    // decode_CNV2: de-interleave subframes 2+3 of the 1800-symbol frame at `window`, LDPC-decode both,
    // CRC-check. Returns true on a fully valid frame (both CRCs pass). (PocketSDR decode_CNV2.)
    bool decode_frame( const uint8_t* window );
    // Parse CNAV-2 subframe-2 (600 LDPC-decoded bits) into eph_ (IS-GPS-800E Table 3.5-2 layout, CNAV
    // scales). Called only on a CRC-valid frame, so the result is trustworthy.
    void parse_sf2( const uint8_t* sf2_bits );

    Satellite_id satellite_id_;
    Ephemeris    eph_; // invalid until Phase 3 (LDPC + ephemeris parse)

    std::vector<uint8_t>                             syms_;   // hard data symbols (0/1), newest at back
    std::array<std::array<uint8_t, SF1_SYMS>, N_TOI> sf1_ {}; // subframe-1 template per TOI
    uint64_t                                         epoch_count_   = 0;
    bool                                             frame_synced_  = false;
    int                                              toi_           = 0; // current frame's TOI
    int                                              rev_           = 1; // polarity (1 normal, 0 reversed)
    uint64_t                                         epoch_at_sync_ = 0; // epoch when boundary last (re)checked
};

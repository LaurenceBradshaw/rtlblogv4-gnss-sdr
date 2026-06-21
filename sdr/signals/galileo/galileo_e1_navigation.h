#pragma once
#include <array>
#include <cstdint>
#include <map>
#include <vector>
#include "almanac.h"
#include "ephemeris.h"
#include "ionospheric.h"
#include "signal.h"

// Decoded Galileo I/NAV ephemeris (word types 1-4 = orbit + clock, word 5 = GST week/TOW).
// Field names/units follow the Galileo OS SIS-ICD (mirrors GPS's Gps_ephemeris layout).
struct Galileo_ephemeris : public Ephemeris
{
    int iod_nav = -1; // issue of data (words 1-4 must share it)
    int sisa    = 0;  // signal-in-space accuracy index
};

// Galileo E1-B I/NAV decoder.
//
// One symbol per tracking epoch (250 sym/s) arrives via process() as the E1-B data prompt
// (the pilot is tracked for carrier/code; this is the data-component prompt). Pipeline:
//   sync to the 10-symbol preamble (0101100000, recurs every 250 symbols) -> collect the
//   240 data symbols of each page part -> de-interleave (8x30) -> un-invert G2 -> soft
//   Viterbi (K=7, rate 1/2, G={171,133}oct, G2 inverted) -> 114 bits -> even/odd page
//   assembly -> CRC-24Q. CRC pass == Galileo nav decoded.
//
// Reference: Galileo OS SIS-ICD; GNSS-SDR galileo_telemetry_decoder_gs.cc, viterbi_decoder,
// galileo_inav_message.cc.
class Galileo_e1b_decoder : public Nav_decoder
{
public:
    Galileo_e1b_decoder( Satellite_id prn, const Signal_params& params );

    void process( double prompt_i, double prompt_i_prev ) override;
    bool bit_sync_found() const override
    {
        return false;
    } // not used for E1 loop control
    bool frame_synced() const override
    {
        return synced_; // I/NAV page-part preamble sync (see synced_)
    }
    bool sw_loop() const override
    {
        return false;
    }

    bool ephemeris_valid() const
    {
        return eph_.valid;
    }
    const Ephemeris& ephemeris() const override
    {
        return eph_;
    }

    const Iono* iono() const override
    {
        static const Iono dummy_iono {};
        return &dummy_iono; // TODO: no iono model in Galileo for now; revisit if that changes - currently GPS only
    }

    // Reduced almanac decoded from I/NAV word types 7-10 (coarse orbit for acquisition aiding). Keyed by PRN
    // within Galileo; the receiver re-keys by sv_key(constellation, prn) when it aggregates.
    const std::map<int, Almanac>& almanac() const override
    {
        return almanac_;
    }

    int get_current_bit_index() const override
    {
        return 0;
    }
    int get_current_ms_tick() const override
    {
        return 0;
    }
    // Milliseconds of SV transmit time since the TOW reference (the even-page-part preamble that
    // eph_.tow / TOW5 refers to). One symbol = one E1 code period (4 ms); tow_ref_symbol_ is the
    // symbol count at that preamble (set when word 5 is decoded). channel.cpp adds eph_.tow to this.
    const double ms_since_tow_update() const override
    {
        return static_cast<double>( symbol_count_ - tow_ref_symbol_ ) * params_.code_period_s * 1000.0;
    }

private:
    void decode_page_part();            // data_buf_ (240 soft symbols) -> de-interleave/Viterbi/CRC
    void extract_word( const int* jk ); // 128-bit I/NAV word -> ephemeris fields

    static constexpr int PREAMBLE_LEN   = 10;
    static constexpr int PAGE_PART_SYMS = 250; // 10 preamble + 240 data
    static constexpr int DATA_SYMS      = 240;
    static constexpr int PAGE_PART_BITS = 114; // Viterbi output bits per page part

    // I/NAV sync pattern 0101100000 in +/-1 (bit 1 -> +1, bit 0 -> -1).
    static constexpr std::array<int, PREAMBLE_LEN> PREAMBLE = { -1, 1, -1, 1, 1, -1, -1, -1, -1, -1 };

    Satellite_id         satellite_id_;
    const Signal_params& params_; // this signal's physics: name (logs), code_period_s, ...
    uint64_t             symbol_count_  = 0;
    uint64_t             tow_ref_symbol_ = 0; // symbol_count_ at the TOW5 reference (even preamble)

    // Preamble search (sign-based) until frame sync. Detections vote by page-part phase
    // (symbol_count % 250); the phase that reaches PREAMBLE_VOTE_TH is the frame boundary.
    static constexpr int            PREAMBLE_VOTE_TH = 4;
    std::array<int, PREAMBLE_LEN>   pre_win_ {};
    int                             pre_fill_ = 0;
    std::array<int, PAGE_PART_SYMS> pre_votes_ {}; // per-phase preamble-detection counts
    std::array<int, PAGE_PART_SYMS> pre_pol_ {};   // per-phase summed polarity (sign of corr)

    // Frame state once synced.
    bool               synced_         = false;
    int                frame_polarity_ = 1;
    int                frame_pos_      = 0; // position within the 250-symbol page part
    std::vector<float> data_buf_;           // 240 collected soft data symbols (polarity-corrected)

    // Page assembly + CRC stats.
    std::vector<int> page_even_; // 114 bits of the stored even page part
    bool             have_even_  = false;
    int              crc_ok_   = 0;
    int              crc_fail_ = 0;

    // Ephemeris assembly: words 1-4 (bits 0-3, IOD-gated) + word 5 GST (bit 4).
    Galileo_ephemeris eph_;
    int               eph_words_  = 0;
    int               eph_iodnav_ = -1;

    // Reduced-almanac assembly: the 3 SVs of one almanac are spread across I/NAV word types 7,8,9,10, all
    // sharing an IODa. Buffer each word's 128 bits + its IODa; once all four are present with a matching IODa,
    // assemble_almanac() decodes the 3 SVs into almanac_. (Mirrors RTKLIB decode_gal_inav_alm.)
    void store_almanac_word( int word_type, const int* jk ); // buffer a word 7-10, assemble when complete
    void assemble_almanac();                                 // decode the 3 SVs once 7-10 share an IODa

    std::array<std::array<int, 128>, 4> alm_words_ {};         // buffered bits, index 0..3 = word type 7..10
    std::array<int, 4>                  alm_ioda_ { -1, -1, -1, -1 };
    int                                 alm_decoded_ioda_ = -1; // last IODa already assembled (avoid re-decode)
    std::map<int, Almanac>              almanac_;                // PRN -> coarse orbit
};

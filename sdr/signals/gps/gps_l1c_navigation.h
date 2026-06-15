#pragma once
#include "ephemeris.h"
#include "ionospheric.h"
#include "signal.h"

// GPS L1C CNAV-2 navigation decoder - SKELETON (codes + tracking are the focus first; the message
// decoder is deferred). CNAV-2 is very different from LNAV: a 1800-symbol frame (18 s) of subframe 1
// (9-bit TOI, BCH) + subframe 2 (600 bits, LDPC(1200,600) + block interleaver) + subframe 3 (274 bits,
// LDPC + CRC). Decoding it needs an LDPC decoder - see [[gps-l1c-plan]] Phase 4.
//
// Until then this accepts epochs but never declares sync and exposes no ephemeris, so an L1C channel
// will acquire + track (Phase 2/3) but produce no nav data.
class Gps_l1c_decoder : public Nav_decoder
{
public:
    explicit Gps_l1c_decoder( Satellite_id prn );

    void             process( double prompt_i, double prompt_i_prev ) override;
    bool             bit_sync_found() const override;
    bool             sw_loop() const override;
    const Ephemeris& ephemeris() const override;
    const Iono*      iono() const override;

    int          get_current_bit_index() const override;
    int          get_current_ms_tick() const override;
    const double ms_since_tow_update() const override;

private:
    Satellite_id satellite_id_;
    Ephemeris    eph_; // always invalid - no decoding yet
};

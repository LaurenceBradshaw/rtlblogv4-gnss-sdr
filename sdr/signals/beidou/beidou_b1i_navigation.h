#pragma once
#include "ephemeris.h"
#include "ionospheric.h"
#include "signal.h"

// BeiDou B1I D1 navigation decoder - SKELETON (parameters/codes are implemented; this
// message decoder is deferred until B1I IQ data is available to verify against).
//
// D1 (MEO/IGSO satellites): 50 bps, each data bit spread by a 20-bit Neuman-Hofman
// secondary code; superframe of 5 subframes; each word protected by BCH(15,11,1) with
// interleaving - substantially different from GPS LNAV's XOR parity.
//
// TODO to make this real:
//   - NH(20) secondary-code wipeoff + 20 ms bit synchronisation.
//   - BCH(15,11,1) decode + de-interleave of each word.
//   - Subframe/superframe assembly + D1 ephemeris bit-field extraction (BeiDou SIS-ICD).
//   - GEO satellites transmit D2 (500 bps, no NH) - a separate decoder path.
class Beidou_b1i_decoder : public Nav_decoder
{
public:
    explicit Beidou_b1i_decoder( Satellite_id prn );

    void             process( double prompt_i, double prompt_i_prev ) override;
    bool             bit_sync_found() const override;
    bool             sw_loop() const override;
    const Ephemeris& ephemeris() const override;
    const Iono*  iono() const override;

    int          get_current_bit_index() const override;
    int          get_current_ms_tick() const override;
    const double ms_since_tow_update() const override;

private:
    Satellite_id satellite_id_;
    Ephemeris    eph_; // always invalid - no decoding yet
};

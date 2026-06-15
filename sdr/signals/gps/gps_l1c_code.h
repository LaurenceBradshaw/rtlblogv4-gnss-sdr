#pragma once
#include <vector>
#include "types.h"

namespace gps
{

// GPS L1C ranging code generators (IS-GPS-800), ported from GNSS-SDRLIB sdrcode.c.
//
// L1C has a data component (L1Cd, carries the CNAV-2 message) and a pilot component (L1Cp, data-free,
// tracked for a clean lock - like Galileo E1-C). Both primaries are 10230-chip WEIL codes at 1.023 Mcps
// (10 ms period): the Legendre sequence (length 10223) is combined with a per-PRN Weil index, then a
// fixed 7-chip sequence is inserted at a per-PRN index to reach 10230. The pilot additionally carries a
// per-PRN 1800-chip OVERLAY (L1Co) at one chip per 10 ms primary period (18 s period).
//
// Modulation is TMBOC(6,1,4/33) for the pilot and BOC(1,1) for the data; we use the BOC(1,1)
// approximation for both (same approach as Galileo's CBOC->BOC(1,1)) - each chip becomes two opposite
// half-chips, doubling the chip count. The overlay is NOT BOC-modulated (it is the secondary sequence
// applied once per primary period; tracking wipes it off after secondary sync).
//
// GPS L1C is defined for PRN 1..63 here (PRN >= 64 are QZSS/other and use a second S2 LFSR we omit).

// L1Cp - the PILOT (data-free) component: acquisition + carrier/code loops track this.
class L1cp_code
{
public:
    static constexpr double PERIOD_TIME   = 10e-3;
    static constexpr int    PRIMARY_CHIPS = 10230;
    static constexpr int    BOC_CHIPS     = 2 * PRIMARY_CHIPS; // 20460 half-chips, post-BOC(1,1)
    static constexpr int    OVERLAY_CHIPS = 1800;              // L1Co secondary, 1 per 10 ms epoch
    static constexpr int    MAX_PRN       = 63;

    static Complex_buf        generate( Satellite_id prn_id, uint32_t sampling_rate ); // BOC, sampled (acq)
    static std::vector<float> chips( Satellite_id prn_id );                            // BOC half-chips (tracking)
    static std::vector<float> primary_chips( Satellite_id prn_id );                    // raw 10230-chip Weil
    static std::vector<float> overlay_chips( Satellite_id prn_id );                    // 1800-chip L1Co secondary
};

// L1Cd - the DATA component (CNAV-2). Despread at the pilot's code phase to recover nav symbols.
class L1cd_code
{
public:
    static constexpr double PERIOD_TIME   = 10e-3;
    static constexpr int    PRIMARY_CHIPS = 10230;
    static constexpr int    BOC_CHIPS     = 2 * PRIMARY_CHIPS;
    static constexpr int    MAX_PRN       = 63;

    static Complex_buf        generate( Satellite_id prn_id, uint32_t sampling_rate );
    static std::vector<float> chips( Satellite_id prn_id );
    static std::vector<float> primary_chips( Satellite_id prn_id );
};

} // namespace gps

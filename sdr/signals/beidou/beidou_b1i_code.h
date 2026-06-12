#pragma once

#include <vector>
#include "types.h"

namespace beidou
{

// BeiDou B1I ranging code (CB1I) generator.
//   2046-chip Gold code at 2.046 Mcps -> 1 ms period. Two 11-stage LFSRs (G1, G2);
//   the per-SV code is G1 xor (two selected G2 taps), truncated to 2046 chips.
//   Algorithm/tap table follow the BeiDou SIS-ICD (cross-checked against
//   GNSS-SDRLIB gencode_B1IB2I). Mirrors the gps::L1ca_code interface.
//
// NOTE: This is the PRIMARY code only. B1I D1 navigation additionally overlays a
// 20-bit Neuman-Hofman secondary code at the 1 kHz primary-code rate (one NH chip per
// ms, repeating every 20 ms = one data bit) - that belongs in tracking/nav, not here.
class B1i_code
{
public:
    static constexpr double PERIOD_TIME  = 1e-3; // 1 ms
    static constexpr int    PERIOD_CHIPS = 2046;
    static constexpr int    MAX_PRN      = 37; // phase-assignment table size (BeiDou SIS-ICD)

    // Samples of one 2046-chip period (Q = 0) at the given sampling rate.
    static Complex_buf generate( Satellite_id prn_id, uint32_t sampling_rate );

    // Raw chip values (2046 +/-1 floats, one per chip).
    static std::vector<float> chips( Satellite_id prn_id );
};

} // namespace beidou

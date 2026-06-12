#pragma once

#include <vector>
#include "types.h"

namespace galileo
{

// Galileo E1-B ranging code generator.
//   4092-chip primary memory code at 1.023 Mcps (4 ms period), from the Galileo OS
//   SIS-ICD tables. Modulated with a BOC(1,1) subcarrier - each chip becomes two
//   opposite-sign half-chips - giving 8184 half-chips at an effective 2.046 Mcps.
//   This is the BOC(1,1) approximation of the true CBOC(6,1,11) (the BOC(6,1)
//   component, ~1/11 of the power, is omitted), matching GNSS-SDRLIB.
//
// E1-B is the DATA channel (carries the I/NAV message) and has no secondary code.
// (E1-C, the pilot, adds a 25-chip CS25 secondary - not generated here.)
//
// The 4092-chip primary codes are MEMORY codes (not LFSR-generated): the hex tables
// in galileo_e1b_hex.h are lifted verbatim from GNSS-SDRLIB (Galileo SIS-ICD).
class E1b_code
{
public:
    static constexpr double PERIOD_TIME   = 4e-3;              // 4 ms
    static constexpr int    PRIMARY_CHIPS = 4092;              // ranging chips, pre-BOC
    static constexpr int    BOC_CHIPS     = 2 * PRIMARY_CHIPS; // 8184 half-chips, post-BOC(1,1)
    static constexpr int    MAX_PRN       = 50;

    // One BOC(1,1)-modulated code period (Q = 0) sampled at sample_rate_hz.
    static Complex_buf generate( Satellite_id prn_id, uint32_t sampling_rate );

    // The BOC(1,1) half-chip sequence (BOC_CHIPS +/-1 floats) - the tracking replica.
    static std::vector<float> chips( Satellite_id prn_id );

    // The raw 4092-chip primary code (+/-1), before BOC - for inspection / testing.
    static std::vector<float> primary_chips( Satellite_id prn_id );
};

// Galileo E1-C ranging code generator (the PILOT channel).
//   Same construction as E1-B (4092-chip memory code -> BOC(1,1) -> 8184 half-chips),
//   but the E1-C primary is additionally overlaid with a 25-chip Neuman-Hofman secondary
//   code (CS25, the same sequence for every SV) at the 1 ms... no - at the primary-code
//   rate: one secondary chip per 4 ms primary period, so the secondary repeats every
//   25 * 4 ms = 100 ms.
//
// The pilot carries NO data, so once the (deterministic, known) secondary is wiped off in
// tracking, the FLL works and coherent integration can extend across epochs -> a clean
// carrier lock that the data-bearing E1-B cannot achieve. Tracking holds the secondary
// out of the replica (chips()/generate() return the PRIMARY only) and applies it per
// epoch after secondary sync; secondary_chips() exposes the 25-chip sequence for that.
class E1c_code
{
public:
    static constexpr double PERIOD_TIME     = 4e-3;
    static constexpr int    PRIMARY_CHIPS   = 4092;
    static constexpr int    BOC_CHIPS       = 2 * PRIMARY_CHIPS; // 8184 half-chips
    static constexpr int    SECONDARY_CHIPS = 25;                // CS25 overlay (100 ms period)
    static constexpr int    MAX_PRN         = 50;

    static Complex_buf        generate( Satellite_id prn_id, uint32_t sampling_rate );
    static std::vector<float> chips( Satellite_id prn_id );
    static std::vector<float> primary_chips( Satellite_id prn_id );

    // The 25-chip CS25 secondary (+/-1), identical for all SVs. One chip per 4 ms epoch.
    static std::vector<float> secondary_chips();
};

} // namespace galileo

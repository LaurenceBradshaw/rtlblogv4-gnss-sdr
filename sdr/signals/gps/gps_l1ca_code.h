#pragma once

#include <vector>
#include "types.h"

namespace gps
{

class L1ca_code
{
public:
    static constexpr double PERIOD_TIME  = 1e-3;
    static constexpr int    PERIOD_CHIPS = 1023;

    // generate samples of one 1 ms period at the given sampling rate
    static Complex_buf generate( Satellite_id prn_id, uint32_t sampling_rate );

    // return raw chip values (1023 +/-1 floats, one per chip)
    static std::vector<float> chips( Satellite_id prn_id );
};

} // namespace gps
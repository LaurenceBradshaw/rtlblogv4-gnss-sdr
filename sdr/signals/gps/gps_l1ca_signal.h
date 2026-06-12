#pragma once
#include "signal.h"

// GPS L1 C/A - the fully implemented reference signal.
//   1023-chip Gold codes, 1.023 Mcps, 1 ms period, BPSK, 50 bps LNAV (20 ms/bit).
// Code generation delegates to gps::L1ca_code; the decoder is Gps_l1ca_decoder.
class Gps_l1ca_signal : public Signal
{
public:
    const Signal_params& params() const override
    {
        return params_;
    }
    std::pair<int, int> sv_range() const override
    {
        return { 1, 32 };
    }

    Complex_buf                  code_samples( Satellite_id sv, double sample_rate_hz ) const override;
    std::vector<float>           code_chips( Satellite_id sv ) const override;
    std::unique_ptr<Nav_decoder> make_nav_decoder( Satellite_id sv ) const override;

private:
    static const Signal_params params_;
};

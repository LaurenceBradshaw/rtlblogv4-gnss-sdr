#pragma once
#include "signal.h"

// BeiDou B1I. The primary code generator (beidou::B1i_code) is implemented; the D1
// navigation decoder (Beidou_b1i_decoder) is still a skeleton, and tracking does not yet
// wipe the NH secondary code - so a B1I channel can acquire/track but not yet decode.
//
// NB: B1I is centred at 1561.098 MHz, OUTSIDE a 1575.42 MHz L1 capture - verifying the
// tracking/nav path needs a B1I IQ recording. The code generator is unit-tested
// independently (tests at the bottom of beidou_b1i_code.cpp).
class Beidou_b1i_signal : public Signal
{
public:
    const Signal_params& params() const override
    {
        return params_;
    }
    std::pair<int, int> sv_range() const override
    {
        return { 1, 37 };
    }

    Complex_buf                  code_samples( Satellite_id sv, double sample_rate_hz ) const override;
    std::vector<float>           code_chips( Satellite_id sv ) const override;
    std::unique_ptr<Nav_decoder> make_nav_decoder( Satellite_id sv ) const override;

private:
    static const Signal_params params_;
};

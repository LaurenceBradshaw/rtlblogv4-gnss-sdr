#pragma once
#include "signal.h"

// GPS L1C. Like Galileo E1, acquisition + the carrier/code loops track the data-free PILOT (L1Cp):
// make_tracker() returns a Pilot_tracker. code_samples/code_chips return the L1Cp primary BOC(1,1)
// code; secondary_code(sv) returns the per-PRN 1800-chip L1Co overlay; data_code_chips(sv) returns the
// L1Cd data code for the nav-symbol prompt. The nav decoder (CNAV-2) is a skeleton for now.
class Gps_l1c_signal : public Signal
{
public:
    const Signal_params& params() const override
    {
        return params_;
    }
    std::pair<int, int> sv_range() const override
    {
        return { 1, 32 }; // same operational GPS PRNs as L1 C/A (codes defined to 63)
    }

    Complex_buf                  code_samples( Satellite_id sv, double sample_rate_hz ) const override;
    std::vector<float>           code_chips( Satellite_id sv ) const override;
    std::unique_ptr<Nav_decoder> make_nav_decoder( Satellite_id sv ) const override;
    std::vector<float>           secondary_code( Satellite_id sv ) const override;
    std::vector<float>           data_code_chips( Satellite_id sv ) const override;
    std::unique_ptr<Tracker>     make_tracker( Satellite_id sv, double sample_rate_hz ) const override;

private:
    static const Signal_params params_;
};

#pragma once
#include "signal.h"

// Galileo E1. Acquisition and the carrier/code loops track the E1-C PILOT (data-free), the
// way GNSS-SDR does (Galileo_E1 DLL_PLL VEML, track_pilot=true): make_tracker() returns a
// Pilot_tracker that pulls in with Costas+FLL, syncs + wipes the CS25 secondary, then runs a
// pure PLL on the pilot. code_samples/code_chips return the E1-C primary BOC(1,1) code;
// secondary_code() returns CS25. (E1-B data demod for nav comes via a data-component prompt.)
class Galileo_e1_signal : public Signal
{
public:
    const Signal_params& params() const override
    {
        return params_;
    }
    std::pair<int, int> sv_range() const override
    {
        return { 1, 36 };
    }

    Complex_buf                  code_samples( Satellite_id sv, double sample_rate_hz ) const override;
    std::vector<float>           code_chips( Satellite_id sv ) const override;
    std::unique_ptr<Nav_decoder> make_nav_decoder( Satellite_id sv ) const override;
    std::vector<float>           secondary_code() const override;
    std::vector<float>           data_code_chips( Satellite_id sv ) const override;
    std::unique_ptr<Tracker>     make_tracker( Satellite_id sv, double sample_rate_hz ) const override;

private:
    static const Signal_params params_;
};

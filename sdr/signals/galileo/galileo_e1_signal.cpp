#include "galileo_e1_signal.h"
#include "galileo_e1_code.h"
#include "galileo_e1_navigation.h"
#include "pilot_tracker.h"

const Signal_params Galileo_e1_signal::params_ = {
    /* constellation     */ Constellation::Galileo,
    /* band              */ Band::E1,
    /* code              */ Code::B,
    /* name              */ "Galileo E1",
    /* carrier_freq_hz   */ 1575.42e6,                        // same L1 band as GPS - reachable from an L1 capture
    /* chip_rate_hz      */ 1.023e6,                          // ranging chip rate (BOC subcarrier handled in tracking)
    /* code_length_chips */ galileo::E1c_code::PRIMARY_CHIPS, // 4092 ranging chips (single source of truth)
    /* code_period_s     */ galileo::E1c_code::PERIOD_TIME,   // 4 ms
    /* nav_bit_ms        */ 4,                                // I/NAV 250 sym/s -> 4 ms/symbol (E1-B data)
    /* modulation        */ Modulation::Boc11,
    /* acq_integrations  */ 4, // 4 ms code -> ~6 dB/epoch, fewer epochs needed
    /* acq_fft_factor    */ 1, // single-period circular FFT (4x cheaper/call)
    /* loop_bw_wide      */ { 5.0, 30.0, 200.0 }, // {dll,pll,fll} Costas+strong-FLL pull-in (atan FLL)
    /* loop_bw_narrow    */ { 3.0, 20.0, 50.0 },  // tight pure-PLL steady state on the data-free pilot
    /* frame_sync_timeout_s */ 45.0, // I/NAV page-vote sync is slower (can take 30 s+); generous margin
};

Complex_buf Galileo_e1_signal::code_samples( Satellite_id sv, double sample_rate_hz ) const
{
    return galileo::E1c_code::generate( sv, static_cast<uint32_t>( sample_rate_hz ) ); // pilot
}

std::vector<float> Galileo_e1_signal::code_chips( Satellite_id sv ) const
{
    return galileo::E1c_code::chips( sv ); // E1-C pilot primary BOC code
}

std::unique_ptr<Nav_decoder> Galileo_e1_signal::make_nav_decoder( Satellite_id sv ) const
{
    return std::make_unique<Galileo_e1b_decoder>( sv, params_ );
}

std::vector<float> Galileo_e1_signal::secondary_code( Satellite_id /*sv*/ ) const
{
    return galileo::E1c_code::secondary_chips(); // CS25 (25 chips, 1 per 4 ms epoch) - same for all SVs
}

std::vector<float> Galileo_e1_signal::data_code_chips( Satellite_id sv ) const
{
    return galileo::E1b_code::chips( sv ); // E1-B data code (nav symbols), at the pilot's phase
}

std::unique_ptr<Tracker> Galileo_e1_signal::make_tracker( Satellite_id sv, double sample_rate_hz ) const
{
    return std::make_unique<Pilot_tracker>(
        code_chips( sv ), sample_rate_hz, params_, secondary_code( sv ), data_code_chips( sv )
    );
}

#include "beidou_b1i_signal.h"
#include "beidou_b1i_code.h"
#include "beidou_b1i_navigation.h"

const Signal_params Beidou_b1i_signal::params_ = {
    /* constellation     */ Constellation::Beidou,
    /* band              */ Band::B1,
    /* code              */ Code::I,
    /* name              */ "BeiDou B1I",
    /* carrier_freq_hz   */ 1561.098e6, // B1I band - NOT in an L1 (1575.42 MHz) capture
    /* chip_rate_hz      */ 2.046e6,
    /* code_length_chips */ beidou::B1i_code::PERIOD_CHIPS, // single source of truth (beidou_b1i_code.h)
    /* code_period_s     */ beidou::B1i_code::PERIOD_TIME,
    /* nav_bit_ms        */ 20,
    /* modulation        */ Modulation::Bpsk,
    /* acq_integrations  */ 10,
    /* acq_fft_factor    */ 2,
    /* loop_bw_wide      */ { 4.0, 40.0, 25.0 }, // {dll,pll,fll} pull-in (cross/dot FLL)
    /* loop_bw_narrow    */ { 2.0, 25.0, 10.0 }, // steady state after frame sync
    /* frame_sync_timeout_s */ 30.0, // B1I D1 subframes (6 s) like GPS; no decoder yet, generous default
};

Complex_buf Beidou_b1i_signal::code_samples( Satellite_id sv, double sample_rate_hz ) const
{
    return beidou::B1i_code::generate( sv, static_cast<uint32_t>( sample_rate_hz ) );
}

std::vector<float> Beidou_b1i_signal::code_chips( Satellite_id sv ) const
{
    return beidou::B1i_code::chips( sv );
}

std::unique_ptr<Nav_decoder> Beidou_b1i_signal::make_nav_decoder( Satellite_id sv ) const
{
    return std::make_unique<Beidou_b1i_decoder>( sv );
}

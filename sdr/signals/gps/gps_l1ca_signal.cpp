#include "gps_l1ca_signal.h"
#include <cstdint>
#include "gps_l1ca_code.h"
#include "gps_l1ca_navigation.h"

const Signal_params Gps_l1ca_signal::params_ = {
    /* constellation     */ Constellation::Gps,
    /* band              */ Band::L1,
    /* code              */ Code::CA,
    /* name              */ "GPS L1 C/A",
    /* carrier_freq_hz   */ 1575.42e6,
    /* chip_rate_hz      */ 1.023e6,
    /* code_length_chips */ gps::L1ca_code::PERIOD_CHIPS, // single source of truth (gps_l1ca_code.h)
    /* code_period_s     */ gps::L1ca_code::PERIOD_TIME,
    /* nav_bit_ms        */ 20,
    /* modulation        */ Modulation::Bpsk,
    /* acq_integrations  */ 10,
    /* acq_fft_factor    */ 2,
    /* loop_bw_wide      */ { 4.0, 40.0, 25.0 }, // {dll,pll,fll} pull-in (cross/dot FLL)
    /* loop_bw_narrow    */ { 2.0, 25.0, 10.0 }, // steady state after frame sync
    /* frame_sync_timeout_s */ 25.0,             // subframe preamble lands by ~13 s; margin for a re-sync
};

Complex_buf Gps_l1ca_signal::code_samples( Satellite_id sv, double sample_rate_hz ) const
{
    return gps::L1ca_code::generate( sv, static_cast<uint32_t>( sample_rate_hz ) );
}

std::vector<float> Gps_l1ca_signal::code_chips( Satellite_id sv ) const
{
    return gps::L1ca_code::chips( sv );
}

std::unique_ptr<Nav_decoder> Gps_l1ca_signal::make_nav_decoder( Satellite_id sv ) const
{
    return std::make_unique<Gps_l1ca_decoder>( sv, params_ );
}

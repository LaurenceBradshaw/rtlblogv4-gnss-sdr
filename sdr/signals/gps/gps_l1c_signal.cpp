#include "gps_l1c_signal.h"
#include "gps_l1c_code.h"
#include "gps_l1c_navigation.h"
#include "pilot_tracker.h"

const Signal_params Gps_l1c_signal::params_ = {
    /* constellation     */ Constellation::Gps,
    /* band              */ Band::L1,
    /* code              */ Code::Cd,
    /* name              */ "GPS L1C",
    /* carrier_freq_hz   */ 1575.42e6,
    /* chip_rate_hz      */ 1.023e6,                       // ranging chip rate (BOC subcarrier handled in tracking)
    /* code_length_chips */ gps::L1cp_code::PRIMARY_CHIPS, // 10230 ranging chips
    /* code_period_s     */ gps::L1cp_code::PERIOD_TIME,   // 10 ms
    /* nav_bit_ms        */ 10,                            // CNAV-2: 100 sym/s -> 10 ms/symbol (= one code period)
    /* modulation        */ Modulation::Boc11,            // BOC(1,1) approximation of TMBOC(6,1,4/33)
    /* acq_integrations  */ 3,  // 10 ms code -> ~10 dB/epoch, few epochs needed
    /* acq_fft_factor    */ 1,  // single-period circular FFT (the 10 ms code is long; halves cost)
    // Tuned for the 10 ms epoch (gnss-sdr has no L1C; values bracketed from its E1 4 ms / L2C 20 ms sets
    // and swept on the Skydel sim). The wide FLL is the key knob: 200 (B*T~2.0) made pull-in intermittently
    // diverge (cos2phi stuck ~0); 30 (B*T~0.3) was too cold to pull in at all; 100 (B*T~1.0) locks 100% of
    // the time, deterministically, cos2phi ~0.92-0.97. The lower narrow set tightens the steady lock.
    /* loop_bw_wide      */ { 3.0, 25.0, 100.0 }, // {dll,pll,fll} pull-in (atan FLL)
    /* loop_bw_narrow    */ { 1.5, 12.0, 20.0 },  // steady state
    /* frame_sync_timeout_s */ 600.0, // no CNAV-2 decode yet -> never frame-syncs; don't drop the channel
};

Complex_buf Gps_l1c_signal::code_samples( Satellite_id sv, double sample_rate_hz ) const
{
    return gps::L1cp_code::generate( sv, static_cast<uint32_t>( sample_rate_hz ) ); // pilot
}

std::vector<float> Gps_l1c_signal::code_chips( Satellite_id sv ) const
{
    return gps::L1cp_code::chips( sv ); // L1Cp pilot primary BOC code
}

std::unique_ptr<Nav_decoder> Gps_l1c_signal::make_nav_decoder( Satellite_id sv ) const
{
    return std::make_unique<Gps_l1c_decoder>( sv );
}

std::vector<float> Gps_l1c_signal::secondary_code( Satellite_id sv ) const
{
    return gps::L1cp_code::overlay_chips( sv ); // per-PRN 1800-chip L1Co overlay (1 per 10 ms epoch)
}

std::vector<float> Gps_l1c_signal::data_code_chips( Satellite_id sv ) const
{
    return gps::L1cd_code::chips( sv ); // L1Cd data code (nav symbols), at the pilot's phase
}

std::unique_ptr<Tracker> Gps_l1c_signal::make_tracker( Satellite_id sv, double sample_rate_hz ) const
{
    return std::make_unique<Pilot_tracker>(
        code_chips( sv ), sample_rate_hz, params_, secondary_code( sv ), data_code_chips( sv )
    );
}

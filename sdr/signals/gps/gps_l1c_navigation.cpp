#include "gps_l1c_navigation.h"

Gps_l1c_decoder::Gps_l1c_decoder( Satellite_id prn )
    : satellite_id_( prn )
{
}

// SKELETON: accepts epochs, performs no CNAV-2 decoding yet, never declares sync.
void Gps_l1c_decoder::process( double /*prompt_i*/, double /*prompt_i_prev*/ ) {}

bool Gps_l1c_decoder::bit_sync_found() const
{
    return false;
}

bool Gps_l1c_decoder::sw_loop() const
{
    return false;
}

const Ephemeris& Gps_l1c_decoder::ephemeris() const
{
    return eph_;
}

const Iono* Gps_l1c_decoder::iono() const
{
    static const Iono dummy_iono {};
    return &dummy_iono; // no L1C iono decode yet
}

int Gps_l1c_decoder::get_current_bit_index() const
{
    return 0;
}

int Gps_l1c_decoder::get_current_ms_tick() const
{
    return 0;
}

const double Gps_l1c_decoder::ms_since_tow_update() const
{
    return 0.0;
}

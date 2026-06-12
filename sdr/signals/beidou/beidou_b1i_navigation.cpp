#include "beidou_b1i_navigation.h"

Beidou_b1i_decoder::Beidou_b1i_decoder( Satellite_id prn )
    : satellite_id_( prn )
{
}

// SKELETON: accepts epochs but performs no decoding yet, and never declares sync - so a
// B1I channel will acquire and track (once tracking handles the NH code) but produce no
// nav data. See the header TODO for the real D1 decode path.
void Beidou_b1i_decoder::process( double /*prompt_i*/, double /*prompt_i_prev*/ ) {}

bool Beidou_b1i_decoder::bit_sync_found() const
{
    return false;
}

bool Beidou_b1i_decoder::sw_loop() const
{
    return false;
}

const Ephemeris& Beidou_b1i_decoder::ephemeris() const
{
    return eph_;
}

const Iono* Beidou_b1i_decoder::iono() const
{
    static const Iono dummy_iono {};
    return &dummy_iono; // TODO: no iono model in BeiDou for now; revisit if that changes - currently GPS only
}

int Beidou_b1i_decoder::get_current_bit_index() const
{
    return 0;
}

int Beidou_b1i_decoder::get_current_ms_tick() const
{
    return 0;
}

const double Beidou_b1i_decoder::ms_since_tow_update() const
{
    return 0.0;
}

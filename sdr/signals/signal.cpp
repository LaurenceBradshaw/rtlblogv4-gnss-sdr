#include "signal.h"
#include <stdexcept>
#include "costas_tracker.h"
#include "beidou_b1i_signal.h"
#include "galileo_e1_signal.h"
#include "gps_l1ca_signal.h"

// Default tracker: data-bearing BPSK strategy (Costas + FLL, 3-tap E/P/L). Signals that
// track a pilot / BOC component override this (Galileo E1 -> Pilot_tracker).
std::unique_ptr<Tracker> Signal::make_tracker( Satellite_id sv, double sample_rate_hz ) const
{
    return std::make_unique<Costas_tracker>( code_chips( sv ), sample_rate_hz, params(), secondary_code() );
}

std::unique_ptr<Signal> make_signal( Constellation constellation, Band band, Code code )
{
    // Only the triples implemented below are valid; everything else throws. New signals are
    // registered by adding a case (and its Signal subclass).
    if( constellation == Constellation::Gps && band == Band::L1 && code == Code::CA )
        return std::make_unique<Gps_l1ca_signal>();
    if( constellation == Constellation::Galileo && band == Band::E1 && code == Code::B )
        return std::make_unique<Galileo_e1_signal>();
    if( constellation == Constellation::Beidou && band == Band::B1 && code == Code::I )
        return std::make_unique<Beidou_b1i_signal>();

    throw std::invalid_argument( "make_signal: unsupported (constellation, band, code) combination" );
}

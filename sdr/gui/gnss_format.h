#pragma once
#include <QColor>
#include "constellations.h" // Constellation
#include "signal.h"         // Code

// Shared GUI formatting for GNSS constellations - one home for the small maps that several widgets
// (satellite list/row, detail window, sky plot, ...) would otherwise each copy.
namespace gui_format
{
// One-letter prefix, GNSS convention: GPS=G, Galileo=E, BeiDou=C.
inline const char* constellation_prefix( Constellation c )
{
    switch( c )
    {
    case Constellation::Gps:
        return "G";
    case Constellation::Galileo:
        return "E";
    case Constellation::Beidou:
        return "C";
    default:
        return "?";
    }
}

// Full name for headings/labels (the core map; re-exported here for convenience).
using ::constellation_name;

// Short label for a signal component, to disambiguate the same SV tracked on >1 code (e.g. GPS L1CA vs
// L1C). One physical satellite can appear as two rows / bars / graphs, one per code.
// CAVEAT - NEEDS BAND when more bands land: Code is unique PER CONSTELLATION today (GPS {CA,Cd}, Galileo
// {B}, BeiDou {I}), so the GUI keys on (constellation, prn, code) and labels by Code alone. But the same
// code component recurs across bands - BeiDou "I" on B1/B2/B3, Galileo "B"/"C" on E1/E6 - so once a
// constellation has the same Code on two bands, BOTH this label (B1I vs B2I) and the keying
// (Satellite_list_widget::key_of, Receiver::history_key, Configured_satellite) must also take Band.
// Channel_snapshot already carries `band`; Configured_satellite would need it added.
inline const char* code_label( Code c )
{
    switch( c )
    {
    case Code::CA:
        return "L1CA";
    case Code::Cd:
        return "L1C";
    case Code::B:
        return "E1";
    case Code::I:
        return "B1I";
    }
    return "?";
}

// Plot/legend colour per constellation (sky-plot dots, etc.).
inline QColor constellation_color( Constellation c )
{
    switch( c )
    {
    case Constellation::Gps:
        return QColor( 80, 200, 120 ); // green
    case Constellation::Galileo:
        return QColor( 90, 150, 240 ); // blue
    case Constellation::Beidou:
        return QColor( 235, 160, 60 ); // orange
    default:
        return QColor( 180, 180, 180 );
    }
}
} // namespace gui_format

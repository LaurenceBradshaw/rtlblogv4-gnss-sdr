#pragma once
#include <QColor>
#include "constellations.h" // Constellation

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

#include "sky_plot_widget.h"
#include <QPainter>
#include <algorithm>
#include <cmath>
#include "gnss_format.h"

namespace
{
using gui_format::constellation_color;
using gui_format::constellation_prefix;

// Polar mapping: elevation 90 deg -> centre, 0 deg -> rim; azimuth clockwise from north (up).
QPointF sky_point( double az_deg, double el_deg, QPointF centre, double radius )
{
    const double r   = ( 90.0 - el_deg ) / 90.0 * radius;
    const double rad = az_deg * M_PI / 180.0;
    return { centre.x() + r * std::sin( rad ), centre.y() - r * std::cos( rad ) };
}
} // namespace

Sky_plot_widget::Sky_plot_widget( QWidget* parent )
    : Gui_panel( parent )
{
    setMinimumSize( 320, 320 );
}

void Sky_plot_widget::update_frame( const Gui_frame& frame )
{
    sats_ = frame.sky;
    update();
}

void Sky_plot_widget::paintEvent( QPaintEvent* )
{
    QPainter p( this );
    p.setRenderHint( QPainter::Antialiasing );
    p.fillRect( rect(), QColor( 20, 20, 20 ) );

    const QPointF centre( width() / 2.0, height() / 2.0 );
    const double  radius = std::min( width(), height() ) / 2.0 - 24.0;

    // Elevation rings (0 / 30 / 60 deg) and the zenith.
    p.setPen( QColor( 90, 90, 90 ) );
    for( int el = 0; el <= 60; el += 30 )
    {
        const double r = ( 90.0 - el ) / 90.0 * radius;
        p.drawEllipse( centre, r, r );
    }

    // Cardinal radials + labels.
    const char* const labels[] = { "N", "E", "S", "W" };
    for( int i = 0; i < 4; ++i )
    {
        const double  az  = i * 90.0;
        const QPointF rim = sky_point( az, 0.0, centre, radius );
        p.setPen( QColor( 90, 90, 90 ) );
        p.drawLine( centre, rim );
        p.setPen( Qt::white );
        const QPointF lbl = sky_point( az, -6.0, centre, radius );
        p.drawText( QRectF( lbl.x() - 8, lbl.y() - 8, 16, 16 ), Qt::AlignCenter, labels[i] );
    }

    // Satellites above the horizon.
    p.setFont( QFont( QStringLiteral( "monospace" ), 8 ) );
    for( const Sky_satellite& s : sats_ )
    {
        if( s.elevation_deg < 0.0 )
        {
            continue;
        }
        const QPointF pt  = sky_point( s.azimuth_deg, s.elevation_deg, centre, radius );
        const QColor  col = constellation_color( s.constellation );

        p.setPen( QPen( col, 2 ) );
        p.setBrush( s.has_lock ? QBrush( col ) : Qt::NoBrush ); // filled if locked, hollow if not
        p.drawEllipse( pt, 7.0, 7.0 );

        p.setPen( Qt::white );
        const QString tag = QString::asprintf( "%s%02d", constellation_prefix( s.constellation ), s.satellite_id );
        p.drawText( QRectF( pt.x() + 9, pt.y() - 8, 44, 16 ), Qt::AlignLeft | Qt::AlignVCenter, tag );
    }

    if( sats_.empty() )
    {
        p.setPen( Qt::gray );
        p.drawText( rect(), Qt::AlignCenter, QStringLiteral( "No satellites (needs a fix)" ) );
    }
}

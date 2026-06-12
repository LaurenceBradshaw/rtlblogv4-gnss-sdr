#include "iq_heatmap_widget.h"
#include <algorithm>
#include <cmath>
#include <QImage>
#include <QPainter>
#include "plasma.h"

Iq_heatmap_widget::Iq_heatmap_widget( QWidget* parent )
    : QWidget( parent )
{
    setMinimumSize( 280, 280 );
}

void Iq_heatmap_widget::set_snapshot( const Tracking_history::Iq_snapshot& snapshot )
{
    snap_ = snapshot;
    update();
}

void Iq_heatmap_widget::paintEvent( QPaintEvent* )
{
    QPainter p( this );
    p.fillRect( rect(), Qt::black ); // black background, matching the Python plot

    constexpr int N    = Tracking_history::IQ_GRID;
    const int     side = std::min( width(), height() );
    const QRect   box( ( width() - side ) / 2, ( height() - side ) / 2, side, side );

    if( snap_.has_data )
    {
        float max_heat = 0.0f;
        for( float h : snap_.heat )
        {
            max_heat = std::max( max_heat, h );
        }
        const double lmax = std::log( 1.0 + max_heat );

        // Build the grid as an image (log-scaled heat -> plasma), then smooth-scale it to fill the box.
        // Low-heat cells fade toward BLACK (the lowest ~12% of intensity) so empty space reads black like
        // the Python plot, rather than plasma's dark-purple floor.
        QImage img( N, N, QImage::Format_RGB32 );
        for( int q = 0; q < N; ++q )
        {
            for( int i = 0; i < N; ++i )
            {
                const float  h = snap_.heat[q * N + i];
                const double t = ( lmax > 0.0 ) ? std::log( 1.0 + h ) / lmax : 0.0;
                QColor       c = gui_format::plasma( t );
                const double b = std::clamp( t / 0.12, 0.0, 1.0 ); // black floor for empty/faint cells
                c.setRgb(
                    static_cast<int>( c.red() * b ), static_cast<int>( c.green() * b ), static_cast<int>( c.blue() * b )
                );
                // +Q points up on screen, so flip the row.
                img.setPixelColor( i, N - 1 - q, c );
            }
        }
        p.drawImage( box, img.scaled( side, side, Qt::IgnoreAspectRatio, Qt::SmoothTransformation ) );
    }
    else
    {
        p.setPen( Qt::gray );
        p.drawText( rect(), Qt::AlignCenter, QStringLiteral( "no data (not tracking)" ) );
    }

    // Axes through the origin + labels.
    p.setPen( QColor( 150, 150, 150, 150 ) );
    p.drawLine( box.center().x(), box.top(), box.center().x(), box.bottom() );
    p.drawLine( box.left(), box.center().y(), box.right(), box.center().y() );
    p.setPen( Qt::white );
    p.drawText( box.adjusted( 4, 4, -4, -4 ), Qt::AlignRight | Qt::AlignVCenter, QStringLiteral( "I" ) );
    p.drawText( box.adjusted( 4, 4, -4, -4 ), Qt::AlignHCenter | Qt::AlignTop, QStringLiteral( "Q" ) );
}

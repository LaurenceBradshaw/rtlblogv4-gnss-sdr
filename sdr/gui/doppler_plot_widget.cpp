#include "doppler_plot_widget.h"
#include <QPainter>
#include <QPolygonF>
#include <algorithm>

Doppler_plot_widget::Doppler_plot_widget( QWidget* parent )
    : QWidget( parent )
{
    setMinimumSize( 360, 220 );
}

void Doppler_plot_widget::set_series( const Tracking_history::Doppler_series& series )
{
    series_ = series;
    update();
}

void Doppler_plot_widget::paintEvent( QPaintEvent* )
{
    QPainter p( this );
    p.setRenderHint( QPainter::Antialiasing );
    p.fillRect( rect(), Qt::black );

    constexpr int margin_l = 56, margin_b = 22, margin_t = 8, margin_r = 10;
    const QRectF  plot( margin_l, margin_t, width() - margin_l - margin_r, height() - margin_t - margin_b );

    if( series_.t_s.size() < 2 )
    {
        p.setPen( Qt::gray );
        p.drawText( rect(), Qt::AlignCenter, QStringLiteral( "no data (not tracking)" ) );
        return;
    }

    // Auto-scale: x to the time span, y to the value range (padded). Guard against a flat line.
    const double t0     = series_.t_s.front();
    const double t1     = series_.t_s.back();
    auto [lo_it, hi_it] = std::minmax_element( series_.hz.begin(), series_.hz.end() );
    double ymin = *lo_it, ymax = *hi_it;
    double pad = std::max( 1.0, ( ymax - ymin ) * 0.1 );
    ymin -= pad;
    ymax += pad;
    const double tspan = std::max( 1e-6, t1 - t0 );
    const double yspan = std::max( 1e-6, ymax - ymin );

    const auto px = [&]( double t ) { return plot.left() + ( t - t0 ) / tspan * plot.width(); };
    const auto py = [&]( double hz ) { return plot.bottom() - ( hz - ymin ) / yspan * plot.height(); };

    // Grid + axis labels (3 y ticks, start/end x ticks).
    p.setPen( QColor( 70, 70, 70 ) );
    p.drawRect( plot );
    p.setPen( QColor( 160, 160, 160 ) );
    for( int i = 0; i <= 2; ++i )
    {
        const double hz = ymin + yspan * i / 2.0;
        const double y  = py( hz );
        p.drawText( QRectF( 0, y - 8, margin_l - 4, 16 ), Qt::AlignRight | Qt::AlignVCenter, QString::number( hz, 'f', 0 ) );
    }
    p.drawText(
        QRectF( plot.left(), plot.bottom() + 4, plot.width(), margin_b - 4 ),
        Qt::AlignLeft,
        QString::number( t0, 'f', 1 ) + " s"
    );
    p.drawText(
        QRectF( plot.left(), plot.bottom() + 4, plot.width(), margin_b - 4 ),
        Qt::AlignRight,
        QString::number( t1, 'f', 1 ) + " s"
    );

    // The trace.
    QPolygonF poly;
    poly.reserve( static_cast<int>( series_.t_s.size() ) );
    for( size_t i = 0; i < series_.t_s.size(); ++i )
    {
        poly << QPointF( px( series_.t_s[i] ), py( series_.hz[i] ) );
    }
    p.setPen( QPen( QColor( 0, 191, 255 ), 1.2 ) ); // deepskyblue
    p.drawPolyline( poly );
}

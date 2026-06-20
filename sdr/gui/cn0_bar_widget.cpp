#include "cn0_bar_widget.h"
#include <QHelpEvent>
#include <QPainter>
#include <QToolTip>
#include <algorithm>

namespace
{
constexpr double CN0_MAX      = 55.0; // dB-Hz top of the scale (typical strong-signal ceiling)
constexpr double CN0_STEP     = 10.0; // gridline spacing
constexpr int    MARGIN_L     = 34;   // room for the Y-axis labels
constexpr int    MARGIN_B     = 32;   // room for the two-line SV + code labels
constexpr int    MARGIN_T     = 8;
constexpr int    MARGIN_R     = 8;
constexpr double BAR_GAP_FRAC = 0.25; // fraction of each slot left as gap
} // namespace

Cn0_bar_widget::Cn0_bar_widget( QWidget* parent )
    : QWidget( parent )
{
    setMinimumWidth( 180 );
}

void Cn0_bar_widget::set_bars( std::vector<Bar> bars )
{
    bars_ = std::move( bars );
    update();
}

void Cn0_bar_widget::paintEvent( QPaintEvent* )
{
    QPainter p( this );
    p.fillRect( rect(), QColor( 25, 25, 25 ) );

    const QRectF plot(
        MARGIN_L, MARGIN_T, std::max( 0, width() - MARGIN_L - MARGIN_R ), std::max( 0, height() - MARGIN_T - MARGIN_B )
    );

    // Y-axis gridlines + C/N0 labels.
    p.setPen( QColor( 120, 120, 120 ) );
    QFont small = p.font();
    small.setPointSizeF( std::max( 7.0, small.pointSizeF() - 1.0 ) );
    p.setFont( small );
    for( double v = 0.0; v <= CN0_MAX + 0.1; v += CN0_STEP )
    {
        const double y = plot.bottom() - ( v / CN0_MAX ) * plot.height();
        p.setPen( QColor( 55, 55, 55 ) );
        p.drawLine( QPointF( plot.left(), y ), QPointF( plot.right(), y ) );
        p.setPen( QColor( 150, 150, 150 ) );
        p.drawText( QRectF( 0, y - 7, MARGIN_L - 4, 14 ), Qt::AlignRight | Qt::AlignVCenter, QString::number( int( v ) ) );
    }

    bar_rects_.assign( bars_.size(), QRectF() );
    if( bars_.empty() )
    {
        p.setPen( QColor( 150, 150, 150 ) );
        p.drawText( plot, Qt::AlignCenter, QStringLiteral( "No tracking satellites" ) );
        return;
    }

    const double slot     = plot.width() / static_cast<double>( bars_.size() );
    const double bar_w    = slot * ( 1.0 - BAR_GAP_FRAC );
    const double gap_half = ( slot - bar_w ) / 2.0;

    for( size_t i = 0; i < bars_.size(); ++i )
    {
        const Bar&   b   = bars_[i];
        const double cn0 = std::clamp( b.cn0_db_hz, 0.0, CN0_MAX );
        const double h   = ( cn0 / CN0_MAX ) * plot.height();
        const double x   = plot.left() + i * slot + gap_half;
        const QRectF bar( x, plot.bottom() - h, bar_w, h );
        bar_rects_[i] = bar;

        QColor fill = b.color;
        if( !b.lock )
        {
            fill.setAlpha( 90 ); // dim the unlocked (acquiring-quality) bars
        }
        p.fillRect( bar, fill );

        // Two-line label under the bar: SV on top, signal component below (so a "G04" + "L1C" pair does
        // not run together at narrow slot widths). Dropped entirely if slots get too narrow to read.
        if( slot >= 16.0 )
        {
            const double lh = ( MARGIN_B - 2 ) / 2.0;
            p.setPen( QColor( 200, 200, 200 ) );
            p.drawText( QRectF( plot.left() + i * slot, plot.bottom() + 2, slot, lh ), Qt::AlignCenter, b.label );
            p.setPen( QColor( 150, 150, 150 ) ); // the code dimmer, as secondary identity
            p.drawText( QRectF( plot.left() + i * slot, plot.bottom() + 2 + lh, slot, lh ), Qt::AlignCenter, b.code );
        }
    }
}

bool Cn0_bar_widget::event( QEvent* event )
{
    if( event->type() == QEvent::ToolTip )
    {
        auto* help = static_cast<QHelpEvent*>( event );
        for( size_t i = 0; i < bar_rects_.size() && i < bars_.size(); ++i )
        {
            // Widen the hit area to the bar's full slot height so the whole column is hoverable.
            QRectF hit = bar_rects_[i];
            hit.setTop( MARGIN_T );
            hit.setBottom( height() - MARGIN_B );
            if( hit.contains( help->pos() ) )
            {
                QToolTip::showText(
                    help->globalPos(),
                    QStringLiteral( "%1 %2: %3 dB-Hz" )
                        .arg( bars_[i].label )
                        .arg( bars_[i].code )
                        .arg( bars_[i].cn0_db_hz, 0, 'f', 1 ),
                    this
                );
                return true;
            }
        }
        QToolTip::hideText();
        event->ignore();
        return true;
    }
    return QWidget::event( event );
}

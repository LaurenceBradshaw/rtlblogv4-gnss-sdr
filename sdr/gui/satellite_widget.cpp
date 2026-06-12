#include "satellite_widget.h"
#include <QFont>
#include <QHBoxLayout>
#include <QLabel>
#include <QMouseEvent>
#include "gnss_format.h"

namespace
{
using gui_format::constellation_prefix;

const char* state_text( Channel_state s )
{
    switch( s )
    {
    case Channel_state::TRACKING:
        return "TRAK";
    case Channel_state::ACQUIRING:
        return "ACQ";
    default:
        return "IDLE";
    }
}

// Shared column layout for the header and every row (monospace, so the columns line up):
//   SV(4)  STATE(6)  C/N0(7)  DOPPLER(12)  LOCK(6)  EPH(4)
constexpr const char* ROW_FORMAT = "%-4s %-6s %7s %12s %-6s %-4s";
} // namespace

Satellite_widget::Satellite_widget( QWidget* parent )
    : QWidget( parent )
{
    setCursor( Qt::PointingHandCursor );
    setToolTip( QStringLiteral( "Click for details" ) );

    auto* layout = new QHBoxLayout( this );
    layout->setContentsMargins( 6, 1, 6, 1 );

    line_ = new QLabel( this );
    line_->setFont( QFont( QStringLiteral( "monospace" ) ) );
    layout->addWidget( line_ );
}

QString Satellite_widget::header_text()
{
    return QString::asprintf( ROW_FORMAT, "SV", "STATE", "C/N0", "DOPPLER", "LOCK", "EPH" );
}

void Satellite_widget::update_snapshot( const Channel_snapshot& s )
{
    constellation_ = s.constellation;
    prn_           = static_cast<int>( s.satellite_id );

    const QString sv       = QString::asprintf( "%s%02d", constellation_prefix( s.constellation ), s.satellite_id );
    const bool    tracking = ( s.state == Channel_state::TRACKING );

    // Only a tracking SV has meaningful carrier metrics; acquiring/idle rows show placeholders.
    const QString cn0  = tracking ? QString::asprintf( "%.1f", s.cn0_db_hz ) : QStringLiteral( "-" );
    const QString dopp = tracking ? QString::asprintf( "%+.0f Hz", s.carrier_doppler_hz ) : QStringLiteral( "-" );
    const QString lock = ( tracking && s.has_lock ) ? QStringLiteral( "true" ) : QStringLiteral( "-" );
    const QString eph  = s.has_observable ? QStringLiteral( "valid" ) : QStringLiteral( "-" );

    line_->setText(
        QString::asprintf(
            ROW_FORMAT,
            sv.toUtf8().constData(),
            state_text( s.state ),
            cn0.toUtf8().constData(),
            dopp.toUtf8().constData(),
            lock.toUtf8().constData(),
            eph.toUtf8().constData()
        )
    );
}

void Satellite_widget::mousePressEvent( QMouseEvent* event )
{
    if( event->button() == Qt::LeftButton )
    {
        emit clicked( constellation_, prn_ );
    }
}

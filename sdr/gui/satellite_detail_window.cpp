#include "satellite_detail_window.h"
#include <QCloseEvent>
#include <QFont>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>
#include "doppler_graph_window.h"
#include "gnss_format.h"
#include "iq_constellation_window.h"

namespace
{
using gui_format::constellation_name;
using gui_format::constellation_prefix;

const char* state_text( Channel_state s )
{
    switch( s )
    {
    case Channel_state::TRACKING:
        return "TRACKING";
    case Channel_state::ACQUIRING:
        return "ACQUIRING";
    default:
        return "IDLE";
    }
}

QString sv_label( Constellation c, int prn )
{
    return QString::asprintf( "%s%02d  (%s PRN %d)", constellation_prefix( c ), prn, constellation_name( c ), prn );
}
} // namespace

Satellite_detail_window::Satellite_detail_window(
    Constellation constellation, int prn, const Receiver& receiver, QWidget* parent
)
    : QWidget( parent, Qt::Window ) // separate top-level window even though parented (for cleanup on exit)
    , constellation_( constellation )
    , prn_( prn )
    , receiver_( receiver )
{
    setAttribute( Qt::WA_DeleteOnClose );
    setWindowTitle( sv_label( constellation, prn ) );
    resize( 360, 560 );

    auto* layout = new QVBoxLayout( this );
    text_        = new QLabel( this );
    text_->setFont( QFont( QStringLiteral( "monospace" ) ) );
    text_->setTextInteractionFlags( Qt::TextSelectableByMouse );
    text_->setAlignment( Qt::AlignTop | Qt::AlignLeft );
    layout->addWidget( text_ );

    // Graph buttons. Always enabled (a not-yet-tracking SV just opens an empty graph that fills in once
    // it locks). More graph types will be added here.
    auto* graphs = new QHBoxLayout;
    auto* iq_btn = new QPushButton( QStringLiteral( "I/Q constellation" ), this );
    connect( iq_btn, &QPushButton::clicked, this, &Satellite_detail_window::open_iq_constellation );
    graphs->addWidget( iq_btn );
    auto* dop_btn = new QPushButton( QStringLiteral( "Doppler" ), this );
    connect( dop_btn, &QPushButton::clicked, this, &Satellite_detail_window::open_doppler );
    graphs->addWidget( dop_btn );
    graphs->addStretch();
    layout->addLayout( graphs );

    Channel_snapshot initial;
    initial.constellation = constellation;
    initial.satellite_id  = static_cast<Satellite_id>( prn );
    update_snapshot( initial );
}

void Satellite_detail_window::open_iq_constellation()
{
    if( iq_window_ ) // already open -> raise it
    {
        iq_window_->raise();
        iq_window_->activateWindow();
        return;
    }
    iq_window_ = new Iq_constellation_window( constellation_, prn_, receiver_, this );
    connect( iq_window_, &QObject::destroyed, this, [this] { iq_window_ = nullptr; } );
    iq_window_->show();
}

void Satellite_detail_window::open_doppler()
{
    if( doppler_window_ )
    {
        doppler_window_->raise();
        doppler_window_->activateWindow();
        return;
    }
    doppler_window_ = new Doppler_graph_window( constellation_, prn_, receiver_, this );
    connect( doppler_window_, &QObject::destroyed, this, [this] { doppler_window_ = nullptr; } );
    doppler_window_->show();
}

void Satellite_detail_window::update_snapshot( const Channel_snapshot& s )
{
    QString t;
    const auto line = [&t]( const char* label, const QString& value ) {
        t += QString::asprintf( "%-16s ", label ) + value + '\n';
    };
    const auto num = [&line]( const char* label, double v, const char* fmt = "%.6g" ) {
        line( label, QString::asprintf( fmt, v ) );
    };

    t += sv_label( s.constellation, s.satellite_id ) + "\n\n";

    line( "State", QString::fromLatin1( state_text( s.state ) ) );
    line( "Carrier lock", s.has_lock ? QStringLiteral( "yes" ) : QStringLiteral( "no" ) );
    line( "Has observable", s.has_observable ? QStringLiteral( "yes" ) : QStringLiteral( "no" ) );

    const bool tracking = ( s.state == Channel_state::TRACKING );
    line( "C/N0 (dB-Hz)", tracking ? QString::asprintf( "%.1f", s.cn0_db_hz ) : QStringLiteral( "-" ) );
    line( "Doppler (Hz)", tracking ? QString::asprintf( "%+.1f", s.carrier_doppler_hz ) : QStringLiteral( "-" ) );
    line( "Carrier acc",
          tracking ? QString::asprintf( "%+.3f Hz/s", s.carrier_acceleration ) : QStringLiteral( "-" ) );
    num( "Wavelength (m)", s.wavelength_m );
    num( "Transmit time (s)", s.transmission_time_s, "%.3f" );

    t += "\n--- Ephemeris ---\n";
    const Ephemeris& e = s.eph;
    if( !e.valid )
    {
        t += "not decoded yet\n";
    }
    else
    {
        line( "Week", QString::number( e.week ) );
        num( "TOW (s)", e.tow, "%.1f" );
        num( "Toe (s)", e.toe, "%.1f" );
        num( "Toc (s)", e.toc, "%.1f" );
        num( "af0 (s)", e.af0 );
        num( "af1 (s/s)", e.af1 );
        num( "af2 (s/s^2)", e.af2 );
        num( "sqrtA (m^.5)", e.sqrt_a, "%.6f" );
        num( "Eccentricity", e.e );
        num( "M0 (rad)", e.m0 );
        num( "delta_n (rad/s)", e.delta_n );
        num( "Omega0 (rad)", e.omega0 );
        num( "omega (rad)", e.omega );
        num( "OmegaDot", e.omegadot );
        num( "i0 (rad)", e.i0 );
        num( "iDot (rad/s)", e.idot );
        num( "Cuc / Cus", e.cuc );
        num( "Crc / Crs", e.crc );
        num( "Cic / Cis", e.cic );
        num( "Group delay (s)", e.group_delay );
    }

    text_->setText( t );
}

void Satellite_detail_window::closeEvent( QCloseEvent* event )
{
    emit closed( constellation_, prn_ ); // the list drops us from its map before WA_DeleteOnClose fires
    QWidget::closeEvent( event );
}

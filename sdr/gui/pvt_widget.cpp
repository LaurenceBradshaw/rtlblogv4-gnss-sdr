#include "pvt_widget.h"
#include <QFont>
#include <QLabel>
#include <QVBoxLayout>
#include "gui_frame.h"
#include "map_widget.h"

Pvt_widget::Pvt_widget( QWidget* parent )
    : Gui_panel( parent )
{
    auto* layout = new QVBoxLayout( this );

    status_   = new QLabel( this );
    geodetic_ = new QLabel( this );
    position_ = new QLabel( this );
    velocity_ = new QLabel( this );
    clock_    = new QLabel( this );
    map_      = new Map_widget( this );

    // Monospace so the columns line up as the numbers change.
    const QFont mono( QStringLiteral( "monospace" ) );
    for( QLabel* label : { geodetic_, position_, velocity_, clock_ } )
    {
        label->setFont( mono );
        label->setTextInteractionFlags( Qt::TextSelectableByMouse );
    }

    layout->addWidget( status_ );
    layout->addWidget( geodetic_ );
    layout->addWidget( position_ );
    layout->addWidget( velocity_ );
    layout->addWidget( clock_ );
    layout->addWidget( map_, 1 ); // map takes the remaining space

    update_frame( Gui_frame {} ); // start in the "no fix" state
}

void Pvt_widget::update_frame( const Gui_frame& frame )
{
    if( !frame.position || !frame.position->valid )
    {
        status_->setText( QStringLiteral( "No fix yet - waiting for satellites..." ) );
        geodetic_->clear();
        position_->clear();
        velocity_->clear();
        clock_->clear();
        map_->set_fix( 0.0, 0.0, false );
        return;
    }

    const Position_solution& p = *frame.position;
    status_->setText( QStringLiteral( "Fix valid" ) );

    if( frame.has_geodetic )
    {
        geodetic_->setText( QString::asprintf(
            "Lat %11.6f deg   Lon %11.6f deg   Alt %9.1f m", frame.lat_deg, frame.lon_deg, frame.alt_m
        ) );
        map_->set_fix( frame.lat_deg, frame.lon_deg, true );
    }

    position_->setText( QString::asprintf(
        "ECEF  x = %12.1f m    y = %12.1f m    z = %12.1f m", p.ecef_x_m, p.ecef_y_m, p.ecef_z_m
    ) );
    velocity_->setText( QString::asprintf(
        "Vel   x = %8.2f m/s  y = %8.2f m/s  z = %8.2f m/s", p.ecef_x_m_s, p.ecef_y_m_s, p.ecef_z_m_s
    ) );
    // Reference constellation is always shown; an inter-system bias only for non-reference
    // constellations that are actually present (so a single-constellation fix shows none).
    QString ref_isb = QString::asprintf( "    ref = %s", constellation_name( p.reference ) );
    for( int c = 0; c < NUM_CONSTELLATIONS; ++c )
    {
        if( p.isb_present[c] )
        {
            ref_isb += QString::asprintf(
                "    %s ISB = %.1f m", constellation_name( static_cast<Constellation>( c ) ), p.isb_m[c]
            );
        }
    }
    clock_->setText(
        QString::asprintf( "Clock bias = %.3f m    drift = %.3f m/s", p.clock_bias_m, p.clock_drift_m_s ) + ref_isb
    );
}

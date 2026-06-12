#include "doppler_graph_window.h"
#include <QComboBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QVBoxLayout>
#include "doppler_plot_widget.h"
#include "gnss_format.h"

Doppler_graph_window::Doppler_graph_window(
    Constellation constellation, int prn, const Receiver& receiver, QWidget* parent
)
    : Graph_window(
          constellation,
          prn,
          receiver,
          QString::asprintf( "%s%02d  Carrier Doppler", gui_format::constellation_prefix( constellation ), prn ),
          parent
      )
{
    resize( 520, 300 );
    auto* layout = new QVBoxLayout( this );

    auto* controls = new QHBoxLayout;
    controls->addWidget( new QLabel( QStringLiteral( "Span:" ), this ) );
    mode_ = new QComboBox( this );
    mode_->addItem( QStringLiteral( "Seconds" ) ); // index 0 -> fast ring
    mode_->addItem( QStringLiteral( "Minutes" ) ); // index 1 -> slow ring
    connect( mode_, QOverload<int>::of( &QComboBox::currentIndexChanged ), this, [this]( int ) { refresh(); } );
    controls->addWidget( mode_ );
    controls->addStretch();
    layout->addLayout( controls );

    plot_ = new Doppler_plot_widget( this );
    layout->addWidget( plot_ );

    start_polling();
}

void Doppler_graph_window::refresh()
{
    if( auto h = history() )
    {
        plot_->set_series( mode_->currentIndex() == 0 ? h->doppler.fast : h->doppler.slow );
    }
}

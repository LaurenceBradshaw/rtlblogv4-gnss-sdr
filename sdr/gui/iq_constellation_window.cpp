#include "iq_constellation_window.h"
#include <QVBoxLayout>
#include "gnss_format.h"
#include "iq_heatmap_widget.h"

Iq_constellation_window::Iq_constellation_window(
    Constellation constellation, int prn, Code code, const Receiver& receiver, QWidget* parent
)
    : Graph_window(
          constellation,
          prn,
          code,
          receiver,
          QString::asprintf(
              "%s%02d %s  I/Q constellation",
              gui_format::constellation_prefix( constellation ),
              prn,
              gui_format::code_label( code )
          ),
          parent
      )
{
    resize( 360, 400 );
    auto* layout = new QVBoxLayout( this );
    heatmap_     = new Iq_heatmap_widget( this );
    layout->addWidget( heatmap_ );

    start_polling();
}

void Iq_constellation_window::refresh()
{
    if( auto h = history() )
    {
        heatmap_->set_snapshot( h->iq );
    }
}

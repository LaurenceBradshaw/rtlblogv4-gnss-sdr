#pragma once
#include <vector>
#include "gui_frame.h" // Sky_satellite
#include "gui_panel.h"

// The "Sky" tab: a polar az/el plot of the tracked satellites as seen from the receiver. Zenith at
// the centre, horizon at the rim, north up, azimuth clockwise; rings at 30/60 deg elevation. Each SV
// is a dot coloured by constellation (filled = carrier lock) and labelled with its PRN.
class Sky_plot_widget : public Gui_panel
{
    Q_OBJECT
public:
    explicit Sky_plot_widget( QWidget* parent = nullptr );

    QString tab_title() const override
    {
        return QStringLiteral( "Sky" );
    }
    void update_frame( const Gui_frame& frame ) override;

protected:
    void paintEvent( QPaintEvent* event ) override;

private:
    std::vector<Sky_satellite> sats_;
};

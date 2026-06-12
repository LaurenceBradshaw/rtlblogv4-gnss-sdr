#pragma once
#include "gui_panel.h"

class QLabel;
class Map_widget;

// Shows the latest PVT solution: geodetic lat/lon/alt, ECEF position and velocity, receiver clock
// bias/drift and the inter-system bias, plus a live map marking the fix.
class Pvt_widget : public Gui_panel
{
    Q_OBJECT
public:
    explicit Pvt_widget( QWidget* parent = nullptr );

    QString tab_title() const override
    {
        return QStringLiteral( "Position" );
    }
    void update_frame( const Gui_frame& frame ) override;

private:
    QLabel*     status_   = nullptr;
    QLabel*     geodetic_ = nullptr;
    QLabel*     position_ = nullptr;
    QLabel*     velocity_ = nullptr;
    QLabel*     clock_    = nullptr;
    Map_widget* map_      = nullptr;
};

#pragma once
#include <QString>
#include <QWidget>

struct Gui_frame;

// Base for every GUI panel (one per tab). A panel is a QWidget that knows how to name its tab and
// refresh itself from a Gui_frame. Main_window holds a list of panels and pushes each new frame to
// all of them, so adding a panel is just: write the subclass, register it in one place. No panel
// knows about any other, and none touch the receiver - they only read the frame they are handed.
class Gui_panel : public QWidget
{
public:
    using QWidget::QWidget;
    ~Gui_panel() override = default;

    // The tab label for this panel.
    virtual QString tab_title() const = 0;
    // Refresh the panel's contents from the latest frame (called on the GUI thread, ~15 Hz).
    virtual void update_frame( const Gui_frame& frame ) = 0;
};

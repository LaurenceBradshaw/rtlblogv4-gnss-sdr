#pragma once
#include <vector>
#include <QMainWindow>
#include "receiver_view.h"

class QTabWidget;
class QTimer;
class QLabel;
class QPushButton;
class Gui_panel;
class Signals_widget;
class Source_widget;
class Receiver;
class Receiver_controller;

// Top-level window: one tab per Gui_panel, all refreshed from the Receiver on a timer. Panels are
// registered in the constructor (one add_panel() each) - the refresh loop and tab wiring are
// generic, so adding a panel touches nothing else. A Start/Stop toolbar drives the receiver via the
// controller (the receiver does not auto-run).
class Main_window : public QMainWindow
{
    Q_OBJECT
public:
    Main_window( const Receiver& receiver, Receiver_controller& controller, QWidget* parent = nullptr );

private:
    void add_panel( Gui_panel* panel ); // take ownership (via Qt parent) + add a tab
    void refresh();                     // poll the receiver, push the frame to every panel + status bar
    void on_start();                    // stage the Signals-tab selection, then start the receiver
    void stage_signals();               // apply the Signals-tab selection now (live pre-Start preview)

    Receiver_view        view_;
    Receiver_controller& controller_;
    Source_widget*       source_panel_  = nullptr; // for staging the source/run params on Start
    Signals_widget*      signals_panel_ = nullptr; // for staging the selection on Start

    QTabWidget*             tabs_          = nullptr;
    QTimer*                 timer_         = nullptr;
    QPushButton*            start_button_  = nullptr;
    QPushButton*            stop_button_   = nullptr;
    QLabel*                 status_label_  = nullptr; // exec / stream / channel times
    QLabel*                 warning_label_ = nullptr; // "behind real-time" - shown only when behind
    std::vector<Gui_panel*> panels_;
};

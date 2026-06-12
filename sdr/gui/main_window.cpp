#include "main_window.h"
#include <QFont>
#include <QLabel>
#include <QPushButton>
#include <QStatusBar>
#include <QTabWidget>
#include <QTimer>
#include <QToolBar>
#include "gui_frame.h"
#include "gui_panel.h"
#include "pvt_widget.h"
#include "receiver_controller.h"
#include "satellite_list_widget.h"
#include "sky_plot_widget.h"

namespace
{
constexpr int    REFRESH_INTERVAL_MS = 66;  // ~15 Hz
constexpr double BEHIND_WARN_S       = 0.2; // show the warning once channels lag exec by this (matches console)
} // namespace

Main_window::Main_window( const Receiver& receiver, Receiver_controller& controller, QWidget* parent )
    : QMainWindow( parent ),
      view_( receiver ),
      controller_( controller )
{
    setWindowTitle( QStringLiteral( "GNSS SDR" ) );
    resize( 960, 600 );

    // Start/Stop toolbar. The receiver does not auto-run; Start launches it (an IQ file replays from
    // the beginning), Stop halts + tears it down so Start can run a fresh pass. refresh() keeps the
    // two buttons' enabled state in sync with the controller.
    auto* toolbar = addToolBar( QStringLiteral( "Controls" ) );
    toolbar->setMovable( false );
    start_button_ = new QPushButton( QStringLiteral( "Start" ), this );
    stop_button_  = new QPushButton( QStringLiteral( "Stop" ), this );
    stop_button_->setEnabled( false );
    connect( start_button_, &QPushButton::clicked, this, [this] { controller_.start(); } );
    connect( stop_button_, &QPushButton::clicked, this, [this] { controller_.stop(); } );
    toolbar->addWidget( start_button_ );
    toolbar->addWidget( stop_button_ );

    tabs_ = new QTabWidget( this );
    setCentralWidget( tabs_ );

    // Register panels - one line each, in tab order.
    add_panel( new Pvt_widget );
    add_panel( new Satellite_list_widget( receiver ) );
    add_panel( new Sky_plot_widget );

    // Persistent status bar (part of QMainWindow, so it stays across tab changes). The exec/stream/
    // channel readout sits on the left; the "behind real-time" warning is appended and only shown when
    // the receiver actually lags.
    status_label_ = new QLabel( this );
    status_label_->setFont( QFont( QStringLiteral( "monospace" ) ) );
    statusBar()->addWidget( status_label_ );

    warning_label_ = new QLabel( this );
    warning_label_->setStyleSheet( QStringLiteral( "color: #e03030; font-weight: bold;" ) );
    warning_label_->setContentsMargins( 24, 0, 0, 0 );
    statusBar()->addWidget( warning_label_ );

    timer_ = new QTimer( this );
    connect( timer_, &QTimer::timeout, this, &Main_window::refresh );
    timer_->start( REFRESH_INTERVAL_MS );
}

void Main_window::add_panel( Gui_panel* panel )
{
    panel->setParent( tabs_ ); // Qt owns it
    panels_.push_back( panel );
    tabs_->addTab( panel, panel->tab_title() );
}

void Main_window::refresh()
{
    // Keep the buttons in step with the controller (which also flips back to stopped on its own when
    // an IQ file drains or the receiver errors out).
    const bool running = controller_.is_running();
    start_button_->setEnabled( !running );
    stop_button_->setEnabled( running );

    const Gui_frame frame = view_.poll();
    for( Gui_panel* panel : panels_ )
    {
        panel->update_frame( frame );
    }

    const Receiver_status& s = frame.status;
    status_label_->setText(
        QString::asprintf( "exec (s): %.1f   |   stream (s): %.1f   |   channel (s): %.1f", s.exec_s, s.stream_s, s.channel_s )
    );

    const double behind = s.exec_s - s.channel_s;
    if( behind > BEHIND_WARN_S )
    {
        warning_label_->setText( QString::asprintf( "[WARNING: channels are %.1fs behind real-time]", behind ) );
    }
    else
    {
        warning_label_->clear();
    }
}

#include "graph_window.h"
#include <QTimer>
#include "receiver.h"

namespace
{
constexpr int REFRESH_MS = 50; // 20 Hz (history is published at ~10 Hz; this just samples it)
}

Graph_window::Graph_window(
    Constellation constellation, int prn, const Receiver& receiver, const QString& title, QWidget* parent
)
    : QWidget( parent, Qt::Window )
    , constellation_( constellation )
    , prn_( prn )
    , receiver_( receiver )
{
    setAttribute( Qt::WA_DeleteOnClose );
    setWindowTitle( title );
    receiver_.subscribe_history( constellation_, prn_, true );
}

Graph_window::~Graph_window()
{
    receiver_.subscribe_history( constellation_, prn_, false );
}

void Graph_window::start_polling()
{
    timer_ = new QTimer( this );
    connect( timer_, &QTimer::timeout, this, &Graph_window::refresh );
    timer_->start( REFRESH_MS );
    refresh();
}

std::optional<Tracking_history::Snapshot> Graph_window::history() const
{
    return receiver_.published_history( constellation_, prn_ );
}

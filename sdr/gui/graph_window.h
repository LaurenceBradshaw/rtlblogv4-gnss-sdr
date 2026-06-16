#pragma once
#include <optional>
#include <QWidget>
#include "constellations.h"    // Constellation
#include "signal.h"            // Code
#include "tracking_history.h"  // Tracking_history::Snapshot

class Receiver;
class QTimer;

// Base for per-SV graph windows. A top-level window bound to one (constellation, PRN). Owns the shared
// lifecycle: subscribe that SV's history on the Receiver while open, unsubscribe on close, and a poll
// timer. Derived classes build their view widget(s) and implement refresh() to pull the published
// history (via history()) and update the view; they must call start_polling() at the END of their
// constructor (so the first/timer-driven refresh() runs only once the derived object is fully built).
class Graph_window : public QWidget
{
    Q_OBJECT
public:
    Graph_window(
        Constellation constellation, int prn, Code code, const Receiver& receiver, const QString& title,
        QWidget* parent = nullptr
    );
    ~Graph_window() override; // unsubscribe

protected:
    virtual void refresh() = 0; // pull history() + update the view (called on a timer)
    void         start_polling();
    std::optional<Tracking_history::Snapshot> history() const; // the published history for this SV

    Constellation   constellation_;
    int             prn_;
    Code            code_;
    const Receiver& receiver_;

private:
    QTimer* timer_ = nullptr;
};

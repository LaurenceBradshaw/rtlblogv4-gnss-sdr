#pragma once
#include <map>
#include <vector>
#include "constellations.h" // Constellation
#include "gui_panel.h"
#include "signal.h" // Code

struct Channel_snapshot;
class Receiver;
class QLabel;
class QVBoxLayout;
class Satellite_widget;
class Satellite_detail_window;
class Cn0_bar_widget;

// The "Satellites" tab: a row per EVERY configured SV (idle/acquiring ones show placeholders), with
// TRACKING satellites sorted to the top. Each row is clickable and opens a per-SV detail window, which
// this widget keeps refreshed while open. Each frame reconciles the channel snapshots into the rows.
class Satellite_list_widget : public Gui_panel
{
    Q_OBJECT
public:
    explicit Satellite_list_widget( const Receiver& receiver, QWidget* parent = nullptr );

    QString tab_title() const override
    {
        return QStringLiteral( "Satellites" );
    }
    void update_frame( const Gui_frame& frame ) override;

private:
    // Keyed by (constellation, prn, code) so an SV tracked on two components (e.g. GPS L1CA + L1C) gets a
    // row, detail window, and C/N0 bar per code rather than colliding into one.
    static int key_of( Constellation constellation, int prn, Code code );
    void       open_detail( Constellation constellation, int prn, Code code ); // row-click handler

    const Receiver&                         receiver_;              // for opening detail/graph windows
    QLabel*                                 summary_     = nullptr; // "Tracking N of M satellites"
    QVBoxLayout*                            rows_layout_ = nullptr; // holds the rows (+ a trailing stretch)
    Cn0_bar_widget*                         cn0_bars_    = nullptr; // C/N0 bar chart (right half)
    std::map<int, Satellite_widget*>        rows_;                  // key -> row (std::map keeps SV order)
    std::map<int, Satellite_detail_window*> detail_windows_;        // key -> open detail window
    std::vector<int>                        current_order_;         // keys in current layout order
};

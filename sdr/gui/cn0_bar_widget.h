#pragma once
#include <vector>
#include <QColor>
#include <QRectF>
#include <QString>
#include <QWidget>

// A C/N0 bar chart for the Satellites panel: one vertical bar per tracking satellite, height
// proportional to its carrier-to-noise density (dB-Hz), coloured by constellation. Hovering a bar
// shows the exact C/N0 value as a tooltip. Bars with carrier lock are drawn solid; tracking-but-
// unlocked ones are drawn dimmed. Pure QPainter, no external plotting dependency.
class Cn0_bar_widget : public QWidget
{
    Q_OBJECT
public:
    struct Bar
    {
        QString label;     // SV label, e.g. "G05" (drawn on the top line under the bar)
        QString code;      // signal component, e.g. "L1C" (drawn on the line below, to avoid overlap)
        double  cn0_db_hz; // current C/N0
        QColor  color;     // constellation colour
        bool    lock;      // carrier lock (solid vs dimmed)
    };

    explicit Cn0_bar_widget( QWidget* parent = nullptr );

    // Replace the displayed set (already in the desired left-to-right order).
    void set_bars( std::vector<Bar> bars );

protected:
    void paintEvent( QPaintEvent* event ) override;
    bool event( QEvent* event ) override; // tooltip hit-testing

private:
    std::vector<Bar>     bars_;
    std::vector<QRectF>  bar_rects_; // screen rects from the last paint, parallel to bars_
};

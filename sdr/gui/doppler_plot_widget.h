#pragma once
#include <QString>
#include <QWidget>
#include "tracking_history.h" // Tracking_history::Doppler_series

// A simple time-series line plot for Doppler (deepskyblue line on black). Given a
// Doppler_series it auto-scales both axes; the owning window picks which series (fast/seconds or
// slow/minutes) to show.
class Doppler_plot_widget : public QWidget
{
    Q_OBJECT
public:
    explicit Doppler_plot_widget( QWidget* parent = nullptr );

    void set_series( const Tracking_history::Doppler_series& series );

protected:
    void paintEvent( QPaintEvent* event ) override;

private:
    Tracking_history::Doppler_series series_;
};

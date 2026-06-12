#pragma once
#include <QWidget>
#include "tracking_history.h" // Tracking_history::Iq_snapshot

// Renders a channel's prompt I/Q constellation as a plasma heatmap (the live equivalent of the Python
// scatter, but density/recency-coloured from a fixed grid). I horizontal, Q vertical, origin centred.
class Iq_heatmap_widget : public QWidget
{
    Q_OBJECT
public:
    explicit Iq_heatmap_widget( QWidget* parent = nullptr );

    void set_snapshot( const Tracking_history::Iq_snapshot& snapshot );

protected:
    void paintEvent( QPaintEvent* event ) override;

private:
    Tracking_history::Iq_snapshot snap_;
};

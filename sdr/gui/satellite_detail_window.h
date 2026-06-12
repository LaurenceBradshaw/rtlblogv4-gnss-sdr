#pragma once
#include <QWidget>
#include "channel.h" // Channel_snapshot, Constellation

class QLabel;
class Receiver;
class Iq_constellation_window;
class Doppler_graph_window;

// A separate top-level window with full detail for ONE satellite: identity, tracking state, carrier
// metrics (C/N0, Doppler, acceleration), lock/ephemeris flags, and the broadcast ephemeris fields.
// The Satellite_list_widget creates one per clicked SV, pushes the live snapshot each frame, and drops
// it from its map when the window closes (the window deletes itself via WA_DeleteOnClose).
class Satellite_detail_window : public QWidget
{
    Q_OBJECT
public:
    Satellite_detail_window( Constellation constellation, int prn, const Receiver& receiver, QWidget* parent = nullptr );

    Constellation constellation() const
    {
        return constellation_;
    }
    int prn() const
    {
        return prn_;
    }

    void update_snapshot( const Channel_snapshot& snapshot );

signals:
    void closed( Constellation constellation, int prn );

protected:
    void closeEvent( QCloseEvent* event ) override;

private:
    void open_iq_constellation(); // graph button -> opens / raises the I/Q heatmap window
    void open_doppler();          // graph button -> opens / raises the Doppler-over-time window

    Constellation            constellation_;
    int                      prn_;
    const Receiver&          receiver_;
    QLabel*                  text_           = nullptr; // monospace block of all fields
    Iq_constellation_window* iq_window_      = nullptr; // open graph windows (if any), so we raise vs dup
    Doppler_graph_window*    doppler_window_ = nullptr;
};

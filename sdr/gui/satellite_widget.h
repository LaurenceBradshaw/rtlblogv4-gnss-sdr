#pragma once
#include <QString>
#include <QWidget>
#include "channel.h" // Channel_snapshot, Constellation

class QLabel;

// One row in the satellite list: a single satellite's live state (PRN, state, C/N0, Doppler, carrier
// lock, ephemeris). Shown for EVERY configured SV - acquiring/idle ones display placeholders. The row
// is clickable; clicking emits clicked() so the list can open a detail window for that SV.
class Satellite_widget : public QWidget
{
    Q_OBJECT
public:
    explicit Satellite_widget( QWidget* parent = nullptr );

    void update_snapshot( const Channel_snapshot& snapshot );

    // The column-header line whose layout matches a row (used by the list's header).
    static QString header_text();

signals:
    void clicked( Constellation constellation, int prn );

protected:
    void mousePressEvent( QMouseEvent* event ) override;

private:
    QLabel*       line_          = nullptr;
    Constellation constellation_ = Constellation::Unknown; // this row's SV identity (for clicked())
    int           prn_           = 0;
};

#pragma once
#include <QCache>
#include <QColor>
#include <QNetworkAccessManager>
#include <QPixmap>
#include <QPoint>
#include <QSet>
#include <QString>
#include <QWidget>
#include <vector>

// A lightweight slippy map: drag to pan, wheel to zoom down to street level, using OpenStreetMap
// raster tiles fetched on demand via Qt6Network and cached to disk (so revisited areas need no
// further requests). Plots a marker at the receiver fix and, in "follow" mode, keeps it centred
// (double-click re-centres + resumes following). Pure Qt Widgets - no Qt Location / QML.
class Map_widget : public QWidget
{
    Q_OBJECT
public:
    explicit Map_widget( QWidget* parent = nullptr );

    // Set the receiver fix (degrees). While following, the map re-centres on it; the first valid
    // fix also zooms in from the initial world view to a local view.
    void set_fix( double lat_deg, double lon_deg, bool valid );

    // A satellite's sub-satellite (nadir) ground point, drawn as a labelled marker on the map.
    struct Sat_marker
    {
        double  lat_deg;
        double  lon_deg;
        QColor  color;
        QString label;  // e.g. "G05"
        bool    filled; // carrier lock -> filled, else hollow
    };
    // Replace the set of satellite markers shown this frame.
    void set_satellites( std::vector<Sat_marker> sats );

protected:
    void paintEvent( QPaintEvent* event ) override;
    void wheelEvent( QWheelEvent* event ) override;
    void mousePressEvent( QMouseEvent* event ) override;
    void mouseMoveEvent( QMouseEvent* event ) override;
    void mouseReleaseEvent( QMouseEvent* event ) override;
    void mouseDoubleClickEvent( QMouseEvent* event ) override;

private:
    QPixmap tile_pixmap( int z, int x, int y );  // cache lookup; fetches on miss
    void    request_tile( int z, int x, int y ); // async network fetch (dedup'd)
    void    screen_to_lonlat( QPointF px, double& lat_deg, double& lon_deg ) const;
    void    set_center_from_world_px( double world_x, double world_y ); // -> center_lat_/lon_
    void    update_follow();                                            // dead-zone + eased follow of the fix

    // View state
    double center_lat_ = 20.0;
    double center_lon_ = 0.0;
    int    zoom_       = 2;
    bool   follow_     = true;
    bool   first_fix_  = true;

    // Fix marker
    bool   has_fix_ = false;
    double fix_lat_ = 0.0;
    double fix_lon_ = 0.0;

    // Satellite ground markers (sub-satellite points), redrawn each frame.
    std::vector<Sat_marker> sats_;

    // Panning
    bool   dragging_ = false;
    QPoint last_pos_;

    // Tile caches + networking
    QNetworkAccessManager    net_;
    QCache<QString, QPixmap> mem_cache_; // key "z/x/y"
    QSet<QString>            in_flight_;
    QString                  cache_dir_;
};

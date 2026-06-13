#include "map_widget.h"
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QMouseEvent>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QPainter>
#include <QStandardPaths>
#include <QWheelEvent>
#include <algorithm>
#include <cmath>
#include "logging.h"

namespace
{
constexpr double TILE            = 256.0;
constexpr int    MIN_ZOOM        = 2;
constexpr int    MAX_ZOOM        = 19;
constexpr double LAT_LIMIT       = 85.05112878; // web-mercator latitude clamp
constexpr int    FIRST_FIX_ZOOM  = 15;          // local view on the first fix
constexpr int    MEM_CACHE_TILES = 512;

// Follow-mode smoothing. The fix may move within FOLLOW_DEADZONE_M (metres) of centre WITHOUT panning -
// this kills the static-receiver jitter; being in metres, its on-screen size scales with zoom (it grows
// ~2x per zoom-in). Beyond it the centre EASES toward keeping the fix at the dead-zone edge (a low-pass
// that smooths the residual motion). EARTH_CIRC_M maps metres to web-mercator pixels.
constexpr double FOLLOW_DEADZONE_M = 5.0;
constexpr double FOLLOW_EASE       = 0.15; // fraction of the way to the target per update (~15 Hz)
constexpr double EARTH_CIRC_M      = 40075016.686;

double clamp_lat( double lat )
{
    return std::max( -LAT_LIMIT, std::min( LAT_LIMIT, lat ) );
}

// Web-mercator <-> fractional tile coordinates at zoom z.
double lon_to_tile_x( double lon_deg, int z )
{
    return ( lon_deg + 180.0 ) / 360.0 * ( 1 << z );
}
double lat_to_tile_y( double lat_deg, int z )
{
    const double lat = lat_deg * M_PI / 180.0;
    return ( 1.0 - std::asinh( std::tan( lat ) ) / M_PI ) / 2.0 * ( 1 << z );
}
double tile_x_to_lon( double x, int z )
{
    return x / ( 1 << z ) * 360.0 - 180.0;
}
double tile_y_to_lat( double y, int z )
{
    const double t = M_PI - 2.0 * M_PI * y / ( 1 << z );
    return 180.0 / M_PI * std::atan( std::sinh( t ) );
}

QString tile_key( int z, int x, int y )
{
    return QStringLiteral( "%1/%2/%3" ).arg( z ).arg( x ).arg( y );
}
} // namespace

Map_widget::Map_widget( QWidget* parent )
    : QWidget( parent )
{
    setMinimumHeight( 240 );
    mem_cache_.setMaxCost( MEM_CACHE_TILES );
    cache_dir_ = QStandardPaths::writableLocation( QStandardPaths::CacheLocation ) + QStringLiteral( "/osm_tiles" );

    logging::log( logging::Level::Info, "Map tile cache directory: " + cache_dir_.toStdString() );
}

void Map_widget::set_fix( double lat_deg, double lon_deg, bool valid )
{
    has_fix_ = valid;
    if( !valid )
    {
        return;
    }
    fix_lat_ = lat_deg;
    fix_lon_ = lon_deg;

    if( first_fix_ )
    {
        first_fix_  = false;
        zoom_       = FIRST_FIX_ZOOM;
        center_lat_ = clamp_lat( lat_deg );
        center_lon_ = lon_deg;
    }
    else if( follow_ )
    {
        update_follow();
    }
    update();
}

void Map_widget::set_satellites( std::vector<Sat_marker> sats )
{
    sats_ = std::move( sats );
    update();
}

void Map_widget::update_follow()
{
    // Work in world pixels at the current zoom.
    const double fx = lon_to_tile_x( fix_lon_, zoom_ ) * TILE;
    const double fy = lat_to_tile_y( fix_lat_, zoom_ ) * TILE;
    double       cx = lon_to_tile_x( center_lon_, zoom_ ) * TILE;
    double       cy = lat_to_tile_y( center_lat_, zoom_ ) * TILE;

    // Dead zone in pixels = metres -> web-mercator pixels at this zoom/latitude (so it scales with zoom).
    const double world_px = TILE * ( 1 << zoom_ );
    const double px_per_m = world_px / ( EARTH_CIRC_M * std::cos( fix_lat_ * M_PI / 180.0 ) );
    const double dz       = FOLLOW_DEADZONE_M * px_per_m;

    // Target centre: move only enough to keep the fix within the dead zone; inside it, hold position.
    const double dx = fx - cx;
    const double dy = fy - cy;
    double       tx = cx;
    double       ty = cy;
    if( std::abs( dx ) > dz )
    {
        tx = fx - std::copysign( dz, dx );
    }
    if( std::abs( dy ) > dz )
    {
        ty = fy - std::copysign( dz, dy );
    }

    // Ease the centre toward the target (low-pass smoothing of the residual motion).
    cx += FOLLOW_EASE * ( tx - cx );
    cy += FOLLOW_EASE * ( ty - cy );
    set_center_from_world_px( cx, cy );
}

void Map_widget::paintEvent( QPaintEvent* )
{
    QPainter p( this );
    p.fillRect( rect(), QColor( 30, 30, 30 ) );

    const int    z        = zoom_;
    const int    n        = 1 << z;
    const double cx_px    = lon_to_tile_x( center_lon_, z ) * TILE;
    const double cy_px    = lat_to_tile_y( center_lat_, z ) * TILE;
    const double origin_x = cx_px - width() / 2.0; // world px at screen (0,0)
    const double origin_y = cy_px - height() / 2.0;

    const int tx0 = static_cast<int>( std::floor( origin_x / TILE ) );
    const int ty0 = static_cast<int>( std::floor( origin_y / TILE ) );
    const int tx1 = static_cast<int>( std::floor( ( origin_x + width() ) / TILE ) );
    const int ty1 = static_cast<int>( std::floor( ( origin_y + height() ) / TILE ) );

    for( int ty = ty0; ty <= ty1; ++ty )
    {
        if( ty < 0 || ty >= n )
        {
            continue; // no vertical wrap
        }
        for( int tx = tx0; tx <= tx1; ++tx )
        {
            const int     wx = ( ( tx % n ) + n ) % n; // wrap horizontally
            const double  sx = tx * TILE - origin_x;
            const double  sy = ty * TILE - origin_y;
            const QPixmap pm = tile_pixmap( z, wx, ty );
            if( !pm.isNull() )
            {
                p.drawPixmap( QPointF( sx, sy ), pm );
            }
            else
            {
                p.fillRect( QRectF( sx, sy, TILE, TILE ), QColor( 50, 50, 50 ) );
            }
        }
    }

    // Satellite sub-satellite (nadir) markers. These sit far from the fix (the orbit is ~20 000 km
    // up), so they are mostly only on-screen when zoomed out; draw only those within the viewport.
    p.setRenderHint( QPainter::Antialiasing, true );
    for( const Sat_marker& s : sats_ )
    {
        const double sx = lon_to_tile_x( s.lon_deg, z ) * TILE - origin_x;
        const double sy = lat_to_tile_y( clamp_lat( s.lat_deg ), z ) * TILE - origin_y;
        if( sx < -20 || sx > width() + 20 || sy < -20 || sy > height() + 20 )
        {
            continue; // off-screen
        }
        p.setPen( QPen( s.color.darker( 150 ), 1.5 ) );
        p.setBrush( s.filled ? QBrush( s.color ) : Qt::NoBrush );
        p.drawEllipse( QPointF( sx, sy ), 4.0, 4.0 );
        p.setPen( s.color.lighter( 130 ) );
        p.drawText( QPointF( sx + 6.0, sy + 4.0 ), s.label );
    }

    if( has_fix_ )
    {
        const double fx = lon_to_tile_x( fix_lon_, z ) * TILE - origin_x;
        const double fy = lat_to_tile_y( fix_lat_, z ) * TILE - origin_y;
        p.setRenderHint( QPainter::Antialiasing, true );
        p.setPen( QPen( Qt::white, 2 ) );
        p.setBrush( QColor( 220, 40, 40 ) );
        p.drawEllipse( QPointF( fx, fy ), 6.0, 6.0 );
    }

    // Attribution (OSM tile usage policy requires it).
    p.setPen( Qt::white );
    p.drawText( rect().adjusted( 4, 4, -4, -4 ), Qt::AlignBottom | Qt::AlignRight, QStringLiteral( "© OpenStreetMap" ) );
}

QPixmap Map_widget::tile_pixmap( int z, int x, int y )
{
    const QString key = tile_key( z, x, y );
    if( QPixmap* cached = mem_cache_.object( key ) )
    {
        return *cached;
    }
    const QString path = cache_dir_ + '/' + key + QStringLiteral( ".png" );
    QPixmap       pm;
    if( pm.load( path ) )
    {
        mem_cache_.insert( key, new QPixmap( pm ) );
        return pm;
    }
    request_tile( z, x, y );
    return {};
}

void Map_widget::request_tile( int z, int x, int y )
{
    const QString key = tile_key( z, x, y );
    if( in_flight_.contains( key ) )
    {
        return;
    }
    in_flight_.insert( key );

    QNetworkRequest req( QUrl( QStringLiteral( "https://tile.openstreetmap.org/%1.png" ).arg( key ) ) );
    req.setHeader( QNetworkRequest::UserAgentHeader, QStringLiteral( "rtlblogv4-gnss-sdr/0.1 (hobby GNSS receiver)" ) );

    QNetworkReply* reply = net_.get( req );
    connect(
        reply,
        &QNetworkReply::finished,
        this,
        [this, reply, key]()
        {
            reply->deleteLater();
            in_flight_.remove( key );
            if( reply->error() != QNetworkReply::NoError )
            {
                return;
            }
            const QByteArray data = reply->readAll();
            QPixmap          pm;
            if( !pm.loadFromData( data ) )
            {
                return;
            }
            mem_cache_.insert( key, new QPixmap( pm ) );

            const QString path = cache_dir_ + '/' + key + QStringLiteral( ".png" );
            QDir().mkpath( QFileInfo( path ).absolutePath() );
            QFile f( path );
            if( f.open( QIODevice::WriteOnly ) )
            {
                f.write( data );
            }
            update();
        }
    );
}

void Map_widget::screen_to_lonlat( QPointF px, double& lat_deg, double& lon_deg ) const
{
    const int    z       = zoom_;
    const double cx_px   = lon_to_tile_x( center_lon_, z ) * TILE;
    const double cy_px   = lat_to_tile_y( center_lat_, z ) * TILE;
    const double world_x = cx_px - width() / 2.0 + px.x();
    const double world_y = cy_px - height() / 2.0 + px.y();
    lon_deg              = tile_x_to_lon( world_x / TILE, z );
    lat_deg              = tile_y_to_lat( world_y / TILE, z );
}

void Map_widget::set_center_from_world_px( double world_x, double world_y )
{
    center_lon_ = tile_x_to_lon( world_x / TILE, zoom_ );
    center_lat_ = clamp_lat( tile_y_to_lat( world_y / TILE, zoom_ ) );
}

void Map_widget::wheelEvent( QWheelEvent* event )
{
    const int step     = ( event->angleDelta().y() > 0 ) ? 1 : -1;
    const int new_zoom = std::max( MIN_ZOOM, std::min( MAX_ZOOM, zoom_ + step ) );
    if( new_zoom == zoom_ )
    {
        return;
    }

    // Keep the geo point under the cursor fixed across the zoom.
    double lat = 0.0;
    double lon = 0.0;
    screen_to_lonlat( event->position(), lat, lon );
    zoom_           = new_zoom;
    const double px = lon_to_tile_x( lon, zoom_ ) * TILE;
    const double py = lat_to_tile_y( lat, zoom_ ) * TILE;
    set_center_from_world_px( px - ( event->position().x() - width() / 2.0 ), py - ( event->position().y() - height() / 2.0 ) );
    follow_ = false;
    update();
}

void Map_widget::mousePressEvent( QMouseEvent* event )
{
    if( event->button() == Qt::LeftButton )
    {
        dragging_ = true;
        last_pos_ = event->pos();
        follow_   = false;
    }
}

void Map_widget::mouseMoveEvent( QMouseEvent* event )
{
    if( !dragging_ )
    {
        return;
    }
    const QPoint delta = event->pos() - last_pos_;
    last_pos_          = event->pos();
    const double cx_px = lon_to_tile_x( center_lon_, zoom_ ) * TILE - delta.x();
    const double cy_px = lat_to_tile_y( center_lat_, zoom_ ) * TILE - delta.y();
    set_center_from_world_px( cx_px, cy_px );
    update();
}

void Map_widget::mouseReleaseEvent( QMouseEvent* event )
{
    if( event->button() == Qt::LeftButton )
    {
        dragging_ = false;
    }
}

void Map_widget::mouseDoubleClickEvent( QMouseEvent* )
{
    // Re-centre on the fix and resume following.
    if( has_fix_ )
    {
        follow_     = true;
        center_lat_ = clamp_lat( fix_lat_ );
        center_lon_ = fix_lon_;
        update();
    }
}

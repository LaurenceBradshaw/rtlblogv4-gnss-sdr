#include "receiver_view.h"
#include <cmath>
#include "geodesy.h"
#include "orbit.h"
#include "receiver.h"

Receiver_view::Receiver_view( const Receiver& receiver )
    : receiver_( receiver )
{
}

Gui_frame Receiver_view::poll() const
{
    Gui_frame frame;
    frame.position = receiver_.latest_position();
    frame.channels = receiver_.channel_snapshots();
    frame.status   = receiver_.status();

    // Before Start (or before the first tick) there are no live channels yet - show every configured
    // satellite as an idle placeholder so the list is fully populated from the moment the GUI opens.
    if( frame.channels.empty() )
    {
        for( const Configured_satellite& sat : receiver_.configured_satellites() )
        {
            Channel_snapshot s;
            s.constellation = sat.constellation;
            s.satellite_id  = static_cast<Satellite_id>( sat.prn );
            s.state         = Channel_state::IDLE;
            frame.channels.push_back( s );
        }
    }

    if( frame.position && frame.position->valid )
    {
        const Ecef user { frame.position->ecef_x_m, frame.position->ecef_y_m, frame.position->ecef_z_m };

        const Geodetic g   = ecef_to_geodetic( user );
        frame.has_geodetic = true;
        frame.lat_deg      = g.lat_rad * 180.0 / M_PI;
        frame.lon_deg      = g.lon_rad * 180.0 / M_PI;
        frame.alt_m        = g.alt_m;

        // Look angles of every tracked SV with a usable ephemeris, for the sky plot. The SV transmit
        // time is a fine enough time base to propagate the (slow) orbit for az/el.
        for( const Channel_snapshot& s : frame.channels )
        {
            if( !s.has_observable )
            {
                continue;
            }
            const Ecef sat = orbit::satellite_ecef_pos( s.eph, s.transmission_time_s, s.constellation );
            double     elevation_rad = 0.0;
            double     azimuth_rad   = 0.0;
            look_angles( user, sat, elevation_rad, azimuth_rad );

            // Sub-satellite ground point (nadir) for the map marker.
            const Geodetic satg = ecef_to_geodetic( sat );

            Sky_satellite sky;
            sky.constellation = s.constellation;
            sky.satellite_id  = s.satellite_id;
            sky.azimuth_deg   = azimuth_rad * 180.0 / M_PI;
            sky.elevation_deg = elevation_rad * 180.0 / M_PI;
            sky.cn0_db_hz     = s.cn0_db_hz;
            sky.has_lock      = s.has_lock;
            sky.sub_lat_deg   = satg.lat_rad * 180.0 / M_PI;
            sky.sub_lon_deg   = satg.lon_rad * 180.0 / M_PI;
            frame.sky.push_back( sky );
        }
    }

    return frame;
}

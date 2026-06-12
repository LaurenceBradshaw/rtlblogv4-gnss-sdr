#pragma once
#include "types.h"

// Receiver-frame geometry helpers (constellation-independent). WGS-84.
// Used by the atmosphere models (need receiver lat/lon/height + satellite look angles)
// and later by elevation-based measurement weighting in the position solver.

namespace WGS84
{
// WGS84 ellipsoid constants
constexpr double A   = 6378137.0;                     // semi-major axis (m)
constexpr double F   = 1.0 / 298.257223563;           // flattening
constexpr double B   = A * ( 1.0 - F );               // semi-minor axis (m)
constexpr double E2  = ( A * A - B * B ) / ( A * A ); // eccentricity^2
constexpr double EP2 = ( A * A - B * B ) / ( B * B ); // second eccentricity^2
} // namespace WGS84

struct Geodetic
{
    double lat_rad; // geodetic latitude  (rad)
    double lon_rad; // geodetic longitude (rad)
    double alt_m;   // height above the WGS-84 ellipsoid (m)
};

// WGS-84 ECEF -> geodetic latitude/longitude/height.
// (Closed-form Bowring)
Geodetic ecef_to_geodetic( const Ecef& p );

// Look angles of `sat` as seen from `user` (both ECEF, metres).
// elevation_rad in [-pi/2, pi/2] (negative = below the local horizon),
// azimuth_rad in [0, 2*pi) measured clockwise from true north.
// (Form the ENU vector to the satellite, then elevation = asin(up/range),
//  azimuth = atan2(east, north).)
void look_angles( const Ecef& user, const Ecef& sat, double& elevation_rad, double& azimuth_rad );

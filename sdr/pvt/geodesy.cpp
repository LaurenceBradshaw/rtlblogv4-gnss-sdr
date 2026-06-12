#include "geodesy.h"
#include <cmath>

Geodetic ecef_to_geodetic( const Ecef& ecef )
{
    double x = ecef.x;
    double y = ecef.y;
    double z = ecef.z;

    double longitude = std::atan2( y, x );
    double p         = std::sqrt( x * x + y * y );

    // Bowring’s method for initial latitude
    double theta    = std::atan2( z * WGS84::A, p * WGS84::B );
    double latitude = std::atan2(
        z + WGS84::EP2 * WGS84::B * std::pow( std::sin( theta ), 3 ),
        p - WGS84::E2 * WGS84::A * std::pow( std::cos( theta ), 3 )
    );

    // radius of curvature in the prime vertical
    double sin_lat = std::sin( latitude );
    double N       = WGS84::A / std::sqrt( 1.0 - WGS84::E2 * sin_lat * sin_lat );

    double altitude = ( p / std::cos( latitude ) ) - N;

    return Geodetic { latitude, longitude, altitude };
}

void look_angles( const Ecef& user, const Ecef& sat, double& elevation_rad, double& azimuth_rad )
{
    Ecef     los { sat.x - user.x, sat.y - user.y, sat.z - user.z };
    Geodetic user_geo = ecef_to_geodetic( user );

    double sin_lat = std::sin( user_geo.lat_rad );
    double cos_lat = std::cos( user_geo.lat_rad );
    double sin_lon = std::sin( user_geo.lon_rad );
    double cos_lon = std::cos( user_geo.lon_rad );

    double east  = -sin_lon * los.x + cos_lon * los.y;
    double north = -sin_lat * cos_lon * los.x - sin_lat * sin_lon * los.y + cos_lat * los.z;
    double up    = cos_lat * cos_lon * los.x + cos_lat * sin_lon * los.y + sin_lat * los.z;

    double los_norm = std::sqrt( east * east + north * north + up * up );

    elevation_rad = std::asin( up / los_norm );
    azimuth_rad   = std::atan2( east, north );
    if( azimuth_rad < 0.0 )
    {
        azimuth_rad += 2.0 * M_PI;
    }
    else if( azimuth_rad >= 2.0 * M_PI )
    {
        azimuth_rad -= 2.0 * M_PI;
    }
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace
{
// Forward WGS-84 geodetic -> ECEF (the standard closed form), used here only to build
// independent inputs for the round-trip and look-angle tests.
Ecef geodetic_to_ecef( double lat_rad, double lon_rad, double alt_m )
{
    const double s = std::sin( lat_rad );
    const double N = WGS84::A / std::sqrt( 1.0 - WGS84::E2 * s * s );
    return {
        ( N + alt_m ) * std::cos( lat_rad ) * std::cos( lon_rad ),
        ( N + alt_m ) * std::cos( lat_rad ) * std::sin( lon_rad ),
        ( N * ( 1.0 - WGS84::E2 ) + alt_m ) * s
    };
}
const double DEG = M_PI / 180.0;
} // namespace

TEST_CASE( "ecef_to_geodetic_exact_anchors", "[pvt][geodesy]" )
{
    // Independent exact points: a point on the equator at the prime meridian sits at (A,0,0)
    // with zero height; 90 deg east is (0,A,0).
    const Geodetic g0 = ecef_to_geodetic( Ecef { WGS84::A, 0.0, 0.0 } );
    REQUIRE( g0.lat_rad == Catch::Approx( 0.0 ).margin( 1e-12 ) );
    REQUIRE( g0.lon_rad == Catch::Approx( 0.0 ).margin( 1e-12 ) );
    REQUIRE( g0.alt_m == Catch::Approx( 0.0 ).margin( 1e-6 ) );

    const Geodetic g90 = ecef_to_geodetic( Ecef { 0.0, WGS84::A, 0.0 } );
    REQUIRE( g90.lat_rad == Catch::Approx( 0.0 ).margin( 1e-12 ) );
    REQUIRE( g90.lon_rad == Catch::Approx( M_PI / 2.0 ) );
    REQUIRE( g90.alt_m == Catch::Approx( 0.0 ).margin( 1e-6 ) );
}

TEST_CASE( "ecef_to_geodetic_round_trip", "[pvt][geodesy]" )
{
    // geodetic -> ecef -> geodetic recovers the input for a spread of realistic sites
    // (incl. the CTTC capture site ~41.275 N / 1.987 E / 80 m).
    struct
    {
        double lat_deg, lon_deg, alt_m;
    } sites[] = {
        { 41.275, 1.987, 80.0 }, { 0.0, 0.0, 0.0 }, { 51.5, -0.1, 35.0 }, { -33.9, 151.2, 58.0 }, { 60.0, -120.0, 1500.0 }, { -70.0, 90.0, 0.0 }
    };
    for( const auto& s : sites )
    {
        const Ecef     p = geodetic_to_ecef( s.lat_deg * DEG, s.lon_deg * DEG, s.alt_m );
        const Geodetic g = ecef_to_geodetic( p );
        INFO( "site lat=" << s.lat_deg << " lon=" << s.lon_deg );
        REQUIRE( g.lat_rad == Catch::Approx( s.lat_deg * DEG ).margin( 1e-9 ) );
        REQUIRE( g.lon_rad == Catch::Approx( s.lon_deg * DEG ).margin( 1e-9 ) );
        REQUIRE( g.alt_m == Catch::Approx( s.alt_m ).margin( 1e-4 ) );
    }
}

TEST_CASE( "look_angles_cardinal_geometry", "[pvt][geodesy]" )
{
    // User on the equator at the prime meridian: local up=+x, east=+y, north=+z.
    const Ecef user { WGS84::A, 0.0, 0.0 };
    double     el = 0.0, az = 0.0;

    // Straight up -> elevation 90 deg.
    look_angles( user, Ecef { WGS84::A + 20.0e6, 0.0, 0.0 }, el, az );
    REQUIRE( el == Catch::Approx( M_PI / 2.0 ) );

    // On the horizon due north (+z) -> elevation 0, azimuth 0.
    look_angles( user, Ecef { WGS84::A, 0.0, 1.0e6 }, el, az );
    REQUIRE( el == Catch::Approx( 0.0 ).margin( 1e-9 ) );
    REQUIRE( az == Catch::Approx( 0.0 ).margin( 1e-9 ) );

    // On the horizon due east (+y) -> elevation 0, azimuth 90 deg.
    look_angles( user, Ecef { WGS84::A, 1.0e6, 0.0 }, el, az );
    REQUIRE( el == Catch::Approx( 0.0 ).margin( 1e-9 ) );
    REQUIRE( az == Catch::Approx( M_PI / 2.0 ) );

    // A satellite below the local horizon (toward -x) -> negative elevation.
    look_angles( user, Ecef { WGS84::A - 1.0e6, 0.0, 0.0 }, el, az );
    REQUIRE( el == Catch::Approx( -M_PI / 2.0 ) );
}
#endif

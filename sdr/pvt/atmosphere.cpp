#include "atmosphere.h"
#include <cmath>
#include "constants.h"

double klobuchar_iono_delay_m(
    const double alpha[4],
    const double beta[4],
    double       user_lat_rad,
    double       user_lon_rad,
    double       elevation_rad,
    double       azimuth_rad,
    double       gps_tow_s
)
{
    // TODO(you): IS-GPS-200 20.3.3.5.2.5 Klobuchar model. Outline:
    //   1. earth-centred angle psi (semicircles) from elevation
    //   2. ionospheric pierce point geodetic lat/lon (phi_i, lambda_i)
    //   3. geomagnetic latitude phi_m
    //   4. local time t = 43200*lambda_i + gps_tow  (wrap to [0, 86400))
    //   5. amplitude AMP = sum alpha_n * phi_m^n   (>= 0)
    //      period PER = sum beta_n * phi_m^n       (>= 72000)
    //   6. phase x = 2*pi*(t - 50400)/PER
    //   7. slant factor F = 1 + 16*(0.53 - el)^3   (el in semicircles)
    //   8. T_iono(s) = F * (5e-9 + AMP*(1 - x^2/2 + x^4/24))   if |x| < 1.57
    //                = F * 5e-9                                  otherwise
    //   return T_iono * c   (metres)

    double psi = 0.0137 / ( elevation_rad / M_PI + 0.11 ) - 0.022; // earth-centred angle (semicircles)

    double phi_i = user_lat_rad / M_PI + psi * std::cos( azimuth_rad ); // ionospheric pierce point latitude (semicircles)
    if( phi_i > 0.416 )
    {
        phi_i = 0.416;
    }
    else if( phi_i < -0.416 )
    {
        phi_i = -0.416;
    }

    double lambda_i = user_lon_rad / M_PI + ( psi * std::sin( azimuth_rad ) ) /
                                                std::cos( phi_i * M_PI ); // ionospheric pierce point longitude (semicircles)

    double phi_m = phi_i + 0.064 * std::cos( ( lambda_i - 1.617 ) * M_PI ); // geomagnetic latitude (semicircles)

    double t = std::fmod( 43200.0 * lambda_i + gps_tow_s, 86400.0 ); // local time (s)
    if( t < 0.0 )
    {
        t += 86400.0;
    }

    double amp = 0.0; // amplitude (s)
    double per = 0.0; // period (s)
    for( int n = 0; n < 4; n++ )
    {
        amp += alpha[n] * std::pow( phi_m, n );
        per += beta[n] * std::pow( phi_m, n );
    }
    if( amp < 0.0 )
    {
        amp = 0.0;
    }
    if( per < 72000.0 )
    {
        per = 72000.0;
    }

    double x = 2.0 * M_PI * ( t - 50400.0 ) / per;                      // phase (rad)
    double F = 1.0 + 16.0 * std::pow( 0.53 - elevation_rad / M_PI, 3 ); // slant factor
    double T_iono_s;
    if( std::abs( x ) < 1.57 )
    {
        T_iono_s = F * ( 5e-9 + amp * ( 1.0 - std::pow( x, 2 ) / 2.0 + std::pow( x, 4 ) / 24.0 ) );
    }
    else
    {
        T_iono_s = F * 5e-9;
    }

    return T_iono_s * constants::SPEED_OF_LIGHT_M_S; // metres
}

double tropo_delay_m( double elevation_rad, double user_alt_m )
{
    // TODO(you): zenith tropo delay (standard atmosphere) * elevation mapping function.
    // A simple Saastamoinen/Black model is plenty for a few-metre-class fix.

    double zenith_delay_m   = 2.3;                             // zenith delay in standard atmosphere (m)
    double mapping_function = 1.0 / std::sin( elevation_rad ); // mapping function (slant/zenith)

    double tropo_delay = zenith_delay_m * mapping_function;

    return tropo_delay;
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace
{
const double DEG = M_PI / 180.0;

// The Klobuchar night-time floor: with zero amplitude coefficients the cosine series vanishes
// and the delay collapses to F * 5 ns * c, where F is the obliquity factor.
double klobuchar_floor_m( double elevation_rad )
{
    const double F = 1.0 + 16.0 * std::pow( 0.53 - elevation_rad / M_PI, 3 );
    return F * 5e-9 * constants::SPEED_OF_LIGHT_M_S;
}
} // namespace

TEST_CASE( "tropo_delay_zenith_and_mapping", "[pvt][atmosphere]" )
{
    // Zenith = the 2.3 m standard zenith delay; the 1/sin(el) mapping is exact off-zenith.
    REQUIRE( tropo_delay_m( M_PI / 2.0, 0.0 ) == Catch::Approx( 2.3 ) );
    REQUIRE( tropo_delay_m( 30.0 * DEG, 0.0 ) == Catch::Approx( 2.3 / std::sin( 30.0 * DEG ) ) ); // 4.6 m
    REQUIRE( tropo_delay_m( 45.0 * DEG, 0.0 ) == Catch::Approx( 2.3 / std::sin( 45.0 * DEG ) ) );

    // Monotonic: lower elevation -> longer slant path -> larger delay.
    REQUIRE( tropo_delay_m( 10.0 * DEG, 0.0 ) > tropo_delay_m( 60.0 * DEG, 0.0 ) );
}

TEST_CASE( "klobuchar_zero_coefficients_give_obliquity_floor", "[pvt][atmosphere]" )
{
    const double a[4] = { 0, 0, 0, 0 };
    const double b[4] = { 0, 0, 0, 0 };
    for( double el_deg : { 10.0, 30.0, 60.0, 90.0 } )
    {
        const double el = el_deg * DEG;
        const double d  = klobuchar_iono_delay_m( a, b, 0.7, 0.0, el, 0.0, 50000.0 );
        INFO( "el=" << el_deg );
        REQUIRE( d == Catch::Approx( klobuchar_floor_m( el ) ) );
    }
}

TEST_CASE( "klobuchar_realistic_is_positive_and_elevation_monotonic", "[pvt][atmosphere]" )
{
    // Representative broadcast Klobuchar coefficients (mid solar activity).
    const double a[4] = { 1.0e-8, 1.49e-8, -5.96e-8, -1.19e-7 };
    const double b[4] = { 8.8e4, 3.28e4, -1.97e5, -3.93e5 };
    const double lat = 41.0 * DEG, lon = 2.0 * DEG, az = 1.0;
    const double tow = 50400.0; // local-noon-ish, near the model's daytime peak

    const double low  = klobuchar_iono_delay_m( a, b, lat, lon, 10.0 * DEG, az, tow );
    const double high = klobuchar_iono_delay_m( a, b, lat, lon, 80.0 * DEG, az, tow );
    REQUIRE( high > 0.0 );
    REQUIRE( low > high );                             // lower elevation -> larger slant delay
    REQUIRE( low < 100.0 );                            // sane upper bound (m)
    REQUIRE( high > klobuchar_floor_m( 80.0 * DEG ) ); // daytime exceeds the night floor
}
#endif

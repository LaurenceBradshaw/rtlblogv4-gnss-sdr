#include "orbit.h"
#include <cmath>

namespace orbit
{

static constexpr double OMEGA_E_DOT      = 7.2921151467e-5; // Earth's rotation rate (rad/s) - GPS & Galileo
static constexpr double C                = 299792458.0;     // Speed of light (m/s)
static constexpr double SECONDS_PER_WEEK = 604800.0;        // seconds in a (GPS/Galileo/BeiDou) week

// Gravitational parameter mu (m^3/s^2) per constellation's reference frame.
static double constellation_mu( Constellation c )
{
    switch( c )
    {
    case Constellation::Galileo:
        return 3.986004418e14; // Galileo Terrestrial Reference Frame
    case Constellation::Beidou:
        return 3.986004418e14; // BeiDou CGCS2000 (same value as Galileo)
    case Constellation::Gps:
    default:
        return 3.986005e14; // GPS WGS-84
    }
}

// Wrap a time difference into +/- half a week (toe/toc are seconds-of-week).
static double wrap_week( double dt )
{
    if( dt < -SECONDS_PER_WEEK / 2 )
    {
        return dt + SECONDS_PER_WEEK;
    }
    if( dt > SECONDS_PER_WEEK / 2 )
    {
        return dt - SECONDS_PER_WEEK;
    }
    return dt;
}

Ecef satellite_ecef_pos( const Ephemeris& eph, double t_system, Constellation c )
{
    // IS-GPS-200L table 20-IV (identical Keplerian propagation for Galileo E1 / BeiDou).
    const double MU = constellation_mu( c );

    double a  = eph.sqrt_a * eph.sqrt_a;
    double n0 = std::sqrt( MU / ( a * a * a ) );

    double tk = wrap_week( t_system - eph.toe );

    double n = n0 + eph.delta_n;
    double M = eph.m0 + n * tk;

    double E = M;
    for( int i = 0; i < 10; ++i )
    {
        E += ( M - E + eph.e * std::sin( E ) ) / ( 1 - eph.e * std::cos( E ) );
    }

    double nu = 2 * std::atan( std::sqrt( ( 1 + eph.e ) / ( 1 - eph.e ) ) * std::tan( E / 2 ) );

    double phi = nu + eph.omega;

    double u = phi + eph.cus * std::sin( 2 * phi ) + eph.cuc * std::cos( 2 * phi );
    double r = a * ( 1 - eph.e * std::cos( E ) ) + eph.crs * std::sin( 2 * phi ) + eph.crc * std::cos( 2 * phi );
    double i = eph.i0 + eph.idot * tk + eph.cis * std::sin( 2 * phi ) + eph.cic * std::cos( 2 * phi );

    double x_prime = r * std::cos( u );
    double y_prime = r * std::sin( u );

    double omega = eph.omega0 + ( eph.omegadot - OMEGA_E_DOT ) * tk - OMEGA_E_DOT * eph.toe;

    return {
        x_prime * std::cos( omega ) - y_prime * std::cos( i ) * std::sin( omega ),
        x_prime * std::sin( omega ) + y_prime * std::cos( i ) * std::cos( omega ),
        y_prime * std::sin( i )
    };
}

Ecef satellite_ecef_vel( const Ephemeris& eph, double t_system, Constellation c )
{
    const double MU = constellation_mu( c );

    // Shared with the position computation.
    double a  = eph.sqrt_a * eph.sqrt_a;
    double n0 = std::sqrt( MU / ( a * a * a ) );

    double tk = wrap_week( t_system - eph.toe );

    double n = n0 + eph.delta_n;
    double M = eph.m0 + n * tk;

    double E = M;
    for( int i = 0; i < 10; ++i )
    {
        E += ( M - E + eph.e * std::sin( E ) ) / ( 1 - eph.e * std::cos( E ) );
    }

    double nu  = 2 * std::atan( std::sqrt( ( 1 + eph.e ) / ( 1 - eph.e ) ) * std::tan( E / 2 ) );
    double phi = nu + eph.omega;

    double u = phi + eph.cus * std::sin( 2 * phi ) + eph.cuc * std::cos( 2 * phi );
    double r = a * ( 1 - eph.e * std::cos( E ) ) + eph.crs * std::sin( 2 * phi ) + eph.crc * std::cos( 2 * phi );
    double i = eph.i0 + eph.idot * tk + eph.cis * std::sin( 2 * phi ) + eph.cic * std::cos( 2 * phi );

    double x_prime = r * std::cos( u );
    double y_prime = r * std::sin( u );
    double omega   = eph.omega0 + ( eph.omegadot - OMEGA_E_DOT ) * tk - OMEGA_E_DOT * eph.toe;

    double x_pos = x_prime * std::cos( omega ) - y_prime * std::cos( i ) * std::sin( omega );
    double y_pos = x_prime * std::sin( omega ) + y_prime * std::cos( i ) * std::cos( omega );

    // Orbital-plane rates.
    double M_dot = n;
    double E_dot = M_dot / ( 1.0 - eph.e * std::cos( E ) );

    double nu_dot  = std::sqrt( 1.0 - eph.e * eph.e ) * E_dot / ( 1.0 - eph.e * std::cos( E ) );
    double phi_dot = nu_dot;

    double u_dot = phi_dot + 2.0 * phi_dot * ( eph.cus * std::cos( 2.0 * phi ) - eph.cuc * std::sin( 2.0 * phi ) );
    double r_dot = a * eph.e * std::sin( E ) * E_dot +
                   2.0 * phi_dot * ( eph.crs * std::cos( 2.0 * phi ) - eph.crc * std::sin( 2.0 * phi ) );
    double i_dot = eph.idot + 2.0 * phi_dot * ( eph.cis * std::cos( 2.0 * phi ) - eph.cic * std::sin( 2.0 * phi ) );

    double x_prime_dot = r_dot * std::cos( u ) - r * u_dot * std::sin( u );
    double y_prime_dot = r_dot * std::sin( u ) + r * u_dot * std::cos( u );

    double omega_dot = eph.omegadot - OMEGA_E_DOT;

    double dx = x_prime_dot * std::cos( omega ) - y_prime_dot * std::cos( i ) * std::sin( omega ) +
                y_prime * std::sin( i ) * i_dot * std::sin( omega ) - omega_dot * y_pos;
    double dy = x_prime_dot * std::sin( omega ) + y_prime_dot * std::cos( i ) * std::cos( omega ) -
                y_prime * std::sin( i ) * i_dot * std::cos( omega ) + omega_dot * x_pos;
    double dz = y_prime_dot * std::sin( i ) + y_prime * std::cos( i ) * i_dot;

    return { dx, dy, dz };
}

double satellite_clock_offset( const Ephemeris& eph, double t_sv, Constellation c )
{
    const double MU = constellation_mu( c );
    const double F  = -2.0 * std::sqrt( MU ) / ( C * C ); // relativistic term (s/m^0.5)

    double delta_t    = wrap_week( t_sv - eph.toc );
    double delta_t_sv = eph.af0 + eph.af1 * delta_t + eph.af2 * delta_t * delta_t - eph.group_delay;
    double t_system   = t_sv - delta_t_sv;

    double a  = eph.sqrt_a * eph.sqrt_a;
    double n0 = std::sqrt( MU / ( a * a * a ) );

    double tk = wrap_week( t_system - eph.toe );

    double n = n0 + eph.delta_n;
    double M = eph.m0 + n * tk;

    double E = M;
    for( int i = 0; i < 10; ++i )
    {
        E += ( M - E + eph.e * std::sin( E ) ) / ( 1 - eph.e * std::cos( E ) );
    }

    double delta_t_rel = F * eph.e * eph.sqrt_a * std::sin( E );

    return delta_t_sv + delta_t_rel;
}

double satellite_clock_drift( const Ephemeris& eph, double t_sv, Constellation c )
{
    const double MU = constellation_mu( c );
    const double F  = -2.0 * std::sqrt( MU ) / ( C * C );

    double delta_t    = wrap_week( t_sv - eph.toc );
    double delta_t_sv = eph.af0 + eph.af1 * delta_t + eph.af2 * delta_t * delta_t - eph.group_delay;
    double t_system   = t_sv - delta_t_sv;

    double a  = eph.sqrt_a * eph.sqrt_a;
    double n0 = std::sqrt( MU / ( a * a * a ) );

    double tk = wrap_week( t_system - eph.toe );

    double n = n0 + eph.delta_n;
    double M = eph.m0 + n * tk;

    double E = M;
    for( int i = 0; i < 10; ++i )
    {
        E += ( M - E + eph.e * std::sin( E ) ) / ( 1 - eph.e * std::cos( E ) );
    }

    // d(delta_t_sv)/dt (polynomial) + d(relativistic)/dt.
    double delta_t_sv_dot = eph.af1 + 2.0 * eph.af2 * delta_t;

    double M_dot = n;
    double E_dot = M_dot / ( 1.0 - eph.e * std::cos( E ) );

    double delta_t_rel_dot = F * eph.e * eph.sqrt_a * std::cos( E ) * E_dot;

    return delta_t_sv_dot + delta_t_rel_dot;
}

} // namespace orbit

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace
{
double norm( const Ecef& p )
{
    return std::sqrt( p.x * p.x + p.y * p.y + p.z * p.z );
}

// A representative (mildly eccentric, all correction terms non-zero) GPS broadcast ephemeris,
// so every branch of the propagation is exercised by the finite-difference checks.
Ephemeris realistic_gps_eph()
{
    Ephemeris e;
    e.valid       = true;
    e.sqrt_a      = 5153.65; // a ~ 26560 km
    e.e           = 0.006;
    e.m0          = 0.3;
    e.delta_n     = 4.5e-9;
    e.omega0      = 1.0;
    e.omega       = -1.7;
    e.omegadot    = -8.0e-9;
    e.i0          = 0.95;
    e.idot        = -2.0e-10;
    e.cuc         = 1.0e-6;
    e.cus         = 8.0e-6;
    e.crc         = 250.0;
    e.crs         = -30.0;
    e.cic         = -1.0e-7;
    e.cis         = 1.0e-7;
    e.toe         = 0.0;
    e.toc         = 0.0;
    e.af0         = 1.0e-4;
    e.af1         = 3.0e-12;
    e.af2         = 0.0;
    e.group_delay = 5.0e-9;
    return e;
}
} // namespace

TEST_CASE( "orbit_circular_radius_matches_semi_major_axis", "[pvt][orbit]" )
{
    // For a zero-eccentricity orbit with no correction terms, the radius is exactly a, so the
    // ECEF magnitude must equal a at any time (a pure rotation of the orbital-plane point).
    Ephemeris e;
    e.valid        = true;
    e.sqrt_a       = 5153.65;
    e.e            = 0.0;
    e.m0           = 0.5;
    e.i0           = 0.96;
    e.omega0       = 0.7;
    const double a = e.sqrt_a * e.sqrt_a;

    for( double t : { 0.0, 600.0, 1800.0, 3600.0 } )
    {
        const Ecef p = orbit::satellite_ecef_pos( e, t, Constellation::Gps );
        REQUIRE( norm( p ) == Catch::Approx( a ).epsilon( 1e-9 ) );
    }
    // Sanity: GPS orbital radius ~26 560 km.
    REQUIRE( a == Catch::Approx( 26.56e6 ).epsilon( 0.01 ) );
}

TEST_CASE( "orbit_velocity_is_position_derivative", "[pvt][orbit]" )
{
    // The analytic velocity must equal the central finite-difference of the position.
    const Ephemeris e = realistic_gps_eph();
    const double    h = 1.0;
    for( double t : { 100.0, 1000.0, 3000.0 } )
    {
        const Ecef pa = orbit::satellite_ecef_pos( e, t + h, Constellation::Gps );
        const Ecef pb = orbit::satellite_ecef_pos( e, t - h, Constellation::Gps );
        const Ecef fd { ( pa.x - pb.x ) / ( 2 * h ), ( pa.y - pb.y ) / ( 2 * h ), ( pa.z - pb.z ) / ( 2 * h ) };
        const Ecef v = orbit::satellite_ecef_vel( e, t, Constellation::Gps );
        INFO( "t=" << t );
        REQUIRE( v.x == Catch::Approx( fd.x ).epsilon( 1e-5 ) );
        REQUIRE( v.y == Catch::Approx( fd.y ).epsilon( 1e-5 ) );
        REQUIRE( v.z == Catch::Approx( fd.z ).epsilon( 1e-5 ) );
    }
    // Sanity: a GPS satellite's ECEF speed sits in a few-km/s band (less than the ~3.9 km/s
    // inertial speed because Earth rotates beneath it).
    const double speed = norm( orbit::satellite_ecef_vel( e, 100.0, Constellation::Gps ) );
    REQUIRE( speed > 2.0e3 );
    REQUIRE( speed < 4.5e3 );
}

TEST_CASE( "orbit_clock_drift_is_offset_derivative", "[pvt][orbit]" )
{
    const Ephemeris e = realistic_gps_eph();
    const double    h = 1.0;
    for( double t : { 100.0, 1000.0, 3000.0 } )
    {
        const double oa    = orbit::satellite_clock_offset( e, t + h, Constellation::Gps );
        const double ob    = orbit::satellite_clock_offset( e, t - h, Constellation::Gps );
        const double fd    = ( oa - ob ) / ( 2 * h );
        const double drift = orbit::satellite_clock_drift( e, t, Constellation::Gps );
        INFO( "t=" << t );
        REQUIRE( drift == Catch::Approx( fd ).epsilon( 1e-4 ).margin( 1e-18 ) );
    }
}

TEST_CASE( "orbit_galileo_radius_plausible", "[pvt][orbit]" )
{
    // Galileo uses a larger orbit (~29 600 km) and its own mu; check the magnitude band.
    Ephemeris e;
    e.valid      = true;
    e.sqrt_a     = 5440.6; // a ~ 29.6e6
    e.e          = 0.0;
    e.m0         = 0.2;
    e.i0         = 0.98;
    const Ecef p = orbit::satellite_ecef_pos( e, 300.0, Constellation::Galileo );
    REQUIRE( norm( p ) == Catch::Approx( 29.6e6 ).epsilon( 0.01 ) );
}
#endif

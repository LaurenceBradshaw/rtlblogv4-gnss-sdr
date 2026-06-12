#pragma once

// Atmospheric pseudorange corrections (metres of extra path delay). Both delays are
// positive: the signal is slowed/bent, so the measured pseudorange is longer than the
// geometric range. The Observation_engine subtracts them from each pseudorange.
//
// Kept free of any constellation type (raw coefficient arrays + angles) so the same
// models serve GPS now and other constellations later. Klobuchar is the L1 delay; for
// other carriers scale by (f_L1 / f)^2 at the call site when that day comes.

// Broadcast (Klobuchar) ionospheric slant delay on L1, in metres.
//   alpha[4], beta[4] : Klobuchar coefficients (GPS SF4 page 18).
//   user_lat_rad, user_lon_rad : receiver geodetic latitude / longitude.
//   elevation_rad, azimuth_rad : satellite look angles from the receiver.
//   gps_tow_s : receiver GPS time of week at reception (s).
// Reference: IS-GPS-200 section 20.3.3.5.2.5 (the half-cosine model). Note the spec works
// in semicircles - convert the angles, and clamp the night-time / large-x cases.
double klobuchar_iono_delay_m(
    const double alpha[4],
    const double beta[4],
    double       user_lat_rad,
    double       user_lon_rad,
    double       elevation_rad,
    double       azimuth_rad,
    double       gps_tow_s
);

// Tropospheric slant delay, in metres, from the satellite elevation and receiver height.
// A standard-atmosphere model is fine here (e.g. Saastamoinen zenith delay + a 1/sin(el)
// style mapping, or Black/Hopfield). No met sensors, so assume standard P/T/humidity.
double tropo_delay_m( double elevation_rad, double user_alt_m );

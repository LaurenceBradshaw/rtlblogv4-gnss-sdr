#pragma once
#include "constellations.h"
#include "ephemeris.h"
#include "types.h"

// Broadcast-ephemeris satellite position / velocity / clock. Shared across constellations that
// use the same Keplerian broadcast model with a 604800 s week (GPS L1 C/A, Galileo E1, BeiDou).
// The only per-constellation differences are the gravitational parameter mu (and the relativistic
// term F derived from it); the single-frequency clock group delay rides in Ephemeris::group_delay.
// t_system / t_sv are in that constellation's own system time (GPST for GPS, GST for Galileo) -
// any small inter-system offset is absorbed by the receiver's inter-system bias in the solver.
namespace orbit
{
Ecef   satellite_ecef_pos( const Ephemeris& eph, double t_system, Constellation c );
Ecef   satellite_ecef_vel( const Ephemeris& eph, double t_system, Constellation c );
double satellite_clock_offset( const Ephemeris& eph, double t_sv, Constellation c );
double satellite_clock_drift( const Ephemeris& eph, double t_sv, Constellation c );
} // namespace orbit

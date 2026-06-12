#pragma once
#include "constellations.h"
#include "types.h"

struct Ephemeris
{
    bool          valid = false;
    Constellation constellation;
    Satellite_id  prn;

    int    week     = 0;
    double tow      = 0.0;       // time of week at subframe start (s)
    double toc      = 0.0;       // clock reference time (s)
    double af0      = 0.0;       // clock bias (s)
    double af1      = 0.0;       // clock drift (s/s)
    double af2      = 0.0;       // clock drift rate (s/s^2)
    double toe      = 0.0;       // reference time of ephemeris (s)
    double sqrt_a   = 0.0;       // sqrt(semi-major axis) (m^0.5)
    double e        = 0.0;       // eccentricity
    double m0       = 0.0;       // mean anomaly at toe (rad)
    double delta_n  = 0.0;       // mean motion correction (rad/s)
    double omega0   = 0.0;       // right ascension at toe (rad)
    double omega    = 0.0;       // argument of perigee (rad)
    double omegadot = 0.0;       // rate of right ascension (rad/s)
    double i0       = 0.0;       // inclination at toe (rad)
    double idot     = 0.0;       // rate of inclination (rad/s)
    double cuc = 0.0, cus = 0.0; // latitude argument corrections (rad)
    double crc = 0.0, crs = 0.0; // radius corrections (m)
    double cic = 0.0, cis = 0.0; // inclination corrections (rad)

    double group_delay = 0.0; // single-frequency clock group delay (s): GPS Tgd / Galileo BGD
};

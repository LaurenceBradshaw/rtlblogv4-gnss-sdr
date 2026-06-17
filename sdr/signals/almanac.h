#pragma once
#include "constellations.h"
#include "ephemeris.h"
#include "types.h"

// A satellite's broadcast ALMANAC: the coarse Keplerian orbit + clock used for acquisition aiding
// (predicting which SVs are up and their Doppler), NOT for the PVT range solve (that needs the precise
// Ephemeris). One physical SV broadcasts the almanac for its WHOLE constellation, subcommutated over the
// broadcast, so any tracked SV slowly fills it in (GPS: subframes 4/5, 25 pages, ~12.5 min; Galileo:
// I/NAV word types 7-10; BeiDou: similar). Constellation-agnostic: the decode produces standard Keplerian
// units (each constellation's references/scales are applied there), so the orbit model just needs the
// constellation (for mu).
struct Almanac
{
    bool          valid         = false;
    Constellation constellation = Constellation::Unknown;
    Satellite_id  prn           = 0;

    int    week     = 0;   // system week the toa belongs to (taken from the decoded ephemeris)
    double toa      = 0.0; // almanac reference time of week (s)
    double e        = 0.0; // eccentricity
    double sqrt_a   = 0.0; // sqrt(semi-major axis) (m^0.5)
    double i0       = 0.0; // inclination (rad)
    double omega0   = 0.0; // right ascension of ascending node (rad)
    double omega    = 0.0; // argument of perigee (rad)
    double m0       = 0.0; // mean anomaly (rad)
    double omegadot = 0.0; // rate of right ascension (rad/s)
    double af0      = 0.0; // clock bias (s)
    double af1      = 0.0; // clock drift (s/s)
    int    health   = 0;   // SV health (0 = healthy)

    // The almanac is just a coarse ephemeris with the perturbation terms zeroed, so it can drive the same
    // orbit model (orbit::satellite_ecef_pos/_vel) for az/el + Doppler prediction. toe := toa; the
    // constellation carries through so the orbit uses the right mu.
    Ephemeris as_ephemeris() const
    {
        Ephemeris e_out;
        e_out.constellation = constellation;
        e_out.prn           = prn;
        e_out.valid         = valid;
        e_out.week          = week;
        e_out.toe = e_out.toc = toa;
        e_out.sqrt_a          = sqrt_a;
        e_out.e               = e;
        e_out.m0              = m0;
        e_out.omega0          = omega0;
        e_out.omega           = omega;
        e_out.omegadot        = omegadot;
        e_out.i0              = i0;
        e_out.af0             = af0;
        e_out.af1             = af1;
        // delta_n / idot / cuc.. cis / group_delay all stay 0 - the almanac carries no perturbations.
        return e_out;
    }
};

#pragma once

// Shared physical constants - one source of truth so a value (and its precision) can't drift between
// translation units. (The Earth-rotation rate was once 7.292115e-5 in one place and 7.2921151467e-5 in
// another - a ~1.5e-7 inconsistency in the Sagnac correction.) These are well-known fundamental constants;
// domain-specific DERIVED values (e.g. per-signal scale factors) stay local to their site.
namespace constants
{
constexpr double SPEED_OF_LIGHT_M_S        = 299792458.0;     // c (m/s) - exact by SI definition
constexpr double EARTH_ROTATION_RATE_RAD_S = 7.2921151467e-5; // WGS-84 / GTRF Earth rotation rate (rad/s)
} // namespace constants

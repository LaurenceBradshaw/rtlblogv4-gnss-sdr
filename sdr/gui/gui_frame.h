#pragma once
#include <optional>
#include <vector>
#include "channel.h"  // Channel_snapshot
#include "position.h" // Position_solution
#include "receiver.h" // Receiver_status

// Immutable per-frame view of everything the GUI renders this tick, assembled by Receiver_view
// from the Receiver's thread-safe accessors. Panels render from this and never touch the receiver
// or its channels directly. To surface new data to the GUI: add a field here and fill it in
// Receiver_view::poll() - no panel needs to know where the data came from.
// One satellite's look angles for the sky plot, derived (in Receiver_view::poll) from its broadcast
// ephemeris + the latest fix. Only satellites with a usable ephemeris and an available fix appear.
struct Sky_satellite
{
    Constellation constellation  = Constellation::Unknown;
    int           satellite_id   = 0;
    double        azimuth_deg    = 0.0; // clockwise from true north
    double        elevation_deg  = 0.0; // above the local horizon
    double        cn0_db_hz      = 0.0;
    bool          has_lock       = false;
};

struct Gui_frame
{
    std::optional<Position_solution> position; // nullopt until the PVT filter is seeded
    std::vector<Channel_snapshot>    channels; // one entry per receiver channel
    std::vector<Sky_satellite>       sky;      // tracked SVs with az/el (needs a fix); empty until fixed

    // Geodetic position derived from `position` (valid only when has_geodetic). Computed once in
    // Receiver_view::poll so every panel (read-out + map) shares it.
    bool   has_geodetic = false;
    double lat_deg      = 0.0;
    double lon_deg      = 0.0;
    double alt_m        = 0.0;

    Receiver_status status; // exec / stream / channel seconds (status bar)
};

#pragma once

// Constellations the receiver knows about. The enum ORDER is also the PRIORITY order used to
// choose the PVT reference constellation - the one whose system time the receiver clock bias is
// measured against. Every other constellation present gets its own inter-system bias in the
// solver (N constellations -> N-1 biases). Keep NUM_CONSTELLATIONS in sync when adding one.
enum class Constellation
{
    Gps,
    Galileo,
    Beidou,
    Unknown // sentinel (e.g. a not-yet-populated snapshot); MUST stay last so it is not counted
};

constexpr int NUM_CONSTELLATIONS = 3; // real constellations (Gps..Beidou); excludes Unknown

// A single integer identifying a satellite across constellations (constellation * 1000 + PRN). Used to
// key per-SV maps (acquisition predictions, GUI rows, graph history, ...). Per-CODE keys build on this as
// sv_key(con,prn) * 10 + code (Code has < 10 values; that *10 is kept where Code is in scope).
inline int sv_key( Constellation c, int prn )
{
    return static_cast<int>( c ) * 1000 + prn;
}

// Human-readable constellation name (for logs + the GUI).
inline const char* constellation_name( Constellation c )
{
    switch( c )
    {
    case Constellation::Gps:
        return "GPS";
    case Constellation::Galileo:
        return "Galileo";
    case Constellation::Beidou:
        return "BeiDou";
    default:
        return "?";
    }
}

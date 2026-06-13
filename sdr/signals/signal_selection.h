#pragma once
#include <set>
#include <string>
#include <vector>
#include "signal.h" // Constellation, Band, Code

// User-facing selection of WHICH signals + PRNs the receiver searches. Kept separate from the signal
// physics so both the CLI (main.cpp) and the GUI tab build the same config. The receiver currently
// knows one signal per constellation, so selection is by constellation token; the full
// (constellation, band, code) triple is retained so finer band/code selection can be added later.
struct Signal_id
{
    Constellation constellation;
    Band          band;
    Code          code;
};

inline bool operator==( const Signal_id& a, const Signal_id& b )
{
    return a.constellation == b.constellation && a.band == b.band && a.code == b.code;
}

// One selected signal plus the PRNs to search within it. An empty prns set means "all PRNs in this
// signal's range" - each signal carries its OWN PRN allowlist, so different constellations can search
// different subsets.
struct Signal_selection
{
    Signal_id     id;
    std::set<int> prns;
};

// The default set searched when none is specified: GPS L1 C/A + Galileo E1-B, all PRNs.
std::vector<Signal_selection> default_signal_selection();

// Every signal the receiver can search (grouped by constellation in enum order). Used by the GUI to
// build the selection tab; grows as new bands/codes are added.
std::vector<Signal_id> all_signals();

// Resolve a Signal_id from a constellation token ("gps", "galileo"/"gal", "beidou"/"bds") and an
// optional signal-component name ("l1ca", "e1b", "b1i"). Case-insensitive. An empty `signal` selects
// the constellation's default (only) signal. Throws std::invalid_argument on an unknown constellation
// or a signal name that doesn't belong to it.
Signal_id signal_from_tokens( const std::string& constellation, const std::string& signal = "" );

// Canonical lower-case constellation token ("gps" / "galileo" / "beidou" / "?").
const char* signal_token( const Signal_id& id );
// Canonical lower-case signal-component name ("l1ca" / "e1b" / "b1i" / "?").
const char* signal_name( const Signal_id& id );

// Build the signal selection from the CLI's two flags (each a flattened token stream, since the parser
// may comma-split a value):
//   signal_tokens : `CONSTELLATION[:COMPONENT]` each, e.g. "gps", "gps:l1ca", "galileo:e1b".
//   prn_tokens    : per-constellation PRN lists `CONSTELLATION:PRNS`, e.g. "gps:1,4,6-32" (the comma
//                   split leaves "gps:1","4","6-32" - an alpha-headed token starts a constellation, a
//                   numeric/range token attaches to it). PRNs apply to ALL of that constellation's
//                   selected components.
// Throws std::invalid_argument on an unknown signal/constellation, a malformed PRN, a PRN token before
// any constellation, or PRNs given for a constellation that is not in signal_tokens.
std::vector<Signal_selection>
parse_selection( const std::vector<std::string>& signal_tokens, const std::vector<std::string>& prn_tokens );

// Highest PRN for a constellation (its signal's sv_range upper bound). Used to resolve open-ended PRN
// ranges. 0 if the constellation has no known signal.
int max_prn_for( Constellation c );

// Parse a PRN range list like "1-4,7,9-20,32" into {1,2,3,4,7,9..20,32}. Open-ended ranges are allowed
// when max_prn>0: "N-" means N..max_prn and "-N" means 1..N. Throws std::invalid_argument on a malformed
// token (incl. an open-ended range when max_prn==0). Empty/whitespace -> empty set ("all PRNs").
std::set<int> parse_prn_list( const std::string& spec, int max_prn = 0 );

// Inverse of parse_prn_list: collapse a PRN set into a compact "1-4,7,9-11" string (empty set -> "").
std::string format_prn_list( const std::set<int>& prns );

// True if `prn` is selected by `filter` (an empty filter selects everything).
inline bool prn_selected( const std::set<int>& filter, int prn )
{
    return filter.empty() || filter.count( prn ) != 0;
}

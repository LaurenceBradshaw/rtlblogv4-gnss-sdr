#include "signal_selection.h"
#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>

namespace
{
std::string to_lower( std::string s )
{
    std::transform( s.begin(), s.end(), s.begin(), []( unsigned char c ) { return static_cast<char>( std::tolower( c ) ); } );
    return s;
}

std::string trim( const std::string& s )
{
    size_t a = 0, b = s.size();
    while( a < b && std::isspace( static_cast<unsigned char>( s[a] ) ) )
    {
        ++a;
    }
    while( b > a && std::isspace( static_cast<unsigned char>( s[b - 1] ) ) )
    {
        --b;
    }
    return s.substr( a, b - a );
}

bool has_alpha( const std::string& s )
{
    return std::any_of( s.begin(), s.end(), []( unsigned char c ) { return std::isalpha( c ) != 0; } );
}

// Split on `sep`, keeping every field (including empties, which callers ignore).
std::vector<std::string> split( const std::string& s, char sep )
{
    std::vector<std::string> out;
    size_t                   i = 0;
    while( i <= s.size() )
    {
        const size_t at = s.find( sep, i );
        const size_t end = ( at == std::string::npos ) ? s.size() : at;
        out.push_back( s.substr( i, end - i ) );
        if( at == std::string::npos )
        {
            break;
        }
        i = at + 1;
    }
    return out;
}

int parse_int( const std::string& s )
{
    size_t pos = 0;
    int    v   = 0;
    try
    {
        v = std::stoi( s, &pos );
    }
    catch( const std::exception& )
    {
        throw std::invalid_argument( "PRN list: not an integer: '" + s + "'" );
    }
    if( pos != s.size() )
    {
        throw std::invalid_argument( "PRN list: not an integer: '" + s + "'" );
    }
    return v;
}

// Parse one PRN token into `out`: "7" (single), "9-20" (range), "9-" (9..max_prn) or "-9" (1..9). The
// open-ended forms require max_prn>0.
void add_prn_token( std::set<int>& out, const std::string& raw, int max_prn )
{
    const std::string tok = trim( raw );
    if( tok.empty() )
    {
        return;
    }
    const size_t dash = tok.find( '-' );
    if( dash == std::string::npos )
    {
        out.insert( parse_int( tok ) );
        return;
    }
    const std::string lo_s = tok.substr( 0, dash );
    const std::string hi_s = tok.substr( dash + 1 );
    const int         lo   = lo_s.empty() ? 1 : parse_int( lo_s ); // "-N" -> from 1
    int               hi;
    if( hi_s.empty() ) // "N-" -> to the constellation's last PRN
    {
        if( max_prn <= 0 )
        {
            throw std::invalid_argument( "PRN range '" + tok + "' is open-ended but no maximum is known here" );
        }
        hi = max_prn;
    }
    else
    {
        hi = parse_int( hi_s );
    }
    if( lo > hi )
    {
        throw std::invalid_argument( "PRN range out of order: '" + tok + "'" );
    }
    for( int p = lo; p <= hi; ++p )
    {
        out.insert( p );
    }
}
} // namespace

std::vector<Signal_selection> default_signal_selection()
{
    return { { { Constellation::Gps, Band::L1, Code::CA }, {} }, { { Constellation::Galileo, Band::E1, Code::B }, {} } };
}

std::vector<Signal_id> all_signals()
{
    return { { Constellation::Gps, Band::L1, Code::CA },
             { Constellation::Gps, Band::L1, Code::Cd }, // GPS L1C
             { Constellation::Galileo, Band::E1, Code::B },
             { Constellation::Beidou, Band::B1, Code::I } };
}

int max_prn_for( Constellation c )
{
    // Cache: make_signal builds a Signal just to read its (constant) sv_range upper bound.
    static std::map<Constellation, int> cache;
    if( const auto it = cache.find( c ); it != cache.end() )
    {
        return it->second;
    }
    int top = 0;
    for( const Signal_id& id : all_signals() )
    {
        if( id.constellation == c )
        {
            top = make_signal( id.constellation, id.band, id.code )->sv_range().second;
            break;
        }
    }
    cache[c] = top;
    return top;
}

Signal_id signal_from_tokens( const std::string& constellation, const std::string& signal )
{
    const std::string c = to_lower( trim( constellation ) );
    const std::string s = to_lower( trim( signal ) );

    Constellation con;
    Signal_id     def; // the constellation's default (currently only) signal
    if( c == "gps" )
    {
        con = Constellation::Gps;
        def = { Constellation::Gps, Band::L1, Code::CA };
    }
    else if( c == "galileo" || c == "gal" )
    {
        con = Constellation::Galileo;
        def = { Constellation::Galileo, Band::E1, Code::B };
    }
    else if( c == "beidou" || c == "bds" )
    {
        con = Constellation::Beidou;
        def = { Constellation::Beidou, Band::B1, Code::I };
    }
    else
    {
        throw std::invalid_argument( "Unknown constellation '" + constellation + "' (expected gps, galileo, beidou)" );
    }

    if( s.empty() )
    {
        return def; // bare constellation -> its default signal (the .constellation is what callers use)
    }
    // Resolve the named component within the constellation.
    if( con == Constellation::Gps )
    {
        if( s == "l1ca" )
            return { Constellation::Gps, Band::L1, Code::CA };
        if( s == "l1c" )
            return { Constellation::Gps, Band::L1, Code::Cd };
        throw std::invalid_argument( "Unknown signal '" + signal + "' for gps (expected l1ca, l1c)" );
    }
    if( con == Constellation::Galileo && s == "e1" )
        return def;
    if( con == Constellation::Beidou && s == "b1i" )
        return def;
    throw std::invalid_argument( "Unknown signal '" + signal + "' for " + c + " (expected " + signal_name( def ) + ")" );
}

const char* signal_token( const Signal_id& id )
{
    switch( id.constellation )
    {
    case Constellation::Gps:
        return "gps";
    case Constellation::Galileo:
        return "galileo";
    case Constellation::Beidou:
        return "beidou";
    default:
        return "?";
    }
}

const char* signal_name( const Signal_id& id )
{
    switch( id.code ) // the component (a constellation can now have several)
    {
    case Code::CA:
        return "l1ca";
    case Code::Cd:
        return "l1c";
    case Code::B:
        return "e1"; // whole-signal token: the receiver tracks the E1-C pilot + decodes E1-B data as one signal
    case Code::I:
        return "b1i";
    default:
        return "?";
    }
}

std::vector<Signal_selection>
parse_selection( const std::vector<std::string>& signal_tokens, const std::vector<std::string>& prn_tokens )
{
    // 1. Which components to search: each token is CONSTELLATION[:COMPONENT]. A bare CONSTELLATION (no
    //    component) selects ALL of that constellation's components; an explicit one selects just it.
    //    De-duplicated, order kept.
    std::vector<Signal_id> comps;
    auto                   add = [&comps]( const Signal_id& id ) {
        if( std::find( comps.begin(), comps.end(), id ) == comps.end() )
        {
            comps.push_back( id );
        }
    };
    for( const std::string& raw : signal_tokens )
    {
        const std::string tok = trim( raw );
        if( tok.empty() )
        {
            continue;
        }
        const size_t      colon = tok.find( ':' );
        const std::string con   = ( colon == std::string::npos ) ? tok : tok.substr( 0, colon );
        const std::string comp  = ( colon == std::string::npos ) ? "" : tok.substr( colon + 1 );
        if( comp.empty() )
        {
            // Bare constellation -> all its components. (signal_from_tokens validates the token.)
            const Constellation c = signal_from_tokens( con ).constellation;
            for( const Signal_id& s : all_signals() )
            {
                if( s.constellation == c )
                {
                    add( s );
                }
            }
        }
        else
        {
            add( signal_from_tokens( con, comp ) );
        }
    }

    // 2. Per-constellation PRN allowlists. Tokens were comma-split by the CLI, so an alpha-headed token
    //    (CONSTELLATION[:firstPRN]) starts a constellation and numeric/range tokens attach to it.
    std::map<Constellation, std::set<int>> prn_by_con;
    bool                                   have_con = false;
    Constellation                          cur_con  = Constellation::Unknown;
    for( const std::string& raw : prn_tokens )
    {
        const std::string tok = trim( raw );
        if( tok.empty() )
        {
            continue;
        }
        const size_t      colon = tok.find( ':' );
        const std::string head  = ( colon == std::string::npos ) ? tok : tok.substr( 0, colon );
        if( has_alpha( head ) )
        {
            cur_con  = signal_from_tokens( head ).constellation;
            have_con = true;
            prn_by_con.try_emplace( cur_con ); // ensure the key exists even with no PRNs yet
            if( colon != std::string::npos )
            {
                for( const std::string& p : split( tok.substr( colon + 1 ), ',' ) )
                {
                    add_prn_token( prn_by_con[cur_con], p, max_prn_for( cur_con ) );
                }
            }
        }
        else
        {
            if( !have_con )
            {
                throw std::invalid_argument( "PRN '" + tok + "' given before any constellation in --prns" );
            }
            add_prn_token( prn_by_con[cur_con], tok, max_prn_for( cur_con ) );
        }
    }

    // 3. Every --prns constellation must actually be searched.
    for( const auto& [con, prns] : prn_by_con )
    {
        const bool searched = std::any_of( comps.begin(), comps.end(), [&]( const Signal_id& id ) { return id.constellation == con; } );
        if( !searched )
        {
            throw std::invalid_argument(
                "PRNs given for '" + std::string( signal_token( { con, Band::L1, Code::CA } ) ) + "' but it is not in --signal"
            );
        }
    }

    // 4. Combine: each component gets its constellation's PRN set (empty = all).
    std::vector<Signal_selection> out;
    for( const Signal_id& id : comps )
    {
        const auto it = prn_by_con.find( id.constellation );
        out.push_back( { id, it != prn_by_con.end() ? it->second : std::set<int> {} } );
    }
    return out;
}

std::set<int> parse_prn_list( const std::string& spec, int max_prn )
{
    std::set<int> out;
    for( const std::string& tok : split( spec, ',' ) )
    {
        add_prn_token( out, tok, max_prn );
    }
    return out;
}

std::string format_prn_list( const std::set<int>& prns )
{
    std::string out;
    auto        flush = [&]( int lo, int hi ) {
        if( !out.empty() )
        {
            out += ',';
        }
        out += std::to_string( lo );
        if( hi != lo )
        {
            out += '-' + std::to_string( hi );
        }
    };
    bool have = false;
    int  lo = 0, hi = 0;
    for( int p : prns ) // std::set iterates ascending
    {
        if( have && p == hi + 1 )
        {
            hi = p; // extend the current run
        }
        else
        {
            if( have )
            {
                flush( lo, hi );
            }
            lo = hi = p;
            have    = true;
        }
    }
    if( have )
    {
        flush( lo, hi );
    }
    return out;
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "parse_prn_list_ranges_and_singletons", "[selection][prn]" )
{
    REQUIRE( parse_prn_list( "1-4,7,9-11,32" ) == std::set<int> { 1, 2, 3, 4, 7, 9, 10, 11, 32 } );
    REQUIRE( parse_prn_list( "  5 , 5 ,3 " ) == std::set<int> { 3, 5 } ); // dedup + whitespace
    REQUIRE( parse_prn_list( "" ).empty() );                              // empty -> all
    REQUIRE( parse_prn_list( "10-10" ) == std::set<int> { 10 } );         // single-element range
    REQUIRE_THROWS_AS( parse_prn_list( "4-1" ), std::invalid_argument );  // reversed range
    REQUIRE_THROWS_AS( parse_prn_list( "1,x" ), std::invalid_argument );  // non-integer
}

TEST_CASE( "parse_prn_list_open_ended_ranges", "[selection][prn]" )
{
    // "N-" means N..max_prn, "-N" means 1..N (max_prn>0 required for the open high end).
    REQUIRE( parse_prn_list( "30-", 32 ) == std::set<int> { 30, 31, 32 } );
    REQUIRE( parse_prn_list( "-3", 32 ) == std::set<int> { 1, 2, 3 } );
    REQUIRE( parse_prn_list( "1,28-", 32 ) == std::set<int> { 1, 28, 29, 30, 31, 32 } );
    REQUIRE( parse_prn_list( "-", 5 ) == std::set<int> { 1, 2, 3, 4, 5 } ); // both ends open = all
    REQUIRE_THROWS_AS( parse_prn_list( "30-" ), std::invalid_argument );    // open high end, no max known
}

TEST_CASE( "format_prn_list_collapses_runs", "[selection][prn]" )
{
    REQUIRE( format_prn_list( { 1, 2, 3, 4, 7, 9, 10, 11, 32 } ) == "1-4,7,9-11,32" );
    REQUIRE( format_prn_list( {} ).empty() );
    REQUIRE( format_prn_list( { 5 } ) == "5" );
    // Round-trips with parse_prn_list.
    for( const std::string& s : { "1-4,7,9-11,32", "1,3,5,7", "10-10", "2-2,4-9" } )
    {
        REQUIRE( parse_prn_list( format_prn_list( parse_prn_list( s ) ) ) == parse_prn_list( s ) );
    }
}

TEST_CASE( "prn_selected_empty_filter_is_all", "[selection][prn]" )
{
    REQUIRE( prn_selected( {}, 17 ) );
    const std::set<int> some { 1, 5, 9 };
    REQUIRE( prn_selected( some, 5 ) );
    REQUIRE_FALSE( prn_selected( some, 6 ) );
}

TEST_CASE( "parse_selection_components_and_per_constellation_prns", "[selection][signals]" )
{
    // --signal gps:l1ca --signal galileo:e1   --prns gps:1,4,6-8 --prns galileo:3,11
    // (the CLI comma-splits the --prns values, so they arrive flattened):
    const auto s = parse_selection(
        { "gps:l1ca", "galileo:e1" }, { "gps:1", "4", "6-8", "galileo:3", "11" }
    );
    REQUIRE( s.size() == 2 );
    REQUIRE( s[0].id == Signal_id { Constellation::Gps, Band::L1, Code::CA } );
    REQUIRE( s[0].prns == std::set<int> { 1, 4, 6, 7, 8 } );
    REQUIRE( s[1].id == Signal_id { Constellation::Galileo, Band::E1, Code::B } );
    REQUIRE( s[1].prns == std::set<int> { 3, 11 } );

    // Bare constellation -> ALL of its components. GPS now has TWO (L1 C/A + L1C), so "gps galileo"
    // expands to 3 selections; no PRNs anywhere -> all PRNs.
    const auto t = parse_selection( { "gps", "galileo" }, {} );
    REQUIRE( t.size() == 3 );
    REQUIRE( t[0].id == Signal_id { Constellation::Gps, Band::L1, Code::CA } );
    REQUIRE( t[1].id == Signal_id { Constellation::Gps, Band::L1, Code::Cd } );
    REQUIRE( t[2].id == Signal_id { Constellation::Galileo, Band::E1, Code::B } );
    REQUIRE( t[0].prns.empty() );

    // PRNs are per-CONSTELLATION: one --prns gps:... applies to ALL of gps's selected components.
    const auto u = parse_selection( { "gps:l1ca", "gps:l1c", "galileo:e1" }, { "gps:1,2,5-7" } );
    REQUIRE( u.size() == 3 );
    REQUIRE( u[0].prns == std::set<int> { 1, 2, 5, 6, 7 } ); // gps L1 C/A
    REQUIRE( u[1].prns == std::set<int> { 1, 2, 5, 6, 7 } ); // gps L1C (same per-constellation set)
    REQUIRE( u[2].prns.empty() );

    // Bare "gps" (=l1ca+l1c) + an explicit "gps:l1ca" de-dups to just the two GPS components.
    REQUIRE( parse_selection( { "gps", "gps:l1ca" }, {} ).size() == 2 );

    // Open-ended PRN range resolves against the constellation's max (GPS = 32).
    const auto o = parse_selection( { "gps" }, { "gps:30-" } );
    REQUIRE( o[0].prns == std::set<int> { 30, 31, 32 } );

    REQUIRE( std::string( signal_token( s[1].id ) ) == "galileo" );
    REQUIRE( std::string( signal_name( s[0].id ) ) == "l1ca" );
    REQUIRE( default_signal_selection().size() == 2 );
}

TEST_CASE( "parse_selection_errors", "[selection][signals]" )
{
    REQUIRE_THROWS_AS( parse_selection( { "glonass" }, {} ), std::invalid_argument );          // unknown constellation
    REQUIRE_THROWS_AS( parse_selection( { "gps:e1" }, {} ), std::invalid_argument );          // component not in constellation
    REQUIRE_THROWS_AS( parse_selection( { "gps" }, { "gps:4-1" } ), std::invalid_argument );   // reversed PRN range
    REQUIRE_THROWS_AS( parse_selection( { "gps" }, { "4" } ), std::invalid_argument );         // PRN before any constellation
    REQUIRE_THROWS_AS( parse_selection( { "gps" }, { "galileo:1,2" } ), std::invalid_argument ); // PRNs for un-searched constellation
    REQUIRE_THROWS_AS( signal_from_tokens( "gps", "l5" ), std::invalid_argument );             // unknown gps component
    REQUIRE( signal_from_tokens( "gps" ).code == Code::CA );                                   // default ok
}
#endif

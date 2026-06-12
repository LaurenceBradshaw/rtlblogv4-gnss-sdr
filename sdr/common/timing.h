#pragma once
#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>

#include "logging.h"

// Lightweight, thread-safe wall-clock profiling for the hot paths.
//
// Channels run on separate threads, so each bucket accumulates atomically. Buckets
// are a fixed enum (not a string map) so the per-call cost is one steady_clock pair
// + two relaxed atomic adds - negligible next to the work being measured.
//
// Usage:
//   {
//       timing::Scoped t( timing::Bucket::Acquisition );
//       acquisition_.integrate( block );
//   }
// and periodically: logging::log(Info, timing::report(wall_s));
namespace timing
{

enum class Bucket : size_t
{
    Acquisition,    // Acquisition_engine::integrate (FFT correlation, all Doppler bins)
    TrackCorrelate, // Costas_tracker::correlate_epoch (per-sample E/P/L correlation)
    TrackLoops,     // Costas_tracker::run_loops (cumsum + PLL/DLL filters)
    Count
};

inline constexpr std::array<const char*, static_cast<size_t>( Bucket::Count )> kNames = {
    "acquisition", "track.correlate", "track.loops"
};

struct Stat
{
    std::atomic<uint64_t> ns { 0 };
    std::atomic<uint64_t> calls { 0 };
};

inline std::array<Stat, static_cast<size_t>( Bucket::Count )>& registry()
{
    static std::array<Stat, static_cast<size_t>( Bucket::Count )> r;
    return r;
}

inline void add( Bucket b, uint64_t ns )
{
    Stat& s = registry()[static_cast<size_t>( b )];
    s.ns.fetch_add( ns, std::memory_order_relaxed );
    s.calls.fetch_add( 1, std::memory_order_relaxed );
}

// RAII scope timer - adds elapsed ns to the bucket on destruction.
class Scoped
{
public:
    explicit Scoped( Bucket b )
        : bucket_( b ),
          start_( std::chrono::steady_clock::now() )
    {
    }
    ~Scoped()
    {
        const auto end = std::chrono::steady_clock::now();
        add( bucket_, static_cast<uint64_t>( std::chrono::duration_cast<std::chrono::nanoseconds>( end - start_ ).count() ) );
    }
    Scoped( const Scoped& )            = delete;
    Scoped& operator=( const Scoped& ) = delete;

private:
    Bucket                                bucket_;
    std::chrono::steady_clock::time_point start_;
};

// Format cumulative stats. cores = total CPU-time in the bucket / wall time, i.e. the
// equivalent number of fully-busy cores (can exceed 1 - work is summed across threads).
inline std::string report( double wall_s )
{
    std::string out = "Timing (cumulative):";
    auto&       r   = registry();
    for( size_t i = 0; i < r.size(); i++ )
    {
        const uint64_t ns    = r[i].ns.load( std::memory_order_relaxed );
        const uint64_t calls = r[i].calls.load( std::memory_order_relaxed );
        const double   ms    = static_cast<double>( ns ) * 1e-6;
        const double   us    = calls ? static_cast<double>( ns ) * 1e-3 / static_cast<double>( calls ) : 0.0;
        const double   cores = wall_s > 0.0 ? ( static_cast<double>( ns ) * 1e-9 ) / wall_s : 0.0;
        out += fmt::format(
            "\n  {:<16} {:>10.1f} ms  {:>9} calls  {:>8.2f} us/call  {:>5.2f} cores", kNames[i], ms, calls, us, cores
        );
    }
    return out;
}

} // namespace timing

#include "tracking_history.h"
#include <algorithm>
#include <cmath>

namespace
{
constexpr double IQ_EXTENT_HEADROOM = 1.6;   // grid edge sits this far past the peak prompt magnitude
constexpr double IQ_EXTENT_DECAY    = 0.999; // per-fade shrink (~20 Hz) so the extent recovers SLOWLY
constexpr double FADE_PERIOD_S      = 0.05;  // fade the whole grid at most this often (stream time)
constexpr double FADE_TAU_S         = 0.5;   // heat e-folding time -> ~1-2 s visible trail

constexpr double FAST_DT  = 0.02; // Doppler fast-ring sample period (s) -> 50 Hz
constexpr size_t FAST_CAP = 600;  // ~12 s of fast samples
constexpr size_t SLOW_CAP = 1000; // ~16 min of per-second means
} // namespace

void Tracking_history::record( double time_s, double doppler_hz, double prompt_i, double prompt_q )
{
    std::lock_guard<std::mutex> lk( mu_ );

    // ---- I/Q heatmap ----
    // Auto-scale the grid extent: ratchet up to ~headroom x the PEAK prompt magnitude immediately, then
    // hold it (only a slow per-fade decay below). Ratcheting to the peak means data never exceeds it, so
    // nothing clamps at the edges.
    const double wanted = std::max( std::abs( prompt_i ), std::abs( prompt_q ) ) * IQ_EXTENT_HEADROOM;
    if( wanted > iq_extent_ )
    {
        iq_extent_ = static_cast<float>( wanted );
    }

    if( !iq_has_data_ )
    {
        last_fade_s_ = time_s;
        iq_has_data_ = true;
    }
    else if( time_s - last_fade_s_ >= FADE_PERIOD_S )
    {
        const float decay = static_cast<float>( std::exp( -( time_s - last_fade_s_ ) / FADE_TAU_S ) );
        for( float& c : iq_heat_ )
        {
            c *= decay;
        }
        iq_extent_   = static_cast<float>( iq_extent_ * IQ_EXTENT_DECAY );
        last_fade_s_ = time_s;
    }

    if( iq_extent_ > 0.0f )
    {
        const auto bin = [this]( double v ) {
            const int b = static_cast<int>( ( v / iq_extent_ + 1.0 ) * 0.5 * IQ_GRID );
            return std::clamp( b, 0, IQ_GRID - 1 );
        };
        iq_heat_[bin( prompt_q ) * IQ_GRID + bin( prompt_i )] += 1.0f;
    }

    // ---- Doppler: fast ring (downsampled to FAST_DT) ----
    if( !fast_started_ || time_s - last_fast_s_ >= FAST_DT )
    {
        fast_.emplace_back( time_s, static_cast<float>( doppler_hz ) );
        if( fast_.size() > FAST_CAP )
        {
            fast_.pop_front();
        }
        last_fast_s_  = time_s;
        fast_started_ = true;
    }

    // ---- Doppler: slow ring (one mean per whole second) ----
    const long sec = static_cast<long>( std::floor( time_s ) );
    if( sec != slow_sec_ )
    {
        if( slow_count_ > 0 )
        {
            slow_.emplace_back( slow_sec_ + 0.5, static_cast<float>( slow_sum_ / slow_count_ ) );
            if( slow_.size() > SLOW_CAP )
            {
                slow_.pop_front();
            }
        }
        slow_sec_   = sec;
        slow_sum_   = 0.0;
        slow_count_ = 0;
    }
    slow_sum_ += doppler_hz;
    ++slow_count_;
}

void Tracking_history::reset()
{
    std::lock_guard<std::mutex> lk( mu_ );
    iq_heat_.fill( 0.0f );
    iq_extent_   = 0.0f;
    last_fade_s_ = 0.0;
    iq_has_data_ = false;

    fast_.clear();
    last_fast_s_  = 0.0;
    fast_started_ = false;
    slow_.clear();
    slow_sec_   = -1;
    slow_sum_   = 0.0;
    slow_count_ = 0;
}

Tracking_history::Snapshot Tracking_history::snapshot() const
{
    std::lock_guard<std::mutex> lk( mu_ );

    Snapshot s;
    s.iq.heat     = iq_heat_;
    s.iq.extent   = iq_extent_;
    s.iq.has_data = iq_has_data_;

    const auto fill = []( const std::deque<std::pair<double, float>>& src, Doppler_series& dst ) {
        dst.t_s.reserve( src.size() );
        dst.hz.reserve( src.size() );
        for( const auto& [t, hz] : src )
        {
            dst.t_s.push_back( t );
            dst.hz.push_back( hz );
        }
    };
    fill( fast_, s.doppler.fast );
    fill( slow_, s.doppler.slow );

    return s;
}

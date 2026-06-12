#pragma once
#include <array>
#include <deque>
#include <mutex>
#include <utility>
#include <vector>

// Per-channel rolling history for the GUI graphs. Written by the tracking worker EVERY epoch (record(),
// cheap), read by the GUI via a locked snapshot copy. Recorded ALWAYS - so a graph opened later already
// has data - and reset() on re-acquisition (a new lock = a new history).
//
// Holds: (1) the prompt I/Q constellation heatmap; (2) Doppler over time as two tiers - a fast ring
// (recent, high rate) and a slow ring (per-second means, long span) so both time scales are cheap.
class Tracking_history
{
public:
    static constexpr int IQ_GRID = 128; // I/Q heatmap resolution (cells per axis)

    // I/Q heatmap, ready to render.
    struct Iq_snapshot
    {
        std::array<float, IQ_GRID * IQ_GRID> heat {};         // cell heat, [row = Q bin][col = I bin]
        float                                extent   = 0.0f; // grid spans [-extent, +extent] in I and Q
        bool                                 has_data = false;
    };

    // A time series of Doppler samples (parallel arrays for easy plotting).
    struct Doppler_series
    {
        std::vector<double> t_s; // channel stream time (s)
        std::vector<float>  hz;  // Doppler (Hz)
    };
    struct Doppler_snapshot
    {
        Doppler_series fast; // ~50 Hz, last ~12 s     (the "seconds" view)
        Doppler_series slow; // 1 Hz per-second means  (the "minutes" view)
    };

    // Everything a graph might want, copied under one lock.
    struct Snapshot
    {
        Iq_snapshot     iq;
        Doppler_snapshot doppler;
    };

    // Called by the owning worker each tracking epoch. time_s = the channel's stream time (s).
    void record( double time_s, double doppler_hz, double prompt_i, double prompt_q );
    // Clear all history (e.g. on re-acquisition).
    void reset();

    Snapshot snapshot() const;

private:
    mutable std::mutex mu_;

    // I/Q heatmap: per-cell heat, periodically time-faded; extent auto-scaled to the prompt magnitude.
    std::array<float, IQ_GRID * IQ_GRID> iq_heat_ {};
    float                                iq_extent_   = 0.0f;
    double                               last_fade_s_ = 0.0;
    bool                                 iq_has_data_ = false;

    // Doppler rings. fast_: every FAST_DT seconds, capped to FAST_CAP. slow_: one per-second mean,
    // capped to SLOW_CAP. (pair = (time_s, Doppler Hz)).
    std::deque<std::pair<double, float>> fast_;
    double                               last_fast_s_ = 0.0;
    bool                                 fast_started_ = false;

    std::deque<std::pair<double, float>> slow_;
    long                                 slow_sec_   = -1; // integer second currently being averaged
    double                               slow_sum_   = 0.0;
    int                                  slow_count_ = 0;
};

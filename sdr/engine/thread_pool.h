#pragma once
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <mutex>
#include <thread>
#include <vector>

// Fixed-size worker pool with two priority lanes.
//
// Workers always drain the High lane before the Low lane, so latency-critical work
// (channel tracking, ~90 us, must run every ms) is never stuck behind slow,
// deferrable work (acquisition, ~35 ms). Tasks are plain std::function<void()>.
//
// The High lane is ORDERED by the per-task `order` key (lowest first), not FIFO: passing each
// channel's next_sample makes the FURTHEST-BEHIND channel (least data processed) run first, so the
// slowest channels catch up instead of being served in arbitrary submission order - which tightens
// the spread between channels (and the run-to-run variance it drives). The Low lane (acquisition)
// stays FIFO and ignores `order`.
//
// Lifetime: the destructor lets workers finish any in-flight + queued tasks, then
// joins. Submit nothing after you intend to shut down (the owner - Scheduler - stops
// dispatching first, so the queue drains cleanly).
class Thread_pool
{
public:
    enum class Priority
    {
        High,
        Low
    };

    explicit Thread_pool( unsigned num_workers );
    ~Thread_pool();

    Thread_pool( const Thread_pool& )            = delete;
    Thread_pool& operator=( const Thread_pool& ) = delete;

    // order: High-lane ordering key (lowest runs first; pass the channel's next_sample). Ignored by
    // the Low lane.
    void submit( Priority priority, uint64_t order, std::function<void()> task );

    unsigned worker_count() const
    {
        return static_cast<unsigned>( workers_.size() );
    }

private:
    void worker_loop();

    std::vector<std::thread>                       workers_;
    std::multimap<uint64_t, std::function<void()>> high_; // ordered: lowest `order` (furthest behind) first
    std::deque<std::function<void()>>              low_;  // FIFO
    std::mutex                                     mutex_;
    std::condition_variable                        cv_;
    bool                                           stop_ = false;
};

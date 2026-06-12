#pragma once
#include <atomic>
#include <thread>
#include <vector>
#include "channel.h"
#include "sample_buffer.h"
#include "thread_pool.h"

// Drives channel processing. A single lightweight thread periodically scans the
// channels and hands the ready ones to the Thread_pool as work-items:
//   - tracking -> High priority (latency-critical, fast)
//   - acquisition -> Low priority (slow, deferrable). Acquisition may use every worker
//     NOT currently demanded by tracking: the per-scan budget is worker_count minus the
//     tracking channels needing a worker right now. So at startup (nothing tracking yet)
//     acquisition saturates the pool, and it backs off only as satellites start tracking.
// Each scan also advances the buffer's oldest index to the slowest channel.
//
// Replaces the old thread-per-channel model (32 OS threads busy-spinning).
class Scheduler
{
public:
    Scheduler( std::vector<Channel*> channels, Thread_pool& pool, Sample_buffer& sample_buffer );
    ~Scheduler();

    Scheduler( const Scheduler& )            = delete;
    Scheduler& operator=( const Scheduler& ) = delete;

    // Slowest channel's stream position (samples), updated each scan. For status.
    Sample_index min_next_sample() const
    {
        return min_next_.load( std::memory_order_relaxed );
    }

private:
    void dispatch_loop();

    std::vector<Channel*> channels_;
    Thread_pool&          pool_;
    Sample_buffer&        sample_buffer_;

    std::atomic<unsigned>     acq_inflight_ { 0 }; // concurrent acquisition tasks in flight
    std::atomic<Sample_index> min_next_ { 0 };

    std::atomic<bool> stop_ { false };
    std::thread       thread_; // declared last - starts after everything is ready
};

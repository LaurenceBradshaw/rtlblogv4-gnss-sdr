#include "scheduler.h"
#include <algorithm>
#include <limits>
#include <utility>

namespace
{
// Scan cadence. Readiness is wall-clock based (the Sample_buffer real-time gate), so
// a timer - not a push signal - is the right trigger. A few ms of latency is nothing
// for a 1 ms epoch rate, and the gate keeps everything bounded.
constexpr auto SCAN_INTERVAL = std::chrono::milliseconds( 1 );
} // namespace

Scheduler::Scheduler( std::vector<Channel*> channels, Thread_pool& pool, Sample_buffer& sample_buffer )
    : channels_( std::move( channels ) ),
      pool_( pool ),
      sample_buffer_( sample_buffer ),
      thread_( &Scheduler::dispatch_loop, this )
{
}

Scheduler::~Scheduler()
{
    stop_.store( true, std::memory_order_relaxed );
    if( thread_.joinable() )
    {
        thread_.join();
    }
}

void Scheduler::dispatch_loop()
{
    while( !stop_.load( std::memory_order_relaxed ) )
    {
        std::this_thread::sleep_for( SCAN_INTERVAL );

        Sample_index min_next = std::numeric_limits<Sample_index>::max();

        // Pass 1: reclaim floor, and count the workers tracking needs RIGHT NOW (channels
        // tracking and either already running or with pending work this scan). Acquisition
        // is then allowed every other worker - see acq_budget below.
        unsigned tracking_demand = 0;
        for( Channel* ch : channels_ )
        {
            // Single-owner flag. acquire so a settled channel's state/next_sample are
            // visible (published by its last worker's release).
            const bool sched = ch->scheduled().load( std::memory_order_acquire );

            // Buffer reclaim floor: while a channel is RUNNING its state is in flux
            // (it may be transitioning acquire->track), so conservatively pin
            // next_sample (the lowest sample it could read; only moves forward within
            // a quantum). When settled, trust buffer_floor() - an acquiring channel
            // jumps to the newest block, so it must not pin the buffer at a stale index.
            min_next = std::min( min_next, sched ? ch->next_sample() : ch->buffer_floor() );

            if( !ch->is_acquiring() && ( sched || ch->has_pending_work() ) )
            {
                ++tracking_demand;
            }
        }

        // Acquisition may use every worker not currently demanded by tracking. Until sats
        // start tracking this is the whole pool (fast startup acquisition); it tightens as
        // tracking comes online. The High-priority lane still lets a tracking task that
        // becomes ready mid-scan jump ahead of queued acquisitions for the next free worker.
        const unsigned workers    = pool_.worker_count();
        const unsigned acq_budget = std::max( 1u, workers > tracking_demand ? workers - tracking_demand : 1u );

        // Pass 2: dispatch ready channels.
        for( Channel* ch : channels_ )
        {
            const bool sched = ch->scheduled().load( std::memory_order_acquire );
            if( sched )
            {
                continue; // already queued or running
            }
            if( !ch->has_pending_work() )
            {
                continue;
            }

            const bool acquiring = ch->is_acquiring();
            if( acquiring )
            {
                // Dynamic budget: only this thread increments, so the load->add window
                // can't overshoot.
                if( acq_inflight_.load( std::memory_order_relaxed ) >= acq_budget )
                {
                    continue;
                }
                acq_inflight_.fetch_add( 1, std::memory_order_relaxed );
                // Publish the jump target now so this channel doesn't pin a stale
                // (backed-off) next_sample once it's scheduled. Safe: we are its sole
                // accessor until scheduled_ is released below.
                ch->refresh_acquire_floor();
            }

            ch->scheduled().store( true, std::memory_order_release );

            pool_.submit(
                acquiring ? Thread_pool::Priority::Low : Thread_pool::Priority::High,
                // High-lane order key: the channel's stream position, so the furthest-behind
                // (lowest next_sample) tracking channel is served first. Ignored for the Low lane.
                static_cast<uint64_t>( ch->next_sample() ),
                [this, ch, acquiring]
                {
                    ch->process();
                    if( acquiring )
                    {
                        acq_inflight_.fetch_sub( 1, std::memory_order_relaxed );
                    }
                    ch->scheduled().store( false, std::memory_order_release );
                }
            );
        }

        if( min_next != std::numeric_limits<Sample_index>::max() )
        {
            sample_buffer_.advance_oldest( min_next );
            min_next_.store( min_next, std::memory_order_relaxed );
        }
    }
}

#include "thread_pool.h"
#include <exception>
#include "logging.h"

Thread_pool::Thread_pool( unsigned num_workers )
{
    if( num_workers == 0 )
    {
        num_workers = 1;
    }
    workers_.reserve( num_workers );

    for( unsigned i = 0; i < num_workers; i++ )
    {
        workers_.emplace_back( [this] { worker_loop(); } );
    }
}

Thread_pool::~Thread_pool()
{
    {
        std::lock_guard<std::mutex> lock( mutex_ );
        stop_ = true;
    }
    cv_.notify_all();
    for( auto& w : workers_ )
    {
        if( w.joinable() )
        {
            w.join();
        }
    }
}

void Thread_pool::submit( Priority priority, uint64_t order, std::function<void()> task )
{
    {
        std::lock_guard<std::mutex> lock( mutex_ );
        if( priority == Priority::High )
        {
            high_.emplace( order, std::move( task ) ); // ordered by `order`: furthest-behind first
        }
        else
        {
            low_.push_back( std::move( task ) ); // FIFO
        }
    }
    cv_.notify_one();
}

void Thread_pool::worker_loop()
{
    for( ;; )
    {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> lock( mutex_ );
            cv_.wait( lock, [this] { return stop_ || !high_.empty() || !low_.empty(); } );

            if( high_.empty() && low_.empty() )
            {
                return; // woken only by stop_ with nothing left -> exit
            }

            // High lane first, taking the lowest-order (furthest-behind) task; else FIFO from Low.
            if( !high_.empty() )
            {
                auto it = high_.begin();
                task    = std::move( it->second );
                high_.erase( it );
            }
            else
            {
                task = std::move( low_.front() );
                low_.pop_front();
            }
        }

        // A throwing task must not take down the worker.
        try
        {
            task();
        }
        catch( const std::exception& e )
        {
            logging::log( logging::Level::Error, std::string( "Thread_pool task threw: " ) + e.what() );
        }
        catch( ... )
        {
            logging::log( logging::Level::Error, "Thread_pool task threw unknown exception" );
        }
    }
}

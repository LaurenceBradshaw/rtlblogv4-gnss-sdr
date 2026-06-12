#include "sample_buffer.h"
#include <cassert>
#include <chrono>
#include <cstring>
#include <stdexcept>
#include <thread>

Sample_buffer::Sample_buffer( size_t capacity, double sample_rate_hz )
    : capacity_( capacity ),
      sample_rate_hz_( sample_rate_hz ),
      buffer_( capacity )
{
    assert( capacity > 0 && ( capacity & ( capacity - 1 ) ) == 0 );
}

// Called from the streaming thread.
//
// Uses two memcpy calls (<= one wrap) instead of a per-sample loop.
// The old per-sample approach computed (base + i) % capacity_ on every
// iteration, which generated a full integer divide on each sample because
// capacity_ is a runtime value - very slow for 8192 samples per call.
//
// capacity_ is always a power of 2 so (index & (capacity_-1)) == (index % capacity_).
//
// Release-store on write_index_ guarantees channel threads that acquire-load
// it will see all the buffer_ writes completed above.
void Sample_buffer::push( const Complex_buf& samples )
{
    const size_t n = samples.size();
    if( n == 0 )
    {
        return;
    }

    // Record stream start on the first push so that is_available() can gate
    // samples to real time.  Written before the release-store so that any
    // thread that acquire-loads stream_started_=true sees stream_start_.
    if( !stream_started_.load( std::memory_order_relaxed ) )
    {
        stream_start_ = std::chrono::steady_clock::now();
        stream_started_.store( true, std::memory_order_release );
    }

    // Block until there is space. The buffer holds >= 2 s of samples so sleeping
    // 1 s is always enough time for channels to drain at least half of it.
    // A simple sleep avoids the CPU waste of a tight spin.
    while( write_index_.load( std::memory_order_relaxed ) - oldest_valid_index_.load( std::memory_order_relaxed ) +
               static_cast<Sample_index>( n ) >
           static_cast<Sample_index>( capacity_ ) )
    {
        std::this_thread::sleep_for( std::chrono::seconds( 1 ) );
    }

    const Sample_index base         = write_index_.load( std::memory_order_relaxed );
    const size_t       mask         = capacity_ - 1;
    const size_t       start_offset = static_cast<size_t>( base & mask );
    const size_t       space_to_end = capacity_ - start_offset;

    if( n <= space_to_end )
    {
        std::memcpy( buffer_.data() + start_offset, samples.data(), n * sizeof( Complex_sample ) );
    }
    else
    {
        std::memcpy( buffer_.data() + start_offset, samples.data(), space_to_end * sizeof( Complex_sample ) );
        std::memcpy( buffer_.data(), samples.data() + space_to_end, ( n - space_to_end ) * sizeof( Complex_sample ) );
    }

    write_index_.store( base + n, std::memory_order_release );
}

size_t Sample_buffer::index_to_offset( Sample_index index ) const
{
    return static_cast<size_t>( index & ( capacity_ - 1 ) );
}

const Complex_sample& Sample_buffer::at( Sample_index index ) const
{
    if( !is_available( index ) )
    {
        throw std::runtime_error( "Sample_buffer underrun" );
    }

    return buffer_[index_to_offset( index )];
}

Sample_block Sample_buffer::at( Sample_index start, size_t n ) const
{
    const Sample_index end = start + n;

    if( !is_valid_range( start, end ) )
    {
        throw std::runtime_error( "Sample_buffer underrun" );
    }

    const size_t start_offset = index_to_offset( start );
    const size_t end_offset   = index_to_offset( end );

    if( end_offset > start_offset )
    {
        return { &buffer_[start_offset], n, nullptr, 0 };
    }

    const size_t first_len  = buffer_.size() - start_offset;
    const size_t second_len = end_offset;
    return { &buffer_[start_offset], first_len, &buffer_[0], second_len };
}

bool Sample_buffer::is_available( Sample_index index ) const
{
    const Sample_index w = write_index_.load( std::memory_order_acquire );
    const Sample_index o = oldest_valid_index_.load( std::memory_order_relaxed );

    if( index < o || index >= w )
    {
        return false;
    }

    // Real-time gate: sample i is visible only once i/sample_rate_hz seconds
    // have elapsed since the first push().  This ensures channels process data
    // at the rate hardware would deliver it, even though the ring buffer may
    // already hold samples far ahead in time.
    //
    // stream_started_ is acquire-loaded so that stream_start_ (written before
    // the release-store in push()) is guaranteed visible here.
    if( stream_started_.load( std::memory_order_acquire ) )
    {
        const auto elapsed = std::chrono::steady_clock::now() - stream_start_;
        const auto elapsed_samples =
            static_cast<Sample_index>( std::chrono::duration<double>( elapsed ).count() * sample_rate_hz_ );

        if( index >= elapsed_samples )
            return false;
    }

    return true;
}

Sample_index Sample_buffer::write_index() const
{
    return write_index_.load( std::memory_order_acquire );
}

Sample_index Sample_buffer::newest_valid_index() const
{
    const Sample_index w = write_index_.load( std::memory_order_acquire );
    const Sample_index o = oldest_valid_index_.load( std::memory_order_relaxed );

    if( w == 0 || w <= o )
    {
        return o; // no data has been written yet
    }

    // Start from the last written sample
    Sample_index newest = w - 1;

    // Apply the same real-time gate used by is_available():
    // sample i is only visible once i/sample_rate_hz seconds have elapsed.
    if( stream_started_.load( std::memory_order_acquire ) )
    {
        const auto elapsed = std::chrono::steady_clock::now() - stream_start_;
        const auto elapsed_samples =
            static_cast<Sample_index>( std::chrono::duration<double>( elapsed ).count() * sample_rate_hz_ );

        if( elapsed_samples == 0 )
        {
            return o; // real-time gate hasn't opened yet
        }

        newest = std::min( newest, elapsed_samples - 1 );
    }

    return newest < o ? o : newest;
}

Sample_index Sample_buffer::oldest_valid_index() const
{
    return oldest_valid_index_.load( std::memory_order_relaxed );
}

void Sample_buffer::advance_oldest( Sample_index index )
{
    Sample_index current = oldest_valid_index_.load( std::memory_order_relaxed );
    while( index > current )
    {
        if( oldest_valid_index_.compare_exchange_weak( current, index, std::memory_order_relaxed ) )
        {
            break;
        }
    }
}

bool Sample_buffer::is_valid_range( Sample_index start, Sample_index end ) const
{
    const Sample_index w = write_index_.load( std::memory_order_acquire );
    const Sample_index o = oldest_valid_index_.load( std::memory_order_relaxed );
    return start >= o && end <= w;
}

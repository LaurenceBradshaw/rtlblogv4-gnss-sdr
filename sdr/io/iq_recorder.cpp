#include "iq_recorder.h"
#include <cctype>
#include <cstring>
#include <stdexcept>
#include "logging.h"

// Recording dumps the raw in-memory sample layout; the file is only replayable as float32 if that holds.
static_assert( sizeof( Complex_sample ) == 2 * sizeof( float ), "Complex_sample must be interleaved float32 I/Q" );

Iq_recorder::Iq_recorder( const std::string& path ) : out_( path, std::ios::binary | std::ios::trunc )
{
    if( !out_ )
    {
        throw std::runtime_error( "Could not open IQ record file: " + path );
    }
    writer_ = std::thread( [this] { writer_loop(); } );
}

Iq_recorder::~Iq_recorder()
{
    {
        std::lock_guard<std::mutex> lock( mtx_ );
        stop_ = true;
    }
    cv_.notify_one();
    if( writer_.joinable() )
    {
        writer_.join(); // the writer drains the queue before returning, so nothing recorded is lost
    }
    out_.flush();
    if( dropped_blocks_ > 0 )
    {
        logging::log(
            logging::Level::Warning,
            fmt::format( "IQ recorder dropped {} block(s) - disk could not keep up with the stream", dropped_blocks_ )
        );
    }
}

void Iq_recorder::record( const Complex_buf& samples )
{
    if( samples.empty() )
    {
        return;
    }
    const size_t bytes = samples.size() * sizeof( Complex_sample );
    {
        std::lock_guard<std::mutex> lock( mtx_ );
        if( queued_bytes_ + bytes > MAX_QUEUE_BYTES )
        {
            ++dropped_blocks_; // bounded: degrade the recording, never block the live pipeline
            return;
        }
        queue_.push_back( samples ); // copy: the writer thread outlives this call
        queued_bytes_ += bytes;
    }
    cv_.notify_one();
}

std::string with_iq_extension( const std::string& path )
{
    if( path.empty() )
    {
        return path; // empty = no recording; leave it
    }
    const auto ends_with = [&]( const char* ext )
    {
        const size_t n = std::strlen( ext );
        if( path.size() < n )
        {
            return false;
        }
        for( size_t i = 0; i < n; ++i )
        {
            if( std::tolower( static_cast<unsigned char>( path[path.size() - n + i] ) ) != ext[i] )
            {
                return false;
            }
        }
        return true;
    };
    if( ends_with( ".f32" ) || ends_with( ".iq" ) )
    {
        return path;
    }
    return path + ".f32";
}

void Iq_recorder::writer_loop()
{
    for( ;; )
    {
        Complex_buf block;
        {
            std::unique_lock<std::mutex> lock( mtx_ );
            cv_.wait( lock, [this] { return stop_ || !queue_.empty(); } );
            if( queue_.empty() )
            {
                return; // woken to stop with an empty queue -> done
            }
            block = std::move( queue_.front() );
            queue_.pop_front();
            queued_bytes_ -= block.size() * sizeof( Complex_sample );
        }
        out_.write(
            reinterpret_cast<const char*>( block.data() ),
            static_cast<std::streamsize>( block.size() * sizeof( Complex_sample ) )
        );
    }
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_test_macros.hpp>

TEST_CASE( "with_iq_extension_defaults_and_preserves", "[io][record]" )
{
    REQUIRE( with_iq_extension( "" ) == "" );                       // empty stays empty (= no recording)
    REQUIRE( with_iq_extension( "capture" ) == "capture.f32" );     // no extension -> default .f32
    REQUIRE( with_iq_extension( "/tmp/run1" ) == "/tmp/run1.f32" );
    REQUIRE( with_iq_extension( "my.capture" ) == "my.capture.f32" ); // a dot that isn't a known ext
    REQUIRE( with_iq_extension( "x.f32" ) == "x.f32" );            // known extensions preserved
    REQUIRE( with_iq_extension( "x.iq" ) == "x.iq" );
    REQUIRE( with_iq_extension( "X.F32" ) == "X.F32" );            // case-insensitive
    REQUIRE( with_iq_extension( "data.IQ" ) == "data.IQ" );
}
#endif

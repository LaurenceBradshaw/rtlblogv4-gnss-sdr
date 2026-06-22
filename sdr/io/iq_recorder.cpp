#include "iq_recorder.h"
#include <cctype>
#include <cstring>
#include <stdexcept>
#include "logging.h"

Iq_recorder::Iq_recorder( const std::string& path, Iq_sample_format format )
    : out_( path, std::ios::binary | std::ios::trunc ),
      format_( format )
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
    std::vector<uint8_t> block = encode_iq( samples, format_ ); // float -> on-disk format (off the writer thread)
    {
        std::lock_guard<std::mutex> lock( mtx_ );
        if( queued_bytes_ + block.size() > MAX_QUEUE_BYTES )
        {
            ++dropped_blocks_; // bounded: degrade the recording, never block the live pipeline
            return;
        }
        queued_bytes_ += block.size();
        queue_.push_back( std::move( block ) );
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
        std::vector<uint8_t> block;
        {
            std::unique_lock<std::mutex> lock( mtx_ );
            cv_.wait( lock, [this] { return stop_ || !queue_.empty(); } );
            if( queue_.empty() )
            {
                return; // woken to stop with an empty queue -> done
            }
            block = std::move( queue_.front() );
            queue_.pop_front();
            queued_bytes_ -= block.size();
        }
        out_.write( reinterpret_cast<const char*>( block.data() ), static_cast<std::streamsize>( block.size() ) );
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

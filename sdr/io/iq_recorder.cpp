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

namespace
{
// The format-tagging extensions (with dot), one per Iq_sample_format - see iq_format_extension.
const char* const MANAGED_EXTS[] = { ".f32", ".i8", ".ui8", ".i16", ".ui16" };

bool ends_with_ci( const std::string& s, const char* suffix )
{
    const size_t n = std::strlen( suffix );
    if( s.size() < n )
    {
        return false;
    }
    for( size_t i = 0; i < n; ++i )
    {
        if( std::tolower( static_cast<unsigned char>( s[s.size() - n + i] ) ) != suffix[i] )
        {
            return false;
        }
    }
    return true;
}

// Length (incl. dot) of a trailing managed extension, or 0 if none.
size_t managed_ext_len( const std::string& s )
{
    for( const char* ext : MANAGED_EXTS )
    {
        if( ends_with_ci( s, ext ) )
        {
            return std::strlen( ext );
        }
    }
    return 0;
}

// True if the final path component has any '.'-extension at all.
bool has_extension( const std::string& s )
{
    const size_t slash = s.find_last_of( "/\\" );
    const size_t dot   = s.find_last_of( '.' );
    return dot != std::string::npos && ( slash == std::string::npos || dot > slash ) && dot + 1 < s.size();
}
} // namespace

std::string with_iq_extension( const std::string& path, Iq_sample_format format )
{
    if( path.empty() )
    {
        return path; // empty = no recording; leave it
    }
    const std::string want = std::string( "." ) + iq_format_extension( format );
    if( const size_t n = managed_ext_len( path ) )
    {
        return path.substr( 0, path.size() - n ) + want; // managed extension -> retag to this format
    }
    if( has_extension( path ) )
    {
        return path; // a deliberate custom extension -> leave it
    }
    return path + want; // no extension -> append the format's
}

std::string swap_iq_extension( const std::string& path, Iq_sample_format format )
{
    if( const size_t n = managed_ext_len( path ) )
    {
        return path.substr( 0, path.size() - n ) + "." + iq_format_extension( format );
    }
    return path; // not a managed extension (custom, or none) -> leave it
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

TEST_CASE( "with_iq_extension_tags_by_format", "[io][record]" )
{
    using F = Iq_sample_format;
    REQUIRE( with_iq_extension( "", F::INT8 ) == "" );                      // empty stays empty (= no recording)
    REQUIRE( with_iq_extension( "capture", F::INT8 ) == "capture.i8" );     // no extension -> append the format's
    REQUIRE( with_iq_extension( "/tmp/run1", F::UINT8 ) == "/tmp/run1.ui8" );
    REQUIRE( with_iq_extension( "capture.f32", F::INT8 ) == "capture.i8" ); // managed extension -> retag to format
    REQUIRE( with_iq_extension( "x.ui8", F::FLOAT32 ) == "x.f32" );
    REQUIRE( with_iq_extension( "X.I16", F::UINT16 ) == "X.ui16" );         // retag is case-insensitive
    REQUIRE( with_iq_extension( "cap.dat", F::INT8 ) == "cap.dat" );        // custom extension -> left alone
    REQUIRE( with_iq_extension( "my.capture", F::INT8 ) == "my.capture" );  // .capture is custom -> left
}

TEST_CASE( "swap_iq_extension_soft_forces_managed_only", "[io][record]" )
{
    using F = Iq_sample_format;
    REQUIRE( swap_iq_extension( "cap.i8", F::UINT8 ) == "cap.ui8" );  // managed -> retag
    REQUIRE( swap_iq_extension( "cap.f32", F::INT16 ) == "cap.i16" );
    REQUIRE( swap_iq_extension( "cap", F::UINT8 ) == "cap" );         // no extension -> leave (soft force only)
    REQUIRE( swap_iq_extension( "cap.dat", F::UINT8 ) == "cap.dat" ); // custom -> leave
}
#endif

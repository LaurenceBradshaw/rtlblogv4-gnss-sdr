#include "iq_format.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace
{
template <typename T>
T read_trivial( const unsigned char* data )
{
    T value {};
    std::memcpy( &value, data, sizeof( T ) );
    return value;
}

// Branchless normalisers - multiply instead of divide, clamp instead of if/return chains. This lets the
// compiler auto-vectorise the decode loop (the branchy originals blocked vectorisation).
constexpr float INV_127     = 1.0f / 127.0f;
constexpr float INV_127_5   = 1.0f / 127.5f;
constexpr float INV_32767   = 1.0f / 32767.0f;
constexpr float INV_32767_5 = 1.0f / 32767.5f;

inline float normalise_int8( int8_t value )
{
    return std::clamp( static_cast<float>( value ) * INV_127, -1.0f, 1.0f );
}
inline float normalise_int16( int16_t value )
{
    return std::clamp( static_cast<float>( value ) * INV_32767, -1.0f, 1.0f );
}
inline float normalise_uint8( uint8_t value )
{
    return ( static_cast<float>( value ) - 127.5f ) * INV_127_5;
}
inline float normalise_uint16( uint16_t value )
{
    return ( static_cast<float>( value ) - 32767.5f ) * INV_32767_5;
}

// Inverses of the normalisers: float [-1,1] -> the on-disk integer (round + clamp to the type range).
inline int8_t denormalise_int8( float f )
{
    return static_cast<int8_t>( std::clamp( std::lround( f * 127.0f ), -128L, 127L ) );
}
inline int16_t denormalise_int16( float f )
{
    return static_cast<int16_t>( std::clamp( std::lround( f * 32767.0f ), -32768L, 32767L ) );
}
inline uint8_t denormalise_uint8( float f )
{
    return static_cast<uint8_t>( std::clamp( std::lround( f * 127.5f + 127.5f ), 0L, 255L ) );
}
inline uint16_t denormalise_uint16( float f )
{
    return static_cast<uint16_t>( std::clamp( std::lround( f * 32767.5f + 32767.5f ), 0L, 65535L ) );
}
} // namespace

int bytes_per_iq_sample( Iq_sample_format format )
{
    switch( format )
    {
    case Iq_sample_format::UINT8:
    case Iq_sample_format::INT8:
        return 2;
    case Iq_sample_format::UINT16:
    case Iq_sample_format::INT16:
        return 4;
    case Iq_sample_format::FLOAT32:
        return 8;
    }
    throw std::runtime_error( "Unsupported IQ sample format" );
}

const char* iq_format_name( Iq_sample_format format )
{
    switch( format )
    {
    case Iq_sample_format::UINT8:
        return "uint8";
    case Iq_sample_format::INT8:
        return "int8";
    case Iq_sample_format::UINT16:
        return "uint16";
    case Iq_sample_format::INT16:
        return "int16";
    case Iq_sample_format::FLOAT32:
        return "float32";
    }
    return "?";
}

void decode_iq( const unsigned char* buf, size_t len, Iq_sample_format format, Complex_buf& output )
{
    const size_t bps       = static_cast<size_t>( bytes_per_iq_sample( format ) );
    const size_t n_samples = len / bps;
    for( size_t i = 0; i < n_samples; ++i )
    {
        const unsigned char* p = buf + i * bps;
        switch( format )
        {
        case Iq_sample_format::UINT8:
            output[i] = Complex_sample( normalise_uint8( p[0] ), normalise_uint8( p[1] ) );
            break;
        case Iq_sample_format::INT8:
            output[i] = Complex_sample(
                normalise_int8( static_cast<int8_t>( p[0] ) ), normalise_int8( static_cast<int8_t>( p[1] ) )
            );
            break;
        case Iq_sample_format::UINT16:
            output[i] = Complex_sample(
                normalise_uint16( read_trivial<uint16_t>( p ) ),
                normalise_uint16( read_trivial<uint16_t>( p + sizeof( uint16_t ) ) )
            );
            break;
        case Iq_sample_format::INT16:
            output[i] = Complex_sample(
                normalise_int16( read_trivial<int16_t>( p ) ),
                normalise_int16( read_trivial<int16_t>( p + sizeof( int16_t ) ) )
            );
            break;
        case Iq_sample_format::FLOAT32:
            output[i] =
                Complex_sample( read_trivial<float>( p ), read_trivial<float>( p + sizeof( float ) ) );
            break;
        }
    }
}

std::vector<uint8_t> encode_iq( const Complex_buf& samples, Iq_sample_format format )
{
    const size_t         bps = static_cast<size_t>( bytes_per_iq_sample( format ) );
    std::vector<uint8_t> out( samples.size() * bps );
    for( size_t i = 0; i < samples.size(); ++i )
    {
        unsigned char* p  = out.data() + i * bps;
        const float    re = samples[i].real();
        const float    im = samples[i].imag();
        switch( format )
        {
        case Iq_sample_format::INT8:
            p[0] = static_cast<uint8_t>( denormalise_int8( re ) );
            p[1] = static_cast<uint8_t>( denormalise_int8( im ) );
            break;
        case Iq_sample_format::UINT8:
            p[0] = denormalise_uint8( re );
            p[1] = denormalise_uint8( im );
            break;
        case Iq_sample_format::INT16:
        {
            const int16_t qi = denormalise_int16( re ), qq = denormalise_int16( im );
            std::memcpy( p, &qi, sizeof( int16_t ) );
            std::memcpy( p + sizeof( int16_t ), &qq, sizeof( int16_t ) );
            break;
        }
        case Iq_sample_format::UINT16:
        {
            const uint16_t qi = denormalise_uint16( re ), qq = denormalise_uint16( im );
            std::memcpy( p, &qi, sizeof( uint16_t ) );
            std::memcpy( p + sizeof( uint16_t ), &qq, sizeof( uint16_t ) );
            break;
        }
        case Iq_sample_format::FLOAT32:
            std::memcpy( p, &re, sizeof( float ) );
            std::memcpy( p + sizeof( float ), &im, sizeof( float ) );
            break;
        }
    }
    return out;
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace
{
const Iq_sample_format ALL_FORMATS[] = { Iq_sample_format::UINT8,  Iq_sample_format::INT8, Iq_sample_format::UINT16,
                                         Iq_sample_format::INT16,  Iq_sample_format::FLOAT32 };

// Worst-case round-trip error for a format = one quantisation step (≈ 1/full-scale); float32 is exact.
double round_trip_tol( Iq_sample_format f )
{
    switch( f )
    {
    case Iq_sample_format::INT8:
    case Iq_sample_format::UINT8:
        return 1.0 / 127.0;
    case Iq_sample_format::INT16:
    case Iq_sample_format::UINT16:
        return 1.0 / 32767.0;
    case Iq_sample_format::FLOAT32:
        return 1e-6;
    }
    return 1e-6;
}
} // namespace

TEST_CASE( "iq_format_sizes_and_names", "[io][format]" )
{
    REQUIRE( bytes_per_iq_sample( Iq_sample_format::INT8 ) == 2 );
    REQUIRE( bytes_per_iq_sample( Iq_sample_format::UINT8 ) == 2 );
    REQUIRE( bytes_per_iq_sample( Iq_sample_format::INT16 ) == 4 );
    REQUIRE( bytes_per_iq_sample( Iq_sample_format::UINT16 ) == 4 );
    REQUIRE( bytes_per_iq_sample( Iq_sample_format::FLOAT32 ) == 8 );

    REQUIRE( std::string( iq_format_name( Iq_sample_format::INT8 ) ) == "int8" );
    REQUIRE( std::string( iq_format_name( Iq_sample_format::UINT8 ) ) == "uint8" );
    REQUIRE( std::string( iq_format_name( Iq_sample_format::INT16 ) ) == "int16" );
    REQUIRE( std::string( iq_format_name( Iq_sample_format::UINT16 ) ) == "uint16" );
    REQUIRE( std::string( iq_format_name( Iq_sample_format::FLOAT32 ) ) == "float32" );
}

TEST_CASE( "iq_format_encode_decode_round_trip", "[io][format]" )
{
    const Complex_buf in = { { 0.5f, -0.25f }, { 1.0f, -1.0f }, { 0.0f, 0.123f }, { -0.75f, 0.9f } };
    for( Iq_sample_format fmt : ALL_FORMATS )
    {
        const std::vector<uint8_t> bytes = encode_iq( in, fmt );
        REQUIRE( bytes.size() == in.size() * static_cast<size_t>( bytes_per_iq_sample( fmt ) ) );

        Complex_buf out( in.size() );
        decode_iq( bytes.data(), bytes.size(), fmt, out );

        const double tol = round_trip_tol( fmt );
        INFO( "format " << iq_format_name( fmt ) );
        for( size_t i = 0; i < in.size(); ++i )
        {
            REQUIRE( out[i].real() == Catch::Approx( in[i].real() ).margin( tol ) );
            REQUIRE( out[i].imag() == Catch::Approx( in[i].imag() ).margin( tol ) );
        }
    }
}

TEST_CASE( "iq_format_record_replay_is_byte_stable", "[io][format]" )
{
    // A recording is encode(pipeline floats). Replaying it (decode) then re-recording (encode) must reproduce
    // the exact bytes - i.e. decode is the exact inverse of encode's output, so a recorded file round-trips
    // losslessly through replay. (Floats outside the encode lattice would not, but a recording is always on it.)
    const Complex_buf in = { { 0.5f, -0.25f }, { 1.0f, -1.0f }, { -0.3f, 0.42f }, { 0.0f, -0.8f } };
    for( Iq_sample_format fmt : ALL_FORMATS )
    {
        const std::vector<uint8_t> b1 = encode_iq( in, fmt );
        Complex_buf                mid( in.size() );
        decode_iq( b1.data(), b1.size(), fmt, mid );
        const std::vector<uint8_t> b2 = encode_iq( mid, fmt );
        INFO( "format " << iq_format_name( fmt ) );
        REQUIRE( b1 == b2 );
    }
}

TEST_CASE( "iq_format_exact_scaling_anchors", "[io][format]" )
{
    // int8: +1.0 -> 127, 0 -> 0, -1.0 -> -127 (full-scale +/-127, mid 0).
    const std::vector<uint8_t> i8 = encode_iq( { { 1.0f, 0.0f }, { -1.0f, 0.0f } }, Iq_sample_format::INT8 );
    REQUIRE( static_cast<int8_t>( i8[0] ) == 127 );
    REQUIRE( static_cast<int8_t>( i8[1] ) == 0 );
    REQUIRE( static_cast<int8_t>( i8[2] ) == -127 );
    // uint8: +1.0 -> 255, -1.0 -> 0 (mid-rail 127.5; clamped to the byte range).
    const std::vector<uint8_t> u8 = encode_iq( { { 1.0f, -1.0f } }, Iq_sample_format::UINT8 );
    REQUIRE( u8[0] == 255 );
    REQUIRE( u8[1] == 0 );
    // Out-of-range inputs clamp to the type range [-128, 127], never wrap. (Normal [-1,1] input only ever
    // reaches [-127, 127]; -128 is reachable solely by clamping an overdriven sample.)
    const std::vector<uint8_t> clamp = encode_iq( { { 2.0f, -2.0f } }, Iq_sample_format::INT8 );
    REQUIRE( static_cast<int8_t>( clamp[0] ) == 127 );
    REQUIRE( static_cast<int8_t>( clamp[1] ) == -128 );
}
#endif

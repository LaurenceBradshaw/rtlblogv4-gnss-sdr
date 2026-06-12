#include "iq_file_device.h"
#include <algorithm>
#include "logging.h"

namespace
{
template <typename T>
T read_trivial( const unsigned char* data )
{
    T value {};
    std::memcpy( &value, data, sizeof( T ) );
    return value;
}

// Branchless normalisers - multiply instead of divide, clamp instead of
// if/return chains. This allows the compiler to auto-vectorize the conversion
// loop in convert_buffer (the branchy originals blocked vectorisation).
static constexpr float INV_127     = 1.0f / 127.0f;
static constexpr float INV_127_5   = 1.0f / 127.5f;
static constexpr float INV_32767   = 1.0f / 32767.0f;
static constexpr float INV_32767_5 = 1.0f / 32767.5f;

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

int bytes_per_complex_sample( Iq_sample_format sample_format )
{
    switch( sample_format )
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

    throw std::runtime_error( "Unsupported IQ file format" );
}

void convert_buffer( const unsigned char* buf, size_t len, Iq_sample_format sample_format, Complex_buf& output )
{
    const size_t bytes_per_sample = bytes_per_complex_sample( sample_format );
    const size_t n_samples        = len / bytes_per_sample;

    for( size_t i = 0; i < n_samples; ++i )
    {
        const unsigned char* sample_ptr = buf + i * bytes_per_sample;

        switch( sample_format )
        {
        case Iq_sample_format::UINT8:
        {
            output[i] = Complex_sample( normalise_uint8( sample_ptr[0] ), normalise_uint8( sample_ptr[1] ) );
            break;
        }
        case Iq_sample_format::INT8:
        {
            output[i] = Complex_sample(
                normalise_int8( static_cast<int8_t>( sample_ptr[0] ) ), normalise_int8( static_cast<int8_t>( sample_ptr[1] ) )
            );
            break;
        }
        case Iq_sample_format::UINT16:
        {
            output[i] = Complex_sample(
                normalise_uint16( read_trivial<uint16_t>( sample_ptr ) ),
                normalise_uint16( read_trivial<uint16_t>( sample_ptr + sizeof( uint16_t ) ) )
            );
            break;
        }
        case Iq_sample_format::INT16:
        {
            output[i] = Complex_sample(
                normalise_int16( read_trivial<int16_t>( sample_ptr ) ),
                normalise_int16( read_trivial<int16_t>( sample_ptr + sizeof( int16_t ) ) )
            );
            break;
        }
        case Iq_sample_format::FLOAT32:
        {
            const auto i_raw = read_trivial<float>( sample_ptr );
            const auto q_raw = read_trivial<float>( sample_ptr + sizeof( float ) );
            output[i]        = Complex_sample( i_raw, q_raw );
            break;
        }
        }
    }
}
} // namespace

Iq_file_device::Iq_file_device( const std::string& path, uint32_t sample_rate_hz, Iq_sample_format sample_format )
    : sample_rate_hz_( sample_rate_hz ),
      sample_format_( sample_format )
{
    if( path.empty() )
    {
        logging::log( logging::Level::Error, "IQ file path is empty" );
        throw std::runtime_error( "IQ file path is empty" );
    }
    file_.open( path, std::ios::binary );
    if( !file_ )
    {
        logging::log( logging::Level::Error, "Failed to open IQ file" );
        throw std::runtime_error( "Failed to open IQ file" );
    }
}

Iq_file_device::~Iq_file_device()
{
    stop_streaming();

    if( file_.is_open() )
    {
        file_.close();
    }
}

void Iq_file_device::stop_streaming()
{
    running_ = false; // signal the thread to stop (idempotent)
    if( worker_.joinable() )
    {
        worker_.join(); // always join if joinable, regardless of who set running_=false
    }
}

void Iq_file_device::streaming_thread()
{
    constexpr size_t BUFFER_SAMPLES   = 8192;
    const int        bytes_per_sample = bytes_per_complex_sample( sample_format_ );

    std::vector<unsigned char> raw_buf( BUFFER_SAMPLES * bytes_per_sample );

    Complex_buf buffer;
    buffer.resize( BUFFER_SAMPLES );

    while( running_ )
    {
        file_.read( reinterpret_cast<char*>( raw_buf.data() ), static_cast<std::streamsize>( raw_buf.size() ) );

        const size_t bytes_read   = static_cast<size_t>( file_.gcount() );
        const size_t samples_read = bytes_read / bytes_per_sample;
        if( samples_read == 0 )
            break;

        buffer.resize( samples_read );
        convert_buffer( raw_buf.data(), bytes_read, sample_format_, buffer );

        // push() blocks internally when the ring buffer is full, providing
        // natural backpressure - no real-time sleep needed.
        callback_( buffer );

        samples_consumed_.fetch_add( samples_read, std::memory_order_relaxed );
    }

    running_ = false;
}

void Iq_file_device::start_streaming( Sample_callback callback )
{
    if( running_ )
    {
        return;
    }

    callback_ = std::move( callback );
    running_  = true;

    worker_ = std::thread( &Iq_file_device::streaming_thread, this );
}

bool Iq_file_device::is_streaming() const
{
    return running_;
}

uint32_t Iq_file_device::sample_rate_hz() const
{
    return sample_rate_hz_;
}

uint64_t Iq_file_device::samples_consumed() const
{
    return samples_consumed_.load( std::memory_order_relaxed );
}
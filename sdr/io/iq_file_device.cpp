#include "iq_file_device.h"
#include <vector>
#include "logging.h"

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
    const int        bytes_per_sample = bytes_per_iq_sample( sample_format_ );

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
        decode_iq( raw_buf.data(), bytes_read, sample_format_, buffer );

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
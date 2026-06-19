#include "rtlsdr_device.h"
#include <rtl-sdr.h>
#include <stdexcept>
#include "logging.h"

// Open the device. Throws std::runtime_error if none is found or the open fails.
Rtlsdr_device::Rtlsdr_device( int device_index )
    : dev_( nullptr )
{
    if( rtlsdr_get_device_count() == 0 )
    {
        throw std::runtime_error( "No RTL-SDR devices found" );
    }
    if( rtlsdr_open( reinterpret_cast<rtlsdr_dev_t**>( &dev_ ), static_cast<uint32_t>( device_index ) ) != 0 )
    {
        throw std::runtime_error( "Failed to open RTL-SDR device" );
    }
}

Rtlsdr_device::~Rtlsdr_device()
{
    stop_streaming();
    if( dev_ != nullptr )
    {
        rtlsdr_close( reinterpret_cast<rtlsdr_dev_t*>( dev_ ) );
    }
}

void Rtlsdr_device::set_centre_freq_hz( uint32_t freq_hz )
{
    if( rtlsdr_set_center_freq( reinterpret_cast<rtlsdr_dev_t*>( dev_ ), freq_hz ) != 0 )
    {
        throw std::runtime_error( "Failed to set RTL-SDR centre frequency" );
    }
}

void Rtlsdr_device::set_sample_rate_hz( uint32_t rate_hz )
{
    if( rtlsdr_set_sample_rate( reinterpret_cast<rtlsdr_dev_t*>( dev_ ), rate_hz ) != 0 )
    {
        throw std::runtime_error( "Failed to set RTL-SDR sample rate" );
    }
}

void Rtlsdr_device::set_gain_tenths_db( int gain_tenths_db )
{
    auto* dev = reinterpret_cast<rtlsdr_dev_t*>( dev_ );
    if( rtlsdr_set_tuner_gain_mode( dev, 1 ) != 0 ) // manual gain (disable tuner AGC)
    {
        throw std::runtime_error( "Failed to set RTL-SDR manual gain mode" );
    }
    if( rtlsdr_set_tuner_gain( dev, gain_tenths_db ) != 0 )
    {
        throw std::runtime_error( "Failed to set RTL-SDR tuner gain" );
    }
}

void Rtlsdr_device::set_agc( bool enable )
{
    auto* dev = reinterpret_cast<rtlsdr_dev_t*>( dev_ );
    rtlsdr_set_tuner_gain_mode( dev, enable ? 0 : 1 );
    rtlsdr_set_agc_mode( dev, enable ? 1 : 0 );
}

bool Rtlsdr_device::is_streaming() const
{
    return streaming_.load( std::memory_order_relaxed );
}

uint32_t Rtlsdr_device::sample_rate_hz() const
{
    return rtlsdr_get_sample_rate( reinterpret_cast<rtlsdr_dev_t*>( dev_ ) );
}

uint32_t Rtlsdr_device::centre_freq_hz() const
{
    return rtlsdr_get_center_freq( reinterpret_cast<rtlsdr_dev_t*>( dev_ ) );
}

int Rtlsdr_device::tuner_gain_tenths_db() const
{
    return rtlsdr_get_tuner_gain( reinterpret_cast<rtlsdr_dev_t*>( dev_ ) );
}

uint64_t Rtlsdr_device::samples_consumed() const
{
    return samples_consumed_.load( std::memory_order_relaxed );
}

void Rtlsdr_device::start_streaming( Sample_callback callback )
{
    if( streaming_ )
    {
        return;
    }
    callback_  = std::move( callback );
    streaming_ = true;

    if( rtlsdr_reset_buffer( reinterpret_cast<rtlsdr_dev_t*>( dev_ ) ) != 0 )
    {
        streaming_ = false;
        throw std::runtime_error( "Failed to reset RTL-SDR buffer" );
    }

    // rtlsdr_read_async() blocks until cancelled, so run it on a worker thread - start_streaming()
    // then returns immediately, matching Iq_file_device.
    worker_ = std::thread( &Rtlsdr_device::streaming_thread, this );
}

void Rtlsdr_device::streaming_thread()
{
    // buf_count = 0, buf_len = 0 -> librtlsdr defaults (15 buffers of 16384 bytes = 8192 complex
    // samples each, ~2 ms per buffer at 4 MSPS). Blocks here until stop_streaming() cancels it.
    rtlsdr_read_async( reinterpret_cast<rtlsdr_dev_t*>( dev_ ), &Rtlsdr_device::raw_callback, this, 0, 0 );
    streaming_ = false;
}

void Rtlsdr_device::stop_streaming()
{
    if( dev_ != nullptr )
    {
        rtlsdr_cancel_async( reinterpret_cast<rtlsdr_dev_t*>( dev_ ) );
    }
    if( worker_.joinable() )
    {
        worker_.join();
    }
}

void Rtlsdr_device::raw_callback( unsigned char* buf, uint32_t len, void* ctx )
{
    auto* self = static_cast<Rtlsdr_device*>( ctx );

    // RTL-SDR delivers interleaved uint8 I/Q in offset binary (0 -> -127.5, 255 -> +127.5);
    // dividing by 127.5 normalises to [-1, 1] (matching Iq_file_device's UINT8 path).
    const size_t n_samples = len / 2;
    Complex_buf  samples( n_samples );
    for( size_t i = 0; i < n_samples; ++i )
    {
        const float i_sample = ( static_cast<float>( buf[2 * i] ) - 127.5f ) / 127.5f;
        const float q_sample = ( static_cast<float>( buf[2 * i + 1] ) - 127.5f ) / 127.5f;
        samples[i]           = Complex_sample( i_sample, q_sample );
    }

    self->callback_( samples );
    self->samples_consumed_.fetch_add( n_samples, std::memory_order_relaxed );
}

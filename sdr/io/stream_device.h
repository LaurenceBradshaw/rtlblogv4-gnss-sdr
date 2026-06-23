#pragma once

#include <functional>
#include <vector>
#include "types.h"

// Abstract IQ sample source: a recorded file (Iq_file_device) or a live front-end (Rtlsdr_device).
// The rest of the pipeline only sees this interface, so it is source-agnostic. start_streaming()
// spawns the producer; samples are delivered through the callback (push to the Sample_buffer).
class Stream_device
{
public:
    using Sample_callback = std::function<void( const Complex_buf& samples )>;

    virtual ~Stream_device() = default;

    // --- streaming (file + hardware) ------------------------------------------------------
    virtual void     start_streaming( Sample_callback callback ) = 0;
    virtual void     stop_streaming()                            = 0;
    virtual bool     is_streaming() const                        = 0; // false at file EOF; always true on a live source
    virtual uint32_t sample_rate_hz() const                      = 0; // fixed for a file, configured on hardware
    virtual uint64_t samples_consumed() const                    = 0;

    // --- live RF front-end tuning (hardware only) -----------------------------------------
    // No-ops by default: a recorded file's centre frequency / sample rate / gain / AGC are fixed at
    // capture time, so a file source ignores these. Live front-ends (Rtlsdr_device) override them.
    // Keeping them on the base lets callers configure any source uniformly without down-casting.
    virtual void set_centre_freq_hz( uint32_t /*freq_hz*/ ) {}
    virtual void set_sample_rate_hz( uint32_t /*rate_hz*/ ) {}
    virtual void set_gain_tenths_db( int /*gain_tenths_db*/ ) {}
    virtual void set_agc( bool /*enable*/ ) {}
    // Bias-tee: feed DC up the coax to power an ACTIVE antenna's LNA (most GPS antennas need this). Off by
    // default - it puts DC on the antenna port, so enable it deliberately. (Equivalent to `rtl_biast -b 1`.)
    virtual void set_bias_tee( bool /*enable*/ ) {}
};

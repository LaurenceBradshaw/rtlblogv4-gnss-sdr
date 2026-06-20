#pragma once
#include <cstddef>
#include <vector>
#include "types.h"

// Streaming FIR low-pass + integer decimation. Anti-alias filters then keeps every `factor`-th sample,
// so a stream at fs comes out at fs/factor without aliasing the surviving band. Stateful: feed it
// arbitrary-sized blocks via process() and it carries the filter history across calls, so the
// decimation phase stays continuous. The filter is a Hamming-windowed sinc (linear phase, unity DC
// gain); the constant group delay is harmless for GNSS (it folds into code phase / the clock bias).
class Fir_decimator
{
public:
    // factor : decimation ratio M (>=1; 1 is a pass-through still applying the low-pass).
    // num_taps : FIR length; <=0 picks a sensible default (8*factor+1, forced odd).
    // cutoff_norm : low-pass cutoff in cycles/sample; <=0 picks 0.45/factor (just under the fs/2M
    //               output Nyquist, leaving a small transition band).
    explicit Fir_decimator( int factor, int num_taps = 0, double cutoff_norm = 0.0 );

    // Filter + decimate a block; returns the (approximately in.size()/factor) output samples. The
    // first (num_taps-1) input samples warm up the filter before any output is produced.
    Complex_buf process( const Complex_buf& in );

    int factor() const
    {
        return factor_;
    }
    const std::vector<float>& taps() const
    {
        return taps_;
    }

private:
    int                factor_;
    std::vector<float> taps_;
    // Filter state carried across process() calls, kept DEINTERLEAVED (I and Q in separate contiguous
    // float arrays) so each output is two real dot-products - which the SSE dot helper vectorises cleanly
    // even in an unoptimised (debug) build.
    std::vector<float> hist_i_;
    std::vector<float> hist_q_;
    size_t             pos_ = 0; // next output window start within hist_i_/hist_q_
};

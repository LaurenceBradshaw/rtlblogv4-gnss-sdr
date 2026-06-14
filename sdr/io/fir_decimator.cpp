#include "fir_decimator.h"
#include <algorithm>
#include <cmath>
#if defined( __SSE2__ )
#include <emmintrin.h>
#endif

namespace
{
constexpr double PI = 3.14159265358979323846;

// Real dot product sum_{k<n} a[k]*b[k]. SSE2 (4-wide) where available - SSE2 is baseline on x86-64, so
// this vectorises even at -O0 (a debug build); a portable scalar path covers other targets.
float dot( const float* a, const float* b, size_t n )
{
#if defined( __SSE2__ )
    __m128 acc = _mm_setzero_ps();
    size_t i   = 0;
    for( ; i + 4 <= n; i += 4 )
    {
        acc = _mm_add_ps( acc, _mm_mul_ps( _mm_loadu_ps( a + i ), _mm_loadu_ps( b + i ) ) );
    }
    float tmp[4];
    _mm_storeu_ps( tmp, acc );
    float s = tmp[0] + tmp[1] + tmp[2] + tmp[3];
    for( ; i < n; ++i ) // tail (n not a multiple of 4)
    {
        s += a[i] * b[i];
    }
    return s;
#else
    float s = 0.0f;
    for( size_t i = 0; i < n; ++i )
    {
        s += a[i] * b[i];
    }
    return s;
#endif
}

double sinc( double t ) // normalised sinc: sin(pi t)/(pi t), sinc(0)=1
{
    if( std::abs( t ) < 1e-12 )
    {
        return 1.0;
    }
    return std::sin( PI * t ) / ( PI * t );
}

// Hamming-windowed-sinc low-pass, cutoff wc (cycles/sample), normalised to unity DC gain.
std::vector<float> design_lowpass( int n, double wc )
{
    std::vector<float> h( n );
    const double       c   = ( n - 1 ) / 2.0;
    double             sum = 0.0;
    for( int i = 0; i < n; ++i )
    {
        const double x      = i - c;
        const double ideal  = 2.0 * wc * sinc( 2.0 * wc * x );                    // ideal LP impulse response
        const double window = 0.54 - 0.46 * std::cos( 2.0 * PI * i / ( n - 1 ) ); // Hamming
        const double v      = ideal * window;
        h[i]                = static_cast<float>( v );
        sum += v;
    }
    if( sum != 0.0 )
    {
        for( float& v : h )
        {
            v = static_cast<float>( v / sum ); // unity DC gain
        }
    }
    return h;
}
} // namespace

Fir_decimator::Fir_decimator( int factor, int num_taps, double cutoff_norm )
    : factor_( std::max( 1, factor ) )
{
    int n = num_taps > 0 ? num_taps : 8 * factor_ + 1;
    if( n % 2 == 0 )
    {
        ++n; // odd -> symmetric, integer group delay
    }
    const double wc = cutoff_norm > 0.0 ? cutoff_norm : 0.45 / factor_;
    taps_           = design_lowpass( n, wc );
}

Complex_buf Fir_decimator::process( const Complex_buf& in )
{
    hist_i_.reserve( hist_i_.size() + in.size() );
    hist_q_.reserve( hist_q_.size() + in.size() );
    for( const Complex_sample& s : in ) // deinterleave I/Q for the SIMD dot-products
    {
        hist_i_.push_back( s.real() );
        hist_q_.push_back( s.imag() );
    }

    const size_t n = taps_.size();
    const float* t = taps_.data();
    Complex_buf  out;
    out.reserve( hist_i_.size() / static_cast<size_t>( factor_ ) + 1 );

    while( pos_ + n <= hist_i_.size() )
    {
        out.emplace_back( dot( t, hist_i_.data() + pos_, n ), dot( t, hist_q_.data() + pos_, n ) );
        pos_ += static_cast<size_t>( factor_ );
    }

    // Drop the samples no future output needs (everything before the next window start), keeping the
    // tail (< num_taps samples) as the filter state for the next block. Phase carries via pos_ -> 0.
    if( pos_ > 0 )
    {
        const auto cut = static_cast<std::ptrdiff_t>( pos_ );
        hist_i_.erase( hist_i_.begin(), hist_i_.begin() + cut );
        hist_q_.erase( hist_q_.begin(), hist_q_.begin() + cut );
        pos_ = 0;
    }
    return out;
}

#ifdef ENABLE_UNIT_TESTS
#include <complex>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace
{
// Steady-state magnitude response to a complex tone exp(j*2*pi*f*n) through the decimator.
double tone_response( int factor, double f_cycles_per_sample )
{
    Fir_decimator dec( factor );
    Complex_buf   in;
    const int     count = 20000;
    for( int n = 0; n < count; ++n )
    {
        const double ph = 2.0 * PI * f_cycles_per_sample * n;
        in.emplace_back( static_cast<float>( std::cos( ph ) ), static_cast<float>( std::sin( ph ) ) );
    }
    const Complex_buf out = dec.process( in );
    double            sum = 0.0;
    int               k   = 0;
    for( size_t i = out.size() / 2; i < out.size(); ++i, ++k ) // average over the settled back half
    {
        sum += std::abs( std::complex<float>( out[i] ) );
    }
    return k > 0 ? sum / k : 0.0;
}
} // namespace

TEST_CASE( "fir_decimator_taps_unity_dc_gain", "[io][decimator]" )
{
    Fir_decimator dec( 4 );
    double        sum = 0.0;
    for( float t : dec.taps() )
    {
        sum += t;
    }
    REQUIRE( sum == Catch::Approx( 1.0 ).margin( 1e-6 ) ); // DC passes at unity
    REQUIRE( dec.taps().size() % 2 == 1 );                 // symmetric / linear phase
}

TEST_CASE( "fir_decimator_output_rate_and_dc", "[io][decimator]" )
{
    const int     M = 5;
    Fir_decimator dec( M );
    Complex_buf   in( 50000, Complex_sample( 1.0f, 0.0f ) ); // DC

    const Complex_buf out = dec.process( in );
    REQUIRE( out.size() == Catch::Approx( in.size() / M ).epsilon( 0.01 ) ); // ~one output per M inputs
    REQUIRE( out.back().real() == Catch::Approx( 1.0 ).margin( 1e-3 ) );     // unity DC gain after warm-up
    REQUIRE( out.back().imag() == Catch::Approx( 0.0 ).margin( 1e-3 ) );
}

TEST_CASE( "fir_decimator_passband_and_stopband", "[io][decimator]" )
{
    const int M = 8; // output Nyquist = 0.5/8 = 0.0625 cycles/sample
    REQUIRE( tone_response( M, 0.01 ) == Catch::Approx( 1.0 ).margin( 0.05 ) ); // passband tone survives
    REQUIRE( tone_response( M, 0.20 ) < 0.05 );                                 // aliasing tone rejected
}

TEST_CASE( "fir_decimator_phase_continuous_across_blocks", "[io][decimator]" )
{
    // One big block vs many odd-sized blocks must yield the same decimated stream.
    Complex_buf in;
    for( int n = 0; n < 12000; ++n )
    {
        const double ph = 2.0 * PI * 0.013 * n;
        in.emplace_back( static_cast<float>( std::cos( ph ) ), static_cast<float>( std::sin( ph ) ) );
    }

    Fir_decimator     one( 6 );
    const Complex_buf whole = one.process( in );

    Fir_decimator chunked( 6 );
    Complex_buf   pieced;
    for( size_t i = 0; i < in.size(); i += 257 ) // odd block size, not a multiple of the factor
    {
        const Complex_buf blk( in.begin() + static_cast<std::ptrdiff_t>( i ),
                               in.begin() + static_cast<std::ptrdiff_t>( std::min( in.size(), i + 257 ) ) );
        const Complex_buf o = chunked.process( blk );
        pieced.insert( pieced.end(), o.begin(), o.end() );
    }

    REQUIRE( pieced.size() == whole.size() );
    for( size_t i = 0; i < whole.size(); ++i )
    {
        REQUIRE( pieced[i].real() == Catch::Approx( whole[i].real() ).margin( 1e-5 ) );
        REQUIRE( pieced[i].imag() == Catch::Approx( whole[i].imag() ).margin( 1e-5 ) );
    }
}
#endif

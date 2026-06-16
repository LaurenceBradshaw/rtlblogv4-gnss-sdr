#include "acquisition.h"
#include <algorithm>
#include <cmath>
#include <limits>

// constructor
// mirrors GNSS-SDRLIB initsdrch() + initacqstruct():
//
//   acq->intg   = ACQINTG_L1CA;
//   acq->step   = ACQSTEP;
//   acq->nfreq  = 2*(ACQHBAND/ACQSTEP)+1;
//   acq->nfft   = 2*nsamp;
//
//   for (i=0; i<nfft; i++) rcode[i] = 0;      // zero-pad
//   rescode(..., nsamp, rcode);                 // fills first nsamp chips
//   cpxcpx(rcode, NULL, 1.0, nfft, xcode);     // real -> complex, Q=0
//   cpxfft(NULL, xcode, nfft);                  // FFT in-place
Acquisition_engine::Acquisition_engine(
    const Complex_buf& prn_code, Satellite_id satellite_id, double sample_rate_hz, const Signal_params& sig,
    const Acquisition_aiding& aiding
)
    : aiding_( aiding ),
      satellite_id_( satellite_id ),
      constellation_( sig.constellation ),
      doppler_center_bins_( 0 ),
      active_half_bins_( 0 ),
      result_ {},
      intg_count_( 0 )
{
    n_                   = static_cast<int>( prn_code.size() );
    m_                   = sig.acq_fft_factor * n_; // nfft (factor 2 = zero-padded, 1 = circular)
    target_integrations_ = sig.acq_integrations;
    chip_rate_           = sig.chip_rate_hz;
    carrier_freq_hz_     = sig.carrier_freq_hz;
    ti_                  = 1.0 / sample_rate_hz;
    ctime_               = sig.code_period_s; // code period - C/N0 normalisation divisor
    nsamp_chip_          = static_cast<int>( std::round( sample_rate_hz / chip_rate_ ) );

    // Doppler search grid. The FFT-shift acquisition (see integrate()) wipes the
    // carrier by circularly shifting the data spectrum, so Doppler bins must align
    // to the FFT bin spacing fs/m_ (~500 Hz for a 2x-padded 1 ms block). The
    // parabolic peak interpolation in check_acquisition() recovers sub-bin accuracy.
    // The grid is allocated for the WIDE half-width; each attempt searches the active
    // (possibly narrower) sub-range the aiding recommends - see integrate().
    doppler_step_hz_ = sample_rate_hz / static_cast<double>( m_ );
    nfreq_           = 2 * static_cast<int>( Acquisition_aiding::SEARCH_HBAND_WIDE_HZ / doppler_step_hz_ ) + 1;
    doppler_freqs_.resize( nfreq_ );
    for( int i = 0; i < nfreq_; i++ )
    {
        doppler_freqs_[i] = ( i - ( nfreq_ - 1 ) / 2 ) * doppler_step_hz_; // bin index -> Hz
    }

    // Build zero-padded code FFT - mirrors initsdrch():
    //   zero-fill nfft elements, resample code into first nsamp, then FFT.
    // Stored pre-conjugated so integrate() multiplies by conj(FFT(code)) directly.
    code_fft_.assign( m_, Complex_sample( 0.0f, 0.0f ) );
    for( int i = 0; i < n_; i++ )
    {
        code_fft_[i] = prn_code[i]; // Q=0, from L1ca_code::generate
    }

    auto* fc   = reinterpret_cast<fftwf_complex*>( code_fft_.data() );
    fft_plan_  = fftwf_plan_dft_1d( m_, fc, fc, FFTW_FORWARD, FFTW_ESTIMATE );
    ifft_plan_ = fftwf_plan_dft_1d( m_, fc, fc, FFTW_BACKWARD, FFTW_ESTIMATE );
    fftwf_execute( fft_plan_ ); // FFT code once: xcode = fft(code)
    for( auto& c : code_fft_ )
    {
        c = std::conj( c ); // pre-conjugate: integrate() needs conj(FFT(code))
    }

    // Pre-allocate scratch buffers (avoid heap allocation per integrate call)
    data_fft_.resize( m_ );
    datax_.resize( m_ );

    power_.assign( static_cast<size_t>( n_ ) * nfreq_, 0.0 );
}

Acquisition_engine::~Acquisition_engine()
{
    fftwf_destroy_plan( fft_plan_ );
    fftwf_destroy_plan( ifft_plan_ );
}

// reset
void Acquisition_engine::reset()
{
    std::fill( power_.begin(), power_.end(), 0.0 );
    intg_count_ = 0;
    result_     = {};
}

// integrate
// One non-coherent integration epoch - FFT-shift PCPS (parallel code-phase search).
//
// Equivalent to GNSS-SDRLIB's pcorrelator(): for each Doppler the correlation is
//   P = | IFFT( FFT(carrier_wiped_data) * conj(FFT(code)) ) |^2
// but instead of carrier-wiping in time and forward-FFT'ing per Doppler bin, we use
// the shift theorem: wiping by f = k*(fs/m) is a circular shift of the spectrum by
// k bins. So we FFT the raw block ONCE, then per Doppler bin rotate that spectrum by
// k, multiply by conj(FFT(code)), IFFT, and accumulate power. This replaces nfreq
// forward FFTs + nfreq carrier-wipe loops with a single forward FFT.
//
// Returns true if check_acquisition() passes (peak ratio > ACQTH).
bool Acquisition_engine::integrate( const Sample_block& block )
{
    const float m2     = static_cast<float>( m_ ) * static_cast<float>( m_ );
    const int   centre = ( nfreq_ - 1 ) / 2;

    // Latch the recenter AND the search half-width for this whole attempt (the power map
    // accumulates across epochs at one grid, so neither may move mid-attempt). The aiding
    // gives a common-mode offset (receiver clock) to search around instead of 0 Hz, plus how
    // wide to search: wide while bootstrapping, narrow once >=2 SVs pin the offset.
    if( intg_count_ == 0 )
    {
        const Acquisition_aiding::Estimate est = aiding_.estimate( constellation_, satellite_id_, carrier_freq_hz_ );
        doppler_center_bins_                   = static_cast<int>( std::lround( est.center_hz / doppler_step_hz_ ) );
        const int half                         = static_cast<int>( std::lround( est.half_width_hz / doppler_step_hz_ ) );
        active_half_bins_                      = std::min( half, ( nfreq_ - 1 ) / 2 ); // <= allocated grid
    }

    // FFT the raw block ONCE (carrier-free). Guard against a block shorter than m_ (zero-pad the tail)
    // so a sizing mismatch can never read out of bounds - the FFT tolerates a few zero samples.
    const size_t avail = block.len1 + block.len2;
    for( int i = 0; i < m_; i++ )
    {
        const size_t ui = static_cast<size_t>( i );
        data_fft_[i]    = ui >= avail            ? Complex_sample( 0.0f, 0.0f )
                          : ui < block.len1      ? block.ptr1[i]
                                                 : block.ptr2[ui - block.len1];
    }
    auto* fd = reinterpret_cast<fftwf_complex*>( data_fft_.data() );
    fftwf_execute_dft( fft_plan_, fd, fd );

    auto* fx = reinterpret_cast<fftwf_complex*>( datax_.data() );

    // Access the complex buffers as raw float[2] (no std::complex operator* / real()
    // / imag() calls - those are out-of-line at -O0; this is the acquisition hot loop).
    const float* dff = reinterpret_cast<const float*>( data_fft_.data() );
    const float* cff = reinterpret_cast<const float*>( code_fft_.data() ); // pre-conjugated
    float*       dx  = reinterpret_cast<float*>( datax_.data() );

    // Only the active sub-range [centre +/- active_half_bins_] is searched; the rest of the
    // allocated grid stays zeroed (from reset()), so check_acquisition's scan ignores it.
    const int fi_lo = centre - active_half_bins_;
    const int fi_hi = centre + active_half_bins_;
    for( int fi = fi_lo; fi <= fi_hi; fi++ )
    {
        // Carrier wipe by f = k*(fs/m) == circular shift of the spectrum by k bins:
        //   X'[j] = data_fft[(j - k) mod m].  Fuse the shift with the conj(code)
        //   multiply. doppler_center_bins_ offsets the whole search onto the aiding
        //   estimate; |k| stays small vs m_, so the wrap still touches only the edges.
        const int k = fi - centre + doppler_center_bins_;

        for( int j = 0; j < m_; j++ )
        {
            int src = j - k;
            if( src < 0 )
                src += m_;
            else if( src >= m_ )
                src -= m_;
            const float dr  = dff[2 * src];
            const float di  = dff[2 * src + 1];
            const float ccr = cff[2 * j];
            const float cci = cff[2 * j + 1];
            dx[2 * j]       = dr * ccr - di * cci;
            dx[2 * j + 1]   = dr * cci + di * ccr;
        }

        // IFFT - cpxifft(iplan, datax, m) in cpxconv
        fftwf_execute_dft( ifft_plan_, fx, fx );

        // Accumulate power (first n_ samples) - cpxconv flagsum=1
        double* P = power_.data() + fi * n_;
        for( int i = 0; i < n_; i++ )
        {
            const float re = dx[2 * i];
            const float im = dx[2 * i + 1];
            P[i] += static_cast<double>( ( re * re + im * im ) / m2 );
        }
    }

    ++intg_count_;
    return check_acquisition();
}

// check_acquisition
// Exact mirror of GNSS-SDRLIB checkacquisition():
//
//   maxP   = maxvd(P, nsamp*nfreq, -1, -1, &maxi)
//   ind2sub(maxi, nsamp, nfreq, &codei, &freqi)
//   exinds = codei - 2*nsampchip;  if < 0:  += nsamp
//   exinde = codei + 2*nsampchip;  if >=n_: -= nsamp
//   meanP  = meanvd(&P[freqi*nsamp], nsamp, exinds, exinde)
//   cn0    = 10*log10(maxP / meanP / ctime)
//   maxP2  = maxvd(&P[freqi*nsamp], nsamp, exinds, exinde, &maxi2)
//   peakr  = maxP / maxP2
//   return peakr > ACQTH
bool Acquisition_engine::check_acquisition()
{
    const int total = n_ * nfreq_;

    // Global maximum (no exclusion: exinds=exinde=-1 -> all included)
    // mirrors: maxvd(P, nsamp*nfreq, -1, -1, &maxi)
    int    maxi = 0;
    double maxP = power_[0];
    for( int i = 1; i < total; i++ )
    {
        if( power_[i] > maxP )
        {
            maxP = power_[i];
            maxi = i;
        }
    }

    // ind2sub(maxi, nsamp, nfreq, &codei, &freqi)
    //   *subx = ind % nx          -> codei = maxi % n_
    //   *suby = ny*ind/(nx*ny)    -> freqi = maxi / n_
    const int codei = maxi % n_;
    const int freqi = maxi / n_;

    // Exclusion zone: +/-2 chips around the peak, with wrap-around
    int exinds = codei - 2 * nsamp_chip_;
    if( exinds < 0 )
    {
        exinds += n_;
    }
    int exinde = codei + 2 * nsamp_chip_;
    if( exinde >= n_ )
    {
        exinde -= n_;
    }

    // Exclusion predicate (mirrors meanvd / maxvd include conditions)
    // non-wrapped (exinds<=exinde): exclude [exinds, exinde]
    // wrapped     (exinds> exinde): exclude [0,exinde] union [exinds,n_-1]
    const auto is_excluded = [&]( int i ) -> bool
    {
        if( exinds <= exinde )
        {
            return i >= exinds && i <= exinde;
        }
        else
        {
            return i >= exinds || i <= exinde;
        }
    };

    // Mean power in the peak Doppler bin, excluding the peak zone
    // mirrors: meanvd(&P[freqi*nsamp], nsamp, exinds, exinde)
    const double* Pfreq = power_.data() + freqi * n_;
    double        sumP  = 0.0;
    int           ne    = 0;
    for( int i = 0; i < n_; i++ )
    {
        if( !is_excluded( i ) )
        {
            sumP += Pfreq[i];
        }
        else
        {
            ++ne;
        }
    }
    const double meanP = sumP / ( n_ - ne );

    // Second peak in the same Doppler bin, excluding the peak zone
    // mirrors: maxvd(&P[freqi*nsamp], nsamp, exinds, exinde, &maxi2)
    // Note: GNSS-SDRLIB initialises max=data[0] before the loop; replicate that.
    double maxP2 = Pfreq[0];
    for( int i = 1; i < n_; i++ )
    {
        if( !is_excluded( i ) && Pfreq[i] > maxP2 )
        {
            maxP2 = Pfreq[i];
        }
    }

    // C/N0 and peak ratio
    // mirrors: cn0 = 10*log10(maxP/meanP/ctime)
    //          peakr = maxP/maxP2
    const double peakr = maxP / maxP2;
    const double cn0   = 10.0 * std::log10( maxP / meanP / ctime_ );

    // Sub-bin Doppler refinement: parabolic interpolation of the power peak across
    // the three Doppler bins at the winning code phase. A 50 Hz grid otherwise
    // leaves +/-25 Hz residual at hand-off, too far for the 2nd-order PLL to pull in
    // without spinning; interpolation tightens it to ~+/-5-10 Hz for free.
    // doppler_freqs_ is relative to the search centre; add the latched recenter offset.
    double doppler = doppler_freqs_[freqi] + doppler_center_bins_ * doppler_step_hz_;
    if( freqi > 0 && freqi < nfreq_ - 1 )
    {
        const double pm1   = power_[( freqi - 1 ) * n_ + codei];
        const double p0    = power_[freqi * n_ + codei];
        const double pp1   = power_[( freqi + 1 ) * n_ + codei];
        const double denom = pm1 - 2.0 * p0 + pp1;
        if( denom < 0.0 ) // a real peak (concave) - guard against flat/noise
        {
            const double delta = 0.5 * ( pm1 - pp1 ) / denom; // in bins, |delta|<=0.5
            doppler += std::clamp( delta, -0.5, 0.5 ) * doppler_step_hz_;
        }
    }

    result_.found      = ( peakr > ACQTH );
    result_.metric     = peakr;
    result_.cn0_db_hz  = cn0;
    result_.code_phase = static_cast<double>( codei );
    result_.doppler_hz = doppler;

    return result_.found;
}

#ifdef ENABLE_UNIT_TESTS
#include <cstdint>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include "cttc_gps_l1_snippet.h"

namespace
{
constexpr double ACQ_FS = 4.0e6; // match the CTTC capture

// Build one acquisition block (m = acq_fft_factor * n complex samples) of a noise-free GPS L1 C/A
// signal: the PRN code (real +/-1) starting at sample offset `code_phase`, on a carrier of doppler_hz.
Complex_buf synth_block( const Complex_buf& code, int m, int code_phase, double doppler_hz )
{
    const int    n  = static_cast<int>( code.size() );
    const double ti = 1.0 / ACQ_FS;
    Complex_buf  blk( m );
    for( int j = 0; j < m; ++j )
    {
        const int    idx = ( ( ( j - code_phase ) % n ) + n ) % n;
        const float  c   = code[idx].real();
        const double ph  = 2.0 * M_PI * doppler_hz * j * ti;
        blk[j]           = Complex_sample( c * std::cos( ph ), c * std::sin( ph ) );
    }
    return blk;
}
} // namespace

// Synthesised (pure-signal) acquisition: a noise-free GPS L1 C/A code on a known carrier at a known
// code phase must be found at exactly that code phase and Doppler. NOTE the receiver's signal-in-Q
// convention flips the reported Doppler sign (an injected +f comes back as -f); see [[iq-convention-flip]].
TEST_CASE( "acquisition_synthetic_gps_signal", "[acquisition][gps]" )
{
    auto              sig  = make_signal( Constellation::Gps, Band::L1, Code::CA );
    const Complex_buf code = sig->code_samples( /*prn=*/1, ACQ_FS );
    const int         m    = sig->params().acq_fft_factor * static_cast<int>( code.size() );

    const int    inject_phase   = 137;
    const double inject_doppler = 2000.0;
    const Complex_buf blk = synth_block( code, m, inject_phase, inject_doppler );

    Acquisition_aiding aiding;
    Acquisition_engine acq( code, 1, ACQ_FS, sig->params(), aiding );
    Sample_block       block { blk.data(), blk.size(), nullptr, 0 };
    for( int i = 0; i < acq.target_integrations(); ++i )
    {
        acq.integrate( block );
    }
    const Acquisition_result r = acq.result();
    REQUIRE( r.found );
    REQUIRE( r.metric > Acquisition_engine::ACQTH );
    REQUIRE( r.code_phase == Catch::Approx( inject_phase ).margin( 1.0 ) );
    REQUIRE( r.doppler_hz == Catch::Approx( -inject_doppler ).margin( 100.0 ) ); // sign-flipped convention
}

// Real-data acquisition: a small hardcoded slice of the CTTC capture (see cttc_gps_l1_snippet.h). A
// present satellite (GPS PRN 1) must acquire at the same Doppler the live receiver reports (-7033 Hz);
// an absent one (PRN 4) must not. One 2 ms block is enough (PRN 1 peak ratio ~6 vs PRN 4 ~1.5).
// TODO: add a Galileo E1 real-data acquisition fixture once the outstanding Galileo issues are sorted.
TEST_CASE( "acquisition_cttc_real_gps_prn1", "[acquisition][gps][cttc]" )
{
    auto      sig = make_signal( Constellation::Gps, Band::L1, Code::CA );
    const int m   = sig->params().acq_fft_factor * static_cast<int>( sig->code_samples( 1, CTTC_SNIPPET_FS ).size() );
    REQUIRE( 2 * m == static_cast<int>( std::size( CTTC_GPS_L1_SNIPPET ) ) ); // one block of I/Q pairs

    Complex_buf samples( m );
    for( int i = 0; i < m; ++i )
    {
        samples[i] = Complex_sample( CTTC_GPS_L1_SNIPPET[2 * i] / 32767.0f, CTTC_GPS_L1_SNIPPET[2 * i + 1] / 32767.0f );
    }
    Sample_block block { samples.data(), samples.size(), nullptr, 0 };

    SECTION( "present satellite acquires" )
    {
        Acquisition_aiding aiding;
        Acquisition_engine acq( sig->code_samples( 1, CTTC_SNIPPET_FS ), 1, CTTC_SNIPPET_FS, sig->params(), aiding );
        REQUIRE( acq.integrate( block ) );
        const Acquisition_result r = acq.result();
        REQUIRE( r.found );
        REQUIRE( r.metric > Acquisition_engine::ACQTH );
        REQUIRE( r.cn0_db_hz > 40.0 );
        REQUIRE( r.doppler_hz == Catch::Approx( -7033.0 ).margin( 250.0 ) );
    }
    SECTION( "absent satellite does not acquire" )
    {
        Acquisition_aiding aiding;
        Acquisition_engine acq( sig->code_samples( 4, CTTC_SNIPPET_FS ), 4, CTTC_SNIPPET_FS, sig->params(), aiding );
        REQUIRE_FALSE( acq.integrate( block ) );
        REQUIRE( acq.result().metric < Acquisition_engine::ACQTH );
    }
}
#endif

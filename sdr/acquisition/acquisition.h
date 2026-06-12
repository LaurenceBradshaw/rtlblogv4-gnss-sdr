#pragma once
#include <fftw3.h>
#include <complex>
#include <vector>
#include "acquisition_aiding.h"
#include "signal.h"
#include "types.h"

// mirrors GNSS-SDRLIB sdracq.c checkacquisition() output fields
struct Acquisition_result
{
    bool   found;
    double metric;     // peak / second-peak ratio  (threshold: ACQTH)
    double cn0_db_hz;  // C/N0 estimate in dB-Hz
    double doppler_hz; // carrier frequency of the winning Doppler bin
    double code_phase; // sample offset of the peak (0 .. nsamp-1)
};

// mirrors GNSS-SDRLIB pcorrelator() + checkacquisition() + sdraccuisition()
class Acquisition_engine
{
public:
    // Doppler search step. GNSS-SDRLIB uses 200 Hz, but that leaves up to +/-100 Hz
    // residual at hand-off - too far for the 2nd-order PLL to pull in without the
    // phase wrapping (slipping) first. 50 Hz keeps the residual within +/-25 Hz so
    // the carrier loop locks immediately instead of slipping during pull-in.
    static constexpr double ACQSTEP   = 200.0;   // Hz - Doppler search step
    static constexpr double ACQTH     = 3.0;     // peak-ratio threshold            (ACQTH)
    // Doppler search half-width is now supplied per-attempt by Acquisition_aiding (wide while
    // bootstrapping, narrow once the common offset is pinned); the grid is allocated for WIDE.

    // prn_code : one code period, Q=0, as returned by Signal::code_samples()
    // sample_rate_hz : front-end sample rate
    // sig : signal physics (chip rate sets the exclusion zone; code period sets the
    //       C/N0 normalisation). PRN-independent - shared across this signal's channels.
    // aiding : receiver-wide Doppler-recenter estimate, read at the start of each attempt.
    Acquisition_engine( const Complex_buf& prn_code, double sample_rate_hz, const Signal_params& sig,
                        const Acquisition_aiding& aiding );
    ~Acquisition_engine();

    // FFTW plans are not copyable
    Acquisition_engine( const Acquisition_engine& )            = delete;
    Acquisition_engine& operator=( const Acquisition_engine& ) = delete;

    // Correlate one 2*nsamp block and accumulate into the internal power map.
    // Mirrors pcorrelator() with flagsum=1, then calls check_acquisition().
    // Returns true when peak/second-peak ratio exceeds ACQTH.
    bool integrate( const Sample_block& block );

    // Clear the accumulated power map and integration counter.
    // Call after a successful acquisition or when sliding to a new search window.
    void reset();

    int intg_count() const
    {
        return intg_count_;
    }
    // Non-coherent epochs to accumulate before the accept/reject decision (per-signal).
    int target_integrations() const
    {
        return target_integrations_;
    }
    const Acquisition_result& result() const
    {
        return result_;
    }

private:
    // mirrors GNSS-SDRLIB checkacquisition()
    bool check_acquisition();

    // FFT / IFFT plans (size m_, reused every integrate() call)
    fftwf_plan fft_plan_;
    fftwf_plan ifft_plan_;

    Complex_buf code_fft_; // conj(FFT(zero-padded code)) (size m_) - xcode in GNSS-SDRLIB
    Complex_buf data_fft_; // FFT of the raw block, computed once per integrate (size m_)
    Complex_buf datax_;    // scratch: shifted spectrum x conj(code), then IFFT (size m_)

    std::vector<double> doppler_freqs_; // nfreq_ Doppler hypotheses (Hz, relative to recenter)
    std::vector<double> power_;         // n_ x nfreq_ accumulated correlation power - P in GNSS-SDRLIB

    const Acquisition_aiding& aiding_;          // shared receiver-wide recenter estimate
    double                    carrier_freq_hz_; // this signal's carrier (Hz) - for the aiding query
    int    doppler_center_bins_; // recenter offset (integer FFT bins), latched per attempt
    int    active_half_bins_;    // searched half-width (bins each side of centre), latched per attempt

    Acquisition_result result_;
    int                intg_count_;
    int                target_integrations_; // non-coherent epochs per attempt (sig.acq_integrations)

    int    n_;               // nsamp - samples per code period
    int    m_;               // nfft = acq_fft_factor * n_ - FFT / correlation size
    int    nfreq_;           // number of Doppler bins
    int    nsamp_chip_;      // samples per chip - sets exclusion zone in checkacquisition
    double chip_rate_;       // spreading-code chip rate (chip/s) - from Signal_params
    double ti_;              // sampling interval (s)
    double ctime_;           // code period (s) - C/N0 normalisation divisor (from Signal_params)
    double doppler_step_hz_; // Doppler bin spacing = fs/m_ (FFT-bin aligned, ~500 Hz)
};

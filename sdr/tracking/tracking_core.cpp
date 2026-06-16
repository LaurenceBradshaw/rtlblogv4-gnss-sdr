#include "tracking_core.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include "logging.h"

// constructor
// mirrors GNSS-SDRLIB inittrkstruct() + inittrkprmstruct():
//
//   trk->codefreq = crate
//   trk->carrfreq = acq.acqfreq       (set later in initialise())
//   trk->loop     = LOOP_L1CA
//   trk->corrn    = sdrini.trkcorrn   (hardcoded to CORRN=1)
//   trk->corrp[0] = sdrini.trkcorrd   (hardcoded to ~0.5-chip spacing)
//   trk->ne=1, trk->nl=2
//   trk->II/QQ/... = calloc(...)
//   prm1/prm2 coefficients from bandwidths
Tracking_core::Tracking_core(
    const std::vector<float>& code_chips,
    double                    sample_rate_hz,
    const Signal_params&      sig,
    const std::vector<float>& secondary,
    const std::vector<float>& data_code
)
    : code_( code_chips ),
      code_len_( static_cast<int>( code_chips.size() ) ),
      code_rate_( static_cast<double>( code_chips.size() ) / sig.code_period_s ),
      epoch_period_( sig.code_period_s ),
      fll_active_( sig.nav_bit_ms > static_cast<int>( std::round( sig.code_period_s * 1000.0 ) ) ),
      rf_freq_( sig.carrier_freq_hz ),
      sample_rate_( sample_rate_hz ),
      ti_( 1.0 / sample_rate_hz ),
      acq_freq_( 0.0 ),
      code_freq_( static_cast<double>( code_chips.size() ) / sig.code_period_s ),
      carrier_freq_( 0.0 ),
      remaining_code_( 0.0 ),
      remaining_carr_( 0.0 ),
      code_nco_( 0.0 ),
      code_err_( 0.0 ),
      carrier_nco_( 0.0 ),
      carrier_acc_( 0.0 ),
      carrier_err_( 0.0 ),
      freq_err_( 0.0 ),
      data_code_( data_code ),
      has_data_( !data_code.empty() ),
      data_prompt_i_( 0.0 ),
      data_prompt_q_( 0.0 ),
      old_data_prompt_i_( 0.0 ),
      old_data_prompt_q_( 0.0 ),
      epoch_count_( 0 ),
      secondary_( secondary ),
      secondary_sync_( false ),
      secondary_index_( 0 ),
      secondary_polarity_( 1 )
{
    // Loop bandwidths come from the Signal (Signal_params::loop_bw_wide/narrow): per-signal because
    // the right values depend on the integration period and C/N0. wide = pull-in (prm1, pre frame
    // sync); narrow = tight steady state (prm2, post sync). Both run at the actual code-period dt.
    prm1_ = make_prm( sig.loop_bw_wide.dll, sig.loop_bw_wide.pll, sig.loop_bw_wide.fll );
    prm2_ = make_prm( sig.loop_bw_narrow.dll, sig.loop_bw_narrow.pll, sig.loop_bw_narrow.fll );

    // E/L spacing in samples: target ~0.5 chip
    // mirrors sdrini.trkcorrp (default 1 sample for RTL-SDR at 2 MHz)
    corr_spacing_ = std::max( 1, static_cast<int>( std::round( sample_rate_hz / ( 2.0 * code_rate_ ) ) ) );

    // Code elements per ranging chip: 1 for BPSK, 2 for BOC half-chips. Lets the pilot place
    // VEML taps at chip-fraction offsets (e.g. 0.6 chip) regardless of the half-chip coding.
    code_elements_per_chip_ = code_rate_ / sig.chip_rate_hz;

    n_taps_ = 3; // default (Costas P/E/L); configure_taps() sets it each epoch
    tap_offset_chips_.fill( 0.0 );

    II_.fill( 0.0 );
    QQ_.fill( 0.0 );
    old_I_.fill( 0.0 );
    old_Q_.fill( 0.0 );
    sum_I_.fill( 0.0 );
    sum_Q_.fill( 0.0 );
    oldsum_I_.fill( 0.0 );
    oldsum_Q_.fill( 0.0 );
}

// make_prm
// mirrors GNSS-SDRLIB inittrkprmstruct():
//   dllw2 = (dllb/0.53)^2   dllaw = 1.414*(dllb/0.53)
//   pllw2 = (pllb/0.53)^2   pllaw = 1.414*(pllb/0.53)
//   fllw  = fllb/0.25
Tracking_loop_prm Tracking_core::make_prm( double dllb, double pllb, double fllb )
{
    Tracking_loop_prm p;
    p.dllb  = dllb;
    p.pllb  = pllb;
    p.fllb  = fllb;
    p.dllw2 = ( dllb / 0.53 ) * ( dllb / 0.53 );
    p.dllaw = 1.414 * ( dllb / 0.53 );
    p.pllw2 = ( pllb / 0.53 ) * ( pllb / 0.53 );
    p.pllaw = 1.414 * ( pllb / 0.53 );
    p.fllw  = fllb / 0.25;
    return p;
}

// initialise
// Called once after a successful acquisition.
// mirrors: trk->carrfreq = acq.acqfreq; trk->codefreq = sdr->crate; trk->remcode/carr = 0
// remcode is 0 because channel.cpp already advanced next_sample_ to the code boundary.
void Tracking_core::initialise( const Acquisition_result& acq )
{
    // Negate the acquisition Doppler at hand-off: the wipe is now conj (e^{-j}), so the same
    // physical signal needs the opposite-sign carrier_freq_ to be cancelled. Acquisition itself
    // (and the aiding bridge it feeds) keep their own unchanged convention.
    acq_freq_       = -acq.doppler_hz;
    carrier_freq_   = -acq.doppler_hz;
    code_freq_      = code_rate_;
    remaining_code_ = 0.0;
    remaining_carr_ = 0.0;
    code_nco_       = 0.0;
    code_err_       = 0.0;
    carrier_nco_    = 0.0;
    carrier_acc_    = 0.0;
    carrier_err_    = 0.0;
    freq_err_       = 0.0;
    epoch_count_    = 0;

    // Reset the lock detector so a re-acquired channel starts measuring afresh.
    m2_sum_            = 0.0;
    m4_sum_            = 0.0;
    nbd_sum_           = 0.0;
    cn0_count_         = 0;
    cn0_db_hz_         = 0.0;
    carrier_lock_test_ = 0.0;
    lock_seeded_       = false;
    lock_fail_count_   = 0;
    locked_            = false;

    secondary_sync_     = false;
    secondary_index_    = 0;
    secondary_polarity_ = 1;
    sec_i_hist_.clear();
    sec_q_hist_.clear();

    data_prompt_i_ = data_prompt_q_ = old_data_prompt_i_ = old_data_prompt_q_ = 0.0;

    II_.fill( 0.0 );
    QQ_.fill( 0.0 );
    old_I_.fill( 0.0 );
    old_Q_.fill( 0.0 );
    sum_I_.fill( 0.0 );
    sum_Q_.fill( 0.0 );
    oldsum_I_.fill( 0.0 );
    oldsum_Q_.fill( 0.0 );
}

// compute_nsamp
// mirrors: sdr->samples_consumed = (int)((sdr->clen - sdr->trk.remcode) / (sdr->trk.codefreq / sdr->f_sf))
int Tracking_core::compute_samples_needed() const
{
    // remaining_code_ is kept in (-(ci/2), code_len_ - ci/2] by the pre-wrap in correlate(),
    // so (code_len_ - remaining_code_) is always >= ci/2 > 0, giving n >= 1 here.
    // The max(1, n) guard is still here as a safety net.
    const int n = static_cast<int>( ( code_len_ - remaining_code_ ) / ( code_freq_ / sample_rate_ ) );
    return std::max( 1, n );
}

// correlate
// Time-domain Early/Prompt/Late correlator.
// mirrors GNSS-SDRLIB correlator() in sdrcmn.c:
//   *remp = mixcarr(data, dtype, ti, n, freq, phi0, dataI, dataQ)
//   *remc = rescode(codein, coden, coff, smax, ti*crate, n, code_e)
//   dot_23(dataI, dataQ, code, code-s[0], code+s[0], n, II, QQ)
//
// Carrier: wipe by conj(e^(+j * 2pi * carrfreq * ti * i + j * remcarr)) = e^(-j...) (standard
// wipeoff). carrier_freq_ is set to -acq.doppler_hz at hand-off (initialise) so this conjugated
// wipe still cancels the signal carrier - acquisition itself is unchanged.
// Code array indexing: Prompt=index 0, Early=index 1 (code-s), Late=index 2 (code+s)
void Tracking_core::correlate( const Sample_block& block, int n )
{
    // Save previous epoch correlations before overwriting
    // mirrors: memcpy(trk->oldI, trk->II, ...); trk->oldremcode = trk->remcode; ...
    old_I_             = II_;
    old_Q_             = QQ_;
    old_data_prompt_i_ = data_prompt_i_;
    old_data_prompt_q_ = data_prompt_q_;

    const double ci = code_freq_ / sample_rate_; // code elements per sample

    // Strategy sets n_taps_ + the per-tap code-element offsets (Costas: 3-tap E/P/L; Pilot:
    // 5-tap VE/E/P/L/VL). Prompt is always tap 0.
    configure_taps( ci );

    // Carrier phasor - mirrors mixcarr() with DTYPEIQ. Done in plain float (not
    // std::complex) so the hot loop has no out-of-line operator*/real()/imag() calls
    // at -O0; auto-vectorizes at -O3. Math is identical: float carrier recurrence,
    // double E/P/L accumulators.
    //   phi step per sample: ps = freq * 2pi * ti;  w = conj(e^{+j*phi})*raw = e^{-j*phi}*raw
    // remaining_carr_ is kept wrapped to [0, 2pi) below so the start cos/sin are exact.
    const double step_phase = 2.0 * M_PI * carrier_freq_ * ti_;
    const float  sr         = static_cast<float>( std::cos( step_phase ) );      // step real
    const float  si         = static_cast<float>( std::sin( step_phase ) );      // step imag
    float        cr         = static_cast<float>( std::cos( remaining_carr_ ) ); // carrier real
    float        cqi        = static_cast<float>( std::sin( remaining_carr_ ) ); // carrier imag

    // std::complex<float> is layout-compatible with float[2]; access as raw floats.
    const float* p1   = reinterpret_cast<const float*>( block.ptr1 );
    const float* p2   = reinterpret_cast<const float*>( block.ptr2 );
    const int    len1 = static_cast<int>( block.len1 );

    double acc_i[MAX_TAPS] = { 0.0 };
    double acc_q[MAX_TAPS] = { 0.0 };
    double data_i = 0.0, data_q = 0.0; // data-component prompt (pilot tracking)

    for( int i = 0; i < n; i++ )
    {
        float rr, ri; // raw sample real/imag
        if( i < len1 )
        {
            rr = p1[2 * i];
            ri = p1[2 * i + 1];
        }
        else
        {
            const int k = i - len1;
            rr          = p2[2 * k];
            ri          = p2[2 * k + 1];
        }

        // Carrier wipe: w = raw * conj(carrier) = raw * e^{-j*phi} (standard wipeoff). The
        // carrier replica is (cr + j*cqi); conjugating it lands the despread SIGNAL in I (real)
        // and the phase ERROR in Q (imag) - the textbook signal-in-I convention.
        const float wr = rr * cr + ri * cqi;
        const float wi = ri * cr - rr * cqi;
        // Advance carrier: carrier *= step
        const float ncr = cr * sr - cqi * si;
        cqi             = cr * si + cqi * sr;
        cr              = ncr;

        // Prompt code phase (code elements) with wrap-around.
        double phase_P = remaining_code_ + i * ci;
        if( phase_P >= code_len_ )
        {
            phase_P -= code_len_;
        }
        else if( phase_P < 0.0 ) // carried phase can be slightly negative (see end of loop)
        {
            phase_P += code_len_;
        }

        // Each tap = prompt phase + its offset; the offsets are small (< 1 period) so a
        // single conditional wrap suffices (for the 3-tap case this is bit-identical to the
        // previous explicit P/E/L code with phase_E = P-s, phase_L = P+s).
        for( int k = 0; k < n_taps_; k++ )
        {
            double ph = phase_P + tap_offset_chips_[k];
            if( ph < 0.0 )
            {
                ph += code_len_;
            }
            else if( ph >= code_len_ )
            {
                ph -= code_len_;
            }
            const float c = code_[static_cast<int>( ph )];
            acc_i[k] += wr * c;
            acc_q[k] += wi * c;
        }

        // Data-component prompt: despread the data code at the prompt phase (pilot tracking).
        if( has_data_ )
        {
            const float cd = data_code_[static_cast<int>( phase_P )];
            data_i += wr * cd;
            data_q += wi * cd;
        }
    }

    for( int k = 0; k < n_taps_; k++ )
    {
        II_[k] = acc_i[k];
        QQ_[k] = acc_q[k];
    }
    data_prompt_i_ = data_i;
    data_prompt_q_ = data_q;

    // Update carried code phase for the next epoch.
    //
    // n = floor((code_len_ - remaining_code_)/ci), so after advancing n samples the
    // phase lands in (code_len_ - ci, code_len_]. Subtract one full period so the
    // carried phase stays small (in (-ci, 0]) and the NEXT compute_samples_needed()
    // always returns ~one full code period (~nsamp), never a 1-sample epoch.
    //
    // The previous conditional wrap only fired when the leftover was < 0.5*ci, so a
    // leftover in [0.5*ci, ci) left remaining_code_ just under code_len_ and produced
    // a 1-sample epoch (near-zero correlation -> discriminator spike that kicked the
    // loop). Subtracting unconditionally removes that failure mode.
    remaining_code_ += n * ci - code_len_;

    remaining_carr_ += n * step_phase;
    remaining_carr_ = std::fmod( remaining_carr_, 2.0 * M_PI );
    if( remaining_carr_ < 0.0 )
    {
        remaining_carr_ += 2.0 * M_PI;
    }
}

// cumsum_corr
// mirrors GNSS-SDRLIB cumsumcorr():
//   II[i] *= polarity;   QQ[i] *= polarity;
//   oldsumI[i] += oldI[i];   oldsumQ[i] += oldQ[i];
//   sumI[i]    += II[i];     sumQ[i]    += QQ[i];
void Tracking_core::cumsum_corr( int polarity )
{
    for( int i = 0; i < n_taps_; i++ )
    {
        II_[i] *= polarity;
        QQ_[i] *= polarity;

        oldsum_I_[i] += old_I_[i];
        oldsum_Q_[i] += old_Q_[i];
        sum_I_[i] += II_[i];
        sum_Q_[i] += QQ_[i];
    }
}

// clear_cumsum
// mirrors GNSS-SDRLIB clearcumsumcorr()
void Tracking_core::clear_cumsum()
{
    oldsum_I_.fill( 0.0 );
    oldsum_Q_.fill( 0.0 );
    sum_I_.fill( 0.0 );
    sum_Q_.fill( 0.0 );
}

// Pure lock-detector math, factored out so it can be unit-tested directly (the estimator now feeds both
// the GUI C/N0 display and the PVT observable gate, so it is worth pinning).
namespace
{
// M2M4 (2nd/4th-moment) C/N0 estimate (dB-Hz) from the window-averaged prompt-power moments
// M2 = <P>, M4 = <P^2> (P = I^2+Q^2). Signal power Pd = sqrt(2*M2^2 - M4), noise Pn = M2 - Pd,
// SNR = Pd/Pn, C/N0 = 10*log10(SNR / Tcoh). Returns 0 when the moments are inconsistent with a
// coherent signal (no positive SNR). Scale-invariant: scaling P by k leaves the SNR (a ratio) fixed.
double m2m4_cn0_db_hz( double m2, double m4, double epoch_period_s )
{
    const double pd2 = 2.0 * m2 * m2 - m4; // signal power squared (negative when no coherent signal)
    if( pd2 <= 0.0 || m2 <= 0.0 )
    {
        return 0.0;
    }
    const double pd = std::sqrt( pd2 );
    const double pn = m2 - pd;
    if( pn <= 0.0 )
    {
        return 0.0;
    }
    return 10.0 * std::log10( ( pd / pn ) / epoch_period_s );
}

// Van Dierendonck carrier lock test cos(2*phi) ~ NBD/NBP, where NBD = sum(I^2 - Q^2) and
// NBP = sum(I^2 + Q^2) over the window. ~+1 phase-locked (energy on I), <=0 carrier lost.
double carrier_lock_cos2phi( double nbd_sum, double nbp_sum )
{
    return ( nbp_sum > 0.0 ) ? nbd_sum / nbp_sum : 0.0;
}
} // namespace

// update_lock_detectors
// Two complementary detectors over a window of CN0_WINDOW_EPOCHS, then a hysteretic lock decision.
//
// (1) M2M4 C/N0 on the prompt power P = I^2 + Q^2 (see m2m4_cn0_db_hz). P is sign-insensitive (data
//     flips don't matter) and the SNR is a ratio (scale-invariant), so this is identical for every
//     signal; only epoch_period_ sets the dB-Hz offset.
// (2) Van Dierendonck carrier lock test cos(2*phi) ~ (<I^2> - <Q^2>) / (<I^2> + <Q^2>): ~+1 when the
//     carrier is phase-locked (all energy on I), <=0 when it is lost. Power-difference form (data-robust;
//     see the header notes on why not the coherent (sum I)^2 one).
// Both are EMA-smoothed. has_lock() is then latched: lock requires BOTH above threshold; loss is declared
// only after LOCK_FAIL_WINDOWS consecutive failing windows; a single good window re-locks.
void Tracking_core::update_lock_detectors( double prompt_i, double prompt_q )
{
    const double power = prompt_i * prompt_i + prompt_q * prompt_q;
    m2_sum_ += power;
    m4_sum_ += power * power;
    nbd_sum_ += prompt_i * prompt_i - prompt_q * prompt_q;

    if( ++cn0_count_ < CN0_WINDOW_EPOCHS )
    {
        return;
    }

    // (1) M2M4 C/N0. (2) Carrier lock test (NBP = sum I^2+Q^2 = m2_sum_; NBD = sum I^2-Q^2 = nbd_sum_).
    const double cn0       = m2m4_cn0_db_hz( m2_sum_ / cn0_count_, m4_sum_ / cn0_count_, epoch_period_ );
    const double lock_test = carrier_lock_cos2phi( nbd_sum_, m2_sum_ );

    // EMA-smooth both (seed on the first window so they converge quickly).
    if( !lock_seeded_ )
    {
        cn0_db_hz_         = cn0;
        carrier_lock_test_ = lock_test;
        lock_seeded_       = true;
    }
    else
    {
        cn0_db_hz_         = 0.7 * cn0_db_hz_ + 0.3 * cn0;
        carrier_lock_test_ = 0.7 * carrier_lock_test_ + 0.3 * lock_test;
    }

    // Hysteretic lock decision.
    const bool window_ok = ( cn0_db_hz_ > CN0_LOCK_THRESHOLD_DBHZ ) && ( carrier_lock_test_ > CARRIER_LOCK_THRESHOLD );
    if( window_ok )
    {
        locked_          = true;
        lock_fail_count_ = 0;
    }
    else if( ++lock_fail_count_ >= LOCK_FAIL_WINDOWS )
    {
        locked_ = false;
    }

    m2_sum_    = 0.0;
    m4_sum_    = 0.0;
    nbd_sum_   = 0.0;
    cn0_count_ = 0;
}

// pll_update
// 3rd-order PLL aided by a 1st-order FLL - the structure GNSS-SDR uses for this
// file (order=3, FLL-assisted). A 3rd-order loop tracks frequency RATE, so it
// rides through data-bit transitions and line-of-sight dynamics that make the
// 2nd-order loop thrash. Discriminators are unchanged (Costas + cross-FLL).
//
// Standard Kaplan 3rd-order coefficients (a3=1.1, b3=2.4, w0=Bn/0.7845), with two
// integrators carried in carrier_acc_ (acceleration) and carrier_nco_ (frequency):
//
//   carrier_acc_ += w0^3*theta*dt
//   carrier_nco_ += b3*w0*dtheta + a3*w0^2*theta*dt + carrier_acc_*dt + fll*thetaf*dt
//   carrfreq      = acqfreq + carrier_nco_
void Tracking_core::pll_update( const Tracking_loop_prm& prm, double dt, bool use_fll, bool pure_pll )
{
    // Standard signal-in-I convention (after the conj wipe): the despread signal is in I (sum_I)
    // and the phase error in Q (sum_Q). atan2(QP, IP) = atan2(error, signal) drives the prompt
    // onto the +I axis. The accumulator sign (carrier_nco_ += ...) is unchanged: swapping the
    // discriminator inputs already accounts for the wipe flip, leaving negative feedback intact.
    const double IP = sum_I_[0], QP = sum_Q_[0];
    const double oldIP = oldsum_I_[0], oldQP = oldsum_Q_[0];

    double carrErr, freqErr = 0.0;

    // PLL discriminator. pure_pll: full 4-quadrant atan2 (GNSS-SDR pll_four_quadrant_atan)
    // for a data-free pilot after secondary wipeoff - resolves the full 360 deg, no half-cycle
    // ambiguity, tighter lock. Otherwise Costas (180 deg-folded), for data-bearing signals.
    if( pure_pll )
    {
        carrErr = std::atan2( QP, IP ) / M_PI;
    }
    else if( IP > 0.0 )
    {
        carrErr = std::atan2( QP, IP ) / M_PI;
    }
    else
    {
        carrErr = std::atan2( -QP, -IP ) / M_PI;
    }

    // FLL discriminator. Two data/secondary-insensitive forms; pick per signal. The conj wipe
    // reverses the prompt's rotation sense, so the frequency-error sign is NEGATED vs the old
    // signal-in-Q code. (The PLL absorbed its flip through the I/Q-input swap; the FLL reads
    // sum_I/sum_Q directly, so it needs the explicit negation.)
    const double Ic = sum_I_[0], Qc = sum_Q_[0];
    const double Ip = oldsum_I_[0], Qp = oldsum_Q_[0];
    if( fll_active_ )
    {
        // Cross/dot of consecutive prompts (GPS). Valid when the coherent sum spans <1
        // data bit; the data sign (d_prev*d_curr) cancels in every term.
        //   cross = I_prev*Q_curr - I_curr*Q_prev,  dot = I_prev*I_curr + Q_prev*Q_curr
        const double cross = Ip * Qc - Ic * Qp;
        const double dot   = Ip * Ic + Qp * Qc;
        freqErr            = -std::atan2( cross, dot );
    }
    else
    {
        // atan(Q/I) difference (Galileo E1, 1 symbol/epoch - GNSS-SDRLIB form, in the signal-in-I
        // convention). atan() folds +/-pi, so a +/-1 data/secondary flip (negates I and Q) leaves
        // Q/I unchanged -> it survives the modulation every epoch without a multi-epoch sum.
        // UNVERIFIED: Galileo is disabled in main.cpp; re-check this sign/form when re-enabling it.
        const double f1 = ( Ic == 0.0 ) ? M_PI / 2 : std::atan( Qc / Ic );
        const double f2 = ( Ip == 0.0 ) ? M_PI / 2 : std::atan( Qp / Ip );
        freqErr         = -( f1 - f2 );
        if( freqErr > M_PI / 2 )
        {
            freqErr = M_PI - freqErr;
        }
        if( freqErr < -M_PI / 2 )
        {
            freqErr = -M_PI - freqErr;
        }
    }

    const double fll_term = use_fll ? prm.fllw * dt * freqErr : 0.0;
    carrier_nco_ += prm.pllaw * ( carrErr - carrier_err_ ) + prm.pllw2 * dt * carrErr + fll_term;

    carrier_freq_ = acq_freq_ + carrier_nco_;
    carrier_err_  = carrErr;
    freq_err_     = freqErr;
}

// dll_update
// 2nd-order DLL with carrier aiding.
// mirrors GNSS-SDRLIB dll() in sdrtrk.c exactly:
//
//   codeErr = (|E| - |L|) / (|E| + |L|)    normalised envelope discriminator
//   codeNco += dllaw*(codeErr-prevErr) + dllw2*dt*codeErr
//   codefreq = crate - codeNco + (carrfreq - f_if - foffset) / (f_cf / crate)
//            = crate - codeNco + carrfreq * (chip_rate / rf_freq)   [f_if=foffset=0]
void Tracking_core::dll_update( const Tracking_loop_prm& prm, double dt )
{
    // Strategy supplies the normalised code discriminator (Costas: E-L; Pilot: VEML).
    const double codeErr = code_error();

    code_nco_ += prm.dllaw * ( codeErr - code_err_ ) + prm.dllw2 * dt * codeErr;

    // Carrier-aided code frequency - mirrors: crate - codeNco + (carrfreq - f_if - foffset)/(f_cf/crate)
    code_freq_ = code_rate_ - code_nco_ + carrier_freq_ * ( code_rate_ / rf_freq_ );
    code_err_  = codeErr;
}

// correlate_epoch
// First half of a tracking epoch: correlate and capture the per-epoch prompts.
//
// This is split from run_loops() so the caller can run navigation in between -
// matching GNSS-SDRLIB's actual per-epoch order, where sdrnavigation() runs
// after correlator() but BEFORE the loop filters, so flagsync/swloop are current:
//
//   correlator()        <- here
//   sdrnavigation()     <- caller, sets flagsync/swloop for THIS epoch
//   cumsumcorr()        <- run_loops()
//   pll()/dll()         <- run_loops(), gated by the just-updated flags
//   clearcumsumcorr()   <- run_loops()
Tracking_output Tracking_core::correlate_epoch( const Sample_block& block )
{
    const int n = compute_samples_needed();

    if( static_cast<size_t>( n ) > block.len1 + block.len2 )
    {
        return { false, 0, 0.0, 0.0 };
    }

    correlate( block, n );

    // Lock detector: feed this epoch's prompt I/Q to the C/N0 + carrier-lock detectors (the same call
    // serves Costas and pilot tracking alike).
    update_lock_detectors( II_[0], QQ_[0] );

    double cur_I, prev_I;
    if( has_data_ && secondary_sync_ )
    {
        // Pilot-aided data demodulation, pilot SECONDARY SYNCED. The data-component prompt and the pilot
        // share the carrier, so the nav symbol is the data prompt PROJECTED onto the pilot's phasor:
        //   symbol = (data . pilot) / |pilot|.
        // This recovers the symbol regardless of which axis the carrier loop settles on (and tracks
        // residual phase wander), unlike reading a fixed I or Q. The pilot's per-epoch secondary chip is
        // divided out so it doesn't flip the symbol; the constant E1-B/E1-C sign is left for the nav
        // decoder's preamble-polarity resolution. (mirrors GNSS-SDR's pilot-referenced data prompt.)
        const double cs25 = !secondary_.empty() ? ( secondary_[secondary_index_] >= 0.0f ? 1.0 : -1.0 ) : 1.0;
        const double mag  = std::sqrt( II_[0] * II_[0] + QQ_[0] * QQ_[0] ) + 1e-12;
        const double omag = std::sqrt( old_I_[0] * old_I_[0] + old_Q_[0] * old_Q_[0] ) + 1e-12;
        cur_I             = cs25 * ( data_prompt_i_ * II_[0] + data_prompt_q_ * QQ_[0] ) / mag;
        prev_I            = cs25 * ( old_data_prompt_i_ * old_I_[0] + old_data_prompt_q_ * old_Q_[0] ) / omag;
    }
    else if( has_data_ )
    {
        // Pilot tracking but the secondary is NOT synced - e.g. GPS L1C, whose 1800-chip L1Co overlay is
        // not synced (its clean Costas lock doesn't need it). The pilot prompt still CARRIES that overlay,
        // so projecting the data onto it would inject the overlay's pseudorandom +/-1 and scramble every
        // symbol. The data component (L1Cd) has NO overlay, so read its carrier-wiped I directly: the
        // carrier loop holds the despread signal on I (signal-in-I conj wipe), and the constant polarity
        // is resolved by the nav decoder's frame-sync. (For Galileo this branch is only the brief
        // pre-CS25-sync window, whose symbols are not used for decode.)
        cur_I  = data_prompt_i_;
        prev_I = old_data_prompt_i_;
    }
    else
    {
        // Directly-tracked (Costas) signal: with the conj (e^{-j}) wipe the data lands in I.
        cur_I  = II_[0];
        prev_I = old_I_[0];
    }

    return { true, n, cur_I, prev_I };
}

double Tracking_core::get_fractional_chip_time() const
{
    double chips_passed = static_cast<double>( code_len_ ) - remaining_code_;

    return chips_passed / code_rate_;
}

// advance_secondary_sync
// Find the secondary (overlay) code phase from the recent COMPLEX prompt history, using a
// DIFFERENTIAL correlation that is robust to an imperfect carrier lock.
//
// Why differential: before the secondary is wiped the carrier loop can only Costas/FLL-track
// the pilot, leaving a residual frequency that slowly rotates the prompt phase. A coherent
// (single-component) correlation against the secondary then cancels out across the window
// (the phase - and hence the projected sign - drifts), so it never reaches threshold even for
// a strong, "locked" SV. Forming D[k] = P[k]*conj(P[k-1]) cancels the common carrier phase and
// leaves the differential secondary pattern s[k]*s[k-1] times a CONSTANT residual-frequency
// phase, which the magnitude below is immune to. (Differential is blind to the overall code
// polarity - that's fine: the nav decoder resolves data sign from the preamble.)
void Tracking_core::advance_secondary_sync()
{
    const int S = static_cast<int>( secondary_.size() );
    sec_i_hist_.push_back( static_cast<float>( II_[0] ) );
    sec_q_hist_.push_back( static_cast<float>( QQ_[0] ) );

    const int W_MIN = 6 * S;  // need several secondary periods of history before searching
    const int W_MAX = 12 * S; // bound the window (and the search cost)
    if( static_cast<int>( sec_q_hist_.size() ) > W_MAX )
    {
        sec_i_hist_.erase( sec_i_hist_.begin() );
        sec_q_hist_.erase( sec_q_hist_.begin() );
    }
    if( static_cast<int>( sec_q_hist_.size() ) < W_MIN )
    {
        return;
    }

    const int W = static_cast<int>( sec_q_hist_.size() );

    // Differential products D[k] = P[k] * conj(P[k-1]) and the total differential energy.
    std::vector<double> d_re( W ), d_im( W );
    double              energy = 0.0;
    for( int k = 1; k < W; k++ )
    {
        const double ik = sec_i_hist_[k], qk = sec_q_hist_[k];
        const double ip = sec_i_hist_[k - 1], qp = sec_q_hist_[k - 1];
        d_re[k] = ik * ip + qk * qp;
        d_im[k] = qk * ip - ik * qp;
        energy += std::sqrt( d_re[k] * d_re[k] + d_im[k] * d_im[k] );
    }

    // Correlate the differential prompt against the differential secondary at every cyclic
    // offset; the peak magnitude reveals the phase.
    double best     = 0.0;
    int    best_off = 0;
    for( int off = 0; off < S; off++ )
    {
        double cr = 0.0, ci = 0.0;
        for( int k = 1; k < W; k++ )
        {
            const float  s0   = secondary_[( k + off ) % S] >= 0.0f ? 1.0f : -1.0f;
            const float  s1   = secondary_[( k - 1 + off ) % S] >= 0.0f ? 1.0f : -1.0f;
            const double dsec = s0 * s1; // differential secondary chip
            cr += d_re[k] * dsec;
            ci += d_im[k] * dsec;
        }
        const double mag = std::sqrt( cr * cr + ci * ci );
        if( mag > best )
        {
            best     = mag;
            best_off = off;
        }
    }

    // Accept when the aligned differential correlation captures a clear majority of the
    // differential energy. (Differential squares the noise, so the bar is a bit below the old
    // coherent 0.5.)
    if( best >= SECONDARY_SYNC_RATIO * energy )
    {
        secondary_sync_     = true;
        secondary_polarity_ = 1; // overall polarity unresolved (differential); nav resolves it
        // sec_*_hist_[W-1] is THIS epoch (secondary phase (W-1+best_off)%S); the next epoch
        // (first one wiped in run_loops) is one chip later.
        secondary_index_ = ( W + best_off ) % S;
        clear_cumsum(); // restart coherent accumulation cleanly post-sync
    }
}

#ifdef ENABLE_UNIT_TESTS
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

namespace
{
// Closed-form prompt-power moments for a coherent signal of power Pd in complex Gaussian noise of
// power Pn: M2 = Pd + Pn, M4 = Pd^2 + 4*Pd*Pn + 2*Pn^2. The M2M4 estimator inverts these exactly.
void moments_for( double pd, double pn, double& m2, double& m4 )
{
    m2 = pd + pn;
    m4 = pd * pd + 4.0 * pd * pn + 2.0 * pn * pn;
}
} // namespace

TEST_CASE( "m2m4_cn0_recovers_known_snr", "[tracking][cn0]" )
{
    // Feed moments synthesised from a known Pd/Pn -> the estimator must return 10*log10(SNR/T).
    struct
    {
        double pd, pn, t;
    } cases[] = { { 100.0, 1.0, 1e-3 }, { 50.0, 2.0, 1e-3 }, { 10.0, 5.0, 4e-3 }, { 1000.0, 7.0, 1e-3 } };
    for( const auto& c : cases )
    {
        double m2 = 0.0, m4 = 0.0;
        moments_for( c.pd, c.pn, m2, m4 );
        const double expected = 10.0 * std::log10( ( c.pd / c.pn ) / c.t );
        INFO( "Pd=" << c.pd << " Pn=" << c.pn << " T=" << c.t );
        REQUIRE( m2m4_cn0_db_hz( m2, m4, c.t ) == Catch::Approx( expected ) );
    }
}

TEST_CASE( "m2m4_cn0_is_scale_invariant", "[tracking][cn0]" )
{
    double m2 = 0.0, m4 = 0.0;
    moments_for( 80.0, 3.0, m2, m4 );
    const double base = m2m4_cn0_db_hz( m2, m4, 1e-3 );
    // Scaling the prompt power by k scales M2 by k and M4 by k^2; the SNR (a ratio) is unchanged.
    for( double k : { 0.01, 4.0, 1000.0 } )
    {
        REQUIRE( m2m4_cn0_db_hz( k * m2, k * k * m4, 1e-3 ) == Catch::Approx( base ) );
    }
}

TEST_CASE( "m2m4_cn0_zero_when_no_coherent_signal", "[tracking][cn0]" )
{
    // Pure noise (Pd=0): M2=Pn, M4=2*Pn^2 -> 2*M2^2 - M4 = 0 -> no positive SNR -> 0 dB-Hz.
    double m2 = 0.0, m4 = 0.0;
    moments_for( 0.0, 5.0, m2, m4 );
    REQUIRE( m2m4_cn0_db_hz( m2, m4, 1e-3 ) == 0.0 );
    // Inconsistent moments (M4 above the coherent bound) also clamp to 0.
    REQUIRE( m2m4_cn0_db_hz( 10.0, 1000.0, 1e-3 ) == 0.0 );
}

TEST_CASE( "carrier_lock_cos2phi_matches_phase", "[tracking][lock]" )
{
    // I = A*cos(phi), Q = A*sin(phi): NBD/NBP = (I^2 - Q^2)/(I^2 + Q^2) = cos(2*phi).
    const double A = 3.0;
    for( double deg : { 0.0, 30.0, 45.0, 60.0, 90.0 } )
    {
        const double phi = deg * M_PI / 180.0;
        const double i   = A * std::cos( phi );
        const double q   = A * std::sin( phi );
        const double nbd = i * i - q * q;
        const double nbp = i * i + q * q;
        INFO( "phi=" << deg << " deg" );
        REQUIRE( carrier_lock_cos2phi( nbd, nbp ) == Catch::Approx( std::cos( 2.0 * phi ) ).margin( 1e-12 ) );
    }
    REQUIRE( carrier_lock_cos2phi( 0.0, 0.0 ) == 0.0 ); // guarded: no power -> 0
}
#endif

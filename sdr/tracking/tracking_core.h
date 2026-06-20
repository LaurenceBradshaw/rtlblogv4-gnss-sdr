#pragma once
#include <array>
#include <vector>
#include "acquisition.h"
#include "signal.h"
#include "tracker.h"
#include "types.h"

// mirrors sdrtrkprm_t: coefficients for 2nd-order PLL + 1st-order FLL + 2nd-order DLL
struct Tracking_loop_prm
{
    double dllb, pllb, fllb;
    double dllw2, dllaw; // DLL coefficients
    double pllw2, pllaw; // PLL coefficients
    double fllw;         // FLL coefficient
};

// Shared carrier+code tracking core: carrier/code NCO, the plain-float N-tap correlate()
// hot loop, the loop-filter math (PLL/FLL/DLL), and the secondary-code helpers. Abstract:
// concrete strategies derive from it and implement run_loops():
//   - Costas_tracker : data-bearing BPSK (GPS L1 C/A, L2C, BeiDou B1I, ...) - Costas + FLL.
//   - Pilot_tracker  : pilot+data BOC (Galileo E1, GPS L1C, ...) - pure PLL, secondary, VEML.
// mirrors GNSS-SDRLIB sdrtrk_t + sdrtracking() + pll() + dll().
class Tracking_core : public Tracker
{
public:
    // Max correlator taps the shared correlate() supports. Costas uses 3 (P/E/L); the double-estimator
    // pilot tracker uses 7: P, code-E/L each on the in-phase AND quadrature subcarrier (subcarrier-phase-
    // independent envelope), plus subcarrier-E/L. n_taps_ selects how many are active.
    static constexpr int MAX_TAPS  = 8;
    static constexpr int LOOP_L1CA = 10; // loop filter interval (epochs)

    // Secondary (overlay) code sync acceptance: differential correlation peak / differential
    // energy must exceed this (see advance_secondary_sync). Below the old coherent 0.5 because
    // the differential product squares the noise.
    static constexpr double SECONDARY_SYNC_RATIO = 0.35;
    // Secondary-sync VERIFICATION (fixes the intermittent Galileo medium-C/N0 false sync). A correct
    // secondary phase lets the pure PLL clean the pilot so cos(2*phi) climbs above SECONDARY_SYNC_MIN_LOCK
    // within ~SECONDARY_SYNC_VERIFY_EPOCHS; a FALSE sync (picked at ~30-35 dB-Hz) leaves it stuck low. If
    // it hasn't cleaned up by then, drop the sync and re-search (a fresh window may catch the true phase).
    // 0.35 cleanly separates a false sync (cos(2*phi) stuck ~0.12) from a correct one (climbs to 0.6-0.99),
    // with margin so a correct-but-noisy lock's transient dips don't trigger a spurious re-search.
    static constexpr double SECONDARY_SYNC_MIN_LOCK      = 0.35; // "cleaned up" cos(2*phi) bar
    static constexpr int    SECONDARY_SYNC_VERIFY_EPOCHS = 1000; // grace after sync before declaring it false
    // Fast degradation check: if the carrier was ALREADY cleanly locked when the secondary synced
    // (cos(2*phi) > _WAS_CLEAN) but the sync then COLLAPSED it by more than _DEGRADE, that sync broke a
    // good lock - drop it quickly (don't wait the full grace). GPS L1C's clean Costas lock (~0.95) doesn't
    // need its long 1800-chip overlay, which false-syncs and scrambles the pilot; a correct Galileo CS25
    // sync, by contrast, leaves a clean lock clean. (For a DIRTY pre-sync lock the slow check above applies.)
    static constexpr double SECONDARY_SYNC_WAS_CLEAN   = 0.85;
    static constexpr double SECONDARY_SYNC_DEGRADE     = 0.40;
    static constexpr int    SECONDARY_SYNC_FAST_EPOCHS = 50;

    // Loop-filter noise bandwidths now live on each Signal (Signal_params::loop_bw_wide/narrow),
    // read by the Tracking_core ctor -> make_prm. Per-signal because the right values depend on the
    // integration period: the wide FLL that pulls Galileo's 4 ms epoch in is unstable on L1C's 10 ms.
    // (GPS L1 C/A and BeiDou B1I use the cross/dot-FLL set; Galileo E1 and GPS L1C the atan-FLL set -
    // selected by fll_active_, which now only picks the FLL discriminator, not the bandwidths.)
    // Extended coherent integration (pilot, once cleanly secondary-locked): coherently sum the
    // secondary-wiped correlators over EXTEND_SYMBOLS epochs, then run the loop once per window. The
    // ~10*log10(N) dB SNR boost steadies the medium-C/N0 lock that otherwise churns / false-syncs. Mirrors
    // gnss-sdr Galileo E1 (extend_correlation_symbols). We keep the prm2 bandwidth (not gnss-sdr's narrower
    // pll_bw_narrow): the SNR boost cleans the lock while the tighter loop keeps the per-epoch phase (hence
    // the cos(2*phi) lock test, and the PVT lock gate) above threshold - a narrow loop loosened it below.
    static constexpr int EXTEND_SYMBOLS = 2; // epochs per extended window (Galileo 4 ms -> 8 ms)

    // Lock detector. Two tests over a window of CN0_WINDOW_EPOCHS (both EMA-smoothed), combined with
    // hysteresis - the standard gnss-sdr scheme:
    //   (1) SNV (signal-to-noise variance) C/N0: Psig = (<|I|>)^2, noise = <I^2+Q^2> - Psig, SNR = Psig/noise.
    //       Scale-invariant, the SAME for every signal (only epoch_period_ feeds the dB-Hz scaling), and -
    //       unlike the old M2M4 - HIGH-SNR ROBUST: it measures noise as a variance instead of M2M4's
    //       M2 - sqrt(2 M2^2 - M4) difference of large near-equal terms (which saturated ~38 dB-Hz on long
    //       coherent integration, e.g. L1Cd's 10 ms code, making C/N0 incomparable across components).
    //   (2) Van Dierendonck carrier lock test cos(2*phi) ~ (<I^2> - <Q^2>) / (<I^2> + <Q^2>): ~+1 when
    //       all energy is on I (phase-locked), <=0 when the carrier is lost. (Power-difference form, not
    //       the coherent (sum I)^2 one, because an epoch is one code period - shorter than a nav bit -
    //       so coherent summing would cancel across bit flips; this form is data-robust.)
    // has_lock() requires BOTH above their thresholds, and only declares LOSS after LOCK_FAIL_WINDOWS
    // consecutive failing windows (so a momentary dip doesn't flap the indicator); it re-locks on the
    // first good window.
    // NOTE: this detects genuine CARRIER loss. It does NOT flag a false / cross-correlation lock - those
    // keep a deceptively HIGH C/N0 (~43 dB-Hz, measured) AND a coherent carrier (so cos(2*phi) also looks
    // healthy) yet never frame-sync; the frame-sync timeout (Channel) is what recovers those.
    static constexpr int    CN0_WINDOW_EPOCHS       = 100;  // per estimate (GPS 100 ms, Galileo 400 ms)
    static constexpr double CN0_LOCK_THRESHOLD_DBHZ = 28.0; // below this = lost lock / noise
    static constexpr double CARRIER_LOCK_THRESHOLD  = 0.7;  // cos(2*phi) gate (~phi < 23 deg); relaxed
                                                            // from 0.8 (~18 deg) - too strict for marginal
                                                            // low-C/N0 field signals (RTL-SDR), flapped the
                                                            // indicator without genuine carrier loss.
    static constexpr int LOCK_FAIL_WINDOWS = 20;            // consecutive bad windows before loss of lock

    // code_chips: raw +/-1 float chip values from Signal::code_chips()
    // sig: signal physics - chip rate (code NCO), code length, and RF carrier (used in
    //      the carrier-aided code-frequency term). Loop bandwidths above are still
    //      GPS-tuned; revisit per signal when adding non-GPS tracking.
    // secondary: overlay/secondary code chips (Signal::secondary_code()), empty for none.
    //      When present (Galileo E1-C CS25), tracking syncs to it and wipes it off so the
    //      FLL is valid. Empty -> every secondary code path is skipped (GPS byte-identical).
    // data_code: data-component code (pilot tracking), empty for none. When present, a prompt
    //      correlates it at the pilot's code phase to recover nav symbols (correlate_epoch
    //      then returns that data prompt instead of the pilot prompt).
    Tracking_core(
        const std::vector<float>& code_chips,
        double                    sample_rate_hz,
        const Signal_params&      sig,
        const std::vector<float>& secondary = {},
        const std::vector<float>& data_code = {}
    );

    // Initialise tracking state from acquisition result.
    // mirrors: sdr->trk.carrfreq = acq.acqfreq; sdr->trk.codefreq = sdr->crate
    void initialise( const Acquisition_result& acq ) override;

    // Re-center the carrier NCO on an external PHYSICAL Doppler (cross-code/-frequency aiding).
    void steer_carrier_doppler( double doppler_hz ) override;

    // Feed the common (PVT) clock-drift fraction forward into the carrier baseline - see Tracker. (Fix 3)
    void apply_common_clock_drift( double eps_fraction ) override;

    // Samples needed for the next epoch.
    // mirrors: samples_consumed = (clen - remcode) / (codefreq / f_sf)
    int compute_samples_needed() const override;

    // One tracking epoch, split so navigation can run between the two halves
    // (mirrors GNSS-SDRLIB's correlator -> sdrnavigation -> loop-filter order):
    //   1. correlate_epoch(block): correlate, return per-epoch prompts + samples consumed.
    //      block must hold >= compute_samples_needed() samples.
    //   2. caller runs navigation, which updates bit_sync / sw_loop for THIS epoch.
    //   3. run_loops(bit_sync, sw_loop): accumulate + run PLL/DLL using the current flags.
    Tracking_output correlate_epoch( const Sample_block& block ) override;

    // Latched lock state: set/cleared each completed window by update_lock_detectors() from the combined
    // C/N0 + carrier-lock-test criteria with consecutive-failure hysteresis (see the detector notes above).
    bool has_lock() const override
    {
        return locked_;
    }
    double get_cn0_db_hz() const override
    {
        return cn0_db_hz_;
    }

    double get_fractional_chip_time() const override;
    double code_phase_offset_s() const override;

    double get_code_freq() const override
    {
        return code_freq_;
    }

    double get_carrier_doppler_hz() const override
    { // TODO: Figure out what this is actually supposed to return from tracking.
        // return carrier_nco_;
        return carrier_freq_;
        // return carrier_freq_ - rf_freq_;
    }

    double get_carrier_acceleration() const override
    {
        return carrier_acc_;
    }

    double get_prompt_i() const override
    {
        return II_[0]; // raw prompt-I correlator output (data on I after the conj wipe)
    }
    double get_prompt_q() const override
    {
        return QQ_[0];
    }

    // run_loops() is the per-strategy step (Costas_tracker / Pilot_tracker); it stays pure
    // virtual here (inherited from Tracker), making Tracking_core abstract.

protected:
    // Per-tap replica policies for correlate_impl (compile-time, fully inlined - no per-sample branch or
    // virtual). Each forms tap k's reference sample at the running prompt code phase phase_P (code elements):
    // Bpsk reads the code array directly; the double-estimator factors the BOC replica into
    // primary[chip] x subcarrier(elem) with an independent subcarrier delay subc_off (see Pilot_tracker).
    struct Bpsk_replica
    {
        const float*  code;
        int           code_len;
        const double* tap_off;
        float         operator()( int k, double phase_P ) const
        {
            double ph = phase_P + tap_off[k];
            if( ph < 0.0 )
            {
                ph += code_len;
            }
            else if( ph >= code_len )
            {
                ph -= code_len;
            }
            return code[static_cast<int>( ph )];
        }
    };
    struct Boc_de_replica
    {
        const float*  code;
        int           code_len;
        const double* tap_off;
        const double* tap_sc;
        double        subc_off;
        float         operator()( int k, double phase_P ) const
        {
            double cph = phase_P + tap_off[k];
            if( cph < 0.0 )
            {
                cph += code_len;
            }
            else if( cph >= code_len )
            {
                cph -= code_len;
            }
            double sph = phase_P + subc_off + tap_sc[k];
            if( sph < 0.0 )
            {
                sph += code_len;
            }
            else if( sph >= code_len )
            {
                sph -= code_len;
            }
            const float prim = code[static_cast<int>( cph ) | 1];              // +primary[chip] = code_[2*chip+1]
            const float sc   = ( static_cast<int>( sph ) & 1 ) ? 1.0f : -1.0f; // subcarrier element parity
            return prim * sc;
        }
    };

    // The shared correlator hot loop, parameterised on a compile-time Replica policy (above) so the single
    // loop body serves both strategies with no per-sample branch. Each subclass implements correlate() as a
    // one-line call to this with its own policy (Costas -> Bpsk_replica; Pilot -> Boc_de_replica).
    template <class Replica>
    void         correlate_impl( const Sample_block& block, int n, Replica replica );
    virtual void correlate( const Sample_block& block, int n ) = 0; // per-strategy entry (-> correlate_impl)

    // mirrors GNSS-SDRLIB sdrtrk.c functions
    void cumsum_corr( int polarity );
    void clear_cumsum();
    // pure_pll: full 4-quadrant atan2(Q,I) PLL discriminator (data-free pilot after secondary
    // wipeoff) instead of the Costas (180 deg-folded) one; pair with use_fll=false.
    void pll_update( const Tracking_loop_prm& prm, double dt, bool use_fll, bool pure_pll );
    void dll_update( const Tracking_loop_prm& prm, double dt );

    // Lock detectors: accumulate this epoch's prompt I/Q and, once a window completes, update the
    // smoothed C/N0 + carrier lock test and the latched lock state. Called once per epoch from
    // correlate_epoch (shared by both trackers).
    void update_lock_detectors( double prompt_i, double prompt_q );

    // Per-strategy hooks, called by the shared correlate()/dll_update():
    //   configure_taps(ci): set n_taps_, tap_offset_chips_[] and (DE) tap_sc_offset_[] (Costas 3-tap E/P/L;
    //     Pilot 7-tap double-estimator). ci = code elements per sample (for sample-based E/L spacing).
    //   code_error():       normalised code discriminator (Costas E-L; Pilot subcarrier-independent envelope).
    virtual void   configure_taps( double ci ) = 0;
    virtual double code_error() const          = 0;

    // Secondary-code phase search: accumulate prompt-Q signs and, once enough history is
    // collected, correlate against the (known) secondary at every cyclic offset to find
    // its phase + polarity. Sets secondary_sync_ when a strong match is found.
    void advance_secondary_sync();

    static Tracking_loop_prm make_prm( double dllb, double pllb, double fllb );

    // Signal physics (from Signal_params) + hardware constants
    std::vector<float> code_;       // code elements (BPSK chips, or BOC half-chips), +/-1
    int                code_len_;   // number of code elements = code_.size()
    double             code_rate_;  // code-element rate (Hz): chip rate for BPSK, half-chip
                                    // rate for BOC (= code_len_ / code_period). Drives the
                                    // code NCO and the carrier-aiding scale.
    double epoch_period_;           // one code period (s) - loop-filter dt before bit sync
    bool   fll_active_;             // selects the FLL discriminator: true (>1 epoch per
                                    // data symbol, e.g. GPS) uses the cross/dot form;
                                    // false (1 symbol/epoch, e.g. Galileo E1) uses the
                                    // data-insensitive atan(I/Q) difference form.
    double rf_freq_;                // nominal RF carrier (Hz) - carrier-aided code term
    double code_elements_per_chip_; // code_rate_ / ranging chip rate (1 BPSK, 2 BOC)
    int    corr_spacing_;           // E/L spacing in samples
    double sample_rate_;
    double ti_; // 1/sample_rate (s)

    // Carrier/code NCO state - mirrors key sdrtrk_t fields. Defaulted here and (re)set by initialise();
    // the ctor only sets code_freq_ (= code_rate_), the rest start at rest.
    double acq_freq_       = 0.0; // initial acquisition Doppler (Hz), used as PLL baseline
    double code_freq_      = 0.0; // current code frequency (Hz) - mirrors trk.codefreq (set in ctor = code_rate_)
    double carrier_freq_   = 0.0; // current carrier frequency (Hz) - mirrors trk.carrfreq
    double remaining_code_ = 0.0; // remaining code phase (code elements) - mirrors trk.remcode
    double remaining_carr_ = 0.0; // remaining carrier phase (rad) - mirrors trk.remcarr
    double code_nco_       = 0.0; // code NCO accumulator
    double code_err_       = 0.0; // last code error (DLL)
    double carrier_nco_    = 0.0; // carrier NCO frequency offset (Hz) - main loop output
    double carrier_acc_    = 0.0; // 3rd-order PLL acceleration integrator (inner state)
    double carrier_err_    = 0.0; // last carrier error (PLL)
    double freq_err_       = 0.0; // last frequency error (FLL)

    // Common clock-drift feedforward (Fix 3): the amount of the receiver-wide PVT clock drift currently folded
    // into acq_freq_ (carrier Hz), and whether the first PVT push has seeded it (acquisition already baked the
    // recenter into acq_freq_ at hand-off, so the first push only records the baseline - no move - to avoid
    // double-counting). Both reset by initialise().
    double clock_ff_        = 0.0;
    bool   clock_ff_seeded_ = false;

    // Lock detector state (updated every epoch in correlate_epoch; see update_lock_detectors).
    double m2_sum_    = 0.0; // running sum of prompt power (I^2+Q^2) over the window (SNV total power + lock NBP)
    double abs_i_sum_ = 0.0; // running sum of |prompt I| over the window (SNV coherent signal amplitude)
    double nbd_sum_   = 0.0; // running sum of (I^2 - Q^2) over the window (carrier-lock numerator)
    int    cn0_count_ = 0;   // epochs accumulated into the current window
    double cn0_db_hz_ = 0.0; // smoothed C/N0 estimate (dB-Hz); 0 until the first window completes

    double carrier_lock_test_ = 0.0;   // smoothed cos(2*phi) estimate; ~+1 locked, <=0 no carrier
    bool   lock_seeded_       = false; // false until the first window seeds the smoothed estimates
    int    lock_fail_count_   = 0;     // consecutive completed windows failing the lock criteria
    bool   locked_            = false; // latched lock state returned by has_lock()

    // Correlator taps. tap_offset_chips_ holds each tap's code-phase offset (code elements); n_taps_ is how
    // many are active. The tap-index meaning is set by configure_taps(): Costas uses 0=Prompt, 1=Early, 2=Late;
    // the double-estimator pilot uses a 7-tap layout (see Pilot_tracker::configure_taps).
    int                          n_taps_ = 3; // default Costas P/E/L; configure_taps() resets it each epoch
    std::array<double, MAX_TAPS> tap_offset_chips_ {};

    // mirrors sdrtrk_t II/QQ/oldI/oldQ/sumI/sumQ/oldsumI/oldsumQ (per-epoch and cumulative correlator I/Q)
    std::array<double, MAX_TAPS> II_ {}, QQ_ {};
    std::array<double, MAX_TAPS> old_I_ {}, old_Q_ {};
    std::array<double, MAX_TAPS> sum_I_ {}, sum_Q_ {};
    std::array<double, MAX_TAPS> oldsum_I_ {}, oldsum_Q_ {};

    // Data-component prompt (pilot tracking): a single prompt correlating data_code_ at the
    // pilot's code phase, giving the nav symbol. Empty data_code_ -> not used.
    std::vector<float> data_code_;
    bool               has_data_;
    double             data_prompt_i_ = 0.0, data_prompt_q_ = 0.0;         // this epoch's data prompt
    double             old_data_prompt_i_ = 0.0, old_data_prompt_q_ = 0.0; // previous (for nav prev-symbol)

    int epoch_count_ = 0; // total epochs since initialise()

    // Secondary (overlay) code state - only used when secondary_ is non-empty (pilot).
    std::vector<float> secondary_;                  // overlay code chips (e.g. CS25), 1 per epoch
    bool               secondary_sync_     = false; // phase acquired -> wipe-off active, FLL valid
    int                secondary_index_    = 0;     // position in secondary_ for the next epoch
    int                secondary_polarity_ = 1;     // +/-1, resolved at sync (carrier sign)
    int                epoch_at_sync_      = 0;     // epoch_count_ when the secondary last synced (sync verify)
    double             cos2phi_at_sync_    = 0.0;   // carrier lock at sync time (fast degradation check)
    std::vector<float> sec_i_hist_;                 // recent prompt-I values (complex pair with sec_q_hist_)
    std::vector<float> sec_q_hist_;

    Tracking_loop_prm prm1_; // loop params before nav frame sync
    Tracking_loop_prm prm2_; // loop params after  nav frame sync
};

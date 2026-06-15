#pragma once
#include <memory>
#include <utility>
#include <vector>
#include "constellations.h"
#include "nav_decoder.h"

// GNSS signal abstraction
//
// Everything constellation-specific in the receiver is reached through a Signal:
// its physical parameters, how to generate its spreading code for a given SV, the
// range of SV ids it defines, and a factory for its navigation-message decoder.
//
// The acquisition/tracking engines and the Channel are constructed FROM a Signal, so
// adding a signal means writing a new Signal subclass + Nav_decoder - it does not touch
// the DSP engines, the scheduler, or the channel state machine.
//
// Implemented: GPS L1 C/A (LNAV) and Galileo E1 (I/NAV, via Pilot_tracker/VEML); both
// decode end to end. BeiDou B1I has code generation; its nav decoder is a skeleton.

class Tracker; // built per (signal, SV) by Signal::make_tracker(); see tracker.h

// A signal is identified by a (constellation, band, code) triple. The constellation is the
// system; the band is the RF carrier (L1/E1/B1 all sit near 1.5 GHz but are distinct); the
// code is the component within that band that carries (or anchors) the service. make_signal()
// validates the triple - only the combinations it knows are buildable.

// RF band / carrier. Named by the owning system's convention; values across systems that
// happen to share a frequency (L1 == E1 == 1575.42 MHz) are still kept distinct.
enum class Band
{
    L1,
    L2,
    L5,
    E1,
    B1
};

// Code / signal component within a band. For data+pilot services the data component (the one
// carrying the nav message) names the signal; its paired pilot is generated alongside it.
//   CA - GPS L1 C/A (coarse/acquisition)
//   B  - Galileo E1-B (I/NAV data; paired with the E1-C pilot)
//   I  - BeiDou B1I (in-phase)
enum class Code
{
    CA, // GPS L1 C/A
    B,  // Galileo E1-B (I/NAV data; paired with E1-C pilot)
    I,  // BeiDou B1I
    Cd  // GPS L1Cd (CNAV-2 data; paired with the L1Cp pilot)
};

// Spreading-symbol modulation. Determines how the code generator lays out chips and
// (later) how the acquisition/tracking replica is built:
//   Bpsk   - one chip per code chip (GPS L1 C/A, BeiDou B1I).
//   Boc11  - BOC(1,1): each chip split into two opposite-sign half-chips (the code
//            generator bakes this in, doubling the effective chip count / rate). Used
//            for Galileo E1, as the BOC(1,1) approximation of CBOC(6,1,11).
enum class Modulation
{
    Bpsk,
    Boc11
};

// Tracking-loop noise bandwidths (Hz) used to derive the filter coefficients
// (Tracking_core::make_prm). Per signal because the right values depend on the
// integration period (code_period_s) and C/N0: e.g. the wide FLL that pulls Galileo's
// 4 ms epoch in is unstable (B*T too high) on GPS L1C's 10 ms epoch.
struct Loop_bw
{
    double dll; // code (DLL) noise bandwidth
    double pll; // carrier phase (PLL) noise bandwidth
    double fll; // carrier frequency (FLL) noise bandwidth
};

// Fundamental, PRN-independent physics of a signal. Drives acquisition/tracking
// sizing and the carrier/code NCOs.
struct Signal_params
{
    Constellation constellation;
    Band          band;              // RF band (selects the carrier)
    Code          code;              // signal component within the band
    const char*   name;              // human-readable, e.g. "GPS L1 C/A"
    double        carrier_freq_hz;   // nominal RF centre (e.g. 1575.42e6)
    double        chip_rate_hz;      // primary-code chip rate (chip/s)
    int           code_length_chips; // primary-code length (chips)
    double        code_period_s;     // primary-code period (s) - one tracking epoch
    int           nav_bit_ms;        // tracking epochs (ms) per nav data symbol
    Modulation    modulation;        // chip modulation (BPSK vs BOC(1,1))

    // Acquisition tuning. Defaults (GPS): 10 non-coherent epochs, 2x zero-padded FFT.
    //   acq_integrations: non-coherent epochs accumulated before the accept/reject decision.
    //     Signals with a longer coherent code period (Galileo E1 = 4 ms vs GPS 1 ms) get more
    //     gain per epoch, so they need fewer epochs.
    //   acq_fft_factor: acquisition FFT size = factor * samples-per-code-period. 2 = the
    //     classic zero-padded linear correlation; 1 = single-period circular correlation
    //     (half the FFT size and half the Doppler bins, ~3.9 dB worst-case Doppler scalloping
    //     - fine for a strong, long (4 ms) code where each epoch has ample margin).
    int acq_integrations;
    int acq_fft_factor;

    // Carrier/code loop bandwidths. "wide" is the pull-in set (prm1, before nav frame sync);
    // "narrow" is the tighter steady-state set (prm2, after sync). See Loop_bw.
    Loop_bw loop_bw_wide;
    Loop_bw loop_bw_narrow;

    // Frame-sync timeout (s): if a channel tracks this long without achieving nav frame sync
    // (Nav_decoder::frame_synced), it locked onto a false / cross-correlation peak and is dropped
    // back to re-acquire (see Channel::process_tracking). Per-signal because frame sync takes very
    // different times: GPS subframe preamble ~13 s; Galileo I/NAV page-vote can be 30 s+.
    double frame_sync_timeout_s;
};

// A GNSS signal definition. One concrete Signal per constellation/band; it is shared
// (by const reference) across all the Channels tracking that signal.
class Signal
{
public:
    virtual ~Signal() = default;

    virtual const Signal_params& params() const = 0;

    // Inclusive [min, max] SV/PRN ids this signal defines codes for.
    virtual std::pair<int, int> sv_range() const = 0;

    // One sampled code period (Q = 0) at sample_rate_hz - the acquisition replica.
    virtual Complex_buf code_samples( Satellite_id sv, double sample_rate_hz ) const = 0;

    // Raw +/-1 chip values (code_length_chips of them) - the tracking replica.
    virtual std::vector<float> code_chips( Satellite_id sv ) const = 0;

    // A fresh navigation decoder matching this signal's message format.
    virtual std::unique_ptr<Nav_decoder> make_nav_decoder( Satellite_id sv ) const = 0;

    // A fresh carrier+code tracker for one SV. Default builds a Costas_tracker (data/BPSK);
    // signals that need a different strategy (e.g. Galileo E1 -> Pilot_tracker) override it.
    virtual std::unique_ptr<Tracker> make_tracker( Satellite_id sv, double sample_rate_hz ) const;

    // Secondary (overlay) code: one chip per primary-code period, deterministic and
    // known. Tracking syncs to it then wipes it off, which lets the FLL work and
    // coherent integration extend across epochs. Empty for signals without one (GPS L1 C/A,
    // BeiDou B1I, Galileo E1-B data). Per-SV: Galileo E1-C is the same CS25 for all PRNs, but
    // GPS L1Cp's overlay differs per PRN - hence the sv argument (the default ignores it).
    virtual std::vector<float> secondary_code( Satellite_id /*sv*/ ) const
    {
        return {};
    }

    // Data-component code, for pilot-tracking signals: when the carrier/code loops track
    // the (data-free) pilot, a separate prompt correlates this code at the pilot's code
    // phase to recover the nav symbols. Empty for signals tracked directly (the prompt is
    // already the data); Galileo E1-C pilot returns the E1-B data code.
    virtual std::vector<float> data_code_chips( Satellite_id sv ) const
    {
        return {};
    }
};

// Construct the Signal for a (constellation, band, code) triple. The single registry of known
// signals; throws std::invalid_argument for any combination that is not implemented.
std::unique_ptr<Signal> make_signal( Constellation constellation, Band band, Code code );

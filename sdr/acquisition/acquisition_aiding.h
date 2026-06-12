#pragma once
#include <mutex>

// Receiver-wide acquisition aiding: a shared estimate of the common-mode carrier-frequency
// offset - the receiver clock / local-oscillator error, identical for every satellite on a
// given band - used to RECENTER each channel's Doppler search so the window tracks the true
// (clock-shifted) satellite cluster instead of sitting on nominal 0 Hz. One instance is shared
// by every Channel (the offset is a property of the front-end, not of any one satellite).
//
// Why it helps: f_obs(SV) = f_clock (common) + f_doppler(SV) (per-SV, +/-~5 kHz). A cheap SDR
// front-end can have f_clock of several kHz (this capture ~ -7.5 kHz), which alone pushes the
// edge satellites outside a 0-centred +/-ACQHBAND window. Recentering on f_clock recovers them.
//
// Two sources of the estimate, lowest priority first:
//   * PRE-PVT bridge (this file): the mean observed Doppler of the satellites acquired so far.
//     Averaging across SVs cancels the per-SV Doppler, leaving f_clock. With a SINGLE SV we
//     move only halfway, as a damage-limiting hedge - that one SV's own Doppler is an unknown
//     +/-5 kHz error, and over-shifting could push a not-yet-found SV out of the window. We
//     only RECENTER here; the search width is NOT narrowed (one/few SVs can't justify it).
//   * POST-PVT: once PVT solves the receiver clock DRIFT, set_clock_fraction() is called with
//     that rigorous common-mode offset. It then overrides the bridge - see from_pvt.
//
// As well as the recenter, the aiding chooses the search HALF-WIDTH (ACQHBAND): WIDE while
// bootstrapping (0 or 1 SV - the recenter is not yet trustworthy, so keep a broad net), then
// NARROW once >=2 SVs pin the common offset. Narrowing shrinks the Doppler bin count, which is
// the dominant acquisition cost (fewer IFFTs per attempt). The engine allocates for WIDE and
// searches the active half-width each attempt.
//
// The offset is stored as a FRACTION of carrier (LO error is ~proportional to carrier), so it
// is band-agnostic: report and query with the channel's own carrier frequency.
class Acquisition_aiding
{
public:
    // Doppler search half-widths (Hz). Even after recentering on the common-mode clock offset,
    // each satellite's OWN Doppler still spreads +/-~5 kHz around that center (LOS velocity up to
    // ~+/-800 m/s near the horizon), so no recentered window can shrink below that ~5 kHz floor
    // without permanently missing high-Doppler (low/rising) satellites:
    //   WIDE        - uncentered (0/1 SV): covers the clock offset itself, sized for the grid.
    //   NARROW      - >=2 SV bridge center (known to a couple kHz): ~5 kHz spread + that margin.
    //   VERY_NARROW - PVT center (rigorous): just the ~5 kHz per-SV Doppler spread.
    // (A genuinely small window needs per-SV Doppler prediction from the almanac, not just the
    // common clock - that is future work.)
    static constexpr double SEARCH_HBAND_WIDE_HZ        = 8000.0;
    static constexpr double SEARCH_HBAND_NARROW_HZ      = 6000.0;
    static constexpr double SEARCH_HBAND_VERY_NARROW_HZ = 5000.0;

    struct Estimate
    {
        double center_hz;     // Doppler recenter target for the queried carrier (0 = no aiding yet)
        double half_width_hz; // search half-width to use this attempt (WIDE or NARROW)
        int    n;             // satellites contributing to the bridge estimate
        bool   from_pvt;      // true once a PVT clock-drift estimate is set
    };

    // Report an acquired satellite's observed Doppler at its carrier. Thread-safe (channels
    // acquire concurrently on the pool); cheap and rare (once per successful acquisition).
    void report( double doppler_hz, double carrier_hz );

    // Recommended recenter + search half-width for a search on carrier_hz. Thread-safe.
    Estimate estimate( double carrier_hz ) const;

    // PVT hook: called with the receiver clock drift expressed as a fraction of carrier (df/f).
    // Once set it takes precedence over the satellite-mean bridge.
    void set_clock_fraction( double fraction );

    // Clear all aiding state back to the un-aided start (for a fresh run / restart).
    void reset()
    {
        std::lock_guard<std::mutex> lk( mu_ );
        sum_fraction_         = 0.0;
        n_                    = 0;
        clock_fraction_       = 0.0;
        clock_fraction_valid_ = false;
    }

private:
    // The recenter fraction the bridge recommends right now (caller MUST hold mu_). The single
    // source of truth shared by estimate() and the report() logs so they can never disagree:
    // >=2 SVs use the full mean Doppler; a lone SV uses HALF of it (hedge against that one SV's
    // unknown own-Doppler). Returns 0 when no SV has been reported yet.
    double bridge_fraction_locked() const;

    mutable std::mutex mu_;
    double             sum_fraction_ = 0.0; // sum of doppler/carrier over acquired SVs (bridge)
    int                n_            = 0;

    // Rigorous common-mode offset from the PVT clock-drift solution. Unset until PVT.
    double clock_fraction_       = 0.0;
    bool   clock_fraction_valid_ = false;
};

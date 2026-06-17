#pragma once
#include <map>
#include "almanac.h"
#include "ephemeris.h"
#include "ionospheric.h"

// Per-epoch navigation-message decoder. One concrete decoder per signal type; the
// Channel owns one polymorphically and feeds it the prompt-I every tracking epoch.
class Nav_decoder
{
public:
    virtual ~Nav_decoder() = default;

    // Called every tracking epoch immediately after correlation.
    //   prompt_i      : this epoch's raw prompt-I
    //   prompt_i_prev : previous epoch's raw prompt-I (for bit-sync sign changes)
    virtual void process( double prompt_i, double prompt_i_prev ) = 0;

    virtual bool bit_sync_found() const = 0; // data-bit boundary acquired
    virtual bool sw_loop() const        = 0; // fire the (slow) loop filter this epoch

    // True once the decoder has achieved FRAME sync (GPS subframe preamble / Galileo I/NAV page sync) -
    // a stronger condition than bit_sync_found(). The Channel uses this to spot a channel that locked
    // onto a false / cross-correlation peak: it still gets bit sync on the garbage but never frame sync,
    // so a timeout drops it back to re-acquire. Defaults to bit-sync for decoders that don't track it.
    virtual bool frame_synced() const
    {
        return bit_sync_found();
    }

    virtual const Ephemeris& ephemeris() const = 0; // current ephemeris (valid or not); never null
    virtual const Iono*      iono() const      = 0; // current broadcast iono (valid or not); never null

    // Almanac entries decoded so far (PRN -> coarse orbit), subcommutated over the broadcast. Empty for
    // decoders that don't decode the almanac (default); GPS L1CA overrides. Constellation-wide, so the
    // Receiver aggregates these across channels for acquisition aiding.
    virtual const std::map<int, Almanac>& almanac() const
    {
        static const std::map<int, Almanac> empty;
        return empty;
    }

    virtual int          get_current_bit_index() const = 0;
    virtual int          get_current_ms_tick() const   = 0;
    virtual const double ms_since_tow_update() const   = 0;

    const bool tow_anchored() const
    {
        return tow_confirmed_;
    }

protected:
    // TOW-continuity anchor (see tow_is_continuous): the last decoded TOW and the epoch it
    // arrived. A subframe is trusted (and applied) only if its TOW agrees, so two consecutive
    // consistent subframes are needed before any reaches the ephemeris - keeps garbage out of PVT.
    bool     tow_anchored_     = false; // a (tentative) anchor exists
    bool     tow_confirmed_    = false; // a second TOW agreed -> anchor is trustworthy
    double   tow_anchor_value_ = 0.0;   // last good TOW (s)
    uint64_t tow_anchor_epoch_ = 0;     // epoch_count_ at that TOW
    uint64_t active_eph_tow_anchor_epoch_ =
        0; // epoch_count_ at that TOW for the current ephemeris (used to gate TOW continuity on the same subframe, so a new
           // subframe can re-anchor if the old one was bad)
};
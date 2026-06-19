#include "acquisition_aiding.h"
#include "logging.h"

void Acquisition_aiding::report( double doppler_hz, double carrier_hz )
{
    std::lock_guard<std::mutex> lk( mu_ );
    sum_fraction_ += doppler_hz / carrier_hz;
    ++n_;

    if( n_ == 1 && !clock_fraction_valid_ )
    {
        // The single Doppler is sum_fraction_*carrier; the recenter only moves HALFWAY to it.
        logging::log(
            logging::Level::Info,
            fmt::format(
                "Acquisition aiding - 1 satellite acquired (Doppler {:.1f} Hz), recenter {:.1f} Hz (halfway hedge)",
                sum_fraction_ * carrier_hz,
                bridge_fraction_locked() * carrier_hz
            )
        );
    }
    if( n_ == 2 && !clock_fraction_valid_ )
    {
        logging::log(
            logging::Level::Info,
            fmt::format(
                "Acquisition aiding - {} satellites acquired, bridge estimate {:.1f} Hz (mean Doppler), search "
                "half-width narrowed to {:.0f} Hz",
                n_,
                bridge_fraction_locked() * carrier_hz,
                SEARCH_HBAND_NARROW_HZ
            )
        );
    }
}

Acquisition_aiding::Estimate Acquisition_aiding::estimate_common_locked( double carrier_hz ) const
{
    if( clock_fraction_valid_ )
    {
        // Post-PVT: rigorous common-mode offset, trusted enough to narrow the search.
        return { clock_fraction_ * carrier_hz, SEARCH_HBAND_VERY_NARROW_HZ, n_, true, true };
    }
    if( n_ == 0 )
    {
        return { 0.0, SEARCH_HBAND_WIDE_HZ, 0, false, true };
    }
    // Narrow only once >=2 SVs pin the common offset; stay wide while bootstrapping.
    const double half_width = ( n_ >= 2 ) ? SEARCH_HBAND_NARROW_HZ : SEARCH_HBAND_WIDE_HZ;
    return { bridge_fraction_locked() * carrier_hz, half_width, n_, false, true };
}

Acquisition_aiding::Estimate Acquisition_aiding::estimate( double carrier_hz ) const
{
    std::lock_guard<std::mutex> lk( mu_ );
    return estimate_common_locked( carrier_hz );
}

Acquisition_aiding::Estimate Acquisition_aiding::estimate( Constellation con, int prn, double carrier_hz ) const
{
    std::lock_guard<std::mutex> lk( mu_ );
    const int key = sv_key( con, prn );

    // (1) MEASURED sibling Doppler (highest precedence): another channel has this SV cleanly locked, so its
    // full observed Doppler/carrier is the best center - scale by this carrier. A live sibling lock proves
    // the SV is up, so it is always searchable here.
    const auto md = sv_doppler_fraction_.find( key );
    if( md != sv_doppler_fraction_.end() )
    {
        return { md->second * carrier_hz, SEARCH_HBAND_ALMANAC_HZ, n_, clock_fraction_valid_, true };
    }

    const auto it = predictions_.find( key );
    if( it == predictions_.end() )
    {
        return estimate_common_locked( carrier_hz ); // no per-SV measurement or prediction yet
    }
    // Per-SV almanac prediction: center on this SV's predicted LOS Doppler PLUS the common clock offset
    // (the prediction is geometry-only). With both the almanac geometry and a clock pin, the ALMANAC
    // half-width (just the prediction error) collapses the search to a couple of bins. searchable reflects
    // the horizon prediction.
    const double clock  = clock_fraction_valid_ ? clock_fraction_ : bridge_fraction_locked();
    const double center = ( it->second.los_fraction + clock ) * carrier_hz;
    return { center, SEARCH_HBAND_ALMANAC_HZ, n_, clock_fraction_valid_, it->second.above_horizon };
}

void Acquisition_aiding::set_prediction( Constellation con, int prn, bool above_horizon, double los_doppler_fraction )
{
    std::lock_guard<std::mutex> lk( mu_ );
    predictions_[sv_key( con, prn )] = { above_horizon, los_doppler_fraction };
}

void Acquisition_aiding::report_sv_doppler( Constellation con, int prn, double doppler_hz, double carrier_hz )
{
    std::lock_guard<std::mutex> lk( mu_ );
    sv_doppler_fraction_[sv_key( con, prn )] = doppler_hz / carrier_hz; // band-agnostic df/f
}

bool Acquisition_aiding::sv_doppler_fraction( Constellation con, int prn, double& fraction_out ) const
{
    std::lock_guard<std::mutex> lk( mu_ );
    const auto it = sv_doppler_fraction_.find( sv_key( con, prn ) );
    if( it == sv_doppler_fraction_.end() )
    {
        return false;
    }
    fraction_out = it->second;
    return true;
}

bool Acquisition_aiding::searchable( Constellation con, int prn ) const
{
    std::lock_guard<std::mutex> lk( mu_ );
    const auto it = predictions_.find( sv_key( con, prn ) );
    return it == predictions_.end() || it->second.above_horizon; // unknown SV stays searchable
}

void Acquisition_aiding::set_clock_fraction( double fraction )
{
    if( !clock_fraction_valid_ )
    {
        logging::log(
            logging::Level::Info,
            fmt::format(
                "Acquisition aiding - PVT clock-drift solution set, recenter {:.3f} ppm, search half-width narrowed to "
                "{:.0f} Hz",
                fraction * 1e6, // report in ppm for readability
                SEARCH_HBAND_VERY_NARROW_HZ
            )
        );
    }
    std::lock_guard<std::mutex> lk( mu_ );
    clock_fraction_       = fraction;
    clock_fraction_valid_ = true;
}

double Acquisition_aiding::bridge_fraction_locked() const
{
    if( n_ == 0 )
    {
        return 0.0;
    }
    const double mean = sum_fraction_ / n_;
    return ( n_ == 1 ) ? 0.5 * mean : mean;
}

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

Acquisition_aiding::Estimate Acquisition_aiding::estimate( double carrier_hz ) const
{
    std::lock_guard<std::mutex> lk( mu_ );
    if( clock_fraction_valid_ )
    {
        // Post-PVT: rigorous common-mode offset, trusted enough to narrow the search.
        return { clock_fraction_ * carrier_hz, SEARCH_HBAND_VERY_NARROW_HZ, n_, true };
    }
    if( n_ == 0 )
    {
        return { 0.0, SEARCH_HBAND_WIDE_HZ, 0, false };
    }
    // Narrow only once >=2 SVs pin the common offset; stay wide while bootstrapping.
    const double half_width = ( n_ >= 2 ) ? SEARCH_HBAND_NARROW_HZ : SEARCH_HBAND_WIDE_HZ;
    return { bridge_fraction_locked() * carrier_hz, half_width, n_, false };
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

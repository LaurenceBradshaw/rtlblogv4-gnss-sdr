#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace utils
{
template <typename T>
double mean( const std::vector<T>& v )
{
    if( v.empty() )
    {
        return 0.0;
    }

    double sum = std::accumulate( v.begin(), v.end(), 0.0 );
    return sum / static_cast<double>( v.size() );
}

template <typename T>
double variance( const std::vector<T>& v, bool sample = true )
{
    if( v.size() < 2 )
    {
        return 0.0;
    }

    double m     = mean( v );
    double accum = 0.0;

    for( auto val : v )
    {
        accum += ( val - m ) * ( val - m );
    }
    return accum / static_cast<double>( v.size() - ( sample ? 1 : 0 ) );
}

template <typename T>
double standard_deviation( const std::vector<T>& v, bool sample = true )
{
    return std::sqrt( variance( v, sample ) );
}

template <typename T>
T median( std::vector<T> v )
{
    if( v.empty() )
    {
        throw std::invalid_argument( "median of empty vector" );
    }

    std::sort( v.begin(), v.end() );
    size_t n = v.size();
    if( n % 2 == 0 )
    {
        return static_cast<T>( ( v[n / 2 - 1] + v[n / 2] ) / 2.0 );
    }
    else
    {
        return v[n / 2];
    }
}

template <typename T>
int median_index( const std::vector<T>& v )
{
    if( v.empty() )
    {
        throw std::invalid_argument( "median_index of empty vector" );
    }
    // For even, pick the lower index (n/2 - 1)
    return ( v.size() % 2 == 0 ) ? static_cast<int>( v.size() / 2 - 1 ) : static_cast<int>( v.size() / 2 );
}

template <typename T>
T mode( const std::vector<T>& v )
{
    if( v.empty() )
    {
        throw std::invalid_argument( "mode of empty vector" );
    }

    std::unordered_map<T, int> freq;
    for( const auto& val : v )
    {
        ++freq[val];
    }

    auto max_it =
        std::max_element( freq.begin(), freq.end(), []( const auto& a, const auto& b ) { return a.second < b.second; } );

    return max_it->first;
}

// Normalise vector to unit length
template <typename T>
std::vector<double> normalise_vector( const std::vector<T>& v )
{
    double norm = std::sqrt( std::inner_product( v.begin(), v.end(), v.begin(), 0.0 ) );
    if( norm == 0.0 )
    {
        return std::vector<double>( v.size(), 0.0 );
    }

    std::vector<double> result( v.size() );
    std::transform( v.begin(), v.end(), result.begin(), [norm]( T val ) { return static_cast<double>( val ) / norm; } );

    return result;
}

// z-score normalisation (centre and scale)
template <typename T>
std::vector<double> zscore_normalise( const std::vector<T>& v )
{
    if( v.empty() )
    {
        return {};
    }

    double m  = mean( v );
    double sd = standard_deviation( v );
    if( sd == 0.0 )
    {
        return std::vector<double>( v.size(), 0.0 );
    }

    std::vector<double> result( v.size() );
    for( size_t i = 0; i < v.size(); ++i )
    {
        result[i] = ( v[i] - m ) / sd;
    }

    return result;
}

inline double deg_to_rad( double degrees )
{
    return degrees * M_PI / 180.0;
}

inline double rad_to_deg( double radians )
{
    return radians * 180.0 / M_PI;
}

} // namespace utils
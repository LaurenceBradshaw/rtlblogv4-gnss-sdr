#pragma once
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <vector>

namespace utils
{
template <typename T>
constexpr bool contains( const std::vector<T>& v, const T& val )
{
    return std::find( std::begin( v ), std::end( v ), val ) != std::end( v );
}

// Return the index of the element. -1 in the case it doesn't exist
template <typename T>
constexpr int find_index( const std::vector<T>& v, const T& val )
{
    auto it = std::find( std::begin( v ), std::end( v ), val );
    return ( it != std::end( v ) ) ? static_cast<int>( std::distance( std::begin( v ), it ) ) : -1;
}

// Return index of first element >= value
template <typename T, typename U>
int lower_bound_index( const std::vector<T>& v, const U& val )
{
    if( v.empty() )
    {
        return -1;
    }
    T    val_t = static_cast<T>( val );
    auto it    = std::lower_bound( v.begin(), v.end(), val_t );
    return ( it != v.end() ) ? static_cast<int>( std::distance( v.begin(), it ) ) : -1;
}

// Return index of first element > value
template <typename T, typename U>
int upper_bound_index( const std::vector<T>& v, const U& val )
{
    if( v.empty() )
    {
        return -1;
    }
    T    val_t = static_cast<T>( val );
    auto it    = std::upper_bound( v.begin(), v.end(), val_t );
    return ( it != v.end() ) ? static_cast<int>( std::distance( v.begin(), it ) ) : -1;
}

// Check if a value exists in sorted vector
template <typename T, typename U>
bool binary_search_contains( const std::vector<T>& v, const U& val )
{
    if( v.empty() )
    {
        return false;
    }
    T val_t = static_cast<T>( val );
    return std::binary_search( v.begin(), v.end(), val_t );
}

// Return index of element closest to given value
template <typename T, typename U>
int nearest_index( const std::vector<T>& v, const U& val )
{
    if( v.empty() )
    {
        return -1;
    }

    using common_t  = std::common_type_t<T, U>;
    common_t target = static_cast<common_t>( val );

    // do lower_bound on T (cheap) using cast of target to T
    T    val_t = static_cast<T>( target );
    auto it    = std::lower_bound( v.begin(), v.end(), val_t );

    if( it == v.begin() )
        return 0;
    if( it == v.end() )
        return static_cast<int>( v.size() - 1 );

    auto     prev_it = std::prev( it );
    common_t d_curr  = std::abs( static_cast<common_t>( *it ) - target );
    common_t d_prev  = std::abs( static_cast<common_t>( *prev_it ) - target );

    return ( d_curr < d_prev ) ? static_cast<int>( std::distance( v.begin(), it ) )
                               : static_cast<int>( std::distance( v.begin(), prev_it ) );
}

template <typename T>
void remove_value( std::vector<T>& v, const T& val )
{
    v.erase( std::remove( v.begin(), v.end(), val ), v.end() );
}

template <typename T, typename MemberType>
std::vector<MemberType> extract_row( const std::vector<T>& data, MemberType T::* member_ptr )
{
    std::vector<MemberType> vec( data.size() );
    for( size_t i = 0; i < data.size(); ++i )
    {
        vec[i] = data[i].*member_ptr;
    }
    return vec;
}

template <typename T>
std::vector<T> clamp( const std::vector<T>& v, const T& low, const T& high )
{
    std::vector<T> result( v.size() );
    std::transform(
        v.begin(),
        v.end(),
        result.begin(),
        [&low, &high]( const T& val ) { return ( val < low ) ? low : ( ( val > high ) ? high : val ); }
    );
    return result;
}

} // namespace utils
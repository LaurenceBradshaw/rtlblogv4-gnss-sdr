#pragma once
#include <Eigen/Dense>
#include <type_traits>

namespace utils
{

template <typename T, typename MemberType>
Eigen::Matrix<typename MemberType::Scalar, Eigen::Dynamic, MemberType::RowsAtCompileTime>
extract_matrix( const std::vector<T>& data, MemberType T::* member_ptr )
{
    constexpr int dim = MemberType::RowsAtCompileTime;

    Eigen::Matrix<typename MemberType::Scalar, Eigen::Dynamic, dim> mat( data.size(), dim );
    for( size_t i = 0; i < data.size(); ++i )
    {
        mat.row( i ) = ( data[i].*member_ptr ).transpose();
    }
    return mat;
}

template <typename T, typename MemberType>
Eigen::Matrix<MemberType, Eigen::Dynamic, 1> extract_vector( const std::vector<T>& data, MemberType T::* member_ptr )
{
    static_assert( std::is_arithmetic<MemberType>::value, "Member must be scalar" );
    Eigen::Matrix<MemberType, Eigen::Dynamic, 1> vec( data.size() );
    for( size_t i = 0; i < data.size(); ++i )
    {
        vec( i ) = data[i].*member_ptr;
    }
    return vec;
}

inline int median_index( const Eigen::VectorXd& v )
{
    if( v.size() == 0 )
    {
        throw std::invalid_argument( "median_index of empty vector" );
    }

    // store value + index pairs
    std::vector<std::pair<double, int>> val_idx;
    val_idx.reserve( v.size() );
    for( int i = 0; i < v.size(); ++i )
    {
        val_idx.emplace_back( v[i], i );
    }

    // sort by value
    std::sort( val_idx.begin(), val_idx.end(), []( const auto& a, const auto& b ) { return a.first < b.first; } );

    int n = static_cast<int>( v.size() );
    if( n % 2 == 0 )
    {
        // For even size, pick the lower middle
        return val_idx[n / 2 - 1].second;
    }
    else
    {
        return val_idx[n / 2].second;
    }
}

// Eigen clamp (element-wise)
template <typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>
clamp( const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& low, const typename Derived::Scalar& high )
{
    return mat.cwiseMax( low ).cwiseMin( high );
}

} // namespace utils
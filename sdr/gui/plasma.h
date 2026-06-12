#pragma once
#include <algorithm>
#include <array>
#include <QColor>

namespace gui_format
{
// matplotlib "plasma" colormap (dark purple -> magenta -> orange -> yellow), approximated by linear
// interpolation between control points. t in [0,1]. Used for the I/Q heatmap.
inline QColor plasma( double t )
{
    static constexpr int                          N  = 9;
    static constexpr std::array<std::array<int, 3>, N> CP = { {
        { 13, 8, 135 },
        { 75, 3, 161 },
        { 125, 3, 168 },
        { 168, 34, 150 },
        { 203, 70, 121 },
        { 229, 107, 93 },
        { 248, 148, 65 },
        { 253, 195, 40 },
        { 240, 249, 33 },
    } };

    t              = std::clamp( t, 0.0, 1.0 );
    const double x = t * ( N - 1 );
    const int    i = std::min( static_cast<int>( x ), N - 2 );
    const double f = x - i;
    const auto&  a = CP[i];
    const auto&  b = CP[i + 1];
    return QColor(
        static_cast<int>( a[0] + f * ( b[0] - a[0] ) ),
        static_cast<int>( a[1] + f * ( b[1] - a[1] ) ),
        static_cast<int>( a[2] + f * ( b[2] - a[2] ) )
    );
}
} // namespace gui_format

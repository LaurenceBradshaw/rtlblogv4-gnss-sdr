#pragma once
#include <cstdint>

// Bit-field extractors shared by the nav-message decoders, grouped by how the frame is stored in memory
// (both MSB-first, 0-based position):
//   PACKED   - 8 bits per byte (uint8_t buffer), e.g. GPS L1 C/A subframes.    packed_uint / packed_int
//              (+ packed_uint2 / packed_int2: concatenate two split fields, for LNAV).
//   UNPACKED - one bit per element (uint8_t or int), e.g. the Viterbi / LDPC   unpacked_uint / unpacked_int
//              decoder output of GPS L1C and Galileo I/NAV.
// Names omit a "bits" stem since they are always called qualified (bits::...).
namespace bits
{

// --- packed: 8 bits / byte -------------------------------------------------------------------------
inline uint32_t packed_uint( const uint8_t* buf, int pos, int len )
{
    uint32_t v = 0;
    for( int i = pos; i < pos + len; i++ )
    {
        v = ( v << 1 ) | ( ( buf[i / 8] >> ( 7 - i % 8 ) ) & 1u );
    }
    return v;
}

inline int32_t packed_int( const uint8_t* buf, int pos, int len ) // sign-extends the result
{
    const uint32_t u = packed_uint( buf, pos, len );
    if( len > 0 && ( u >> ( len - 1 ) ) )
    {
        return static_cast<int32_t>( u | ( ~0u << len ) );
    }
    return static_cast<int32_t>( u );
}

inline uint32_t packed_uint2( const uint8_t* b, int p1, int l1, int p2, int l2 ) // concatenate two split fields
{
    return ( packed_uint( b, p1, l1 ) << l2 ) | packed_uint( b, p2, l2 );
}

inline int32_t packed_int2( const uint8_t* b, int p1, int l1, int p2, int l2 ) // signed two-field concatenation
{
    if( packed_uint( b, p1, 1 ) )
    {
        return static_cast<int32_t>( ( packed_int( b, p1, l1 ) << l2 ) | packed_uint( b, p2, l2 ) );
    }
    return static_cast<int32_t>( packed_uint2( b, p1, l1, p2, l2 ) );
}

// --- unpacked: one bit per element (T = uint8_t or int) --------------------------------------------
template <class T>
uint64_t unpacked_uint( const T* bits, int start, int len )
{
    uint64_t v = 0;
    for( int i = 0; i < len; i++ )
    {
        v = ( v << 1 ) | ( static_cast<uint64_t>( bits[start + i] ) & 1u );
    }
    return v;
}

template <class T>
int64_t unpacked_int( const T* bits, int start, int len ) // two's-complement sign-extend
{
    uint64_t v = unpacked_uint( bits, start, len );
    if( len < 64 && ( v >> ( len - 1 ) ) )
    {
        v |= ~0ull << len;
    }
    return static_cast<int64_t>( v );
}

} // namespace bits

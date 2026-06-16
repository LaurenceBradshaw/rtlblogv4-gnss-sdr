#pragma once
#include <cstdint>

// GPS L1C CNAV-2 LDPC decoder (rate-1/2 binary LDPC, IS-GPS-800E). Self-contained log-domain
// sum-product (belief-propagation) decoder over the sparse parity-check matrix built from the
// IS-GPS-800E sub-matrix tables (gps_l1c_ldpc_tables.h). Mirrors PocketSDR src/sdr_ldpc.c
// (gen_B_LDPC_H + probability-propagation decode) but with our own decoder so there's no dependency
// on Neal's LDPC-codes library. See [[reference-repos]], [[gps-l1c-plan]] Phase 4b.
namespace gps::l1c::ldpc
{
// Decode one CNAV-2 subframe from its hard symbols (0/1, after de-interleaving). The first m bits of
// the systematic codeword are written to `msg`. Returns true iff every parity check is satisfied
// (a valid codeword) - the caller still confirms with the subframe CRC.
//   SF2: n=1200 symbols -> m=600 message bits.   SF3: n=548 symbols -> m=274 message bits.
bool decode_sf2( const uint8_t* syms, uint8_t* msg ); // msg[600]
bool decode_sf3( const uint8_t* syms, uint8_t* msg ); // msg[274]

} // namespace gps::l1c::ldpc

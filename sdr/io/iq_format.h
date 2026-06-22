#pragma once
#include <cstdint>
#include <vector>
#include "types.h" // Complex_buf

// IQ sample-format codec, independent of the source (file or live RTL-SDR): the on-disk/on-wire byte layout
// of complex I/Q samples and the conversion to/from the in-memory float Complex_buf. Used by the file reader
// (decode) and the recorder (encode), so neither has to depend on the other.
enum class Iq_sample_format
{
    UINT8,
    INT8,
    UINT16,
    INT16,
    FLOAT32,
};

// On-disk size (bytes) of one complex I/Q sample.
int bytes_per_iq_sample( Iq_sample_format format );

// The CLI/--format token for a format ("uint8" / "int8" / "uint16" / "int16" / "float32").
const char* iq_format_name( Iq_sample_format format );

// Decode `len` bytes of packed samples in `format` -> float [-1,1] I/Q. `output` must be pre-sized to the
// sample count (len / bytes_per_iq_sample). MSB-first integers are normalised to [-1,1]; float32 passes through.
void decode_iq( const unsigned char* buf, size_t len, Iq_sample_format format, Complex_buf& output );

// Encode float samples -> packed bytes in `format`: the exact INVERSE of decode_iq (round + clamp). Used by
// the recorder so a recording replays as a normal capture of that format.
std::vector<uint8_t> encode_iq( const Complex_buf& samples, Iq_sample_format format );

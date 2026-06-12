#pragma once

#include <complex>
#include <cstdint>
#include <vector>

using Complex_sample = std::complex<float>;
using Complex_buf    = std::vector<Complex_sample>;

struct Sample_block
{
    const Complex_sample* ptr1;
    size_t                len1;

    const Complex_sample* ptr2;
    size_t                len2;
};

using Sample_index = uint64_t;
using Satellite_id = uint32_t;

struct Ecef
{
    double x;
    double y;
    double z;
};
#pragma once
#include <atomic>
#include <chrono>
#include "types.h"

class Sample_buffer
{
public:
    // capacity must be a power of 2.
    // sample_rate_hz is used by is_available() to gate samples to real time:
    // sample i is visible only once i/sample_rate_hz seconds have elapsed
    // since the first push() - mirroring what hardware would deliver.
    explicit Sample_buffer( size_t capacity, double sample_rate_hz );

    // Called from the streaming thread
    void push( const Complex_buf& samples );

    // Called from channel threads (const, read-only)
    const Complex_sample& at( Sample_index index ) const;
    Sample_block          at( Sample_index start, size_t n ) const;

    bool is_available( Sample_index index ) const;

    Sample_index write_index() const;
    Sample_index newest_valid_index() const;
    Sample_index oldest_valid_index() const;

    void advance_oldest( Sample_index index );

private:
    size_t      capacity_;
    double      sample_rate_hz_;
    Complex_buf buffer_;

    std::atomic<Sample_index> write_index_ { 0 };
    std::atomic<Sample_index> oldest_valid_index_ { 0 };

    // Set on first push() with release, read in is_available() with acquire.
    std::atomic<bool>                     stream_started_ { false };
    std::chrono::steady_clock::time_point stream_start_; // written before stream_started_ publish

    size_t index_to_offset( Sample_index index ) const;
    bool   is_valid_range( Sample_index start, Sample_index end ) const;
};

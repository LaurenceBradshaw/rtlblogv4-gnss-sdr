#pragma once

#include <atomic>
#include <fstream>
#include <thread>
#include "stream_device.h"

enum class Iq_sample_format
{
    UINT8,
    INT8,
    UINT16,
    INT16,
    FLOAT32,
};

class Iq_file_device : public Stream_device
{
public:
    Iq_file_device( const std::string& path, uint32_t sample_rate_hz, Iq_sample_format sample_format );

    ~Iq_file_device();

    void start_streaming( Sample_callback callback ) override;
    void stop_streaming() override;
    bool is_streaming() const override;

    uint32_t sample_rate_hz() const override;

    uint64_t samples_consumed() const override;

private:
    void streaming_thread();

private:
    std::ifstream    file_;
    Iq_sample_format sample_format_;

    uint32_t sample_rate_hz_;

    Sample_callback   callback_;
    std::thread       worker_;
    std::atomic<bool> running_ { false };

    std::atomic<uint64_t> samples_consumed_ { 0 };
};
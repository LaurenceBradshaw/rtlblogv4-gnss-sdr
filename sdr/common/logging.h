#pragma once

#include <fmt/core.h>
#include <unistd.h>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <ostream>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>

namespace logging
{

enum class Level
{
    Info,
    Warning,
    Error,
    Debug
};

enum class Colour
{
    Reset,
    Red,
    Yellow,
    Green,
    Cyan
};

inline std::string to_string( Colour colour )
{
    switch( colour )
    {
    case Colour::Reset:
        return "\033[0m";
    case Colour::Red:
        return "\033[31m";
    case Colour::Yellow:
        return "\033[33m";
    case Colour::Green:
        return "\033[32m";
    case Colour::Cyan:
        return "\033[36m";
    default:
        assert( false );
        return "\033[0m";
    }
}

inline std::string to_string( Level level )
{
    switch( level )
    {
    case Level::Info:
        return "INFO";
    case Level::Warning:
        return "WARNING";
    case Level::Error:
        return "ERROR";
    case Level::Debug:
        return "DEBUG";
    default:
        assert( false );
        return "UNKNOWN";
    }
}

inline Colour to_colour( Level level )
{
    switch( level )
    {
    case Level::Info:
        return Colour::Reset;
    case Level::Warning:
        return Colour::Yellow;
    case Level::Error:
        return Colour::Red;
    case Level::Debug:
        return Colour::Cyan;
    default:
        assert( false );
        return Colour::Reset;
    }
}

std::filesystem::path get_log_dir();

void log( Level level, const std::string& str );
void stop_log();

} // namespace logging

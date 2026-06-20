#include <cxxopts.hpp>
#include <iostream>
#include <string>
#include "logging.h"
#include "receiver.h"
#include "signal_selection.h"
#ifdef SDR_GUI
#include "gui_app.h"
#endif

static Iq_sample_format parse_iq_format( const std::string& value )
{
    if( value == "uint8" )
        return Iq_sample_format::UINT8;
    if( value == "int8" )
        return Iq_sample_format::INT8;
    if( value == "uint16" )
        return Iq_sample_format::UINT16;
    if( value == "int16" )
        return Iq_sample_format::INT16;
    if( value == "float32" )
        return Iq_sample_format::FLOAT32;
    throw std::invalid_argument( "Unsupported IQ format: " + value );
}

// Excluded from the unit-test binary, which provides its own main() via Catch2.
#ifndef ENABLE_UNIT_TESTS
int main( int argc, char** argv )
{
    cxxopts::Options options( "sdr", "GNSS software receiver" );
    options.add_options()
        ( "file",        "IQ file path (file source)",
          cxxopts::value<std::string>() )
        ( "sample-rate", "Sample rate in Hz",
          cxxopts::value<uint32_t>()->default_value( "2048000" ) )
        ( "format",      "IQ sample format (file source): uint8 int8 uint16 int16 float32",
          cxxopts::value<std::string>()->default_value( "int16" ) )
        ( "rtlsdr",      "Use a live RTL-SDR front-end instead of a file source",
          cxxopts::value<bool>()->default_value( "false" ) )
        ( "device",      "RTL-SDR device index",
          cxxopts::value<int>()->default_value( "0" ) )
        ( "gain",        "RTL-SDR tuner gain in dB (negative = hardware AGC)",
          cxxopts::value<double>()->default_value( "-1" ) )
        ( "gui",         "Launch the graphical interface (requires a GUI build)",
          cxxopts::value<bool>()->default_value( "false" ) )
        ( "decimate",    "FIR-decimate the source by this integer factor before processing (1 = none). "
                         "Processing rate becomes sample-rate/factor. e.g. a 25 MHz capture with "
                         "--decimate 12 runs at ~2.083 MHz",
          cxxopts::value<uint32_t>()->default_value( "1" ) )
        ( "signal",      "Signal to search, repeatable: CONSTELLATION[:COMPONENT], where "
                         "CONSTELLATION=gps|galileo|beidou and COMPONENT=l1ca|l1c|e1|b1i. Omit the component "
                         "to search ALL of that constellation's components. e.g. --signal gps:l1ca "
                         "--signal galileo  (default: gps galileo)",
          cxxopts::value<std::vector<std::string>>() )
        ( "prns",        "PRN allowlist for a constellation (repeatable): CONSTELLATION:LIST, applied to "
                         "all of that constellation's selected components. e.g. --prns gps:1,4,6-32 "
                         "--prns galileo:3,11  (default: all PRNs)",
          cxxopts::value<std::vector<std::string>>() )
        ( "h,help",      "Show usage" );

    auto result = options.parse( argc, argv );

    if( result.count( "h" ) )
    {
        std::cout << options.help() << "\n";
        return 0;
    }

    const bool gui = result["gui"].as<bool>();

    Receiver_config config;
    config.use_rtlsdr     = result["rtlsdr"].as<bool>();
    config.sample_rate_hz = result["sample-rate"].as<uint32_t>();
    config.device_index   = result["device"].as<int>();
    config.gain_db        = result["gain"].as<double>();
    config.decimation     = std::max( 1u, result["decimate"].as<uint32_t>() );

    // The GUI's Source tab supplies the file / source params, so --file is optional under --gui.
    if( !gui && !config.use_rtlsdr && !result.count( "file" ) )
    {
        logging::log( logging::Level::Error, "IQ file path is required (--file), or use --rtlsdr" );
        return 1;
    }

    try
    {
        if( !config.use_rtlsdr && result.count( "file" ) )
        {
            config.file_path = result["file"].as<std::string>();
            config.format    = parse_iq_format( result["format"].as<std::string>() );
        }
        if( result.count( "signal" ) || result.count( "prns" ) )
        {
            const std::vector<std::string> sig =
                result.count( "signal" ) ? result["signal"].as<std::vector<std::string>>() : std::vector<std::string> {};
            const std::vector<std::string> prns =
                result.count( "prns" ) ? result["prns"].as<std::vector<std::string>>() : std::vector<std::string> {};
            config.selected_signals = parse_selection( sig, prns );
        }

        Receiver receiver( std::move( config ) );

        if( gui )
        {
#ifdef SDR_GUI
            return Gui_app::run( receiver, argc, argv );
#else
            logging::log( logging::Level::Error, "Built without the GUI - reconfigure with -DSDR_GUI=ON" );
            return 1;
#endif
        }

        receiver.run();
    }
    catch( const std::exception& e )
    {
        logging::log( logging::Level::Error, "Fatal - " + std::string( e.what() ) );
        return 1;
    }

    return 0;
}
#endif // ENABLE_UNIT_TESTS

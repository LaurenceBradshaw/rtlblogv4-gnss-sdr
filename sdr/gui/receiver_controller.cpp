#include "receiver_controller.h"
#include <exception>
#include <string>
#include "logging.h"
#include "receiver.h"

Receiver_controller::Receiver_controller( Receiver& receiver )
    : receiver_( receiver )
{
}

Receiver_controller::~Receiver_controller()
{
    stop();
}

void Receiver_controller::start()
{
    if( running_.load( std::memory_order_relaxed ) )
    {
        return;
    }
    // Reap a previous run that finished on its own (file drained / error) but was never joined -
    // assigning over a joinable std::thread would call std::terminate.
    if( worker_.joinable() )
    {
        worker_.join();
    }

    running_.store( true, std::memory_order_relaxed );
    worker_ = std::thread( [this] {
        try
        {
            receiver_.run(); // blocks: setup -> loop -> teardown
        }
        catch( const std::exception& e )
        {
            logging::log( logging::Level::Error, std::string( "Receiver stopped - " ) + e.what() );
        }
        running_.store( false, std::memory_order_relaxed ); // ended (stop / file drained / error)
    } );
}

void Receiver_controller::stop()
{
    receiver_.stop(); // ask run() to leave its loop
    if( worker_.joinable() )
    {
        worker_.join();
    }
    running_.store( false, std::memory_order_relaxed );
}

void Receiver_controller::apply_signal_selection( std::vector<Signal_selection> selection )
{
    receiver_.set_signal_selection( std::move( selection ) );
}

void Receiver_controller::apply_source_params( Source_params params )
{
    receiver_.set_source_params( std::move( params ) );
}

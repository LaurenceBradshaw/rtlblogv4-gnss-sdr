#pragma once
#include <atomic>
#include <thread>

class Receiver;

// Drives a Receiver on a worker thread for the GUI. The receiver does NOT auto-run: the Start button
// calls start() (spawns the run loop), Stop calls stop() (asks run() to exit, then joins). Restartable
// - each start() runs the receiver afresh (Receiver::setup() rebuilds the whole pipeline and resets the
// clock/EKF; an IQ file replays from the beginning). is_running() drives the button enabled state.
class Receiver_controller
{
public:
    explicit Receiver_controller( Receiver& receiver );
    ~Receiver_controller(); // stops + joins any running worker

    Receiver_controller( const Receiver_controller& )            = delete;
    Receiver_controller& operator=( const Receiver_controller& ) = delete;

    void start(); // run the receiver on a worker thread (no-op if already running)
    void stop();  // ask run() to exit and join the worker (no-op if not running)

    bool is_running() const
    {
        return running_.load( std::memory_order_relaxed );
    }

private:
    Receiver&         receiver_;
    std::thread       worker_;
    std::atomic<bool> running_ { false };
};

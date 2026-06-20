#pragma once

class Receiver;

// Runs the Qt application: shows the main window, drives the receiver on a worker thread, and tears
// it down cleanly when the window closes. One application per process, so this is a free function.
namespace Gui_app
{
// Blocks on the Qt event loop; returns the application exit code. The receiver is run on its own
// thread and stopped/joined before returning.
int run( Receiver& receiver, int argc, char** argv );
} // namespace Gui_app

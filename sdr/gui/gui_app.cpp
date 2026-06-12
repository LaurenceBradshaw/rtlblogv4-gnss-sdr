#include "gui_app.h"
#include <QApplication>
#include "main_window.h"
#include "receiver.h"
#include "receiver_controller.h"

int Gui_app::run( Receiver& receiver, int argc, char** argv )
{
    QApplication app( argc, argv );

    // The receiver does NOT auto-run - the user presses Start once they have the tab they want to
    // watch open. The controller owns the worker thread; the window's Start/Stop buttons drive it.
    Receiver_controller controller( receiver );

    Main_window window( receiver, controller );
    window.show();

    const int rc = app.exec();

    // Window closed: stop the receiver and join its worker before anything is destroyed.
    controller.stop();
    return rc;
}

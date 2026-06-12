#pragma once
#include "gui_frame.h"

class Receiver;

// The single point of contact between the GUI and the receiver: reads the Receiver's thread-safe
// published state and assembles a Gui_frame. Panels depend on Gui_frame, not on Receiver, so the
// receiver internals can change without touching any widget.
class Receiver_view
{
public:
    explicit Receiver_view( const Receiver& receiver );

    // Assemble the current frame. Cheap (a couple of locked copies); called once per GUI tick.
    Gui_frame poll() const;

private:
    const Receiver& receiver_;
};

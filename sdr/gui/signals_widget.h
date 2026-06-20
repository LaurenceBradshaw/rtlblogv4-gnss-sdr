#pragma once
#include <QString>
#include <vector>
#include "gui_panel.h"
#include "signal_selection.h"

class QCheckBox;
class QLabel;
class QLineEdit;
class Receiver;

// The "Signals" tab: choose which signals + PRNs the receiver searches. One PRN text box PER
// CONSTELLATION ("1,4,6-32"; blank = all - PRNs pick which satellites, a constellation-level choice),
// with a checkbox per signal component beneath it. The selection is applied on Start (Main_window reads
// selection() and stages it via the controller). Initialised from the receiver's current selection so
// it reflects any --signal / --prns CLI args.
class Signals_widget : public Gui_panel
{
    Q_OBJECT
public:
    explicit Signals_widget( const Receiver& receiver, QWidget* parent = nullptr );

    QString tab_title() const override
    {
        return QStringLiteral( "Signals" );
    }
    void update_frame( const Gui_frame& frame ) override; // static config - no per-frame work

    // Build the selection from the checked components, each with its constellation's PRN set. Throws
    // std::invalid_argument (message names the constellation) if a constellation with checked components
    // has a malformed PRN box. A blank PRN box means "all PRNs" for that constellation.
    std::vector<Signal_selection> selection() const;

    // Validate every PRN box (marks malformed ones with a red border + shows the first error in the
    // status line, or clears it). Returns true if all are valid; also caches that for is_valid(). Run
    // live on each keystroke and by the Start/tab-change paths before they read selection().
    bool validate();
    // The cached result of the last validate() - true if every PRN box currently parses.
    bool is_valid() const
    {
        return valid_;
    }

    // Lock the inputs while the receiver runs (changes apply on the next Start), and surface messages.
    void set_editable( bool on );
    void show_error( const QString& msg );
    void clear_status();

private:
    struct Component
    {
        Signal_id  id;
        QCheckBox* enabled = nullptr;
    };
    struct Group // one constellation: a shared PRN box + its components
    {
        Constellation          constellation = Constellation::Unknown;
        QLineEdit*             prns          = nullptr;
        std::vector<Component> components;
    };
    std::vector<Group> groups_;
    QLabel*            status_ = nullptr;
    bool               valid_  = true; // cached: all PRN boxes parse (updated by validate())
};

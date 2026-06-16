#pragma once
#include <QString>
#include "gui_panel.h"
#include "receiver.h" // Source_params, Iq_sample_format

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QRadioButton;
class QSpinBox;
class Receiver;

// The "Source" tab: choose the input - an IQ file (path / format / sample rate) or a live RTL-SDR
// front-end (device index / gain) - plus the decimation factor (both sources). Lets a run be set up
// entirely from the GUI, no CLI args. Applied on Start (Main_window reads source_params() and stages it
// via the controller). Initialised from the receiver's current config so it reflects any CLI args.
class Source_widget : public Gui_panel
{
    Q_OBJECT
public:
    explicit Source_widget( const Receiver& receiver, QWidget* parent = nullptr );

    QString tab_title() const override
    {
        return QStringLiteral( "Source" );
    }
    void update_frame( const Gui_frame& /*frame*/ ) override {} // static config - no per-frame work

    // Build the source params from the widgets (the sample rate is taken as-is; validate() first).
    Source_params source_params() const;

    // Validate the inputs (file present + a positive integer sample rate); shows the first error in the
    // status line (or clears it) and caches the result for is_valid(). Run live on edits and before Start.
    bool validate();
    bool is_valid() const
    {
        return valid_;
    }

    void set_editable( bool on ); // lock the inputs while the receiver runs
    void show_error( const QString& msg );
    void clear_status();

private:
    void update_source_visibility(); // show the file vs RTL-SDR group per the radio
    void update_gain_enabled();      // grey out the gain spin box while hardware AGC is selected
    void update_effective_rate();    // refresh the "processing rate = sample rate / decimation" read-out

    QRadioButton*   file_radio_   = nullptr;
    QRadioButton*   rtlsdr_radio_ = nullptr;
    QGroupBox*      file_group_   = nullptr;
    QGroupBox*      rtlsdr_group_ = nullptr;
    QLineEdit*      file_path_    = nullptr;
    QComboBox*      format_       = nullptr;
    QLineEdit*      sample_rate_  = nullptr; // native source rate (Hz); both sources
    QSpinBox*       decimation_   = nullptr; // FIR decimation factor; both sources
    QLabel*         effective_rate_ = nullptr; // live "processing rate after decimation" read-out
    QSpinBox*       device_index_ = nullptr;
    QDoubleSpinBox* gain_         = nullptr;
    QCheckBox*      agc_          = nullptr; // hardware AGC -> gain_db < 0
    QLabel*         status_       = nullptr;
    bool            valid_        = true;
};

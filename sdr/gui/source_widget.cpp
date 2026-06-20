#include "source_widget.h"
#include <algorithm>
#include <array>
#include <QButtonGroup>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QVBoxLayout>
#include "receiver.h"

namespace
{
// The selectable sample formats, in a fixed display order, each paired with the Iq_sample_format it maps
// to (stored as the combo item's data so the index order is free to change).
struct Format_option
{
    Iq_sample_format value;
    const char*      label;
};
const std::array<Format_option, 5> FORMATS = { {
    { Iq_sample_format::INT8, "int8 (8+8 bit)" },
    { Iq_sample_format::UINT8, "uint8" },
    { Iq_sample_format::INT16, "int16" },
    { Iq_sample_format::UINT16, "uint16" },
    { Iq_sample_format::FLOAT32, "float32" },
} };
} // namespace

Source_widget::Source_widget( const Receiver& receiver, QWidget* parent )
    : Gui_panel( parent )
{
    const Receiver_config& cfg = receiver.config();

    auto* root = new QVBoxLayout( this );

    // --- source selector ---
    file_radio_   = new QRadioButton( QStringLiteral( "IQ file" ), this );
    rtlsdr_radio_ = new QRadioButton( QStringLiteral( "RTL-SDR (live)" ), this );
    auto* src_group = new QButtonGroup( this );
    src_group->addButton( file_radio_ );
    src_group->addButton( rtlsdr_radio_ );
    auto* src_row = new QHBoxLayout;
    src_row->addWidget( new QLabel( QStringLiteral( "Source:" ), this ) );
    src_row->addWidget( file_radio_ );
    src_row->addWidget( rtlsdr_radio_ );
    src_row->addStretch();
    root->addLayout( src_row );

    // --- IQ file group ---
    file_group_   = new QGroupBox( QStringLiteral( "IQ file" ), this );
    auto* file_form = new QFormLayout( file_group_ );
    file_path_    = new QLineEdit( QString::fromStdString( cfg.file_path ), file_group_);
    auto* browse  = new QPushButton( QStringLiteral( "Browse..." ), file_group_ );
    auto* path_row = new QHBoxLayout;
    path_row->addWidget( file_path_ );
    path_row->addWidget( browse );
    file_form->addRow( QStringLiteral( "File:" ), path_row );
    format_ = new QComboBox( file_group_ );
    for( const Format_option& f : FORMATS )
    {
        format_->addItem( QString::fromUtf8( f.label ), static_cast<int>( f.value ) );
    }
    format_->setCurrentIndex( format_->findData( static_cast<int>( cfg.format ) ) );
    file_form->addRow( QStringLiteral( "Format:" ), format_ );
    root->addWidget( file_group_ );

    // --- RTL-SDR group ---
    rtlsdr_group_ = new QGroupBox( QStringLiteral( "RTL-SDR" ), this );
    auto* rtl_form = new QFormLayout( rtlsdr_group_ );
    device_index_ = new QSpinBox( rtlsdr_group_ );
    device_index_->setRange( 0, 32 );
    device_index_->setValue( cfg.device_index );
    rtl_form->addRow( QStringLiteral( "Device index:" ), device_index_ );
    gain_ = new QDoubleSpinBox( rtlsdr_group_ );
    gain_->setRange( 0.0, 60.0 );
    gain_->setSingleStep( 0.5 );
    gain_->setSuffix( QStringLiteral( " dB" ) );
    agc_ = new QCheckBox( QStringLiteral( "Hardware AGC" ), rtlsdr_group_ );
    const bool agc_on = cfg.gain_db < 0.0;
    agc_->setChecked( agc_on );
    gain_->setValue( agc_on ? 30.0 : cfg.gain_db );
    auto* gain_row = new QHBoxLayout;
    gain_row->addWidget( gain_ );
    gain_row->addWidget( agc_ );
    rtl_form->addRow( QStringLiteral( "Gain:" ), gain_row );
    root->addWidget( rtlsdr_group_ );

    // --- common: native sample rate + decimation ---
    auto* common      = new QFormLayout;
    sample_rate_      = new QLineEdit( QString::number( cfg.sample_rate_hz ), this );
    common->addRow( QStringLiteral( "Sample rate (Hz):" ), sample_rate_ );
    decimation_ = new QSpinBox( this );
    decimation_->setRange( 1, 64 );
    decimation_->setValue( static_cast<int>( std::max( 1u, cfg.decimation ) ) );
    decimation_->setToolTip( QStringLiteral( "FIR-decimate the source by this factor before processing (1 = none)." ) );
    common->addRow( QStringLiteral( "Decimation:" ), decimation_ );
    effective_rate_ = new QLabel( this );
    effective_rate_->setStyleSheet( QStringLiteral( "color: gray;" ) );
    common->addRow( QString(), effective_rate_ ); // sits under the decimation field, no label column
    hatch_ = new QCheckBox( QStringLiteral( "Hatch carrier-smoothing" ), this );
    hatch_->setChecked( cfg.hatch_enabled );
    hatch_->setToolTip( QStringLiteral(
        "Carrier-smooth the code pseudorange (lower noise, tighter fix). Reset on lock loss / re-acquire." ) );
    common->addRow( QStringLiteral( "Smoothing:" ), hatch_ );
    root->addLayout( common );

    status_ = new QLabel( this );
    status_->setWordWrap( true );
    root->addWidget( status_ );
    root->addStretch();

    file_radio_->setChecked( !cfg.use_rtlsdr );
    rtlsdr_radio_->setChecked( cfg.use_rtlsdr );

    // Wiring: browse a file, keep the file/RTL-SDR groups in step with the radio, grey the gain under AGC,
    // and re-validate live so the Start button enables/disables as inputs change.
    connect( browse, &QPushButton::clicked, this, [this] {
        const QString f = QFileDialog::getOpenFileName(
            this, QStringLiteral( "Choose IQ file" ), file_path_->text(),
            QStringLiteral( "IQ data (*.iq *.dat *.bin);;All files (*)" ) );
        if( !f.isEmpty() )
        {
            file_path_->setText( f );
            validate();
        }
    } );
    connect( file_radio_, &QRadioButton::toggled, this, [this] {
        update_source_visibility();
        validate();
    } );
    connect( agc_, &QCheckBox::toggled, this, [this] { update_gain_enabled(); } );
    connect( file_path_, &QLineEdit::textChanged, this, [this]( const QString& ) { validate(); } );
    connect( sample_rate_, &QLineEdit::textChanged, this, [this]( const QString& ) {
        validate();
        update_effective_rate();
    } );
    connect( decimation_, qOverload<int>( &QSpinBox::valueChanged ), this, [this]( int ) { update_effective_rate(); } );

    update_source_visibility();
    update_gain_enabled();
    update_effective_rate();
    validate();
}

void Source_widget::update_effective_rate()
{
    bool           rate_ok = false;
    const unsigned rate    = sample_rate_->text().trimmed().toUInt( &rate_ok );
    if( !rate_ok || rate == 0 )
    {
        effective_rate_->setText( QStringLiteral( "Processing rate: —" ) ); // em dash when rate invalid
        return;
    }
    const unsigned decim = static_cast<unsigned>( decimation_->value() );
    const double   proc  = static_cast<double>( rate ) / static_cast<double>( decim );
    effective_rate_->setText(
        ( decim == 1 ) ? QStringLiteral( "Processing rate: %1 MHz (no decimation)" ).arg( proc / 1e6, 0, 'f', 4 )
                       : QStringLiteral( "Processing rate: %1 MHz  (= %2 Hz / %3)" )
                             .arg( proc / 1e6, 0, 'f', 4 )
                             .arg( rate )
                             .arg( decim ) );
}

void Source_widget::update_source_visibility()
{
    const bool rtl = rtlsdr_radio_->isChecked();
    file_group_->setVisible( !rtl );
    rtlsdr_group_->setVisible( rtl );
}

void Source_widget::update_gain_enabled()
{
    gain_->setEnabled( !agc_->isChecked() );
}

Source_params Source_widget::source_params() const
{
    Source_params p;
    p.use_rtlsdr     = rtlsdr_radio_->isChecked();
    p.file_path      = file_path_->text().trimmed().toStdString();
    p.format         = static_cast<Iq_sample_format>( format_->currentData().toInt() );
    p.sample_rate_hz = sample_rate_->text().trimmed().toUInt(); // validate() guarantees this parses > 0
    p.device_index   = device_index_->value();
    p.gain_db        = agc_->isChecked() ? -1.0 : gain_->value();
    p.decimation     = static_cast<uint32_t>( decimation_->value() );
    p.hatch_enabled  = hatch_->isChecked();
    return p;
}

bool Source_widget::validate()
{
    clear_status();
    bool ok = true;

    if( file_radio_->isChecked() )
    {
        const QString path = file_path_->text().trimmed();
        if( path.isEmpty() )
        {
            show_error( QStringLiteral( "Choose an IQ file (or switch the source to RTL-SDR)." ) );
            ok = false;
        }
        else if( !QFileInfo::exists( path ) )
        {
            show_error( QStringLiteral( "File not found: %1" ).arg( path ) );
            ok = false;
        }
    }

    bool           rate_ok = false;
    const unsigned rate    = sample_rate_->text().trimmed().toUInt( &rate_ok );
    const bool     rate_valid = rate_ok && rate > 0;
    sample_rate_->setStyleSheet( rate_valid ? QString() : QStringLiteral( "border: 1px solid red;" ) );
    if( !rate_valid )
    {
        if( ok ) // don't clobber a more specific file error
        {
            show_error( QStringLiteral( "Sample rate must be a positive integer (Hz)." ) );
        }
        ok = false;
    }

    valid_ = ok;
    return ok;
}

void Source_widget::set_editable( bool on )
{
    file_radio_->setEnabled( on );
    rtlsdr_radio_->setEnabled( on );
    file_group_->setEnabled( on );
    rtlsdr_group_->setEnabled( on );
    sample_rate_->setEnabled( on );
    decimation_->setEnabled( on );
    hatch_->setEnabled( on );
    if( on )
    {
        update_gain_enabled(); // restore the AGC-driven gain enable state
    }
}

void Source_widget::show_error( const QString& msg )
{
    status_->setStyleSheet( QStringLiteral( "color: red;" ) );
    status_->setText( msg );
}

void Source_widget::clear_status()
{
    status_->clear();
    status_->setStyleSheet( QString() );
}

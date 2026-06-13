#include "signals_widget.h"
#include <QCheckBox>
#include <QFont>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>
#include <algorithm>
#include <stdexcept>
#include "constellations.h" // constellation_name
#include "gnss_format.h"    // constellation_color
#include "receiver.h"

Signals_widget::Signals_widget( const Receiver& receiver, QWidget* parent )
    : Gui_panel( parent )
{
    const std::vector<Signal_selection>& current = receiver.signal_selection();
    // Any selection entry for a constellation carries that constellation's (shared) PRN set.
    auto con_prns = [&current]( Constellation c ) -> const std::set<int>*
    {
        for( const Signal_selection& s : current )
        {
            if( s.id.constellation == c )
            {
                return &s.prns;
            }
        }
        return nullptr;
    };
    auto is_selected = [&current]( const Signal_id& id )
    {
        for( const Signal_selection& s : current )
        {
            if( s.id == id )
            {
                return true;
            }
        }
        return false;
    };

    auto* layout = new QVBoxLayout( this );

    auto* intro = new QLabel(
        QStringLiteral(
            "Signals to search. PRNs (per constellation, blank = all) pick which satellites; "
            "the checkboxes pick which signals."
        ),
        this
    );
    intro->setWordWrap( true );
    layout->addWidget( intro );

    // Group all selectable signals by constellation; one PRN box per group, a checkbox per component.
    for( const Signal_id& id : all_signals() )
    {
        Group* group = groups_.empty() || groups_.back().constellation != id.constellation ? nullptr : &groups_.back();
        if( group == nullptr )
        {
            // New constellation heading row: name (bold, coloured) + its shared PRN box.
            auto* headRow = new QHBoxLayout;
            auto* label   = new QLabel( constellation_name( id.constellation ), this );
            QFont f       = label->font();
            f.setBold( true );
            label->setFont( f );
            label->setStyleSheet(
                QStringLiteral( "color: %1;" ).arg( gui_format::constellation_color( id.constellation ).name() )
            );
            auto* prns = new QLineEdit( this );
            prns->setPlaceholderText( QStringLiteral( "PRNs e.g. 1,4,6-32 (blank = all)" ) );
            if( const std::set<int>* p = con_prns( id.constellation ) )
            {
                prns->setText( QString::fromStdString( format_prn_list( *p ) ) );
            }
            headRow->addWidget( label );
            headRow->addWidget( prns, 1 );
            layout->addLayout( headRow );

            groups_.push_back( { id.constellation, prns, {} } );
            group = &groups_.back();
        }

        auto* check = new QCheckBox( signal_name( id ), this );
        check->setChecked( is_selected( id ) );
        auto* row = new QHBoxLayout;
        row->setContentsMargins( 20, 0, 0, 0 ); // indent under the heading
        row->addWidget( check );
        row->addStretch( 1 );
        layout->addLayout( row );

        group->components.push_back( { id, check } );
    }

    layout->addStretch( 1 );

    status_ = new QLabel( this );
    status_->setWordWrap( true );
    layout->addWidget( status_ );

    // Live validation: re-check every PRN box as the user types.
    for( const Group& g : groups_ )
    {
        connect( g.prns, &QLineEdit::textChanged, this, [this]( const QString& ) { validate(); } );
    }
    validate(); // seed valid_ + clear any styling from the initial (CLI-seeded) values
}

bool Signals_widget::validate()
{
    QString first_error;
    for( const Group& g : groups_ )
    {
        bool ok = true;
        try
        {
            parse_prn_list( g.prns->text().toStdString(), max_prn_for( g.constellation ) ); // empty (= all) is valid
        }
        catch( const std::exception& e )
        {
            ok = false;
            if( first_error.isEmpty() )
            {
                first_error =
                    QString::fromStdString( std::string( constellation_name( g.constellation ) ) + " PRNs: " + e.what() );
            }
        }
        g.prns->setStyleSheet( ok ? QString() : QStringLiteral( "border: 1px solid #e03030;" ) );
    }
    valid_ = first_error.isEmpty();
    if( valid_ )
    {
        clear_status();
    }
    else
    {
        show_error( first_error );
    }
    return valid_;
}

void Signals_widget::update_frame( const Gui_frame& )
{
    // Static configuration; nothing to refresh per frame.
}

std::vector<Signal_selection> Signals_widget::selection() const
{
    std::vector<Signal_selection> out;
    for( const Group& g : groups_ )
    {
        // Only parse a constellation's PRN box if it actually has a checked component.
        const bool any = std::any_of(
            g.components.begin(), g.components.end(), []( const Component& c ) { return c.enabled->isChecked(); }
        );
        if( !any )
        {
            continue;
        }
        std::set<int> prns;
        try
        {
            prns = parse_prn_list( g.prns->text().toStdString(), max_prn_for( g.constellation ) );
        }
        catch( const std::exception& e )
        {
            throw std::invalid_argument( std::string( constellation_name( g.constellation ) ) + " PRNs: " + e.what() );
        }
        for( const Component& c : g.components )
        {
            if( c.enabled->isChecked() )
            {
                out.push_back( { c.id, prns } );
            }
        }
    }
    return out;
}

void Signals_widget::set_editable( bool on )
{
    for( const Group& g : groups_ )
    {
        g.prns->setEnabled( on );
        for( const Component& c : g.components )
        {
            c.enabled->setEnabled( on );
        }
    }
}

void Signals_widget::show_error( const QString& msg )
{
    status_->setStyleSheet( QStringLiteral( "color: #e03030; font-weight: bold;" ) );
    status_->setText( msg );
}

void Signals_widget::clear_status()
{
    status_->clear();
}

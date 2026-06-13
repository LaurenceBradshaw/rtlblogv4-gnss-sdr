#include "satellite_list_widget.h"
#include <QFont>
#include <QHBoxLayout>
#include <QLabel>
#include <QScrollArea>
#include <QVBoxLayout>
#include <set>
#include "cn0_bar_widget.h"
#include "gnss_format.h"
#include "gui_frame.h"
#include "satellite_detail_window.h"
#include "satellite_widget.h"

Satellite_list_widget::Satellite_list_widget( const Receiver& receiver, QWidget* parent )
    : Gui_panel( parent ),
      receiver_( receiver )
{
    auto* outer = new QVBoxLayout( this );

    summary_ = new QLabel( this );
    outer->addWidget( summary_ );

    // Body is split: the text list on the left half, the C/N0 bar chart on the right half.
    auto* body = new QHBoxLayout;
    outer->addLayout( body, 1 );

    // Left column: header + scrollable rows.
    auto* left   = new QVBoxLayout;
    auto* header = new QLabel( Satellite_widget::header_text(), this );
    header->setFont( QFont( QStringLiteral( "monospace" ) ) );
    header->setContentsMargins( 6, 0, 6, 0 );
    left->addWidget( header );

    // Rows live in a scroll area so the full constellation doesn't overflow the tab.
    auto* rows_container = new QWidget;
    rows_layout_         = new QVBoxLayout( rows_container );
    rows_layout_->setContentsMargins( 0, 0, 0, 0 );
    rows_layout_->setSpacing( 0 );
    rows_layout_->addStretch(); // trailing stretch; rows are inserted before it

    auto* scroll = new QScrollArea( this );
    scroll->setWidgetResizable( true );
    scroll->setWidget( rows_container );
    left->addWidget( scroll );

    body->addLayout( left, 1 );

    // Right column: the C/N0 bar chart (equal width to the list).
    cn0_bars_ = new Cn0_bar_widget( this );
    body->addWidget( cn0_bars_, 1 );

    update_frame( Gui_frame {} );
}

int Satellite_list_widget::key_of( Constellation constellation, int prn )
{
    return static_cast<int>( constellation ) * 1000 + prn;
}

void Satellite_list_widget::update_frame( const Gui_frame& frame )
{
    std::set<int>                      tracking_keys;
    std::set<int>                      present_keys; // SVs in this frame; rows for any others are removed
    int                                tracking_count = 0;
    std::map<int, Cn0_bar_widget::Bar> bars; // keyed -> SV-sorted left-to-right

    for( const Channel_snapshot& s : frame.channels )
    {
        const int key = key_of( s.constellation, static_cast<int>( s.satellite_id ) );
        present_keys.insert( key );

        Satellite_widget*& row = rows_[key];
        if( !row )
        {
            row = new Satellite_widget;
            connect( row, &Satellite_widget::clicked, this, &Satellite_list_widget::open_detail );
        }
        row->update_snapshot( s );

        if( s.state == Channel_state::TRACKING )
        {
            tracking_keys.insert( key );
            ++tracking_count;
            bars[key] = Cn0_bar_widget::Bar {
                QString::asprintf( "%s%02d", gui_format::constellation_prefix( s.constellation ), s.satellite_id ),
                s.cn0_db_hz,
                gui_format::constellation_color( s.constellation ),
                s.has_lock
            };
        }

        // Keep an open detail window for this SV refreshed with the live snapshot.
        if( auto it = detail_windows_.find( key ); it != detail_windows_.end() )
        {
            it->second->update_snapshot( s );
        }
    }

    // Drop rows for SVs no longer in the frame (e.g. the signal selection changed on a restart). While
    // running, every configured channel is always present, so this only fires on a configuration change.
    for( auto it = rows_.begin(); it != rows_.end(); )
    {
        if( present_keys.count( it->first ) == 0 )
        {
            rows_layout_->removeWidget( it->second );
            it->second->deleteLater();
            it = rows_.erase( it );
        }
        else
        {
            ++it;
        }
    }

    // Desired order: TRACKING satellites first, each group in SV order (std::map iterates sorted).
    std::vector<int> desired;
    desired.reserve( rows_.size() );
    for( const auto& [key, row] : rows_ )
    {
        if( tracking_keys.count( key ) )
        {
            desired.push_back( key );
        }
    }
    for( const auto& [key, row] : rows_ )
    {
        if( !tracking_keys.count( key ) )
        {
            desired.push_back( key );
        }
    }

    if( desired != current_order_ )
    {
        for( const auto& [key, row] : rows_ )
        {
            rows_layout_->removeWidget( row );
        }
        int index = 0;
        for( int key : desired )
        {
            rows_layout_->insertWidget( index++, rows_[key] );
        }
        current_order_ = std::move( desired );
    }

    summary_->setText( QString::asprintf( "Tracking %d of %zu satellites", tracking_count, rows_.size() ) );

    std::vector<Cn0_bar_widget::Bar> bar_list;
    bar_list.reserve( bars.size() );
    for( auto& [key, bar] : bars )
    {
        bar_list.push_back( std::move( bar ) );
    }
    cn0_bars_->set_bars( std::move( bar_list ) );
}

void Satellite_list_widget::open_detail( Constellation constellation, int prn )
{
    const int key = key_of( constellation, prn );

    if( auto it = detail_windows_.find( key ); it != detail_windows_.end() )
    {
        it->second->raise(); // already open - bring it to the front
        it->second->activateWindow();
        return;
    }

    // Parent to this widget's top-level window (so it is cleaned up on exit) but flagged as its own
    // window. It refreshes from update_frame and removes itself from the map when closed.
    auto* win            = new Satellite_detail_window( constellation, prn, receiver_, window() );
    detail_windows_[key] = win;
    connect(
        win,
        &Satellite_detail_window::closed,
        this,
        [this]( Constellation c, int p ) { detail_windows_.erase( key_of( c, p ) ); }
    );
    win->show();
}

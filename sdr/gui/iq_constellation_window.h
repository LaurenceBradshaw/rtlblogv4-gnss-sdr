#pragma once
#include "graph_window.h"

class Iq_heatmap_widget;

// A graph window showing one SV's live prompt I/Q constellation as a plasma heatmap.
class Iq_constellation_window : public Graph_window
{
    Q_OBJECT
public:
    Iq_constellation_window(
        Constellation constellation, int prn, Code code, const Receiver& receiver, QWidget* parent = nullptr
    );

protected:
    void refresh() override;

private:
    Iq_heatmap_widget* heatmap_ = nullptr;
};

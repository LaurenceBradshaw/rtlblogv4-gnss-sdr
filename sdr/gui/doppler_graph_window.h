#pragma once
#include "graph_window.h"

class Doppler_plot_widget;
class QComboBox;

// A graph window showing one SV's carrier Doppler over time, with a seconds/minutes toggle: the
// "seconds" view draws the fast ring (recent, high rate); "minutes" draws the slow per-second-mean ring.
class Doppler_graph_window : public Graph_window
{
    Q_OBJECT
public:
    Doppler_graph_window(
        Constellation constellation, int prn, Code code, const Receiver& receiver, QWidget* parent = nullptr
    );

protected:
    void refresh() override;

private:
    Doppler_plot_widget* plot_ = nullptr;
    QComboBox*           mode_ = nullptr; // 0 = seconds (fast), 1 = minutes (slow)
};

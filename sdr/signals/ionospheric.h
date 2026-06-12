#pragma once

// Broadcast ionosphere-correction data. It is a SYSTEM-WIDE model (not per-satellite): the receiver
// caches one and applies it to every satellite. Different constellations broadcast different models -
// GPS the Klobuchar 8-coefficient model, Galileo NeQuick-G (effective ionisation ai0..ai2 + flags),
// BeiDou Klobuchar / BDGIM - so `model` says which fields are populated. Only Klobuchar is decoded
// today; when another is added, give it its own fields here and dispatch on `model` at the use site.
struct Iono
{
    enum class Model
    {
        None,
        Klobuchar, // GPS LNAV SF4 page 18 (also BeiDou D1/D2)
    };

    Model model = Model::None;
    bool  valid = false;

    // Klobuchar coefficients (populated when model == Klobuchar).
    double alpha[4]     = {}; // amplitude coeffs (s, s/sc, s/sc^2, s/sc^3)
    double beta[4]      = {}; // period coeffs    (s, s/sc, s/sc^2, s/sc^3)
    int    leap_seconds = 0;  // dt_LS: GPS-UTC offset (s)
};

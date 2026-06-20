#!/usr/bin/env python3
"""Characterise the 2.5 MHz band-limit effect on the BOC subcarrier autocorrelation / double-estimator
subcarrier discriminator, for GPS L1C (TMBOC, 4/33 chips on BOC(6,1)) vs Galileo E1 (effectively BOC(1,1)
at <=2.5 MHz, since its CBOC BOC(6,1) is a small coherent amplitude on every chip and filters away cleanly).

It builds the TRUE transmitted code at a high rate, low-pass filters to the receiver's 2.5 MHz Nyquist
(1.25 MHz), correlates against the receiver's BOC(1,1) replica while sweeping the SUBCARRIER delay (code held
aligned), and reports the subcarrier correlation R(d), the early-late discriminator D(d), and where D crosses
zero. A clean BOC(1,1) should give ONE stable zero at d=0 (+/- the 1-element periodic lobes); false zeros at
fractional offsets are the +/-0.4-element subcarrier-loop equilibria seen on L1C.

Units: 1 chip = 2 "elements" (T_c); the sub-chip width T_s = 1 element = half a chip. Offsets are in elements.
"""
import numpy as np

SPC          = 48                      # samples per chip (>=24 to resolve BOC(6,1)'s 12 half-cycles/chip)
CHIP_HZ      = 1.023e6
FS           = SPC * CHIP_HZ           # high-rate model sample rate
NYQ_HZ       = 1.25e6                  # 2.5 MHz receiver -> 1.25 MHz Nyquist (the decimation anti-alias edge)
N_CHIPS      = 33 * 200                # 200 TMBOC blocks
TMBOC_POS    = {0, 4, 6, 29}           # L1Cp BOC(6,1) chip positions within each 33-chip block (IS-GPS-800)
RNG          = np.random.default_rng(1)


def subcarrier(boc_n, spc):
    """One chip of a BOC(boc_n,1) subcarrier sampled at spc samples/chip: square wave, boc_n full periods/chip."""
    t = (np.arange(spc) + 0.5) / spc            # mid-sample phase within the chip, [0,1)
    return np.sign(np.sin(2 * np.pi * boc_n * t))


SC_BOC1 = subcarrier(1, SPC)   # BOC(1,1) chip waveform  (= [-1..,+1..])
SC_BOC6 = subcarrier(6, SPC)   # BOC(6,1) chip waveform


def build(tmboc):
    """Return (true band-limited signal, replica code, replica BOC(1,1) subcarrier) at the high rate."""
    prn = RNG.choice([-1.0, 1.0], size=N_CHIPS)
    true = np.empty(N_CHIPS * SPC)
    rep_code = np.empty(N_CHIPS * SPC)   # +/-1 primary code, held flat across the chip (for the subcarrier sweep)
    rep_sc = np.empty(N_CHIPS * SPC)     # replica BOC(1,1) subcarrier
    for c in range(N_CHIPS):
        sl = slice(c * SPC, (c + 1) * SPC)
        use6 = tmboc and (c % 33) in TMBOC_POS
        true[sl] = prn[c] * (SC_BOC6 if use6 else SC_BOC1)
        rep_code[sl] = prn[c]
        rep_sc[sl] = SC_BOC1
    # band-limit the TRUE signal to the 2.5 MHz receiver Nyquist
    F = np.fft.rfft(true)
    f = np.fft.rfftfreq(true.size, d=1.0 / FS)
    F[f > NYQ_HZ] = 0.0
    true_bl = np.fft.irfft(F, n=true.size)
    return true_bl, rep_code, rep_sc


def acf(true_bl, rep_code, rep_sc, offsets_el):
    """Subcarrier correlation R(d): code aligned, replica subcarrier shifted by d elements (= d*SPC/2 samples)."""
    base = true_bl * rep_code            # de-spread by the (aligned) code -> leaves subcarrier(true) * code^2
    R = []
    for d in offsets_el:
        shift = int(round(d * SPC / 2))  # 1 element = SPC/2 samples
        R.append(np.dot(base, np.roll(rep_sc, shift)) / base.size)
    return np.array(R)


def zeros(off, D):
    out = []
    for i in range(len(D) - 1):
        if D[i] == 0 or D[i] * D[i + 1] < 0:
            x = off[i] - D[i] * (off[i + 1] - off[i]) / (D[i + 1] - D[i])
            out.append(x)
    return out


def code_acf(true_bl, rep_code, rep_sc, offsets_el):
    """Envelope correlation vs CODE delay: subcarrier held aligned (prompt), replica PRIMARY shifted by d
    elements. This is what the DE code/envelope DLL sees. (Shift rep_code by whole elements = SPC/2 samples.)"""
    base = true_bl * rep_sc              # de-rotate the subcarrier (aligned) -> leaves primary(true)*subcarrier^2
    R = []
    for d in offsets_el:
        shift = int(round(d * SPC / 2))
        R.append(np.dot(base, np.roll(rep_code, shift)) / base.size)
    return np.array(R)


off = np.arange(-2.5, 2.5001, 0.02)
s_el = 0.3                                  # subcarrier E/L spacing used in the tracker (elements)
s_code = 0.818                              # code E/L spacing in the tracker (corr_spacing*ci ~ 0.82 el @2.5MHz)
for name, tmboc in (("Galileo (BOC1,1)", False), ("GPS L1C (TMBOC 4/33)", True)):
    tb, rc, rs = build(tmboc)
    # ---- subcarrier loop ----
    Re = acf(tb, rc, rs, off + s_el); Rl = acf(tb, rc, rs, off - s_el)
    Dsub = np.abs(Re) - np.abs(Rl)
    zsub = [round(x, 3) for x in zeros(off, Dsub) if -0.6 < x < 0.6]
    # ---- code/envelope loop ----
    Rc = code_acf(tb, rc, rs, off)
    Rce = code_acf(tb, rc, rs, off + s_code); Rcl = code_acf(tb, rc, rs, off - s_code)
    Dcode = np.abs(Rce) - np.abs(Rcl)
    zcode = [round(x, 3) for x in zeros(off, Dcode) if -1.2 < x < 1.2]
    Rc = Rc / np.max(np.abs(Rc))
    print(f"\n=== {name} ===")
    print(f"  subcarrier disc zeros in (-0.6,0.6) el: {zsub}")
    print(f"  CODE/envelope disc zeros in (-1.2,1.2) el (spacing {s_code} el): {zcode}")
    print("  code envelope R(d) profile (d el : R):")
    for d in np.arange(-1.4, 1.401, 0.2):
        i = np.argmin(np.abs(off - d))
        bar = "#" * int(abs(Rc[i]) * 30)
        print(f"   {d:+.1f}: {Rc[i]:+.3f} {bar}")

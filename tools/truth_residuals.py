#!/usr/bin/env python3
"""Truth-residual analyzer for the TRUTH_LLH harness.

Run the receiver with TRUTH_LLH="lat,lon,alt" set; it emits TRUTHDIAG lines
(resid = pseudorange - true_geometric_range = c*rx_clock + per-SV error, pre-dedup).
This script removes the per-epoch common term (the receiver clock, via the median
across SVs in that epoch) and reports each SV/code's residual BIAS (mean) and
JITTER (std) - the clock-independent per-SV ranging error. That is the metric for
validating observable work (sub-sample DLL, code biases), NOT the noisy EKF position.

Usage:  ./bin/sdr ... 2>&1 | tools/truth_residuals.py
        tools/truth_residuals.py run.log
"""
import re, sys, statistics
from collections import defaultdict

CODE = {0: "CA", 1: "B", 2: "I", 3: "Cd"}
CON = {0: "GPS", 1: "GAL", 2: "BDS"}

line_re = re.compile(
    r"TRUTHDIAG rx=(\d+) con=(\d+) prn=\s*(\d+) code=(\d+) resid=(-?[\d.]+) range=([\d.]+)"
)

def main():
    src = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
    epochs = defaultdict(list)  # rx -> [(con,prn,code,resid)]
    for raw in src:
        line = re.sub(r"\x1b\[[0-9]*m", "", raw)
        m = line_re.search(line)
        if m:
            rx, con, prn, code, resid, _ = m.groups()
            epochs[int(rx)].append((int(con), int(prn), int(code), float(resid)))

    if not epochs:
        print("no TRUTHDIAG lines found (run with TRUTH_LLH set)", file=sys.stderr)
        return 1

    # Per epoch, subtract the common term PER CONSTELLATION (GPS and Galileo have different clock
    # references / inter-system bias), so each constellation's receiver clock is removed cleanly.
    per_sv = defaultdict(list)  # (con,prn,code) -> [cleaned residual]
    for rx, rows in epochs.items():
        med_con = {}
        for con in {r[0] for r in rows}:
            med_con[con] = statistics.median(r[3] for r in rows if r[0] == con)
        for con, prn, code, resid in rows:
            per_sv[(con, prn, code)].append(resid - med_con[con])

    print(f"epochs={len(epochs)}  (per-epoch common rx-clock removed via median)")
    print(f"{'SV':>8} {'code':>4} {'n':>5} {'bias(m)':>9} {'jitter(m)':>10}")
    for (con, prn, code) in sorted(per_sv):
        v = per_sv[(con, prn, code)]
        bias = statistics.mean(v)
        jit = statistics.pstdev(v) if len(v) > 1 else 0.0
        print(f"{CON.get(con,con):>5}{prn:>3} {CODE.get(code,code):>4} {len(v):>5} {bias:>9.1f} {jit:>10.1f}")

    # Same-SV cross-code differential (e.g. GPS CA vs Cd) - the inter-signal bias.
    biases = {(c, p, cd): statistics.mean(v) for (c, p, cd), v in per_sv.items()}
    pairs = [(c, p) for (c, p, cd) in biases if cd == 0 and (c, p, 3) in biases]
    if pairs:
        print("\nsame-SV CA-vs-Cd bias differential (m):")
        for c, p in sorted(pairs):
            print(f"  {CON.get(c,c)}{p}: {biases[(c,p,3)]-biases[(c,p,0)]:+.1f}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

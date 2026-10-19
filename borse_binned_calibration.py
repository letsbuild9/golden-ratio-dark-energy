"""Calibration of the binned-w(z) comparison (B1 in borse_v2b_tests.py).

Bin-averaging a smooth w(z) uniformly in z is only an approximation to what the piecewise-constant fit returns.
We therefore apply the identical averaging to each data set's own w0waCDM posterior-mean model (official DESI
chains) and compare: deviations shared by DESI's own CPL model are a feature of the averaging, not of the theory.

Run:  .venv/bin/python borse_binned_calibration.py   (output: rerun/binned_calibration.txt)
"""
import numpy as np
from scipy.integrate import quad
from scipy.stats import chi2 as chi2_dist, norm
import borse_desi_chains as DC

PHI = (1 + 5 ** 0.5) / 2
edges = [(0.0, 0.7), (0.7, 1.4), (1.4, 2.1)]


def binavg(w0, wa):
    return np.array([quad(lambda z: w0 + wa * z / (1 + z), a, b)[0] / (b - a) for a, b in edges])


def binned_stats(folder):
    ws, cols = [], []
    for k in range(1, 5):
        path = f"{DC.CHAINS}/{folder}/chain.{k}.txt"
        hdr = open(path).readline().lstrip("#").split()
        d = np.loadtxt(path)
        d = d[int(0.3 * len(d)):]
        ws.append(d[:, hdr.index("weight")])
        cols.append(np.column_stack([d[:, hdr.index(f"w_bin_{i}")] for i in range(3)]))
    w = np.concatenate(ws)
    S = np.vstack(cols)
    return np.average(S, axis=0, weights=w), np.cov(S.T, aweights=w)


def jsig(t, mu, cov):
    c = (t - mu) @ np.linalg.inv(cov) @ (t - mu)
    return c, norm.isf(chi2_dist.sf(c, 3) / 2)


print(f"{'data':12s} {'model':22s} {'bin avgs (0-0.7, 0.7-1.4, 1.4-2.1)':38s} {'chi2(3)':>8s} {'sigma':>6s}")
for b, c, label in (("binned_cmb", "cmb", "DESI+CMB"), ("binned_cmb_pantheonplus", "cmb_pantheonplus", "+Pantheon+"),
                    ("binned_cmb_union3", "cmb_union3", "+Union3"), ("binned_cmb_desy5sn", "cmb_desy5sn", "+DESY5")):
    mu, cov = binned_stats(b)
    w, a, bb = DC.load(c)
    w0m, wam = np.average(a, weights=w), np.average(bb, weights=w)
    for name, (w0, wa) in (("phi-ladder", (-PHI / 2, -1 / PHI)), (f"DESI CPL mean ({w0m:.2f},{wam:.2f})", (w0m, wam))):
        t = binavg(w0, wa)
        c2, s = jsig(t, mu, cov)
        print(f"{label:12s} {name:22s} {', '.join(f'{v:6.3f}' for v in t):38s} {c2:8.2f} {s:6.2f}")
    print(f"{'':12s} {'DESI binned (data)':22s} {', '.join(f'{v:6.3f}' for v in mu)}")

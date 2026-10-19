"""Test 1: where does the prediction sit in DESI's OFFICIAL w0-wa posteriors (full CMB likelihood)?

Chains: DESI DR2 BAO cosmology release (data.desi.lbl.gov/public/papers/y3/bao-cosmo-params),
model base_w_wa, CMB = Planck PR3 low-l TT + low-l EE + Planck PR4 CamSpec TTTEEE + Planck PR4/ACT DR6 lensing
(the default "CMB" combination of arXiv:2503.14738).

Method: weighted 2D histogram of (w0, wa) samples (30% burn-in removed), Gaussian-smoothed; the credible level
of a point is the posterior mass in regions denser than the point (i.e. the contour through it), converted to
a two-sided Gaussian sigma.  A Gaussian (mean + covariance) estimate is printed as a cross-check.

Run:  .venv/bin/python borse_desi_chains.py
"""

import os
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import norm, chi2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CHAINS = os.path.join(HERE, "data", "desi_chains")
OUT = os.path.join(HERE, "figures")
PHI = (1 + 5 ** 0.5) / 2

COMBOS = {"cmb": "DESI DR2 + CMB", "cmb_pantheonplus": "DESI DR2 + CMB + Pantheon+",
          "cmb_union3": "DESI DR2 + CMB + Union3", "cmb_desy5sn": "DESI DR2 + CMB + DESY5"}
POINTS = {"Borse d=2 (-phi/2, -1/phi)": (-PHI / 2, -1 / PHI),
          "Borse d=3 (-1+phi^-3, -1/phi)": (-1 + PHI ** -3, -1 / PHI),
          "non-crossing (-phi/2, -0.191)": (-PHI / 2, -(1 - PHI / 2)),
          "LCDM (-1, 0)": (-1.0, 0.0)}


def load(folder, burn=0.3):
    cols, ws, w0, wa = None, [], [], []
    for k in range(1, 5):
        path = os.path.join(CHAINS, folder, f"chain.{k}.txt")
        with open(path) as fh:
            header = fh.readline().lstrip("#").split()
        d = np.loadtxt(path)
        d = d[int(burn * len(d)):]
        ws.append(d[:, header.index("weight")])
        w0.append(d[:, header.index("w")])
        wa.append(d[:, header.index("wa")])
    return np.concatenate(ws), np.concatenate(w0), np.concatenate(wa)


def density_grid(w, x, y, n=220, smooth=2.0):
    xr = (x.min() - 0.05, x.max() + 0.05)
    yr = (y.min() - 0.1, y.max() + 0.1)
    H, xe, ye = np.histogram2d(x, y, bins=n, range=[xr, yr], weights=w)
    H = gaussian_filter(H, smooth)
    H /= H.sum()
    return H, xe, ye


def level_of(H, xe, ye, pt):
    i = np.searchsorted(xe, pt[0]) - 1
    j = np.searchsorted(ye, pt[1]) - 1
    if not (0 <= i < H.shape[0] and 0 <= j < H.shape[1]):
        return 1.0
    dens = H[i, j]
    return H[H > dens].sum()


def to_sigma(mass):
    mass = min(mass, 1 - 1e-15)
    return norm.isf((1 - mass) / 2)


def contour_levels(H, masses):
    flat = np.sort(H.ravel())[::-1]
    cum = np.cumsum(flat)
    return [flat[np.searchsorted(cum, m)] for m in masses]


def main():
    os.makedirs(OUT, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6), sharey=True)
    for ax, (folder, label) in zip(axes, COMBOS.items()):
        w, x, y = load(folder)
        H, xe, ye = density_grid(w, x, y)
        mean = np.array([np.average(x, weights=w), np.average(y, weights=w)])
        cov = np.cov(np.vstack([x, y]), aweights=w)
        icov = np.linalg.inv(cov)
        print("\n" + "=" * 84 + f"\n{label}   ({len(x)} samples after burn-in)\n" + "=" * 84)
        print(f"posterior mean w0 = {mean[0]:.3f} +- {np.sqrt(cov[0, 0]):.3f},  wa = {mean[1]:.3f} +- {np.sqrt(cov[1, 1]):.3f},"
              f"  corr = {cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]):.2f}")
        print(f"{'point':32s} {'contour level':>14s} {'sigma (2D)':>11s} {'Gaussian sigma':>15s}")
        for name, pt in POINTS.items():
            m = level_of(H, xe, ye, pt)
            d = np.array(pt) - mean
            g = norm.isf(chi2.sf(d @ icov @ d, 2) / 2)
            print(f"{name:32s} {100 * m:13.2f}% {to_sigma(m):11.2f} {g:15.2f}")
        xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
        lv = contour_levels(H, [0.6827, 0.9545, 0.9973])
        ax.contourf(xc, yc, H.T, levels=[lv[2], lv[1], lv[0], H.max()], colors=["#dbe6f5", "#a9c4e8", "#3b6fb6"])
        ax.contour(xc, yc, H.T, levels=[lv[2], lv[1], lv[0]], colors="k", linewidths=0.6)
        ax.plot(*POINTS["Borse d=2 (-phi/2, -1/phi)"], "r*", ms=15, label="$\\varphi$-ladder")
        ax.plot(*POINTS["LCDM (-1, 0)"], "kx", ms=9, label="$\\Lambda$CDM")
        ax.set(title=label, xlabel="$w_0$", xlim=(-1.15, -0.35), ylim=(-2.4, 0.4))
    axes[0].set_ylabel("$w_a$")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Official DESI DR2 posteriors (full Planck + ACT lensing CMB): 68 / 95 / 99.7% regions")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig8_desi_official_chains.png"), dpi=130)
    print("\nSaved figures/fig8_desi_official_chains.png")


if __name__ == "__main__":
    main()

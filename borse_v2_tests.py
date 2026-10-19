"""Version-2 robustness tests (REVISION_PLAN.md, Phase 2: T1-T6, T8, T11, T12).

T6  Table 1 extended with DES-Dovekie, plus absolute goodness of fit (chi2/dof, p-value).
T3  Look-elsewhere test repeated with DES-Dovekie.
T1  Free-x fit: ladder curve w0 = -1 + x^2/2, wa = -x with x free (compressed-CMB pipeline, 5 combinations).
T2  Control constants on the same curve (1/sqrt2, 2/3, 1/e, ln2, pi/5, 0.6, ...).
T1b Posterior density along the same curve in the official DESI chains (full CMB).
T8  Position of the phi-ladder point in alternative official DESI chains (other CMB likelihoods, no CMB, SN cut).
T4  Savage-Dickey Bayes factors from the official chains (DESI priors: w0 in U[-3,1], wa in U[-3,2], w0+wa<0).
T5  DESI pivot-redshift test and fixed-redshift w(z=0.3), w(z=0.5).
T12 Background history implied by the phi-ladder w(a).
T11 Effective-gravity mu0 for a constant growth index gamma = 1/phi.

Run:  .venv/bin/python borse_v2_tests.py      (about 15-20 min; output in rerun/v2_results.txt, figures/fig9_x_profile.png)
"""

import os
import numpy as np
from scipy.stats import chi2 as chi2_dist, norm
from scipy.optimize import brentq, least_squares
from scipy.integrate import quad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from borse_fit import BAO, CMB, SN, MODELS, PHI, OMEGA_NU, OMEGA_R, fit, hubble, sigma_from_dchi2, OUT
import borse_desi_chains as DC

X = 1 / PHI
LADDER = lambda x: (-1 + x ** 2 / 2, -x)
PRIOR_AREA = 4 * 5 - 4.5  # w0 in [-3,1] x wa in [-3,2] minus the region w0 + wa >= 0 (area 4.5)


def head(t):
    print("\n" + "=" * 96 + f"\n{t}\n" + "=" * 96, flush=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    bao, cmb = BAO(), CMB()
    sns = {n: SN(n) for n in ("PantheonPlus", "Union3", "DESY5", "DESDovekie")}
    combos = {"DESI+CMB": [bao, cmb]}
    combos.update({f"DESI+CMB+{n}": [bao, cmb, s] for n, s in sns.items()})

    # ------------------------------------------------------------- T6 table + goodness of fit
    head("T6. Model comparison incl. DES-Dovekie, with absolute goodness of fit")
    print(f"{'data':24s} {'model':10s} {'k':>2s} {'chi2':>9s} {'dof':>5s} {'chi2/dof':>8s} {'p-value':>8s} {'AIC':>9s} {'BIC':>9s} {'nsig vs free':>12s}")
    best = {}
    for cname, likes in combos.items():
        N = sum(l.n for l in likes)
        nsn = sum(1 for l in likes if isinstance(l, SN))
        res = {}
        for m in ("LCDM", "Borse d=2", "Borse d=3", None):
            k = 5 if m is None else 3
            c2, p = fit(likes, MODELS[m] if m else None)
            res[m or "w0waCDM"] = (c2, p, k)
        best[cname] = res
        free = res["w0waCDM"][0]
        for m, (c2, p, k) in res.items():
            dof = N - k - nsn  # analytic marginalisation of the SN absolute magnitude removes one dof per SN sample
            ns = "--" if m == "w0waCDM" else f"{sigma_from_dchi2(c2 - free, 2):.2f}"
            print(f"{cname:24s} {m:10s} {k:2d} {c2:9.2f} {dof:5d} {c2 / dof:8.3f} {chi2_dist.sf(c2, dof):8.3f} "
                  f"{c2 + 2 * k:9.2f} {c2 + k * np.log(N):9.2f} {ns:>12s}")
        print(f"{'':24s} free best fit: w0 = {res['w0waCDM'][1][3]:.3f}, wa = {res['w0waCDM'][1][4]:.3f}")

    # ------------------------------------------------------------- T3 look-elsewhere incl. Dovekie
    head("T3. Look-elsewhere test (403 ladder-family pairs), now incl. DES-Dovekie")
    fam = sorted({round(s * c * X ** d + o, 6) for d in range(1, 5) for c in (1, .5) for s in (1, -1) for o in (0, 1, -1)})
    pairs = [(a, b) for a in fam if -1.3 <= a <= -0.3 for b in fam if -2.5 <= b <= 0.5 and a + b < 0]
    within = {}
    for n in ("PantheonPlus", "Union3", "DESY5", "DESDovekie"):
        cname = f"DESI+CMB+{n}"
        path = os.path.join("data", f"ladder_{cname}.npy")
        if n == "DESDovekie" or not os.path.exists(path):
            chis = np.array([fit(combos[cname], pr, quick=True)[0] for pr in pairs])
            np.save(path, chis)
        chis = np.load(path)
        free = best[cname]["w0waCDM"][0]
        borse = best[cname]["Borse d=2"][0]
        within[n] = (chis - free) < 2.30
        print(f"{cname:24s} phi-ladder rank {int(np.sum(chis < borse - 1e-6)) + 1:3d}/{len(pairs)};"
              f"  pairs within 1 sigma of free fit: {within[n].sum():3d}")
    for label, keys in (("Pantheon+, Union3, DESY5", ("PantheonPlus", "Union3", "DESY5")),
                        ("Pantheon+, Union3, DES-Dovekie", ("PantheonPlus", "Union3", "DESDovekie"))):
        allm = np.logical_and.reduce([within[k] for k in keys])
        print(f"within 1 sigma for all of [{label}]: {[(round(pairs[i][0], 3), round(pairs[i][1], 3)) for i in np.where(allm)[0]]}")

    # ------------------------------------------------------------- T1/T2 free-x profile
    head("T1/T2. Free x on the ladder curve w0 = -1 + x^2/2, wa = -x (compressed-CMB pipeline)")
    xs = np.round(np.arange(0.30, 1.0001, 0.01), 3)
    named = {"1/e": 1 / np.e, "0.5": 0.5, "0.6": 0.6, "1/phi": X, "pi/5": np.pi / 5, "2/3": 2 / 3,
             "ln2": np.log(2), "1/sqrt2": 1 / np.sqrt(2), "0.75": 0.75}
    prof = {}
    for cname, likes in combos.items():
        chis = np.array([fit(likes, LADDER(x), quick=True)[0] for x in xs])
        prof[cname] = chis
        i = int(np.argmin(chis))
        lo = xs[:i + 1][np.where(chis[:i + 1] - chis[i] > 1)[0][-1]] if np.any(chis[:i + 1] - chis[i] > 1) else np.nan
        hi = xs[i:][np.where(chis[i:] - chis[i] > 1)[0][0]] if np.any(chis[i:] - chis[i] > 1) else np.nan
        # refine the minimum with a parabola through the three lowest grid points
        j = min(max(i, 1), len(xs) - 2)
        c = np.polyfit(xs[j - 1:j + 2], chis[j - 1:j + 2], 2)
        xb = -c[1] / (2 * c[0])
        print(f"\n{cname}: best x = {xb:.3f}  (68% interval on 0.01 grid: [{lo:.2f}, {hi:.2f}]);"
              f"  1/phi is at dchi2 = {np.interp(X, xs, chis) - chis.min():.2f}")
        row = []
        for nm, v in named.items():
            row.append(f"{nm}={np.interp(v, xs, chis) - chis.min():.2f}")
        print("   dchi2 at named constants: " + ", ".join(row))
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for cname, chis in prof.items():
        ax.plot(xs, chis - chis.min(), label=cname)
    for nm, v in named.items():
        ax.axvline(v, color="r" if nm == "1/phi" else "grey", ls="-" if nm == "1/phi" else ":", lw=1.2 if nm == "1/phi" else 0.7)
        ax.text(v, 9.3, nm, rotation=90, fontsize=7, ha="right", va="top", color="r" if nm == "1/phi" else "grey")
    ax.axhline(1, color="k", lw=0.6)
    ax.axhline(4, color="k", lw=0.6, ls="--")
    ax.set(xlabel="$x$ in $w_0=-1+x^2/2$, $w_a=-x$", ylabel="$\\Delta\\chi^2$", ylim=(0, 10),
           title="Data-preferred value of $x$ (compressed-CMB fits)")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig9_x_profile.png"), dpi=130)

    # ------------------------------------------------------------- official chains: T1b, T8, T4, T5
    head("T1b/T8/T4/T5. Official DESI DR2 chains")
    chainsets = {
        "cmb": "DESI+CMB (CamSpec) [baseline]", "cmb_pantheonplus": "+Pantheon+ [baseline]",
        "cmb_union3": "+Union3 [baseline]", "cmb_desy5sn": "+DESY5 [baseline]",
        "plik": "DESI+CMB (Planck PR3 plik)", "plik_desy5": "+DESY5 (plik)", "hillipop_desy5": "+DESY5 (HiLLiPoP/LoLLiPoP)",
        "cmb_desy5_zmin01": "+DESY5 with z>0.1", "bao_only": "DESI BAO only (no CMB)", "bao_desy5": "DESI BAO + DESY5 (no CMB)",
    }
    xg = np.linspace(0.3, 1.0, 141)
    print(f"{'chains':32s} {'level':>7s} {'nsig':>5s} | {'x_best (curve)':>15s} | {'BF(phi:free)':>12s} {'BF(LCDM:free)':>13s} | "
          f"{'z_p':>5s} {'w_p data':>14s} {'w_p th':>7s} {'dev':>5s} | {'w(0.3) data':>14s} {'dev':>5s} | {'w(0.5) data':>14s} {'dev':>5s}")
    w03, w05 = (lambda z: LADDER(X)[0] + LADDER(X)[1] * z / (1 + z))(0.3), (lambda z: LADDER(X)[0] + LADDER(X)[1] * z / (1 + z))(0.5)
    for folder, label in chainsets.items():
        if not os.path.exists(os.path.join(DC.CHAINS, folder, "chain.4.txt")):
            print(f"{label:32s} (chains not downloaded)")
            continue
        w, a, b = DC.load(folder)
        H, xe, ye = DC.density_grid(w, a, b)
        lev = DC.level_of(H, xe, ye, LADDER(X))
        # density along the ladder curve
        dens = np.array([H[min(max(np.searchsorted(xe, LADDER(v)[0]) - 1, 0), H.shape[0] - 1),
                           min(max(np.searchsorted(ye, LADDER(v)[1]) - 1, 0), H.shape[1] - 1)] for v in xg])
        k = int(np.argmax(dens))
        inside = xg[dens >= dens[k] * np.exp(-0.5)]
        xbest = f"{xg[k]:.3f} [{inside.min():.2f},{inside.max():.2f}]" if dens[k] > 0 else "n/a"
        # Savage-Dickey: posterior pdf at the point / prior pdf, for smoothing 1.5, 2, 3 bins
        bfs, bfl = [], []
        for sm in (1.5, 2.0, 3.0):
            Hs, xs_, ys_ = DC.density_grid(w, a, b, smooth=sm)
            area = (xs_[1] - xs_[0]) * (ys_[1] - ys_[0])
            pdf = lambda pt: Hs[np.searchsorted(xs_, pt[0]) - 1, np.searchsorted(ys_, pt[1]) - 1] / area \
                if (xs_[0] < pt[0] < xs_[-1] and ys_[0] < pt[1] < ys_[-1]) else 0.0
            bfs.append(pdf(LADDER(X)) * PRIOR_AREA)
            bfl.append(pdf((-1.0, 0.0)) * PRIOR_AREA)
        bf_phi = f"{min(bfs):.1f}-{max(bfs):.1f}"
        bf_l = f"{min(bfl):.2f}-{max(bfl):.2f}" if max(bfl) > 0 else "<~0 (not sampled)"
        # pivot and fixed-z
        c = np.cov(np.vstack([a, b]), aweights=w)
        ap = 1 + c[0, 1] / c[1, 1]
        wp = a + b * (1 - ap)
        mwp, swp = np.average(wp, weights=w), np.sqrt(np.average((wp - np.average(wp, weights=w)) ** 2, weights=w))
        wpt = LADDER(X)[0] + LADDER(X)[1] * (1 - ap)
        out = []
        for z, wt in ((0.3, w03), (0.5, w05)):
            wz = a + b * z / (1 + z)
            m, s = np.average(wz, weights=w), np.sqrt(np.average((wz - np.average(wz, weights=w)) ** 2, weights=w))
            out.append(f"{m:7.3f}+-{s:5.3f} {(wt - m) / s:5.2f}")
        print(f"{label:32s} {100 * lev:6.1f}% {DC.to_sigma(lev):5.2f} | {xbest:>15s} | {bf_phi:>12s} {bf_l:>13s} | "
              f"{1 / ap - 1:5.2f} {mwp:7.3f}+-{swp:5.3f} {wpt:7.3f} {(wpt - mwp) / swp:5.2f} | {out[0]} | {out[1]}")
    print(f"\nTheory: w(z=0.3) = {w03:.4f}, w(z=0.5) = {w05:.4f}.  Prior area for Bayes factors = {PRIOR_AREA}"
          " (DESI DR2: w0 in U[-3,1], wa in U[-3,2], w0+wa<0). BF > 1 favours the fixed point over free (w0,wa);"
          " ranges show sensitivity to histogram smoothing. BF(phi:LCDM) = BF(phi:free)/BF(LCDM:free).")

    # ------------------------------------------------------------- T12 background history
    head("T12. Background history implied by the phi-ladder w(a) (CPL, flat; Om from the DESY5 phi-ladder fit)")
    p = best["DESI+CMB+DESY5"]["Borse d=2"][1]
    Om = (p[0] + OMEGA_NU) / p[2] ** 2
    w0, wa = LADDER(X)
    w = lambda a: w0 + wa * (1 - a)
    rde = lambda a: a ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
    E2 = lambda a: Om * a ** -3 + (1 - Om) * rde(a)
    q = lambda a: 0.5 * (Om * a ** -3 + (1 - Om) * rde(a) * (1 + 3 * w(a))) / E2(a)
    apk = brentq(lambda a: 1 + w(a), 0.3, 1.0)
    aacc = brentq(q, 0.3, 1.0)
    aeq = brentq(lambda a: Om * a ** -3 - (1 - Om) * rde(a), 0.3, 1.0)
    aend = brentq(q, 1.0, 6.0)
    aeq2 = brentq(lambda a: Om * a ** -3 - (1 - Om) * rde(a), 1.5, 20)
    aw1 = 1 - (1 - w0) / wa
    lb = quad(lambda a: 1 / (a * np.sqrt(E2(a))), apk, 1)[0] * 977.8 / (100 * p[2])
    print(f"Om = {Om:.4f}, H0 = {100 * p[2]:.2f}")
    print(f"dark-energy density maximum at z = {1 / apk - 1:.4f} (equals z_cross); lookback time {lb:.2f} Gyr")
    print(f"dark energy overtakes matter at z = {1 / aeq - 1:.3f}; acceleration begins at z = {1 / aacc - 1:.3f}")
    print(f"quintom-B check: w0 = {w0:.4f} > -1 and w0 + wa = {w0 + wa:.4f} < -1 -> {w0 > -1 and w0 + wa < -1}")
    print(f"CPL extrapolated to the future: acceleration ends at a = {aend:.2f}; matter dominates again at a = {aeq2:.2f};"
          f" w exceeds +1 at a = {aw1:.2f} (extrapolation invalid)")

    # ------------------------------------------------------------- T11 mu0 for constant gamma = 1/phi
    head("T11. Effective-gravity mu0 for a constant growth index gamma = 1/phi (phi-ladder background, DESY5 fit)")
    from borse_model import DragModel
    gr = DragModel(p, 0.0)
    z = gr.z
    N = -np.log(1 + z)
    Omz = gr.Om
    e = 1e-4
    dlnH = (np.log(hubble((1 + z) * np.exp(-e) - 1, *p)) - np.log(hubble((1 + z) * np.exp(e) - 1, *p))) / (2 * e)
    mug = (np.gradient(gr.f, N) + gr.f ** 2 + (2 + dlnH) * gr.f) / (1.5 * Omz)
    ode = p[2] ** 2 - p[0] - OMEGA_NU - OMEGA_R
    a_ = 1 / (1 + z)
    Ode = ode * a_ ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a_)) / (hubble(z, *p) / 100) ** 2
    sel = (z < 3) & (z > 0.02)
    for lab, f in (("gamma = 1/phi (constant)", Omz ** X), ("gamma = gamma_GR + F(4)", gr.f * Omz ** (X ** 4 / 2))):
        mu = (np.gradient(f, N) + f ** 2 + (2 + dlnH) * f) / (1.5 * Omz) / mug
        r = mu[sel] - 1
        mu0 = least_squares(lambda m: r - m[0] * Ode[sel] / Ode[0], [0.0]).x[0]
        vals = ", ".join(f"z={zz}: {mu[np.argmin(abs(z - zz))] - 1:+.3f}" for zz in (0, 0.5, 1, 2))
        print(f"{lab:26s}: mu - 1 at {vals};  mu0 = {mu0:+.3f}; max residual {np.max(abs(r - mu0 * Ode[sel] / Ode[0])):.4f}")


if __name__ == "__main__":
    main()

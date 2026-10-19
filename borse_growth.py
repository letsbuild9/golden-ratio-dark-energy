"""Growth-index test: does gamma = 1/phi fit measured f*sigma8(z)?

Growth model: f(z) = Omega_m(z)^gamma,  dlnD/dlna = f,  f*sigma8(z) = f(z) sigma8_0 D(z)/D(0).
Backgrounds: best fits to DESI DR2 + CMB + DESY5 from borse_fit.py (LCDM and Borse d=2).

Two tests
  A. Shape only:  sigma8_0 free  -> tests how f*sigma8 evolves with z, independent of the CMB.
  B. CMB-anchored: amplitude fixed at z=30 by Planck 2018 (sigma8 = 0.8111 +- 0.006 in LCDM),
     then grown to z=0 by each model -> tests the S8 claim as well.

Data (f*sigma8 only; the BOSS/eBOSS sub-covariances are the marginals of their DM, DH, fs8 fits):
  6dFGS z=0.067 (Beutler+2012), SDSS MGS z=0.15 (Howlett+2015; 0.53 +- 0.16 as tabulated in eBOSS DR16 Table 3),
  BOSS DR12 z=0.38, 0.51 + eBOSS LRG z=0.698 (joint covariance), eBOSS ELG z=0.85 (Gaussian approx),
  eBOSS QSO z=1.48   (SDSS DR16 cosmology, arXiv:2007.08991; files from Cobaya bao_data)

Run:  .venv/bin/python borse_growth.py          (SDSS-era data)
      .venv/bin/python borse_growth.py --desi   (DESI DR1 full-shape + peculiar velocities)
"""

import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, minimize_scalar, brentq
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from borse_fit import BAO, CMB, SN, MODELS, OMEGA_NU, OMEGA_R, PHI, fit, hubble, DATA, OUT

GAMMA_BORSE = 1 / PHI
Z_INIT = 30.0
SIG8_PLANCK, SIG8_PLANCK_ERR = 0.8111, 0.0060
PLANCK_BG = (0.1200 + 0.02237, 0.02237, 0.6736, -1.0, 0.0)  # Planck 2018 TT,TE,EE+lowE+lensing (same as sigma8 = 0.8111)


# ------------------------------------------------------------------ data
def load_fs8():
    z, v, blocks = [], [], []
    for name in ("LRG", "QSO"):
        rows = [l.split() for l in open(os.path.join(DATA, f"sdss_DR16_BAOplus_{name}_FSBAO_DMDHfs8.dat"))]
        cov = np.loadtxt(os.path.join(DATA, f"sdss_DR16_BAOplus_{name}_FSBAO_DMDHfs8_covtot.txt"))
        idx = [i for i, r in enumerate(rows) if r[2] == "f_sigma8"]
        z += [float(rows[i][0]) for i in idx]
        v += [float(rows[i][1]) for i in idx]
        blocks.append(cov[np.ix_(idx, idx)])
    extra = [(0.067, 0.423, 0.055), (0.15, 0.53, 0.16), (0.85, 0.315, 0.095)]
    z += [e[0] for e in extra]
    v += [e[1] for e in extra]
    blocks += [np.array([[e[2] ** 2]]) for e in extra]
    n = len(z)
    cov = np.zeros((n, n))
    i = 0
    for b in blocks:
        k = len(b)
        cov[i:i + k, i:i + k] = b
        i += k
    order = np.argsort(z)
    z, v = np.array(z)[order], np.array(v)[order]
    cov = cov[np.ix_(order, order)]
    return z, v, cov


# ------------------------------------------------------------------ growth
def omega_m_of_z(z, p):
    return (p[0] + OMEGA_NU) * (1 + z) ** 3 / (hubble(z, *p) / 100) ** 2


def growth(p, gamma=None):
    """Return (z grid, D normalised to D=a at Z_INIT, f). gamma=None -> exact linear GR growth."""
    lna = np.linspace(-np.log(1 + Z_INIT), 0, 600)
    z = np.exp(-lna) - 1
    if gamma is not None:
        sol = solve_ivp(lambda x, y: [omega_m_of_z(np.exp(-x) - 1, p) ** gamma],
                        [lna[0], 0], [lna[0]], t_eval=lna, rtol=1e-9, atol=1e-12)
        D = np.exp(sol.y[0])
        f = omega_m_of_z(z, p) ** gamma
        return z, D, f

    def dlnH(x):
        e = 1e-4
        return (np.log(hubble(np.exp(-x - e) - 1, *p)) - np.log(hubble(np.exp(-x + e) - 1, *p))) / (2 * e)

    def rhs(x, y):
        zz = np.exp(-x) - 1
        return [y[1], -(2 + dlnH(x)) * y[1] + 1.5 * omega_m_of_z(zz, p) * y[0]]

    a0 = np.exp(lna[0])
    sol = solve_ivp(rhs, [lna[0], 0], [a0, a0], t_eval=lna, rtol=1e-9, atol=1e-12)
    D, dD = sol.y
    return z, D, dD / D


class Model:
    def __init__(self, p, gamma):
        self.p, self.gamma = p, gamma
        zg, D, f = growth(p, gamma)
        self.zg, self.D, self.f = zg[::-1], D[::-1], f[::-1]  # ascending z for interp
        self.D0 = np.interp(0, self.zg, self.D)

    def fs8(self, z, sigma8_0):
        return np.interp(z, self.zg, self.f) * sigma8_0 * np.interp(z, self.zg, self.D) / self.D0

    def sigma8_anchored(self):
        """sigma8 today when the amplitude at Z_INIT is fixed by Planck LCDM."""
        zp, Dp, _ = growth(PLANCK_BG)
        ratio_pl = Dp[-1] / Dp[0]           # D(0)/D(z_init) in Planck LCDM
        ratio = self.D0 / self.D[-1]        # same for this model (D[-1] is z_init)
        return SIG8_PLANCK * ratio / ratio_pl


def chi2_fs8(model, s8, data):
    z, v, icov = data
    r = v - model.fs8(z, s8)
    return r @ icov @ r


def fit_shape(p, gamma, data):
    m = Model(p, gamma)
    r = minimize_scalar(lambda s: chi2_fs8(m, s, data), bounds=(0.4, 1.2), method="bounded")
    return r.fun, r.x


def fit_anchored(p, gamma, data):
    m = Model(p, gamma)
    s_cmb = m.sigma8_anchored()
    rel = SIG8_PLANCK_ERR / SIG8_PLANCK
    f = lambda A: chi2_fs8(m, s_cmb * A, data) + ((A - 1) / rel) ** 2
    r = minimize_scalar(f, bounds=(0.8, 1.2), method="bounded")
    return r.fun, s_cmb * r.x


def profile_gamma(p, data, fitter):
    gs = np.linspace(0.35, 0.95, 61)
    chis = np.array([fitter(p, g, data)[0] for g in gs])
    i = np.argmin(chis)
    best = minimize_scalar(lambda g: fitter(p, g, data)[0], bounds=(gs[max(i - 2, 0)], gs[min(i + 2, 60)]),
                           method="bounded")
    g0, c0 = best.x, best.fun
    lo = brentq(lambda g: fitter(p, g, data)[0] - c0 - 1, 0.2, g0)
    hi = brentq(lambda g: fitter(p, g, data)[0] - c0 - 1, g0, 1.3)
    return gs, chis, g0, c0, g0 - lo, hi - g0


def load_fs8_desi():
    """DESI DR1: ShapeFit f*sigma_s8 in six bins (arXiv:2411.12021, App. A, marginal errors; treated as
    f*sigma8, valid when the fiducial sound horizon is close to the true one) plus the DR1 peculiar-velocity
    consensus f*sigma8(0.07) = 0.450 +- 0.055 (arXiv:2512.03230)."""
    import json
    d = json.load(open(os.path.join(DATA, "desi_dr1_fs8.json")))
    z = np.array([0.07] + [e["z"] for e in d])
    v = np.array([0.450] + [e["fs8"] for e in d])
    err = np.array([0.055] + [e["err"] for e in d])
    return z, v, np.diag(err ** 2)


DATASETS = {"sdss": ("6dFGS, SDSS, BOSS, eBOSS", load_fs8, "fig5_growth_index.png"),
            "desi": ("DESI DR1 full-shape + peculiar velocities", load_fs8_desi, "fig7_growth_index_desi.png")}


def main(which="sdss"):
    os.makedirs(OUT, exist_ok=True)
    dlabel, loader, figname = DATASETS[which]
    z, v, cov = loader()
    data = (z, v, np.linalg.inv(cov))
    print(f"f*sigma8 data ({dlabel}): {len(z)} points at z = {', '.join(f'{x:.3f}' for x in z)}")

    likes = [BAO(), CMB(), SN("DESY5")]
    bgs = {}
    for name in ("LCDM", "Borse d=2"):
        c2, p = fit(likes, MODELS[name])
        bgs[name] = p
        print(f"Background {name:9s}: H0 = {100 * p[2]:.2f}, Om = {(p[0] + OMEGA_NU) / p[2] ** 2:.4f}"
              f"  (fit to DESI DR2 + CMB + DESY5, chi2 = {c2:.2f})")

    profiles = {}
    for label, fitter in (("A. shape only (sigma8 free)", fit_shape),
                          ("B. CMB-anchored amplitude", fit_anchored)):
        print("\n" + "=" * 86 + f"\n{label}\n" + "=" * 86)
        print(f"{'background':10s} {'growth':16s} {'chi2':>7s} {'sigma8_0':>9s} {'S8':>7s}")
        for bname, p in bgs.items():
            om = (p[0] + OMEGA_NU) / p[2] ** 2
            for gname, g in (("GR (exact)", None), ("gamma = 0.55", 0.55), ("gamma = 1/phi", GAMMA_BORSE)):
                c2, s8 = fitter(p, g, data)
                print(f"{bname:10s} {gname:16s} {c2:7.2f} {s8:9.4f} {s8 * np.sqrt(om / 0.3):7.4f}")
            gs, chis, g0, c0, em, ep = profile_gamma(p, data, fitter)
            profiles[(label, bname)] = (gs, chis)
            d_b = fitter(p, GAMMA_BORSE, data)[0] - c0
            d_gr = fitter(p, 0.55, data)[0] - c0
            print(f"{bname:10s} {'gamma free':16s} {c0:7.2f}   best gamma = {g0:.3f} -{em:.3f} +{ep:.3f}"
                  f"   | 1/phi off by {np.sqrt(max(d_b, 0)):.2f} sigma, 0.55 off by {np.sqrt(max(d_gr, 0)):.2f} sigma")
        print(f"{'':10s} {'(dof = ' + str(len(z) - 1) + ' for fixed gamma)':s}")

    # figures
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    zz = np.linspace(0, 1.7, 200)
    ax[0].errorbar(z, v, yerr=np.sqrt(np.diag(cov)), fmt="ko", capsize=3, label=dlabel)
    for bname, gname, g, sty in (("LCDM", "LCDM + GR", None, "k-"), ("Borse d=2", "$\\varphi$-ladder bg + GR", None, "b--"),
                                 ("Borse d=2", "$\\varphi$-ladder bg + $\\gamma=1/\\varphi$", GAMMA_BORSE, "r-")):
        m = Model(bgs[bname], g)
        ax[0].plot(zz, m.fs8(zz, fit_anchored(bgs[bname], g, data)[1]), sty, label=gname + " (CMB-anchored)")
    ax[0].set(xlabel="z", ylabel="$f\\sigma_8$", title="Growth of structure")
    ax[0].legend(fontsize=8)
    for (label, bname), (gs, chis) in profiles.items():
        ax[1].plot(gs, chis - chis.min(), ls="-" if label.startswith("A") else "--",
                   label=f"{label[:2]} " + {"LCDM": "$\\Lambda$CDM bg", "Borse d=2": "$\\varphi$-ladder bg"}[bname])
    ax[1].axvline(GAMMA_BORSE, color="r", ls=":", label="$1/\\varphi$")
    ax[1].axvline(0.55, color="k", ls=":", label="GR 0.55")
    ax[1].axhline(1, color="grey", lw=0.7)
    ax[1].axhline(4, color="grey", lw=0.7)
    ax[1].set(xlabel="$\\gamma$", ylabel="$\\Delta\\chi^2$", ylim=(0, 10), title="Profile in growth index")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, figname), dpi=130)
    print(f"\nSaved figures/{figname}")


if __name__ == "__main__":
    import sys
    main("desi" if "--desi" in sys.argv else "sdss")

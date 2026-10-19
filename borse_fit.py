"""Likelihood fit of the Borse dark-energy prediction against real data.

Models (flat, CPL dark energy w = w0 + wa(1-a)):
  LCDM        w0 = -1,        wa = 0          (3 free: omega_bc, omega_b, h)
  Borse d=2   w0 = -phi/2,    wa = -1/phi     (3 free, DE fixed by theory)
  Borse d=3   w0 = -1+phi^-3, wa = -1/phi     (3 free, alternative branch in the brief)
  w0waCDM     w0, wa free                     (5 free)

Data:
  DESI DR2 BAO (Cobaya bao_data, arXiv:2503.14738)
  CMB distance priors R, l_A, omega_b from Planck 2018 (Chen, Huang & Wang, arXiv:1808.05724)
  Supernovae: Pantheon+, Union3, DESY5 (Cobaya sn_data), absolute magnitude marginalised analytically

Run:  .venv/bin/python borse_fit.py          (best fits, ~1 min)
      .venv/bin/python borse_fit.py --grid   (also w0-wa likelihood contours, ~15 min)
"""

import os
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2_dist, norm
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "figures")

PHI = (1 + 5 ** 0.5) / 2
C_KMS = 299792.458
OMEGA_G = 2.473e-5
# Radiation and neutrinos (CAMB conventions): N_eff = 3.046 with one massive eigenstate (0.06 eV) carrying N_eff/3.
T_CMB = 2.7255
NEFF = 3.046
N_MASSIVE = NEFF / 3
NU_FACTOR = 7 / 8 * (4 / 11) ** (4 / 3)
OMEGA_R = OMEGA_G * (1 + NU_FACTOR * (NEFF - N_MASSIVE))  # photons + massless neutrinos, today
Y_NU0 = 0.06 / (8.617333262e-5 * T_CMB * (4 / 11) ** (1 / 3))  # m_nu / (k T_nu0)
A_NU, P_NU = 180 * 1.2020569031595942 / (7 * np.pi ** 4), 1.83


def omega_nu(z):
    """Massive-neutrino density (in units of rho_crit/h^2) at redshift z: fitting formula of
    Komatsu et al. 2011 (ApJS 192, 18, eq. 26), accurate to ~0.1% from the relativistic to the matter regime."""
    y = Y_NU0 / (1 + z)
    return OMEGA_G * NU_FACTOR * N_MASSIVE * (1 + z) ** 4 * (1 + (A_NU * y) ** P_NU) ** (1 / P_NU)


OMEGA_NU = float(omega_nu(0.0))  # today: behaves as matter

MODELS = {
    "LCDM": (-1.0, 0.0),
    "Borse d=2": (-PHI / 2, -1 / PHI),
    "Borse d=3": (-1 + PHI ** -3, -1 / PHI),
    "w0waCDM": None,
}

# ------------------------------------------------------------------ cosmology
Z_GRID = np.concatenate([np.linspace(0, 0.1, 2001)[:-1], np.linspace(0.1, 2.5, 4801)])


def hubble(z, obc, ob, h, w0, wa):
    ode = h ** 2 - obc - OMEGA_NU - OMEGA_R
    a = 1 / (1 + z)
    de = a ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
    return 100 * np.sqrt(OMEGA_R * (1 + z) ** 4 + obc * (1 + z) ** 3 + omega_nu(z) + ode * de)


def z_star_hs(ob, om):  # Hu & Sugiyama 1996, eq. E-1 (om = total matter density omega_m)
    g1 = 0.0783 * ob ** -0.238 / (1 + 39.5 * ob ** 0.763)
    g2 = 0.560 / (1 + 21.1 * ob ** 1.81)
    return 1048 * (1 + 0.00124 * ob ** -0.738) * (1 + g1 * om ** g2)


_ZCAL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "zstar_calibration.json")


def z_star(ob, obc):
    """Photon-decoupling redshift: Hu & Sugiyama fit times a CAMB-calibrated quadratic correction
    (borse_calibrate_zstar.py). The Chen et al. distance priors use z* from the Planck chains (exact recombination)."""
    zs = z_star_hs(ob, obc + OMEGA_NU)
    if os.path.exists(_ZCAL):
        c = json.load(open(_ZCAL))["coef"]
        xb, xm = ob - 0.0224, obc - 0.143
        zs *= 1 + c[0] + c[1] * xb + c[2] * xm + c[3] * xb ** 2 + c[4] * xm ** 2 + c[5] * xb * xm
    return zs


def r_drag(obc, ob):  # DESI 2024 VI eq. 2.5 (arXiv:2404.03002): omega_m = Omega_m h^2 incl. the 0.06 eV neutrino
    om = obc + OMEGA_NU
    return 147.05 * (ob / 0.02236) ** -0.13 * (om / 0.1432) ** -0.23 * (3.046 / 3.04) ** -0.1


class Cosmo:
    def __init__(self, p):
        self.obc, self.ob, self.h, self.w0, self.wa = p
        H = hubble(Z_GRID, *p)
        self.H = H
        self.DM = cumulative_trapezoid(C_KMS / H, Z_GRID, initial=0)

    def dm(self, z):
        return np.interp(z, Z_GRID, self.DM)

    def dh(self, z):
        return C_KMS / np.interp(z, Z_GRID, self.H)

    def cmb(self):
        p = (self.obc, self.ob, self.h, self.w0, self.wa)
        om = self.obc + OMEGA_NU
        zs = z_star(self.ob, self.obc)
        dm = self.DM[-1] + quad(lambda z: C_KMS / hubble(z, *p), Z_GRID[-1], zs, epsrel=1e-8)[0]
        Rb = lambda z: 3 * self.ob / (4 * OMEGA_G) / (1 + z)
        rs = quad(lambda lz: C_KMS / np.sqrt(3 * (1 + Rb(np.exp(lz) - 1))) / hubble(np.exp(lz) - 1, *p)
                  * np.exp(lz), np.log(1 + zs), np.log(1e8), epsrel=1e-8, limit=200)[0]
        R = np.sqrt(om) * 100 * dm / C_KMS
        return np.array([R, np.pi * dm / rs, self.ob])


# ------------------------------------------------------------------ likelihoods
class BAO:
    def __init__(self):
        rows = [l.split() for l in open(os.path.join(DATA, "desi_dr2_mean.txt")) if not l.startswith("#")]
        self.z = np.array([float(r[0]) for r in rows])
        self.val = np.array([float(r[1]) for r in rows])
        self.kind = [r[2] for r in rows]
        self.icov = np.linalg.inv(np.loadtxt(os.path.join(DATA, "desi_dr2_cov.txt")))
        self.n = len(rows)

    def chi2(self, c):
        rd = r_drag(c.obc, c.ob)
        th = []
        for z, k in zip(self.z, self.kind):
            dm, dh = c.dm(z), c.dh(z)
            th.append({"DM_over_rs": dm, "DH_over_rs": dh, "DV_over_rs": (z * dm ** 2 * dh) ** (1 / 3)}[k] / rd)
        r = self.val - np.array(th)
        return r @ self.icov @ r


class CMB:
    mean = np.array([1.7502, 301.471, 0.02236])
    sig = np.array([0.0046, 0.0895, 0.00015])
    corr = np.array([[1, 0.46, -0.66], [0.46, 1, -0.33], [-0.66, -0.33, 1]])
    icov = np.linalg.inv(corr * np.outer(sig, sig))
    n = 3

    def chi2(self, c):
        r = c.cmb() - self.mean
        return r @ self.icov @ r


class SN:
    """Generic SN likelihood with the absolute magnitude marginalised analytically."""

    INVERSE = {"DESDovekie"}

    def __init__(self, name):
        self.name = name
        cache = os.path.join(DATA, f"{name}.npz")
        if os.path.exists(cache):
            d = np.load(cache)
            self.zhd, self.zhel, self.mu, icov = d["zhd"], d["zhel"], d["mu"], d["icov"]
        else:
            self.zhd, self.zhel, self.mu, mat = getattr(self, "_load_" + name)()
            # DES-Dovekie is distributed as an inverse covariance; the others as covariances
            icov = mat if name in self.INVERSE else np.linalg.inv(mat)
            np.savez(cache, zhd=self.zhd, zhel=self.zhel, mu=self.mu, icov=icov)
        self.icov = icov
        self.e = icov.sum()
        self.colsum = icov.sum(axis=0)
        self.n = len(self.mu)

    @staticmethod
    def _read_cov(path):
        v = np.array(open(path).read().split(), dtype=float)
        n = int(v[0])
        return v[1:].reshape(n, n)

    def _load_PantheonPlus(self):
        d = np.genfromtxt(os.path.join(DATA, "PantheonPlus/Pantheon+SH0ES.dat"), names=True, dtype=None,
                          encoding=None)
        cov = self._read_cov(os.path.join(DATA, "PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov"))
        keep = d["zHD"] > 0.01
        return d["zHD"][keep], d["zHEL"][keep], d["m_b_corr"][keep], cov[np.ix_(keep, keep)]

    def _load_DESY5(self):
        d = np.genfromtxt(os.path.join(DATA, "DESY5/DES-SN5YR_HD.csv"), delimiter=",", names=True)
        cov = self._read_cov(os.path.join(DATA, "DESY5/covsys_000.txt")) + np.diag(d["MUERR_FINAL"] ** 2)
        return d["zHD"], d["zHEL"], d["MU"], cov

    def _load_DESDovekie(self):
        """DES-Dovekie (Popovic et al., arXiv:2511.07517), as packaged in Cobaya sn_data: MU column and the
        upper triangle of the inverse total covariance, unpacked exactly as in cobaya.likelihoods.base_classes.sn."""
        d = np.genfromtxt(os.path.join(DATA, "DES-Dovekie/DES-Dovekie_HD.csv"), delimiter=",", names=True)
        z = np.load(os.path.join(DATA, "DES-Dovekie/covtot_inv_000.npz"))
        n = int(z[z.files[0]][0])
        icov = np.zeros((n, n))
        icov[np.triu_indices(n)] = z[z.files[1]]
        il = np.tril_indices(n, -1)
        icov[il] = icov.T[il]
        return d["zHD"], d["zHEL"], d["MU"], icov

    def _load_Union3(self):
        # header lists one more column than the rows carry; take zcmb, zhel, mb by position
        zcmb, zhel, mb = np.loadtxt(os.path.join(DATA, "Union3/lcparam_full.txt"), usecols=(1, 2, 4), unpack=True)
        cov = self._read_cov(os.path.join(DATA, "Union3/mag_covmat.txt"))
        return zcmb, zhel, mb, cov

    def chi2(self, c):
        mu_th = 5 * np.log10((1 + self.zhel) * c.dm(self.zhd)) + 25
        r = self.mu - mu_th
        icr = self.icov @ r
        return r @ icr - (self.colsum @ r) ** 2 / self.e


# ------------------------------------------------------------------ fitting
def total_chi2(p, likes):
    obc, ob, h, w0, wa = p
    if not (0.08 < obc < 0.2 and 0.015 < ob < 0.03 and 0.5 < h < 0.9 and w0 + wa < 0 and -3 < w0 < 1 and -5 < wa < 3):
        return 1e10
    c = Cosmo(p)
    try:
        return sum(l.chi2(c) for l in likes)
    except (ValueError, FloatingPointError):
        return 1e10


def fit(likes, de=None, quick=False):
    """Minimise chi2. de=(w0, wa) fixes dark energy; None frees it. quick: one start (grids)."""
    if de is None:
        starts = [[0.142, 0.02237, 0.68, -0.8, -0.6], [0.142, 0.02237, 0.66, -0.7, -1.0],
                  [0.142, 0.02237, 0.67, -0.95, -0.2]]
        f = lambda q: total_chi2(q, likes)
    else:
        starts = [[0.142, 0.02237, 0.68], [0.140, 0.02240, 0.66]][:1 if quick else 2]
        f = lambda q: total_chi2([*q, *de], likes)
    best = None
    for s in starts:
        r = minimize(f, s, method="Nelder-Mead",
                     options=dict(xatol=1e-6, fatol=1e-6, maxiter=20000, maxfev=20000, adaptive=True))
        r = minimize(f, r.x, method="Nelder-Mead",
                     options=dict(xatol=1e-7, fatol=1e-7, maxiter=20000, maxfev=20000, adaptive=True))
        if best is None or r.fun < best.fun:
            best = r
    p = list(best.x) if de is None else [*best.x, *de]
    return best.fun, p


def sigma_from_dchi2(d, dof):
    return norm.isf(chi2_dist.sf(max(d, 0), dof) / 2)


def main(do_grid):
    os.makedirs(OUT, exist_ok=True)
    bao, cmb = BAO(), CMB()
    sns = {n: SN(n) for n in ("PantheonPlus", "Union3", "DESY5")}
    combos = {"DESI+CMB": [bao, cmb]}
    combos.update({f"DESI+CMB+{n}": [bao, cmb, s] for n, s in sns.items()})

    summary = {}
    for cname, likes in combos.items():
        N = sum(l.n for l in likes)
        print("\n" + "=" * 92 + f"\n{cname}   (N = {N} data points)\n" + "=" * 92)
        print(f"{'model':11s} {'k':>2s} {'chi2':>9s} {'dchi2 vs LCDM':>14s} {'AIC':>9s} {'BIC':>9s}"
              f" {'H0':>6s} {'Om':>6s} {'w0':>7s} {'wa':>7s}")
        res = {}
        for m, de in MODELS.items():
            k = 5 if de is None else 3
            c2, p = fit(likes, de)
            res[m] = (c2, p, k)
        for m, (c2, p, k) in res.items():
            om = (p[0] + OMEGA_NU) / p[2] ** 2
            print(f"{m:11s} {k:2d} {c2:9.2f} {c2 - res['LCDM'][0]:14.2f} {c2 + 2 * k:9.2f} {c2 + k * np.log(N):9.2f}"
                  f" {100 * p[2]:6.2f} {om:6.4f} {p[3]:7.3f} {p[4]:7.3f}")
        free = res["w0waCDM"][0]
        print("\nLikelihood-ratio tests (2 dof: w0, wa):")
        for m in ("LCDM", "Borse d=2", "Borse d=3"):
            d = res[m][0] - free
            print(f"  {m:10s} is disfavoured relative to the best free (w0, wa) by "
                  f"dchi2 = {d:6.2f}  ->  {sigma_from_dchi2(d, 2):4.2f} sigma")
        summary[cname] = res

    if do_grid:
        ladder_test(combos, summary)
        grid_contours(combos, summary)
    return summary


def ladder_test(combos, summary):
    """Look-elsewhere test: fit every (w0, wa) pair built from the brief's own ladder
    family +-{1, 1/2} x^d (+0, +-1), d = 1..4, and see how special the chosen pair is."""
    x = 1 / PHI
    fam = sorted({round(s * c * x ** d + o, 6) for d in range(1, 5) for c in (1, .5) for s in (1, -1)
                  for o in (0, 1, -1)})
    pairs = [(a, b) for a in fam if -1.3 <= a <= -0.3 for b in fam if -2.5 <= b <= 0.5 and a + b < 0]
    print("\n" + "=" * 92 + f"\nLook-elsewhere: {len(pairs)} (w0, wa) pairs from the brief's ladder family\n" + "=" * 92)
    for cname, likes in combos.items():
        if cname == "DESI+CMB":
            continue
        chis = np.array([fit(likes, pr, quick=True)[0] for pr in pairs])
        lcdm, borse = summary[cname]["LCDM"][0], summary[cname]["Borse d=2"][0]
        order = np.argsort(chis)
        rank = int(np.sum(chis < borse - 1e-6)) + 1
        n_beat = int(np.sum(chis < lcdm - 4))
        print(f"{cname:24s} Borse d=2 ranks {rank:3d}/{len(pairs)};  {n_beat:3d} pairs beat LCDM by dchi2 > 4;"
              f"  best pair w0={pairs[order[0]][0]:.3f}, wa={pairs[order[0]][1]:.3f} (chi2 {chis[order[0]]:.2f})")
        np.save(os.path.join(DATA, f"ladder_{cname}.npy"), chis)


def grid_contours(combos, summary):
    w0s, was = np.linspace(-1.15, -0.45, 16), np.linspace(-2.0, 0.4, 16)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
    for ax, cname in zip(axes, [c for c in combos if c != "DESI+CMB"]):
        likes = combos[cname]
        chi = np.zeros((len(was), len(w0s)))
        for i, wa in enumerate(was):
            for j, w0 in enumerate(w0s):
                chi[i, j] = fit(likes, (w0, wa), quick=True)[0]
            print(f"  {cname}: row {i + 1}/{len(was)}", flush=True)
        dchi = chi - min(chi.min(), summary[cname]["w0waCDM"][0])
        np.save(os.path.join(DATA, f"grid_{cname}.npy"), chi)
        ax.contourf(w0s, was, dchi, levels=[0, 2.30, 6.18], colors=["#3b6fb6", "#a9c4e8"])
        ax.contour(w0s, was, dchi, levels=[2.30, 6.18], colors="k", linewidths=0.8)
        b = summary[cname]["w0waCDM"][1]
        ax.plot(b[3], b[4], "k+", ms=10, label="best fit")
        ax.plot(*MODELS["Borse d=2"], "r*", ms=15, label="$\\varphi$-ladder (d=2)")
        ax.plot(*MODELS["Borse d=3"], "m*", ms=12, label="d=3 branch")
        ax.plot(-1, 0, "kx", ms=9, label="$\\Lambda$CDM")
        ax.set(title=cname, xlabel="$w_0$")
    axes[0].set_ylabel("$w_a$")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Profile likelihood in the $w_0$-$w_a$ plane (68% / 95% regions)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig4_w0wa_likelihood.png"), dpi=130)
    print("Saved figures/fig4_w0wa_likelihood.png")


if __name__ == "__main__":
    main("--grid" in sys.argv)

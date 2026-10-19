"""Version-2 additional tests requested in review.

E1  Joint early-dark-energy + dark-energy check. The hypothesis also predicts f_EDE = x^4/2 = 0.073, while the
    w0-wa fits assume no EDE. (a) EDE-agnostic: BAO (+SN) fits with r_d free (no CMB) - only h*r_d enters, so any
    early-universe change is absorbed. (b) Indicative only: DESI+CMB(+SN) fits with the axion-like EDE fluid
    (n=3, w_n=1/2, f=0.073) included in H(z) for z_c = 1e3, 3.5e3, 1e4; r_d and r_s(z*) rescaled by the EDE/no-EDE
    ratio of the sound-horizon integral. CAVEAT: compressed CMB distance priors assume standard early physics,
    so (b) is not a substitute for a full Boltzmann + CMB likelihood analysis.
B1  DESI DR2 binned w(z) (official chains, model w-binned-3uniform: w piecewise constant in z = 0-0.7, 0.7-1.4,
    1.4-2.1, 2.1-10; w = -1 beyond) versus the bin-averaged phi-ladder w(z).
Q1  Energy exchange between dark matter and dark energy: sign of the apparent phantom crossing and whether a
    constant-w dark energy with Q = -xi H rho_DE or Q = -xi H rho_c (xi > 0: energy flows DE -> DM) reproduces the
    phi-ladder crossing z = 1/sqrt5 and effective wa = -1/phi (background level only).

Run:  .venv/bin/python borse_v2b_tests.py   (output: rerun/v2b_results.txt)
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp
from scipy.optimize import minimize, brentq
from scipy.stats import chi2 as chi2_dist, norm

import borse_fit as BF
from borse_fit import BAO, CMB, SN, MODELS, PHI, OMEGA_NU, OMEGA_R, OMEGA_G, omega_nu, z_star, r_drag, sigma_from_dchi2, Z_GRID, C_KMS
import borse_desi_chains as DC

X = 1 / PHI
LADDER = (-PHI / 2, -1 / PHI)
F_EDE, W_N = X ** 4 / 2, 0.5


def head(t):
    print("\n" + "=" * 96 + f"\n{t}\n" + "=" * 96, flush=True)


# ------------------------------------------------------------------------------------------- E1 helpers
def make_h(p, f=0.0, zc=3500.0):
    obc, ob, h, w0, wa = p
    bg = lambda z: OMEGA_R * (1 + z) ** 4 + obc * (1 + z) ** 3 + omega_nu(z)
    rc = f / (1 - f) * bg(zc) if f else 0.0
    ede = lambda z: 2 * rc / (1 + ((1 + zc) / (1 + z)) ** (3 * (1 + W_N))) if f else 0.0 * z
    ode = h ** 2 - obc - OMEGA_NU - OMEGA_R - (ede(0.0) if f else 0.0)

    def H(z):
        a = 1 / (1 + z)
        de = a ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
        return 100 * np.sqrt(bg(z) + ode * de + ede(z))
    return H


def rs_int(H, z, ob):
    Rb = lambda zz: 3 * ob / (4 * OMEGA_G) / (1 + zz)
    return quad(lambda lz: C_KMS / np.sqrt(3 * (1 + Rb(np.exp(lz) - 1))) / H(np.exp(lz) - 1) * np.exp(lz),
                np.log(1 + z), np.log(1e8), epsrel=1e-7, limit=200)[0]


class CosmoEDE:
    """Background with optional EDE; r_d and r_s(z*) from the calibrated no-EDE values times the EDE ratio."""

    def __init__(self, p, f=0.0, zc=3500.0, rd_free=None):
        self.obc, self.ob, self.h, self.w0, self.wa = p
        self.p = p
        H = make_h(p, f, zc)
        self.Hf = H
        self.Hg = H(Z_GRID)
        self.DM = cumulative_trapezoid(C_KMS / self.Hg, Z_GRID, initial=0)
        self.zs = z_star(self.ob, self.obc)
        if f:
            H0f = make_h(p, 0.0)
            self.ratio_d = rs_int(H, 1060.0, self.ob) / rs_int(H0f, 1060.0, self.ob)
            self.rs_star = rs_int(H, self.zs, self.ob)
        else:
            self.ratio_d = 1.0
            self.rs_star = None
        self.rd = rd_free if rd_free is not None else r_drag(self.obc, self.ob) * self.ratio_d

    def dm(self, z):
        return np.interp(z, Z_GRID, self.DM)

    def dh(self, z):
        return C_KMS / np.interp(z, Z_GRID, self.Hg)

    def cmb(self):
        dm = self.DM[-1] + quad(lambda z: C_KMS / self.Hf(z), Z_GRID[-1], self.zs, epsrel=1e-8)[0]
        rs = self.rs_star if self.rs_star is not None else rs_int(self.Hf, self.zs, self.ob)
        R = np.sqrt(self.obc + OMEGA_NU) * 100 * dm / C_KMS
        return np.array([R, np.pi * dm / rs, self.ob])


def bao_chi2(bao, c):
    th = []
    for z, k in zip(bao.z, bao.kind):
        dm, dh = c.dm(z), c.dh(z)
        th.append({"DM_over_rs": dm, "DH_over_rs": dh, "DV_over_rs": (z * dm ** 2 * dh) ** (1 / 3)}[k] / c.rd)
    r = bao.val - np.array(th)
    return r @ bao.icov @ r


def total(likes, c):
    s = 0.0
    for l in likes:
        s += bao_chi2(l, c) if isinstance(l, BAO) else l.chi2(c)
    return s


def fit_ede(likes, de, f, zc):
    def fun(q):
        obc, ob, h = q[:3]
        w0, wa = (q[3], q[4]) if de is None else de
        if not (0.08 < obc < 0.2 and 0.015 < ob < 0.03 and 0.5 < h < 0.9 and w0 + wa < 0 and -3 < w0 < 1 and -5 < wa < 3):
            return 1e10
        try:
            return total(likes, CosmoEDE([obc, ob, h, w0, wa], f, zc))
        except Exception:
            return 1e10
    starts = ([[0.145, 0.0224, 0.70, -0.8, -0.6], [0.145, 0.0224, 0.69, -0.7, -1.0]] if de is None
              else [[0.145, 0.0224, 0.70]])
    best = None
    for s in starts:
        r = minimize(fun, s, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=6000, maxfev=6000, adaptive=True))
        r = minimize(fun, r.x, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=6000, maxfev=6000, adaptive=True))
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, best.x


def fit_rdfree(likes, de):
    """No CMB: fixed omega_b, h (only Omega_m and h*r_d matter); free omega_bc and r_d (and w0, wa if de is None)."""
    ob, h = 0.02237, 0.6736

    def fun(q):
        obc, rd = q[:2]
        w0, wa = (q[2], q[3]) if de is None else de
        if not (0.05 < obc < 0.3 and 100 < rd < 200 and w0 + wa < 0 and -3 < w0 < 1 and -5 < wa < 3):
            return 1e10
        return total(likes, CosmoEDE([obc, ob, h, w0, wa], rd_free=rd))
    starts = [[0.14, 147.0, -0.8, -0.6], [0.14, 147.0, -0.6, -1.2], [0.14, 147.0, -0.95, -0.1]] if de is None else [[0.14, 147.0], [0.16, 140.0]]
    best = None
    for s in starts:
        r = minimize(fun, s, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-7, maxiter=20000, maxfev=20000, adaptive=True))
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, best.x


def main():
    bao, cmb = BAO(), CMB()
    sns = {n: SN(n) for n in ("PantheonPlus", "Union3", "DESY5", "DESDovekie")}

    # ---------------------------------------------------------------- E1a: r_d free, no CMB
    head("E1a. EDE-agnostic check: DESI BAO (+SN) with r_d free, no CMB (only h*r_d enters)")
    print(f"{'data':22s} {'chi2 phi':>9s} {'chi2 LCDM':>9s} {'chi2 free':>9s} {'phi vs free':>11s} {'LCDM vs free':>12s}  free (w0, wa)    Om(phi)")
    sets = {"DESI BAO": [bao]}
    sets.update({f"BAO+{n}": [bao, s] for n, s in sns.items()})
    for name, likes in sets.items():
        cp, xp = fit_rdfree(likes, LADDER)
        cl, _ = fit_rdfree(likes, (-1.0, 0.0))
        cf, xf = fit_rdfree(likes, None)
        print(f"{name:22s} {cp:9.2f} {cl:9.2f} {cf:9.2f} {sigma_from_dchi2(cp - cf, 2):10.2f}s {sigma_from_dchi2(cl - cf, 2):11.2f}s"
              f"  ({xf[2]:.3f}, {xf[3]:.3f})  {(xp[0] + OMEGA_NU) / 0.6736 ** 2:.3f}")

    # ---------------------------------------------------------------- E1b: EDE included (indicative)
    head("E1b. Indicative: DESI+CMB(+SN) with f_EDE = 0.073 in H(z) (compressed CMB priors assume standard early physics!)")
    combos = {"DESI+CMB": [bao, cmb]}
    combos.update({f"DESI+CMB+{n}": [bao, cmb, s] for n, s in sns.items()})
    print(f"{'z_c':>6s} {'data':24s} {'phi vs free':>11s} {'LCDM vs free':>12s} {'H0(phi)':>8s} {'r_d ratio':>9s}  free (w0, wa)")
    for zc in (1000.0, 3500.0, 10000.0):
        for name, likes in combos.items():
            cp, xp = fit_ede(likes, LADDER, F_EDE, zc)
            cl, _ = fit_ede(likes, (-1.0, 0.0), F_EDE, zc)
            cf, xf = fit_ede(likes, None, F_EDE, zc)
            ratio = CosmoEDE([*xp[:3], *LADDER], F_EDE, zc).ratio_d
            print(f"{zc:6.0f} {name:24s} {sigma_from_dchi2(cp - cf, 2):10.2f}s {sigma_from_dchi2(cl - cf, 2):11.2f}s"
                  f" {100 * xp[2]:8.2f} {ratio:9.4f}  ({xf[3]:.3f}, {xf[4]:.3f})", flush=True)

    # ---------------------------------------------------------------- B1: binned w(z)
    head("B1. DESI DR2 binned w(z) (official chains, w piecewise constant) vs bin-averaged phi-ladder w(z)")
    edges = [(0.0, 0.7), (0.7, 1.4), (1.4, 2.1), (2.1, 10.0)]
    wz = lambda z: LADDER[0] + LADDER[1] * z / (1 + z)
    th = np.array([quad(wz, a, b)[0] / (b - a) for a, b in edges])
    print("theory bin averages (uniform in z): " + ", ".join(f"[{a}-{b}] {t:.3f}" for (a, b), t in zip(edges, th)))
    for folder, label in (("binned_cmb", "DESI+CMB"), ("binned_cmb_pantheonplus", "+Pantheon+"),
                          ("binned_cmb_union3", "+Union3"), ("binned_cmb_desy5sn", "+DESY5")):
        ws, cols = [], []
        for k in range(1, 5):
            path = f"{DC.CHAINS}/{folder}/chain.{k}.txt"
            hdr = open(path).readline().lstrip("#").split()
            d = np.loadtxt(path)
            d = d[int(0.3 * len(d)):]
            ws.append(d[:, hdr.index("weight")])
            cols.append(np.column_stack([d[:, hdr.index(f"w_bin_{i}")] for i in range(4)]))
        w = np.concatenate(ws)
        S = np.vstack(cols)
        mu = np.average(S, axis=0, weights=w)
        cov = np.cov(S.T, aweights=w)
        sd = np.sqrt(np.diag(cov))
        dev = (th - mu) / sd
        d3 = (th - mu)[:3]
        c3 = cov[:3, :3]
        chi3 = d3 @ np.linalg.inv(c3) @ d3
        print(f"{label:12s} " + "  ".join(f"bin{i}: {m:6.3f}+-{s:5.3f} ({x:+.2f}s)" for i, (m, s, x) in enumerate(zip(mu, sd, dev)))
              + f" | first 3 bins jointly: chi2 = {chi3:.2f} (3 dof) -> {norm.isf(chi2_dist.sf(chi3, 3) / 2):.2f} sigma")

    # ---------------------------------------------------------------- Q1: energy exchange
    head("Q1. Energy exchange (background only). Inferred rho_DE,eff = rho_DE + rho_c - rho_c0 a^-3")
    print("Analytic: 1 + w_eff = [(1+w) rho_DE + D] / (rho_DE + D), D = rho_c - rho_c0 a^-3. An apparent phantom past needs D < 0,")
    print("i.e. dark matter LESS dense in the past than a^-3 extrapolation: energy must flow DE -> DM (Q < 0 in rho_c' + 3H rho_c = -Q).")
    Om, w = 0.316, LADDER[0]
    zgrid = np.linspace(0, 3, 601)

    def solve(xi, form):
        def rhs(N, y):
            rc, rde = y
            q = xi * (rde if form == "rho_DE" else rc)  # energy rate DE -> DM, in units of H
            return [-3 * rc + q, -3 * (1 + w) * rde - q]
        N = -np.log(1 + zgrid)
        sol = solve_ivp(rhs, [0, N[-1]], [Om, 1 - Om], t_eval=N, rtol=1e-10, atol=1e-12)
        rc, rde = sol.y
        a = np.exp(N)
        D = rc - Om * a ** -3
        eff = rde + D
        weff = -1 - np.gradient(np.log(np.abs(eff)), N) / 3
        return weff, rde, eff

    def crossing(xi, form):
        weff, rde, eff = solve(xi, form)
        idx = np.where(np.diff(np.sign(weff + 1)))[0]
        return zgrid[idx[0]] if len(idx) else np.inf
    for form in ("rho_DE", "rho_c"):
        print(f"\nQ = -xi H {form} (xi > 0: DE -> DM), dark energy with constant w = {w:.4f}:")
        for xi in (0.0, 0.1, 0.3, 0.5, 1.0):
            weff, rde, eff = solve(xi, form)
            zc_ = crossing(xi, form)
            print(f"   xi = {xi:4.2f}: w_eff(z=0,0.5,1,2) = " + ", ".join(f"{np.interp(z, zgrid, weff):+.3f}" for z in (0, 0.5, 1, 2))
                  + f";  crossing at z = {zc_:.3f}" + ("" if np.all(rde > 0) and np.all(eff > 0) else "  [density sign problem]"))
        try:
            xs = brentq(lambda xi: crossing(xi, form) - 1 / np.sqrt(5), 1e-3, 5.0)
            weff, _, _ = solve(xs, form)
            sel = zgrid <= 1.5
            a = 1 / (1 + zgrid[sel])
            A = np.column_stack([np.ones(sel.sum()), 1 - a])
            w0e, wae = np.linalg.lstsq(A, weff[sel], rcond=None)[0]
            print(f"   -> xi giving crossing at z = 1/sqrt5: {xs:.3f}; effective CPL fit over 0<z<1.5: w0 = {w0e:.3f}, wa = {wae:.3f}"
                  f" (phi-ladder: {LADDER[0]:.3f}, {LADDER[1]:.3f})")
        except ValueError:
            print("   -> no xi in (0, 5] reproduces the crossing at z = 1/sqrt5")
    print("\nGrowth (mu0) for interacting models depends on the perturbation prescription (momentum transfer); not computed here.")


if __name__ == "__main__":
    main()

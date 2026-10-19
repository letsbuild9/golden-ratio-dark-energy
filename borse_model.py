"""Part 3 model: Borse dark energy + pure momentum exchange with dark matter ("dark scattering").

Physics
  Background: flat universe, CPL dark energy with w0 = -phi/2, wa = -1/phi (P1, P2).
  Perturbations: elastic dark matter - dark energy scattering (Simpson 2010, arXiv:1007.1034).
    Only momentum is exchanged, so the background is unchanged. Its action-level form is the
    pure-momentum-transfer ("Type 3") coupling of Pourtsidou, Skordis & Copeland 2013 (arXiv:1307.0458).
  Sub-horizon, quasi-static, smooth dark energy. In e-folds N = ln a:
    dc'' + (2 + dlnH/dN + Drag) dc' = 1.5 (Om_c dc + Om_b db)
    db'' + (2 + dlnH/dN)        db' = 1.5 (Om_c dc + Om_b db)
    Drag(a) = xi (1 + w(a)) Om_DE(a) H(a)/H0
  Dimensionless coupling xi = (sigma_D / m_c) c rho_crit0 / H0.
  Note: for z > z_cross, 1 + w < 0, so the drag reverses sign (anti-friction).

Questions answered
  Q1. What growth index does the brief's own xi = (7 phi - 10)/4 produce?
  Q2. What xi is needed for the growth suppression of gamma = 1/phi (prediction P5)? Is the shift constant in z?
  Q3. What do current f*sigma8 data say about xi?

Run:  .venv/bin/python borse_model.py
"""

import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from borse_fit import BAO, CMB, SN, MODELS, OMEGA_NU, OMEGA_R, PHI, fit, hubble, OUT, sigma_from_dchi2
from borse_growth import load_fs8, growth, PLANCK_BG, SIG8_PLANCK, SIG8_PLANCK_ERR, Z_INIT

X = 1 / PHI
XI_BRIEF = (7 * PHI - 10) / 4
G_SI, MPC_M, C_MS = 6.674e-11, 3.0857e22, 2.99792458e8


class DragModel:
    def __init__(self, p, xi):
        self.p, self.xi = p, xi
        obc, ob, h, w0, wa = p
        self.oc, self.ob = obc - ob + OMEGA_NU, ob  # neutrino mass lumped with CDM
        self.ode = h ** 2 - obc - OMEGA_NU - OMEGA_R
        self.w0, self.wa = w0, wa
        N = np.linspace(-np.log(1 + Z_INIT), 0, 800)
        a0 = np.exp(N[0])
        sol = solve_ivp(self._rhs, [N[0], 0], [a0, a0, a0, a0], t_eval=N, rtol=1e-9, atol=1e-12)
        dc, dcp, db, dbp = sol.y
        a = np.exp(N)
        E2 = (hubble(1 / a - 1, *p) / 100) ** 2
        wc, wb = self.oc * a ** -3 / E2, self.ob * a ** -3 / E2
        self.Om = wc + wb
        dm = (wc * dc + wb * db) / self.Om
        dmp = (wc * dcp + wb * dbp) / self.Om
        self.z = (1 / a - 1)[::-1]           # ascending z
        self.D = dm[::-1]
        self.f = (dmp / dm)[::-1]
        self.Om = self.Om[::-1]
        self.D0 = self.D[0]

    def w(self, a):
        return self.w0 + self.wa * (1 - a)

    def _rhs(self, N, y):
        dc, dcp, db, dbp = y
        a = np.exp(N)
        z = 1 / a - 1
        H = hubble(z, *self.p)
        E2 = (H / 100) ** 2
        e = 1e-4
        dlnH = (np.log(hubble(1 / np.exp(N + e) - 1, *self.p)) - np.log(hubble(1 / np.exp(N - e) - 1, *self.p))) / (2 * e)
        de = a ** (-3 * (1 + self.w0 + self.wa)) * np.exp(-3 * self.wa * (1 - a))
        Om_de = self.ode * de / E2
        drag = self.xi * (1 + self.w(a)) * Om_de * H / (100 * self.p[2])
        src = 1.5 * (self.oc * dc + self.ob * db) * a ** -3 / E2
        return [dcp, -(2 + dlnH + drag) * dcp + src, dbp, -(2 + dlnH) * dbp + src]

    def gamma_eff(self, z):
        return np.interp(z, self.z, np.log(self.f) / np.log(self.Om))

    def fs8(self, z, s8):
        return np.interp(z, self.z, self.f) * s8 * np.interp(z, self.z, self.D) / self.D0

    def sigma8_anchored(self):
        _, Dp, _ = growth(PLANCK_BG)
        return SIG8_PLANCK * (self.D0 / self.D[-1]) / (Dp[-1] / Dp[0])


def xi_physical(xi, h):
    """Convert dimensionless xi to sigma_D/m_c in cm^2/g and barn/GeV."""
    H0 = 100 * h * 1e3 / MPC_M
    rho_c = 3 * H0 ** 2 / (8 * np.pi * G_SI)
    s_m = xi * H0 / (C_MS * rho_c)            # m^2/kg
    return s_m * 10, s_m * 10 / 0.5609         # cm^2/g, b/GeV  (1 b/GeV = 0.5609 cm^2/g)


def suppression(m):
    return m.D0 / m.D[-1]  # D(z=0) / D(z_init)


def growth_ratio(z_asc, f):  # exp( integral f dlna ) from z_init to 0
    N = -np.log(1 + z_asc)
    return np.exp(np.trapezoid(f[::-1], N[::-1]))


def anchored_chi2(m, data):
    z, v, icov = data
    s = m.sigma8_anchored()
    rel = SIG8_PLANCK_ERR / SIG8_PLANCK
    r = minimize_scalar(lambda A: (lambda d: d @ icov @ d)(v - m.fs8(z, s * A)) + ((A - 1) / rel) ** 2,
                        bounds=(0.8, 1.2), method="bounded")
    return r.fun, s * r.x


def main():
    os.makedirs(OUT, exist_ok=True)
    _, p = fit([BAO(), CMB(), SN("DESY5")], MODELS["Borse d=2"])
    h = p[2]
    om0 = (p[0] + OMEGA_NU) / h ** 2
    print(f"Background: Borse d=2 best fit to DESI DR2 + CMB + DESY5: H0 = {100 * h:.2f}, Om = {om0:.4f}")
    zs = [0.0, 0.5, 1.0, 1.5]

    gr = DragModel(p, 0.0)
    print("\nGR (xi = 0): gamma_eff at z = " + ", ".join(f"{z}: {gr.gamma_eff(z):.4f}" for z in zs))

    # Q1: brief's xi
    br = DragModel(p, XI_BRIEF)
    print(f"\nQ1. Brief's xi = (7 phi - 10)/4 = {XI_BRIEF:.4f}")
    print("    gamma_eff at z = " + ", ".join(f"{z}: {br.gamma_eff(z):.4f}" for z in zs))
    print("    Delta gamma    = " + ", ".join(f"{z}: {br.gamma_eff(z) - gr.gamma_eff(z):+.4f}" for z in zs))

    # Q2: can any physical (xi >= 0) drag give the suppression implied by gamma = 1/phi?
    tgt_ratio = suppression(gr) * growth_ratio(gr.z, gr.Om ** X) / growth_ratio(gr.z, gr.f)
    print(f"\nQ2. Target: D(0)/D(z_init) lowered by {tgt_ratio / suppression(gr) - 1:+.2%} (gamma = 1/phi)")
    print("    target Delta gamma = 1/phi - gamma_GR: " + ", ".join(f"z={z}: {X - gr.gamma_eff(z):+.5f}" for z in zs))
    for xi in (0.0, XI_BRIEF, 1.0, 3.0, 10.0):
        m = DragModel(p, xi)
        print(f"    xi = {xi:6.3f}: growth change {suppression(m) / suppression(gr) - 1:+.2%},  Delta gamma"
              + ", ".join(f" z={z}: {m.gamma_eff(z) - gr.gamma_eff(z):+.4f}" for z in zs))
    print("    -> with phantom crossing, positive drag ENHANCES total growth; no xi >= 0 reaches the target.")

    # Q3: data
    z, v, cov = load_fs8()
    data = (z, v, np.linalg.inv(cov))
    xis = np.linspace(0, 6, 25)
    chis = np.array([anchored_chi2(DragModel(p, xi), data)[0] for xi in xis])
    i = int(np.argmin(chis))
    print(f"\nQ3. f*sigma8 (7 points, CMB-anchored), crossing background:")
    for name, xi in (("xi = 0 (GR)", 0.0), ("xi = brief", XI_BRIEF), ("xi = 3", 3.0)):
        c, s8 = anchored_chi2(DragModel(p, xi), data)
        print(f"    {name:12s} chi2 = {c:5.2f}   sigma8 = {s8:.4f}   S8 = {s8 * np.sqrt(om0 / 0.3):.4f}")
    above = (chis - chis.min() > 3.84) & (xis > xis[i])
    print(f"    best xi on grid = {xis[i]:.2f};  95% upper limit xi < {xis[np.argmax(above)] if above.any() else '>6'}")

    # Q4: non-crossing alternative with the same w0. A CPL history never crosses -1 iff w0 + wa >= -1,
    # i.e. wa >= -(1 + w0) = -F(2); the boundary wa = -(1 + w0) is the non-crossing case closest to the data.
    wa_th = -(1 + p[3])
    print(f"\nQ4. Non-crossing alternative: w0 = -phi/2 with wa >= {wa_th:.3f} (boundary point used), w >= -1 always")
    bao, cmb = BAO(), CMB()
    for n in ("PantheonPlus", "Union3", "DESY5"):
        L = [bao, cmb, SN(n)]
        cf, ct = fit(L, None)[0], fit(L, (p[3], wa_th))[0]
        print(f"    DESI+CMB+{n:12s}: non-crossing boundary point disfavoured at {sigma_from_dchi2(ct - cf, 2):.2f} sigma")
    _, pt = fit([bao, cmb, SN("DESY5")], (p[3], wa_th))
    gt = Clipped(pt, 0.0)
    tgt_t = suppression(gt) * growth_ratio(gt.z, gt.Om ** X) / growth_ratio(gt.z, gt.f)
    xi_t = brentq(lambda xi: suppression(Clipped(pt, xi)) - tgt_t, 0.0, 200.0)
    mt = Clipped(pt, xi_t)
    print(f"    xi needed for gamma = 1/phi suppression: {xi_t:.2f} ({xi_t / XI_BRIEF:.0f}x the brief's),"
          f" sigma_D/m_c ~ {xi_physical(xi_t, pt[2])[1]:.0f} b/GeV")
    print("    Delta gamma: " + ", ".join(f"z={z}: {mt.gamma_eff(z) - gt.gamma_eff(z):+.5f}" for z in zs)
          + "   (target 1/phi - gamma_GR: " + ", ".join(f"{X - gt.gamma_eff(z):+.5f}" for z in zs) + ")")

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    zz = np.linspace(0, 2, 200)
    for m, ref, lab, sty in ((br, gr, f"crossing bg, $\\xi=(7\\varphi-10)/4$", "b--"),
                             (DragModel(p, 3.0), gr, "crossing bg, $\\xi$=3", "b-"),
                             (mt, gt, f"non-crossing bg, $\\xi$={xi_t:.1f}", "g-")):
        ax[0].plot(zz, m.gamma_eff(zz) - ref.gamma_eff(zz), sty, label=lab)
    ax[0].plot(zz, X - gr.gamma_eff(zz), "r:", label="target $1/\\varphi - \\gamma_{\\rm GR}$")
    ax[0].axhline(0, color="k", lw=0.7)
    ax[0].axvline(1 / 5 ** 0.5, color="grey", lw=0.7, ls=":")
    ax[0].set(xlabel="z", ylabel="$\\Delta\\gamma$ vs GR", title="Growth-index shift from momentum drag")
    ax[0].legend(fontsize=8)
    ax[1].plot(xis, chis - chis.min(), "k-")
    ax[1].axvline(XI_BRIEF, color="b", ls="--", label="$\\xi=(7\\varphi-10)/4$")
    ax[1].axhline(3.84, color="grey", lw=0.7)
    ax[1].set(xlabel="$\\xi$", ylabel="$\\Delta\\chi^2$", title="$f\\sigma_8$ data vs drag (crossing bg)", ylim=(0, 10))
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig6_drag_model.png"), dpi=130)
    print("\nSaved figures/fig6_drag_model.png")


class Clipped(DragModel):
    """Same model, but w(a) never below -1 (non-crossing)."""

    def w(self, a):
        return np.maximum(-1.0, self.w0 + self.wa * (1 - a))


if __name__ == "__main__":
    main()

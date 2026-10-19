"""Numerical checks of the cosmology claims in the Borse framework.

Run:  .venv/bin/python borse_checks.py
Writes figures to ./figures and prints a summary.

Checks
  1. Late dark energy: w0 = -phi/2, wa = -1/phi  -> phantom crossing z, distance to DESI DR2.
  2. Growth index gamma = 1/phi -> S8 relative to Planck LCDM.
  3. Early dark energy fraction f = x^4/2 -> sound horizon r_s and H0 (fluid EDE model).

Simplifications (deliberate, so every number is traceable):
  - Flat universe, Planck 2018 base parameters, massless neutrinos at early times
    (0.06 eV neutrino mass folded into late-time matter).
  - Recombination / drag redshifts held fixed at Planck values when EDE is added.
  - EDE is the standard fluid approximation (Poulin et al. 2018): constant density
    before z_c, decaying as a^{-3(1+w_n)} after, with w_n = 1/2 (n = 3 potential).
  - A real test needs a Boltzmann code (CLASS/CAMB variant with EDE) plus an MCMC fit
    to Planck/ACT + BAO + SNe; this script only tells you whether the claim is in range.
"""

import os
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHI = (1 + 5 ** 0.5) / 2
X = 1 / PHI
C_KMS = 299792.458
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

# Planck 2018 TT,TE,EE+lowE+lensing (base LCDM)
PL = dict(
    omega_b=0.02237, omega_c=0.1200, h=0.6736, z_star=1089.92, z_drag=1059.94,
    r_star=144.43, r_drag=147.09, theta100=1.04110, S8=0.832,
)
from borse_fit import OMEGA_NU, OMEGA_G, OMEGA_R, omega_nu  # same radiation / massive-neutrino treatment as the fits

# DESI DR2 BAO + CMB + SNe (arXiv:2503.14738); errors symmetrised
DESI = {
    "DESI+CMB+Pantheon+": (-0.838, 0.055, -0.62, 0.205),
    "DESI+CMB+Union3": (-0.667, 0.088, -1.09, 0.29),
    "DESI+CMB+DESY5": (-0.752, 0.057, -0.86, 0.215),
}


def section(title):
    print("\n" + "=" * 72 + f"\n{title}\n" + "=" * 72)


# ---------------------------------------------------------------- 1. late DE
def check_dark_energy():
    section("1. Late dark energy (CPL: w = w0 + wa z/(1+z))")
    w0, wa = -PHI / 2, -1 / PHI
    w = lambda z: w0 + wa * z / (1 + z)
    z_cross = brentq(lambda z: w(z) + 1, 1e-6, 50)
    print(f"Borse w0 = {w0:.4f}, wa = {wa:.4f}")
    print(f"Phantom crossing z = {z_cross:.4f}   (claim 1/sqrt5 = {5 ** -0.5:.4f})")

    print("\nDistance from DESI DR2 best fits (Gaussian approx, w0-wa correlation unknown):")
    print(f"{'dataset':22s} {'dw0/sig':>8s} {'dwa/sig':>8s} {'nsig rho=0':>11s} {'nsig rho=-0.9':>14s}")
    for name, (m0, s0, ma, sa) in DESI.items():
        d0, da = (w0 - m0) / s0, (wa - ma) / sa
        out = []
        for rho in (0.0, -0.9):
            chi2 = (d0 ** 2 + da ** 2 - 2 * rho * d0 * da) / (1 - rho ** 2)
            out.append(np.sqrt(chi2))  # Mahalanobis distance, 2 dof
        print(f"{name:22s} {d0:8.2f} {da:8.2f} {out[0]:11.2f} {out[1]:14.2f}")

    z = np.linspace(0, 3, 400)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(z, w(z), lw=2, label="$\\varphi$-ladder: $w_0=-\\varphi/2$, $w_a=-1/\\varphi$")
    colors = dict(zip(DESI, ("C1", "C2", "C3")))
    for name, (m0, _, ma, _) in DESI.items():
        ax[0].plot(z, m0 + ma * z / (1 + z), ls="--", lw=1, color=colors[name], label=name)
    ax[0].axhline(-1, color="k", lw=0.8)
    ax[0].axvline(z_cross, color="grey", ls=":", lw=1)
    ax[0].set(xlabel="z", ylabel="w(z)", title=f"Phantom crossing at z = {z_cross:.4f}")
    ax[0].legend(fontsize=8)

    for name, (m0, s0, ma, sa) in DESI.items():
        ax[1].errorbar(m0, ma, xerr=s0, yerr=sa, fmt="o", capsize=3, color=colors[name], label=name)
    ax[1].plot(w0, wa, "r*", ms=16, label="$\\varphi$-ladder")
    ax[1].plot(-1, 0, "kx", ms=10, label="$\\Lambda$CDM")
    ax[1].set(xlabel="$w_0$", ylabel="$w_a$", title="vs DESI DR2 (1$\\sigma$, uncorrelated bars)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_dark_energy.png"), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- 2. growth
def growth_today(Om, gamma=None):
    """D(a=1) normalised to D = a deep in matter domination."""
    Om_a = lambda a: Om * a ** -3 / (Om * a ** -3 + 1 - Om)
    ai = 1e-3
    if gamma is not None:  # dlnD/dlna = Omega_m(a)^gamma
        sol = solve_ivp(lambda lna, y: [Om_a(np.exp(lna)) ** gamma], [np.log(ai), 0], [np.log(ai)],
                        rtol=1e-9, atol=1e-12)
        return np.exp(sol.y[0, -1])
    # exact LCDM linear growth: D'' + (2 + dlnH/dlna) D' = 1.5 Om(a) D, primes = d/dlna
    def rhs(lna, y):
        a = np.exp(lna)
        dlnH = -1.5 * Om_a(a)
        return [y[1], -(2 + dlnH) * y[1] + 1.5 * Om_a(a) * y[0]]
    sol = solve_ivp(rhs, [np.log(ai), 0], [ai, ai], rtol=1e-9, atol=1e-12)
    return sol.y[0, -1]


def check_growth():
    section("2. Growth index gamma = 1/phi -> S8")
    Om = (PL["omega_b"] + PL["omega_c"] + OMEGA_NU) / PL["h"] ** 2
    D_exact = growth_today(Om)
    rows = []
    for g in (0.55, X):
        ratio = growth_today(Om, g) / D_exact
        rows.append((g, ratio, PL["S8"] * ratio))
    print(f"Omega_m = {Om:.4f}, Planck S8 (early-universe normalised) = {PL['S8']}")
    for g, r, s8 in rows:
        print(f"gamma = {g:.4f}:  D/D_LCDM = {r:.4f}  ->  S8 = {s8:.3f}")
    print("Claim: S8 ~ 0.778.  Weak lensing: DES Y3 ~0.776, KiDS-1000 ~0.759, KiDS-Legacy ~0.815.")
    print("Note: gamma alone ignores the w0/wa background, which also changes growth.")

    gs = np.linspace(0.50, 0.70, 60)
    s8s = [PL["S8"] * growth_today(Om, g) / D_exact for g in gs]
    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.plot(gs, s8s, lw=2)
    ax.axvline(X, color="r", ls=":", label=f"$\\gamma=1/\\phi$ -> S8 = {rows[1][2]:.3f}")
    ax.axhspan(0.776 - 0.017, 0.776 + 0.017, alpha=0.2, label="DES Y3 S8")
    ax.set(xlabel="growth index $\\gamma$", ylabel="$S_8$", title="S8 vs growth index (Planck normalisation)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_growth.png"), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- 3. EDE
def make_H(h, omega_m, f_ede=0.0, z_c=3500.0, w_n=0.5):
    omega_lam = h ** 2 - omega_m - OMEGA_R
    rho_bg = lambda z: OMEGA_R * (1 + z) ** 4 + (omega_m - OMEGA_NU) * (1 + z) ** 3 + omega_nu(z)
    rho_c = f_ede / (1 - f_ede) * rho_bg(z_c)  # EDE = fraction f of total at z_c

    def H(z):
        ede = 2 * rho_c / (1 + ((1 + z_c) / (1 + z)) ** (3 * (1 + w_n))) if f_ede else 0.0
        return 100 * np.sqrt(rho_bg(z) + omega_lam + ede)
    return H


def sound_horizon(H, z, omega_b):
    R = lambda zz: 3 * omega_b / (4 * OMEGA_G) / (1 + zz)
    integrand = lambda lnz1: (C_KMS / np.sqrt(3 * (1 + R(np.exp(lnz1) - 1)))
                              / H(np.exp(lnz1) - 1) * np.exp(lnz1))
    return quad(integrand, np.log(1 + z), np.log(1e8), limit=400, epsabs=0, epsrel=1e-8)[0]


def comoving_distance(H, z):
    return quad(lambda zz: C_KMS / H(zz), 0, z, limit=400, epsrel=1e-8)[0]


def check_ede():
    section("3. Early dark energy f = x^4/2 -> sound horizon and H0")
    ob, om = PL["omega_b"], PL["omega_b"] + PL["omega_c"] + OMEGA_NU
    f = X ** 4 / 2
    H0f = make_H(PL["h"], om)
    rs0, rd0 = sound_horizon(H0f, PL["z_star"], ob), sound_horizon(H0f, PL["z_drag"], ob)
    theta = rs0 / comoving_distance(H0f, PL["z_star"])
    print("Baseline LCDM (validation of this code against Planck 2018):")
    print(f"  r_s(z*) = {rs0:.2f} Mpc (Planck {PL['r_star']}), r_d = {rd0:.2f} Mpc (Planck {PL['r_drag']})")
    print(f"  100 theta* = {100 * theta:.5f} (Planck {PL['theta100']})")

    def h_from_theta(rs):  # keep theta*, omega_b, omega_c fixed; solve for h
        return brentq(lambda h: rs / comoving_distance(make_H(h, om), PL["z_star"]) - theta, 0.5, 0.95)

    print(f"\nEDE fraction f = x^4/2 = {f:.5f}")
    print(f"Claim: r_s = 136.36 Mpc, H0 = 72.70.  Note 147.09*(1-f) = {PL['r_drag'] * (1 - f):.2f},"
          f" 67.36/(1-f) = {67.36 / (1 - f):.2f}  -> the claim assumes r_s shrinks by exactly f.")
    print(f"\n{'z_c':>7s} {'r_d [Mpc]':>10s} {'dr_d/r_d':>9s} {'H0 (theta* fixed)':>18s} {'H0 (H0*r_d fixed)':>18s}")
    zcs = [1000, 2000, 3500, 5000, 10000]
    for zc in zcs:
        Hf = make_H(PL["h"], om, f, zc)
        rs, rd = sound_horizon(Hf, PL["z_star"], ob), sound_horizon(Hf, PL["z_drag"], ob)
        print(f"{zc:7d} {rd:10.2f} {rd / rd0 - 1:9.2%} {100 * h_from_theta(rs):18.2f} {67.36 * rd0 / rd:18.2f}")

    f_need = brentq(lambda ff: sound_horizon(make_H(PL["h"], om, ff, 3500), PL["z_drag"], ob) - 136.36,
                    1e-4, 0.6)
    print(f"\nEDE fraction needed at z_c=3500 to reach r_d = 136.36 Mpc: f = {f_need:.3f}")
    print("Caveat: in full CMB fits EDE also raises omega_c (~+10%), which pulls H0 back down.")

    fs = np.linspace(0, 0.25, 40)
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for zc in (1000, 3500, 10000):
        rds = [sound_horizon(make_H(PL["h"], om, ff, zc), PL["z_drag"], ob) for ff in fs]
        ax.plot(fs, rds, label=f"fluid EDE, $z_c$={zc}")
    ax.plot(fs, rd0 * (1 - fs), "k--", label="earlier claim $r_d(1-f)$")
    ax.axvline(f, color="r", ls=":", label=f"$f=x^4/2$={f:.3f}")
    ax.axhline(136.36, color="grey", lw=0.8)
    ax.set(xlabel="EDE fraction f at $z_c$", ylabel="$r_d$ [Mpc]", title="Sound horizon vs EDE fraction")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig3_ede.png"), dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- 4. brief
def check_brief_consistency():
    section("4. Internal consistency of the Scientific Brief (pages 4-8)")
    f = X ** 4 / 2
    print(f"x^4/2 = {f:.5f},  (5 - 3 phi)/2 = {(5 - 3 * PHI) / 2:.5f}   (identity holds)")
    print(f"r_s formula: 147.09 * (1 - phi^-4/2) = {147.09 * (1 - f):.2f} Mpc   (brief prints 147.09: typo)")
    print(f"H0 formula: 67.40 / (1 - phi^-4/2) = {67.40 / (1 - f):.2f}   (this is H0*r_s = const, not theta* fixed)")
    print(f"xi = (7 phi - 10)/4 = {(7 * PHI - 10) / 4:.4f}   (new constant: 7, 10, 4 are not from the d-ladder)")
    print(f"d=3 branch: w0 = -1 + phi^-3 = {-1 + X ** 3:.4f}  vs d=2 branch: w0 = {-PHI / 2:.4f}  -> two predictions for one number")

    print("\nPhantom-crossing redshift of each DESI DR2 best fit (CPL):")
    for name, (m0, _, ma, _) in DESI.items():
        u = (-1 - m0) / ma
        print(f"  {name:22s} z_cross = {u / (1 - u):.3f}")
    print(f"  Borse                  z_cross = {5 ** -0.5:.3f}   (fixed by w0, wa: not an independent test)")

    print("\nBrief's 'statistical agreement' column, recomputed:")
    rows = [
        ("w0 vs Pantheon+", -PHI / 2, -0.838, 0.055),
        ("wa vs Pantheon+", -X, -0.62, 0.205),
        ("gamma vs Nguyen+2023 RSD", X, 0.633, 0.025),
        ("S8 (0.778 claimed) vs DES Y3", 0.778, 0.776, 0.017),
        ("S8 (0.812 from gamma) vs DES Y3", 0.812, 0.776, 0.017),
        ("S8 (0.812 from gamma) vs KiDS-Legacy", 0.812, 0.815, 0.021),
        ("H0 vs SH0ES 2022", 67.40 / (1 - f), 73.04, 1.04),
    ]
    for name, pred, obs, sig in rows:
        print(f"  {name:38s} pred {pred:8.4f}  obs {obs:7.3f} +- {sig:5.3f}  -> {abs(pred - obs) / sig:4.2f} sigma")


# ---------------------------------------------------------------- 5. look-elsewhere
def check_numerology():
    section("5. Look-elsewhere test: how easy is it to hit these numbers with phi?")
    # Family the brief itself uses: (p + q*phi)/r, e.g. (5-3phi)/2 and (7phi-10)/4
    vals = {round((p + q * PHI) / r, 10)
            for p in range(-10, 11) for q in range(-10, 11) for r in (1, 2, 3, 4) if q != 0}
    vals = np.array(sorted(vals))
    ladder = {round(s * c * X ** d + o, 10)
              for d in range(1, 5) for c in (1, 0.5) for s in (1, -1) for o in (0, 1, -1)}
    ladder = np.array(sorted(ladder))
    targets = [
        ("w0 (Pantheon+)", -0.838, 0.055), ("w0 (DESY5)", -0.752, 0.057), ("w0 (Union3)", -0.667, 0.088),
        ("wa (Pantheon+)", -0.62, 0.205), ("gamma (RSD)", 0.633, 0.025),
        ("S8/S8_Planck (DES Y3)", 0.776 / 0.832, 0.017 / 0.832), ("EDE f (typical fits)", 0.09, 0.03),
    ]
    print(f"Family A: (p + q phi)/r, |p|,|q| <= 10, r <= 4  -> {len(vals)} distinct numbers")
    print(f"Family B: brief's ladder +-{{1, 1/2}} x^d (+0, +-1), d = 1..4 -> {len(ladder)} distinct numbers")
    print(f"{'target':24s} {'1-sigma window':>18s} {'hits A':>7s} {'hits B':>7s}  example from A")
    for name, mu, sig in targets:
        hits = vals[np.abs(vals - mu) <= sig]
        hb = ladder[np.abs(ladder - mu) <= sig]
        ex = min(hits, key=lambda v: abs(v - mu)) if len(hits) else float("nan")
        print(f"{name:24s} {mu - sig:8.3f} .. {mu + sig:6.3f} {len(hits):7d} {len(hb):7d}  {ex:.4f}")
    print("Many hits => matching a measurement with some phi-expression is expected by chance,")
    print("so a match is evidence only if the expression is fixed BEFORE looking at the data.")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    check_dark_energy()
    check_growth()
    check_ede()
    check_brief_consistency()
    check_numerology()
    print(f"\nFigures written to {OUT}")

"""Energy-exchange detail (Open problem 4): constant w = -phi/2 dark energy with Q = -xi H rho_DE or -xi H rho_c
(xi > 0: energy flows DE -> DM). Finds xi placing the apparent crossing at z = 1/sqrt5 and compares the inferred
w_eff(z) with the phi-ladder w(z); reports where the inferred dark-energy density turns non-positive.

Run:  .venv/bin/python borse_energy_exchange.py   (output: rerun/energy_exchange_detail.txt)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PHI = (1 + 5 ** 0.5) / 2
w, Om = -PHI / 2, 0.316
zg = np.linspace(0, 3, 3001)


def solve(xi, form):
    def rhs(N, y):
        rc, rde = y
        q = xi * (rde if form == "rho_DE" else rc)
        return [-3 * rc + q, -3 * (1 + w) * rde - q]
    N = -np.log(1 + zg)
    s = solve_ivp(rhs, [0, N[-1]], [Om, 1 - Om], t_eval=N, rtol=1e-11, atol=1e-13)
    rc, rde = s.y
    a = np.exp(N)
    eff = rde + rc - Om * a ** -3
    weff = -1 - np.gradient(np.log(np.abs(eff)), N) / 3
    return weff, rde, eff


def zcross(xi, form):
    weff, _, _ = solve(xi, form)
    i = np.where(np.diff(np.sign(weff + 1)))[0]
    return zg[i[0]] if len(i) else np.inf


wl = lambda z: -PHI / 2 - (1 / PHI) * z / (1 + z)
for form in ("rho_DE", "rho_c"):
    xi = brentq(lambda x: zcross(x, form) - 1 / 5 ** .5, 1e-3, 5)
    weff, rde, eff = solve(xi, form)
    neg = zg[np.where(eff <= 0)[0][0]] if np.any(eff <= 0) else None
    print(f"Q = -xi H {form}: xi = {xi:.3f} puts the crossing at z = 1/sqrt5. w_eff vs phi-ladder:")
    for z in (0, 0.3, 0.5, 0.7, 1.0):
        print(f"   z={z:3.1f}: interacting w_eff = {np.interp(z, zg, weff):+.3f}   phi-ladder w = {wl(z):+.3f}")
    print(f"   inferred dark-energy density becomes non-positive at z = {neg:.2f}" if neg is not None else "   inferred density stays positive to z=3")

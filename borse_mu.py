"""Effective-gravity mapping and growth extras on the phi-ladder background (DESI DR2 + CMB + DESY5 fit).

- gamma_GR(z) from exact GR growth, and the predicted gamma for both growth rules;
- mu(a) implied by gamma = gamma_GR + F(4) and by gamma = 1/phi, fitted to mu(a) = 1 + mu0 Omega_DE(a)/Omega_DE,0;
- CMB-anchored f*sigma8 chi2 (SDSS-era data) for GR, gamma = 1/phi and gamma_GR + F(4).

Run:  .venv/bin/python borse_mu.py   (output: rerun/mu_and_extras.txt)
"""

import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.integrate import cumulative_trapezoid

from borse_model import DragModel, F4, fit, BAO, CMB, SN, MODELS, hubble, OMEGA_NU, OMEGA_R
import borse_growth as G

_, p = fit([BAO(), CMB(), SN("DESY5")], MODELS["Borse d=2"])
X = G.GAMMA_BORSE
gr = DragModel(p, 0.0)
print("gamma_GR exact at z=0,0.5,1:", [f"{float(gr.gamma_eff(z)):.6f}" for z in (0, 0.5, 1)],
      "-> predicted gamma (GR+F4):", [f"{float(gr.gamma_eff(z)) + F4:.6f}" for z in (0, 0.5, 1)])
z = gr.z
N = -np.log(1 + z)
Om = gr.Om
e = 1e-4
dlnH = (np.log(hubble((1 + z) * np.exp(-e) - 1, *p)) - np.log(hubble((1 + z) * np.exp(e) - 1, *p))) / (2 * e)
mug = (np.gradient(gr.f, N) + gr.f ** 2 + (2 + dlnH) * gr.f) / (1.5 * Om)
ode = p[2] ** 2 - p[0] - OMEGA_NU - OMEGA_R
a = 1 / (1 + z)
Ode = ode * a ** (-3 * (1 + p[3] + p[4])) * np.exp(-3 * p[4] * (1 - a)) / (hubble(z, *p) / 100) ** 2
sel = (z < 3) & (z > 0.02)
for lab, f in (("GR+F(4)", gr.f * Om ** F4), ("gamma=1/phi", Om ** X)):
    mu = (np.gradient(f, N) + f ** 2 + (2 + dlnH) * f) / (1.5 * Om) / mug
    r = mu[sel] - 1
    mu0 = least_squares(lambda m: r - m[0] * Ode[sel] / Ode[0], [0.0]).x[0]
    print(f"[{lab}] mu-1 at z=0,0.5,1,2: [{', '.join(f'{mu[np.argmin(abs(z - zz))] - 1:.3f}' for zz in (0, 0.5, 1, 2))}]")
    print(f"[{lab}] mu0 = {mu0:+.3f}; max residual = {np.max(abs(r - mu0 * Ode[sel] / Ode[0])):.4f}")

zz, v, cov = G.load_fs8()
data = (zz, v, np.linalg.inv(cov))


def shifted(p, delta):
    m = G.Model(p, None)
    om = G.omega_m_of_z(m.zg, p)
    fn = m.f * om ** delta
    L = -np.log(1 + m.zg[::-1])
    lnD = np.log(m.D[-1]) + cumulative_trapezoid(fn[::-1], L, initial=0)
    m.f, m.D = fn, np.exp(lnD[::-1])
    m.D0 = np.interp(0, m.zg, m.D)
    return m


def anch(m):
    s = m.sigma8_anchored()
    rel = G.SIG8_PLANCK_ERR / G.SIG8_PLANCK
    return minimize_scalar(lambda A: G.chi2_fs8(m, s * A, data) + ((A - 1) / rel) ** 2, bounds=(.8, 1.2), method="bounded").fun


print(f"SDSS anchored chi2: GR {anch(G.Model(p, None)):.2f}  1/phi {anch(G.Model(p, X)):.2f}  GR+F4 {anch(shifted(p, F4)):.2f}")

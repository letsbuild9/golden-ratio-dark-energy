"""Calibrate the photon-decoupling redshift z* against CAMB (exact recombination).

The Chen, Huang & Wang (2019) distance priors were built with z* from the Planck chains, i.e. from a full
recombination calculation. The Hu & Sugiyama (1996) fitting formula is ~0.18% high, which shifts l_A by
several of its standard deviations. This script fits a quadratic correction z*_CAMB / z*_HS - 1 in
(omega_b, omega_cb) on a grid covering all fits in this project, and saves it to data/zstar_calibration.json.

Run:  .venv/bin/python borse_calibrate_zstar.py
"""

import json
import os
import numpy as np
import camb

from borse_fit import z_star_hs, OMEGA_NU

rows = []
for ob in np.linspace(0.0210, 0.0238, 8):
    for ocb in np.linspace(0.128, 0.158, 8):
        P = camb.CAMBparams()
        P.set_cosmology(H0=67.0, ombh2=ob, omch2=ocb - ob, mnu=0.06, num_massive_neutrinos=1, nnu=3.046, TCMB=2.7255)
        zc = camb.get_background(P).get_derived_params()["zstar"]
        rows.append((ob, ocb, zc, z_star_hs(ob, ocb + OMEGA_NU)))
rows = np.array(rows)
xb, xm = rows[:, 0] - 0.0224, rows[:, 1] - 0.143
A = np.column_stack([np.ones(len(rows)), xb, xm, xb ** 2, xm ** 2, xb * xm])
ratio = rows[:, 2] / rows[:, 3] - 1
coef, *_ = np.linalg.lstsq(A, ratio, rcond=None)
resid = ratio - A @ coef
out = {"coef": coef.tolist(), "grid": {"omega_b": [0.0210, 0.0238], "omega_cb": [0.128, 0.158]},
       "max_abs_residual_relative": float(np.abs(resid).max()), "camb_version": camb.__version__}
json.dump(out, open(os.path.join("data", "zstar_calibration.json"), "w"), indent=1)
print(f"z*_CAMB/z*_HS - 1 ranges {ratio.min():.5f} .. {ratio.max():.5f}; quadratic correction residual max "
      f"{np.abs(resid).max():.2e} (relative), i.e. < {np.abs(resid).max() * 1090:.4f} in z*")
print("saved data/zstar_calibration.json")

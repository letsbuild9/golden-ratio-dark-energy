"""T7: cross-check our background and growth code against CAMB (PPF dark energy) for the phi-ladder best fit.

Compares H(z), D_M(z), the drag-epoch sound horizon r_d, the CMB shift parameter R and acoustic scale l_A,
and the shape of f*sigma8(z) (GR growth) between our pipeline and CAMB.

Run:  .venv/bin/python borse_camb_check.py   (output: rerun/camb_check.txt)
"""

import numpy as np
import camb

from borse_fit import BAO, CMB, SN, MODELS, OMEGA_NU, Cosmo, fit, hubble, r_drag, z_star
from borse_model import DragModel

likes = [BAO(), CMB(), SN("DESY5")]
_, p = fit(likes, MODELS["Borse d=2"])
obc, ob, h, w0, wa = p
print(f"phi-ladder best fit (DESI+CMB+DESY5): omega_bc={obc:.5f}, omega_b={ob:.5f}, h={h:.4f}, w0={w0:.4f}, wa={wa:.4f}")

pars = camb.CAMBparams()
pars.set_cosmology(H0=100 * h, ombh2=ob, omch2=obc - ob, mnu=0.06, num_massive_neutrinos=1, nnu=3.046, TCMB=2.7255)
pars.set_dark_energy(w=w0, wa=wa, dark_energy_model="ppf")
pars.InitPower.set_params(As=2.1e-9, ns=0.965)
zs = [0.0, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
pars.set_matter_power(redshifts=zs[::-1], kmax=2.0)
res = camb.get_results(pars)
der = res.get_derived_params()

c = Cosmo(p)
print(f"\n{'quantity':22s} {'ours':>12s} {'CAMB':>12s} {'diff':>9s}")
for z in (0.3, 0.5, 1.0, 2.0):
    ours_H, camb_H = hubble(z, *p), res.hubble_parameter(z)
    ours_D, camb_D = c.dm(z), res.comoving_radial_distance(z)
    print(f"{'H(z=' + str(z) + ')':22s} {ours_H:12.4f} {camb_H:12.4f} {100 * (ours_H / camb_H - 1):8.3f}%")
    print(f"{'D_M(z=' + str(z) + ')':22s} {ours_D:12.3f} {camb_D:12.3f} {100 * (ours_D / camb_D - 1):8.3f}%")
rd = r_drag(obc, ob)
print(f"{'r_d (DESI fit vs CAMB)':22s} {rd:12.3f} {der['rdrag']:12.3f} {100 * (rd / der['rdrag'] - 1):8.3f}%")
zst = z_star(ob, obc)
print(f"{'z* (calibrated vs CAMB)':22s} {zst:12.2f} {der['zstar']:12.2f} {100 * (zst / der['zstar'] - 1):8.3f}%")
R, lA, _ = c.cmb()
om = (obc + OMEGA_NU) / h ** 2
DMs = res.comoving_radial_distance(der["zstar"])
R_c = np.sqrt(om) * 100 * h * DMs / 299792.458
lA_c = np.pi * DMs / der["rstar"]
print(f"{'R (shift parameter)':22s} {R:12.5f} {R_c:12.5f} {100 * (R / R_c - 1):8.3f}%")
print(f"{'l_A (acoustic scale)':22s} {lA:12.4f} {lA_c:12.4f} {100 * (lA / lA_c - 1):8.3f}%")

# growth shape: f*sigma8(z) / f*sigma8(0), GR, scale-independent part
fs8_c = np.array(res.get_fsigma8())[::-1]  # CAMB returns in order of decreasing z as requested -> reverse to increasing z
m = DragModel(p, 0.0)
fs8_o = np.array([m.fs8(z, 1.0) for z in zs])
print("\nGR growth shape f*sigma8(z)/f*sigma8(0):")
for z, a_, b_ in zip(zs, fs8_o / fs8_o[0], fs8_c / fs8_c[0]):
    print(f"  z={z:3.1f}: ours {a_:.4f}  CAMB {b_:.4f}  diff {100 * (a_ / b_ - 1):+.3f}%")

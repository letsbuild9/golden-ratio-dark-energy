# Pre-registered predictions: golden-ratio phase-space hypothesis (version 2)

**Author:** Siddesh I. Borse
**Locked on:** {LOCK_DATE}
**Integrity:** the SHA-256 hash of this file is recorded in `PREDICTIONS.sha256`. Any later edit changes the hash. The public Zenodo record fixes the date.
**History:** this is the first public version. An earlier local draft, never published, is superseded by this file.

## Data status at the time of locking

| Dataset | Status |
|---|---|
| DESI DR2 BAO (arXiv:2503.14738; observations 14 May 2021 – 9 April 2024) | Public. **Used to choose these values, so it is not a test of them.** |
| DESI analyses including observations **after 9 April 2024** | Not public; no release date announced. DESI's public spectroscopic releases so far are the EDR and DR1 (March 2025). **First DESI test of P1–P4.** |
| DES-Dovekie supernova recalibration (arXiv:2511.07517) | Public; used in the analysis |
| Euclid DR1 | Data release from November 2026; cosmology results not expected before mid-2027 (secondary sources). Decisive growth measurements need the full survey. **Test of P5.** |
| CMB early-dark-energy limits | 95% upper limits from 0.065 to 0.12 depending on the data combination (Qu et al. 2024; Poulin et al. 2025; Reeves et al. 2025). **P6 is under pressure.** |

## The rule

- **The constant:** x = 1/φ = (√5 − 1)/2 ≈ 0.618034, the positive root of x² + x − 1 = 0 (the fixed point of x = 1/(1 + x)).
- **The ladder:** each dark-sector parameter that differs from ΛCDM is set by a ladder term **O_d = α_d x^d**, where α_d = 1 for linear or volumetric modes (d = 1, 3) and α_d = ½ for transverse or plasma modes (d = 2, 4). For 1 + w₀, wₐ and f_EDE (all zero in ΛCDM) the ladder term is the departure itself; for growth it is the value of γ, not its shift from the general-relativity value γ_GR ≈ 0.55.
- **The assignment:**
  - late-time dark energy: 1 + w₀ = O₂;
  - growth: γ = O₁;
  - early dark energy: f_EDE = O₄;
  - wₐ = −x. Equivalently, the crossing of w = −1 lies at 1 − a = x/2. This is an equivalent postulate, not a derivation.

## Predictions

| ID | Quantity | Exact | Value |
|---|---|---|---|
| **P1** | Dark-energy equation of state today, w₀ (CPL) | −φ/2 | **−0.80902** |
| **P2** | Dark-energy evolution, wₐ (CPL, a local parametrization) | −1/φ (equivalently wₐ = 1 + 2w₀, or crossing at 1 − a = x/2) | **−0.61803** |
| **P3** | Redshift where w = −1 (dark-energy density maximum) | 1/√5 = 1/(2φ − 1) | **0.44721** |
| **P4** | w at fixed redshifts | w(0.3) = −(7 + 19√5)/52; w(0.5) = −(1 + 5√5)/12 | **−0.9516; −1.0150** |
| **P5** | Growth index γ (constant for 0 ≤ z ≤ 1) | 1/φ | **0.618** |
| **P6** | Peak early-dark-energy fraction (axion-like EDE, n = 3, z_c free) | x⁴/2 = (5 − 3φ)/2 | **0.07295** |
| **P7** | Clustering amplitude S8 (derived from P1, P2, P5; CMB-anchored: lower σ₈, higher Ω_m) | — | **≈ 0.81** (0.812 in the best-fit background; approximate) |

**Explicitly not predicted:** H₀.

## Falsification criteria (fixed now, before the data)

- **P1–P4 are falsified** if, in the first DESI analysis including observations after 9 April 2024, either of the following holds:
  - the point (−0.809, −0.618) lies outside the 3σ region (Δχ² > 11.83) of the w₀–wₐ constraint, from BAO + CMB alone or from BAO + CMB + supernovae for at least two of the principal SN compilations;
  - w(0.3) or w(0.5) differs from the prediction by more than 3σ in at least two of these combinations.
- **P1–P4 are supported** if the point stays within 2σ in all combinations while ΛCDM is excluded at more than 3σ.
- **P5 is falsified** by a measurement with σ(γ) ≤ 0.03 lying more than 3σ from 0.618.
- **P6 is falsified** if CMB data exclude f_EDE ≥ 0.073 at 99.7% for the n = 3 axion-like model with z_c free.

**Overlap caveat.** Every future DESI analysis also contains the DR2 data used to choose P1–P4, so part of any future agreement is built in. Fully independent tests come from other probes: Euclid growth (P5), new supernova samples, and CMB experiments (P6).

## Status at lock time (honest caveats)

1. **P1–P4 were chosen after DESI DR2.** Current agreement is therefore consistency, not confirmation:
   - official DESI DR2 chains (full CMB): 0.5σ with Pantheon+, 1.0σ with Union3, 0.5σ with DESY5, and 1.45σ with no supernovae;
   - compressed-CMB fits: 0.3–1.0σ from the best free fit, including DES-Dovekie.
2. **The look-elsewhere test:** of 403 pairs from the same algebraic family, (−0.809, −0.618) is the only one within 1σ for all three supernova samples.
3. **Robustness:**
   - without CMB data and with r_d free (independent of early-universe physics), the point is 0.1–1.1σ from the best fit;
   - including the predicted f_EDE = 0.073 (indicative, compressed CMB) it is 0.2–1.3σ;
   - against DESI's model-independent binned w(z), the predicted curve deviates by 1.0–1.75σ, comparable to DESI's own best-fit CPL models (0.3–1.5σ).
4. **Fitting x freely gives x = 0.55–0.69 depending on the data.** 1/φ is the only one of nine simple constants inside every 68% range, but neighboring values (0.6, π/5) fit almost as well. The data cannot single out φ.
5. **Growth:** current CMB-anchored fσ8 data (SDSS and DESI DR1) place γ = 1/φ 1.1–1.3σ from the best fit, a mild, non-significant preference for general relativity.
6. **Postulates:** the α_d weights and the assignment of d to each observable are postulates, and wₐ = −1/φ is not derived.

## Reproducibility

- **Analysis code:**
  - `borse_fit.py`
  - `borse_growth.py`
  - `borse_model.py`
  - `borse_mu.py`
  - `borse_desi_chains.py`
  - `borse_v2_tests.py`
  - `borse_checks.py`
  - `borse_camb_check.py`
  - `borse_calibrate_zstar.py`
  - `borse_v2b_tests.py`
  - `borse_binned_calibration.py`
  - `borse_energy_exchange.py`
- **Verification:**
  - `borse_symbolic.py` (exact symbolic checks);
  - `verify_paper.py` (checks every number in the paper against the outputs);
  - `VERIFICATION.md` (the source of every external value).

# Verification ledger (version 2)

**Paper:** *Golden-Ratio Non-Resonance in Phase Space as an Organizing Hypothesis for Dark-Sector Anomalies* (S. I. Borse)
**Prediction record:** `PREDICTIONS.md`, version 2, locked 2026-10-19; SHA-256 in `PREDICTIONS.sha256`
**Public record:** Zenodo, DOI [10.5281/zenodo.23481070](https://doi.org/10.5281/zenodo.23481070), published 2026-10-19 18:05 UTC; file checksums on Zenodo match the local files.

| Class | How it was checked | Result |
|---|---|---|
| A. Mathematics | Exact symbolic proofs with SymPy (`borse_symbolic.py` → `rerun/symbolic_proofs.txt`) | **31 / 31 pass** |
| B. Numbers we computed | Clean re-run of every script from raw data; `verify_paper.py` formats each output as printed and finds it in `paper.tex` (`rerun/verify_paper_output.txt`) | **122 / 122 pass**; pipeline deterministic |
| C. Our code against professional software | CAMB 2.0.5 cross-check (`borse_camb_check.py`) | Agreement 0.000–0.03% (details below) |
| D. Numbers quoted from other papers | Found verbatim in each original PDF, or in the official data products | All found (table below) |
| E. References | arXiv API records (titles, authors, journal refs) and Crossref DOIs | All verified, except 4 classic pre-DOI items (listed) |
| F. Text | Machine spell check; cross-document consistency (paper ↔ PREDICTIONS ↔ Zenodo description) | No misspellings; all prediction values identical |

## How to reproduce

1. `.venv/bin/python borse_calibrate_zstar.py` (needs `camb`; writes `data/zstar_calibration.json`)
2. `borse_checks.py`
3. `borse_fit.py --grid`
4. `borse_growth.py` and `borse_growth.py --desi`
5. `borse_model.py`
6. `borse_mu.py`
7. `borse_desi_chains.py`
8. `borse_v2_tests.py`
9. `borse_v2b_tests.py`
10. `borse_binned_calibration.py`
11. `borse_energy_exchange.py`
12. `borse_camb_check.py`
13. `borse_symbolic.py`
14. `verify_paper.py`

Outputs go to `rerun/`. Data sources are listed in each script header.

---

## A. Mathematics (exact)

- **Golden ratio:**
  - x = 1/φ;
  - (1 − x)/x = x;
  - x² + x − 1 = 0, with exactly one positive root;
  - x = 1/(1 + x);
  - 1 + 1/φ = φ.
- **Lemma 2:**
  - x² = 1 − x;
  - x³ = 2x − 1 = √5 − 2;
  - x⁴ = 2 − 3x;
  - 1 − x²/2 = φ/2;
  - x⁴/2 = (5 − 3φ)/2;
  - x/(2 − x) = 1/√5 = 1/(2φ − 1).
- **Dark energy:**
  - w₀ = −φ/2 and wₐ = −1/φ;
  - wₐ = 1 + 2w₀ (and, on the ladder curve, wₐ − (1 + 2w₀) = −(x² + x − 1));
  - 1 + w₀ + wₐ = (1 − 3x)/2, and the ρ_DE exponents;
  - (1 + w₀)/(−wₐ) = x/2;
  - w = −1 exactly at z = 1/√5;
  - the quintom-B conditions (w₀ > −1, w₀ + wₐ < −1);
  - the fixed-z values w(0.3) = −(7 + 19√5)/52 and w(0.5) = −(1 + 5√5)/12;
  - the non-crossing bound;
  - the d = 3 branch.
- **Remark:** 6/11 + x⁴/2 − 1/φ = (5/2)(89/55 − φ).
- **Proposition 3:** Fibonacci ratios and the contraction bound, checked exactly for n ≤ 25.
- **Contraction (Sec. 2.2):** f(f(u)) = (1 + u)/(2 + u), with derivative 1/(2 + u)² ≤ 1/4 on u ≥ 0, so f∘f is a contraction with constant 1/4; f′(0) = −1, so f itself is not a contraction on [0, ∞).

## B. Computed numbers

- **Every number in the paper's tables and text** that comes from our analysis is checked by `verify_paper.py` (111 items):
  - Tables 1–5;
  - every quoted range (σ values, Δχ², H₀, Ω_m, x ranges, Bayes factors, pivot values);
  - the μ mapping;
  - the growth fits and S8;
  - the drag-model and energy-exchange numbers;
  - the early-dark-energy numbers;
  - the CAMB and Planck validation;
  - the z* calibration residual;
  - the lock hash.
- **Determinism:** two identical runs give byte-identical output.

## C. Code validation

- **Against CAMB 2.0.5** (φ-ladder best fit): H(z) and D_M(z) to ≤ 0.001%, R 0.001%, ℓ_A 0.002%, z* 0.000%, r_d 0.03% (the accuracy of the DESI fitting formula), f σ₈ shape ≤ 0.23% to z = 2.
- **Against Planck 2018** (TT,TE,EE + lowE + lensing parameters): r_s(z*) = 144.43, r_d = 147.09 and 100θ* = 1.04110, identical to the published values.

## D. External values (verified in the original sources)

| Value | Source | Location |
|---|---|---|
| ΛCDM rejected at 2.8 / 3.8 / 4.2σ (with SNe); 3.1σ for DESI + CMB; w₀–wₐ posteriors; PPF used when w crosses −1; priors w₀ ∈ U[−3, 1], wₐ ∈ U[−3, 2], w₀ + wₐ < 0; DR2 observations 14 May 2021 – 9 Apr 2024 | DESI DR2, arXiv:2503.14738 | Main text and parameter table |
| r_d fit 147.05 (ω_m/0.1432)^−0.23 (N_eff/3.04)^−0.1 (ω_b/0.02236)^−0.13 | DESI 2024 VI, arXiv:2404.03002 | Eq. 2.5 |
| Binned w(z) model: edges 0.7, 1.4, 2.1, 10, 15; w_bin_4 = −1; fixed smoothing parameter `smoothing_fac` = 0.02 (functional form not documented in the release, so the paper states only the parameter); priors U[−3, 1] | DESI chain release | `chain.updated.yaml` |
| R, l_A, ω_b and correlation matrix; **z* taken from the Planck chains** | Chen, Huang & Wang, arXiv:1808.05724 | Table I and Sec. II |
| r*, r_drag, 100θ*, σ₈, H₀, ω_b, ω_c | Planck 2018 VI, arXiv:1807.06209 | Tables 1–2 |
| z* fit coefficients | Hu & Sugiyama 1996 | Eq. E-1 |
| Massive-neutrino density f(y) ≈ [1 + (Ay)^p]^{1/p}, A = 0.3173, p = 1.83 | Komatsu et al. 2011, ApJS 192, 18 | Eq. 26 |
| Drag ∝ (1 + w)ρ_DE σ_D | Simpson, arXiv:1007.1034 | Eq. 2 |
| α ≈ 6/11 for ΛCDM; 3/(5 − w/(1 − w)) leading term | Wang & Steinhardt 1998 | Eqs. 13 and B11 |
| V ∝ (1 − cos[ϕ/f])ⁿ; w_n = (n − 1)/(n + 1) | Poulin et al., arXiv:1811.04083 | Main text |
| f_EDE < 0.070 (Planck + lensing), < 0.091 (adding DESI) | Qu et al., arXiv:2404.16805 | Abstract |
| f_EDE < 0.065 (NPIPE-LS), < 0.092 (P-ACT-LS), < 0.109 (P-ACT-LBS), 95%, no SH0ES prior. Labels: L = CMB lensing, S = Pantheon+, B = DESI DR2 BAO, so the first two do **not** include DESI | Poulin et al., arXiv:2505.08051 | Results text and Appendix Table 1 |
| f_EDE < 0.12 (ACT DR6 combination) | Reeves et al., arXiv:2510.06114 | Results table |
| μ₀ = 0.11 (+0.44/−0.54); 0.05 ± 0.22; Σ₀ = 0.008 ± 0.045 (ΛCDM); μ₀ = −0.24 (+0.32/−0.28), Σ₀ = 0.006 ± 0.043 (w₀wₐ) | Ishak et al., arXiv:2411.12026 | Abstract and Table 2 |
| Quintom-B: phantom → quintessence; region w₀ > −1, w₀ + wₐ < −1 | Qiu, Cai et al., arXiv:2511.19994 | Fig. 1 caption and text |
| Interaction ⇒ effective w < −1 | Das, Corasaniti & Khoury, arXiv:astro-ph/0510628 | Abstract |
| Chaos: closed FRW + scalar field; two-field cosmology | Page 1984; Calzetta & El Hasi 1993; Cornish & Shellard 1998; Easther & Maeda 1999 | Crossref and arXiv records, abstracts |
| DESI DR2 Lyα full-shape: fσ8 excluded (10% bias, ≫5σ) | Herbold et al., arXiv:2607.27411 | Abstract and text |
| Unified-tracer DR2 reanalysis: w₀wₐ "extremely consistent" with the baseline | Sanders et al., arXiv:2608.27830 | Text |
| DES-Dovekie data (1820 SNe, packed inverse covariance) | Popovic et al., arXiv:2511.07517; Cobaya `sn_data` | Loader matches `cobaya/likelihoods/base_classes/sn.py` |
| DESI DR2 + CMB (Planck, ACT, SPT) + DES-Dovekie: w₀ = −0.803 ± 0.054, wₐ = −0.72 ± 0.21; 3.2σ (reduced from 4.2σ for DES-SN5YR) | Popovic et al., arXiv:2511.07517 | Abstract |
| DESI + CMB + Pantheon+: w₀ = −0.838 ± 0.055, wₐ = −0.62 (+0.22/−0.19) | DESI DR2, arXiv:2503.14738 | Eq. 26 and parameter table |
| H₀ = 73.04 ± 1.04 | Riess et al. 2022 | Abstract |
| S8 = 0.776 ± 0.017 | DES Y3, arXiv:2105.13549 | Abstract |
| S8 = 0.815 (+0.016/−0.021) | KiDS-Legacy, A&A 703, A158 | Abstract |
| fσ8(0.067) = 0.423 ± 0.055 | Beutler et al. 2012 | Abstract |
| MGS 0.53 ± 0.16; BOSS/eBOSS values | eBOSS DR16, arXiv:2007.08991 | Table 3 |
| fσ8(0.07) = 0.450 ± 0.055; γ_L = 0.580 ± 0.110 | Turner et al., arXiv:2512.03230 | Table 3 and abstract |
| DESI DR1 ShapeFit fσ_s8 (6 bins) | Adame et al., arXiv:2411.12021 | Appendix A (parsed programmatically) |
| γ = 0.633 (+0.025/−0.024) | Nguyen, Huterer & Wen, PRL 131, 111001 | Abstract |
| K_c ≈ 0.971635 | Greene 1979; MacKay 1983 | Secondary literature |
| DESI public spectroscopic releases: EDR (2023) and DR1 (March 2025); no announced date for post-DR2 cosmology | data.desi.lbl.gov release page | Official page |
| DESI data-use acknowledgment text | data.desi.lbl.gov acknowledgments page | Reproduced verbatim |
| SDSS-IV acknowledgment text | sdss4.org/collaboration/citing-sdss | Reproduced verbatim, with the full member list |

## E. References

- **Verified via arXiv records and/or Crossref DOIs:** all cited papers with an arXiv version or a DOI.
- **Classic pre-DOI items, checked only against standard bibliographic knowledge:**
  - Kolmogorov 1954;
  - Moser 1962;
  - Arnold 1964 (Sov. Math. Dokl.);
  - the Arnold and Sakurai textbooks.

---

## Corrections log

### First pass
1. **Thawing coefficient.** wₐ ≈ −1.58(1 + w₀) was unverifiable; it was removed and replaced by the provable non-crossing bound.
2. **MGS error bar:** ±0.175 → ±0.16.
3. **Planck reference parameters** made consistent.
4. **Overclaim removed:** "golden tori proven most robust" relabeled as evidence.
5. **Speculation removed:** an unverified explanation of the Union3 shift was deleted.
6. **μ₀ error:** +0.45 → +0.44.
7. **Pantheon+ count:** "SNe" → "light curves".
8. **References:** several corrected.
9. **Notation:** φ clash with the EDE field resolved.

### Second pass
1. **r_d formula:** uses ω_m (all matter), as in DESI's source.
2. **N_eff:** made consistent (3.046).
3. **Rounding:** 0.630 → 0.629.
4. **Spelling:** standardized to American.
5. **Figure labels:** made to match the paper's terms.
6. **Ω_m rounding** corrected.
7. **Layout:** LaTeX overflows fixed.

### Third pass (version 2)
1. **CMB decoupling redshift.** The Chen et al. priors use z* from the Planck chains, but our code used the Hu–Sugiyama fit, which is 0.18% high. This biased ℓ_A by 0.2% (about 6σ of the prior). Fixed with a CAMB-calibrated correction (residual < 4 × 10⁻⁷).
2. **Massive neutrino.** It was treated as pure matter and also as radiation, biasing r_s by 0.11%. Now uses the Komatsu et al. fitting formula at all redshifts.
3. **Result of fixes 1–2:** Planck values are reproduced exactly, and CAMB agreement is 0.001–0.03%. All fits were re-run; the conclusions are unchanged.
4. **New data:** DES-Dovekie supernovae added.
5. **New tests:**
   - free-x fit;
   - control constants;
   - Bayes factors (with their peak-normalized context);
   - pivot and fixed-z tests;
   - robustness across alternative DESI chains;
   - goodness of fit;
   - the joint early-dark-energy check;
   - DESI binned w(z) (calibrated against DESI's own CPL models);
   - the energy-exchange test.
6. **Sign error corrected (Open problem 4).** An apparent phantom past requires energy flow dark energy → dark matter, not the reverse. Confirmed analytically and by Das, Corasaniti & Khoury.
7. **Growth rule changed to γ = 1/φ**, following the author's α_d rule (Brief 2). The earlier γ_GR + F(4) form is kept as a remark, and the P6 consistency relation is dropped.
8. **Test definition corrected.** The first DESI test is now "data after 9 April 2024": DESI DR3 may contain no new observations, and no date is announced.
9. **Bridge rewritten.** The flat single-field universe is non-chaotic (a 2D autonomous system), so the bridge now invokes the closed single-field and two-field cases, with verified references.
10. **Wording:** wₐ restated as an equivalent postulate (crossing at 1 − a = x/2); the quantum remark shortened; CPL described as a local parametrization.
11. **Disclosures and acknowledgments:** DESI, SDSS and Planck acknowledgments added; AI disclosure now names both Gemini and Claude.
12. **Checker bugs:** a placeholder symbolic check that always passed was found and removed.

### Fourth pass (referee-style review of version 2)
1. **Missing comparison added (M1).** The published DES-Dovekie constraint (w₀ = −0.803 ± 0.054, wₐ = −0.72 ± 0.21) is now quoted in the introduction and Sec. 4.2. The prediction lies 0.1σ and 0.5σ from it, one parameter at a time; no joint level is given because no correlation is published.
2. **Hypothesis wording (M2).** "Each departure from ΛCDM is a ladder term" was inconsistent with γ = 1/φ, which is an absolute value of the growth index, not a shift. The paper, PREDICTIONS.md and the Zenodo description now say that dark-sector parameters are *set by* ladder values, and state explicitly that the ladder value is the departure itself for 1 + w₀, wₐ and f_EDE but the value of γ for growth. The prediction record was re-locked the same day (2026-10-19); new SHA-256 in `PREDICTIONS.sha256`. No predicted value changed.
3. **Drag test retargeted (M3).** `borse_model.py` used the old γ_GR + F(4) target. It now uses γ = 1/φ. Target suppression 2.64% (was 2.63%); non-crossing coupling ξ ≈ 2.3, σ_D/m_c ≈ 36 b/GeV (was 2.4 and 38); Δγ = 0.066–0.072 against a required 0.059–0.069. Conclusion unchanged.
4. **Mathematical statement (M4).** "Contraction-mapping (Banach) fixed point" was imprecise: f has Lipschitz constant 1 on [0, ∞). The paper now says the bound contracts distances to the fixed point and that f∘f is a contraction with constant 1/4 (3 new exact checks).
5. **Causal overreach removed (M5).** The pivot sentence claiming to explain why single datasets look like ΛCDM was replaced by the verifiable statement that the predicted crossing lies inside the pivot range 0.31–0.53.
6. **Bayes-factor wording (M6).** "Essentially zero for DESY5" → the ΛCDM point lies outside the sampled region. "For most samples, over ΛCDM" → over ΛCDM in every combination, with conservative lower bounds B(φ:ΛCDM) ≳ 30, 110 and 960 (min B_φ / max B_ΛCDM).
7. **Abstract shortened (M7)** from 2,510 to 1,733 characters (arXiv limit 1,920).
8. **Factual error fixed (minor 2).** The Poulin et al. 2025 limits 0.065 and 0.092 were described as "with DESI DR2 and Pantheon+"; the source's NPIPE-LS and P-ACT-LS combinations do not include DESI BAO. Corrected, and the P-ACT-LBS limit (0.109, with DESI) added.
9. **Binned w(z) sentence (minor 1).** "Both CPL curves are more negative below z = 0.7" is false for DESI's own DESI+CMB fit (−0.842 vs −0.840); now qualified, and "is a limitation" softened to "appears to be".
10. **Precision of statements (minors 3–7).** Noble numbers defined; KAM non-degeneracy condition stated; "requires at least two fields" → a single canonical field cannot cross; "is again chaotic" → "can be chaotic"; "reproduce exactly" → "to the quoted precision", with the published values; S8 "± 0.01" (no stated origin) removed.
11. **Style (minors 8–14).** H₀ units added; "the author's stated rule" removed; four missing scripts added to Data availability; Table 3 sign convention stated; DESI binned smoothing parameter mentioned; full SDSS-IV acknowledgment; the binned, joint early-dark-energy and drag tests moved to Appendices A–C, with a summary subsection (Sec. 4.5) in the main text.

## Items that cannot be fully verified, and how the paper handles them

- **DESI fσ_s8 values** are used as fσ8. The paper states this requires the fiducial sound horizon to be close to the true one.
- **Compressed CMB priors** assume standard early physics. The joint early-dark-energy test with the CMB is therefore labeled indicative; the r_d-free, no-CMB test is the robust one.
- **Binned w(z) comparison:** uniform bin-averaging is approximate, so it is calibrated against DESI's own CPL models.
- **Energy exchange:** tested at background level only; growth (μ₀) depends on the momentum-transfer prescription and is not computed.
- **Euclid timeline:** partly from secondary sources.
- **The bridge (Open problem 1)** is explicitly a conjecture.

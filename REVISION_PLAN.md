# Master revision plan (version 2): do everything once, in this order

Status (19 October 2026): **Phases 0–5 and K1–K3 done.**
- D1 = B, D2 = re-lock, D3 = no, D4 = T9 dropped (superseded by T1/T2), T10 skipped (no verified forecasts), D5 = included.
- Review additions were included: the joint early-dark-energy check, binned w(z), energy-exchange sign fix and test, the corrected bridge, and the S8 note.
- Referee-style review of v2 applied (7 major, 14 minor items; see VERIFICATION.md, fourth pass). Prediction record re-locked the same day with corrected hypothesis wording; no predicted value changed.
- K4 done: Zenodo DOI 10.5281/zenodo.23481070, published 2026-10-19 18:05 UTC. K5 done: DOI in the paper (2 places), PDF and tarball rebuilt, 122/122 checks.
- D6 done: GitHub https://github.com/letsbuild9/golden-ratio-dark-energy (initial commit 22bd36e).
- K6 in progress: arXiv account created; astro-ph.CO submission awaits endorsement (code JEEE4W). Remaining: revoke the Zenodo token, push later commits, arXiv upload after endorsement, D7.
Rule: no paper edits until Phases 0–2 are complete. Every new number must be produced by a script, so it can be re-verified automatically.

---

## Phase 0: Decisions (needed from the author before anything starts)

| ID | Decision | Options | Recommendation |
|---|---|---|---|
| D1 | Growth rule | **A:** Δγ = F(4) (keeps the P6 cross-check) · **B:** γ = x with the α_d rule (author's own scheme) | B, slightly |
| D2 | Re-lock the predictions as version 2 before the Zenodo upload | Yes / no | Yes |
| D3 | Attach the unpublished version-1 draft to Zenodo as "superseded"? | Yes / no | No (keep one clean version) |
| D4 | Optional tests T9 (mock data) and T10 (forecasts) | Include / skip | Skip T10 unless verified forecast numbers are found |
| D5 | Optional paper item: strengthen "Beyond cosmology" with accepted explanations (folding funnel; Schurger et al. 2012) | Include / skip | Include (short) |
| D6 | Public GitHub repository with tests, linked to Zenodo | Yes / no / later | Later (after arXiv) |
| D7 | Independent human review before arXiv (a cosmologist) | Arrange / skip | Arrange (also helps with endorsement) |

---

## Phase 1: Data and facts to gather (verify each fact in its original source)

| ID | Item | Source | Verify |
|---|---|---|---|
| G1 | DES-Dovekie supernova data and covariance | DES-SN5YR GitHub (`4_DISTANCES_COVMAT/DES-Dovekie_*`, `STAT+SYS.npz`) | File format and redshift column, as for DESY5 |
| G2 | DESI chains for alternative CMB likelihoods (plik; HiLLiPoP/LoLLiPoP), DESY5 with z > 0.1, and DESI-only / BAO-only (no CMB) | data.desi.lbl.gov `bao-cosmo-params/cobaya/base_w_wa/` | Folder names and dataset definitions in the README |
| G3 | DESI's w₀–wₐ priors (for the Bayes factor) | `chain.updated.yaml` in the chains | Exact prior ranges and the w₀ + wₐ < 0 cut |
| G4 | DESI data timeline: which observations DR3 contains, and when results using data after 9 April 2024 are expected | DESI official site / papers | Confirmed in a primary source, or stated as unknown |
| G5 | Unified-tracer reanalysis of DR2 BAO (arXiv:2608.27830): does it shift w₀–wₐ? | arXiv | Values from the paper text |
| G6 | DESI DR2 Lyα full-shape (arXiv:2607.27411): fσ8 excluded due to a 10% bias | arXiv | Quote from the paper |
| G7 | "Quintom-B" terminology (phantom → quintessence crossing) | A review, e.g. Cai, Saridakis et al. (2010) | Definition in the original |
| G8 | Prior art: golden-ratio cosmology papers and "parameter-free" dark-energy proposals | arXiv / ADS search | Read the abstracts and cite the relevant ones |
| G9 | Data-use acknowledgment texts (DESI, SDSS-IV, Planck/ESA, DES, Pantheon+, Union3) | Each collaboration's policy page | Exact wording |
| G10 | Classic reference page numbers (Hurwitz 1891, Kolmogorov 1954, Arnold 1963/64, Moser 1962, Greene 1979, MacKay 1983, Chevallier–Polarski 2001, Wang–Steinhardt 1998, Hu–Sugiyama 1996) | Publisher / ADS | Volume, page and year |
| G11 | (If D5) Schurger et al. 2012 (PNAS) and folding-funnel references | PubMed / ADS | Correct citation |
| G12 | (If T10) Euclid γ forecast and DESI final-survey forecast | Euclid Collaboration forecast papers | Numbers from the original, otherwise drop T10 |
| G13 | Newest early-dark-energy limits (ACT DR6, SPT-3G) | arXiv | Exact limits and datasets |

---

## Phase 2: Tests and analyses (each as a script with saved output)

| ID | Test | Method | Output for the paper | Report even if unfavorable? |
|---|---|---|---|---|
| T1 | **Measure x from data** | Fit the one-parameter ladder (w₀ = −1 + α₂x², wₐ = −x) with the compressed pipeline (all 4 combinations), and profile along the same curve in the official DESI chains | x ± σ per dataset; distance from 0.618 | Yes |
| T2 | **Control constants** | Same fit with x fixed to 1/√2, 2/3, 1/e, 0.6, 0.65, …, plus a continuous χ²(x) curve | Plot and table: is 0.618 special? | Yes |
| T3 | **DES-Dovekie** | Repeat Table 1, the look-elsewhere test and T1 with Dovekie in place of DESY5 | Extra rows and columns | Yes |
| T4 | **Bayes factor** | Savage–Dickey density ratio at (−φ/2, −1/φ) from the official chains, with DESI's priors (G3) | Bayes factor and its prior dependence | Yes |
| T5 | **Pivot / fixed-z test** | DESI pivot table (already computed), plus theory w(z) at fixed z = 0.3 and 0.5 versus the chains | Table; fixed-z values become locked predictions | Yes |
| T6 | **Absolute goodness of fit** | χ²/degrees of freedom and p-value for every Table 1 model | Extra column or sentence | Yes |
| T7 | **CAMB cross-validation** | Install `camb`; compare H(z), D_M, D_H, r_d and fσ8 (PPF dark energy) for the φ-ladder best fit against our code | Agreement level (target < 0.1% for distances) | Yes; fix the code if it disagrees |
| T8 | **Robustness** | Official-chain position of the point for alternative CMB likelihoods, DESY5 with z > 0.1, and DESI without CMB | Table or sentence | Yes |
| T9 | (D4) Mock-data p-value for the look-elsewhere effect | Simulate data from ΛCDM and w₀wₐ; count chance matches within the family | p-value | Yes |
| T10 | (D4) Forecasts | Only with verified forecast errors (G12) | Expected σ separation | Skip if unverified |
| T11 | **Growth rule per D1** | If B: recompute μ₀ for constant γ = x, and update the growth tables and the S8 statement. If A: no new numbers | Updated μ₀ and growth text | Yes |
| T12 | **Background history as a script** | Dark-energy density peak (= z_cross), dark-energy/matter equality, acceleration onset, acceleration end (a = 1.67), CPL breakdown (w > 1 at a = 3.93), lookback time to z_cross | Numbers for the "meaning of z_cross" and "transient acceleration" remarks | Yes |
| T13 | **New symbolic proofs** | wₐ = 1 + 2w₀ ⟺ x² + x − 1 = 0; z_cross = 1/(2φ − 1); α_d ladder values (x, x²/2, x³, x⁴/2) | Added to `rerun/symbolic_proofs.txt` | — |
| T14 | **Recheck after G5** | If the unified-tracer reanalysis changes w₀–wₐ, note it; re-run the tests only if DESI publishes new chains | Sentence | Yes |

---

## Phase 3: Paper changes (only after Phase 2)

**Front matter**
- **P1 Title and abstract:** update all numbers. Add the free-x result, the pivot result and the Bayes factor. Adjust the growth sentence per D1.
- **P2 Introduction:** add a "Related work / prior art" paragraph (G8). Correct the data-timeline statement (G4).

**Section 2: mathematics**
- **P3** Add wₐ = 1 + 2w₀ as a lemma, and the form z_cross = 1/(2φ − 1).

**Section 3: hypothesis**
- **P4** Replace "½ not derived" with the author's α_d rule (Brief 2), labeled as a postulate. State the d = 3 "volumetric" value x³.
- **P5** Growth rule per D1, with the alternative form mentioned as a numerically near-equivalent variant.
- **P6** Remarks with T12 numbers:
  - z_cross is when the dark-energy density peaked;
  - dark energy is diluting now, so there is no Big Rip (quintom-B, per G7);
  - acceleration is transient, but only qualitatively, with the CPL-breakdown caveat;
  - the correct timeline: dark energy overtook matter at z ≈ 0.33 and acceleration began at z ≈ 0.72.

**Section 4: tests**
- **P7** Add sections or tables for T1, T2, T3, T4, T5, T6, T7 and T8 (and T9 if done).
- **P8** Mention the DESI DR2 Lyα full-shape (G6) and the unified-tracer reanalysis (G5).
- **P9** Update the early-dark-energy limits (G13).

**Section 5: predictions**
- **P10** Rewrite per version 2:
  - fixed-z w predictions;
  - falsification defined on DESI data after 9 April 2024 plus independent probes;
  - an explicit note that future DESI results overlap with DR2.
- **P11** Remove most of "Notes added after the lock" (the corrections are now made before the lock). Keep a short history note only if D3 = yes.

**Sections 6–7 and back matter**
- **P12** Open problems: update (α_d still a postulate; mechanism; bridge).
- **P13** Optional (D5): beyond cosmology, with the accepted explanations and what a testable prediction would need.
- **P14** Data acknowledgments (G9).
- **P15** AI disclosure: name both Gemini (origin of the framework) and Claude (analysis and drafting).
- **P16** References: add all new ones; fix the classic page numbers (G10).
- **P17** Data availability: new scripts and the verification ledger.

---

## Phase 4: Prediction record, version 2

- **L1 New PREDICTIONS.md:**
  - version 2, dated with the **actual upload day**;
  - the growth rule per D1;
  - fixed-z w(z = 0.3) and w(z = 0.5);
  - the α_d rule;
  - wₐ = 1 + 2w₀ as an equivalent statement;
  - corrected rounding (0.629);
  - the second-hand 0.061 limit removed;
  - DR2 agreement stated from the official chains (Union3 at 70%);
  - the test definition from P10;
  - the overlap caveat.
- **L2** A new SHA-256 hash in `PREDICTIONS.sha256`.
- **L3** The paper's lock section updated (date and hash).

---

## Phase 5: Verification (all must pass before upload)

| ID | Check | Pass criterion |
|---|---|---|
| V1 | Clean re-run of every script from raw data (caches deleted) | Completes; outputs saved in `rerun/` |
| V2 | Determinism | A repeated run gives identical output |
| V3 | `verify_paper.py` extended to every new number (T1–T9, T11–T12) | 100% pass |
| V4 | Symbolic proofs (T13 added) | 100% pass |
| V5 | Source checks: every new external number found in its original PDF | All found |
| V6 | arXiv-record check of every new reference | Titles, authors and journals match |
| V7 | Spell check and consistent spelling (American) | No real misspellings |
| V8 | Cross-document consistency: paper ↔ PREDICTIONS v2 ↔ Zenodo form ↔ upload-script description ↔ figure labels | No conflicts |
| V9 | LaTeX build | No errors, no overfull lines, no "??" |
| V10 | Visual check of every changed page and figure | OK |
| V11 | VERIFICATION.md updated (new tests, new sources, corrections log) | Complete |
| V12 | PREDICTIONS hash check (`shasum -c`) and an identical Zenodo copy | OK |

---

## Phase 6: Packaging and publishing

- **K1** Zenodo form and upload script: title (version 2, real date), publication date = upload day, description, expected hash.
- **K2** Rebuild `borse_analysis_v1.zip` (consider renaming it v2) with all scripts, `rerun/` outputs, the ledger and `verify_paper.py`.
- **K3** Rebuild `paper/arxiv_submission.tar.gz`.
- **K4** Author: Zenodo login → token → run the script → review the draft → **Publish**.
- **K5** Insert the DOI in the paper (2 places) and the ledger → rebuild → re-run V3 and V9.
- **K6** arXiv: account, endorsement (astro-ph.CO), upload the tarball.
- **K7** (D6) GitHub repository linked to Zenodo, if chosen.

---

## Definition of done

1. All Phase 2 tests run, with outputs saved and their results reported whether favorable or not.
2. V1–V12 all pass.
3. Predictions version 2 is public on Zenodo with the real date, and the DOI is in the paper.
4. No number in the paper is typed by hand without a script output behind it.

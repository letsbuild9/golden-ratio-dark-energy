# Golden-ratio non-resonance and the dark sector

Code, outputs and paper for

> S. I. Borse, *Golden-Ratio Non-Resonance in Phase Space as an Organizing Hypothesis for Dark-Sector Anomalies: Pre-registered Predictions for the Dark-Energy Equation of State, the Growth Index and Early Dark Energy* (2026).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23481070.svg)](https://doi.org/10.5281/zenodo.23481070)

**Pre-registration:** the predictions were time-stamped on Zenodo on 19 October 2026, before any DESI data beyond DR2: DOI [10.5281/zenodo.23481070](https://doi.org/10.5281/zenodo.23481070). `PREDICTIONS.md` here is byte-identical to the Zenodo copy (SHA-256 in `PREDICTIONS.sha256`) and must never be edited.

## The predictions

With x = 1/φ ≈ 0.618034 (the fixed point of x = 1/(1 + x)):

| | Quantity | Value |
|---|---|---|
| P1 | w₀ | −φ/2 = −0.80902 |
| P2 | wₐ | −1/φ = −0.61803 |
| P3 | redshift where w = −1 | 1/√5 = 0.44721 |
| P4 | w(0.3), w(0.5) | −0.9516, −1.0150 |
| P5 | growth index γ | 1/φ = 0.618 |
| P6 | early-dark-energy fraction | x⁴/2 = 0.07295 |
| P7 | S8 (derived) | ≈ 0.81 |

Falsification criteria, data status and caveats are in [`PREDICTIONS.md`](PREDICTIONS.md). H₀ is not predicted.

## Repository layout

| Path | Contents |
|---|---|
| `paper/` | `paper.tex`, `paper.pdf` and figures |
| `borse_*.py` | analysis scripts (likelihoods, fits, official DESI chains, growth, mechanism tests, CAMB check, symbolic proofs) |
| `verify_paper.py` | checks that every number in the paper matches the outputs in `rerun/` |
| `rerun/` | outputs of a clean run of every script |
| `figures/` | figures produced by the scripts |
| `VERIFICATION.md` | ledger: how every number was checked, the source of every external value, corrections log |
| `PREDICTIONS.md`, `PREDICTIONS.sha256` | the locked prediction record |
| `lock_predictions.py`, `zenodo_upload/` | how the record was locked and uploaded |
| `data/` | small inputs only (see below) |

## Reproducing the results

Python 3.12; exact package versions in `requirements.txt`.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

### Data not included in this repository

The large third-party data sets are not redistributed here. Download them from the official sources:

- **Supernovae** (Pantheon+, Union3, DESY5, DES-Dovekie): from the Cobaya [`sn_data`](https://github.com/CobayaSampler/sn_data) repository, into `data/PantheonPlus/`, `data/Union3/`, `data/DESY5/` and `data/DES-Dovekie/` (same file names as in that repository).
- **Official DESI DR2 chains** (about 1.4 GB): from the DESI release at <https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params>, into `data/desi_chains/<combination>/` using the folder names in `borse_desi_chains.py`, `borse_v2_tests.py` and `borse_v2b_tests.py` (for example `cmb`, `cmb_pantheonplus`, `cmb_union3`, `cmb_desy5sn`, `binned_cmb`).

Included: DESI DR2 BAO means and covariance and the SDSS DR16 f σ8 files (from Cobaya [`bao_data`](https://github.com/CobayaSampler/bao_data)), the DESI DR1 full-shape f σ8 values, the z* calibration and the data-use acknowledgment texts.

### Run order

```text
borse_calibrate_zstar.py      (needs camb; writes data/zstar_calibration.json)
borse_checks.py
borse_fit.py --grid
borse_growth.py ; borse_growth.py --desi
borse_model.py
borse_mu.py
borse_desi_chains.py
borse_v2_tests.py             (15–20 min)
borse_v2b_tests.py
borse_binned_calibration.py
borse_energy_exchange.py
borse_camb_check.py
borse_symbolic.py
verify_paper.py               (expect: all checks passed)
```

## Citing

Please cite the Zenodo record (DOI above); GitHub's "Cite this repository" button uses `CITATION.cff`.

## License

- **Code:** MIT (see `LICENSE`).
- **Paper, prediction record and other text:** CC BY 4.0.
- **Third-party data** remain under their providers' terms; see the acknowledgments in the paper.

## Use of AI tools

The original framework was developed in conversations with Gemini (Google); numerical analysis, verification and drafting were assisted by Claude (Anthropic). The author is responsible for all content. See the paper's acknowledgements.

## Author

Siddesh I. Borse, independent researcher. ORCID [0009-0002-3306-7175](https://orcid.org/0009-0002-3306-7175).

"""Lock the prediction record with TODAY's date (never backdated) and propagate the hash everywhere.

Steps:
  1. PREDICTIONS_template.md -> PREDICTIONS.md with "Locked on: <today>" (ISO date)
  2. SHA-256 -> PREDICTIONS.sha256
  3. Copy both into zenodo_upload/
  4. Update the paper's \\lockdate and \\lockhash macros
  5. Regenerate zenodo_upload/ZENODO_FORM.txt with the date and hash

Run on the day of the Zenodo upload:  .venv/bin/python lock_predictions.py
Then rebuild the paper and the upload packages (see REVISION_PLAN.md, Phase 6).
"""

import datetime
import hashlib
import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
today = datetime.date.today()
iso, human = today.isoformat(), today.strftime("%-d %B %Y")

tpl = open(os.path.join(HERE, "PREDICTIONS_template.md")).read()
assert "{LOCK_DATE}" in tpl
pred = tpl.replace("{LOCK_DATE}", iso)
open(os.path.join(HERE, "PREDICTIONS.md"), "w").write(pred)
sha = hashlib.sha256(pred.encode()).hexdigest()
open(os.path.join(HERE, "PREDICTIONS.sha256"), "w").write(f"{sha}  PREDICTIONS.md\n")
os.makedirs(os.path.join(HERE, "zenodo_upload"), exist_ok=True)
for f in ("PREDICTIONS.md", "PREDICTIONS.sha256"):
    shutil.copy(os.path.join(HERE, f), os.path.join(HERE, "zenodo_upload", f))

tex_path = os.path.join(HERE, "paper", "paper.tex")
tex = open(tex_path).read()
tex = re.sub(r"\\newcommand\{\\lockdate\}\{[^}]*\}", lambda m: "\\newcommand{\\lockdate}{" + human + "}", tex)
tex = re.sub(r"\\newcommand\{\\lockhash\}\{[^}]*\}", lambda m: "\\newcommand{\\lockhash}{" + sha + "}", tex)
open(tex_path, "w").write(tex)

form = f"""ZENODO FORM (only needed if you upload by hand instead of running upload_to_zenodo.py)
=====================================================================================

FILES:            PREDICTIONS.md, PREDICTIONS.sha256, borse_analysis_v2.zip
DOI:              "Do you already have a DOI?" -> No -> "Get a DOI now!"
RESOURCE TYPE:    Publication -> Technical note
TITLE:            Pre-registered cosmological predictions of a golden-ratio phase-space hypothesis (locked {human})
CREATOR:          Borse, Siddesh I. | ORCID 0009-0002-3306-7175 | Independent researcher
PUBLICATION DATE: {iso}   (must be the actual upload day: do not backdate)
VERSION:          2.0
LANGUAGE:         English
LICENSE:          Creative Commons Attribution 4.0 International (CC BY 4.0)
KEYWORDS:         cosmology; dark energy; DESI; growth index; early dark energy; pre-registration; golden ratio
ACCESS:           Public, no embargo

DESCRIPTION: copy the DESCRIPTION text from upload_to_zenodo.py, replacing {{HASH}} with:
{sha}

THEN: Preview -> check -> Publish. Published files can never be changed (corrections become a new version).
"""
open(os.path.join(HERE, "zenodo_upload", "ZENODO_FORM.txt"), "w").write(form)
print(f"Locked {iso}; SHA-256 {sha}")
print("Updated: PREDICTIONS.md, PREDICTIONS.sha256, zenodo_upload/{PREDICTIONS.*, ZENODO_FORM.txt}, paper/paper.tex macros")

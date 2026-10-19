"""Create a Zenodo DRAFT with the prediction record, reserve its DOI, and upload the files.

It does NOT publish. You review the draft on zenodo.org and click "Publish" yourself.
The lock date and hash are read from PREDICTIONS.md / PREDICTIONS.sha256 (written by lock_predictions.py),
so nothing here needs editing after a re-lock.

Usage (from the Borse folder):
    .venv/bin/python zenodo_upload/upload_to_zenodo.py --dry-run     # show what will be sent, no network
    .venv/bin/python zenodo_upload/upload_to_zenodo.py --sandbox     # optional test on sandbox.zenodo.org
    .venv/bin/python zenodo_upload/upload_to_zenodo.py               # real draft on zenodo.org

The script asks for your Zenodo personal access token (typing is hidden), or reads it from the
ZENODO_TOKEN environment variable. The token is only sent to Zenodo and is never printed or saved.
"""

import argparse
import datetime
import getpass
import hashlib
import json
import os
import re
import ssl
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = ["PREDICTIONS.md", "PREDICTIONS.sha256", "borse_analysis_v2.zip"]

DESCRIPTION = """
<p>This record time-stamps quantitative cosmological predictions <b>before</b> the data that can test them:
DESI analyses including observations after 9 April 2024, Euclid growth measurements, and next-generation CMB limits on early dark energy.</p>
<p><b>Hypothesis.</b> Let x = 1/&phi; &asymp; 0.618034, the fixed point of x = 1/(1 + x). Each dark-sector parameter that differs from &Lambda;CDM is set by a ladder term
&alpha;<sub>d</sub> x<sup>d</sup>, with &alpha;<sub>d</sub> = 1 for linear or volumetric modes (d = 1, 3) and &frac12; for transverse or plasma modes (d = 2, 4).
For 1 + w<sub>0</sub>, w<sub>a</sub> and f<sub>EDE</sub> the ladder term is the departure itself; for growth it is the value of &gamma;.</p>
<p><b>Predictions:</b></p>
<ul>
<li>P1: w<sub>0</sub> = &minus;&phi;/2 = &minus;0.80902</li>
<li>P2: w<sub>a</sub> = &minus;1/&phi; = &minus;0.61803 (equivalently w<sub>a</sub> = 1 + 2w<sub>0</sub>)</li>
<li>P3: w = &minus;1 (dark-energy density maximum) at z = 1/&radic;5 = 0.44721</li>
<li>P4: w(0.3) = &minus;0.9516 and w(0.5) = &minus;1.0150</li>
<li>P5: growth index &gamma; = 1/&phi; = 0.618</li>
<li>P6: early-dark-energy fraction f<sub>EDE</sub> = x<sup>4</sup>/2 = 0.07295</li>
<li>P7 (derived): S<sub>8</sub> &asymp; 0.81</li>
</ul>
<p>H<sub>0</sub> is explicitly not predicted.</p>
<p>PREDICTIONS.md lists the falsification criteria fixed in advance, the data status at lock time, an overlap caveat
(future DESI analyses contain the DR2 data used to choose P1&ndash;P4) and honest caveats: current data measure x &asymp; 0.55&ndash;0.69,
consistent with 1/&phi; but unable to single it out; P6 is under pressure from CMB limits.</p>
<p>The SHA-256 hash of PREDICTIONS.md is <code>{HASH}</code> (also in PREDICTIONS.sha256).</p>
<p>borse_analysis_v2.zip contains the analysis code and outputs:</p>
<ul>
<li>fits to DESI DR2 BAO + Planck distance priors + Pantheon+/Union3/DESY5/DES-Dovekie supernovae;</li>
<li>the prediction's position in DESI's official DR2 posteriors and robustness checks across CMB choices;</li>
<li>the free-x fit and control constants;</li>
<li>Bayes factors;</li>
<li>pivot and fixed-redshift tests;</li>
<li>robustness tests: DESI's binned w(z), BAO without CMB (r<sub>d</sub> free), and the predicted early dark energy included;</li>
<li>momentum-exchange and energy-exchange mechanism tests;</li>
<li>the look-elsewhere test;</li>
<li>growth-index fits to SDSS and DESI DR1 f&sigma;<sub>8</sub>;</li>
<li>a CAMB cross-check;</li>
<li>exact symbolic proofs;</li>
<li>an automated check that every number in the paper matches the outputs;</li>
<li>a verification ledger.</li>
</ul>
""".strip()


def lock_info():
    pred = open(os.path.join(HERE, "PREDICTIONS.md")).read()
    m = re.search(r"\*\*Locked on:\*\* (\d{4}-\d{2}-\d{2})", pred)
    if not m:
        sys.exit("PREDICTIONS.md has no lock date; run lock_predictions.py first.")
    sha_rec = open(os.path.join(HERE, "PREDICTIONS.sha256")).read().split()[0]
    return m.group(1), sha_rec


def metadata(orcid, date, sha):
    creator = {"name": "Borse, Siddesh I.", "affiliation": "Independent researcher", "orcid": "0009-0002-3306-7175"}
    if orcid:
        creator["orcid"] = orcid
    human = datetime.date.fromisoformat(date).strftime("%-d %B %Y")
    return {
        "upload_type": "publication",
        "publication_type": "technicalnote",
        "title": f"Pre-registered cosmological predictions of a golden-ratio phase-space hypothesis (locked {human})",
        "creators": [creator],
        "description": DESCRIPTION.replace("{HASH}", sha),
        "publication_date": date,
        "version": "2.0",
        "language": "eng",
        "access_right": "open",
        "license": "cc-by-4.0",
        "keywords": ["cosmology", "dark energy", "DESI", "growth index", "early dark energy",
                     "pre-registration", "golden ratio"],
        "prereserve_doi": True,
    }


def api(method, url, token, body=None, raw=None, ctype="application/json"):
    data = raw if raw is not None else (json.dumps(body).encode() if body is not None else None)
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {token}")
    if data is not None:
        req.add_header("Content-Type", ctype)
    try:
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        ctx = None
    try:
        with urllib.request.urlopen(req, context=ctx) as r:
            txt = r.read().decode()
            return json.loads(txt) if txt else {}
    except urllib.error.HTTPError as e:
        msg = e.read().decode()[:600]
        sys.exit(f"Zenodo returned HTTP {e.code} for {method} {url.split('?')[0]}:\n{msg}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org (needs a separate sandbox account/token)")
    ap.add_argument("--dry-run", action="store_true", help="print what would be sent; no network calls")
    ap.add_argument("--orcid", default="", help="your ORCID iD, e.g. 0000-0002-1825-0097 (optional)")
    a = ap.parse_args()

    paths = [os.path.join(HERE, f) for f in FILES]
    for p in paths:
        if not os.path.exists(p):
            sys.exit(f"Missing file: {p}")
    date, sha_rec = lock_info()
    sha = hashlib.sha256(open(paths[0], "rb").read()).hexdigest()
    if sha != sha_rec:
        sys.exit(f"PREDICTIONS.md hash {sha} does not match PREDICTIONS.sha256 ({sha_rec}); the locked file was changed. Aborting.")
    today = datetime.date.today().isoformat()
    if date != today and not a.dry_run:
        sys.exit(f"Lock date {date} is not today ({today}). Re-run lock_predictions.py today so the record is not backdated.")
    print(f"PREDICTIONS.md hash verified: {sha} (locked {date})")

    md = metadata(a.orcid, date, sha)
    if a.dry_run:
        print(json.dumps({k: v for k, v in md.items() if k != "description"}, indent=2))
        print("Files:", ", ".join(f"{f} ({os.path.getsize(p) / 1024:.0f} KB)" for f, p in zip(FILES, paths)))
        return

    base = "https://sandbox.zenodo.org" if a.sandbox else "https://zenodo.org"
    token = os.environ.get("ZENODO_TOKEN") or getpass.getpass(f"Zenodo personal access token for {base} (hidden): ").strip()
    if not token:
        sys.exit("No token given.")

    dep = api("POST", f"{base}/api/deposit/depositions", token, body={})
    dep_id, bucket = dep["id"], dep["links"]["bucket"]
    print(f"Draft created: id {dep_id}")
    for f, p in zip(FILES, paths):
        with open(p, "rb") as fh:
            api("PUT", f"{bucket}/{f}", token, raw=fh.read(), ctype="application/octet-stream")
        print(f"  uploaded {f}")
    dep = api("PUT", f"{base}/api/deposit/depositions/{dep_id}", token, body={"metadata": md})
    doi = dep.get("metadata", {}).get("prereserve_doi", {}).get("doi", "(not returned)")
    print("\nDONE. The draft is NOT published yet.")
    print(f"  Reserved DOI : {doi}")
    print(f"  Review here  : {base}/uploads/{dep_id}")
    print("Open that link, check everything, then click 'Publish'.")


if __name__ == "__main__":
    main()

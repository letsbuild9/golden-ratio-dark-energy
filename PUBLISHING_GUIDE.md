# Publishing guide (independent researcher)

**Status (20 October 2026):** Step 0 done (ORCID 0009-0002-3306-7175). Step 1 done (Zenodo DOI 10.5281/zenodo.23481070, published 2026-10-19 18:05 UTC). Step 2 done (DOI and GitHub link in the paper). Step 3 done (https://github.com/letsbuild9/golden-ratio-dark-energy; later commits still to push). Step 4 in progress: arXiv account `siddeshborse` created, submission waiting for endorsement (code JEEE4W); emails in `outreach/`. Zenodo token revoked (20 Oct). To do: push to GitHub; delete the GitHub token after the push.

The order matters. **Zenodo comes first**, because its public date is what proves the predictions came before the data. Then GitHub, then arXiv.

| Step | What | Who does it | Time |
|---|---|---|---|
| 0 | ORCID iD | You | 5 min |
| 1 | Zenodo: time-stamp the predictions and get a DOI | You (Claude can guide in the browser) | 20 min |
| 2 | Put the DOI in the paper and rebuild | Claude | 5 min |
| 3 | GitHub: public code repository | You + Claude | 20 min |
| 4 | arXiv: preprint (needs an endorsement) | You (Claude can guide) | 1–14 days |
| 5 | Optional: journal submission | You | weeks |

**Safety rules for every step:** you type all passwords, tokens and 2-factor codes yourself, and you click every final **Publish** / **Submit** button yourself. Claude never sees or stores your credentials. Anything published on Zenodo or arXiv is **permanent** (it can get a new version, never be deleted).

---

## Step 0 — ORCID iD (the standard researcher ID)

1. Go to <https://orcid.org/register> and register with `siddesh.borase1997@gmail.com`.
2. Fill in: name *Siddesh I. Borse*, employment *Independent researcher*.
3. Note your iD (format `0000-000X-XXXX-XXXX`). It goes on Zenodo, arXiv and GitHub, and links all your work to you.

## Step 1 — Zenodo (do this on the lock date: 19 October 2026)

> If you do this on a later day, first run `.venv/bin/python lock_predictions.py` (it re-stamps the date; never backdate), then tell Claude so the paper is rebuilt.

1. **Log in:** <https://zenodo.org> → *Log in* → *Log in with ORCID* (or GitHub). Verify your email if asked.
2. **Link ORCID** (if you used another login): top-right menu → *Settings* → *Linked accounts*.
3. **Create a token:** top-right menu → *Applications* → *Personal access tokens* → *New token*.
   - Name: `borse-upload`
   - Scopes: tick **deposit:write** and **deposit:actions** only.
   - Copy the token (shown once). Do not paste it into the chat.
4. **Create the draft** (Terminal, in the `Borse` folder):
   ```bash
   .venv/bin/python zenodo_upload/upload_to_zenodo.py
   ```
   Paste the token when asked (typing is hidden). The script verifies the hash, uploads `PREDICTIONS.md`, `PREDICTIONS.sha256` and `borse_analysis_v2.zip`, reserves a DOI, and **does not publish**.
5. **Review the draft:** Zenodo → *Uploads* → open the draft. Check:
   - title says *locked 19 October 2026*; publication date is today;
   - creator *Borse, Siddesh I.*, add your ORCID iD to the creator;
   - 3 files present; license CC BY 4.0; resource type *Publication → Technical note*;
   - the reserved DOI (`10.5281/zenodo.NNNNNNN`). Copy it.
6. Click **Publish** (yourself). This is the time-stamp.
7. **Revoke the token** (*Applications* → delete it). Good practice; it is no longer needed.
8. Send Claude the DOI.

## Step 2 — DOI into the paper (Claude)

Claude replaces the two `[to be added]` placeholders in `paper/paper.tex` with the DOI, adds it to `VERIFICATION.md`, rebuilds `paper.pdf` and `arxiv_submission.tar.gz`, and re-runs `verify_paper.py`. The Zenodo record itself is never edited.

## Step 3 — GitHub (public code)

Standard practice: code under an MIT license, a README that says how to reproduce, and a `CITATION.cff` so GitHub shows a "Cite this repository" button.

1. **Account:** <https://github.com/signup> (if you don't have one). Turn on 2-factor authentication.
2. **Claude prepares the folder:** `README.md`, `LICENSE` (MIT, code), `CITATION.cff` (with ORCID and the Zenodo DOI), `.gitignore`, and a local git commit.
   - **Not uploaded:** `.venv/` and the large downloaded data (`data/desi_chains/` is 1.4 GB, plus the supernova files). These belong to DESI and the SN collaborations; the README gives their official download links instead.
3. **Create the repo:** github.com → *New repository* → name `golden-ratio-dark-energy` (or your choice) → Public → **do not** add a README/license (we have them) → *Create*.
4. **Push:** easiest with **GitHub Desktop** (<https://desktop.github.com>): *File → Add local repository* → choose `Borse` → *Publish repository*. Or give Claude the repo URL and Claude runs `git push` — your browser/keychain handles the login.
5. Optional: *Releases → Create release* `v2.0` to match the Zenodo record.

## Step 4 — arXiv

1. **Account:** <https://arxiv.org/user/register>. Use your real name exactly as on the paper, affiliation *Independent researcher*, and link your ORCID in the profile.
2. **Start a submission:** *Start new submission* → category **astro-ph.CO** (Cosmology and Nongalactic Astrophysics); optional cross-list: none needed.
3. **Endorsement.** New submitters usually need one. arXiv shows an endorsement code. Send it, politely and briefly, to a researcher who has published in astro-ph.CO (for example, an author you cite whose work is close to yours), with the PDF and the Zenodo DOI. Expect some to say no; one yes is enough. Claude can draft the email.
4. **Upload:** `paper/arxiv_submission.tar.gz` (paper.tex + 7 figures). arXiv compiles it with pdflatex. **Check the generated PDF preview page by page** before continuing (our build used Tectonic; only standard packages are used, so it should compile).
5. **Metadata:**
   - Title: as in the paper.
   - Authors: `Siddesh I. Borse`
   - Abstract: Claude will give you a plain-text version (under 1,920 characters).
   - Comments: `21 pages, 6 figures. Predictions pre-registered on Zenodo, DOI 10.5281/zenodo.NNNNNNN; code at <GitHub URL>`
   - License: **CC BY 4.0**.
6. **Submit** (yourself). Submissions before 14:00 US Eastern on a weekday are announced that evening. Moderators may reclassify or hold it; that is normal.

## Step 5 — Optional: journal

The **Open Journal of Astrophysics** (free, peer-reviewed, arXiv-overlay, accepts independent authors) is a natural fit once the arXiv number exists. Alternatives: *Physical Review D*, *JCAP* (publication charges may apply). Declare AI assistance as in the paper's acknowledgements.

---

## After publishing

- When DESI releases results with data after 9 April 2024, compare them against the criteria in `PREDICTIONS.md` and report the outcome, favorable or not.
- Any correction becomes a **new version** on Zenodo and arXiv, never an edit of the locked record.

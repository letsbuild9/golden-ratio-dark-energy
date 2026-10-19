# arXiv submission kit

Everything to copy into the arXiv forms. You type passwords and click **Submit** yourself.

## 1. Register (once)

<https://arxiv.org/user/register>

| Field | Value |
|---|---|
| Email | siddesh.borase1997@gmail.com |
| Username / password | your choice (you type them) |
| First name / Last name | Siddesh I. / Borse |
| Affiliation | Independent researcher |
| Country | India |
| Career status | Other (or the closest option) |
| Default category | Astrophysics → **Cosmology and Nongalactic Astrophysics (astro-ph.CO)** |

Then verify the email arXiv sends. In your arXiv profile, link ORCID **0009-0002-3306-7175**.

## 2. Start the submission

*Start New Submission* → accept the submission agreement → license **CC BY 4.0** → archive **astro-ph**, category **astro-ph.CO**.

If arXiv says you need an **endorsement**, it shows a 6-character **endorsement code**. Note it and go to section 5. (You can keep the submission saved and finish it after you are endorsed.)

## 3. Upload

File: `paper/arxiv_submission.tar.gz` (paper.tex + 7 PNG figures, about 0.75 MB). arXiv compiles it with pdflatex. The source is pure ASCII and uses only standard packages; it compiles cleanly from a fresh folder here (Tectonic). **Read arXiv's generated PDF page by page** before continuing: 21 pages, title page shows your ORCID, the DOI appears in Sec. 5 and in Data availability.

## 4. Metadata

**Title**
```
Golden-Ratio Non-Resonance in Phase Space as an Organizing Hypothesis for Dark-Sector Anomalies: Pre-registered Predictions for the Dark-Energy Equation of State, the Growth Index and Early Dark Energy
```

**Authors**
```
Siddesh I. Borse
```

**Abstract** (1,747 characters; limit 1,920)
```
We propose that the departures of the dark sector from $\Lambda$CDM are set by one number, the inverse golden ratio $x = 1/\varphi = (\sqrt5-1)/2$. The proposal has two strictly separated parts. The mathematical part is proven, apart from one labeled numerical result: distinct phase-space states never merge; $1/\varphi$ is the global attractor of the self-similar partition $x = 1/(1+x)$ and, by Hurwitz's theorem, the most badly approximable number; and golden KAM tori are found numerically to be the most robust barriers in model four-dimensional phase spaces. The physical part is a hypothesis: dark-sector parameters take ladder values $\alpha_d x^d$, with $\alpha_d = 1$ or $1/2$. Without free parameters it gives $w_0 = -\varphi/2$, $w_a = -1/\varphi$, a crossing of $w = -1$ at $z = 1/\sqrt5$, a growth index $\gamma = 1/\varphi$ and an early-dark-energy fraction $f_{\rm EDE} = x^4/2 = 0.073$. With DESI DR2 BAO, CMB data and four supernova samples, $(w_0, w_a) = (-0.809, -0.618)$ lies $0.3$--$1.0\sigma$ from the free fit, has the lowest AIC and BIC of the models tested, and lies at $0.5$--$1.0\sigma$ in DESI's official posteriors with supernovae. Of 403 pairs from the same algebraic family, it alone is within $1\sigma$ for all three main supernova samples. The agreement holds without CMB data and with the predicted early dark energy. Fitting $x$ freely gives $x = 0.55$--$0.69$: consistent with $1/\varphi$, but neighboring values fit almost as well. The growth index is $1.1$--$1.3\sigma$ from CMB-anchored $f\sigma_8$ data. An earlier claim of $H_0 = 72.7$ km/s/Mpc is withdrawn, and momentum exchange is shown not to produce the predicted growth suppression. All predictions are time-stamped before any DESI data beyond DR2.
```

**Comments**
```
21 pages, 6 figures. Predictions pre-registered on Zenodo (DOI: 10.5281/zenodo.23481070). Code and outputs: https://github.com/letsbuild9/golden-ratio-dark-energy
```

**Report number / Journal reference / DOI:** leave empty (the Zenodo DOI is for the prediction record, not this paper).

**MSC / ACM classes:** leave empty.

Preview → check → **Submit** (you). Submissions before 14:00 US Eastern (23:30 IST) on a weekday are announced that evening (US time).

## 5. Endorsement

arXiv requires new authors in a category to be endorsed by someone who has published there. Being endorsed is not peer review; it only confirms the paper is relevant to the category and of a scholarly standard.

**Finding an endorser**

- Pick 3–5 active astro-ph.CO researchers whose work you cite or that is close to yours (dark-energy phenomenology, DESI analyses, growth index, early dark energy).
- On any of their arXiv abstract pages, the link **"Which authors of this paper are endorsers?"** shows who is eligible.
- Write to **one at a time**, politely; wait about a week before the next. Expect some not to reply.

**Email template** (fill in the brackets; attach the PDF)

> Subject: arXiv endorsement request (astro-ph.CO): pre-registered dark-energy predictions
>
> Dear Dr. [Name],
>
> I am an independent researcher and am preparing my first submission to arXiv astro-ph.CO. I am writing to ask whether you would consider endorsing me. I came across your work on [specific paper/topic], which is closely related.
>
> The paper (attached, 21 pages) tests a parameter-free hypothesis for the dark-energy equation of state and the growth index against DESI DR2 BAO, CMB distance priors, four supernova samples and DESI's official posterior chains. It reports favorable and unfavorable results, quantifies the look-elsewhere effect, and withdraws earlier claims that did not survive testing. The predictions were pre-registered on Zenodo before any post-DR2 DESI data (DOI 10.5281/zenodo.23481070), and all code is public (https://github.com/letsbuild9/golden-ratio-dark-energy).
>
> If you are willing, my endorsement code is JEEE4W; it can be entered at https://arxiv.org/auth/endorse. I fully understand if you are unable to, and I would be grateful for any comments on the paper either way.
>
> Kind regards,
> Siddesh I. Borse
> Independent researcher, ORCID 0009-0002-3306-7175
> siddesh.borase1997@gmail.com

**Be ready for:** moderators can reclassify a paper (for example to gr-qc or physics.gen-ph) or put it on hold. That is normal for a first, single-author, independent submission; accept a reclassification rather than arguing.

## 6. Your endorsement status (19 October 2026)

**Code: `JEEE4W`** (endorsers enter it at https://arxiv.org/auth/endorse). Eligibility below was checked on arXiv's "Which authors can endorse?" pages for astro-ph.CO.

| Order | Endorser | Why them | Opening line for the email |
|---|---|---|---|
| 1 | **Brodie Popovic** (lead author, DES-Dovekie, arXiv:2511.07517) | You use the DES-Dovekie data, and your prediction lies 0.1σ / 0.5σ from their published w₀, wₐ | "I came across your DES-Dovekie reanalysis (arXiv:2511.07517); my paper uses the Dovekie sample, and the predicted (w₀, wₐ) lies within 0.5σ of your published constraint." |
| 2 | **Vivian Poulin** (arXiv:2505.08051; originator of the axion-like EDE model) | Prediction P6 uses his n = 3 EDE model and quotes his 2025 limits | "My paper predicts an axion-like early-dark-energy fraction in your n = 3 model and compares it with your 2025 ACT DR6 / DESI DR2 limits (arXiv:2505.08051)." |
| 3 | **Ryan Turner** (lead author, DESI DR1 peculiar velocities, arXiv:2512.03230) | You use his fσ8(0.07) and reproduce his growth index γ = 0.58 ± 0.11 | "My growth-index analysis uses your DESI DR1 peculiar-velocity measurement (arXiv:2512.03230) and reproduces your γ = 0.58 ± 0.11 as a pipeline check." |
| Backups | Rodrigo Calderón, Théo Simon (2505.08051); Niall Jeffrey, Felipe Andrade-Oliveira (2511.07517); Anthony Carr, Caitlin Ross (2512.03230) | Co-authors of the same papers, all eligible | — |

Not eligible for astro-ph.CO: the quintom-review authors (2511.19994).

Find each person's current email on their recent arXiv papers or their institution's page. Write to one person at a time and wait about a week before the next. Keep it short, attach the PDF, and do not send follow-ups more than once.

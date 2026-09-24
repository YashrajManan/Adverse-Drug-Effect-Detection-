# FDA Adverse Drug Event Analysis (Pharmacovigilance)

Disproportionality signal detection (PRR/ROR/chi-square) and a report-seriousness classifier on the
real, raw **FAERS (FDA Adverse Event Reporting System)** quarterly ASCII extracts for all of 2024
(Q1–Q4) — downloaded directly from FDA's own export servers, not the openFDA REST API. Full,
independent dual-language implementation: **Python** (data engineering, PRR/ROR/chi-square computed by
hand from raw 2×2 contingency tables, LogisticRegression + GradientBoosting classifier) and **R**
(`pvda` package for PRR/ROR, independent `glm()` logistic regression) — both run on identical
deduplicated case data, cross-checked against each other.

## The question

Using a full year of real, raw adverse-event reports, which drug-reaction pairs are reported together
disproportionately often — a candidate drug-safety signal? And can a report's own structure (age, sex,
how many drugs/roles/reactions it involves) predict whether it will be classified as medically serious,
with no clinical judgment data available in the report's fields?

## Data

Real FAERS 2024 quarterly ASCII extracts (`fis.fda.gov/content/Exports/faers_ascii_2024q{1-4}.zip`) —
seven pipe(`$`)-delimited relational tables per quarter, joined by `primaryid`/`caseid`. This project
uses DEMO (demographics), DRUG (every drug on the report + its role), REAC (coded reactions), and OUTC
(outcome codes). Real scale before deduplication: ~1.6M DEMO rows, 7.7M DRUG rows, 5.8M REAC rows.
After deduplicating each case to its latest `caseversion`: **1,484,350 distinct real 2024 case reports**.

## Method

1. **Deduplicate** — keep only each case's max-`caseversion` row before any counting (an amended case
   can appear more than once across quarters).
2. **Label seriousness** — a case is "serious" if it has ≥1 row in OUTC (FAERS' own regulatory
   definition), not a trusted pre-computed flag. 59.5% of deduplicated cases are serious.
3. **Signal detection** — for the top-150 most-reported drugs × top-150 most-reported reactions, build
   a 2×2 contingency table (a = both present, b = drug only, c = reaction only, d = neither) over the
   full 1,484,350-case universe, then compute PRR, ROR, and chi-square. **Deliberate scope choice: all
   four drug roles (primary suspect / secondary suspect / concomitant / interacting) are included**, not
   just suspect-only (standard convention) — a disclosed departure, reasoned through in
   `PROJECT_NARRATIVE.md`.
4. **Seriousness classifier** — logistic regression (Python + R) and gradient boosting (Python) predict
   `serious` from age, sex, the four role counts, reaction count, and total drug count.

## Real results

**Top signal: PAXLOVID × "Disease recurrence" — PRR = 203.5, ROR = 344.0** (matches the well-documented
"COVID-19 rebound" effect). **Twelve independent DMARD/biologic drugs** (leflunomide, sulfasalazine,
several TNF-inhibitors, JAK inhibitors, etc.) all pair with **"Hepatic enzyme increased"** as a top
signal (PRR ≈ 35–70) — consistent with the liver-function-monitoring requirements on these drugs'
actual labels. Several other top signals are **confounding by indication** (e.g. Dupixent ×
"Dermatitis atopic," Vedolizumab × "Colitis ulcerative") — the drug's own indication re-reported as an
"event," a known limitation of disproportionality analysis, not a modeling error.

**Classifier:** Python LogReg AUC = 0.658, independent R `glm()` AUC = 0.659 (near-exact cross-language
agreement), GradientBoosting AUC = 0.704 — real signal above chance, well short of diagnostic-grade.

**A real cross-language bug was found and fixed during this build:** R's PRR/ROR denominator initially
used the wrong case count (the classifier's age/sex-filtered subset, 863,150, instead of the true full
deduplicated universe, 1,484,350), which silently deflated every R-side PRR/ROR relative to Python's.
Caught by comparing the same pair's numbers across both languages, not by either script erroring — full
account, and the fix, in `PROJECT_NARRATIVE.md`.

Full results: `results_py/`, `results_R/` (`signal_scores_*top200.csv` = top 200 pairs by PRR;
`feature_importance_*`, `model_metrics_*` = classifier outputs).

## Files

- `fda_adverse_event_analysis.ipynb` — Python: fetch, dedup, seriousness label, PRR/ROR/chi-square,
  classifier.
- `fda_adverse_event_analysis.R` — R twin: reads the Python-exported bridge files, runs `pvda`-based
  PRR/ROR and an independent `glm()` classifier.
- `results_py/`, `results_R/` — real output CSVs from both languages.

## Caveats (real, not boilerplate)

Spontaneous reports have no reliable exposure denominator and no confirmed causality per report — PRR/
ROR flag disproportionate co-reporting, a hypothesis worth clinical follow-up, never proof of a real
adverse drug reaction on their own. Reporting itself is biased (media attention, a drug's newness, and
litigation all inflate report volume independent of true harm). This project's all-drug-roles scope is
a deliberate departure from suspect-only convention and likely understates true suspect-drug signal
strength rather than inflating it (reasoned through, not directly measured against a suspect-only
rerun). Chi-square is near-automatically "significant" at FAERS' report volume — PRR/ROR magnitude and
the minimum case-count floor (`a ≥ 3`) are what actually separate real signal from noise here, not the
p-value.

# FDA Adverse Drug Event Analysis (Pharmacovigilance)

Disproportionality signal detection (PRR/ROR/chi-square) and a report-seriousness classifier on the
real, raw **FAERS (FDA Adverse Event Reporting System)** quarterly ASCII extracts for all of 2024
(Q1–Q4) — downloaded directly from FDA's own export servers, not the openFDA REST API. Full,
independent dual-language implementation: **Python** (data engineering, PRR/ROR/chi-square computed by
hand from raw 2×2 contingency tables, LogisticRegression + GradientBoosting classifier) and **R**
(`pvda` package for PRR/ROR, independent `glm()` logistic regression) — both run on identical
deduplicated case data, cross-checked against each other.

## Aim

Using a full year of real, raw adverse-event reports, which drug-reaction pairs are reported together
disproportionately often — a candidate drug-safety signal? And can a report's own structure (age, sex,
how many drugs/roles/reactions it involves) predict whether it will be classified as medically serious,
with no clinical judgment data available in the report's fields?

## Objective

Build a real, hand-computed pharmacovigilance signal-detection pipeline (PRR, ROR, chi-square) directly
from FDA's raw quarterly FAERS extracts, label report seriousness from FAERS' own regulatory outcome
codes rather than a pre-computed flag, train a seriousness classifier from purely structural report
features, and cross-validate every real result independently in Python and R.

## Data fetch

Real FAERS 2024 quarterly ASCII extracts (`fis.fda.gov/content/Exports/faers_ascii_2024q{1-4}.zip`) —
seven pipe(`$`)-delimited relational tables per quarter, joined by `primaryid`/`caseid`, downloaded
directly from FDA's own export servers (not the openFDA REST API). This project uses DEMO
(demographics), DRUG (every drug on the report + its role), REAC (coded reactions), and OUTC (outcome
codes). Real scale before deduplication: ~1.6M DEMO rows, 7.7M DRUG rows, 5.8M REAC rows.

## Data describe

After deduplicating each case to its latest `caseversion` (an amended case can appear more than once
across quarters, so only the max-`caseversion` row per case is kept before any counting):
**1,484,350 distinct real 2024 case reports**. Seriousness is labeled from FAERS' own regulatory
definition — a case is "serious" if it has ≥1 row in OUTC, not a trusted pre-computed flag — giving
59.5% of deduplicated cases labeled serious.

## Methods / Workflow — what we did

1. **Deduplicate** — keep only each case's max-`caseversion` row before any counting.
2. **Label seriousness** — from OUTC row presence (FAERS' own regulatory definition), not a
   pre-computed flag.
3. **Signal detection** — for the top-150 most-reported drugs × top-150 most-reported reactions, build
   a 2×2 contingency table (a = both present, b = drug only, c = reaction only, d = neither) over the
   full 1,484,350-case universe, then compute PRR, ROR, and chi-square by hand from the raw counts.
   **Deliberate scope choice: all four drug roles (primary suspect / secondary suspect / concomitant /
   interacting) are included**, not just suspect-only (standard convention) — a disclosed departure,
   reasoned through in `PROJECT_NARRATIVE.md`.
4. **Seriousness classifier** — logistic regression (Python + R) and gradient boosting (Python) predict
   `serious` from age, sex, the four role counts, reaction count, and total drug count.
5. **Cross-language validation** — repeat PRR/ROR (via R's `pvda` package) and the logistic classifier
   (via R's `glm()`) independently in R on identical case data, comparing every real number against the
   Python results rather than trusting either language alone.

## Results

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

## Biology interpretation of results

The PAXLOVID × "Disease recurrence" signal (PRR = 203.5) is a real, independently corroborated
pharmacological finding — it matches the well-documented "COVID-19 rebound" phenomenon reported in
clinical follow-up studies, giving external validation that the hand-computed PRR/ROR pipeline
recovers a genuine, known signal rather than noise. The twelve-drug DMARD/biologic cluster pairing
with "Hepatic enzyme increased" is mechanistically coherent: these drug classes (TNF-inhibitors, JAK
inhibitors, and related immunomodulators) carry real, label-mandated liver-function-monitoring
requirements, so elevated hepatic-enzyme reporting across independently-acting drugs converging on the
same reaction term is consistent with a shared class-level hepatotoxicity monitoring practice rather
than coincidence. By contrast, signals like Dupixent × "Dermatitis atopic" and Vedolizumab × "Colitis
ulcerative" are a textbook pharmacovigilance artifact called **confounding by indication** — the
disease a drug treats gets re-reported as if it were a caused "event," inflating disproportionality
scores for reasons that have nothing to do with drug harm. Distinguishing genuine safety signal
(Paxlovid, the DMARD cluster) from indication-confounding (Dupixent, Vedolizumab) in the same ranked
output is the actual clinical-reasoning work of pharmacovigilance, and both kinds of signal are
reported honestly here rather than only keeping the "clean-looking" ones. The classifier's modest AUC
(0.658-0.704) is itself a biologically sensible result: a report's own structural metadata (age, sex,
drug/reaction counts) carries real but limited information about medical seriousness, because
seriousness is fundamentally a clinical judgment that these structural fields only partially proxy for.

## Learning through project

A cross-language discrepancy caught by comparing the *same real number* across two independent
implementations — not by either script crashing — is one of the most reliable ways to find a real bug:
here, R's PRR/ROR was silently computed against the wrong denominator (the classifier's filtered
subset rather than the true full case universe), and it was the mismatch against Python's numbers, not
an error message, that surfaced it. A disproportionality-analysis scope choice (here: including all
four drug roles rather than suspect-only) is not a neutral technical default — it's a real
methodological decision that shifts results in a knowable, disclosed direction (likely understating
true suspect-drug signal strength), and naming that direction explicitly is what keeps a departure from
convention honest rather than hidden. Recovering a well-documented real-world signal (Paxlovid/COVID
rebound) from a from-scratch PRR/ROR implementation is meaningful external validation that the pipeline
itself is computing something real, independent of any of the project's own internal cross-checks.
Reporting confounded signals (confounding by indication) alongside genuine ones, and naming which is
which, is more scientifically honest than filtering the output down to only the "clean" results.

## Limitations

Spontaneous reports have no reliable exposure denominator and no confirmed causality per report — PRR/
ROR flag disproportionate co-reporting, a hypothesis worth clinical follow-up, never proof of a real
adverse drug reaction on their own. Reporting itself is biased (media attention, a drug's newness, and
litigation all inflate report volume independent of true harm). This project's all-drug-roles scope is
a deliberate departure from suspect-only convention and likely understates true suspect-drug signal
strength rather than inflating it (reasoned through, not directly measured against a suspect-only
rerun). Chi-square is near-automatically "significant" at FAERS' report volume — PRR/ROR magnitude and
the minimum case-count floor (`a ≥ 3`) are what actually separate real signal from noise here, not the
p-value.

## Reproduce

**Python:** run `fda_adverse_event_analysis.ipynb` top to bottom (installs its own dependencies).
Requires downloading the real FAERS 2024 quarterly ASCII extracts from FDA's export servers.

**R:** run the Python notebook first (produces the bridge CSVs), then run
`fda_adverse_event_analysis.R`. Requires the `pvda` package for PRR/ROR.

## Tech

`Python` (pandas, scikit-learn: LogisticRegression + GradientBoosting) · `R` (`pvda`, `glm()`) ·
FAERS raw ASCII quarterly extracts · hand-computed PRR/ROR/chi-square from 2×2 contingency tables

## Files

```
fda_adverse_event_analysis.ipynb   # Python: fetch, dedup, seriousness label, PRR/ROR/chi-square, classifier
fda_adverse_event_analysis.R       # R twin: reads the Python-exported bridge files, pvda-based PRR/ROR + independent glm() classifier
results_py/, results_R/            # real output CSVs from both languages (signal_scores_*, feature_importance_*, model_metrics_*)
```

## License

All rights reserved — see `LICENSE`. This repository is public for portfolio/demonstration purposes
only; no permission is granted to copy, modify, or reuse any part of it.

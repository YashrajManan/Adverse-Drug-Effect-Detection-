# =============================================================================
# #22B -- FDA Adverse Drug Event Analysis -- R twin
# Reads the deduplicated bridge files the Python notebook (fda_adverse_event_analysis.ipynb)
# exports to data_py/, rather than re-implementing the multi-million-row fetch+dedup
# logic a second time in R -- a disclosed cross-language bridge (see PROJECT_NARRATIVE.md).
#
# This R twin focuses on the signal-detection step using a real, dedicated
# pharmacovigilance package -- `pvda`, which implements PRR/ROR natively.
#
# NOTE: the original plan named PhViD, the standard R disproportionality-analysis
# package. PhViD was archived/removed from CRAN on 2024-05-18 (unmaintained) and can
# no longer be installed on a current R version -- discovered live when the install
# failed during the actual build. `pvda` is the actively maintained CRAN replacement
# used instead; its function signatures are genuinely different from PhViD's.
# =============================================================================

# ---- packages + paths ----
# install.packages(c("tidyverse", "pvda", "pROC"))  # run once if not already installed
library(tidyverse)
library(pvda)
library(pROC)

DATA_DIR <- "data_py"
RESULTS_DIR <- "results_R"
dir.create(RESULTS_DIR, showWarnings = FALSE)

# =============================================================================
# STEP 1 -- Load the Python-exported bridge file
# model_features.csv is the deduplicated, case-level feature table (age, sex,
# drug-role counts, reaction count, the regulation-grounded `serious` label)
# that the Python notebook already built from the raw FAERS extracts.
# =============================================================================
model_features <- read.csv(file.path(DATA_DIR, "model_features.csv"), stringsAsFactors = FALSE)

cat("rows:", nrow(model_features), "\n")
cat("fraction serious:", mean(model_features$serious), "\n")

# =============================================================================
# STEP 2 -- Build the 2x2 cell counts (a, b, c, d) for every co-occurring
#           drug x reaction pair
# case_drug_long.csv (caseid, drug) and case_reac_long.csv (caseid, reaction)
# are two more Python-exported bridge files, restricted to the same top-150
# drug/reaction shortlists Python used. Inner-joining them on caseid gives
# every (drug, reaction) pair that co-occurs in at least one case; counting
# gives `a` (cases with both). Each drug's/reaction's own total case count
# (a+b, a+c) comes from separately counting the two long tables.
#
# total_cases is read from total_cases.csv -- the TRUE full deduplicated case
# universe Python's own Step 5 computed (demo_dedup's case count), NOT
# nrow(model_features), which is narrowed by dropna(age, sex) for the
# classifier and is the wrong (smaller) denominator for signal detection --
# a real bug found and fixed during this project's build (see
# PROJECT_NARRATIVE.md Section 7).
# =============================================================================
case_drug_long <- read.csv(file.path(DATA_DIR, "case_drug_long.csv"), stringsAsFactors = FALSE)
case_reac_long <- read.csv(file.path(DATA_DIR, "case_reac_long.csv"), stringsAsFactors = FALSE)

total_cases <- read.csv(file.path(DATA_DIR, "total_cases.csv"))$total_cases

drug_totals <- case_drug_long %>% count(drug, name = "n_drug")
reac_totals <- case_reac_long %>% count(reaction, name = "n_event")

drug_reaction_pairs <- inner_join(case_drug_long, case_reac_long, by = "caseid")
pair_counts <- drug_reaction_pairs %>% count(drug, reaction, name = "a")

MIN_CASES <- 3
pair_counts <- pair_counts %>%
  filter(a >= MIN_CASES) %>%
  left_join(drug_totals, by = "drug") %>%
  left_join(reac_totals, by = "reaction") %>%
  mutate(
    b = n_drug - a,
    c = n_event - a,
    d = total_cases - a - b - c
  )

# =============================================================================
# STEP 3 -- Run PRR/ROR signal detection via pvda
# pvda::prr()'s real internal formula (confirmed from its source,
# R/lower_level_disprop_analysis.R) is obs / (n_drug * (n_event_prr / n_tot_prr)),
# which only reduces to the standard Evans et al. PRR = [a/(a+b)]/[c/(c+d)] when
# n_event_prr and n_tot_prr are the BACKGROUND counts (c, c+d) excluding the
# drug -- not the pair's marginal totals. pair_counts already has c and d from
# Step 2. pvda::ror(a, b, c, d) takes the classic 2x2 cells directly.
# =============================================================================
prr_result <- pvda::prr(obs = pair_counts$a, n_drug = pair_counts$a + pair_counts$b,
                         n_event_prr = pair_counts$c, n_tot_prr = pair_counts$c + pair_counts$d)
ror_result <- pvda::ror(a = pair_counts$a, b = pair_counts$b, c = pair_counts$c, d = pair_counts$d)

signal_df_R <- pair_counts %>%
  mutate(PRR = prr_result$prr, ROR = ror_result$ror) %>%
  arrange(desc(PRR))

write.csv(signal_df_R, file.path(RESULTS_DIR, "signal_scores_R.csv"), row.names = FALSE)

# =============================================================================
# STEP 4 -- Logistic regression on the seriousness label
# Same feature set as Python's LogReg (age, sex_encoded, 4 drug-role counts,
# n_reactions, n_drugs_total), fit independently via base R glm() -- a
# cross-language cross-check of the classifier's AUC and coefficients.
# Note: n_drugs_total is an exact linear combination of the four role counts
# (n_role_PS + n_role_SS + n_role_C + n_role_I), so glm()'s unregularized MLE
# fit aliases its coefficient to NA (perfect collinearity) -- a real, expected
# difference from Python's L2-regularized sklearn fit, not a bug (see
# PROJECT_NARRATIVE.md Section 7).
# =============================================================================
model_features$sex_encoded <- ifelse(model_features$sex == "M", 0, ifelse(model_features$sex == "F", 1, 2))

set.seed(42)
n <- nrow(model_features)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))
train_df <- model_features[train_idx, ]
test_df <- model_features[-train_idx, ]

glm_fit <- glm(serious ~ age + sex_encoded + n_role_PS + n_role_SS + n_role_C + n_role_I +
                 n_reactions + n_drugs_total, data = train_df, family = binomial)

test_pred <- predict(glm_fit, newdata = test_df, type = "response")
roc_obj <- roc(test_df$serious, test_pred)
cat("R logistic regression AUC:", auc(roc_obj), "\n")

# =============================================================================
# STEP 5 -- Export R-side results
# =============================================================================
coef_df <- data.frame(feature = names(coef(glm_fit)), coefficient = coef(glm_fit))
write.csv(coef_df, file.path(RESULTS_DIR, "feature_importance_logreg_R.csv"), row.names = FALSE)

auc_df <- data.frame(model = "LogReg_R", auc = as.numeric(auc(roc_obj)))
write.csv(auc_df, file.path(RESULTS_DIR, "model_metrics_R.csv"), row.names = FALSE)

cat("done\n")

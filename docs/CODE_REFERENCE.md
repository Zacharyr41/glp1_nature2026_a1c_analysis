# Code reference — Part II (GLP-1 outcome analysis)

Every script, what it reads, what it writes, and the exact model it fits. Line numbers refer to
the committed files. Where several scripts share a pattern it is described once in
`ARCHITECTURE.md` §3 and only the differences are listed here.

---

## A. Shared modules

### `persistence.py` — the persistence-of-therapy rule (327 lines)

| Symbol | Line | Role |
|---|---|---|
| `MENTION_VALUES = (1.0, 2.0)`, `STOP_VALUE = 2.0` | 122, 125 | Values of `glp1_event_for_adherance` that count as evidence / explicit stop |
| `REQUIRED_COLUMNS` | 128 | `patient_id`, `glp1_event_for_adherance`, `glp1_days_from_baseline` |
| `normalize_patient_id(values)` | 133 | Stripped strings; the only place IDs are canonicalised (fixes a dtype-mismatch defect from review) |
| `mention_timeline(df)` | 159 | `{pid: [(value, day), …]}` sorted, non-negative days only |
| `censor_day_for_patient(pairs, gap_days)` | 189 | Per-patient last day within persistence: walk mention days, stop at the first inter-mention interval `> gap` (deliberately not `> gap+1`, lines 60–73), coverage ends `prev + gap`; capped at `stop_day − 1` |
| `censor_days(source, gap_days)` | 248 | Map over patients |
| `adherence_flags(df, gap_days_list)` | 272 | Adds `adherence_{g}` per row: 1 iff the most recent mention at or before the row is within `g` days **and** the row precedes any stop day |

The rule, stated once (lines 19–33): an observation on day *d* is within persistence for gap *g*
when some mention day *m ≤ d* has *d − m ≤ g*, and *d < stop_day*. step1 uses `adherence_flags`
(per row); step0a uses `censor_days` (per patient) because the attrition table counts calendar
coverage, not retained observations. The two readings agree given the day-0 anchor that step1
guarantees. The study data never contain value 2, so the stop branch is unexercised.

### `model_spec.py` (150 lines)

- `A1C_ORDER = ["Normal Glycemia", "Prediabetes", "Type 2 Diabetes", "Poorly Controlled Diabetes"]`, `REF_CATEGORY = A1C_ORDER[0]`.
- `load_spline_df(config_json, key="best_df")` — raises `FileNotFoundError` / `KeyError` / `ValueError`; **never defaults**, because the df *is* the model specification (ten scripts used to fall back to 3 silently).
- `enforce_a1c_order(df, column, order, ref_category, require_ref, context)` — ordered categorical; raises if the column is absent, no value matches, or the reference level is missing.

### `covariates.py` (97 lines)

`filter_estimable(candidates, df, *, exclude, context, logger)` keeps a covariate iff it is not
excluded, present, and has `nunique(dropna=True) > 1`; logs every drop with the reason. Replaces
twelve inline list comprehensions that dropped silently.

### `analysis_config.py` (87 lines)

`CONF_LEVEL` from env `CI_CONF_LEVEL` (default 0.95, validated in (0,1));
`z_critical(level=None) = norm.ppf(0.5 + level/2)`; `conf_level_label()`. Every script sets
`Z_CRIT = z_critical()` at import, replacing 61 hard-coded `1.96`s.

### `gap_grids.sh`

`GAPS_PRIMARY=(30 60 90 120 150 180 365 730)`, `GAPS_WITH_548` (adds 548),
`GAPS_STRUCTURED` (= primary), `GAP_PRIMARY=120`, `MAX_DAYS=730`. Sourced by all four runners.

---

## B. Cohort construction

### `step1_prepare_analysis_dataset.py` — weight (905 lines)

**CLI**: `--input-csv` (default `root_data/step8g_with_unstructured_flags.csv`), `--outdir`
(`output/step1_prepare_analysis_dataset`), `--max-days` (730), `--adherence-gaps` (nargs, default `90`), `--log-level`.

**`load_and_prepare`** (line 335), in order:

1. `_fill_exclusion_flag` for `pregnant_during_glp1`, `bariatric_surgery_flag` (absent → 0); drop rows with 1.
2. Cohort filter: `baseline_glp1_brand_final` matches `\bozempic\b|\bwegovy\b|\bmounjaro\b|\bzepbound\b` **or** `baseline_glp1_ingredient_final` matches `\bsemaglutide\b|\btirzepatide\b`; fallback `baseline_glp1` with either pattern (lines 39–50, 353–372). `GLP1_BRAND_NAMES` / `GLP1_INGREDIENT_NAMES` are the single source for these patterns.
3. Require non-null `baseline_weight_final` **and** `baseline_a1c_final` (both pipelines need both baselines).
4. `event_date = date` (else `medication_date`); compute `pct_weight_change` if absent; drop rows with missing outcome.
5. `_fill_covariates` (Unknown / −1 fills), `_derive_age_groups` (decades; 20–39 vs 40+; 20–49 vs 50+).
6. Ordered categorical `baseline_a1c_category` with `Unknown`, then drop `Unknown`.
7. `_ensure_glp1_evidence_and_metformin`: `glp1_evidence_therapy = is_glp1_event_med | is_glp1_event_clinical | event_type=="medication"`; `metformin_with_glp1_baseline` from `medication_history` JSON within ±180 d of baseline, computed once per patient, **fails closed** on unparseable input (lines 155–232).
8. `days_from_baseline = (event_date − baseline_glp1_date).days`, keep `0 ≤ d ≤ max_days`.

**`run_for_gaps`** (line 702): `_ensure_baseline_rows` (synthetic day-0 row per patient lacking one, `pct_weight_change = 0.0`, `glp1_event_for_adherance = 1`, `baseline_carried_to_day0 = 1`, lines 430–528) → helper distance columns → `persistence.adherence_flags(df, gaps)` → for each gap, drop rows at or after the patient's first `adherence_{g} == 0` day (`_rows_before_first_lapse`) → write `analysis_ready_gap{g}.csv` with the `keep_cols` list (lines 759–808). Raises if the marker column is missing from the output.

### `step1_prepare_analysis_dataset_a1c.py` — HbA1c (603 lines)

Imports the filters and helpers from the weight module. Differences:

- Requires `baseline_a1c_final` and at least one of `a1c_value` / `abs_a1c_change`; derives the other (`abs_a1c_change = a1c_value − baseline_a1c_final`).
- Missing-outcome rows are not dropped at load; instead patients are dropped unless they have ≥1 row with `days_from_baseline > 0` and a non-missing `a1c_value` (lines 159–169).
- `_ensure_baseline_rows_a1c` is the structural twin but creates no rows on the study data (HbA1c appears only on measurement dates; 38 % of patients have an observed day-0 value, the rest contribute no day-0 observation). `baseline_carried_to_day0` is emitted and sums to zero.
- Output `analysis_ready_a1c_gap{g}.csv`; see `DATA_SCHEMA.md` §2 for the column differences.

### `structured_only/step0_prefilter_raw.py` (110 lines)

Run on the raw CSV before step1 for the structured-only sensitivity cohort: where
`glp1_has_structured == 0` set `glp1_event_for_adherance` to 0; where `weight_has_structured == 0`
null `weight_in_pounds_final` and `pct_weight_change`; where `a1c_has_structured == 0` null
`a1c_value`, `abs_a1c_change`, `A1C_SOURCE`. Missing flag columns default to 0 (everything treated
as unstructured). Reads the whole file into memory.

### `structured_only/gap120/step0_filter_to_structured.py` (340 lines)

Alternative (post-step1) structured-only filter: keep rows with the flag set or `days_from_baseline == 0`,
then patients with ≥1 post-baseline flagged row. Not used by the production runner.

---

## C. Descriptives

### `step0_analysis_population_table.py` (435 lines)

Table 3. Reads a weight analysis-ready CSV plus the matching A1c file
(`usecols=["patient_id","abs_a1c_change"]` and `["patient_id","baseline_a1c_final"]`), optionally
`root_data/step8f.csv` for brand/ingredient. One row per patient (first row after sorting by
`baseline_glp1_date`). Rows: Mean/Median/Mode/SD(ddof=1)/Range of `age`, `baseline_weight_final`,
`height_in_inches_final`, `baseline_bmi_final`, `baseline_a1c_final`; counts of
`baseline_bmi_final_category`, `age_group`, `gender`, `race`, brand, ingredient; achievement flags
`min(pct_weight_change) ≤ −5/−10/−15` and `min(abs_a1c_change) ≤ −0.5…−2.5`. Output
`step0_baseline_table_gap{g}.csv` sorted alphabetically by `variable, level` (so clinical ordering
is lost; the published table was re-ordered downstream). CLI: `--input-csv`, `--outdir`,
`--adherence-gap-days` (90), `--log-level`.

### `step0a_samplesize_analysis.py` (336 lines)

Supplementary Table 1. For each `analysis_ready_gap*.csv` in `--analysis-dir`, takes the cohort
(`patient_id`), rebuilds each patient's mention timeline from the **untrimmed source**
(`--source-csv`, chunked 1.5 M rows; rows kept when `pct_weight_change` is non-missing and
`0 ≤ date − baseline_glp1_date ≤ 730`), adds a synthetic `(1, 0)` mention for patients without an
observed day-0 row, and counts patients with `persistence.censor_days(...) ≥ 30·m` for months
0–18. `--censor-from gapfiles` reproduces the submitted (double-censored) construction. Output
`samplesize_by_month.csv` (`month_number, days_from_baseline, gap_30 … gap_730`).

---

## D. Trajectory models

### `step2_select_spline_df.py` / `_a1c.py`

Spline df selection. Formula `"{y} ~ bs(days_from_baseline, df={df}) + <covariates>"` for df in
`--df-grid` (default `3,4,5,6`), covariates = those of `age_group, gender, baseline_a1c_category,
baseline_bmi_final_category, race, weight_change_med` (A1c also `metformin_with_glp1_baseline`)
with more than one level; GEE Gaussian / Independence / `groups=patient_id`; `res.qic()` →
pick minimum **QICu**. Writes `df_selection*.csv` and `model_config*.json` (`best_df`, `best_qicu`).
Production chose 3 for both outcomes. The A1c version has no gap routing.

### `step3_fit_gee_baseline.py` / `_a1c.py`

Main-effects GEE at `best_df` (read with plain `json.load`, not `load_spline_df`):
`"{y} ~ bs(days_from_baseline, df) + age_group + gender + baseline_a1c_category +
baseline_bmi_final_category + race + weight_change_med + metformin_with_glp1_baseline"` (single-
level covariates dropped with a log line; `--drop-covariate` to drop by name). Outputs
`coefficients.csv` (`term, estimate, std_error, ci_lower, ci_upper, p_value`), `model_summary.txt`,
`model_config_used.json`. The README's "`_nomet`" refers to running with metformin excluded; there
is no flag or directory named `_nomet` — later scripts simply omit the term from their lists.

### `conf_int/gap_120/step4_predictive_plots.py` (weight) / `_a1c.py` — Fig. 3a–b, ED Table 3

Per baseline-A1c category (skip if `n_obs < --min-nobs` 100):
`"pct_weight_change ~ bs(days_from_baseline, df) + age_group + gender + baseline_bmi_final_category + race [+ weight_change_med]"`
(A1c: `abs_a1c_change ~ …`; `baseline_a1c_category` is in the weight candidate list but is always
dropped as single-valued inside a stratum). Prediction grid 0..min(548, max day) step 14 at modal
covariates; delta-method CI. Writes per-category coefficient/covariance/overview CSVs,
uncentered and centered trajectory PNGs, grouped overlays (Fig. 3a/b), and
`tables/predicted_means_by_day_*.csv` (ED Table 3 uses days 84, 183, 280, 365). A `supported`
mask (≥40 obs and ≥30 patients within ±28 d, ≥10 % of stratum, CI width ≤ 5 pp or ≤ 3× median)
is written but not applied to the plots. Fixed y-limits (−10, 4) weight, (−3, 1) A1c. Env knobs
`PLOT_CI_*`, `PLOT_SHOW_CI`, `PLOT_TRAJECTORY_CI_STYLE`.

### `conf_int/gap_120/step4_observed_summary_plots.py`

Model-free companions: 90-day-bin and ±30-day-month means ± `z·sd/√n` of `pct_weight_change`,
`weight_in_pounds_final`, `a1c_value`, `abs_a1c_change` by category, anchored to the first bin.
CLI `--weight-csv`, `--a1c-csv`, `--outdir`, `--max-days` (548; runner passes 730), `--bin-width` (90).

### `step5_forest_contrasts_weight.py` / `_a1c.py` — ED Table 2 (GEE block)

One combined model with interaction:
`"{y} ~ bs(days_from_baseline, df) * baseline_a1c_category + age_group + gender + baseline_bmi_final_category + race [+ weight_change_med]"`.
`enforce_a1c_order` fixes the reference. Predictions at `--time-days` (weight
`90,180,270,365,450,548,630,730`; A1c `90,180,365,548`) per category at modal covariates;
contrast `(X_cat − X_ref)·β` with delta-method SE from the robust covariance. Outputs
`gee_combined_coefficients.csv`, `forest_predictions.csv`, `forest_contrasts_vs_ref.csv`, forest
PNGs per day and grouped 12 m / 18 m. Optional `step7_adherence_counts/adherence_counts.csv`
(`scope, day, n_unique_window`) only annotates titles.

### `conf_int/gap_120/step6_stratified_by_covariates_weight.py` / `_a1c.py` — Fig. 3c–d

For each stratifier in `age_group, age_group_20_39_vs_40_plus, age_group_20_49_vs_50_plus,
gender, race, metformin_with_glp1_baseline, baseline_bmi_final_category, glp1_user_group` (missing
ones skipped) and each subgroup: `"{y} ~ bs(days_from_baseline, df) + <covariates minus the
stratifier>"` with candidates `baseline_a1c_category, baseline_bmi_final_category, age_group,
gender, race, weight_change_med`. Grouped overlays per stratifier (Fig. 3c sex, 3d age).
`stratified_summary_counts.csv`. A failing stratum is logged and skipped.

### `step6b_stratified_contrasts_weight.py` / `_a1c.py`

Per stratifier, separate GEE per subgroup; difference of predicted means at `--time-days`
(default 365) vs a reference subgroup (youngest age bin, `20-39`/`20-49`, else largest N);
`se_diff = sqrt(se_ref² + se_sub²)`. Output `contrasts_{strat_var}.csv`. `--min-nobs` 100,
`--max-days` 548.

### `step6c_stratified_forest_plots.py` (A1c) / `_weight.py`

Per stratifier × subgroup: (a) interaction model `bs(days) * baseline_a1c_category + covariates`
with contrasts vs Normal Glycemia, (b) main-effects model `bs(days) + baseline_a1c_category +
covariates` → `subgroup_predictions_main.csv`; per-day forests, grouped 12 m / 18 m forests
(`day ∈ {365}` and `{548, 547, 540}`), trajectories (weight: gender only). Known defects in the
A1c version: `globals().get("a1c_cats", [])` is always empty (line 1005), a `zip(y_ticks,
texts_map.items())` label bug (line 1092), and the all-days PNG is overwritten by the grouped PNG.
These affect supplementary figures only.

### `conf_int/gap_120/step6cc_3waystrat_covariates_weight.py` / `_a1c.py`

Sex × age (<40 / 40+, from `age` or `age_group_20_39_vs_40_plus`) strata:
`"{y} ~ bs(days_from_baseline, df=3) + baseline_a1c_category + <remaining covariates>"`.
Two defects: `getattr(args, "spline-df", 3)` never resolves the `--spline-df` argument, so df is
always 3 (which equals the selected value); the A1c script contains a duplicated nested stratum
loop (lines 367–536) that refits everything per outer iteration and leaves legend counts at 0.
Figures only; no CSV.

### `step6d_groups_by_glp1.py` (928 lines)

Semaglutide-only vs tirzepatide-only within each glycemic category:
`"{y} ~ bs(days_from_baseline, df) * glp1_group2 * baseline_a1c_category + covariates"`, rows
truncated to `[0, --truncate-days 548]`. `glp1_user_group` is merged from `root_data/step8f.csv`
when absent (it is absent from step1 outputs); labels containing `sema`/`tirz` and not `switch`
map to the two groups. Trajectories, differences and forests at days 90–548; all curves centred
at the day-0 prediction (CIs shifted, not re-derived).

---

## E. Time-to-event

### `step8_survival_time_to_weight_loss.py` / `step8_survival_time_to_a1c_drop.py`

Per patient and threshold (`WEIGHT_THRESHOLDS = [−5, −10, −15]` %;
`A1C_THRESHOLDS = [0.5, 1.0, 1.5, 2.0]` pp): event = first day with `pct_weight_change ≤ thr` /
`abs_a1c_change ≤ −thr`, else censored at the last observed day; `MAX_FOLLOWUP_DAYS = 540`
converts any time beyond 540 to a censoring at 540. Writes the event table and a product-limit
KM summary with Greenwood / log(−log) median CIs. Required columns: `patient_id`,
`days_from_baseline`, outcome, `baseline_a1c_category`.

### `conf_int/gap_120/step8_survival_plots_and_cox.py` — Fig. 3e, ED Table 2 (Cox block)

`KaplanMeierFitter` per category (x capped 548, lifelines' own 95 % band);
`CoxPHFitter().fit(X, "time_days", "event")` with `pd.get_dummies(drop_first=True)` over
`baseline_a1c_category` (ordered, reference Normal Glycemia), `age_group`, `gender`,
`baseline_bmi_final_category`, `race`, plus numeric `weight_change_med`;
`proportional_hazard_test(..., time_transform="rank")`; HR CI `exp(coef ± z·se)`. Outputs
`cox_*_by_baseline_a1c.csv`, `models/cox_*_phreg.csv` (SAS-style for step8b), KM PNGs at 600 dpi.
Quirk: covariates for **both** outcomes come from the A1c analysis CSV; weight-cohort patients
absent from it get all-zero dummies (treated as reference).

### `step8b_cox_threshold_summary_table.py`

Joins the step8 KM summaries (`n`, `n_events`) with the PHREG CSVs into
`step8b_{weight,a1c}_threshold_cox_summary.csv` (`Outcome, Threshold, Baseline glycemic status,
e / n, HR (LCI-UCI), p`). Hard-codes `output/gap_{g}/…` paths, so the runners build a temporary
symlinked `output/` tree before calling it. `--out-markdown` names the output directory; no
markdown is written.

---

## F. Note-derived assessments (`add_unstructured/`)

Common conventions: no CLI; env `AU_DATADIR` (prepared data), `AU_OUTROOT` (output root),
`AU_STEP1` (weight analysis-ready CSV); `BASE_COV_FORMULA = "age + C(gender) + C(race) +
C(baseline_a1c_category) + C(baseline_bmi_final_category)"`; age median-imputed; GEE Gaussian /
Independence; `covid_era` and (PHQ-9) `antidepressant_baseline` as extra 0/1 terms.

### `prepare_assessment_data.py` (262 lines)

Reads the source CSV with `usecols` = ids, covariates and the six `{d}_*` columns per domain;
window `−180 ≤ days ≤ 365`; a row is an observation of domain *d* iff `{d}_present` and `{d}_value`
non-null; de-duplicate per patient-day (mean value); drop PHQ-2-only rows; derive `post`,
`time_months = days/30.44`, `time_post`, `has_both_periods`, `covid_era` (baseline date in
2020-03-11..2022-05-11), and for PHQ-9 `antidepressant_baseline` from a second chunked pass over
`medication_history` (26 keyword generics/brands, window −180/+90 d). Writes `{d}_prepared.csv`
and `domain_summary.json`. No unit conversion or plausibility bounds live here; those are upstream.

### `run_for_gap.py` (199 lines)

`--gap g [--name add_unstructured] [--force]`; env `CONF_INT_DIR`. Cutoff per patient = max
`days_from_baseline` in `output/step1_prepare_analysis_dataset/analysis_ready_gap{g}.csv`
(hard-coded main-pipeline path); keeps prepared rows with `days < 0 or days ≤ cutoff` for cohort
patients; then runs the six analysis scripts as subprocesses with `AU_*` set.

### `run_trajectory_plots.py` — Fig. 3f

Rows with `has_both_periods == 1`; `"{val} ~ bs(days_from_baseline, df=3, include_intercept=False)
+ BASE_COV_FORMULA [+ antidepressant_baseline] [+ covid_era]"`; grid −180..365 step 7; elevated
overlays for PHQ-9 ≥5/≥10, pain ≥4/≥7 (pre-period mean), alcohol ≥12. CI band suppressed where
<3 patients within ±21 d. Figures only.

### `run_baseline_anchor_analysis.py` — ED Table 4 (change from baseline)

`baseline_score` = mean of values within ±30 d; follow-up rows 30–365 d; `change_score = value −
baseline_score`; per-window (1–3, 3–6, 6–9, 9–12 mo) one-sample t and Wilcoxon; LMM
`smf.mixedlm("change_score ~ time_months + covariates", groups=patient_id, reml=True)`; GEE with a
synthetic day-0 anchor row (`change_score = 0`) and `bs(days, df=3)`; grid 0..365 step 7; report
`baseline_anchored_report.md`. The ED Table 4 point estimates themselves come from
`forest_point_estimates.py`.

### `forest_point_estimates.py` — ED Table 4

ITS-style estimate (full −180..365 GEE on `has_both_periods` rows) and CFB estimate (as above) at
days 91/183/274/365: `delta = ŷ(t) − ŷ(0)`, `se = sqrt(se_t² + se_0²)`. Subgroups: PHQ-9 All/≥5/≥10,
pain All/≥4/≥7, alcohol All/≥8.6 ("top 25 %"), waist and muscle All. `covid_era` is **not**
included here, unlike the other scripts. Output `point_estimates_3_6_9_12mo.csv`.

### `run_its_analysis.py`, `run_elevated_analysis.py`

Mixed-effects interrupted time series `"{val} ~ time_months + post + time_post + covariates"`
(random intercept per patient, REML; PHQ-9 and pain only) plus paired pre/post tests for all five
domains; the elevated variant restricts to pre-period mean ≥ threshold and drops the demographic
fixed effects. Supplementary results. The elevated script's comparison panel looks for a stale
path and labels and is never produced.

### `run_time_varying_covar.py`

Adds concurrent `pct_weight_change` (nearest weight observation via `merge_asof`, not the ±60 d
the docstring claims) as a time-varying covariate; base vs weight-adjusted GEE and CFB; alcohol
"top 25 %" is the run-time 75th percentile here. Outputs coefficient and point-estimate tables
and an xlsx.

### `select_spline_df_assessment.py`

QICu sweep over df 2–6 with no covariates per domain and gap; informational — the analysis
scripts hard-code df = 3.

### `freq/population_description.py` — Supp. Table 2; `freq/mention_frequency.py`

Module-level scripts (run on import). Population description of the PHQ-9 / pain / waist
populations against the gap-120 weight and A1c cohorts; mention/encounter frequency counts per
persistence definition. The production runner's `exec_module` path override is ineffective
(module top level reassigns `DATADIR`/`OUTDIR`), so they read and write `output/submitted_analysis/`.

---

## G. Runners

| Runner | Grid | Notes |
|---|---|---|
| `rerun_conf_int_clean_full.sh` | `GAPS_WITH_548` (all-data, note-derived), `GAPS_STRUCTURED` | Produced the manuscript numbers. `--resume`, `--keep-going`; `set -euo pipefail`; completion markers; `CI_CONF_LEVEL` inherited from the environment |
| `conf_int/run_all_gaps_all_data.sh` | `GAPS_PRIMARY` | Phase-3 subset into `output/conf_int/gap_g/all_data/`; does not pre-check configs |
| `conf_int/run_all_gaps_structured_only.sh` | `GAPS_STRUCTURED` | Phase-4 chain into `output/conf_int/gap_g/structured_only/`; **line 391 references an unset `$FORCE` under `set -u`, which aborts at step8b** |
| `add_unstructured/run_all_gaps.sh` | `GAPS_WITH_548` | Does not set `CONF_INT_DIR`, so writes to `output/submitted_analysis/` |

## H. Tests

`python3 tests/test_persistence_and_step0a.py` (needs pandas, numpy, scipy; no data). Pins:
integer patient IDs survive normalisation; the removed rule variants are absent; randomised
equivalence of `censor_day_for_patient` and `adherence_flags`; stop-day strictness; the abutting
`gap + 1` boundary; the day-0 marker; `z_critical`; fail-closed metformin helper; exclusion-flag
default. All pass on the committed code.

# Data schema — Part II (GLP-1 outcome analysis)

The pipeline consumes one wide CSV and produces per-gap analysis-ready CSVs, model configs,
coefficient tables, event tables and figures. None of the data is distributed; the schemas
below are reconstructed from how the code reads and writes each file.

## 1. Source CSV — `root_data/merged/step8g_with_unstructured_flags_with_assessments_weightcleaned.csv`

One row per **patient × event-day** (a measurement, a medication event, or a note-derived
assessment), long format, ~1.24 M rows. Columns are grouped by role. "Req." marks columns whose
absence stops the pipeline; the rest are optional with a documented fallback.

### 1.1 Identity and time

| Column | Type | Meaning | Req. |
|---|---|---|---|
| `patient_id` | string (UUID) | Patient key; normalised to stripped string by `persistence.normalize_patient_id` | yes |
| `date` | date | Event date of the row; becomes `event_date` | yes (or `medication_date`) |
| `medication_date` | date | Fallback event date when `date` is entirely null | |
| `baseline_glp1_date` | date | Index date = first qualifying injectable GLP-1 RA event (day 0) | yes |
| `days_from_baseline` | int | Present in source for assessment rows; recomputed by step1 as `(event_date − baseline_glp1_date).days` | |
| `event_type` | string | Row type; `"medication"` counts as GLP-1 evidence | |

### 1.2 Exposure and persistence

| Column | Type | Meaning | Req. |
|---|---|---|---|
| `baseline_glp1_brand_final` | string | Baseline brand; cohort filter matches `ozempic|wegovy|mounjaro|zepbound` as whole words | one of the three |
| `baseline_glp1_ingredient_final` | string | Baseline ingredient; matches `semaglutide|tirzepatide` | one of the three |
| `baseline_glp1` | string | Legacy agent label; fallback filter | one of the three |
| `glp1_event_for_adherance` | {0,1,2} | 1 = evidence of ongoing GLP-1 therapy on this row, 2 = explicit stop (never present in study data), 0 = none | yes |
| `glp1_days_from_baseline` | int | Day index of that GLP-1 evidence | yes |
| `is_glp1_event_med`, `is_glp1_event_clinical` | 0/1 | Structured prescription vs note-derived evidence; OR-ed into `glp1_evidence_therapy` | |
| `medication_history` | JSON string | List of `{"medication_name", "medication_date"}`; used for `metformin_with_glp1_baseline` (±180 d) and `antidepressant_baseline` (−180/+90 d) | |
| `glp1_user_group` | string | e.g. `"sema-only"`, `"tirz-only"`, `"switch…"`; used by step6d and the assessment scripts; **not propagated by step1** (step6d merges it from `root_data/step8f.csv`) | |
| `glp1_has_structured`, `glp1_has_unstructured` | 0/1 | Provenance flags; the structured-only prefilter zeroes `glp1_event_for_adherance` where `glp1_has_structured == 0` | |

### 1.3 Outcomes

| Column | Type | Meaning | Req. |
|---|---|---|---|
| `weight_in_pounds_final` | float, lb | Weight on this row | weight pipeline |
| `baseline_weight_final` | float, lb | Baseline weight (closest in −60/+14 d window) | **both** pipelines (A1c cohort also requires it) |
| `pct_weight_change` | float, % | `100·(weight − baseline)/baseline`; recomputed if absent | |
| `a1c_value` | float, % | HbA1c on this row | A1c pipeline (or `abs_a1c_change`) |
| `baseline_a1c_final` | float, % | Baseline HbA1c (lab or EAG-derived, per Methods) | A1c pipeline |
| `abs_a1c_change` | float, pp | `a1c_value − baseline_a1c_final`; recomputed if absent | |
| `A1C_SOURCE` | string | Source of the A1c value; nulled by the structured prefilter | |
| `weight_has_structured`, `weight_has_unstructured`, `a1c_has_structured`, `a1c_has_unstructured` | 0/1 | Provenance flags used by `structured_only/` | |

### 1.4 Covariates

| Column | Type / levels | Meaning |
|---|---|---|
| `baseline_a1c_category` | `Normal Glycemia` (<5.7 %), `Prediabetes` (5.7–<6.5), `Type 2 Diabetes` (6.5–<9.0), `Poorly Controlled Diabetes` (≥9.0), `Unknown` (dropped) | Effect modifier in every model; reference = Normal Glycemia |
| `baseline_bmi_final_category` | `Underweight`, `Normal`, `Overweight`, `Obese I`, `Obese II`, `Obese III`, `Unknown` | Covariate / stratifier |
| `age` | float, years (−1 when missing after step1 fill) | Source of `age_group` (decades), `age_group_20_39_vs_40_plus`, `age_group_20_49_vs_50_plus`, and step6cc's `<40/40+` |
| `gender` | `M`, `F`, `Unknown` | Covariate / stratifier |
| `race` | string, `Unknown` filled | Covariate; Supp. Table 2 collapses to Caucasian vs other |
| `weight_change_med` | 0/1 | Concomitant weight-affecting medication (RxNorm set, per Methods); covariate in every reported model |
| `baseline_height_final`, `height_in_inches_final`, `baseline_bmi_final`, `BMI_final` | float | Descriptives only |

### 1.5 Exclusion flags

| Column | Meaning |
|---|---|
| `pregnant_during_glp1` | 1 → row excluded (NaN → 0; absent → 0) |
| `bariatric_surgery_flag` | 1 → row excluded |

Oral-semaglutide exclusion and the age ≥ 18 restriction (Fig. 1) happen **upstream**; they are
not in this code.

### 1.6 Note-derived assessment columns (six per domain)

For each `d ∈ {phq9, pain_score, waist_circumference, alcohol, muscle_strength}`:

| Column | Meaning |
|---|---|
| `{d}_present` | 0/1 — an extraction of this domain exists on this patient-day |
| `{d}_value` | numeric value (PHQ-9 0–27; pain 0–10; waist inches; alcohol drinks/day; MRC 0–5), already unit-harmonised and plausibility-filtered upstream |
| `{d}_min`, `{d}_max` | min/max across extractions that day |
| `{d}_n_rows` | number of underlying NLP extractions |
| `{d}_source_names` | labels of the source instruments; PHQ rows matching `phq-2|phq2` without `phq-9|phq9` are dropped |

## 2. Analysis-ready CSVs — `analysis_ready_gap{g}.csv` (weight), `analysis_ready_a1c_gap{g}.csv`

Long format, post-baseline rows only (`0 ≤ days_from_baseline ≤ 730`), censored at each
patient's first non-persistent observation under gap `g`, one file per gap.

| Column | Weight | A1c | Notes |
|---|---|---|---|
| `patient_id`, `days_from_baseline` | ✓ | ✓ | |
| `baseline_carried_to_day0` | ✓ (1 for synthetic day-0 rows) | ✓ (always 0) | |
| `glp1_event_for_adherance`, `glp1_days_from_baseline` | ✓ | ✓ | persistence inputs |
| `days_since_last_glp1`, `days_since_last_glp1_evidence`, `days_to_{prev,next,nearest}_glp1`, `days_to_{prev,next,nearest}_glp1_evidence` | ✓ | ✓ | helper distances, unused by models |
| `adherence_{g}` | ✓ | ✓ | 1 for every retained row by construction |
| `pct_weight_change` | ✓ | — | outcome |
| `baseline_a1c_final`, `a1c_value`, `abs_a1c_change` | — | ✓ | outcome |
| `baseline_a1c_category`, `baseline_bmi_final_category`, `age_group`, `age_group_20_39_vs_40_plus`, `age_group_20_49_vs_50_plus`, `gender`, `race`, `weight_change_med`, `metformin_with_glp1_baseline` | ✓ | ✓ | covariates |
| `baseline_glp1`, `baseline_glp1_brand_final`, `baseline_glp1_ingredient_final`, `baseline_glp1_date`, `glp1_evidence_therapy` | ✓ | ✓ | exposure |
| `baseline_weight_final`, `weight_in_pounds_final` | ✓ | ✓ | |
| `age`, `baseline_height_final`, `height_in_inches_final`, `baseline_bmi_final`, `BMI_final` | ✓ | — | descriptives |
| `a1c_has_*`, `weight_has_*`, `glp1_has_*` | ✓ | ✓ | provenance |

Not emitted: `event_date`, `glp1_user_group`, `medication_history`.

## 3. Intermediate and output artefacts

| File | Producer | Schema |
|---|---|---|
| `model_config.json` / `model_config_a1c.json` | step2 | `{"best_df": int, "best_qicu": float}`; `df_selection*.csv` has `df, QIC, QICu[, error]` |
| `coefficients.csv` | step3 | `term, estimate, std_error, ci_lower, ci_upper, p_value` |
| `gee_combined_coefficients.csv` | step5 | `term, Coef, CI Lower, CI Upper, p-value` |
| `forest_predictions.csv` | step5, step6c | `day, a1c_group, pred, ci_low, ci_high, se, n_obs_window, n_unique_window` |
| `forest_contrasts_vs_ref.csv` | step5, step6c | `day, ref, a1c_group, diff, ci_low, ci_high, se` |
| `tables/predicted_means_by_day_*.csv` | step4 pred | `a1c_group, day, pred_mean_uncentered, ci_low_uncentered, ci_high_uncentered, pred_mean_centered, ci_low_centered, ci_high_centered, supported` → ED Table 3 (days 84, 183, 280, 365) |
| `summary_*_by_baseline_a1c_{90day,months}.csv` | step4 observed | binned observed means/SDs |
| `contrasts_{strat_var}.csv` | step6b | `strat_var, ref, subgroup, day, diff, ci_low, ci_high, se, n_people_ref, n_obs_ref, n_people_sub, n_obs_sub` |
| `stratified_summary_counts.csv` | step6 | `strat_var, subgroup, n_people, n_obs, scale` |
| `pooled_coefficients.csv`, `{cat}/{group}/predicted_trajectory.csv`, `forests/forest_day_{d}.csv` | step6d | sema-only vs tirz-only |
| `step8_weight_time_to_threshold_events.csv` | step8 weight | `patient_id, baseline_a1c_category, threshold_pct ∈ {−5,−10,−15}, time_days (≤540), event` |
| `step8_a1c_time_to_threshold_events.csv` | step8 A1c | same with `threshold_abs ∈ {0.5,1.0,1.5,2.0}` |
| `step8_*_time_to_threshold_summary.csv` | step8 | `baseline_a1c_category, threshold, n, n_events, median_time_days, median_time_days_lower_ci, median_time_days_upper_ci` |
| `cox_*_by_baseline_a1c.csv` | step8 KM+Cox | lifelines `cph.summary` + `HR, HR_lower_95, HR_upper_95, ph_global_p_for_a1c_block` |
| `models/cox_*_phreg.csv` | step8 KM+Cox | SAS-PHREG-like: `section, metric, value, Parameter, Level, HazardRatio, HRLowerCL, HRUpperCL, WaldChiSq, PrChiSq, covariate` |
| `step8b_{weight,a1c}_threshold_cox_summary.csv` | step8b | `Outcome, Threshold, Baseline glycemic status, e / n, HR (LCI-UCI), p` → ED Table 2 Cox block |
| `step0_baseline_table_gap{g}.csv` | step0 | `variable, level, Total, Normal Glycemia, Prediabetes, Type 2 Diabetes, Poorly Controlled Diabetes` → Table 3 |
| `samplesize_by_month.csv` | step0a | `month_number, days_from_baseline, gap_30 … gap_730` → Supp. Table 1 |
| `{d}_prepared.csv` | prepare_assessment_data / run_for_gap | `patient_id, days_from_baseline, {d}_value, {d}_min, {d}_max, {d}_n_rows, {d}_source_names, age, gender, race, baseline_a1c_category, baseline_bmi_final_category, glp1_user_group, baseline_glp1_date, date, post, time_months, time_post, has_both_periods, covid_era [, antidepressant_baseline]` |
| `point_estimates_3_6_9_12mo.csv` | forest_point_estimates | `analysis ∈ {ITS, CFB}, subgroup, timepoint, day, n, estimate, se, ci_lo, ci_hi, domain` → ED Table 4 |
| `tvc_coefs_{d}.csv`, `point_estimates_3_6_9_12mo_tvc.csv` | run_time_varying_covar | weight-adjusted sensitivity |
| `descriptive_summary.csv`, `model_coef_{d}.csv` | run_its_analysis | paired pre/post and ITS coefficients |
| `add_unstructured_population_description.{csv,xlsx}` | freq/population_description | Supp. Table 2 |

## 4. Output directory layout (as produced by `rerun_conf_int_clean_full.sh`)

```
output/conf_int_clean/
├─ rerun_master.log, logs/
├─ step1_weight/analysis_ready_gap{30..730}.csv          step1_a1c/analysis_ready_a1c_gap{…}.csv
├─ prefiltered_structured_only.csv
├─ structured_only_step1_weight/…                         structured_only_step1_a1c/…
├─ 1_no_adherence/  (and 1_no_adherence_full/, identical)
│  ├─ data/{phq9,pain_score,waist_circumference,alcohol,muscle_strength}_prepared.csv, domain_summary.json
│  ├─ ITS/{figures,tables,models}/, GEE/figures/, CFB/{figures,tables}/, Forest/{figures,tables}/, time_var_covar/{figures,tables}/
├─ gap_{g}/
│  ├─ add_unstructured/{data,ITS,GEE,CFB,Forest,time_var_covar}/
│  ├─ all_data/{step4_predictive_plots, step4_predictive_plots_a1c, step4_observed_summary_plots,
│  │            step6_stratified_by_covariates_{weight,a1c}, step6cc_3way_by_age_sex_{weight,a1c},
│  │            step8_survival_time_to_{weight_loss,a1c_drop}, step8_survival_plots_and_cox}/
│  └─ structured_only/{step0_…, step2_…, step3_…, step4_…, step5_…, step6_…, step6b_…, step6c_…, step6cc_…, step6d_…, step8_…, step8b}/
└─ spline_selection/
```

Phase 3 additionally **reads** `output/gap_{g}/step2_select_spline_df/model_config.json` and
`…_a1c/model_config_a1c.json` from the main pipeline, and `step6d` / `step0` fall back to
`root_data/step8f.csv` (an earlier wide file with `glp1_user_group` and standardised brand
columns) when those columns are missing.

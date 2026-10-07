# Paper coverage and modularity — Part II (GLP-1 outcome analysis)

## 1. Methods-to-code map

Every analytic decision described in the Methods, with where it lives (or does not).

| Methods statement | Implemented in | Status |
|---|---|---|
| Data source, de-identification, upstream normalisation to standard terminologies | — | upstream platform |
| Agent-assisted note extraction, delta tracking, verification, ontology mapping, concept discovery, KG | — | upstream platform (Part I validates its outputs) |
| Active adults, GLP-1 exposure via prescriptions + note mentions, 1 M-visit sample | — | upstream (Fig. 1 top boxes) |
| Age ≥18, missing age/sex excluded; oral semaglutide excluded | — | upstream; **not in this code** |
| Pregnancy during window, bariatric history excluded | `step1_*` lines 343–351 | here (flags consumed, not derived) |
| Injectable semaglutide / tirzepatide cohort (brand + ingredient) | `step1_prepare_analysis_dataset.py` 39–50, 353–372 | here |
| Baseline window −60/+14 d, nearest value, ties on/before | — | upstream (`baseline_*_final` columns arrive computed) |
| Day-0 observation for every patient | `_ensure_baseline_rows` | here (weight only; see README §1) |
| Weight plausibility (45–200 % of median, 40 % spikes), BMI recalculation, HbA1c spikes, CDC outlier screen, EAG conversion | — | upstream; the input is already `…_weightcleaned.csv` |
| HbA1c categories 5.7 / 6.5 / 9.0 | — | upstream (`baseline_a1c_category` arrives labelled); order and reference enforced by `model_spec` |
| BMI categories | — | upstream |
| Concomitant weight-affecting medication covariate (RxNorm set) | — | upstream (`weight_change_med` arrives as a flag) |
| Metformin near baseline | `_has_metformin_near_baseline` | here (candidate covariate; excluded from reported models) |
| Persistence: backward gap rule, eight thresholds, g = 120 primary | `persistence.py`, `gap_grids.sh` | here |
| Follow-up ≤18 months | 730 / 548 / 540 caps (README §5) | here |
| GEE, Gaussian, independence working correlation, robust SE, cluster = patient, B-spline 3 df chosen by QICu | step2, step3, step4, step5, step6* | here |
| Single pooled model with category as covariate (ED Table 2 GEE differences) | `step3_fit_gee_baseline*` (main effects) and `step5` (interaction contrasts) | here |
| Separate models per category / sex / age for Fig. 3a–d, ED Table 3 | `conf_int/gap_120/step4_predictive_plots*`, `step6_stratified_by_covariates_*` | here |
| KM curves, Cox PH with PH test; thresholds 5/10/15 % and 0.5/1.0/1.5 pt | `step8_*`, `conf_int/gap_120/step8_survival_plots_and_cox.py` | here (2.0 pt also computed; omitted from the paper as non-estimable) |
| Structured-only sensitivity | `structured_only/step0_prefilter_raw.py` + full chain in runner phase 4 | here |
| Note-derived assessments: ±30 d baseline, 30–365 d follow-up, GEE 3 df, LMM ITS, COVID indicator, antidepressant adjustment, time-varying weight | `add_unstructured/*` | here |
| Unit harmonisation and plausibility bounds for assessments (pain ≤10, PHQ-9 ≤27, MRC ≤5, 1st–99th percentile) | — | upstream; **not in this code** |
| Elevated-baseline subgroups (PHQ-9 ≥5/≥10, pain ≥4/≥7, alcohol top quartile) | `forest_point_estimates.py`, `run_trajectory_plots.py`, `run_time_varying_covar.py` | here, with **three different alcohol definitions** (8.6, 12, run-time p75) |
| Two-sided tests, CIs, no multiplicity adjustment | `analysis_config` + each script | here |

## 2. Display-item map

| Item | Script(s) | Output file | Reproducible without data? |
|---|---|---|---|
| Fig. 1 | counts from step1 logs, `step0a` | — | no (figure script not released) |
| Table 3 | `step0_analysis_population_table.py` | `step0_baseline_table_gap120.csv` | no |
| Fig. 3a, b | `conf_int/gap_120/step4_predictive_plots{,_a1c}.py` | `grouped/grouped_trajectory_centered_outcome_{weight,a1c}.png` | no |
| Fig. 3c, d | `conf_int/gap_120/step6_stratified_by_covariates_weight.py` | `by_gender/plots/grouped/…`, `by_age_group_20_39_vs_40_plus/…` | no |
| Fig. 3e | `conf_int/gap_120/step8_survival_plots_and_cox.py` | `km_*_by_baseline_a1c.png` | no |
| Fig. 3f | `add_unstructured/run_trajectory_plots.py` (and CFB variants) | `GEE/figures/trajectory_gee_{phq9,pain_score,waist_circumference}.png` | no |
| ED Table 2 (GEE) | `step5_forest_contrasts_{weight,a1c}.py` | `forest_contrasts_vs_ref.csv` | no |
| ED Table 2 (Cox) | `step8_survival_plots_and_cox.py` → `step8b` | `step8b_*_threshold_cox_summary.csv` | no |
| ED Table 3 | `step4_predictive_plots*` | `tables/predicted_means_by_day_*.csv` at days 84/183/280/365 | no |
| ED Table 4 | `forest_point_estimates.py` (+ `run_time_varying_covar.py`) | `Forest/tables/point_estimates_3_6_9_12mo.csv` | no |
| Supp. Table 1 | `step0a_samplesize_analysis.py` | `samplesize_by_month.csv` | no |
| Supp. Table 2 | `freq/population_description.py` | `add_unstructured_population_description.csv` | no |
| Supplementary figure sets (stratified forests, 3-way, sema vs tirz, ITS) | step6b/6c/6cc/6d, `run_its_analysis.py`, `run_elevated_analysis.py` | various | no |

Not included at all: the xlsx/PDF table formatters, multi-panel figure assemblers, the
patient-flow figure script (README scope note), and the `build_*_pdf.py` builders the runner
looks for in phase 5.

## 3. What is not available

1. **The data.** The single input CSV is licensed; nothing in the repo can be executed end to end
   without it. The regression tests are the only runnable artefact.
2. **The cleaning layer.** Everything the Methods describe under "Data preprocessing and
   plausibility criteria", "Weight", "BMI", "HbA1c", "Anthropometric cleaning and CDC outlier
   detection", "Baseline definition" and "Derived variables" (categories, EAG conversion) is
   upstream. The CSV arrives cleaned and labelled.
3. **Cohort entry.** Age, oral semaglutide, "active adult", and the 1 M-visit sampling are
   upstream.
4. **`root_data/step8f.csv`**, an earlier wide file that step0 and step6d fall back to for
   `glp1_user_group` and brand columns.
5. **`step7`** (adherence counts; only annotates step5 titles).
6. **Deliverable formatting** scripts.

## 4. Modularity assessment

| Layer | Modular? | Evidence |
|---|---|---|
| Persistence rule | **Yes.** One module, pure functions over three columns, tested, two documented readings. | `persistence.py` |
| Model specification guards | **Yes.** `load_spline_df`, `enforce_a1c_order`, `filter_estimable`, `z_critical` are small, importable, and used by most modelling scripts. | `model_spec.py`, `covariates.py`, `analysis_config.py` |
| Cohort construction | Partly. `load_and_prepare` is one function but mixes exclusion, cohort filter, outcome derivation, covariate fills, date handling and windowing; the A1c module reuses it by importing private helpers from the weight module. Agent names are a single constant. | `step1_*` |
| Trajectory modelling | **No.** The GEE recipe is re-implemented in ~14 scripts (`_pred_ci`, modal covariates, grid, centering, plotting), each with its own `A1C_ORDER` literal, y-limits, legend text and outdir routing. step3 keeps a private copy of the covariate filter. | `step4`, `step5`, `step6*` |
| Time-to-event | Partly. Event construction is isolated in `step8_*`; KM/Cox in one script; step8b depends on a hard-coded path tree. | `step8*` |
| Note-derived assessments | **No.** Seven standalone scripts with no CLI, env-variable configuration, duplicated covariate formulas, and three definitions of "top 25 % alcohol". | `add_unstructured/` |
| Orchestration | Partly. One canonical runner with markers and fail-fast; three others overlap; phase 3 depends on a main-pipeline directory outside `conf_int_clean/`. | `*.sh` |
| Configuration | **No central config.** Column names, category labels, thresholds, horizons and paths are literals in each file; `gap_grids.sh` and `CI_CONF_LEVEL` are the only shared settings. | |

Net: the review response (v2) made the *rules* modular (persistence, df, reference level,
estimability, confidence) while leaving the *scripts* as a flat collection of near-duplicates.
Reusing this for a new question is a find-and-replace exercise across ~40 files rather than a
configuration change; see `REUSE_GUIDE.md` for the exact edit points.

## 5. Defects worth knowing before trusting a re-run

Verified against the committed code (none affects the primary reported estimates; see
`CHANGE_LIST.md`/`VERIFICATION.md` references in the README for the authors' own audit).

- `conf_int/run_all_gaps_structured_only.sh:391` uses `$FORCE`, never set, under `set -u` — the
  script aborts at step8b of the first gap. Use `rerun_conf_int_clean_full.sh` instead.
- `rerun_conf_int_clean_full.sh` phase 5 overrides `DATADIR`/`OUTDIR` on the freq modules before
  `exec_module`, which re-executes the module top level and resets them, so `mention_frequency.py`
  and `population_description.py` read/write `output/submitted_analysis/`, not `conf_int_clean/`.
- `step6cc_*`: `getattr(args, "spline-df", 3)` ignores `--spline-df` (df is always 3); the A1c
  script has a duplicated nested loop that refits every stratum once per outer stratum.
- `step8_survival_plots_and_cox.py`: Cox covariates for the weight outcome are taken from the A1c
  cohort file; weight-only patients get reference-level dummies. `_format_covariate_name_to_phreg`
  splits multi-underscore dummy names at the first underscore.
- `run_for_gap.py` and `freq/mention_frequency.py` always take cutoffs from
  `output/step1_prepare_analysis_dataset/` regardless of `CONF_INT_DIR`.
- `run_elevated_analysis.py`'s comparison panel reads a stale path and labels and is never produced.
- `run_time_varying_covar.py` matches the nearest weight in time, not within ±60 d as its docstring says.
- Support masks in `step4_predictive_plots*` and `step6_*` are computed but not applied to the plots.

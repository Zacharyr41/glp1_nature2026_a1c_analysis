# Reuse guide — same pipeline, different clinical question, or MIMIC

The pipeline is a generic "exposure-anchored longitudinal trajectory + time-to-threshold" design:
index date, baseline window, persistence-based censoring, GEE splines by a baseline effect
modifier, KM/Cox to thresholds, plus a note-derived-assessment module. Everything GLP-1-specific
is in column names, label strings and constants. This guide lists, in order, what to produce,
what to edit, and what you have to write from scratch.

## 1. The contract: produce one wide CSV

The entire pipeline starts from one long-format file (`DATA_SCHEMA.md` §1). Reproducing the paper
on new data means producing that file. The minimum viable column set:

```
patient_id, date, baseline_glp1_date,
baseline_glp1_brand_final | baseline_glp1_ingredient_final | baseline_glp1,   # cohort filter
glp1_event_for_adherance, glp1_days_from_baseline,                            # persistence
weight_in_pounds_final, baseline_weight_final,                                # outcome 1
a1c_value, baseline_a1c_final,                                                # outcome 2
baseline_a1c_category, baseline_bmi_final_category, age, gender, race, weight_change_med,
pregnant_during_glp1, bariatric_surgery_flag,
[*_has_structured, *_has_unstructured]                                        # only for the sensitivity cohort
[{d}_present, {d}_value, {d}_min, {d}_max, {d}_n_rows, {d}_source_names]      # only for note-derived assessments
[medication_history]                                                          # only for metformin / antidepressant flags
```

Everything the Methods describe as cleaning (plausibility ranges, spike removal, BMI
recalculation, EAG → HbA1c, baseline-window selection, category assignment) must be done **before**
this file exists. None of it is in the repo; §4 lists what to write.

## 2. Adapting to a different clinical question

Think of the design as five slots. Each slot maps to a small set of literals.

### 2.1 Exposure / index event

| What | Where |
|---|---|
| Agent names used by the cohort filter | `step1_prepare_analysis_dataset.py` `GLP1_BRAND_NAMES`, `GLP1_INGREDIENT_NAMES` (lines 39–40); the A1c module imports the patterns |
| Index date column | `baseline_glp1_date` throughout (step1, step0, step0a, `prepare_assessment_data.py`) |
| Evidence-of-ongoing-therapy rows | `glp1_event_for_adherance` ∈ {0,1,2} and `glp1_days_from_baseline`; keep the names and the rule, just populate them from your exposure events (fills, orders, note mentions) |
| Persistence thresholds | `gap_grids.sh`; `GAPS_DEFAULT` in both step1 files |
| Agent sub-groups (step6d) | `GLP1_TWO_GROUPS`, the `sema`/`tirz`/`switch` substring mapping (lines 85, 466–476) |

For a non-drug exposure (surgery, diagnosis, device) set `glp1_event_for_adherance = 1` on the
index row only and use a single very large gap; the persistence machinery then reduces to
"follow-up until last observation".

### 2.2 Outcomes

| What | Where |
|---|---|
| Continuous outcome 1 (`pct_weight_change`) | step1 weight: `_pct_weight_change`; formula strings in step2/3/4/5/6*/6cc/6d weight scripts; y-limits (−10, 4) in step4/step6 |
| Continuous outcome 2 (`abs_a1c_change`) | step1 A1c lines 63–67; `_a1c` scripts; y-limits (−3, 1) |
| Event thresholds | `WEIGHT_THRESHOLDS`, `A1C_THRESHOLDS` in `step8_survival_time_to_*.py`; file-name tags (`5pct_loss`, `0p5_reduction`) in `step8_survival_plots_and_cox.py` and `step8b` |
| Sign convention | "improvement is negative" is assumed in step8 (`<= thr`) and in step0's achievement flags |
| Follow-up caps | 730 (`MAX_DAYS`), 548 (`--max-days`, `--truncate-days`), 540 (`MAX_FOLLOWUP_DAYS`) |

The simplest route is to keep the column names and write your outcome into `pct_weight_change`
(percent change from baseline) and `abs_a1c_change` (absolute change); then only labels and
thresholds need editing.

### 2.3 Effect modifier (the "baseline glycaemic category")

| What | Where |
|---|---|
| Level names and reference | `model_spec.A1C_ORDER` **and** local copies in step0, step4 (with `Unknown`), step5, step6b A1c, step6, step8 Cox, `population_description.py` |
| Legend cut-points text | `step4_predictive_plots*` lines 246–252 |
| Column name `baseline_a1c_category` | ~230 occurrences; keep the name, change the levels |

### 2.4 Covariates

Candidate lists are literals in every modelling script: `age_group, gender,
baseline_a1c_category, baseline_bmi_final_category, race, weight_change_med
[, metformin_with_glp1_baseline]`. `filter_estimable` will drop anything absent, so the quickest
adaptation is to supply columns with those names (e.g. put your comorbidity flag in
`weight_change_med`). Age grouping functions are in step1 `_derive_age_groups`.

### 2.5 Note-derived assessments

`DOMAINS` in `prepare_assessment_data.py` (line 40) and the `DOMAINS` dicts in each
`add_unstructured/run_*.py` (labels, scale ranges, elevated thresholds, spline df, which get an
ITS model). The `{d}_present/_value/...` naming convention is the only coupling to the CSV.
Change `ANTIDEP_KEYWORDS`, the COVID dates, and the −180/+365 window as needed.

### 2.6 Practical order of work

1. Build the CSV (§1). Run `step1_*` for one gap and inspect `analysis_ready_gap*.csv`.
2. Run step2 → step3 → `conf_int/gap_120/step4_predictive_plots*` for that gap; this is Fig. 3a–b.
3. Run step5 for ED Table 2, step8 + `step8_survival_plots_and_cox.py` for Fig. 3e / ED Table 2.
4. Only then run the full runner across gaps.
5. Set `CI_CONF_LEVEL` once; never edit `1.96`.

## 3. Adapting to MIMIC

MIMIC-IV (with `mimic-iv-note`) can feed this pipeline, but the paper's design is an **outpatient,
multi-year, note-dense** cohort and MIMIC is an **inpatient / ED** database. Expect fewer repeat
measurements per patient, shorter observable follow-up, and a different meaning of "persistence".

### 3.1 Table mapping for the wide CSV

| CSV column | MIMIC-IV source | How |
|---|---|---|
| `patient_id` | `hosp.patients.subject_id` | cast to string |
| `date` | `hosp.labevents.charttime`, `icu.chartevents.charttime`, `hosp.omr.chartdate`, `hosp.prescriptions.starttime` | one row per measurement / order |
| `baseline_glp1_date` | first `hosp.prescriptions` (or `hosp.emar`) row whose `drug` matches your exposure; or `hosp.diagnoses_icd` for a diagnosis index | per patient |
| `baseline_glp1_brand_final` / `_ingredient_final` | `prescriptions.drug`, `gsn`/`ndc` → RxNorm | string-match or RxNav lookup |
| `glp1_event_for_adherance`, `glp1_days_from_baseline` | each `prescriptions`/`emar` row for the exposure → 1; a documented discontinuation → 2 | days = `starttime − baseline` |
| `weight_in_pounds_final` | `icu.chartevents` itemids 226512 (admission weight, kg), 224639 (daily weight, kg); `hosp.omr` `result_name = 'Weight (Lbs)'` | kg × 2.20462; apply the Methods' plausibility rules yourself |
| `baseline_weight_final` | nearest weight in −60/+14 d of index, ties on/before | per patient |
| `a1c_value` | `hosp.labevents` itemid 50852 (`%`) | |
| `baseline_a1c_final` | nearest in −60/+14 d | |
| `baseline_a1c_category` | from `baseline_a1c_final` at 5.7 / 6.5 / 9.0 | |
| `baseline_bmi_final_category` | `hosp.omr` `BMI (kg/m2)` or weight/height² (height itemid 226730, cm) | CDC classes |
| `age` | `anchor_age + (year(index) − anchor_year)`; ages > 89 are 91 | |
| `gender` | `patients.gender` | M/F |
| `race` | `admissions.race` | collapse to your scheme |
| `weight_change_med` | `prescriptions.drug` against your RxNorm list | per patient-window |
| `pregnant_during_glp1` | ICD-10 O00–O9A / ICD-9 630–679 in `diagnoses_icd` within the window | |
| `bariatric_surgery_flag` | ICD-10-PCS 0DB64Z3, 0D164ZA, … or ICD-9-CM 43.89, 44.31, 44.38, 44.39 in `procedures_icd` | |
| `medication_history` | JSON of `prescriptions` rows per patient | needed only for the metformin / antidepressant flags |
| `*_has_structured` / `*_has_unstructured` | 1/0 by source table (labs/chartevents vs note extraction) | |
| `{d}_present`, `{d}_value`, … | **from your own note extraction** over `discharge.csv` / `radiology.csv` (MIMIC-III: `NOTEEVENTS`) | see §4 |
| `glp1_user_group` | derived from the exposure drugs | only for step6d |

Useful public helpers: the MIT-LCP `mimic-code` repository has concept SQL for weight,
height, BMI and comorbidities; use it rather than re-deriving itemids.

### 3.2 Design adjustments you will need to make

- **Index event.** GLP-1 initiation in MIMIC is rare and inpatient; a more natural MIMIC
  question is e.g. "first insulin order", "first loop-diuretic order", "sepsis onset", with
  outcomes like creatinine, lactate, weight, or a note-derived score. The pipeline does not care.
- **Persistence.** Inpatient `emar` administrations are daily, so a 120-day gap is meaningless;
  use gaps of 1–3 days within an admission, or treat the admission as the follow-up window.
- **Time scale.** `days_from_baseline` can be hours ÷ 24 as floats; the splines and KM do not
  require integers. Replace the 730/548/540 caps and the 14-day grid step accordingly.
- **Day-0 anchoring** (`_ensure_baseline_rows`) assumes a baseline value exists within the window;
  with inpatient data most patients will have an observed day-0 row anyway.
- **Repeat measures.** GEE with independence working correlation and robust SE is fine for ICU
  density, but the spline df selection (QICu, grid 3–6) may prefer higher df; let step2 choose.
- **De-identification.** MIMIC dates are shifted per patient into 2100–2200; `covid_era` is
  therefore meaningless and should be removed from the assessment formulas.
- **Race / sex.** `admissions.race` has ~30 levels; collapse before step1 or `filter_estimable`
  will keep sparse levels and Cox dummies will explode.

### 3.3 The note-derived module on MIMIC

`prepare_assessment_data.py` expects the per-domain columns already extracted. For MIMIC you
write the extractor (Part I's `docs/REUSE_GUIDE.md` §3 covers it). Discharge summaries do carry
pain scores, PHQ-type screens occasionally, weight/waist rarely, and alcohol use (drinks/day) in
social history frequently; MRC muscle strength almost never. Choose domains accordingly.

## 4. What you must write yourself (not in either repository)

Ordered by how early you need it.

1. **Cohort and index-date builder** — identify exposure, first qualifying event, active-patient
   criteria, age/sex completeness, oral-formulation exclusion.
2. **Baseline-window selector** — nearest value in −60/+14 d with tie rule, for each baseline variable.
3. **Plausibility cleaning** — the full Methods list: height 48–90 in, weight 60–700 lb, BMI 10–80,
   45–200 % of patient median, 40 % spike rule, 60 % of median baseline rule, BMI recalculation,
   HbA1c 3–20 % with the ±365-day spike rule, EAG → HbA1c `(EAG + 46.7)/28.7`, CDC percentile screen.
4. **Category assignment** — HbA1c and BMI classes, age groups (step1 does age groups; the rest is upstream).
5. **Concomitant-medication flag** — RxNorm set for weight-affecting drugs → `weight_change_med`.
6. **Exposure evidence rows** — `glp1_event_for_adherance` from fills, orders and note mentions;
   explicit stops as 2.
7. **Provenance flags** — `*_has_structured` / `*_has_unstructured` per outcome per row.
8. **Note extraction + verification + ontology mapping** for anything note-derived
   (Part I guide), producing `{d}_present/_value/_min/_max/_n_rows/_source_names` and the
   `is_glp1_event_clinical` evidence rows.
9. **Assessment harmonisation** — units (waist to inches), bounds (pain ≤10, PHQ-9 ≤27, MRC ≤5),
   1st–99th percentile trimming.
10. **`root_data/step8f.csv` equivalent** — only if you use step6d or step0's brand fallback; otherwise add `glp1_user_group` to the main CSV.
11. **Deliverable formatting** — xlsx/PDF tables, multi-panel figure assembly, flow diagram.
12. **A configuration layer** (recommended) — the repo has none; a single YAML mapping of
    `{outcome columns, category levels, thresholds, caps, covariates}` would replace the
    per-script literals listed in §2.

## 5. Minimal smoke test without real data

Build a synthetic CSV with ~200 patients, 10–30 rows each, the §1 columns, and run:

```bash
python3 code/step1_prepare_analysis_dataset.py --input-csv synthetic.csv --outdir out/step1 --adherence-gaps 120
python3 code/step1_prepare_analysis_dataset_a1c.py --input-csv synthetic.csv --outdir out/step1_a1c --adherence-gaps 120
python3 code/step2_select_spline_df.py --input-csv out/step1/analysis_ready_gap120.csv --outdir out/step2 --adherence-gap-days 120
python3 code/step3_fit_gee_baseline.py --input-csv out/step1/analysis_ready_gap120.csv --config-json output/gap_120/step2_select_spline_df/model_config.json --adherence-gap-days 120
python3 code/conf_int/gap_120/step4_predictive_plots.py --input-csv out/step1/analysis_ready_gap120.csv --config-json output/gap_120/step2_select_spline_df/model_config.json --outdir output/gap_120/step4 --adherence-gap-days 120
python3 code/step8_survival_time_to_weight_loss.py --input-csv out/step1/analysis_ready_gap120.csv --adherence-gap-days 120
```

Note the gap routing: step2 writes to `output/gap_120/step2_select_spline_df/` regardless of
`--outdir` when `--adherence-gap-days` is given and the outdir lacks `/gap_`.

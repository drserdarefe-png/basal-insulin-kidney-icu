# Basal insulin prescribing across kidney function in critically ill adults

Analysis code for the manuscript:

> Efe S, Özgür A. *Basal Insulin Prescribing across Kidney Function in Critically Ill Adults: Two Observational Cohorts.* Submitted to the Clinical Journal of the American Society of Nephrology (CJASN), 2026.

The study describes basal glargine prescribing across estimated glomerular filtration rate (eGFR) in MIMIC-IV v3.1 and eICU-CRD v2.0, compares observed doses with three published kidney function-based dosing schedules, and tests whether observational estimators (augmented inverse-probability weighting, within-patient first differencing) can recover the glucose-lowering effect of insulin. A simulation study with a known effect benchmarks the estimators.

## Data access

**No patient-level data are included in this repository.** MIMIC-IV and eICU-CRD are available from PhysioNet to credentialed users who complete the required training and sign the data use agreement:

- MIMIC-IV v3.1: https://doi.org/10.13026/kpb9-mt58
- eICU-CRD v2.0: https://doi.org/10.13026/C2WM1R

The SQL scripts query the PhysioNet BigQuery datasets (`physionet-data.mimiciv_3_1_*`, `physionet-data.eicu_crd`). Run them in your own Google Cloud project after linking your PhysioNet account.

## Repository structure

| Path | Content |
|---|---|
| `EXCLUSION_REGISTRY.md` | Versioned exclusion-criteria registry (v1.3). Every exclusion line in the SQL refers to an ID here (`-- EXC-xx`). |
| `sql/` | MIMIC-IV cohort, exposure, covariate and dose-event extraction (00 to 07). |
| `sql/eicu/` | eICU-CRD glargine order extraction, dose parsing and center-level summaries. |
| `ref/` | Reference tables: insulin pharmacokinetics, enteral formula and dextrose carbohydrate content, corticosteroid potency and routes. Item counts are aggregate. |
| `python/` | Analyses (see below). |
| `rerun_27.11.2026.sh` | Re-runs the primary analyses on the final data set and checks the expected counts (7,954 glargine events, 1,117 patients). |

Code comments are partly in Turkish.

## Analysis map

| Manuscript analysis | Script |
|---|---|
| Dose determinants, SHAP, cross-validated R² | `python/07_gps_diagnostic.py` |
| Doubly robust (AIPW) dose-response by eGFR band | `python/09_aipw_basal.py`, `python/09_aipw_basal_v1_1.py` |
| Within-patient first-difference models and sensitivity | `python/10_within_patient.py`, `python/10b_within_patient_sens.py` |
| Supplemental tables, hypoglycemia by eGFR and dose | `python/12_paper2_supp.py`, `python/13_steroid_arm.py` to `python/19_ci_sepsis.py` |
| Simulation (three data-generating worlds) | `python/20_three_worlds.py` |
| eICU-CRD prescribing, policy comparison, Mundlak, ICC | `python/21_eicu_analysis.py` to `python/24_eicu_closing.py` |
| Observed versus policy-implied dose figure | `python/25_fig4_v9.py` to `python/27_fig4_f10.py` |
| Exploratory scripts not reported as primary results | `python/06_positivity_check.py`, `python/07_isf_two_stage.py`, `python/08_dose_optimizer.py` |

## Requirements

Python 3.11 or later. See `requirements.txt`. The SQL outputs are exported to Parquet and read locally with DuckDB (`--parquet 'data/basal_events/*.parquet'`).

## Citation

If you use this code, please cite the article above together with the MIMIC-IV, eICU-CRD and PhysioNet references listed in it.

## License

MIT License. See `LICENSE`.

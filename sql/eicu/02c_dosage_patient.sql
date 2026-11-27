-- eICU IS3c: HASTA-duzeyi parse edilebilirlik (>=1 sayisal glarjin dozu)
WITH g AS (
  SELECT patientunitstayid, dosage
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
)
SELECT
  COUNT(DISTINCT patientunitstayid) AS n_glargine_stays,
  COUNT(DISTINCT IF(REGEXP_CONTAINS(dosage, r"[0-9]"), patientunitstayid, NULL)) AS n_with_numeric_dose,
  ROUND(100*COUNT(DISTINCT IF(REGEXP_CONTAINS(dosage, r"[0-9]"), patientunitstayid, NULL))
        / COUNT(DISTINCT patientunitstayid), 1) AS pct_parseable_patient
FROM g

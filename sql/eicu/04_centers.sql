-- eICU IS5: merkez varyasyonu (glarjin) -- EN DEGERLI
WITH gd AS (
  SELECT m.patientunitstayid, p.hospitalid,
    SAFE_CAST(REGEXP_EXTRACT(m.dosage, r"[0-9]+(?:\.[0-9]+)?") AS FLOAT64) AS dose_u
  FROM `physionet-data.eicu_crd.medication` m
  JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid)
  WHERE REGEXP_CONTAINS(LOWER(m.drugname), r"glargine|lantus|toujeo|basaglar|semglee")
),
per_hosp AS (
  SELECT hospitalid,
    COUNT(DISTINCT patientunitstayid) AS n_stays,
    APPROX_QUANTILES(IF(dose_u BETWEEN 1 AND 150, dose_u, NULL), 2)[OFFSET(1)] AS median_dose
  FROM gd GROUP BY hospitalid
)
SELECT
  COUNT(*)                                          AS n_hosp_with_glargine,
  COUNTIF(n_stays >= 20)                            AS n_hosp_ge20,
  APPROX_QUANTILES(IF(n_stays>=20, median_dose, NULL), 2)[OFFSET(1)]  AS med_of_hosp_medians,
  MIN(IF(n_stays>=20, median_dose, NULL))           AS min_hosp_median,
  APPROX_QUANTILES(IF(n_stays>=20, median_dose, NULL), 4)[OFFSET(1)]  AS p25_hosp_median,
  APPROX_QUANTILES(IF(n_stays>=20, median_dose, NULL), 4)[OFFSET(3)]  AS p75_hosp_median,
  MAX(IF(n_stays>=20, median_dose, NULL))           AS max_hosp_median
FROM per_hosp

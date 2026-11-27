-- eICU Kontrol2: glarjin doz-eksikligi MERKEZ basina (kumelenme + merkez kaybi)
WITH g AS (
  SELECT patientunitstayid, MAX(IF(REGEXP_CONTAINS(dosage, r"[0-9]"), 1, 0)) AS has_num
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
  GROUP BY patientunitstayid
),
gh AS (SELECT g.has_num, p.hospitalid FROM g JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid)),
per_hosp AS (
  SELECT hospitalid, COUNT(*) AS n_glargine, SUM(has_num) AS n_parseable,
         ROUND(100*(1 - SUM(has_num)/COUNT(*)), 1) AS pct_nullonly
  FROM gh GROUP BY hospitalid
)
SELECT
  COUNT(*) AS n_hosp_total,
  COUNTIF(n_glargine >= 20)  AS n_hosp_ge20_raw,
  COUNTIF(n_parseable >= 20) AS n_hosp_ge20_after_drop,
  COUNTIF(pct_nullonly = 100) AS n_hosp_ALL_null,
  COUNTIF(pct_nullonly >= 50)  AS n_hosp_majority_null,
  ROUND(APPROX_QUANTILES(pct_nullonly, 100)[OFFSET(50)],1) AS pct_nullonly_median_hosp,
  ROUND(APPROX_QUANTILES(pct_nullonly, 4)[OFFSET(1)],1)    AS pct_nullonly_p25,
  ROUND(APPROX_QUANTILES(pct_nullonly, 4)[OFFSET(3)],1)    AS pct_nullonly_p75,
  ROUND(MIN(pct_nullonly),1) AS min_null, ROUND(MAX(pct_nullonly),1) AS max_null
FROM per_hosp

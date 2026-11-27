-- eICU Kontrol3b: ardisik glarjin orderlar arasi saat farki (18-30h penceresi)
WITH g AS (
  SELECT DISTINCT patientunitstayid, drugstartoffset
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
    AND drugstartoffset IS NOT NULL
),
gaps AS (
  SELECT (LEAD(drugstartoffset) OVER (PARTITION BY patientunitstayid ORDER BY drugstartoffset) - drugstartoffset) AS gap_min
  FROM g
)
SELECT COUNT(*) AS n_gaps,
  APPROX_QUANTILES(gap_min,4)[OFFSET(1)] AS gap_q25_min,
  APPROX_QUANTILES(gap_min,4)[OFFSET(2)] AS gap_median_min,
  APPROX_QUANTILES(gap_min,4)[OFFSET(3)] AS gap_q75_min,
  ROUND(100*COUNTIF(gap_min BETWEEN 1080 AND 1800)/COUNT(*),1) AS pct_gap_18_30h,
  ROUND(100*COUNTIF(gap_min BETWEEN 720 AND 2160)/COUNT(*),1)  AS pct_gap_12_36h
FROM gaps WHERE gap_min IS NOT NULL

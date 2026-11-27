-- eICU Kontrol3a: stay basina glarjin order/zaman sayisi + offset kapsami
WITH g AS (
  SELECT patientunitstayid, drugstartoffset, drugorderoffset
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
),
per_stay AS (
  SELECT patientunitstayid,
    COUNT(*) AS n_orders,
    COUNT(DISTINCT drugstartoffset) AS n_distinct_times
  FROM g WHERE drugstartoffset IS NOT NULL GROUP BY patientunitstayid
)
SELECT
  (SELECT COUNT(DISTINCT patientunitstayid) FROM g) AS n_glargine_stays_all,
  (SELECT ROUND(100*COUNTIF(drugstartoffset IS NOT NULL)/COUNT(*),1) FROM g) AS pct_rows_startoffset,
  (SELECT ROUND(100*COUNTIF(drugorderoffset IS NOT NULL)/COUNT(*),1) FROM g) AS pct_rows_orderoffset,
  COUNT(*) AS n_stays_with_time,
  COUNTIF(n_distinct_times >= 2) AS n_stays_ge2_times,
  ROUND(100*COUNTIF(n_distinct_times >= 2)/COUNT(*),1) AS pct_ge2_times,
  APPROX_QUANTILES(n_distinct_times, 4)[OFFSET(1)] AS q25_times,
  APPROX_QUANTILES(n_distinct_times, 4)[OFFSET(2)] AS median_times,
  APPROX_QUANTILES(n_distinct_times, 4)[OFFSET(3)] AS q75_times,
  APPROX_QUANTILES(n_distinct_times,100)[OFFSET(95)] AS p95_times
FROM per_stay

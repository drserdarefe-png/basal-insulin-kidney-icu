-- eICU K2-EK: merkez-duzeyi doz-tamlik esigi
WITH g AS (
  SELECT patientunitstayid, MAX(IF(REGEXP_CONTAINS(dosage, r"[0-9]"),1,0)) AS has_num
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
  GROUP BY patientunitstayid),
pc AS (SELECT p.hospitalid, COUNT(*) AS n_glargine, SUM(has_num) AS n_complete,
              ROUND(100*SUM(has_num)/COUNT(*),1) AS pct_complete
       FROM g JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid) GROUP BY 1)
SELECT
  COUNTIF(n_glargine>=20)                        AS n_centers_ge20,
  COUNTIF(n_glargine>=20 AND pct_complete>=50)   AS n_ge20_and_complete50,
  COUNTIF(n_glargine>=20 AND pct_complete>=80)   AS n_ge20_and_complete80,
  ROUND(APPROX_QUANTILES(IF(n_glargine>=20,pct_complete,NULL),4)[OFFSET(1)],1) AS complete_p25,
  ROUND(APPROX_QUANTILES(IF(n_glargine>=20,pct_complete,NULL),2)[OFFSET(1)],1) AS complete_median,
  ROUND(APPROX_QUANTILES(IF(n_glargine>=20,pct_complete,NULL),4)[OFFSET(3)],1) AS complete_p75
FROM pc

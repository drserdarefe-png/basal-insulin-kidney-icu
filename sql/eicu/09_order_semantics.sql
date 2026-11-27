-- eICU K3-EK: order=doz mu titrasyon mu -- ardisik fark MEDYANI merkez basina
WITH gstay AS (
  SELECT DISTINCT patientunitstayid FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")),
cc AS (SELECT p.hospitalid, COUNT(*) AS n_stays
       FROM gstay JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid) GROUP BY 1),
g AS (SELECT DISTINCT patientunitstayid, drugstartoffset
      FROM `physionet-data.eicu_crd.medication`
      WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
        AND drugstartoffset IS NOT NULL),
gaps AS (SELECT patientunitstayid,
           (LEAD(drugstartoffset) OVER (PARTITION BY patientunitstayid ORDER BY drugstartoffset) - drugstartoffset) AS gap
         FROM g),
gh AS (SELECT p.hospitalid, gaps.gap FROM gaps
       JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid) WHERE gaps.gap IS NOT NULL),
pc AS (SELECT hospitalid, APPROX_QUANTILES(gap,2)[OFFSET(1)] AS median_gap FROM gh GROUP BY 1)
SELECT
  COUNTIF(cc.n_stays >= 20)                                          AS n_centers_ge20,
  COUNTIF(cc.n_stays >= 20 AND pc.median_gap BETWEEN 1380 AND 1500)  AS n_daily_1380_1500,
  COUNTIF(cc.n_stays >= 20 AND pc.median_gap BETWEEN 1500 AND 2000)  AS n_mid_1500_2000,
  COUNTIF(cc.n_stays >= 20 AND pc.median_gap > 2000)                 AS n_orderlevel_gt2000,
  COUNTIF(cc.n_stays >= 20 AND pc.median_gap < 1380)                 AS n_sub_daily_lt1380,
  ROUND(APPROX_QUANTILES(IF(cc.n_stays>=20, pc.median_gap, NULL),4)[OFFSET(1)],0) AS gap_p25_center,
  ROUND(APPROX_QUANTILES(IF(cc.n_stays>=20, pc.median_gap, NULL),2)[OFFSET(1)],0) AS gap_median_center,
  ROUND(APPROX_QUANTILES(IF(cc.n_stays>=20, pc.median_gap, NULL),4)[OFFSET(3)],0) AS gap_p75_center
FROM pc JOIN cc USING (hospitalid)

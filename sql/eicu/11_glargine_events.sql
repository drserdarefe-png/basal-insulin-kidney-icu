-- eICU KOL: glarjin ORDER-olayi tablosu (kovaryatlar order zamanina yoneltilmis)
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.eicu_glargine_events` AS
WITH ord AS (
  SELECT patientunitstayid, drugstartoffset,
    SAFE_CAST(REGEXP_EXTRACT(LOWER(dosage), r"^\s*([0-9]+(?:\.[0-9]+)?)\s*unit(?:\(s\)|s)?\s*$") AS FLOAT64) AS dose_u
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
    AND drugstartoffset IS NOT NULL
),
ord2 AS (SELECT patientunitstayid, drugstartoffset, MAX(dose_u) AS dose_u FROM ord GROUP BY 1,2),
cr AS (
  SELECT patientunitstayid, labresultoffset, labresult AS scr
  FROM `physionet-data.eicu_crd.lab`
  WHERE LOWER(labname) = "creatinine" AND labresult BETWEEN 0.1 AND 20
),
joined AS (
  SELECT o.patientunitstayid, o.drugstartoffset, o.dose_u, c.scr
  FROM ord2 o
  LEFT JOIN cr c
    ON o.patientunitstayid = c.patientunitstayid
   AND c.labresultoffset <= o.drugstartoffset
   AND o.drugstartoffset - c.labresultoffset <= 2880
  QUALIFY ROW_NUMBER() OVER (PARTITION BY o.patientunitstayid, o.drugstartoffset ORDER BY c.labresultoffset DESC) = 1
),
feat AS (
  SELECT j.patientunitstayid, j.drugstartoffset, j.dose_u, j.scr,
    p.hospitalid, p.admissionweight AS weight, p.gender,
    IF(p.age = "> 89", 91.4, SAFE_CAST(p.age AS FLOAT64)) AS age,
    (p.age = "> 89") AS age_topcoded,
    (dm.patientunitstayid IS NOT NULL) AS dm
  FROM joined j
  JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid)
  LEFT JOIN (SELECT DISTINCT patientunitstayid FROM `physionet-data.eicu_crd.pasthistory`
             WHERE LOWER(pasthistorypath) LIKE "%diabetes%" OR LOWER(pasthistoryvalue) LIKE "%diabetes%") dm
    USING (patientunitstayid)
  WHERE p.admissionweight BETWEEN 30 AND 300
),
eg AS (
  SELECT *, SAFE_DIVIDE(dose_u, weight) AS dose_u_kg,
    CASE WHEN scr IS NULL OR age IS NULL OR gender NOT IN ("Male","Female") THEN NULL ELSE
      142.0 * POW(LEAST(scr/IF(gender="Female",0.7,0.9),1.0), IF(gender="Female",-0.241,-0.302))
            * POW(GREATEST(scr/IF(gender="Female",0.7,0.9),1.0), -1.200)
            * POW(0.9938, age) * IF(gender="Female",1.012,1.0) END AS egfr
  FROM feat
)
SELECT *,
  CASE WHEN egfr IS NULL THEN NULL WHEN egfr<15 THEN "lt15" WHEN egfr<45 THEN "b15_45"
       WHEN egfr<90 THEN "b45_90" ELSE "ge90" END AS egfr_band,
  ROW_NUMBER() OVER w AS order_seq,
  (ROW_NUMBER() OVER w = 1) AS is_first_dose,
  LAG(dose_u_kg) OVER w AS prev_dose_u_kg,
  (drugstartoffset - LAG(drugstartoffset) OVER w)/60.0 AS hrs_since_prev
FROM eg
WINDOW w AS (PARTITION BY patientunitstayid ORDER BY drugstartoffset)

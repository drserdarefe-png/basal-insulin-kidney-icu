-- eICU IS4: glarjin hastalarinda 4 degisken kapsami
WITH g AS (
  SELECT DISTINCT patientunitstayid
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
),
cr AS (
  SELECT DISTINCT patientunitstayid FROM `physionet-data.eicu_crd.lab`
  WHERE LOWER(labname) LIKE "%creatinine%"
),
dx AS (
  SELECT DISTINCT patientunitstayid FROM `physionet-data.eicu_crd.diagnosis`
  WHERE LOWER(diagnosisstring) LIKE "%diabetes%"
),
ph AS (
  SELECT DISTINCT patientunitstayid FROM `physionet-data.eicu_crd.pasthistory`
  WHERE LOWER(pasthistorypath) LIKE "%diabetes%" OR LOWER(pasthistoryvalue) LIKE "%diabetes%"
)
SELECT
  COUNT(*)                                                          AS n_glargine,
  ROUND(100*COUNTIF(p.admissionweight IS NOT NULL)/COUNT(*),1)      AS pct_weight,
  ROUND(100*COUNTIF(cr.patientunitstayid IS NOT NULL)/COUNT(*),1)   AS pct_creatinine,
  ROUND(100*COUNTIF(p.age IS NOT NULL AND p.age != "")/COUNT(*),1)  AS pct_age,
  COUNTIF(p.age = "> 89")                                           AS n_age_over89,
  ROUND(100*COUNTIF(p.gender IS NOT NULL AND p.gender != "")/COUNT(*),1) AS pct_gender,
  ROUND(100*COUNTIF(dx.patientunitstayid IS NOT NULL)/COUNT(*),1)   AS pct_dm_diagnosis,
  ROUND(100*COUNTIF(ph.patientunitstayid IS NOT NULL)/COUNT(*),1)   AS pct_dm_pasthistory,
  ROUND(100*COUNTIF(dx.patientunitstayid IS NOT NULL OR ph.patientunitstayid IS NOT NULL)/COUNT(*),1) AS pct_dm_either
FROM g
JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid)
LEFT JOIN cr USING (patientunitstayid)
LEFT JOIN dx USING (patientunitstayid)
LEFT JOIN ph USING (patientunitstayid)

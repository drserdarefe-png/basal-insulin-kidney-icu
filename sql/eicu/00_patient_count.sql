-- eICU IS1: kohort buyuklugu + kapsam
SELECT
  COUNT(*)                              AS n_rows,
  COUNT(DISTINCT patientunitstayid)     AS n_unitstays,
  COUNT(DISTINCT patienthealthsystemstayid) AS n_hospstays,
  COUNT(DISTINCT uniquepid)             AS n_patients,
  COUNT(DISTINCT hospitalid)            AS n_hospitals
FROM `physionet-data.eicu_crd.patient`

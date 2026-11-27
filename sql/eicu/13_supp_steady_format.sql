-- eICU A6/A7 supp: steady-state (ikinci kreatinin) + doz format tier
WITH ev AS (
  SELECT * FROM `uludag-mimic-insulin.insulin.eicu_glargine_events`
  WHERE dose_u_kg IS NOT NULL AND egfr IS NOT NULL
),
cr AS (
  SELECT patientunitstayid, labresultoffset, labresult
  FROM `physionet-data.eicu_crd.lab`
  WHERE LOWER(labname) = "creatinine" AND labresult BETWEEN 0.1 AND 20
),
crpair AS (
  SELECT e.patientunitstayid, e.drugstartoffset,
    ARRAY_AGG(c.labresult ORDER BY c.labresultoffset DESC LIMIT 2) AS scrs
  FROM ev e JOIN cr c
    ON c.patientunitstayid = e.patientunitstayid
   AND c.labresultoffset <= e.drugstartoffset
   AND e.drugstartoffset - c.labresultoffset <= 2880
  GROUP BY 1, 2
),
fmt AS (
  SELECT patientunitstayid, drugstartoffset,
    MAX(IF(REGEXP_CONTAINS(LOWER(dosage), r"^\s*[0-9]+(\.[0-9]+)?\s*units?\s*$"), 1, 0)) AS tier1_exact
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
    AND drugstartoffset IS NOT NULL
  GROUP BY 1, 2
)
SELECT e.patientunitstayid, e.hospitalid, e.dose_u_kg, e.egfr, e.egfr_band,
  e.prev_dose_u_kg, e.hrs_since_prev, e.age, e.weight, CAST(e.dm AS INT64) AS dm, e.is_first_dose,
  ARRAY_LENGTH(cp.scrs) AS n_scr,
  IF(ARRAY_LENGTH(cp.scrs) = 2 AND cp.scrs[OFFSET(0)] > 0,
     ABS(cp.scrs[OFFSET(0)] - cp.scrs[OFFSET(1)]) / cp.scrs[OFFSET(0)], NULL) AS dcr_frac,
  fmt.tier1_exact
FROM ev e
LEFT JOIN crpair cp USING (patientunitstayid, drugstartoffset)
LEFT JOIN fmt USING (patientunitstayid, drugstartoffset)

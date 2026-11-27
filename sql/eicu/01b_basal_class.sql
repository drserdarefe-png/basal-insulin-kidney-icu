-- eICU IS2b: bazal SINIF bazinda DISTINCT hasta (varyant cakismasi duzeltildi)
WITH med AS (
  SELECT patientunitstayid, LOWER(drugname) AS dn
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|levemir|detemir|toujeo|degludec|tresiba|basaglar|semglee|nph|humulin n|novolin n")
),
cls AS (
  SELECT patientunitstayid,
    CASE
      WHEN REGEXP_CONTAINS(dn, r"glargine|lantus|toujeo|basaglar|semglee") THEN "glargine"
      WHEN REGEXP_CONTAINS(dn, r"detemir|levemir")                          THEN "detemir"
      WHEN REGEXP_CONTAINS(dn, r"degludec|tresiba")                         THEN "degludec"
      WHEN REGEXP_CONTAINS(dn, r"nph|humulin n|novolin n")                  THEN "nph"
      ELSE "other" END AS basal_class
  FROM med
)
SELECT c.basal_class,
       COUNT(DISTINCT c.patientunitstayid) AS n_stays,
       COUNT(DISTINCT p.hospitalid)        AS n_hosp
FROM cls c JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid)
GROUP BY c.basal_class
ORDER BY n_stays DESC

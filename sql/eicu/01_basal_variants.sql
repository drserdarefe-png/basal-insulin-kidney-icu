-- eICU IS2: bazal insulin drugname VARYANTLARI (medication + infusiondrug)
-- drugname SERBEST METIN. Varyant listesi ref tablosunun temeli.
WITH base AS (
  SELECT "medication" AS src, patientunitstayid, LOWER(drugname) AS dn
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|levemir|detemir|toujeo|degludec|tresiba|basaglar|semglee|nph|humulin n|novolin n")
  UNION ALL
  SELECT "infusiondrug" AS src, patientunitstayid, LOWER(drugname) AS dn
  FROM `physionet-data.eicu_crd.infusiondrug`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|levemir|detemir|toujeo|degludec|tresiba|basaglar|semglee|nph|humulin n|novolin n")
)
SELECT b.src, b.dn AS drugname_lc,
       COUNT(*)                          AS n_rows,
       COUNT(DISTINCT b.patientunitstayid) AS n_stays,
       COUNT(DISTINCT p.hospitalid)      AS n_hosp
FROM base b
JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid)
GROUP BY b.src, b.dn
ORDER BY n_stays DESC

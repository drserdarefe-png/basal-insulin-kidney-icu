-- eICU K1-EK: NPH marka-adi (insulin SARTI OLMADAN) -- dairesel kacis testi
SELECT LOWER(drugname) AS dn, COUNT(*) AS n, COUNT(DISTINCT patientunitstayid) AS n_stay
FROM `physionet-data.eicu_crd.medication`
WHERE REGEXP_CONTAINS(LOWER(drugname), r"humulin|novolin|lantus|levemir|toujeo|tresiba|basaglar|semglee|isophane|nph|70/30|75/25|ryzodeg|apidra|admelog")
GROUP BY dn ORDER BY n_stay DESC

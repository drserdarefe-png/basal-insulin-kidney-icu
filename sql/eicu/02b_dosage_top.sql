SELECT dosage, COUNT(*) AS n, COUNT(DISTINCT patientunitstayid) AS n_stays
FROM `physionet-data.eicu_crd.medication`
WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
GROUP BY dosage ORDER BY n DESC LIMIT 30

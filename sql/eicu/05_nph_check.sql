-- eICU Kontrol1a: genisletilmis NPH/mix regex
SELECT LOWER(drugname) AS dn,
       COUNT(*) AS n_rows,
       COUNT(DISTINCT m.patientunitstayid) AS n_stays,
       COUNT(DISTINCT p.hospitalid) AS n_hosp
FROM `physionet-data.eicu_crd.medication` m
JOIN `physionet-data.eicu_crd.patient` p USING (patientunitstayid)
WHERE REGEXP_CONTAINS(LOWER(drugname), r"isophane|humulin.?n|novolin.?n|nph|70.?30|75.?25|humalog mix|novolog mix|ryzodeg")
GROUP BY dn ORDER BY n_stays DESC

-- eICU Kontrol1b: insulin gecen TUM varyantlar (ne kaciriyoruz)
SELECT LOWER(drugname) AS dn,
       COUNT(DISTINCT patientunitstayid) AS n_stays
FROM `physionet-data.eicu_crd.medication`
WHERE LOWER(drugname) LIKE "%insulin%"
GROUP BY dn ORDER BY n_stays DESC LIMIT 100

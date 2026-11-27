-- eICU IS3: glarjin dosage ham degerleri + parse edilebilirlik (KRITIK)
WITH g AS (
  SELECT patientunitstayid, dosage
  FROM `physionet-data.eicu_crd.medication`
  WHERE REGEXP_CONTAINS(LOWER(drugname), r"glargine|lantus|toujeo|basaglar|semglee")
),
classified AS (
  SELECT dosage,
    CASE
      WHEN dosage IS NULL OR TRIM(dosage) = "" THEN "z_null_empty"
      WHEN REGEXP_CONTAINS(LOWER(dosage), r"/\s*kg")                       THEN "per_kg"
      WHEN REGEXP_CONTAINS(LOWER(dosage), r"^\s*[0-9]+(\.[0-9]+)?\s*units?\s*$") THEN "N_units"
      WHEN REGEXP_CONTAINS(LOWER(dosage), r"^\s*[0-9]+(\.[0-9]+)?\s*u\s*$")      THEN "N_u"
      WHEN REGEXP_CONTAINS(dosage, r"[0-9]")                               THEN "has_number_other"
      ELSE "z_unparseable_nonnum" END AS parse_class
  FROM g
)
SELECT parse_class,
       COUNT(*)                          AS n,
       ROUND(100*COUNT(*)/SUM(COUNT(*)) OVER (), 1) AS pct
FROM classified
GROUP BY parse_class
ORDER BY n DESC

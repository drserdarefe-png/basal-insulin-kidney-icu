-- =====================================================================
-- 01b — EXC-03: INSULIN POMPASI / BELGELENMEMIS MARUZIYET      v2.0
-- Registry v1.1 EXC-03: ICD (Z96.41/V45.85) VE ICU dozu=0 VE glu_max>180
-- Uygulanmis: PATCH-05 (ref.insulin_pk TEK KAYNAK; label LIKE ayrismasi YOK)
--             BULGU-07 (dongusel bagimlilik yok; pump_icd_any passthrough;
--                       ADIM 2: cohort_final EXC-03 uygulanarak kurulur)
-- CIKTI: pump_status, attrition, cohort_final, grid
-- Cerceve: kriter "pompa var mi" DEGIL, "insulin maruziyeti belgelenmis mi".
-- =====================================================================

-- ---------------------------------------------------------------------
-- ADIM 0 — SINIFLAMA + PREVALANS
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.pump_status` AS
WITH pump_icd AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE icd_code IN ('Z9641', 'V4585')       -- pompa varligi (ICD-10 / ICD-9)
),
ins_doc AS (
  -- PATCH-05: ref TEK KAYNAK. label LIKE KULLANILMAZ (03 ile ayrisma riski).
  SELECT ie.stay_id,
         SUM(ie.amount) AS total_units,
         COUNT(*)       AS n_ins_events
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN `uludag-mimic-insulin.ref.insulin_pk` pk ON pk.itemid = ie.itemid
  WHERE ie.statusdescription NOT IN ('Rewritten')
    AND ie.amountuom IN ('units','U') AND ie.amount > 0
    AND pk.verified = TRUE
  GROUP BY 1
),
glu_hi AS (
  SELECT stay_id, MAX(valuenum) AS glu_max
  FROM `physionet-data.mimiciv_3_1_icu.chartevents`
  WHERE itemid IN (225664, 220621, 226537) AND valuenum BETWEEN 10 AND 2000
  GROUP BY 1
),
dm AS (
  SELECT hadm_id, LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^E1[0-4]|^250')) AS dm_any
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` GROUP BY 1
)
SELECT
  c.stay_id, c.hadm_id,
  (p.hadm_id IS NOT NULL)    AS pump_icd_any,   -- BULGU-07(3): passthrough
  COALESCE(i.total_units, 0) AS ins_units,
  COALESCE(g.glu_max, 0)     AS glu_max,
  COALESCE(d.dm_any, FALSE)  AS dm_any,
  CASE
    -- A: hayalet doz imzasi -> DISLA (EXC-03)
    WHEN p.hadm_id IS NOT NULL AND COALESCE(i.total_units, 0) = 0
     AND COALESCE(g.glu_max, 0) > 180             THEN 'A_pump_undocumented'
    -- B: pompa cikarilmis, doz belgeleniyor -> TUT + kovaryat (agirlikli T1DM)
    WHEN p.hadm_id IS NOT NULL AND COALESCE(i.total_units, 0) > 0
                                                  THEN 'B_prior_pump_converted'
    -- C: pompa yok ama ayni hedef-kirlenmesi riski -> ANA ANALIZDE TUT,
    --    duyarlilik analizinde cikar. 06'da dose_bin='0_no_dose' kirletir.
    WHEN p.hadm_id IS NULL AND COALESCE(i.total_units, 0) = 0
     AND COALESCE(g.glu_max, 0) > 180
     AND COALESCE(d.dm_any, FALSE)                THEN 'C_suspect_undocumented'
    WHEN COALESCE(i.total_units, 0) = 0           THEN 'D_true_no_dose'
    ELSE 'E_normal'
  END AS pump_class
FROM `uludag-mimic-insulin.insulin.cohort_raw` c
LEFT JOIN pump_icd p ON p.hadm_id = c.hadm_id
LEFT JOIN ins_doc  i ON i.stay_id = c.stay_id
LEFT JOIN glu_hi   g ON g.stay_id = c.stay_id
LEFT JOIN dm       d ON d.hadm_id = c.hadm_id;

-- ---- PREVALANS RAPORU -> qc/01b_prevalence.txt ----------------------
-- SELECT pump_class, COUNT(*) n,
--        ROUND(100.0*COUNT(*)/SUM(COUNT(*)) OVER (), 3) AS pct,
--        ROUND(AVG(glu_max),0)   AS mean_glu_max,
--        ROUND(AVG(ins_units),1) AS mean_units
-- FROM `uludag-mimic-insulin.insulin.pump_status`
-- GROUP BY 1 ORDER BY n DESC;
--
-- KARAR ESIGI (registry v1.2'ye yazilacak):
--   A < %0.5 -> sadece A disla, tek satir duyarlilik yeter
--   A > %2   -> CONSORT'ta ayri kutu
--   C > %5   -> 06'da dose_bin='0_no_dose' KIRLI, ayri bayrak

-- ---------------------------------------------------------------------
-- ADIM 1 — ATTRITION (CONSORT dogrudan buradan)
-- ---------------------------------------------------------------------
CREATE OR REPLACE VIEW `uludag-mimic-insulin.insulin.attrition` AS
SELECT
  COUNT(*)                     AS n_all_icu_stays,
  COUNTIF(c.exc_01_dka)        AS n_exc_01_dka,
  COUNTIF(c.exc_02_pregnancy)  AS n_exc_02_pregnancy,
  COUNTIF(ps.pump_class = 'A_pump_undocumented') AS n_exc_03_pump,
  COUNTIF(c.exc_04_los)        AS n_exc_04_los,
  COUNTIF(c.exc_05_few_glu)    AS n_exc_05_few_glu,
  COUNTIF(c.exc_06_alf)        AS n_exc_06_alf,
  COUNTIF(c.exc_07_liver_tx)   AS n_exc_07_liver_tx,
  COUNTIF(c.exc_08_panc)       AS n_exc_08_panc,
  COUNTIF(c.exc_09_no_nut_doc) AS n_exc_09_no_nut_doc,
  COUNTIF(c.exc_10_bad_uom)    AS n_exc_10_bad_uom,
  COUNTIF(c.exc_11_readmit)    AS n_exc_11_readmit,
  COUNTIF(NOT (c.exc_01_dka OR c.exc_02_pregnancy
            OR ps.pump_class = 'A_pump_undocumented'
            OR c.exc_04_los OR c.exc_05_few_glu OR c.exc_06_alf
            OR c.exc_07_liver_tx OR c.exc_08_panc OR c.exc_09_no_nut_doc
            OR c.exc_10_bad_uom OR c.exc_11_readmit)) AS n_analytic
FROM `uludag-mimic-insulin.insulin.cohort_raw` c
JOIN `uludag-mimic-insulin.insulin.pump_status` ps USING (stay_id);

-- ---------------------------------------------------------------------
-- ADIM 2 — COHORT_FINAL (BULGU-07(4): EXC-03 uygulanarak)
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.cohort_final`
CLUSTER BY stay_id AS
SELECT
  c.subject_id, c.hadm_id, c.stay_id, c.intime, c.outtime, c.los,
  ps.pump_class,
  ps.pump_class = 'B_prior_pump_converted' AS cov_prior_pump,     -- 05 kovaryati
  ps.pump_class = 'C_suspect_undocumented' AS sens_suspect_undoc, -- duyarlilik
  ps.pump_icd_any
FROM `uludag-mimic-insulin.insulin.cohort_raw` c
JOIN `uludag-mimic-insulin.insulin.pump_status` ps USING (stay_id)
WHERE NOT (c.exc_01_dka OR c.exc_02_pregnancy
        OR ps.pump_class = 'A_pump_undocumented'
        OR c.exc_04_los OR c.exc_05_few_glu OR c.exc_06_alf
        OR c.exc_07_liver_tx OR c.exc_08_panc OR c.exc_09_no_nut_doc
        OR c.exc_10_bad_uom OR c.exc_11_readmit);

-- ---------------------------------------------------------------------
-- ADIM 3 — SAATLIK GRID (ilk 14 gune kirp: uzun kuyruk seyrek + maliyet)
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.grid`
CLUSTER BY stay_id AS
SELECT
  c.subject_id, c.hadm_id, c.stay_id, hr,
  TIMESTAMP_ADD(c.intime, INTERVAL hr HOUR) AS ts
FROM `uludag-mimic-insulin.insulin.cohort_final` c,
UNNEST(GENERATE_ARRAY(0, LEAST(CAST(FLOOR(c.los * 24) AS INT64), 14*24))) AS hr;

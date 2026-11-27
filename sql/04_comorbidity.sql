-- =====================================================================
-- 04 — KOMORBIDITE BLOGU (effect modifier'lar)              v2.0
-- Uygulanmis: BULGU-01 (MELD kirpma sirasi *10 SONRASINA + UNOS 2016)
--             BULGU-04 (cr_age_hr = MIN, MAX degil)
--             BULGU-05 (CKD-EPI 2021 sabiti 142, 141 degil)
--             PATCH-02 (Child-Pugh YALNIZ sirotikte tanimli)
--             PATCH-04 (MELD kreatinin RRT'de 4.0'a sabit + tanim beyani)
--             PATCH-10 (steroid: prescriptions tabanli; ICU modulunde YOK)
-- =====================================================================

-- =====================================================================
-- 4A.1 BASELINE KREATININ -- MIMIC'te YOK. Uc proxy, hepsi saklanir.
-- =====================================================================
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.baseline_cr` AS
WITH pre_adm AS (        -- proxy A: yatistan onceki 365 gun medyani (EN IYI)
  SELECT c.subject_id, c.stay_id,
         APPROX_QUANTILES(le.valuenum, 2)[OFFSET(1)] AS cr_base_prior,
         COUNT(*) AS n_prior
  FROM `uludag-mimic-insulin.insulin.cohort_final` c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.subject_id = c.subject_id AND le.itemid = 50912
   AND le.charttime BETWEEN TIMESTAMP_SUB(c.intime, INTERVAL 365 DAY)
                        AND TIMESTAMP_SUB(c.intime, INTERVAL 1 DAY)
   AND le.valuenum BETWEEN 0.1 AND 20
  GROUP BY 1,2
),
in_stay_min AS (         -- proxy B: ilk 7 gunun NADIRI
  SELECT c.subject_id, c.stay_id, MIN(le.valuenum) AS cr_base_nadir
  FROM `uludag-mimic-insulin.insulin.cohort_final` c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.subject_id = c.subject_id AND le.itemid = 50912
   AND le.charttime BETWEEN c.intime AND TIMESTAMP_ADD(c.intime, INTERVAL 7 DAY)
   AND le.valuenum BETWEEN 0.1 AND 20
  GROUP BY 1,2
),
demo AS (
  SELECT c.stay_id, c.subject_id, p.gender,
         p.anchor_age + (EXTRACT(YEAR FROM c.intime) - p.anchor_year) AS age
  FROM `uludag-mimic-insulin.insulin.cohort_final` c
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` p USING (subject_id)
)
SELECT
  d.stay_id, d.subject_id, d.gender, d.age,
  pa.cr_base_prior, pa.n_prior, ism.cr_base_nadir,
  -- proxy C: MDRD geri-hesabi (eGFR=75 varsayimi) -- ikisi de yoksa
  CASE WHEN pa.cr_base_prior IS NULL AND ism.cr_base_nadir IS NULL THEN
    POW(75.0 / (175 * POW(d.age, -0.203) * IF(d.gender='F', 0.742, 1.0)), 1/-1.154)
  END AS cr_base_mdrd_backcalc,
  COALESCE(pa.cr_base_prior, ism.cr_base_nadir) AS cr_baseline,
  CASE WHEN pa.cr_base_prior  IS NOT NULL THEN 'prior_median'
       WHEN ism.cr_base_nadir IS NOT NULL THEN 'stay_nadir'
       ELSE 'mdrd_backcalc' END AS cr_baseline_method
FROM demo d
LEFT JOIN pre_adm pa USING (stay_id)
LEFT JOIN in_stay_min ism USING (stay_id);
-- DUYARLILIK ANALIZI ZORUNLU: cr_baseline_method'a gore stratifiye tekrarla.

-- =====================================================================
-- 4A.2 KRONIK BOBREK HASTALIGI -- ICD *VE* baseline kreatinin KESISIMI
--      ICD tek basina AKI/KBY ayirmiyor.
-- =====================================================================
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.renal_static` AS
WITH icd AS (
  SELECT hadm_id,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(N18|585)'))     AS icd_ckd,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(N186|5856|Z992|V451|V560)')) AS icd_esrd,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(N17|584)'))     AS icd_aki
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` GROUP BY 1
),
egfr AS (
  SELECT b.stay_id, b.cr_baseline, b.age, b.gender,
    -- BULGU-05: CKD-EPI 2021 IRKSIZ (Inker NEJM 2021) -> sabit 142 (141 DEGIL)
    142
    * POW(LEAST(b.cr_baseline / IF(b.gender='F', 0.7, 0.9), 1),
          IF(b.gender='F', -0.241, -0.302))
    * POW(GREATEST(b.cr_baseline / IF(b.gender='F', 0.7, 0.9), 1), -1.200)
    * POW(0.9938, b.age)
    * IF(b.gender='F', 1.012, 1.0) AS egfr_baseline
  FROM `uludag-mimic-insulin.insulin.baseline_cr` b
  WHERE b.cr_baseline IS NOT NULL
)
SELECT
  c.stay_id, c.hadm_id,
  e.cr_baseline, e.egfr_baseline,
  CASE WHEN e.egfr_baseline >= 90 THEN 'G1'
       WHEN e.egfr_baseline >= 60 THEN 'G2'
       WHEN e.egfr_baseline >= 45 THEN 'G3a'
       WHEN e.egfr_baseline >= 30 THEN 'G3b'
       WHEN e.egfr_baseline >= 15 THEN 'G4'
       WHEN e.egfr_baseline <  15 THEN 'G5'
       ELSE 'unknown' END AS ckd_stage_baseline,
  (COALESCE(i.icd_ckd, FALSE) AND e.egfr_baseline < 60)     AS ckd_confirmed,
  (COALESCE(i.icd_ckd, FALSE) AND e.egfr_baseline >= 60)    AS ckd_icd_only_discordant,
  (NOT COALESCE(i.icd_ckd, FALSE) AND e.egfr_baseline < 60) AS ckd_lab_only_discordant,
  COALESCE(i.icd_esrd, FALSE) AS esrd,
  COALESCE(i.icd_aki,  FALSE) AS icd_aki_present
FROM `uludag-mimic-insulin.insulin.cohort_final` c
LEFT JOIN icd i USING (hadm_id)
LEFT JOIN egfr e USING (stay_id);

-- =====================================================================
-- 4A.3 SAATLIK RENAL
-- =====================================================================
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.renal_hourly`
CLUSTER BY stay_id AS
WITH cr_ts AS (
  SELECT c.stay_id, le.charttime, le.valuenum AS cr
  FROM `uludag-mimic-insulin.insulin.cohort_final` c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.subject_id = c.subject_id AND le.itemid = 50912
   AND le.charttime BETWEEN c.intime AND c.outtime
   AND le.valuenum BETWEEN 0.1 AND 20
),
cr_grid AS (
  SELECT g.stay_id, g.hr, g.ts,
    ARRAY_AGG(cr_ts.cr ORDER BY cr_ts.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS cr_last,
    -- BULGU-04: EN TAZE kreatininin yasi -> MIN (MAX en ESKIyi verirdi)
    MIN(TIMESTAMP_DIFF(g.ts, cr_ts.charttime, HOUR)) AS cr_age_hr
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN cr_ts
    ON cr_ts.stay_id = g.stay_id
   AND cr_ts.charttime <= g.ts
   AND cr_ts.charttime > TIMESTAMP_SUB(g.ts, INTERVAL 48 HOUR)
  GROUP BY 1,2,3
),
rrt AS (
  SELECT g.stay_id, g.hr,
    MAX(IF(cr.stay_id IS NOT NULL, 1, 0)) AS on_crrt
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.crrt` cr
    ON cr.stay_id = g.stay_id
   AND cr.charttime BETWEEN TIMESTAMP_SUB(g.ts, INTERVAL 1 HOUR) AND g.ts
  GROUP BY 1,2
),
kdigo AS (
  SELECT g.stay_id, g.hr, MAX(k.aki_stage) AS kdigo_stage
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.kdigo_stages` k
    ON k.stay_id = g.stay_id
   AND k.charttime BETWEEN TIMESTAMP_SUB(g.ts, INTERVAL 6 HOUR) AND g.ts
  GROUP BY 1,2
),
bun AS (          -- INT-03: uremik insulin direnci
  SELECT g.stay_id, g.hr,
    ARRAY_AGG(le.valuenum ORDER BY le.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS bun
  FROM `uludag-mimic-insulin.insulin.grid` g
  JOIN `uludag-mimic-insulin.insulin.cohort_final` c USING (stay_id)
  LEFT JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.subject_id = c.subject_id AND le.itemid = 51006
   AND le.charttime <= g.ts AND le.charttime > TIMESTAMP_SUB(g.ts, INTERVAL 48 HOUR)
  GROUP BY 1,2
)
SELECT
  cg.stay_id, cg.hr, cg.ts, cg.cr_last, cg.cr_age_hr,
  CASE WHEN cg.cr_last IS NOT NULL THEN
    142                                        -- BULGU-05
    * POW(LEAST(cg.cr_last / IF(d.gender='F',0.7,0.9), 1), IF(d.gender='F',-0.241,-0.302))
    * POW(GREATEST(cg.cr_last / IF(d.gender='F',0.7,0.9), 1), -1.200)
    * POW(0.9938, d.age) * IF(d.gender='F', 1.012, 1.0)
  END AS egfr_now,
  COALESCE(k.kdigo_stage, 0) AS kdigo_stage,
  COALESCE(r.on_crrt, 0)     AS on_crrt,
  b.bun,
  -- eGFR AKI'de gecersiz: kreatinin hizla degisiyorsa kararli hal YOK
  IF(COALESCE(k.kdigo_stage,0) > 0, 1, 0) AS egfr_unsteady_flag
FROM cr_grid cg
JOIN `uludag-mimic-insulin.insulin.baseline_cr` d USING (stay_id)
LEFT JOIN kdigo k USING (stay_id, hr)
LEFT JOIN rrt   r USING (stay_id, hr)
LEFT JOIN bun   b USING (stay_id, hr);

-- =====================================================================
-- 4B HEPATIK
-- MELD-Na TANIMI: UNOS 2016 (OPTN Policy 9). Kim 2008 NEJM formu DEGIL.
--   Bilincli secim, SAP'ta beyan edildi. Kim 2008 duyarlilik analizinde.
--   Fark dusuk-MELD bolgesinde: kompanse sirotik Kim=14.5 vs UNOS=11.0
--   Yuvarlama KRITIK: raw=11.28 -> yuvarlanmis 11 (duzeltme YOK) = 11.0
--                                  yuvarlanmamis 11.28 > 11        = 17.9
--   OPTN tam sayi MELD kullanir -> ROUND() ACIK olmali.
-- =====================================================================
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.hepatic_static` AS
WITH icd AS (
  SELECT hadm_id,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(K703|K717|K74|K766|K765)')
            OR REGEXP_CONTAINS(icd_code, r'^(5712|5715|5716|5723)'))    AS cirrhosis,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(I85|I864|K766)')
            OR REGEXP_CONTAINS(icd_code, r'^(4560|4561|4562)'))         AS portal_htn,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(R18|K7031)')
            OR REGEXP_CONTAINS(icd_code, r'^(7895|5715)'))              AS ascites_icd,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(K7282|K729)')
            OR REGEXP_CONTAINS(icd_code, r'^(5722)'))                   AS he_icd,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^(E1[0-4])|^250'))           AS dm_any,
    LOGICAL_OR(REGEXP_CONTAINS(icd_code, r'^E10|^2500[13]'))            AS dm_t1
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` GROUP BY 1
),
he_rx AS (
  SELECT DISTINCT hadm_id, TRUE AS he_rx
  FROM `physionet-data.mimiciv_3_1_hosp.prescriptions`
  WHERE REGEXP_CONTAINS(LOWER(drug), r'lactulose|rifaximin')
),
-- PATCH-04: OPTN kurali -- RRT'de kreatinin 4.0'a SABITLENIR
rrt_24h AS (
  SELECT DISTINCT c.stay_id, TRUE AS on_rrt_24h
  FROM `uludag-mimic-insulin.insulin.cohort_final` c
  JOIN `physionet-data.mimiciv_3_1_derived.crrt` cr
    ON cr.stay_id = c.stay_id
   AND cr.charttime BETWEEN c.intime AND TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
),
labs AS (
  SELECT c.stay_id,
    MAX(IF(le.itemid = 50885, le.valuenum, NULL)) AS bilirubin,
    MIN(IF(le.itemid = 50862, le.valuenum, NULL)) AS albumin,
    MAX(IF(le.itemid = 51237, le.valuenum, NULL)) AS inr,
    MAX(IF(le.itemid = 50912, le.valuenum, NULL)) AS cr_meld,
    MIN(IF(le.itemid IN (50983, 52623), le.valuenum, NULL)) AS sodium
  FROM `uludag-mimic-insulin.insulin.cohort_final` c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.subject_id = c.subject_id
   AND le.charttime BETWEEN c.intime AND TIMESTAMP_ADD(c.intime, INTERVAL 24 HOUR)
  GROUP BY 1
),
meld_calc AS (
  SELECT l.stay_id,
    -- PATCH-04: RRT -> cr = 4.0
    IF(COALESCE(r.on_rrt_24h, FALSE), 4.0, l.cr_meld) AS cr_meld_adj,
    COALESCE(r.on_rrt_24h, FALSE) AS meld_cr_capped_rrt,
    l.bilirubin, l.inr, l.sodium
  FROM labs l
  LEFT JOIN rrt_24h r USING (stay_id)
),
meld_raw AS (
  SELECT stay_id, meld_cr_capped_rrt, sodium,
    CASE WHEN bilirubin IS NULL OR inr IS NULL OR cr_meld_adj IS NULL THEN NULL
    ELSE ROUND(                                    -- OPTN: TAM SAYI MELD
      LEAST(GREATEST(
        ( 0.957 * LN(GREATEST(LEAST(cr_meld_adj, 4.0), 1.0))
        + 0.378 * LN(GREATEST(bilirubin, 1.0))
        + 1.120 * LN(GREATEST(inr, 1.0))
        + 0.643 ) * 10                             -- BULGU-01: *10 ONCE
      , 6.0), 40.0)                                --           kirpma SONRA
    ) END AS meld_int
  FROM meld_calc
)
SELECT
  c.stay_id, c.hadm_id,
  COALESCE(i.cirrhosis, FALSE)   AS cirrhosis,
  COALESCE(i.portal_htn, FALSE)  AS portal_htn,
  COALESCE(i.ascites_icd, FALSE) AS ascites,
  COALESCE(i.he_icd, FALSE) OR COALESCE(h.he_rx, FALSE) AS hepatic_enceph,
  COALESCE(i.dm_any, FALSE) AS dm_any,
  COALESCE(i.dm_t1, FALSE)  AS dm_t1,
  l.bilirubin, l.albumin, l.inr, l.sodium,
  m.meld_cr_capped_rrt,
  -- PATCH-02: Child-Pugh YALNIZ sirotikte tanimli.
  -- Non-sirotik septik hastada bili/INR/albumin akut hastalikla bozulur;
  -- skor hesaplanabilir ama ANLAMSIZDIR -> NULL.
  CASE WHEN NOT COALESCE(i.cirrhosis, FALSE) THEN NULL
       WHEN l.bilirubin IS NULL OR l.albumin IS NULL OR l.inr IS NULL THEN NULL
  ELSE
      (CASE WHEN l.bilirubin < 2 THEN 1 WHEN l.bilirubin <= 3 THEN 2 ELSE 3 END)
    + (CASE WHEN l.albumin > 3.5 THEN 1 WHEN l.albumin >= 2.8 THEN 2 ELSE 3 END)
    + (CASE WHEN l.inr < 1.7 THEN 1 WHEN l.inr <= 2.3 THEN 2 ELSE 3 END)
    + (CASE WHEN COALESCE(i.ascites_icd,FALSE) THEN 2 ELSE 1 END)
    + (CASE WHEN COALESCE(i.he_icd,FALSE) OR COALESCE(h.he_rx,FALSE) THEN 2 ELSE 1 END)
  END AS child_pugh_score,
  -- UNOS 2016 MELD-Na: Na duzeltmesi YALNIZ MELD > 11 iken
  CASE WHEN m.meld_int IS NULL THEN NULL
       WHEN m.meld_int > 11 THEN
         LEAST(GREATEST(
           m.meld_int
           + 1.32  * (137 - LEAST(GREATEST(COALESCE(m.sodium,137), 125), 137))
           - 0.033 * m.meld_int
                   * (137 - LEAST(GREATEST(COALESCE(m.sodium,137), 125), 137))
         , 6.0), 40.0)
       ELSE m.meld_int
  END AS meld_na
FROM `uludag-mimic-insulin.insulin.cohort_final` c
LEFT JOIN icd i USING (hadm_id)
LEFT JOIN he_rx h USING (hadm_id)
LEFT JOIN labs l USING (stay_id)
LEFT JOIN meld_raw m USING (stay_id);

-- =====================================================================
-- 4C STEROID -- PATCH-10: prescriptions tabanli
--   00.5 BOS DONDU: MIMIC-IV ICU modulunde (inputevents/d_items) steroid YOK.
--   v1.0'daki itemid'ler UYDURMA idi (221319; 223258 gercekte INSULIN).
--   Neden emar degil: glukokortikoid etkisi saatler-gunler olceginde ve duz;
--   order cozunurlugu YETERLI. SINIRLILIK: order != uygulama (beyan edilir).
-- =====================================================================
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.steroid_orders` AS
SELECT
  p.hadm_id, p.starttime, p.stoptime,
  p.drug, p.route, sp.generic, sp.potency_pred_eq,
  p.dose_val_rx, p.dose_unit_rx, p.doses_per_24_hrs,
  SAFE_CAST(REGEXP_EXTRACT(p.dose_val_rx, r'^([0-9]+\.?[0-9]*)') AS FLOAT64)
    * COALESCE(p.doses_per_24_hrs, 1)
    * sp.potency_pred_eq AS pred_eq_per_day
FROM `physionet-data.mimiciv_3_1_hosp.prescriptions` p
JOIN `uludag-mimic-insulin.ref.steroid_potency` sp
  ON REGEXP_CONTAINS(LOWER(TRIM(p.drug)), sp.drug_regex)
JOIN `uludag-mimic-insulin.ref.steroid_routes` sr
  ON sr.route = p.route AND sr.systemic = TRUE       -- ROTA BEYAZ LISTESI
WHERE sp.verified = TRUE
  AND p.dose_unit_rx = 'mg'
  -- INVESTIGATIONAL DISLA: korlenmis calismada PLASEBO olabilir
  AND NOT REGEXP_CONTAINS(LOWER(p.drug), r'^inv-')
  -- DEPO IM formlari: yavas salinim, order penceresi anlamsiz
  AND NOT REGEXP_CONTAINS(LOWER(p.drug), r'acetate')
  AND p.starttime IS NOT NULL AND p.stoptime IS NOT NULL
  AND p.stoptime > p.starttime
  -- KK-7: ILAC-SPESIFIK TEK-DOZ TAVANI (dose_val_rx aykirilari).
  --   00.5b/04_kk7c: metilprednizolon dose_val_rx kuyrugu fizyolojik-disi
  --   (8695/8600/8050 mg TEK doz; cogunda doses_per_24_hrs=NULL). Puls tedavi
  --   tavani (metilpred 1000, deksa 40, pred/prednizolon 100, HC/kortizon 500,
  --   fludrokortizon 1 mg) asan ORDER dislanir. Bu bir VERI-KALITE filtresidir
  --   (registry dislama kriteri DEGIL; hasta duzeyinde kohort etkilemez).
  --   dose_val_rx parse edilemiyorsa (NULL) satir KORUNUR: pred_eq zaten NULL,
  --   yalniz on_steroid=1 katkisi kalir; mevcut davranis degismez.
  AND (
    SAFE_CAST(REGEXP_EXTRACT(p.dose_val_rx, r'^([0-9]+\.?[0-9]*)') AS FLOAT64) IS NULL
    OR SAFE_CAST(REGEXP_EXTRACT(p.dose_val_rx, r'^([0-9]+\.?[0-9]*)') AS FLOAT64) <=
       CASE sp.generic
         WHEN 'methylprednisolone' THEN 1000
         WHEN 'dexamethasone'      THEN 40
         WHEN 'prednisone'         THEN 100
         WHEN 'prednisolone'       THEN 100
         WHEN 'hydrocortisone'     THEN 500
         WHEN 'cortisone'          THEN 500
         WHEN 'fludrocortisone'    THEN 1
         ELSE 1e15
       END
  );

CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.steroid_hourly`
CLUSTER BY stay_id AS
SELECT
  g.stay_id, g.hr, g.ts,
  COALESCE(SUM(s.pred_eq_per_day), 0) AS pred_eq_24h,
  MAX(IF(s.hadm_id IS NOT NULL, 1, 0)) AS on_steroid
FROM `uludag-mimic-insulin.insulin.grid` g
JOIN `uludag-mimic-insulin.insulin.cohort_final` c USING (stay_id)
LEFT JOIN `uludag-mimic-insulin.insulin.steroid_orders` s
  ON s.hadm_id = c.hadm_id
 AND g.ts >= s.starttime AND g.ts < s.stoptime  -- her ikisi de DATETIME (grid.ts + hosp.prescriptions)
GROUP BY 1,2,3;

-- =====================================================================
-- 4D VAZOPRESOR + AGIRLIK + steroid birlestir
-- =====================================================================
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.other_hourly`
CLUSTER BY stay_id AS
WITH vaso AS (
  SELECT g.stay_id, g.hr,
    MAX(IF(v.norepinephrine_equivalent_dose > 0, 1, 0)) AS on_vasopressor,
    MAX(v.norepinephrine_equivalent_dose) AS ne_eq
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.norepinephrine_equivalent_dose` v
    ON v.stay_id = g.stay_id
   AND v.starttime <= g.ts AND v.endtime > g.ts
  GROUP BY 1,2
),
wt AS (
  SELECT g.stay_id, g.hr,
    ARRAY_AGG(w.weight ORDER BY w.starttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS weight_kg
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.weight_durations` w
    ON w.stay_id = g.stay_id AND w.starttime <= g.ts
  GROUP BY 1,2
)
SELECT g.stay_id, g.hr, g.ts,
       s.pred_eq_24h, s.on_steroid,
       v.on_vasopressor, v.ne_eq, wt.weight_kg
FROM `uludag-mimic-insulin.insulin.grid` g
LEFT JOIN `uludag-mimic-insulin.insulin.steroid_hourly` s USING (stay_id, hr)
LEFT JOIN vaso v USING (stay_id, hr)
LEFT JOIN wt      USING (stay_id, hr);

-- =====================================================================
-- KALITE KONTROLLERI
-- =====================================================================
-- KK-7 (PATCH-10): steroid yakalama. 00.5b beklentisi:
--   methylprednisolone ~12.9k, dexamethasone ~17k, prednisone ~10.8k,
--   hydrocortisone ~9.6k hadm
-- SELECT generic, COUNT(DISTINCT hadm_id) n_hadm,
--        APPROX_QUANTILES(pred_eq_per_day, 4) AS pred_eq_q
-- FROM `uludag-mimic-insulin.insulin.steroid_orders` GROUP BY 1 ORDER BY n_hadm DESC;
--   SAGLAMA (DUZELTILDI): pred_eq_per_day PRED-ESDEGERI skalasindadir; eski
--     "medyan 20-60 mg" esigi HATALIYDI. Puls/yuksek-doz rejimlerde medyan
--     50-100 mg pred-eq NORMALDIR (4 mg deksa x6.67 ~=27; metilpred 80 mg
--     x1.25 =100; puls 500-1000 mg metilpred cok daha yuksek). Beklenen
--     medyanlar: prednizon ~20, deksa ~40-80, metilpred ~50-100, HC ~25.
--     ASIL sorun dose_val_rx KUYRUGU (tek-doz aykirilari) idi -> yukaridaki
--     tek-doz tavani ile filtrelendi (bkz. qc/04_steroid_outliers.txt).
--     Cok kucuk medyan (<5 mg) hala lokal/oftalmik kalem sizintisina isaret eder.

-- KK-7b: rota beyaz listesi calisti mi? Sadece systemic=TRUE cikmali
-- SELECT route, COUNT(*) n FROM `uludag-mimic-insulin.insulin.steroid_orders`
-- GROUP BY 1 ORDER BY n DESC;

-- KK-10 (PATCH-02): Child-Pugh non-sirotikte NULL olmali
-- SELECT cirrhosis, COUNT(*) n, COUNTIF(child_pugh_score IS NOT NULL) AS n_cp
-- FROM `uludag-mimic-insulin.insulin.hepatic_static` GROUP BY 1;
--   cirrhosis=false -> n_cp = 0 OLMAK ZORUNDA

-- KK-11 (BULGU-01 + PATCH-04): MELD makul mu?
-- SELECT cirrhosis, meld_cr_capped_rrt, COUNT(*) n,
--        APPROX_QUANTILES(meld_na, 4) AS meld_q
-- FROM `uludag-mimic-insulin.insulin.hepatic_static` GROUP BY 1,2;
--   Medyan 10-25 bekleniyor. Hepsi 54-58 ise BULGU-01 uygulanmamis.

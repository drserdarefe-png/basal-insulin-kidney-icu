-- =====================================================================
-- 05 — OZELLIK TABLOSU + HEDEF DEGISKENLER                  v2.0
-- Uygulanmis: PATCH-09 (glukoz zamansal leakage: asimetrik zaman damgasi)
-- Kural 4 : Hedef DOZ DEGIL -> doz verilince GLUKOZ NE OLUR
-- Kural 3 : split SUBJECT duzeyinde
-- Kural 17: OZELLIK = kullanilabilirlik ani | HEDEF = ornekleme ani
-- =====================================================================

-- ---------------------------------------------------------------------
-- 5.1 GLUKOZ -- PATCH-09
--   00.4b: 220621 (serum), labevents kopyasi. 629.379 kayittan 606.337'si
--          (%96.3) +-30dk icinde eslesiyor, 563.808'i (%89.6) BIREBIR AYNI
--          deger. Eslesme charttime uzerinden -> sonuc kan alinma saatine
--          GERI TARIHLENIYOR.
--   00.4c: gercek kullanilabilirlik gecikmesi (storetime - charttime):
--          POC (225664)        medyan 13 dk  (D10 -8, D90 57; %21 negatif
--                              -> charttime saate yuvarlaniyor)
--          SERUM (220621)      medyan 66 dk  (D10 46, D90 107; %0.18 negatif)
--          WHOLE_BLOOD (226537) medyan 4 dk  (D10 2, D90 21)
--   -> Eski kod charttime ile forward-fill ediyordu: model, karar aninda
--      HENUZ VAR OLMAYAN serum glukozunu goruyordu. Tahmin ufku 2 saat,
--      sizan bilgi 1 saat.
--   GREATEST (COALESCE degil): n_null_store=0 -> COALESCE gereksiz. Ve
--      storetime'i dogrudan kullanmak POC'ta TERS leakage yaratir (%21'inde
--      deger olculmeden once bilinmis sayilirdi). GREATEST tutucu.
--   CIFT SAYIM YOK: labevents (50809/50931) EKLENMEZ; 220621 zaten kopyasi.
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.glu`
CLUSTER BY stay_id AS
WITH raw AS (
  SELECT ce.stay_id,
         ce.charttime,                      -- HEDEF: fizyolojik ornekleme ani
         -- OZELLIK: kullanilabilirlik ani. |lag|>24h -> veri hatasi, charttime'a dus
         CASE WHEN ABS(TIMESTAMP_DIFF(ce.storetime, ce.charttime, MINUTE)) > 1440
              THEN ce.charttime
              ELSE GREATEST(ce.charttime, ce.storetime)
         END AS available_time,
         IF(ABS(TIMESTAMP_DIFF(ce.storetime, ce.charttime, MINUTE)) > 1440, 1, 0)
              AS store_lag_outlier,
         ce.valuenum AS glu,
         CASE ce.itemid WHEN 225664 THEN 'POC'
                        WHEN 220621 THEN 'SERUM'
                        WHEN 226537 THEN 'WHOLE_BLOOD' END AS glu_src
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  WHERE ce.itemid IN (225664, 220621, 226537)
    AND ce.valuenum BETWEEN 10 AND 2000     -- sentinel: 999999, 1653550, -124
),
grid_glu AS (
  SELECT g.stay_id, g.hr, g.ts,
    ARRAY_AGG(STRUCT(r.glu, r.glu_src, r.store_lag_outlier)
              ORDER BY r.available_time DESC LIMIT 1)[SAFE_OFFSET(0)] AS last_meas,
    MIN(TIMESTAMP_DIFF(g.ts, r.available_time, MINUTE)) AS glu_age_min
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN raw r
    ON r.stay_id = g.stay_id
   AND r.available_time <= g.ts                        -- <<< PATCH-09 CEKIRDEGI
   AND r.available_time > TIMESTAMP_SUB(g.ts, INTERVAL 4 HOUR)   -- FF siniri
  GROUP BY 1,2,3
),
-- fwd / fwd4: charttime KULLANILIR (PATCH-09 asimetrisi).
-- Hedef = hastanin t+2h'teki fizyolojik glukozu; ne zaman ogrenildigi ONEMSIZ.
fwd AS (
  SELECT g.stay_id, g.hr,
    ARRAY_AGG(r.glu ORDER BY ABS(TIMESTAMP_DIFF(
        r.charttime, TIMESTAMP_ADD(g.ts, INTERVAL 2 HOUR), MINUTE))
      LIMIT 1)[SAFE_OFFSET(0)] AS glu_t2h
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN raw r
    ON r.stay_id = g.stay_id
   AND r.charttime BETWEEN TIMESTAMP_ADD(g.ts, INTERVAL 90 MINUTE)
                       AND TIMESTAMP_ADD(g.ts, INTERVAL 150 MINUTE)
  GROUP BY 1,2
),
fwd4 AS (
  SELECT g.stay_id, g.hr,
    MIN(r.glu)   AS glu_nadir_0_4h,
    COUNT(r.glu) AS n_meas_0_4h
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN raw r
    ON r.stay_id = g.stay_id
   AND r.charttime > g.ts
   AND r.charttime <= TIMESTAMP_ADD(g.ts, INTERVAL 4 HOUR)
  GROUP BY 1,2
)
SELECT
  gg.stay_id, gg.hr, gg.ts,
  gg.last_meas.glu               AS glu_now,
  gg.last_meas.glu_src           AS glu_src,
  gg.last_meas.store_lag_outlier AS glu_store_lag_outlier,
  gg.glu_age_min,
  f.glu_t2h,
  f4.glu_nadir_0_4h, f4.n_meas_0_4h,
  -- HEDEF 2: guvenlik
  IF(f4.glu_nadir_0_4h < 70, 1, 0) AS hypo_l1_0_4h,
  IF(f4.glu_nadir_0_4h < 54, 1, 0) AS hypo_l2_0_4h,
  -- Kural 8: olcum sikligi = klinisyen endisesi proxy'si -> OZELLIK DEGIL.
  -- Yalniz duyarlilik analizinde raporlanir.
  LAG(gg.glu_age_min, 1) OVER (PARTITION BY gg.stay_id ORDER BY gg.hr) AS glu_age_prev
FROM grid_glu gg
LEFT JOIN fwd  f  USING (stay_id, hr)
LEFT JOIN fwd4 f4 USING (stay_id, hr);

-- ---------------------------------------------------------------------
-- 5.2 BIRLESTIRME
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.features`
CLUSTER BY stay_id AS
SELECT
  g.subject_id, g.stay_id, g.hr, g.ts,

  -- ===== MARUZIYET =====
  i.dose_u_this_hr, i.basal_u_this_hr, i.iob_bolus_u, i.iob_basal_u,
  i.dose_24h_cum, g.hr - i.last_dose_hr AS hrs_since_dose,
  SAFE_DIVIDE(i.dose_u_this_hr, o.weight_kg) AS dose_u_per_kg,

  -- ===== GLUKOZ DURUMU =====
  gl.glu_now, gl.glu_src, gl.glu_age_min, gl.glu_store_lag_outlier,
  gl.glu_now - LAG(gl.glu_now, 1) OVER w AS d_glu_1h,
  gl.glu_now - LAG(gl.glu_now, 3) OVER w AS d_glu_3h,
  AVG(gl.glu_now) OVER (PARTITION BY g.stay_id ORDER BY g.hr
        ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS glu_mean_24h,
  STDDEV(gl.glu_now) OVER (PARTITION BY g.stay_id ORDER BY g.hr
        ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS glu_sd_24h,

  -- ===== NUTRISYON (effect modifier) =====
  ch.cho_g_hr, ch.cho_en, ch.cho_pn, ch.cho_maint, ch.cho_citrate,
  ch.cho_diluent_abx, ch.n_diluent_abx,
  ch.cho_d50_rescue, ch.n_d50_rescue,          -- KURAL 7: cho_g_hr'da YOK
  ch.d_cho_2h, ch.d_cho_4h, ch.cho_24h_cum,
  ch.cho_interrupted_4h, ch.fasting_hr, ch.nutrition_route,
  ch.tpn_unknown_flag, ch.en_diluted_flag,
  SAFE_DIVIDE(ch.cho_g_hr, o.weight_kg) AS cho_g_hr_per_kg,

  -- ===== RENAL (effect modifier) =====
  rs.egfr_baseline, rs.ckd_stage_baseline, rs.ckd_confirmed, rs.esrd,
  rs.ckd_icd_only_discordant, rs.ckd_lab_only_discordant,
  rh.egfr_now, rh.egfr_unsteady_flag, rh.cr_last, rh.cr_age_hr,
  rh.kdigo_stage, rh.on_crrt, rh.bun,

  -- ===== HEPATIK (effect modifier) =====
  hs.cirrhosis, hs.portal_htn, hs.ascites, hs.hepatic_enceph,
  hs.child_pugh_score, hs.meld_na, hs.meld_cr_capped_rrt,
  hs.bilirubin, hs.albumin, hs.inr,

  -- ===== DIGER =====
  hs.dm_any, hs.dm_t1, o.pred_eq_24h, o.on_steroid,
  o.on_vasopressor, o.ne_eq, o.weight_kg,
  bc.age, bc.gender, bc.cr_baseline_method,
  cf.cov_prior_pump, cf.sens_suspect_undoc, cf.pump_class,

  -- ===== ON-TANIMLI ETKILESIMLER (registry INT-01..05) =====
  IF(hs.cirrhosis AND ch.cho_interrupted_4h = 1, 1, 0)  AS int01_cirr_x_interrupt,
  i.dose_u_this_hr * IF(rh.egfr_now < 15, 1, 0)         AS int02_dose_x_egfr15,
  i.dose_u_this_hr * SAFE_DIVIDE(rh.bun, 20.0)          AS int03_dose_x_bun,
  o.pred_eq_24h * ch.cho_g_hr                           AS int04_steroid_x_cho,
  IF(gl.glu_src = 'POC' AND o.on_vasopressor = 1, 1, 0) AS int05_poc_x_vaso,

  -- ===== HEDEFLER =====
  gl.glu_t2h,
  gl.glu_t2h - gl.glu_now AS d_glu_t2h,        -- ana hedef (Asama 1 ISF)
  gl.glu_nadir_0_4h, gl.hypo_l1_0_4h, gl.hypo_l2_0_4h, gl.n_meas_0_4h,

  -- ===== SPLIT -- SUBJECT DUZEYINDE (Kural 3) =====
  MOD(ABS(FARM_FINGERPRINT(CAST(g.subject_id AS STRING))), 100) AS split_bucket

FROM `uludag-mimic-insulin.insulin.grid` g
LEFT JOIN `uludag-mimic-insulin.insulin.ins`  i  USING (stay_id, hr)
LEFT JOIN `uludag-mimic-insulin.insulin.glu`  gl USING (stay_id, hr)
LEFT JOIN `uludag-mimic-insulin.insulin.cho`  ch USING (stay_id, hr)
LEFT JOIN `uludag-mimic-insulin.insulin.renal_hourly` rh USING (stay_id, hr)
LEFT JOIN `uludag-mimic-insulin.insulin.other_hourly` o  USING (stay_id, hr)
LEFT JOIN `uludag-mimic-insulin.insulin.renal_static`   rs USING (stay_id)
LEFT JOIN `uludag-mimic-insulin.insulin.hepatic_static` hs USING (stay_id)
LEFT JOIN `uludag-mimic-insulin.insulin.baseline_cr`    bc USING (stay_id)
LEFT JOIN `uludag-mimic-insulin.insulin.cohort_final`   cf USING (stay_id)
WINDOW w AS (PARTITION BY g.stay_id ORDER BY g.hr);

-- ---------------------------------------------------------------------
-- 5.3 MODEL SETI
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.model_set` AS
SELECT *,
  CASE WHEN split_bucket < 60 THEN 'train'
       WHEN split_bucket < 80 THEN 'valid'
       ELSE 'test' END AS split
FROM `uludag-mimic-insulin.insulin.features`
WHERE glu_now IS NOT NULL
  AND glu_t2h IS NOT NULL
  AND glu_age_min <= 240;

-- =====================================================================
-- KALITE KONTROLLERI
-- =====================================================================
-- KK-LEAKAGE (Kural 3): SONUC 0 OLMAK ZORUNDA. Degilse DUR.
-- SELECT COUNT(*) AS leaked_subjects FROM (
--   SELECT subject_id FROM `uludag-mimic-insulin.insulin.model_set`
--   GROUP BY 1 HAVING COUNT(DISTINCT split) > 1
-- );

-- KK-5 (PATCH-09): SERUM olcumlerinin YASI, POC'tan ~1 saat BUYUK olmali.
--   Fark yoksa duzeltme uygulanmamis demektir.
-- SELECT glu_src, COUNT(*) n,
--        APPROX_QUANTILES(glu_age_min, 4) AS age_min_q,
--        AVG(glu_store_lag_outlier) AS pct_lag_outlier
-- FROM `uludag-mimic-insulin.insulin.glu` WHERE glu_now IS NOT NULL GROUP BY 1;

-- KK-12: model_set boyutu ve hedef dagilimi
-- SELECT split, COUNT(*) n, COUNT(DISTINCT subject_id) n_subj,
--        AVG(hypo_l1_0_4h) AS hypo_rate,
--        APPROX_QUANTILES(d_glu_t2h, 4) AS d_glu_q,
--        COUNTIF(dose_u_this_hr > 0) AS n_dosed
-- FROM `uludag-mimic-insulin.insulin.model_set` GROUP BY 1;
--   hypo_rate ~%3-8 bekleniyor. %0 ise fwd4 bozuk.

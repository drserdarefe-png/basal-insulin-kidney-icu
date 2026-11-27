-- =====================================================================
-- 07 — OLAY DUZEYINDE BAZAL VERI SETI (Makale 2)                v1.0
-- GRID BAGIMLILIGI YOK: kovaryatlar/hedefler OLAY ZAMANINA (dose_time)
--   yoneltildi. 14-gun grid kirpmasi UYGULANMAZ (qc/06_basal_unmatched.json:
--   kirpma 2.804 olay/343 hasta + YUKSEK-DOZ kuyrugu -- 20 U vs 16 U -- kaybettiriyordu).
-- Zaman: tum kaynaklar DATETIME (Kural 19: acik cast gerekmez; DATETIME_* kullan).
-- Kural 17: OZELLIK glu -> GREATEST(charttime,storetime) | HEDEF glu -> charttime.
-- Kural 8: n_meas_* OZELLIK DEGIL -> yalniz duyarlilik.
-- =====================================================================

-- ADIM 1: olay indeksi (kararli event_id). Glargine BIRINCIL, NPH ikincil.
--   pre-admission (event_hr<0, 4 olay) DISLA.
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.basal_event_index` AS
SELECT *,
  -- ayni yatis + ayni tip (glarjin analizinde = onceki glarjin) onceki doz
  LAG(units_total) OVER w AS prev_units_total,
  DATETIME_DIFF(dose_time, LAG(dose_time) OVER w, MINUTE) / 60.0 AS hrs_since_prev_glargine
FROM (
  SELECT
    GENERATE_UUID() AS event_id,
    e.stay_id, c.subject_id, c.hadm_id, c.intime,
    e.starttime AS dose_time, e.units_total,
    CASE WHEN e.label = 'Insulin - Glargine' THEN 'glargine' ELSE 'nph' END AS basal_type,
    DATETIME_DIFF(e.starttime, c.intime, HOUR) AS event_hr,
    c.los AS los_days
  FROM `uludag-mimic-insulin.insulin.ins_events` e
  JOIN `uludag-mimic-insulin.insulin.cohort_final` c USING (stay_id)
  WHERE e.is_basal = TRUE
    AND e.label IN ('Insulin - Glargine', 'Insulin - NPH')
    AND DATETIME_DIFF(e.starttime, c.intime, HOUR) >= 0
)
WINDOW w AS (PARTITION BY stay_id, basal_type ORDER BY dose_time);

-- ADIM 2: olay-zamanli kovaryatlar + hedefler + split
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.basal_events` AS
WITH idx AS (SELECT * FROM `uludag-mimic-insulin.insulin.basal_event_index`),
-- ===== glukoz ham (Kural 17 available_time) -- kohort filtreli (Kural 11) =====
glu_raw AS (
  SELECT ce.stay_id, ce.charttime, ce.valuenum AS glu,
    CASE WHEN ABS(DATETIME_DIFF(ce.storetime, ce.charttime, MINUTE)) > 1440 THEN ce.charttime
         ELSE GREATEST(ce.charttime, ce.storetime) END AS available_time
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  JOIN `uludag-mimic-insulin.insulin.cohort_final` c USING (stay_id)
  WHERE ce.itemid IN (225664, 220621, 226537) AND ce.valuenum BETWEEN 10 AND 2000
),
-- ===== RENAL @ dose_time =====
cr AS (
  SELECT bi.event_id,
    ARRAY_AGG(le.valuenum ORDER BY le.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS cr_last,
    MIN(DATETIME_DIFF(bi.dose_time, le.charttime, HOUR)) AS cr_age_hr
  FROM idx bi
  LEFT JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.subject_id = bi.subject_id AND le.itemid = 50912 AND le.valuenum BETWEEN 0.1 AND 20
   AND le.charttime <= bi.dose_time AND le.charttime > DATETIME_SUB(bi.dose_time, INTERVAL 48 HOUR)
  GROUP BY 1
),
bun AS (
  SELECT bi.event_id,
    ARRAY_AGG(le.valuenum ORDER BY le.charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS bun
  FROM idx bi
  LEFT JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le
    ON le.subject_id = bi.subject_id AND le.itemid = 51006
   AND le.charttime <= bi.dose_time AND le.charttime > DATETIME_SUB(bi.dose_time, INTERVAL 48 HOUR)
  GROUP BY 1
),
kdigo AS (
  SELECT bi.event_id, MAX(k.aki_stage) AS kdigo_stage
  FROM idx bi
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.kdigo_stages` k
    ON k.stay_id = bi.stay_id
   AND k.charttime <= bi.dose_time AND k.charttime > DATETIME_SUB(bi.dose_time, INTERVAL 6 HOUR)
  GROUP BY 1
),
crrt AS (
  SELECT bi.event_id, MAX(IF(cr.stay_id IS NOT NULL, 1, 0)) AS on_crrt
  FROM idx bi
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.crrt` cr
    ON cr.stay_id = bi.stay_id
   AND cr.charttime <= bi.dose_time AND cr.charttime > DATETIME_SUB(bi.dose_time, INTERVAL 1 HOUR)
  GROUP BY 1
),
-- ===== OTHER @ dose_time =====
wt AS (
  SELECT bi.event_id,
    ARRAY_AGG(w.weight ORDER BY w.starttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS weight_kg
  FROM idx bi
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.weight_durations` w
    ON w.stay_id = bi.stay_id AND w.starttime <= bi.dose_time
  GROUP BY 1
),
-- STEROID: generic + zamanlama-duyarli (pred_eq_24h yerine; theta_steroid~1 fix,
--   bkz. qc/07_steroid_timing.json). dur_hr: deksa 48, metilpred/pred 24, HC 10.
--   night_active = etki gece penceresine (doz+4h) uzaniyor mu (hrs_since <= dur-4).
steroid AS (
  SELECT bi.event_id,
    COALESCE(SUM(s.pred_eq_per_day), 0) AS pred_eq_active,
    MAX(IF(s.hadm_id IS NOT NULL, 1, 0)) AS on_steroid,
    MIN(DATETIME_DIFF(bi.dose_time, s.starttime, HOUR)) AS hrs_since_steroid_start,
    MAX(IF(DATETIME_DIFF(bi.dose_time, s.starttime, HOUR) <=
        (CASE s.generic WHEN 'dexamethasone' THEN 48 WHEN 'methylprednisolone' THEN 24
              WHEN 'prednisone' THEN 24 WHEN 'prednisolone' THEN 24
              WHEN 'hydrocortisone' THEN 10 WHEN 'cortisone' THEN 10
              WHEN 'fludrocortisone' THEN 24 ELSE 24 END) - 4, 1, 0)) AS steroid_night_active,
    ARRAY_AGG(s.generic IGNORE NULLS ORDER BY s.pred_eq_per_day DESC LIMIT 1)[SAFE_OFFSET(0)]
      AS steroid_generic
  FROM idx bi
  LEFT JOIN `uludag-mimic-insulin.insulin.steroid_orders` s
    ON s.hadm_id = bi.hadm_id
   AND bi.dose_time >= s.starttime AND bi.dose_time < s.stoptime
  GROUP BY 1
),
vaso AS (
  SELECT bi.event_id, MAX(IF(v.norepinephrine_equivalent_dose > 0, 1, 0)) AS on_vasopressor
  FROM idx bi
  LEFT JOIN `physionet-data.mimiciv_3_1_derived.norepinephrine_equivalent_dose` v
    ON v.stay_id = bi.stay_id AND v.starttime <= bi.dose_time AND v.endtime > bi.dose_time
  GROUP BY 1
),
-- ===== GLUKOZ DURUMU @ dose_time =====
glu_now AS (
  SELECT bi.event_id,
    ARRAY_AGG(r.glu ORDER BY r.available_time DESC LIMIT 1)[SAFE_OFFSET(0)] AS glu_now
  FROM idx bi
  LEFT JOIN glu_raw r ON r.stay_id = bi.stay_id
   AND r.available_time <= bi.dose_time
   AND r.available_time > DATETIME_SUB(bi.dose_time, INTERVAL 24 HOUR)
  GROUP BY 1
),
glu_stat AS (
  SELECT bi.event_id, AVG(r.glu) AS glu_mean_24h, STDDEV(r.glu) AS glu_sd_24h,
         MIN(r.glu) AS glu_min_24h
  FROM idx bi
  LEFT JOIN glu_raw r ON r.stay_id = bi.stay_id
   AND r.charttime <= bi.dose_time AND r.charttime > DATETIME_SUB(bi.dose_time, INTERVAL 24 HOUR)
  GROUP BY 1
),
-- ===== BOLUS insulin @ dose_time (bolus/bazal ayrimi) =====
bolus AS (SELECT stay_id, starttime, units_total, tau_min
          FROM `uludag-mimic-insulin.insulin.ins_events` WHERE is_basal = FALSE),
dose_bolus AS (
  SELECT bi.event_id,
    COALESCE(SUM(b.units_total), 0) AS dose_24h_cum,
    COALESCE(SUM(b.units_total *
      EXP(-DATETIME_DIFF(bi.dose_time, b.starttime, MINUTE) / NULLIF(b.tau_min, 0))), 0) AS iob_bolus_u
  FROM idx bi
  LEFT JOIN bolus b ON b.stay_id = bi.stay_id
   AND b.starttime <= bi.dose_time AND b.starttime > DATETIME_SUB(bi.dose_time, INTERVAL 24 HOUR)
  GROUP BY 1
),
-- ===== CHO @ dose_time (grid-saati; horizon otesi NULL + bayrak) =====
cho_ev AS (
  SELECT bi.event_id, ch.cho_g_hr, ch.cho_24h_cum,
    IF(ch.stay_id IS NOT NULL, 1, 0) AS cho_from_grid
  FROM idx bi
  LEFT JOIN `uludag-mimic-insulin.insulin.cho` ch
    ON ch.stay_id = bi.stay_id
   AND bi.dose_time >= ch.ts AND bi.dose_time < DATETIME_ADD(ch.ts, INTERVAL 1 HOUR)
),
-- ===== HEDEFLER (charttime; Kural 17) doz+0..24h ve gece doz+4..12h =====
tgt AS (
  SELECT bi.event_id,
    MIN(r.glu) AS glu_nadir_0_24h,
    COUNT(r.glu) AS n_meas_0_24h,
    SAFE_DIVIDE(COUNTIF(r.glu BETWEEN 140 AND 180), COUNT(r.glu)) AS tir_0_24h,
    MIN(IF(r.charttime >  DATETIME_ADD(bi.dose_time, INTERVAL 4 HOUR)
       AND r.charttime <= DATETIME_ADD(bi.dose_time, INTERVAL 12 HOUR), r.glu, NULL)) AS glu_nadir_night,
    COUNTIF(r.charttime >  DATETIME_ADD(bi.dose_time, INTERVAL 4 HOUR)
        AND r.charttime <= DATETIME_ADD(bi.dose_time, INTERVAL 12 HOUR)) AS n_meas_night
  FROM idx bi
  LEFT JOIN glu_raw r ON r.stay_id = bi.stay_id
   AND r.charttime >  bi.dose_time
   AND r.charttime <= DATETIME_ADD(bi.dose_time, INTERVAL 24 HOUR)
  GROUP BY 1
)
SELECT
  bi.event_id, bi.subject_id, bi.stay_id, bi.hadm_id, bi.dose_time,
  bi.basal_type, bi.units_total, bi.event_hr, bi.los_days,
  -- RENAL
  cr.cr_last, cr.cr_age_hr, bn.bun, kd.kdigo_stage, ct.on_crrt,
  CASE WHEN cr.cr_last IS NOT NULL THEN
    142 * POW(LEAST(cr.cr_last / IF(d.gender='F',0.7,0.9), 1), IF(d.gender='F',-0.241,-0.302))
        * POW(GREATEST(cr.cr_last / IF(d.gender='F',0.7,0.9), 1), -1.200)
        * POW(0.9938, d.age) * IF(d.gender='F', 1.012, 1.0)
  END AS egfr_now,
  -- OTHER (steroid: generic + zamanlama-duyarli)
  wt.weight_kg, st.on_steroid, st.pred_eq_active, st.hrs_since_steroid_start,
  st.steroid_night_active, st.steroid_generic, vs.on_vasopressor,
  -- HEPATIK (statik)
  hs.cirrhosis, hs.child_pugh_score, hs.meld_na, hs.dm_any, hs.dm_t1,
  -- NUTRISYON (grid-saati; bayrak)
  che.cho_g_hr, che.cho_24h_cum, che.cho_from_grid,
  -- GLUKOZ DURUMU
  gn.glu_now, gst.glu_mean_24h, gst.glu_sd_24h, gst.glu_min_24h,
  db.dose_24h_cum, db.iob_bolus_u,
  d.age, d.gender,
  -- DOZ GECMISI (ayni yatis onceki glarjin)
  SAFE_DIVIDE(bi.prev_units_total, wt.weight_kg) AS prev_glargine_dose_u_kg,
  bi.hrs_since_prev_glargine,
  -- doz U/kg (positivity/GPS icin)
  SAFE_DIVIDE(bi.units_total, wt.weight_kg) AS dose_u_kg,
  -- HEDEFLER
  tg.glu_nadir_0_24h,
  IF(tg.glu_nadir_0_24h < 70, 1, 0) AS hypo_24h_l1,
  IF(tg.glu_nadir_0_24h < 54, 1, 0) AS hypo_24h_l2,
  tg.tir_0_24h,
  tg.glu_nadir_night,                              -- BIRINCIL HEDEF
  IF(tg.glu_nadir_night < 70, 1, 0) AS hypo_night_l1,
  IF(tg.glu_nadir_night < 54, 1, 0) AS hypo_night_l2,
  tg.n_meas_0_24h,                                  -- Kural 8: duyarlilik
  tg.n_meas_night,                                  -- Kural 8: duyarlilik
  -- SPLIT (model_set ile AYNI bucket -- Kural 3)
  MOD(ABS(FARM_FINGERPRINT(CAST(bi.subject_id AS STRING))), 100) AS split_bucket,
  CASE WHEN MOD(ABS(FARM_FINGERPRINT(CAST(bi.subject_id AS STRING))), 100) < 60 THEN 'train'
       WHEN MOD(ABS(FARM_FINGERPRINT(CAST(bi.subject_id AS STRING))), 100) < 80 THEN 'valid'
       ELSE 'test' END AS split
FROM idx bi
LEFT JOIN cr        USING (event_id)
LEFT JOIN bun       bn USING (event_id)
LEFT JOIN kdigo     kd USING (event_id)
LEFT JOIN crrt      ct USING (event_id)
LEFT JOIN wt        USING (event_id)
LEFT JOIN steroid   st USING (event_id)
LEFT JOIN vaso      vs USING (event_id)
LEFT JOIN glu_now   gn USING (event_id)
LEFT JOIN glu_stat  gst USING (event_id)
LEFT JOIN dose_bolus db USING (event_id)
LEFT JOIN cho_ev    che USING (event_id)
LEFT JOIN tgt       tg USING (event_id)
LEFT JOIN `uludag-mimic-insulin.insulin.hepatic_static` hs USING (stay_id)
LEFT JOIN `uludag-mimic-insulin.insulin.baseline_cr`    d  USING (stay_id);

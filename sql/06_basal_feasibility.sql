-- =====================================================================
-- 06 — BAZAL INSULIN FIZIBILITE TANISI (kapsam/ikinci-asama karari icin)
-- Kaynak: insulin.ins_events (is_basal=TRUE) x cohort_final.
-- Zaman: grid.ts, ins_events.starttime HEPSI DATETIME (Kural 19: cast yok).
-- Ciktilar: qc/06_basal_feasibility.json, qc/06_basal_x_egfr.json
-- =====================================================================

-- ---------------------------------------------------------------------
-- SORGU 1: label bazli fizibilite (qc/06_basal_feasibility.json)
-- ---------------------------------------------------------------------
WITH ev AS (
  SELECT e.label, c.subject_id, e.stay_id, e.starttime, e.units_total
  FROM `uludag-mimic-insulin.insulin.ins_events` e
  JOIN `uludag-mimic-insulin.insulin.cohort_final` c USING (stay_id)
  WHERE e.is_basal = TRUE
),
per_pat AS (SELECT label, subject_id, COUNT(*) n_ev FROM ev GROUP BY 1, 2),
lab_perpat AS (SELECT label, APPROX_QUANTILES(n_ev, 4) events_per_pat_q FROM per_pat GROUP BY 1),
lab_hour AS (SELECT label, EXTRACT(HOUR FROM starttime) hod, COUNT(*) n FROM ev GROUP BY 1, 2),
lab_hour_arr AS (SELECT label, ARRAY_AGG(STRUCT(hod, n) ORDER BY hod) hour_hist FROM lab_hour GROUP BY 1),
lab_ev AS (
  SELECT label, COUNT(*) n_events, COUNT(DISTINCT stay_id) n_stays,
         COUNT(DISTINCT subject_id) n_subjects,
         APPROX_QUANTILES(units_total, 10) units_total_deciles
  FROM ev GROUP BY 1
)
SELECT le.label, le.n_events, le.n_stays, le.n_subjects,
       lp.events_per_pat_q, le.units_total_deciles, lha.hour_hist
FROM lab_ev le
JOIN lab_perpat lp USING (label)
JOIN lab_hour_arr lha USING (label)
ORDER BY le.n_events DESC;

-- ---------------------------------------------------------------------
-- SORGU 2: glarjin doz x eGFR (qc/06_basal_x_egfr.json)
--   doz saatindeki egfr_now (renal_hourly) + weight_kg (other_hourly)
-- ---------------------------------------------------------------------
WITH glar AS (
  SELECT e.stay_id, c.subject_id, e.starttime, e.units_total
  FROM `uludag-mimic-insulin.insulin.ins_events` e
  JOIN `uludag-mimic-insulin.insulin.cohort_final` c USING (stay_id)
  WHERE e.is_basal = TRUE AND e.label = 'Insulin - Glargine'
),
matched AS (
  SELECT g.stay_id, gl.subject_id, gl.units_total, rh.egfr_now, oh.weight_kg
  FROM glar gl
  JOIN `uludag-mimic-insulin.insulin.grid` g
    ON g.stay_id = gl.stay_id
   AND gl.starttime >= g.ts AND gl.starttime < DATETIME_ADD(g.ts, INTERVAL 1 HOUR)
  LEFT JOIN `uludag-mimic-insulin.insulin.renal_hourly` rh ON rh.stay_id = g.stay_id AND rh.hr = g.hr
  LEFT JOIN `uludag-mimic-insulin.insulin.other_hourly` oh ON oh.stay_id = g.stay_id AND oh.hr = g.hr
)
SELECT
  CASE WHEN egfr_now <  15 THEN '1_egfr<15'
       WHEN egfr_now <  45 THEN '2_egfr_15_45'
       WHEN egfr_now <  90 THEN '3_egfr_45_90'
       WHEN egfr_now >= 90 THEN '4_egfr>=90'
       ELSE '9_unknown' END AS egfr_band,
  COUNT(*) n_events, COUNT(DISTINCT subject_id) n_subjects,
  APPROX_QUANTILES(units_total, 10) units_total_deciles,
  APPROX_QUANTILES(SAFE_DIVIDE(units_total, weight_kg), 10) units_per_kg_deciles
FROM matched
GROUP BY 1 ORDER BY egfr_band;

-- ---------------------------------------------------------------------
-- SORGU 3: eslesme kaybi tanisi (qc/06_basal_unmatched.json)
--   matched / no_grid / no_egfr. event_hr = intime'dan bu yana saat (grid
--   horizon 336h=14gun testi). BULGU: no_grid'in %99.86'si event_hr>336
--   (uzun yatis, grid 14-gun kirpmasi); no_egfr yalniz 53 olay (ihmal).
-- ---------------------------------------------------------------------
-- WITH glar AS (
--   SELECT e.stay_id, e.units_total,
--          DATETIME_DIFF(e.starttime, c.intime, HOUR) AS event_hr,
--          DATETIME_DIFF(c.outtime, c.intime, HOUR)/24.0 AS los_days, e.starttime
--   FROM `uludag-mimic-insulin.insulin.ins_events` e
--   JOIN `uludag-mimic-insulin.insulin.cohort_final` c USING (stay_id)
--   WHERE e.is_basal = TRUE AND e.label = 'Insulin - Glargine'
-- ),
-- gmatch AS (
--   SELECT gl.stay_id, gl.units_total, gl.event_hr, gl.los_days,
--          g.hr AS grid_hr, rh.egfr_now
--   FROM glar gl
--   LEFT JOIN `uludag-mimic-insulin.insulin.grid` g ON g.stay_id = gl.stay_id
--     AND gl.starttime >= g.ts AND gl.starttime < DATETIME_ADD(g.ts, INTERVAL 1 HOUR)
--   LEFT JOIN `uludag-mimic-insulin.insulin.renal_hourly` rh
--     ON rh.stay_id = g.stay_id AND rh.hr = g.hr
-- )
-- SELECT CASE WHEN grid_hr IS NULL THEN 'no_grid'
--             WHEN egfr_now IS NULL THEN 'no_egfr' ELSE 'matched' END AS match_status,
--        COUNT(*) n_events, COUNT(DISTINCT stay_id) n_stays,
--        APPROX_QUANTILES(units_total,4) units_total_q,
--        APPROX_QUANTILES(event_hr,4) event_hr_q,
--        APPROX_QUANTILES(los_days,4) los_days_q
-- FROM gmatch GROUP BY 1 ORDER BY match_status;

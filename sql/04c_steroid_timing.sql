-- =====================================================================
-- 04c — STEROID ZAMANLAMA TANISI (Makale 2, gece hipoglisemi maruziyeti)
-- Baglam: bolus modelinde theta_steroid=1.00 [0.86-1.16] (sifir etki, klinik
--   imkansiz). pred_eq_24h kumulatif+generic-agnostik. Steroid etkisi ZAMANLAMAYA
--   ve GENERICE bagli (etki suresi: deksa 36-54h, metilpred/pred 18-36h, HC 8-12h).
-- Ciktilar: qc/07_steroid_timing.json
-- =====================================================================

-- SORGU 1: order baslama-saati histogrami (sabah piki?)
SELECT generic, EXTRACT(HOUR FROM starttime) AS hod, COUNT(*) n
FROM `uludag-mimic-insulin.insulin.steroid_orders`
GROUP BY 1,2 ORDER BY generic, hod;

-- SORGU 2: glarjin doz olaylarinda AKTIF steroid + gece penceresi kapsama
-- night_active = hrs_since_start <= dur-4 (etki gece penceresi basina uzaniyor mu)
WITH gl AS (
  SELECT event_id, hadm_id, dose_time
  FROM `uludag-mimic-insulin.insulin.basal_events` WHERE basal_type='glargine'
),
active AS (
  SELECT g.event_id, s.generic, s.pred_eq_per_day,
    DATETIME_DIFF(g.dose_time, s.starttime, HOUR) AS hrs_since_start,
    CASE s.generic
      WHEN 'dexamethasone' THEN 48 WHEN 'methylprednisolone' THEN 24
      WHEN 'prednisone' THEN 24 WHEN 'prednisolone' THEN 24
      WHEN 'hydrocortisone' THEN 10 WHEN 'cortisone' THEN 10
      WHEN 'fludrocortisone' THEN 24 ELSE 24 END AS dur_hr
  FROM gl g
  JOIN `uludag-mimic-insulin.insulin.steroid_orders` s
    ON s.hadm_id = g.hadm_id
   AND g.dose_time >= s.starttime AND g.dose_time < s.stoptime
)
SELECT generic, COUNT(*) n_active_orders, COUNT(DISTINCT event_id) n_glargine_events,
  APPROX_QUANTILES(hrs_since_start, 4) hrs_since_start_q,
  APPROX_QUANTILES(pred_eq_per_day, 4) pred_eq_active_q,
  ROUND(AVG(IF(hrs_since_start <= dur_hr - 4, 1, 0)), 3) AS frac_night_active
FROM active GROUP BY 1 ORDER BY n_glargine_events DESC;

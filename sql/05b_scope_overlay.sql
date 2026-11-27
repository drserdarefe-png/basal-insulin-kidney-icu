-- =====================================================================
-- 05b — KAPSAM KARARI TANISI (model_set uzerine NPO + infuzyon bindirmesi)
-- Amac: kapsam duyarliligi. model_set'i DEGISTIRMEZ; yeni tablo model_set_scope.
-- Kural 19 notu: grid.ts, poe.ordertime, ins_events.starttime/endtime HEPSI
--   DATETIME. Cross-modul TIMESTAMP() cast'i GEREKMEZ ve zararlidir
--   (bkz. KK-7 DATETIME fix, commit 87777fd). Bu boru hattinda icu tablolari
--   da DATETIME (grid, ins_events pipeline-turevi).
-- NPO aralik tavani 48h AMPIRIK: baglantili (discontinue-link) altkumenin
--   sure decileleri D50~17.6h, D90~28.5h (max 259h aykiri) -> 48h guvenli.
-- =====================================================================
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.model_set_scope` AS
WITH cohort_hadm AS (
  SELECT DISTINCT hadm_id FROM `uludag-mimic-insulin.insulin.cohort_final`
),
npo_raw AS (                       -- NPO order baslangiclari (D/C transaction'i DEGIL)
  SELECT p.hadm_id, p.poe_id, p.ordertime AS npo_start, p.discontinued_by_poe_id
  FROM `physionet-data.mimiciv_3_1_hosp.poe` p
  JOIN cohort_hadm ch USING (hadm_id)
  WHERE p.order_type = 'Nutrition'
    AND p.order_subtype = 'NPO/Diet for Procedure'
    AND p.transaction_type IN ('New', 'Change')
),
diet AS (                          -- Diet Order (sonraki diet zamani icin)
  SELECT p.hadm_id, p.ordertime AS diet_time
  FROM `physionet-data.mimiciv_3_1_hosp.poe` p
  JOIN cohort_hadm ch USING (hadm_id)
  WHERE p.order_type = 'Nutrition' AND p.order_subtype = 'Diet Order'
),
next_diet AS (
  SELECT nr.poe_id, MIN(dt.diet_time) AS next_diet_time
  FROM npo_raw nr
  JOIN diet dt ON dt.hadm_id = nr.hadm_id AND dt.diet_time > nr.npo_start
  GROUP BY 1
),
npo_int AS (                       -- aralik ucu: min(disc-link, sonraki diet, 48h tavani)
  SELECT n.hadm_id, n.npo_start,
    LEAST(
      DATETIME_ADD(n.npo_start, INTERVAL 48 HOUR),
      COALESCE(d.ordertime,       DATETIME_ADD(n.npo_start, INTERVAL 48 HOUR)),
      COALESCE(nd.next_diet_time, DATETIME_ADD(n.npo_start, INTERVAL 48 HOUR))
    ) AS npo_end
  FROM npo_raw n
  LEFT JOIN `physionet-data.mimiciv_3_1_hosp.poe` d
    ON d.poe_id = n.discontinued_by_poe_id AND d.ordertime > n.npo_start
  LEFT JOIN next_diet nd ON nd.poe_id = n.poe_id
),
npo_hourly AS (
  SELECT g.stay_id, g.hr, 1 AS npo_active
  FROM `uludag-mimic-insulin.insulin.grid` g
  JOIN npo_int n
    ON n.hadm_id = g.hadm_id
   AND g.ts >= n.npo_start AND g.ts < n.npo_end        -- DATETIME vs DATETIME
  GROUP BY 1, 2
),
inf_hourly AS (                    -- o saatte VEYA son 2 saatte aktif infuzyon
  SELECT g.stay_id, g.hr, 1 AS insulin_infusion_active
  FROM `uludag-mimic-insulin.insulin.grid` g
  JOIN `uludag-mimic-insulin.insulin.ins_events` e
    ON e.stay_id = g.stay_id
   AND e.route_mode = 'infusion'
   AND e.starttime <= g.ts
   AND e.endtime  > DATETIME_SUB(g.ts, INTERVAL 2 HOUR)
  GROUP BY 1, 2
)
SELECT m.*,
  COALESCE(nh.npo_active, 0)              AS npo_active,
  COALESCE(ih.insulin_infusion_active, 0) AS insulin_infusion_active
FROM `uludag-mimic-insulin.insulin.model_set` m
LEFT JOIN npo_hourly nh USING (stay_id, hr)
LEFT JOIN inf_hourly ih USING (stay_id, hr);

-- ---------------------------------------------------------------------
-- KAPSAM KOMPOZISYONU (qc/05_scope3.json)
-- bq query --format=prettyjson < (asagisi) > qc/05_scope3.json
-- ---------------------------------------------------------------------
-- SELECT nutrition_route, npo_active, insulin_infusion_active,
--        COUNT(*) n, COUNT(DISTINCT subject_id) n_subj,
--        COUNTIF(dose_u_this_hr > 0) n_dosed,
--        ROUND(AVG(hypo_l1_0_4h),4) hypo_rate, SUM(hypo_l1_0_4h) hypo_events,
--        APPROX_QUANTILES(d_glu_t2h, 4) d_glu_q
-- FROM `uludag-mimic-insulin.insulin.model_set_scope`
-- GROUP BY 1,2,3 ORDER BY n DESC;

-- =====================================================================
-- 03 — INSULIN OLAYLARI + IOB                               v2.0
-- Uygulanmis: PATCH-05 (ref.insulin_pk TEK KAYNAK; inline tablo SILINDI)
--             BULGU-06 (Regular icin rota-bagimli tau: IV inf 45, SC 150)
-- 00.1: Regular n=337.446, %47 infuzyon (rate NOT NULL) -> rota ayrimi kritik
-- Kural 14: rota SATIR duzeyinde
-- =====================================================================

CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.ins_events` AS
SELECT
  ie.stay_id, ie.itemid, pk.label, pk.is_basal,
  ie.starttime, ie.endtime,
  ie.amount AS units_total,
  ie.rate   AS units_per_hr,
  CASE WHEN ie.rate IS NOT NULL AND ie.rateuom LIKE '%hour%'
       THEN 'infusion' ELSE 'bolus' END AS route_mode,
  -- BULGU-06: Regular IV infuzyon farkli farmakokinetik
  CASE WHEN ie.rate IS NOT NULL AND ie.rateuom LIKE '%hour%'
            AND pk.tau_infusion_min IS NOT NULL
       THEN pk.tau_infusion_min
       ELSE pk.tau_min
  END AS tau_min,
  ie.ordercategorydescription
FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
JOIN `uludag-mimic-insulin.ref.insulin_pk` pk USING (itemid)
WHERE ie.statusdescription NOT IN ('Rewritten')
  AND ie.amount > 0
  AND ie.amountuom IN ('units','U')
  AND pk.verified = TRUE
  -- AYKIRI TAVANI (Kural 15, qc/03_dose_outliers.txt): fizyolojik-disi tek-olay.
  --   BASAL-SPESIFIK: basal (glarjin/NPH SC) tavani 150 U (~2 U/kg x 75 kg;
  --     60-80 U dirençli/obezde MESRU). Bolus (rapid/regular) 50 U. Infuzyon 50 U/h.
  --   50 U bolus-tavani basal glarjini TARAFLI kirpiyordu (yuksek doz=yuksek ihtiyac
  --     =dirençli/obez -> weight ayagini zayiflatir). Asil hedef Regular inf/bolus
  --     (max 3013/1600 U) 50'de kalir.
  AND NOT (
    CASE WHEN ie.rate IS NOT NULL AND ie.rateuom LIKE '%hour%' THEN ie.rate > 50
         WHEN pk.is_basal THEN ie.amount > 150
         ELSE ie.amount > 50
    END
  );

-- ---------------------------------------------------------------------
-- 3.1 IOB -- ustel bozunum, 12h geriye bakan pencere
--     Bolus  : tek noktada (starttime)
--     Infuzyon: saatlik dilimlere ayrilir, her dilim ayri bolus gibi
--     Bazal (glargine/NPH/karisimlar) AYRI HAVUZ:
--       glargine ~24h DUZ profil -> ustel model fizyolojik DEGIL
--       NPH PIKLI -> ustel yaklasim ZAYIF
--     Ikisini bolus havuzuna katmak ANA MARUZIYET DEGISKENINI bozar.
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.ins`
CLUSTER BY stay_id AS
WITH inf_slices AS (
  SELECT e.stay_id, e.tau_min, e.is_basal,
         TIMESTAMP_ADD(e.starttime, INTERVAL s HOUR) AS give_time,
         e.units_per_hr AS units
  FROM `uludag-mimic-insulin.insulin.ins_events` e,
  UNNEST(GENERATE_ARRAY(0,
    GREATEST(CAST(FLOOR(TIMESTAMP_DIFF(e.endtime, e.starttime, HOUR)) AS INT64) - 1, 0))) AS s
  WHERE e.route_mode = 'infusion'
),
bolus AS (
  SELECT stay_id, tau_min, is_basal, starttime AS give_time, units_total AS units
  FROM `uludag-mimic-insulin.insulin.ins_events`
  WHERE route_mode = 'bolus'
),
give AS (SELECT * FROM inf_slices UNION ALL SELECT * FROM bolus),

iob AS (
  SELECT
    g.stay_id, g.hr, g.ts,
    COALESCE(SUM(IF(NOT gv.is_basal,
      gv.units * EXP(-TIMESTAMP_DIFF(g.ts, gv.give_time, MINUTE) / gv.tau_min), 0)), 0)
      AS iob_bolus_u,
    COALESCE(SUM(IF(gv.is_basal,
      gv.units * EXP(-TIMESTAMP_DIFF(g.ts, gv.give_time, MINUTE) / gv.tau_min), 0)), 0)
      AS iob_basal_u
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN give gv
    ON gv.stay_id = g.stay_id
   AND gv.give_time <= g.ts
   AND gv.give_time >  TIMESTAMP_SUB(g.ts, INTERVAL 12 HOUR)
  GROUP BY 1,2,3
),
-- bu saatte fiilen verilen doz = MARUZIYET degiskeni
dose_now AS (
  SELECT g.stay_id, g.hr,
    COALESCE(SUM(IF(NOT gv.is_basal, gv.units, 0)), 0) AS dose_u_this_hr,
    COALESCE(SUM(IF(gv.is_basal,     gv.units, 0)), 0) AS basal_u_this_hr
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN give gv
    ON gv.stay_id = g.stay_id
   AND gv.give_time >= g.ts
   AND gv.give_time <  TIMESTAMP_ADD(g.ts, INTERVAL 1 HOUR)
  GROUP BY 1,2
)
SELECT
  i.*, d.dose_u_this_hr, d.basal_u_this_hr,
  SUM(d.dose_u_this_hr) OVER (PARTITION BY i.stay_id ORDER BY i.hr
      ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS dose_24h_cum,
  MAX(IF(d.dose_u_this_hr > 0, i.hr, NULL)) OVER (
      PARTITION BY i.stay_id ORDER BY i.hr ROWS UNBOUNDED PRECEDING) AS last_dose_hr
FROM iob i
JOIN dose_now d USING (stay_id, hr);

-- ---------------------------------------------------------------------
-- KK-8: IOB medyani doz medyaninin ~1.5-3 kati olmali (tau birikimi)
-- ---------------------------------------------------------------------
-- SELECT APPROX_QUANTILES(iob_bolus_u, 4) AS iob_q,
--        APPROX_QUANTILES(dose_u_this_hr, 4) AS dose_q,
--        APPROX_QUANTILES(iob_basal_u, 4) AS basal_q
-- FROM `uludag-mimic-insulin.insulin.ins` WHERE dose_u_this_hr > 0;

-- KK-9 (BULGU-06): rota-bagimli tau uygulanmis mi?
-- SELECT route_mode, label, COUNT(*) n, ANY_VALUE(tau_min) AS tau
-- FROM `uludag-mimic-insulin.insulin.ins_events` GROUP BY 1,2 ORDER BY n DESC;
-- Regular/infusion -> tau=45 | Regular/bolus -> tau=150 OLMALI

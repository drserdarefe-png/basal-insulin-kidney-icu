-- =====================================================================
-- 02 — KARBONHIDRAT AKISI (g/saat)                          v2.0
-- Uygulanmis: PATCH-01 (bolus saat kovasi), PATCH-08 (yeniden yazim:
--             rota SATIR duzeyinde, aykiri filtreleri, Dextrose PN gram)
-- Kural 14: rota itemid'in degil SATIRIN ozelligi (rate IS NULL)
-- Kural 7 : D50 kurtarma cho_g_hr TOPLAMINA GIRMEZ (sonuc, girdi degil)
-- Kural 2 : verified=FALSE olan ref satiri boru hattina GIRMEZ
-- =====================================================================

-- ---------------------------------------------------------------------
-- 2.1 SUREKLI INFUZYONLAR (rate NOT NULL)
--     mL/min: 5 olay, mL/kg/hour: 1 olay -> ihmal, dipnot (00.3c)
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.cho_infusion` AS
WITH en AS (
  SELECT ie.stay_id, ie.starttime, ie.endtime, 'EN' AS src,
         ie.rate * f.cho_g_per_l_eff / 1000.0 AS cho_g_hr,  -- _eff: sulandirma dahil
         ie.rate AS ml_hr,
         f.strength_frac
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN `uludag-mimic-insulin.ref.enteral_formula_cho` f ON f.itemid = ie.itemid
  WHERE ie.statusdescription NOT IN ('Rewritten')
    AND ie.rateuom = 'mL/hour'
    AND ie.rate > 0 AND ie.rate <= 250          -- aykiri: 21600/60000 mL/h var
    AND f.verified = TRUE
),
dex AS (
  SELECT ie.stay_id, ie.starttime, ie.endtime,
         CASE WHEN d.fluid_name LIKE 'ACD-A%' THEN 'CRRT_CITRATE' ELSE 'MAINT' END AS src,
         ie.rate * d.cho_g_per_l / 1000.0 AS cho_g_hr,
         ie.rate AS ml_hr,
         CAST(NULL AS FLOAT64) AS strength_frac
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN `uludag-mimic-insulin.ref.dextrose_sources` d ON d.itemid = ie.itemid
  WHERE ie.statusdescription NOT IN ('Rewritten')
    AND ie.rateuom = 'mL/hour'
    AND ie.rate > 0 AND ie.rate <= 250          -- aykiri: -118.5 ve 30000 var
    AND d.special_handling = 'rate_or_bolus'
    AND d.cho_g_per_l IS NOT NULL
    AND d.verified = TRUE
)
SELECT * FROM en UNION ALL SELECT * FROM dex;

-- ---------------------------------------------------------------------
-- 2.2 BOLUSLAR (rate IS NULL)  -- PATCH-08 CEKIRDEGI
--     220949 rate NULL: 453.573 antibiyotik sulandiricisi (49.463 hasta)
--                      + 1.998 sivi bolusu.  Eski kodda HER IKI CTE'den
--                      de dusuyordu.
--     220952: 16.333 D50 kurtarma (rateuom bos -> tamami bolus)
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.cho_bolus` AS
SELECT
  ie.stay_id,
  ie.starttime AS give_time,
  CASE
    WHEN d.special_handling = 'd50_rescue'        THEN 'D50_RESCUE'
    WHEN ie.ordercategoryname LIKE '%Antibiotic%' THEN 'DILUENT_ABX'
    WHEN ie.ordercategorydescription = 'Bolus'    THEN 'FLUID_BOLUS'
    ELSE 'BOLUS_OTHER'
  END AS src,
  (CASE ie.amountuom WHEN 'L' THEN ie.amount * 1000 ELSE ie.amount END)
    * d.cho_g_per_l / 1000.0 AS cho_g,           -- GRAM (hiz degil)
  ie.ordercategoryname
FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
JOIN `uludag-mimic-insulin.ref.dextrose_sources` d ON d.itemid = ie.itemid
WHERE ie.statusdescription NOT IN ('Rewritten')
  AND ie.rate IS NULL                            -- <<< ROTA AYRIMI (Kural 14)
  AND ie.amountuom IN ('mL','ml','L','cm3')
  AND ie.amount > 0
  AND (CASE ie.amountuom WHEN 'L' THEN ie.amount*1000 ELSE ie.amount END) <= 1000
  AND d.cho_g_per_l IS NOT NULL
  AND d.verified = TRUE;

-- ---------------------------------------------------------------------
-- 2.3 PARENTERAL DEKSTROZ (gram bazli) -- BULGU-C
--     787 hasta (TPN'in %41'i). Medyan 208 g / 22.8 h = 9.1 g/h
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.cho_pn` AS
SELECT
  ie.stay_id, ie.starttime, ie.endtime, 'PN' AS src,
  ie.amount / (TIMESTAMP_DIFF(ie.endtime, ie.starttime, MINUTE) / 60.0) AS cho_g_hr,
  ie.amount AS pn_grams_total,
  TIMESTAMP_DIFF(ie.endtime, ie.starttime, MINUTE) / 60.0 AS pn_dur_h
FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
WHERE ie.itemid = 225947
  AND ie.statusdescription NOT IN ('Rewritten')
  AND ie.amountuom = 'grams'                     -- 1 adet 'pg' var
  AND ie.amount > 0 AND ie.amount <= 1000        -- 2100 g = 7100 kcal imkansiz
  AND TIMESTAMP_DIFF(ie.endtime, ie.starttime, MINUTE) / 60.0 >= 1.0;

-- ---------------------------------------------------------------------
-- 2.4 ICERIGI BILINMEYEN PN (%59) -- sadece bayrak (Kural 6)
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.pn_unknown` AS
SELECT ie.stay_id, ie.starttime, ie.endtime, ie.rate AS tpn_ml_hr
FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
WHERE ie.itemid IN (225916, 225917, 225920)
  AND ie.statusdescription NOT IN ('Rewritten')
  AND ie.rate > 0
  AND ie.stay_id NOT IN (SELECT stay_id FROM `uludag-mimic-insulin.insulin.cho_pn`);

-- ---------------------------------------------------------------------
-- 2.5 GRID'E BINDIR -- UC AYRI JOIN MANTIGI (Kural 12: ucu icin ayri KK)
-- ---------------------------------------------------------------------
CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.cho`
CLUSTER BY stay_id AS
WITH inf_grid AS (                    -- interval containment
  SELECT g.stay_id, g.hr, g.ts,
    COALESCE(SUM(IF(c.src='EN',           c.cho_g_hr, 0)), 0) AS cho_en,
    COALESCE(SUM(IF(c.src='MAINT',        c.cho_g_hr, 0)), 0) AS cho_maint,
    COALESCE(SUM(IF(c.src='CRRT_CITRATE', c.cho_g_hr, 0)), 0) AS cho_citrate,
    COALESCE(SUM(IF(c.src='EN', c.ml_hr, 0)), 0)              AS en_ml_hr,
    MIN(c.strength_frac)                                      AS min_strength_frac
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `uludag-mimic-insulin.insulin.cho_infusion` c
    ON c.stay_id = g.stay_id AND g.ts >= c.starttime AND g.ts < c.endtime
  GROUP BY 1,2,3
),
pn_grid AS (                          -- interval containment
  SELECT g.stay_id, g.hr, COALESCE(SUM(p.cho_g_hr), 0) AS cho_pn
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `uludag-mimic-insulin.insulin.cho_pn` p
    ON p.stay_id = g.stay_id AND g.ts >= p.starttime AND g.ts < p.endtime
  GROUP BY 1,2
),
pnunk_grid AS (
  SELECT g.stay_id, g.hr,
    MAX(IF(u.stay_id IS NOT NULL, 1, 0)) AS tpn_unknown_flag,
    COALESCE(SUM(u.tpn_ml_hr), 0)        AS tpn_ml_hr
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `uludag-mimic-insulin.insulin.pn_unknown` u
    ON u.stay_id = g.stay_id AND g.ts >= u.starttime AND g.ts < u.endtime
  GROUP BY 1,2
),
bolus_grid AS (                       -- SAAT KOVASI (PATCH-01)
  SELECT g.stay_id, g.hr,
    COALESCE(SUM(IF(b.src='DILUENT_ABX', b.cho_g, 0)), 0) AS cho_diluent_abx,
    COALESCE(SUM(IF(b.src='FLUID_BOLUS', b.cho_g, 0)), 0) AS cho_fluid_bolus,
    COALESCE(SUM(IF(b.src='BOLUS_OTHER', b.cho_g, 0)), 0) AS cho_bolus_other,
    COALESCE(SUM(IF(b.src='D50_RESCUE',  b.cho_g, 0)), 0) AS cho_d50_rescue,
    COUNTIF(b.src='D50_RESCUE')                           AS n_d50_rescue,
    COUNTIF(b.src='DILUENT_ABX')                          AS n_diluent_abx
  FROM `uludag-mimic-insulin.insulin.grid` g
  LEFT JOIN `uludag-mimic-insulin.insulin.cho_bolus` b
    ON b.stay_id = g.stay_id
   AND b.give_time >= g.ts
   AND b.give_time <  TIMESTAMP_ADD(g.ts, INTERVAL 1 HOUR)
  GROUP BY 1,2
),
tot AS (
  SELECT i.*, p.cho_pn, u.tpn_unknown_flag, u.tpn_ml_hr,
         bg.cho_diluent_abx, bg.cho_fluid_bolus, bg.cho_bolus_other,
         bg.cho_d50_rescue, bg.n_d50_rescue, bg.n_diluent_abx,
    -- KURAL 7: D50 kurtarma TOPLAMDA YOK
    i.cho_en + i.cho_maint + i.cho_citrate + p.cho_pn
      + bg.cho_diluent_abx + bg.cho_fluid_bolus + bg.cho_bolus_other AS cho_g_hr
  FROM inf_grid i
  JOIN pn_grid    p  USING (stay_id, hr)
  JOIN pnunk_grid u  USING (stay_id, hr)
  JOIN bolus_grid bg USING (stay_id, hr)
)
SELECT *,
  cho_g_hr - LAG(cho_g_hr, 2) OVER w AS d_cho_2h,
  cho_g_hr - LAG(cho_g_hr, 4) OVER w AS d_cho_4h,
  SUM(cho_g_hr) OVER (PARTITION BY stay_id ORDER BY hr
                      ROWS BETWEEN 23 PRECEDING AND CURRENT ROW) AS cho_24h_cum,
  IF(MAX(cho_g_hr) OVER (PARTITION BY stay_id ORDER BY hr
        ROWS BETWEEN 3 PRECEDING AND CURRENT ROW) < 2.0, 1, 0) AS cho_interrupted_4h,
  hr - MAX(IF(cho_g_hr >= 2.0, hr, NULL)) OVER (
        PARTITION BY stay_id ORDER BY hr ROWS UNBOUNDED PRECEDING) AS fasting_hr,
  IF(min_strength_frac < 1.0, 1, 0) AS en_diluted_flag,     -- BULGU-A duyarlilik
  CASE
    WHEN cho_en >= 2 AND (cho_pn > 0 OR tpn_unknown_flag = 1) THEN 'EN+PN'
    WHEN cho_en >= 2                                          THEN 'EN'
    WHEN cho_pn > 0 OR tpn_unknown_flag = 1                   THEN 'PN'
    WHEN cho_maint > 0 OR cho_diluent_abx > 0                 THEN 'MAINT_ONLY'
    ELSE 'none'
  END AS nutrition_route
FROM tot
WINDOW w AS (PARTITION BY stay_id ORDER BY hr);

-- =====================================================================
-- 2.6 KALITE KONTROLLERI -- Kural 12: her join mantigi icin AYRI KK
-- =====================================================================
-- KK-1: EN medyani. BAND 4-12 g/h, medyan ~7. (ESKI 8-15 YANLISTI:
--       gercek urun x hiz hesabi agirlikli medyan ~7 g/h verdi -> yanlis
--       durdurma yaratirdi.) Disindaysa DUR.
-- SELECT nutrition_route, COUNT(*) n,
--        APPROX_QUANTILES(cho_g_hr, 10) AS cho_deciles,
--        AVG(cho_interrupted_4h) AS pct_interrupted,
--        AVG(en_diluted_flag)    AS pct_diluted
-- FROM `uludag-mimic-insulin.insulin.cho` GROUP BY 1 ORDER BY n DESC;

-- KK-2 (PATCH-01): bolus yakalama orani = 1.0 olmali (eski kodda ~0.25)
-- SELECT
--   (SELECT COUNT(*) FROM `uludag-mimic-insulin.insulin.cho_bolus`
--     WHERE src='D50_RESCUE') AS n_bolus_events,
--   (SELECT SUM(n_d50_rescue) FROM `uludag-mimic-insulin.insulin.cho`) AS n_on_grid;

-- KK-3 (PATCH-08): antibiyotik sulandiricisi grid'e bindi mi?
--       Beklenen ~49.463 hasta, hasta basina ~9 olay
-- SELECT COUNT(DISTINCT stay_id) AS n_stays_with_diluent,
--        SUM(n_diluent_abx)      AS n_diluent_on_grid,
--        APPROX_QUANTILES(cho_diluent_abx, 4) AS diluent_g_hr_q
-- FROM `uludag-mimic-insulin.insulin.cho` WHERE n_diluent_abx > 0;

-- KK-4 (BULGU-C): PN g/h makul mu? Beklenen medyan ~9 g/h, 787 hasta
-- SELECT APPROX_QUANTILES(cho_g_hr, 4) AS pn_g_hr_q,
--        COUNT(DISTINCT stay_id) AS n_pn_stays
-- FROM `uludag-mimic-insulin.insulin.cho_pn`;

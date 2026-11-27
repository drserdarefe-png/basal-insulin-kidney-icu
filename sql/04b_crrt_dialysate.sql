-- =====================================================================
-- 04b — CRRT DIYALIZAT/REPLASMAN TANISI (z_crrt yorumu icin)
-- crrt tablosunda dialysate_fluid/replacement_fluid/rate SUTUNLARI VAR.
-- Ciktilar: qc/04_crrt_dialysate.json
-- BULGU: fluid etiketi yalniz POTASYUM tiyeri (Prismasate K2/K4/K0);
--   glukozlu (B22GK4/0) vs glukozsuz (BGK0/0) urun kodu AYIRT EDILEMEZ.
-- =====================================================================

-- SORGU 1: fluid urun dagilimi
WITH c AS (
  SELECT cr.* FROM `physionet-data.mimiciv_3_1_derived.crrt` cr
  JOIN `uludag-mimic-insulin.insulin.cohort_final` f USING (stay_id)
)
SELECT 'dialysate_fluid' AS src, dialysate_fluid AS fluid,
       COUNT(*) n, COUNT(DISTINCT stay_id) n_stays
FROM c WHERE dialysate_fluid IS NOT NULL GROUP BY 1,2
UNION ALL
SELECT 'replacement_fluid', replacement_fluid, COUNT(*), COUNT(DISTINCT stay_id)
FROM c WHERE replacement_fluid IS NOT NULL GROUP BY 1,2
ORDER BY n DESC;

-- SORGU 2: rate quantilleri (mL/h)
-- WITH c AS (SELECT cr.* FROM derived.crrt cr JOIN cohort_final f USING(stay_id))
-- SELECT APPROX_QUANTILES(dialysate_rate,10), APPROX_QUANTILES(replacement_rate,10),
--        APPROX_QUANTILES(prefilter_replacement_rate,10),
--        APPROX_QUANTILES(postfilter_replacement_rate,10), APPROX_QUANTILES(citrate,4)
-- FROM c;

-- SORGU 3: eGFR<15 x on_crrt 2x2 (hasta duzeyi, 07 Z ile ayni: egfr=median, crrt=max)
-- WITH pat AS (SELECT subject_id, APPROX_QUANTILES(egfr_now,2)[OFFSET(1)] egfr_med,
--                     MAX(on_crrt) crrt_any FROM model_set GROUP BY 1)
-- SELECT IF(egfr_med<15,'egfr_lt15','egfr_ge15'), IF(crrt_any=1,'crrt','no_crrt'),
--        COUNT(*) FROM pat GROUP BY 1,2;

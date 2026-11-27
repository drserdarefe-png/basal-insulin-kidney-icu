-- =====================================================================
-- 00 — ITEMID KESIF SORGULARI
-- MIMIC-IV v3.1 | proje: uludag-mimic-insulin
-- Ciktiyi ref/*.csv dosyalarina isle, verified=TRUE yap (SEN onaylarsin).
-- Her sorgu ONCE --dry_run ile calistirilir (Kural 11).
-- NOT: 00 kohorttan ONCE calisir; kohort filtresi uygulanamaz.
-- =====================================================================

-- ---------- 00.1  INSULIN ----------
SELECT d.itemid, d.label, d.abbreviation, d.category, d.unitname,
       COUNT(*) AS n_events,
       COUNT(DISTINCT ie.stay_id) AS n_stays,
       COUNTIF(ie.rate IS NOT NULL) AS n_with_rate,
       COUNTIF(ie.rate IS NULL)     AS n_without_rate,
       APPROX_QUANTILES(ie.amount, 4) AS amount_q,
       ARRAY_AGG(DISTINCT ie.amountuom IGNORE NULLS) AS uoms,
       ARRAY_AGG(DISTINCT ie.ordercategorydescription IGNORE NULLS) AS ordercat
FROM `physionet-data.mimiciv_3_1_icu.d_items` d
JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie USING (itemid)
WHERE LOWER(d.label) LIKE '%insulin%'
GROUP BY 1,2,3,4,5
ORDER BY n_events DESC;

-- ---------- 00.2  ENTERAL NUTRISYON ----------
SELECT d.itemid, d.label, d.category, d.unitname,
       COUNT(*) AS n_events, COUNT(DISTINCT ie.stay_id) AS n_stays,
       APPROX_QUANTILES(ie.rate, 4) AS rate_q,
       ARRAY_AGG(DISTINCT ie.rateuom IGNORE NULLS) AS rateuoms
FROM `physionet-data.mimiciv_3_1_icu.d_items` d
JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie USING (itemid)
WHERE d.category LIKE '%Nutrition%'
   OR REGEXP_CONTAINS(LOWER(d.label),
        r'nutren|jevity|osmolite|vital|glucerna|impact|peptamen|isosource|'
        r'promote|pulmocare|nepro|replete|boost|ensure|tube feed|enteral|'
        r'fibersource|two ?cal|novasource|diabetisource|perative|oxepa')
GROUP BY 1,2,3,4
ORDER BY n_events DESC;

-- ---------- 00.3  DEKSTROZ ICEREN HER SEY ----------
SELECT d.itemid, d.label, d.category, d.unitname,
       COUNT(*) AS n_events, COUNT(DISTINCT ie.stay_id) AS n_stays,
       ARRAY_AGG(DISTINCT ie.rateuom IGNORE NULLS) AS rateuoms,
       ARRAY_AGG(DISTINCT ie.amountuom IGNORE NULLS) AS amtuoms
FROM `physionet-data.mimiciv_3_1_icu.d_items` d
JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie USING (itemid)
WHERE REGEXP_CONTAINS(LOWER(d.label),
        r'dextrose|d5|d10|d20|d25|d50|glucose|tpn|ppn|parenteral|'
        r'amino ?acid|clinimix|kabiven|acd-a|acd a|citrate|prismasate|'
        r'prismocitrate|regiocit|dialysate|replacement fluid')
GROUP BY 1,2,3,4
ORDER BY n_events DESC;

-- ---------- 00.4  GLUKOZ OLCUM KAYNAKLARI ----------
SELECT 'chartevents' AS src, d.itemid, d.label, COUNT(*) n,
       APPROX_QUANTILES(ce.valuenum, 4) AS q
FROM `physionet-data.mimiciv_3_1_icu.d_items` d
JOIN `physionet-data.mimiciv_3_1_icu.chartevents` ce USING (itemid)
WHERE d.itemid IN (225664, 220621, 226537, 228388)
GROUP BY 1,2,3
UNION ALL
SELECT 'labevents', d.itemid, d.label, COUNT(*),
       APPROX_QUANTILES(le.valuenum, 4)
FROM `physionet-data.mimiciv_3_1_hosp.d_labitems` d
JOIN `physionet-data.mimiciv_3_1_hosp.labevents` le USING (itemid)
WHERE d.itemid IN (50809, 50931, 52027, 52569)
GROUP BY 1,2,3
ORDER BY n DESC;

-- ---------- 00.5  STEROID (BULGU-02: potensler duzeltilecek) ----------
SELECT d.itemid, d.label, d.unitname, COUNT(*) n
FROM `physionet-data.mimiciv_3_1_icu.d_items` d
JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie USING (itemid)
WHERE REGEXP_CONTAINS(LOWER(d.label),
        r'hydrocortisone|methylprednisolone|solu-?medrol|solu-?cortef|'
        r'dexamethasone|prednisone|prednisolone')
GROUP BY 1,2,3 ORDER BY n DESC;

-- ---------- 00.6  PATCH-05 AYRISMA DENETIMI ----------
-- ref/insulin_pk.csv BigQuery'ye yuklendikten SONRA calistir.
-- BOS DONMELI. Donmuyorsa ref eksik -> 01b ile 03 ayrisir (sessiz hata).
SELECT d.itemid, d.label, COUNT(*) AS n_events
FROM `physionet-data.mimiciv_3_1_icu.d_items` d
JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie USING (itemid)
LEFT JOIN `uludag-mimic-insulin.ref.insulin_pk` pk ON pk.itemid = d.itemid
WHERE LOWER(d.label) LIKE '%insulin%'
  AND pk.itemid IS NULL
GROUP BY 1,2
ORDER BY n_events DESC;

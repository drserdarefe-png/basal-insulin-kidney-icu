-- =====================================================================
-- 01 — KOHORT + SAATLIK GRID                              v2.0
-- MIMIC-IV v3.1 | proje: uludag-mimic-insulin
-- Dislamalar EXCLUSION_REGISTRY.md v1.1 (EXC-xx) ile birebir.
-- Uygulanmis: BULGU-03 (EXC-01 OR mantigi + anyon acigi + kapsam),
--             PATCH-03 (DKA triadi +-6h zamansal bag)
-- NOT: cohort_final ve grid burada DEGIL, 01b'de kurulur (EXC-03 sonrasi).
-- =====================================================================

CREATE OR REPLACE TABLE `uludag-mimic-insulin.insulin.cohort_raw` AS
WITH base AS (
  SELECT
    ic.subject_id, ic.hadm_id, ic.stay_id,
    ic.intime, ic.outtime, ic.los,
    ROW_NUMBER() OVER (PARTITION BY ic.subject_id ORDER BY ic.intime) AS icu_seq
  FROM `physionet-data.mimiciv_3_1_icu.icustays` ic
),
cohort_hadm AS (
  SELECT DISTINCT hadm_id FROM base          -- tarama kapsami (Kural 11)
),

-- ---- ICD tabanli bayraklar ------------------------------------------
dx AS (
  SELECT d.hadm_id,
    -- EXC-01 DKA/HHS  (ICD kolu)
    LOGICAL_OR(REGEXP_CONTAINS(d.icd_code, r'^(E10[01]|E11[01]|E13[01]|E1[034][01])')
            OR REGEXP_CONTAINS(d.icd_code, r'^(2501|2502|2503)'))            AS f_dka,
    -- EXC-02 gebelik
    LOGICAL_OR(REGEXP_CONTAINS(d.icd_code, r'^(O[0-9]|O[AB]|Z34|Z33)')
            OR REGEXP_CONTAINS(d.icd_code, r'^(63[0-9]|64[0-9]|65[0-9]|V22|V23)')) AS f_preg,
    -- EXC-06 akut KC yetmezligi
    LOGICAL_OR(REGEXP_CONTAINS(d.icd_code, r'^(K720|K721|K711)')
            OR REGEXP_CONTAINS(d.icd_code, r'^(5720|5722|5701)'))            AS f_alf,
    -- EXC-08 pankreatik yetmezlik
    LOGICAL_OR(REGEXP_CONTAINS(d.icd_code, r'^(K863|K861|K903)')
            OR REGEXP_CONTAINS(d.icd_code, r'^(5773|5771)'))                 AS f_panc_insuf
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d
  JOIN cohort_hadm USING (hadm_id)
  GROUP BY 1
),
px AS (
  SELECT p.hadm_id,
    LOGICAL_OR(REGEXP_CONTAINS(p.icd_code, r'^(0FY0)') OR p.icd_code IN ('5059','5051'))
      AS f_liver_tx,
    LOGICAL_OR(REGEXP_CONTAINS(p.icd_code, r'^(0FTG)') OR REGEXP_CONTAINS(p.icd_code, r'^(526)'))
      AS f_pancreatectomy
  FROM `physionet-data.mimiciv_3_1_hosp.procedures_icd` p
  JOIN cohort_hadm USING (hadm_id)
  GROUP BY 1
),

-- ---- EXC-01 LAB KOLU: v1.2'de KALDIRILDI (asagidaki yorum) --------

-- ---- EXC-05 glukoz olcum sayisi -------------------------------------
glu_n AS (
  SELECT ce.stay_id, COUNT(*) AS n_glu
  FROM `physionet-data.mimiciv_3_1_icu.chartevents` ce
  WHERE ce.itemid IN (225664, 220621, 226537)
    AND ce.valuenum BETWEEN 10 AND 2000
  GROUP BY 1
),

-- ---- EXC-09 nutrisyon dokumantasyonu --------------------------------
nut_doc AS (
  SELECT ie.stay_id, COUNT(*) AS n_nut
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN `physionet-data.mimiciv_3_1_icu.d_items` d USING (itemid)
  WHERE d.category LIKE '%Nutrition%'
  GROUP BY 1
),

-- ---- EXC-10 insulin birim butunlugu (ref tabanli, PATCH-05) ---------
ins_uom AS (
  SELECT ie.stay_id,
         LOGICAL_OR(ie.amountuom IS NULL OR ie.amountuom NOT IN ('units','U')) AS f_bad_uom
  FROM `physionet-data.mimiciv_3_1_icu.inputevents` ie
  JOIN `uludag-mimic-insulin.ref.insulin_pk` pk ON pk.itemid = ie.itemid
  WHERE pk.verified = TRUE
  GROUP BY 1
)

SELECT
  b.subject_id, b.hadm_id, b.stay_id, b.intime, b.outtime, b.los, b.icu_seq,
  -- her dislama ayri sutun -> attrition tablosu buradan
  -- EXC-01: registry = ICD **VEYA** lab triadi (BULGU-03: eski kod AND idi)
  -- REGISTRY v1.2 AMENDMANI: lab kolu KALDIRILDI. Yalniz ICD.
  --   Gerekce (qc/01c_lab_only_who.json): lab-only 136 hastanin 80'inde
  --   (%58.8) DIYABET YOK; ortalama maks laktat 9.6 mmol/L; %72'sinde
  --   laktat>4. Anyon acigini keton degil LAKTAT yukseltiyor.
  --   'glu>600' zaten DKA degil HHS esigi (DKA tipik 350-500; SGLT2i
  --   alaninda oglisemik DKA <250). Triad, agir sok + stres
  --   hiperglisemisini yakaliyordu = CALISMA POPULASYONUNUN KENDISI.
  --   Yanlis pozitifler rastgele degil, sistematik olarak EN AGIR
  --   hiperglisemik hastalar -> secim yanliligi.
  --   Keton kurtaramaz: BHB n_hadm=1051 (%1.1); idrar ketonu bol ama
  --   YBU'de spesifik degil (aclik ketozu).
  --   Dislamada SPESIFIKLIK > duyarlilik.
  COALESCE(dx.f_dka, FALSE) AS exc_01_dka,
  COALESCE(dx.f_preg, FALSE)                                 AS exc_02_pregnancy,
  b.los < 1                                                  AS exc_04_los,
  COALESCE(gn.n_glu, 0) < 4                                  AS exc_05_few_glu,
  COALESCE(dx.f_alf, FALSE)                                  AS exc_06_alf,
  COALESCE(px.f_liver_tx, FALSE)                             AS exc_07_liver_tx,
  COALESCE(px.f_pancreatectomy, FALSE)
    OR COALESCE(dx.f_panc_insuf, FALSE)                      AS exc_08_panc,
  COALESCE(nd.n_nut, 0) = 0                                  AS exc_09_no_nut_doc,
  COALESCE(iu.f_bad_uom, FALSE)                              AS exc_10_bad_uom,
  b.icu_seq > 1                                              AS exc_11_readmit,
  -- teshis: ICD kolu mu lab kolu mu tetikledi
  COALESCE(dx.f_dka, FALSE)     AS dka_by_icd,
  FALSE AS dka_by_lab   -- v1.2: retired
FROM base b
LEFT JOIN dx      ON dx.hadm_id = b.hadm_id
LEFT JOIN px      ON px.hadm_id = b.hadm_id
LEFT JOIN glu_n   gn ON gn.stay_id = b.stay_id
LEFT JOIN nut_doc nd ON nd.stay_id = b.stay_id
LEFT JOIN ins_uom iu ON iu.stay_id = b.stay_id;

-- ---------------------------------------------------------------------
-- KK: DKA kollarinin ayrisimi (PATCH-03 etkisi CONSORT dipnotuna)
-- ---------------------------------------------------------------------
-- SELECT
--   COUNTIF(dka_by_icd) AS n_icd_only_or_both,
--   COUNTIF(dka_by_lab) AS n_lab_triad,
--   COUNTIF(dka_by_icd AND dka_by_lab) AS n_both,
--   COUNTIF(exc_01_dka) AS n_exc_01_total
-- FROM `uludag-mimic-insulin.insulin.cohort_raw`;

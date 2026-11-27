# Insülin Doz-Yanıt Çalışması — Dışlama Kriteri Registry
Sürüm: v1.3 — 16,07,26 (Makale 2 taslak sabitleme)
Sorumlu: S. Efe
Durum: **DONDURULMADI** — İmza bekliyor (bkz. altta imza satırı)

## Veri kaynağı ve tanım beyanları (v1.3 sabit)
- Veri: **MIMIC-IV v3.1** (`physionet-data.mimiciv_3_1_{icu,hosp,derived}`), sürüm sabit.
- Makale 2 analiz seti: **7.147 glarjin doz olayı / 1.074 hasta** (gece-glukozlu),
  kohort akışı 94.458 ICU yatışından (bkz. qc/paper2_tables/T1_cohort.md).
- **MELD-Na = UNOS 2016** (OPTN Policy 9), tam sayı; Na düzeltmesi yalnız MELD>11.
- **eGFR = CKD-EPI 2021 ırksız** (Inker NEJM 2021), sabit 142.
- **İnsülin tek-olay doz tavanı = 50 U** (bolus units_total>50 veya infüzyon
  rate>50 U/h dışlanır; Kural 15, qc/03_dose_outliers.txt). Basal-spesifik daha
  yüksek tavan S. Efe kararına açık (glarjin 50-180 U kuyruğu düştü).

Bu dosya TEK doğruluk kaynağıdır (CLAUDE.md Kural 1). SQL'deki her dışlama
satırı buradaki bir ID'ye referans verir (`-- EXC-xx`). Satır SİLİNMEZ;
`status` alanı `retired` yapılır ve gerekçe yazılır.

## Veri kaynağı
MIMIC-IV **v3.1** — `physionet-data.mimiciv_3_1_{icu,hosp,derived}`
Erişim: 15,07,26 | CITI Record ID 76867340 | n_all_icu_stays = 94.458
Sürüm analiz boyunca sabittir (Kural 13).

## Dışlama kriterleri

| ID | Kriter | Operasyonel tanım | Gerekçe | n | Status |
|----|--------|-------------------|---------|---|--------|
| EXC-01 | DKA / HHS | ICD E10.1x, E11.0x, E11.1x, E13.1x, E1[034][01], 2501-2503. **v1.2: lab kolu KALDIRILDI** (aşağıya bkz.) | Farklı doz rejimi (IV protokol); hedef değişkenin anlamı bozulur | 2.174 | active |
| EXC-02 | Gebelik | ICD O00–O9A, Z33/Z34, 63x–65x, V22/V23 | Ayrı glisemik hedef | — | active |
| EXC-03 | İnsülin pompası / belgelenmemiş maruziyet | `diagnoses_icd` Z96.41/V45.85 **VE** ICU insülin dozu=0 **VE** glu_max>180 (01b `A_pump_undocumented`). *v1.0: "chartevents/notes pompa kaydı" — retired* | Doz kaydı eksik/dış sistem → hedef kirlenmesi | ? | active |
| EXC-04 | ICU LOS < 24 h | `icustays.los < 1` | Yetersiz zaman serisi | — | active |
| EXC-05 | Glukoz ölçümü < 4 | chartevents 225664/220621/226537, 10–2000 aralığı | Doz-yanıt eğrisi kurulamaz | — | active |
| EXC-06 | Akut karaciğer yetmezliği | ICD K72.0x, K72.1x, K71.1x, 5720/5722/5701 | İnsülin değil dekstroz alıyor; ters nedensellik | — | active |
| EXC-07 | KC transplant peri-op | prosedür ICD 0FY0x, 5059/5051 | Ekstrem, atipik glisemi dinamiği | — | active |
| EXC-08 | Pankreatektomi / pankreatik yetmezlik | prosedür 0FTGx, 526x; tanı K863/K861/K903, 5773/5771 | Endojen insülin = 0, ayrı popülasyon | — | active |
| EXC-09 | Nütrisyon rotası hiç belgelenmemiş | `d_items.category LIKE '%Nutrition%'` kaydı yok | Kayıt kalitesi proxy'si; CHO akışı hesaplanamaz | — | active |
| EXC-10 | İnsülin birim bütünlüğü | `amountuom` NULL veya ∉ {units, U} | Hedef kirlenmesi | — | active |
| EXC-11 | Tekrar yatış (2.+) | `icustays` sıralaması | Bağımlılık; duyarlılık analizinde geri alınabilir | — | active |

## DIŞLANMAYANLAR (bilinçli karar — çalışmanın konusu)
- Siroz / kronik KC hastalığı → **stratum**
- KBY (tüm evreler), CRRT, IHD → **stratum**
- Nütrisyon almayan hasta → **stratum** (`nutrition_route = 'none'`)
- Steroid alan hasta → kovaryat + etkileşim (INT-04)
- Tip 1 / Tip 2 DM → kovaryat
- `B_prior_pump_converted` (pompa ICD var, YBÜ dozu belgeleniyor) → **TUT**,
  `cov_prior_pump` kovaryatı. Ağırlıklı T1DM; `z_dm_t1` ile örtüşür ama aynı değil.
- `C_suspect_undocumented` (pompa ICD yok, doz=0, glu>180, DM var) → **ANA
  ANALİZDE TUT**. EXC-03 ile aynı hedef-kirlenmesi riskini taşır ancak pompa
  kanıtı yok. Bu asimetri bilinçlidir. Duyarlılık analizinde çıkarılır
  (`sens_suspect_undoc`). 06'da `dose_bin='0_no_dose'` hücresini kirletebilir
  → ayrı bayrakla raporlanır.

## Ön-tanımlı (pre-specified) etkileşim hipotezleri
| ID | Hipotez | Yön | Durum |
|----|---------|-----|-------|
| INT-01 | `sirroz × CHO_kesintisi(>4h)` → hipoglisemi riski | multiplikatif ↑ | **RETIRED (v1.3)** |
| INT-02 | `eGFR<15 × doz` → ISF büyüklüğü | ↑ (klirens ↓) | **bazal analize taşındı (v1.3)** |
| INT-03 | `üremi(BUN yüksek) × doz` → ISF büyüklüğü | ↓ (periferik direnç) | active |
| INT-04 | `steroid × CHO_g_hr` → glukoz artışı | ↑ | active |
| INT-05 | `parmak_ucu_ölçüm × vazopresör` → ölçüm sapması | bias | active |

**INT-01 RETIRED (v1.3) — iki gerekçe:**
1. **Positivity <3 seviye** üç kapsamda da (full/fed/clean): `cirr&interrupt` hücresi
   Kural 10 eşiğini (≥3 destekli doz seviyesi) geçemedi (bkz. qc/scope_*/model_card_support.md).
2. **Maruziyet geçersiz:** `cho_interrupted_4h` bayrağı, beslenmesi kesilen (NPO) hastayı
   güvenilir ayırmıyor — bayrak yiyen hastada positivity/örtüşme 0,847 vs aç hastada 0,708;
   maruziyet tanımı NPO'yu temsil etmiyor. Kesinti kolu güç yok + tanım kirli.

**INT-02 → bazal analize taşındı (v1.3):** Makale 2 (glarjin) `dose_u_kg × egfr_band`
etkileşimi olarak; T3 (kılavuz-pratik açığı) + AIPW/within-patient eGFR-bandı
sonuçlarında raporlanıyor. Bolus ISF modeli birincil olmaktan çıkarıldığı için
INT-02'nin bolus formu askıya alındı.

Analizden ÖNCE dondurulur. Post-hoc bulunanlar ayrı "exploratory" başlığında.

---
## İMZA (dondurma)
- [ ] v1.3 FROZEN — S. Efe imza, tarih: __________

## Tanım beyanları
- **MELD-Na: UNOS 2016** (OPTN Policy 9), Kim 2008 NEJM formu DEĞİL.
  Tam sayı MELD; Na düzeltmesi yalnız MELD>11 iken. Kim 2008 duyarlılık
  analizinde. Fark düşük-MELD bölgesinde: kompanse sirotik Kim=14,5 vs UNOS=11,0.
- **eGFR: CKD-EPI 2021 ırksız** (Inker NEJM 2021), sabit 142.
- **MELD kreatinin: RRT'de 4,0'a sabitlenir** (OPTN). İlk 24h CRRT penceresi.
- **Karışım insülinleri** (70/30, Humalog 75/25) `is_basal=TRUE`, τ=300 dk.
  426 hasta, marjinal etki, klinik yargı.

---

## Sürüm notları

### v1.2 (16,07,26/02:00) — EXC-01 lab kolu retired

v1.1: `ICD VEYA (glu>600 & pH<7.30 & AG>16)`
v1.2: **yalnız ICD**

Kanıt: `qc/01_dka.json`, `qc/01b_ketones.json`, `qc/01c_lab_only_who.json`

- ICD–lab uyumu **%7,1** (2.174 ICD'nin 154'ü lab triadını karşılıyor)
- Lab-only 136 hastanın **80'inde (%58,8) diyabet ICD'si YOK** → DKA olamaz
- Ortalama maks laktat **9,6 mmol/L**; %72'sinde laktat >4 → laktik asidoz
- Anyon açığını keton değil **laktat** yükseltiyor; triad keton bileşeni
  olmadan DKA için spesifik değil
- `glu>600` DKA değil **HHS** eşiği (DKA tipik 350–500; SGLT2i alanında
  öglisemik DKA <250)
- Yanlış pozitifler sistematik olarak **en ağır hiperglisemik hastalar** =
  çalışma popülasyonunun kendisi → seçim yanlılığı
- Keton ile düzeltme uygulanamaz: BHB `n_hadm=1.051` (%1,1); idrar ketonu
  bol (159.282 hadm) ama YBÜ'de spesifik değil (açlık ketozu, NPO)

Etki: n_exc_01 **2.310 → 2.174** (136 hasta geri kazanıldı)
Sınırlılık: kodlanmamış hafif DKA kohortta kalır. Bu hastalar IV protokolde
farklı doz deseni gösterir, 06 positivity'de görünür. Dışlamada
**spesifiklik > duyarlılık**.

### v1.1 (15,07,26/22:45) — EXC-03 op-tanımı revize
v1.0 "chartevents/notes pompa kaydı" → ICD (Z96.41/V45.85) + belgelenmemiş-
maruziyet imzası. Gerekçe: MIMIC-IV'te güvenilir yapılandırılmış pompa
itemid'i yok; notes NLP tekrar-üretilebilirliği düşük; EXC-03'ün asıl amacı
"doz kaydı eksik/dış sistem" ve ICD+imza bunu doğrudan hedefler, `pump_class`
A–E ile denetlenebilir kılar.

### v1.0 (15,07,26/00:00) — ilk sürüm

---

## Süreç notu
v1.0 ve v1.1 tanımları MIMIC-IV yapısı görülmeden yazıldı. Registry **henüz
dondurulmadı**, bu yüzden v1.2 düzeltmesi meşru. **Dondurmadan sonra** aynı
hamle protokol amendmanı olur, ayrı raporlanır, makalede beyan edilir
(WORKFLOW.md rol tablosu).

Registry'nin sahibi **S. Efe**'dir. Kod yazan aktör (Claude Code / Science)
registry'yi değiştiremez, yalnız **önerir**. Dondurma imzası:

- [ ] 01b ADIM 0 prevalansı → A sınıfı % registry'ye
- [ ] ref/*.csv `verified=TRUE` onayı
- [ ] v1.3 FROZEN — imza, tarih

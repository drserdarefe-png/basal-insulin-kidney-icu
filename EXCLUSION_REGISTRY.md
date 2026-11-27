# Insülin Doz-Yanıt Çalışması: Dışlama Kriteri Registry
Sürüm: v1.3 (16.07.2026; kriter tanımları bu tarihten sonra değişmedi)
Sorumlu: S. Efe
Durum: **SON SÜRÜM, KAPATILDI** (27.11.2026, S. Efe). Ayrıntı: en alttaki süreç notu.

## Veri kaynağı ve tanım beyanları (v1.3 sabit)
- Veri: **MIMIC-IV v3.1** (`physionet-data.mimiciv_3_1_{icu,hosp,derived}`), sürüm sabit.
- Analiz seti (makaledeki son sayılar): 94.458 ICU yatışından 7.876 uygun yatış;
  **1.117 hasta, 7.954 glarjin doz olayı**; gece glukozu olan **7.560 olay (1.090 hasta)**;
  çift sağlam (AIPW) modellerde **7.517 olay (1.085 hasta)** (makale Tablo 1, Ek Tablo 1).
- **MELD-Na = UNOS 2016** (OPTN Policy 9), tam sayı; Na düzeltmesi yalnız MELD>11.
- **eGFR = CKD-EPI 2021 ırksız** (Inker NEJM 2021), sabit 142.
- **İnsülin tek-olay doz tavanı (sınıfa özgü, son hal):** infüzyon >50 U/h, bazal
  olmayan bolus >50 U ve bazal insülin >150 U (75 kg'da yaklaşık 2 U/kg) dışlanır.
  İlk uygulanan tek 50 U tavanı, obez ve insüline dirençli hastalardaki gerçek yüksek
  bazal dozları kestiği için bu şekilde değiştirildi (aşağıda sürüm notu).

Bu dosya dışlama kriterleri için tek doğruluk kaynağıdır. SQL'deki her dışlama
satırı buradaki bir ID'ye referans verir (`-- EXC-xx`). Satır SİLİNMEZ;
`status` alanı `retired` yapılır ve gerekçe yazılır.

## Veri kaynağı
MIMIC-IV **v3.1**: `physionet-data.mimiciv_3_1_{icu,hosp,derived}`
Erişim: 15,07,26 | CITI Record ID 76867340 | n_all_icu_stays = 94.458
Sürüm analiz boyunca sabittir.

## Dışlama kriterleri

| ID | Kriter | Operasyonel tanım | Gerekçe | n | Status |
|----|--------|-------------------|---------|---|--------|
| EXC-01 | DKA / HHS | ICD E10.1x, E11.0x, E11.1x, E13.1x, E1[034][01], 2501-2503. **v1.2: lab kolu KALDIRILDI** (aşağıya bkz.) | Farklı doz rejimi (IV protokol); hedef değişkenin anlamı bozulur | 2.174 | active |
| EXC-02 | Gebelik | ICD O00–O9A, Z33/Z34, 63x–65x, V22/V23 | Ayrı glisemik hedef | - | active |
| EXC-03 | İnsülin pompası / belgelenmemiş maruziyet | `diagnoses_icd` Z96.41/V45.85 **VE** ICU insülin dozu=0 **VE** glu_max>180 (01b `A_pump_undocumented`). *v1.0: "chartevents/notes pompa kaydı", retired* | Doz kaydı eksik/dış sistem → hedef kirlenmesi | ? | active |
| EXC-04 | ICU LOS < 24 h | `icustays.los < 1` | Yetersiz zaman serisi | - | active |
| EXC-05 | Glukoz ölçümü < 4 | chartevents 225664/220621/226537, 10–2000 aralığı | Doz-yanıt eğrisi kurulamaz | - | active |
| EXC-06 | Akut karaciğer yetmezliği | ICD K72.0x, K72.1x, K71.1x, 5720/5722/5701 | İnsülin değil dekstroz alıyor; ters nedensellik | - | active |
| EXC-07 | KC transplant peri-op | prosedür ICD 0FY0x, 5059/5051 | Ekstrem, atipik glisemi dinamiği | - | active |
| EXC-08 | Pankreatektomi / pankreatik yetmezlik | prosedür 0FTGx, 526x; tanı K863/K861/K903, 5773/5771 | Endojen insülin = 0, ayrı popülasyon | - | active |
| EXC-09 | Nütrisyon rotası hiç belgelenmemiş | `d_items.category LIKE '%Nutrition%'` kaydı yok | Kayıt kalitesi proxy'si; CHO akışı hesaplanamaz | - | active |
| EXC-10 | İnsülin birim bütünlüğü | `amountuom` NULL veya ∉ {units, U} | Hedef kirlenmesi | - | active |
| EXC-11 | Tekrar yatış (2.+) | `icustays` sıralaması | Bağımlılık; duyarlılık analizinde geri alınabilir | - | active |

## DIŞLANMAYANLAR (bilinçli karar: çalışmanın konusu)
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

**INT-01 RETIRED (v1.3), iki gerekçe:**
1. **Positivity <3 seviye** üç kapsamda da (full/fed/clean): `cirr&interrupt` hücresi
   önceden belirlenen eşiği (≥3 destekli doz seviyesi) geçemedi.
2. **Maruziyet geçersiz:** `cho_interrupted_4h` bayrağı, beslenmesi kesilen (NPO) hastayı
   güvenilir ayırmıyor: bayrak yiyen hastada positivity/örtüşme 0,847 vs aç hastada 0,708;
   maruziyet tanımı NPO'yu temsil etmiyor. Kesinti kolu güç yok + tanım kirli.

**INT-02 → bazal analize taşındı (v1.3):** Makale 2 (glarjin) `dose_u_kg × egfr_band`
etkileşimi olarak; T3 (kılavuz-pratik açığı) + AIPW/within-patient eGFR-bandı
sonuçlarında raporlanıyor. Bolus ISF modeli birincil olmaktan çıkarıldığı için
INT-02'nin bolus formu askıya alındı.

Analizden ÖNCE dondurulur. Post-hoc bulunanlar ayrı "exploratory" başlığında.

---
## Kapanış
- [x] v1.3 son sürüm olarak kapatıldı: S. Efe, 27.11.2026

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

### Doz tavanı notu (dışlama kriteri değil, veri temizliği)
Tek 50 U tavanı yerine sınıfa özgü tavan getirildi (bazal 150 U). Bu değişiklik dört bazal olayı
(iki glarjin, iki NPH) dışlar ve 50 U üstündeki 386 glarjin olayını analize geri kazandırır.
AIPW, hasta-içi ve doz-belirleyici analizler son veriyle yeniden çalıştırıldı (27.11.2026,
`rerun_27.11.2026.sh`); makalede yalnız bu sonuçlar raporlanıyor. EXC-01 ile EXC-11 arasındaki
kriterler bu değişiklikten etkilenmedi.

### v1.2 (16.07.2026 02:00): EXC-01 lab kolu retired

v1.1: `ICD VEYA (glu>600 & pH<7.30 & AG>16)`
v1.2: **yalnız ICD**

Kanıt: yerel kalite kontrol çıktıları (hasta düzeyinde veri içerdiği için depoda yok).

- ICD–lab uyumu **%7,1** (2.174 ICD'nin 154'ü lab triadını karşılıyor)
- Lab-only 136 hastanın **80'inde (%58,8) diyabet ICD'si YOK** → DKA olamaz
- Ortalama maks laktat **9,6 mmol/L**; %72'sinde laktat >4 → laktik asidoz
- Anyon gap'i keton değil **laktat** yükseltiyor; triad keton bileşeni
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

### v1.1 (15.07.2026 22:45): EXC-03 operasyonel tanımı revize
v1.0 "chartevents/notes pompa kaydı" → ICD (Z96.41/V45.85) + belgelenmemiş-
maruziyet imzası. Gerekçe: MIMIC-IV'te güvenilir yapılandırılmış pompa
itemid'i yok; notes NLP tekrar-üretilebilirliği düşük; EXC-03'ün asıl amacı
"doz kaydı eksik/dış sistem" ve ICD+imza bunu doğrudan hedefler, `pump_class`
A–E ile denetlenebilir kılar.

### v1.0 (15.07.2026): ilk sürüm

---

## Süreç notu
v1.0 ve v1.1 tanımları MIMIC-IV yapısı görülmeden yazıldı; v1.2 düzeltmesi sonuç analizlerinden
önce yapıldı. Kriter tanımları v1.3'ten (16.07.2026) sonra değişmedi. Registry o tarihte resmi
imzayla dondurulmadı; 27.11.2026'da, makale gönderiminden önce son sürüm olarak kapatıldı.
Bu tarihten sonraki her değişiklik protokol değişikliği sayılır, ayrı sürüm notuyla yazılır ve
makalede beyan edilir.

Registry'nin sahibi **S. Efe**'dir. Analiz kodunun yazımında yapay zekâ araçları kullanıldı
(makalede beyan edildi); bu araçlar registry'yi değiştiremez, yalnız öneri yapar.

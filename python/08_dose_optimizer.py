#!/usr/bin/env python3
"""
08 — KISITLI DOZ OPTIMIZASYONU                              v2.0
================================================================
ML doz TAHMIN ETMEZ (Kural 4). ML doz->glukoz haritasini verir; dozu
BURASI secer. Boylece behavior cloning'den kacinilir ve karar kurali
DENETLENEBILIR olur.

Amac fonksiyonu -- ASIMETRIK (NICE-SUGAR NEJM 2009):
  siki kontrol (81-108) mortaliteyi ARTIRDI, mekanizma hipoglisemi.
  -> hedef bant 140-180; hipoglisemiden kacinma BIRINCIL, dusurme ikincil.
  -> 200 mg/dL'nin cezasi, 65 mg/dL'nin cezasindan COK daha kucuk.

Karar kurali (sirayla):
  4. guvenlik freni (modelden ONCE, sert kural)
  1. common support disi -> 'insufficient_data', oneri YOK
  2. P(glu<70 | 0-4h) > tau -> doz reddedilir (sert kisit)
  3. kalan dozlar icinde beklenen kaybi minimize eden EN DUSUK doz
"""
import numpy as np
from dataclasses import dataclass

TARGET_LO, TARGET_HI = 140.0, 180.0
HYPO_TAU    = 0.02     # kabul edilen maks P(glu<70). Kalibrasyonla ayarla.
HYPO_TAU_L2 = 0.005    # P(glu<54)


def loss_glucose(g: np.ndarray) -> np.ndarray:
    """Asimetrik, parcali. Birim mg/dL. Hipoglisemi tarafi ~kubik."""
    l = np.zeros_like(g, dtype=float)
    hi = g > TARGET_HI
    l[hi] = ((g[hi] - TARGET_HI) / 40.0) ** 2          # yumusak
    vhi = g > 250
    l[vhi] += ((g[vhi] - 250) / 50.0) ** 2             # >250 hizlanir
    lo = g < TARGET_LO
    l[lo] = ((TARGET_LO - g[lo]) / 40.0) ** 2
    l1 = g < 70
    l[l1] += 25.0 + ((70 - g[l1]) / 5.0) ** 3          # L1 hipoglisemi: DIK
    l2 = g < 54
    l[l2] += 200.0                                      # L2: neredeyse yasak
    return l


@dataclass
class Prediction:
    """Doz->glukoz modelinin ciktisi: nokta DEGIL DAGILIM."""
    dose_grid: np.ndarray        # denenen dozlar (U)
    glu_samples: np.ndarray      # (n_dose, n_mc) t+2h glukoz posterior
    nadir_samples: np.ndarray    # (n_dose, n_mc) 0-4h nadir glukoz
    in_support: np.ndarray       # (n_dose,) bool -- 06'dan gelir


def predict_from_isf(glu_now, isf_samples, cho_g_hr, gamma, dose_grid,
                     sigma, rng, n_mc=2000):
    """07'deki hiyerarsik modelden basit ileri-simulasyon."""
    idx = rng.choice(len(isf_samples), n_mc)
    isf = isf_samples[idx]
    base_drift = gamma["c_cho"] * (cho_g_hr / 10.0) * 10.0
    noise = rng.standard_t(4, size=(len(dose_grid), n_mc)) * sigma
    glu = glu_now + base_drift - dose_grid[:, None] * isf[None, :] + noise
    nadir = glu - np.abs(rng.normal(0.12 * glu, 0.08 * np.abs(glu)))
    return glu, nadir


def choose_dose(pred: Prediction, glu_now: float, cho_interrupted: bool,
                cirrhosis: bool = False):
    n = len(pred.dose_grid)
    p_l1 = (pred.nadir_samples < 70).mean(axis=1)
    p_l2 = (pred.nadir_samples < 54).mean(axis=1)
    exp_loss = loss_glucose(pred.glu_samples.reshape(-1)).reshape(n, -1).mean(axis=1)

    # --- 4. GUVENLIK FRENI (modelden ONCE, sert kural) ---
    if glu_now < 100:
        return dict(dose=0.0, reason="safety_brake_glu<100", exp_loss=None)
    if cho_interrupted and cirrhosis:
        # INT-01: sirroz + beslenme kesintisi = glikojen deposu YOK
        return dict(dose=0.0, reason="safety_brake_cirrhosis_x_NPO", exp_loss=None)
    if cho_interrupted and glu_now < 150:
        return dict(dose=0.0, reason="safety_brake_NPO", exp_loss=None)

    # --- 1. COMMON SUPPORT (06'nin model_card_support.md ciktisi) ---
    ok = pred.in_support.copy()
    if not ok.any():
        return dict(dose=None, reason="insufficient_data", exp_loss=None)

    # --- 2. SERT HIPOGLISEMI KISITI ---
    ok &= (p_l1 <= HYPO_TAU) & (p_l2 <= HYPO_TAU_L2)
    if not ok.any():
        return dict(dose=0.0, reason="all_doses_violate_hypo_constraint", exp_loss=None)

    # --- 3. MIN KAYIP; berabere kalirsa EN DUSUK doz (tutucu) ---
    masked = np.where(ok, exp_loss, np.inf)
    best = masked.min()
    cand = np.where(masked <= best * 1.02)[0]      # %2 tolerans
    j = cand[0]                                     # dose_grid artan sirali
    return dict(dose=float(pred.dose_grid[j]), reason="optimized",
                exp_loss=float(exp_loss[j]),
                p_hypo_l1=float(p_l1[j]), p_hypo_l2=float(p_l2[j]),
                glu_pred_median=float(np.median(pred.glu_samples[j])),
                glu_pred_ci=(float(np.percentile(pred.glu_samples[j], 2.5)),
                             float(np.percentile(pred.glu_samples[j], 97.5))))


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    dose_grid = np.arange(0, 12.5, 0.5)
    senaryolar = [
        ("Normal renal/hepatik, EN devam, glu 220",
         dict(glu_now=220, isf_mu=35, isf_sd=10, cho=12, interrupted=False, cirr=False)),
        ("eGFR 12 (ISF yuksek), EN devam, glu 220",
         dict(glu_now=220, isf_mu=58, isf_sd=22, cho=12, interrupted=False, cirr=False)),
        ("Child-Pugh C sirroz, EN KESILDI, glu 190",
         dict(glu_now=190, isf_mu=45, isf_sd=25, cho=0, interrupted=True, cirr=True)),
        ("Steroid, TPN, glu 280",
         dict(glu_now=280, isf_mu=25, isf_sd=8, cho=20, interrupted=False, cirr=False)),
    ]
    for ad, p in senaryolar:
        isf_s = np.exp(rng.normal(np.log(p["isf_mu"]), p["isf_sd"]/p["isf_mu"], 4000))
        glu, nad = predict_from_isf(p["glu_now"], isf_s, p["cho"],
                                    {"c_cho": 0.9}, dose_grid, 22.0, rng)
        pred = Prediction(dose_grid, glu, nad,
                          in_support=dose_grid <= (8 if p["isf_mu"] > 50 else 12))
        r = choose_dose(pred, p["glu_now"], p["interrupted"], p["cirr"])
        print(f"\n{ad}\n  -> doz: {r['dose']} U   ({r['reason']})")
        if r.get("p_hypo_l1") is not None:
            print(f"     P(glu<70)={r['p_hypo_l1']:.3f}  tahmin t+2h={r['glu_pred_median']:.0f} "
                  f"[{r['glu_pred_ci'][0]:.0f}-{r['glu_pred_ci'][1]:.0f}]")

#!/usr/bin/env python3
"""
09 — AIPW SUREKLI-DOZ YANIT EGRISI (Makale 2, ADIM 1)              v1.0
======================================================================
Hedef deneme emulasyonu: bazal glarjin doz_u_kg -> gece glukoz nadiri.
Cift-saglam (AIPW/DR) surekli-tedavi tahmini (Kennedy 2017 tarzi kernel-DR),
2-kat CAPRAZ-UYDURMA (DML) ile. eGFR bandina gore AYRI egri.

  GPS (doz yogunlugu f(D|X)): ILK-DOZ ve TITRASYON kollari AYRI (atama
    mekanizmasi farkli; titrasyonda onceki doz belirleyici).
  Sonuc modeli mu(D,X): tum glarjin (LightGBM).
  AIPW(d) = mean_i[ mu(d,X_i) + Kh(D_i-d)/f(D_i|X_i) * (Y_i - mu(D_i,X_i)) ]
    agirlik budandi (truncation) -- pozitivite/ust-kuyruk stabilizasyonu.

Y = glu_nadir_night (doz+4..12h min glukoz). Dusuk = hipoglisemi riski.
KULLANIM: python python/09_aipw_basal.py --parquet 'data/basal_events/*.parquet'
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightgbm as lgb

X_CONF = ["glu_now", "glu_mean_24h", "glu_sd_24h", "glu_min_24h", "weight_kg",
          "egfr_now", "on_crrt", "kdigo_stage", "bun", "pred_eq_24h", "on_steroid",
          "on_vasopressor", "cirrhosis", "child_pugh_score", "dm_t1", "dm_any",
          "cho_g_hr", "cho_24h_cum", "iob_bolus_u", "dose_24h_cum", "age", "event_hr"]
X_TITR = ["prev_glargine_dose_u_kg", "hrs_since_prev_glargine"]  # titrasyon ek
DOSE_GRID = np.round(np.arange(0.05, 0.55, 0.05), 3)   # U/kg destekli aralik
BW = 0.06          # kernel bant genisligi (U/kg)
WT_TRUNC = 20.0    # agirlik ust budama


def gauss(x, m, s):
    s = np.clip(s, 1e-6, None)
    return np.exp(-0.5 * ((x - m) / s) ** 2) / (s * np.sqrt(2 * np.pi))


def boot_ci(subj, vals, B=500, seed=0):
    """SUBJECT duzeyinde (cluster) bootstrap: E[Y(d)] band-ortalamasi icin %95 CI.
    Hasta bazli yeniden ornekleme (olay bagimliligini hesaba katar)."""
    g = pd.DataFrame({"s": subj, "v": vals}).groupby("s")["v"].agg(["sum", "count"])
    S, C = g["sum"].values, g["count"].values
    if len(S) < 5:
        m = float(vals.mean()) if len(vals) else float("nan")
        return m, m
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(S), size=(B, len(S)))
    bs = S[idx].sum(axis=1) / np.clip(C[idx].sum(axis=1), 1, None)
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def load(parquet):
    df = duckdb.connect().execute(f"""
        SELECT * FROM read_parquet('{parquet}')
        WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL
          AND glu_nadir_night IS NOT NULL AND weight_kg IS NOT NULL
    """).df()
    for b in ("cirrhosis", "dm_t1", "dm_any", "on_crrt", "on_vasopressor", "on_steroid"):
        if b in df: df[b] = df[b].astype("float")
    # aykiri kontrol (gorev-1): dose_u_kg p99.5 kirp, dose_24h_cum p99 winsorize
    df = df[df.dose_u_kg <= df.dose_u_kg.quantile(0.995)].reset_index(drop=True)
    cap = df.dose_24h_cum.quantile(0.99)
    df["dose_24h_cum"] = df.dose_24h_cum.clip(upper=cap)
    df["is_first"] = df.prev_glargine_dose_u_kg.isna()
    df["egfr_band"] = pd.cut(df.egfr_now, [-1, 15, 45, 90, 1e9],
                             labels=["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"])
    return df


def fit_gps_density(dtr, dte, arm_col):
    """Kol-spesifik GPS: f(D|X). dtr'de fit, dte'de her satirin f(D_i|X_i)."""
    fte = np.full(len(dte), np.nan)
    for arm in (True, False):
        tr = dtr[dtr[arm_col] == arm]; te_mask = (dte[arm_col] == arm).values
        if tr[X_CONF].shape[0] < 30 or te_mask.sum() == 0:
            continue
        m = lgb.LGBMRegressor(n_estimators=300, learning_rate=0.05, num_leaves=15,
                              subsample=0.8, random_state=0, verbose=-1)
        m.fit(tr[X_CONF], tr.dose_u_kg)
        resid_sd = float(np.std(tr.dose_u_kg - m.predict(tr[X_CONF])) + 1e-6)
        pred_te = m.predict(dte[X_CONF][te_mask])
        fte[te_mask] = gauss(dte.dose_u_kg.values[te_mask], pred_te, resid_sd)
    return fte


def fit_outcome(dtr):
    feats = X_CONF + X_TITR + ["dose_u_kg"]
    m = lgb.LGBMRegressor(n_estimators=500, learning_rate=0.03, num_leaves=31,
                          subsample=0.8, colsample_bytree=0.8, random_state=0, verbose=-1)
    m.fit(dtr[feats], dtr.glu_nadir_night)
    return m, feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    df = load(a.parquet)
    print(f"[AIPW] n={len(df):,} glarjin olay (gece-glukozlu), {df.subject_id.nunique()} hasta; "
          f"ilk-doz={int(df.is_first.sum())}, titrasyon={int((~df.is_first).sum())}")

    # ---- 2-kat capraz-uydurma (subject bazli) ----
    rng = np.random.default_rng(0)
    subs = df.subject_id.unique()
    fold_of = {s: int(i % 2) for i, s in enumerate(rng.permutation(subs))}
    df["fold"] = df.subject_id.map(fold_of)

    psi = {d: np.full(len(df), np.nan) for d in DOSE_GRID}   # AIPW pseudo per doz
    for k in (0, 1):
        tr = df[df.fold != k]; te = df[df.fold == k]; idx = te.index.values
        fD = fit_gps_density(tr, te, "is_first")             # f(D_i|X_i) kol-spesifik
        mo, feats = fit_outcome(tr)
        muD = mo.predict(te[feats])                          # mu(D_i,X_i)
        for d in DOSE_GRID:
            Xte = te.copy(); Xte["dose_u_kg"] = d
            mud = mo.predict(Xte[feats])                     # mu(d,X_i)
            Kh = gauss(te.dose_u_kg.values, d, BW)
            w = np.where(fD > 1e-8, Kh / fD, 0.0)
            w = np.clip(w, 0, WT_TRUNC)
            psi[d][idx] = mud + w * (te.glu_nadir_night.values - muD)

    # ---- eGFR bandina gore egri ----
    curves = {}
    subj_all = df.subject_id.values
    for band in ["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]:
        mask = (df.egfr_band == band).values
        pts = []
        for d in DOSE_GRID:
            ok = mask & ~np.isnan(psi[d])
            v = psi[d][ok]; s = subj_all[ok]
            lo, hi = boot_ci(s, v)
            pts.append(dict(dose_u_kg=float(d), n=int(len(v)),
                            E_nadir=round(float(np.mean(v)), 1),
                            ci_lo=round(lo, 1), ci_hi=round(hi, 1)))
        curves[band] = dict(n_events=int(mask.sum()),
                            n_subjects=int(df.subject_id[mask].nunique()), curve=pts)
        print(f"[AIPW] {band:12} n={int(mask.sum()):5}  "
              + "  ".join(f"{p['dose_u_kg']:.2f}:{p['E_nadir']:.0f}[{p['ci_lo']:.0f},{p['ci_hi']:.0f}]"
                          for p in pts[::3]))

    res = dict(outcome="glu_nadir_night (doz+4..12h min glukoz, mg/dL)",
               method="AIPW/DR surekli-doz, 2-kat capraz-uydurma, kol-spesifik GPS; "
                      "CI = SUBJECT-duzeyi cluster bootstrap (500x, %95 persentil)",
               dose_grid_u_kg=DOSE_GRID.tolist(), bandwidth=BW, weight_trunc=WT_TRUNC,
               n_first=int(df.is_first.sum()), n_titration=int((~df.is_first).sum()),
               curves_by_egfr=curves)
    (out / "07_aipw_curve.json").write_text(json.dumps(res, indent=2))

    # ---- figur ----
    fig, ax = plt.subplots(figsize=(9, 6))
    for band, col in zip(["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"],
                         ["#c33", "#e90", "#3a3", "#39c"]):
        c = curves[band]["curve"]
        d = [p["dose_u_kg"] for p in c]; e = [p["E_nadir"] for p in c]
        elo = [p["E_nadir"] - p["ci_lo"] for p in c]; ehi = [p["ci_hi"] - p["E_nadir"] for p in c]
        ax.errorbar(d, e, yerr=[elo, ehi], marker="o", lw=2, color=col, capsize=2,
                    label=f"{band} (n={curves[band]['n_events']})")
    ax.axhline(70, color="k", ls="--", lw=1, label="hipoglisemi (70)")
    ax.set_xlabel("bazal glarjin doz (U/kg)"); ax.set_ylabel("E[gece glukoz nadiri] (mg/dL)")
    ax.set_title("AIPW doz-yanit: bazal doz -> gece nadiri, eGFR bandina gore")
    ax.legend(); fig.tight_layout(); fig.savefig(out / "07_aipw_curve.png", dpi=120); plt.close(fig)
    print(f"\nCiktilar -> {out.resolve()}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
07-GPS — GENELLESTIRILMIS EGILIM SKORU TANISI (Makale 2 fizibilite)
====================================================================
p(doz_u_kg | X) LightGBM ile. Model KURMAZ; karistirma yapisini TANIR.
  (a) CV R2: doz ne kadar ongorulebilir?  >0.6 guclu karistirma |
      <0.3 atama ~rastgele -> nedensel iddia GUCLENIR.
  (b) SHAP: dozu en cok ne belirliyor? egfr_now listenin neresinde?
  (c) Ortusme: tahmini GPS dagilimi eGFR bantlari arasinda ortusuyor mu?
Ciktilar: qc/07_gps_diagnostic.json, qc/07_gps_shap.png
KULLANIM: python python/07_gps_diagnostic.py --parquet 'data/basal_events/*.parquet'
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightgbm as lgb
import shap
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score

# Karistirici X (doz aninda). Hedefler / id / split / doz DISLANIR.
CONFOUNDERS = [
    "egfr_now", "on_crrt", "kdigo_stage", "bun", "cr_last", "cr_age_hr",
    "weight_kg", "pred_eq_24h", "on_steroid", "on_vasopressor",
    "cirrhosis", "child_pugh_score", "meld_na", "dm_t1", "dm_any",
    "glu_now", "glu_mean_24h", "glu_sd_24h", "glu_min_24h",
    "cho_g_hr", "cho_24h_cum", "dose_24h_cum", "iob_bolus_u",
    "prev_glargine_dose_u_kg", "hrs_since_prev_glargine", "age", "event_hr",
]
INERTIA = ["prev_glargine_dose_u_kg", "hrs_since_prev_glargine"]  # doz atalet terimleri


def load(parquet, target, where=None):
    conds = [f"{target} IS NOT NULL"]
    if where:
        conds.append(f"({where})")
    df = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{parquet}') WHERE {' AND '.join(conds)}").df()
    for b in ("cirrhosis", "dm_t1", "dm_any", "on_crrt", "on_vasopressor", "on_steroid"):
        if b in df:
            df[b] = df[b].astype("float")
    return df


def cv_r2(X, y, groups, seed=42):
    gkf = GroupKFold(n_splits=5)
    preds = np.zeros(len(y)); scores = []
    for tr, te in gkf.split(X, y, groups):
        m = lgb.LGBMRegressor(n_estimators=400, learning_rate=0.03,
                              num_leaves=31, subsample=0.8, colsample_bytree=0.8,
                              random_state=seed, verbose=-1)
        m.fit(X.iloc[tr], y.iloc[tr])
        p = m.predict(X.iloc[te]); preds[te] = p
        scores.append(r2_score(y.iloc[te], p))
    return float(np.mean(scores)), float(np.std(scores)), preds


def hist_overlap(a, b, lo, hi, bins=30):
    """Iki dagilimin histogram ortusme katsayisi (0-1)."""
    pa, edges = np.histogram(a, bins=bins, range=(lo, hi), density=True)
    pb, _ = np.histogram(b, bins=bins, range=(lo, hi), density=True)
    w = edges[1] - edges[0]
    return float(np.sum(np.minimum(pa, pb) * w))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    ap.add_argument("--target", default="dose_u_kg", help="tedavi (doz) sutunu")
    ap.add_argument("--confounders", default=None, help="virgullu; yoksa basal preset")
    ap.add_argument("--where", default=None, help="altkume filtresi (SQL WHERE parcasi)")
    ap.add_argument("--out-name", default="07_gps_diagnostic", help="cikti dosya taban adi")
    ap.add_argument("--log1p", action="store_true", help="hedefe log1p uygula")
    ap.add_argument("--clip-q", type=float, default=None,
                    help="hedefin bu quantile ustunu KIRP (satir dislar; or. 0.995)")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    conf = [c.strip() for c in a.confounders.split(",")] if a.confounders else CONFOUNDERS

    df = load(a.parquet, a.target, a.where)
    n0 = len(df)
    if a.clip_q:                          # hedef ust-kuyruk aykirilarini disla
        thr = df[a.target].quantile(a.clip_q)
        df = df[df[a.target] <= thr].reset_index(drop=True)
    feats = [c for c in conf if c in df.columns]
    X = df[feats].astype("float32")
    y = df[a.target].astype("float")
    if a.log1p:
        y = np.log1p(y)
    groups = df["subject_id"]
    print(f"[GPS] target={a.target}{'(log1p)' if a.log1p else ''} where={a.where or '-'} "
          f"clip_q={a.clip_q} n={len(df):,}/{n0:,} satir, "
          f"{df.subject_id.nunique():,} hasta, {len(feats)} kovaryat")

    # (a) CV R2 -- tam + (varsa) atalet HARIC
    r2_full, sd_full, gps_pred = cv_r2(X, y, groups)
    feats_np = [c for c in feats if c not in INERTIA]
    if feats_np != feats:
        r2_np, _, _ = cv_r2(df[feats_np].astype("float32"), y, groups)
    else:
        r2_np = None
    print(f"[GPS] CV R2 tam = {r2_full:.3f} +-{sd_full:.3f}"
          + (f" | atalet HARIC = {r2_np:.3f}" if r2_np is not None else ""))

    # (b) SHAP -- tum veride fit
    m = lgb.LGBMRegressor(n_estimators=400, learning_rate=0.03, num_leaves=31,
                          subsample=0.8, colsample_bytree=0.8, random_state=42, verbose=-1)
    m.fit(X, y)
    sv = shap.TreeExplainer(m).shap_values(X)
    mabs = np.abs(sv).mean(axis=0)
    rank = pd.DataFrame({"feature": feats, "mean_abs_shap": mabs}) \
        .sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    egfr_rank = int(rank.index[rank.feature == "egfr_now"][0]) + 1 if "egfr_now" in feats else None
    print("[GPS] SHAP siralama (ilk 8):")
    for i, r in rank.head(8).iterrows():
        print(f"   {i+1:2d}. {r.feature:<26} {r.mean_abs_shap:.4f}"
              + ("   <-- egfr_now" if r.feature == "egfr_now" else ""))

    # (c) Ortusme: tahmini GPS (gps_pred) eGFR bantlari arasinda
    band = pd.cut(df.egfr_now, [-1, 15, 45, 90, 1e9],
                  labels=["<15", "15-45", "45-90", ">=90"])
    lo, hi = float(np.nanpercentile(gps_pred, 1)), float(np.nanpercentile(gps_pred, 99))
    bands = ["<15", "15-45", "45-90", ">=90"]
    ov = {}
    for i in range(len(bands)):
        for j in range(i + 1, len(bands)):
            gi = gps_pred[(band == bands[i]).values]; gj = gps_pred[(band == bands[j]).values]
            if len(gi) > 5 and len(gj) > 5:
                ov[f"{bands[i]}|{bands[j]}"] = round(hist_overlap(gi, gj, lo, hi), 3)
    min_ov = min(ov.values()) if ov else None
    print(f"[GPS] GPS ortusme (min ikili) = {min_ov}")

    diag = dict(
        n_events=int(len(df)), n_subjects=int(df.subject_id.nunique()),
        n_features=len(feats),
        target=a.target, where=a.where,
        cv_r2_full=round(r2_full, 4), cv_r2_full_sd=round(sd_full, 4),
        cv_r2_no_inertia=(round(r2_np, 4) if r2_np is not None else None),
        r2_interpretation=("guclu_karistirma(>0.6)" if r2_full > 0.6
                           else "orta(0.3-0.6)" if r2_full >= 0.3 else "zayif~rastgele(<0.3)"),
        egfr_now_rank=egfr_rank,
        egfr_now_mean_abs_shap=(round(float(rank.loc[rank.feature == "egfr_now",
                                "mean_abs_shap"].iloc[0]), 4) if egfr_rank else None),
        shap_ranking=[dict(rank=i + 1, feature=r.feature,
                           mean_abs_shap=round(float(r.mean_abs_shap), 4))
                      for i, r in rank.iterrows()],
        gps_overlap_pairwise=ov, gps_overlap_min=min_ov,
    )
    (out / f"{a.out_name}.json").write_text(json.dumps(diag, indent=2))

    # ---- figur ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    top = rank.head(15)[::-1]
    ax[0].barh(top.feature, top.mean_abs_shap, color="#39c")
    if "egfr_now" in top.feature.values:
        yi = list(top.feature).index("egfr_now")
        ax[0].get_yticklabels()[yi].set_color("#c33")
        ax[0].get_yticklabels()[yi].set_fontweight("bold")
    ax[0].set_title(f"SHAP onem ({a.target}); CV R2={r2_full:.2f}")
    ax[0].set_xlabel("ortalama |SHAP|")
    for bnd, col in zip(bands, ["#c33", "#e90", "#3a3", "#39c"]):
        g = gps_pred[(band == bnd).values]
        if len(g) > 5:
            ax[1].hist(g, bins=30, range=(lo, hi), density=True, histtype="step",
                       lw=2, color=col, label=f"eGFR {bnd} (n={len(g)})")
    ax[1].set_title(f"GPS (tahmini doz_u_kg) ortusme; min={min_ov}")
    ax[1].set_xlabel("tahmini doz_u_kg"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / f"{a.out_name}_shap.png", dpi=120); plt.close(fig)
    print(f"\nCiktilar -> {out.resolve()}")


if __name__ == "__main__":
    main()

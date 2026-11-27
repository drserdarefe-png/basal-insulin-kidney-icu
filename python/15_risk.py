#!/usr/bin/env python3
"""
15 — RISK TABLOSU + BETIMLEYICI RISK SKORU (Makale 2). Okuyucuya verilen sey.
Kimliklendirilemez ama BETIMLENEBILIR: kalibre bir hipo OLASILIGI.
KRITIK: TAHMIN modeli, NEDENSEL DEGIL. Katsayilar 'etki' diye YORUMLANMAZ.
Cikti: qc/paper2_tables/T15_risk_table.json/.md, T16_reader_guide.md,
       qc/paper2_figs/F5_calibration.png
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from statsmodels.stats.proportion import proportion_confint
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss
from sklearn.calibration import calibration_curve
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIN_SUBJ = 20


def load(parquet):
    df = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{parquet}') WHERE basal_type='glargine' "
        "AND dose_u_kg IS NOT NULL AND glu_nadir_night IS NOT NULL").df()
    for b in ("dm_t1", "dm_any", "cirrhosis"):
        df[b] = df[b].astype(float)
    df["dm_status"] = np.where(df.dm_t1 == 1, "T1DM", np.where(df.dm_any == 1, "T2DM", "stres"))
    df["eb"] = pd.cut(df.egfr_now, [-1, 15, 45, 90, 1e9],
                      labels=["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]).astype(str)
    df["dtert"] = pd.qcut(df.dose_u_kg, 3, labels=["D1", "D2", "D3"]).astype(str)
    df["snight"] = df.steroid_night_active.fillna(0).astype(int)
    return df


def rate_ci(d):
    n = len(d); k = int(d.hypo_night_l1.sum()); ns = d.subject_id.nunique()
    if ns < MIN_SUBJ:
        return dict(n_events=n, n_subjects=ns, rate="yetersiz veri (<20 hasta)")
    lo, hi = proportion_confint(k, n, method="wilson")
    return dict(n_events=n, n_subjects=ns, hypo_rate=round(k / n, 4),
                ci95=[round(lo, 4), round(hi, 4)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True); ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    T = Path(a.outdir) / "paper2_tables"; F = Path(a.outdir) / "paper2_figs"
    T.mkdir(parents=True, exist_ok=True); F.mkdir(parents=True, exist_ok=True)
    df = load(a.parquet)

    # ---- T15 MARJINAL risk tablolari (doz tertiline gore stratifiye) ----
    dims = {"egfr_band": "eb", "cirrhosis": "cirrhosis", "steroid_night_active": "snight",
            "dm_status": "dm_status"}
    t15 = {"NOT": "Gece hipoglisemisi (hypo_night_l1) orani + Wilson %95 CI. Hucre esigi "
                  ">=20 hasta (Kural 10). Tam capraz cok seyrek -> MARJINAL, doz tertiline gore.",
           "overall": rate_ci(df), "by_dose_tertile": {dt: rate_ci(df[df.dtert == dt]) for dt in ["D1", "D2", "D3"]},
           "marginals_x_dose_tertile": {}}
    for nm, col in dims.items():
        t15["marginals_x_dose_tertile"][nm] = {}
        for dt in ["D1", "D2", "D3"]:
            sub = df[df.dtert == dt]
            t15["marginals_x_dose_tertile"][nm][dt] = {
                str(v): rate_ci(sub[sub[col] == v])
                for v in sorted(sub[col].dropna().unique(), key=str)}
    (T / "T15_risk_table.json").write_text(json.dumps(t15, indent=2))
    # md ozeti
    md = ["# T15 — Gece hipoglisemi risk tablosu (betimleyici)\n",
          f"Genel: {t15['overall']['hypo_rate']*100:.1f}% "
          f"[{t15['overall']['ci95'][0]*100:.1f}-{t15['overall']['ci95'][1]*100:.1f}]  "
          f"(n={t15['overall']['n_events']}, {t15['overall']['n_subjects']} hasta)\n",
          "## Doz tertiline gore"]
    for dt in ["D1", "D2", "D3"]:
        r = t15["by_dose_tertile"][dt]
        md.append(f"- {dt}: {r.get('hypo_rate', r.get('rate'))}")
    (T / "T15_risk_table.md").write_text("\n".join(md), encoding="utf-8")

    # ---- Betimleyici RISK SKORU (F5) -- TAHMIN, nedensel DEGIL ----
    feat_num = ["dose_u_kg", "glu_mean_24h", "weight_kg", "prev_glargine_dose_u_kg"]
    X = df.copy()
    X["prev_glargine_dose_u_kg"] = X.prev_glargine_dose_u_kg.fillna(0.0)
    for c in ["glu_mean_24h", "weight_kg"]:
        X[c] = X[c].fillna(X[c].median())
    Xd = pd.get_dummies(X[["eb", "dm_status"]].astype(str), drop_first=True)
    Xall = pd.concat([X[feat_num + ["cirrhosis", "snight"]].reset_index(drop=True),
                      Xd.reset_index(drop=True)], axis=1).astype(float)
    y = X.hypo_night_l1.values
    tr = (X.split == "train").values; te = (X.split.isin(["valid", "test"])).values
    m = LogisticRegression(max_iter=1000)   # kalibrasyon icin dengeleme YOK
    m.fit(Xall[tr], y[tr])
    p = m.predict_proba(Xall[te])[:, 1]
    prauc = float(average_precision_score(y[te], p)); brier = float(brier_score_loss(y[te], p))
    frac, mean = calibration_curve(y[te], p, n_bins=8, strategy="quantile")
    risk = dict(model="lojistik (TAHMIN, NEDENSEL DEGIL -- katsayi yorumlanmaz)",
                features=list(Xall.columns), split="subject bucket (model_set ile ayni)",
                n_train=int(tr.sum()), n_test=int(te.sum()),
                base_rate_test=round(float(y[te].mean()), 4),
                PR_AUC=round(prauc, 3), Brier=round(brier, 4))
    (T / "T15b_risk_score.json").write_text(json.dumps(risk, indent=2))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].plot([0, max(mean)*1.1], [0, max(mean)*1.1], "k--", label="ideal")
    ax[0].plot(mean, frac, "o-", color="0.2")
    ax[0].set_xlabel("tahmini olasilik"); ax[0].set_ylabel("gozlenen oran")
    ax[0].set_title(f"F5 Kalibrasyon (Brier={brier:.3f})"); ax[0].legend()
    from sklearn.metrics import precision_recall_curve
    pr, rc, _ = precision_recall_curve(y[te], p)
    ax[1].plot(rc, pr, color="0.2"); ax[1].axhline(y[te].mean(), color="k", ls=":", label="baz oran")
    ax[1].set_xlabel("recall"); ax[1].set_ylabel("precision")
    ax[1].set_title(f"F5 PR egrisi (PR-AUC={prauc:.3f})"); ax[1].legend()
    fig.tight_layout(); fig.savefig(F / "F5_calibration.png", dpi=300); plt.close(fig)

    # ---- T16 okuyucu kilavuzu ----
    g = ["# T16 — Okuyucuya ne diyoruz, ne demiyoruz\n",
         "## DIYORUZ (betimleyici)",
         "- Klinisyen bazal glarjin dozunu **eGFR'ye gore ayarlamiyor** (eGFR<15'te 2.27x Aronoff, "
         "KW p=0.76, SHAP rank dusuk). Sirroz (SHAP 24/24) ve steroid (23/24) icin de degil.",
         "- Buna karsin **agirliga** (units_total 12->18->20, p<0.001) ve **DM fenotipine** "
         "(T1/T2/stres p<0.001) gore ayarliyorlar -> ayarlayabiliyorlar, komorbidite icin degil.",
         "- Su profillerde gece hipoglisemisi **daha sik goruluyor** (T15): eGFR<15, yuksek doz "
         "tertili, steroid-sonmus gece.",
         f"- Bu risk **kalibre bir olasilikla** ifade edilebilir (F5: PR-AUC={prauc:.2f}, Brier={brier:.3f}).\n",
         "## DEMIYORUZ (cunku doz-sonuc kimliklendirilemiyor)",
         "- 'Optimal doz X'tir' — HAYIR.",
         "- 'Doz azaltilmalidir' — HAYIR.",
         "- 'Bu risk doz YUZUNDEN' — HAYIR (bkz. qc/paper2_identification.md; AIPW+within-patient "
         "zamanla-degisen karistirma).\n",
         "## Alan baglantisi",
         "GLUCOSE ('klinisyeni gectik') ve Insulin4RL ayni zamanla-degisen karistirmayi tasir ama "
         "BEYAN ETMEZ. Katkimiz: ayni veri/yontem, ama **sinirin gosterilmesi** + kalibre risk."]
    (T / "T16_reader_guide.md").write_text("\n".join(g), encoding="utf-8")

    print(f"[15] T15 overall hypo={t15['overall']['hypo_rate']}; "
          f"D1/D2/D3={[t15['by_dose_tertile'][d].get('hypo_rate') for d in ['D1','D2','D3']]}")
    print(f"[15] risk skoru: PR-AUC={prauc:.3f} Brier={brier:.4f} base={y[te].mean():.4f} "
          f"(n_test={te.sum()})")


if __name__ == "__main__":
    main()

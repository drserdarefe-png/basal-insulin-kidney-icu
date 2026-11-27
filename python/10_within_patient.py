#!/usr/bin/env python3
"""
10 — HASTA-ICI DELTA-DOZ TASARIMI (Makale 2 BIRINCIL)             v1.0
=====================================================================
KARAR: hastalar-arasi ISF ~100 kat yayilim (sigma_log 1.14) + olculmemis
"gercek insulin ihtiyaci" karistiricisi -> AIPW rezidu-biasli. Ayni hastada
ARDISIK glarjin doz ciftlerini FARKLA -> zaman-degismez hasta karistiricilari
(gercek ihtiyac dahil) DUSER.

Birim: (doz_t, doz_{t-1}) ayni yatis, hrs_since_prev in [18,30] (gunluk doz).
  d_dose_u_kg   = dose_t - dose_{t-1}
  d_nadir_night = nadir_t - nadir_{t-1}
MODEL A (birincil): d_nadir ~ d_dose + d_kovaryatlar; cluster-robust SE (subject).
MODEL B: A + prev_night_nadir (zamanla-degisen karistirma / ortalamaya donus).
MODEL C: d_dose x egfr_band (INT-02 bazal karsiligi).
Beklenen beta_d_dose NEGATIF (doz artinca gece nadiri duser).
DUR: |d_dose|>0.05 cift < 500 -> varyasyon yok.
KULLANIM: python python/10_within_patient.py --parquet 'data/basal_events/*.parquet'
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import statsmodels.formula.api as smf

LAGCOLS = ["dose_u_kg", "glu_nadir_night", "cho_24h_cum", "pred_eq_active",
           "glu_mean_24h", "weight_kg", "egfr_now", "on_crrt", "iob_bolus_u", "glu_now"]
RHS = "d_dose + d_cho24 + d_pred_eq + d_glu_mean + d_weight + d_egfr + d_on_crrt + d_iob"


def build_pairs(parquet):
    df = duckdb.connect().execute(f"""
        SELECT * FROM read_parquet('{parquet}')
        WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL
    """).df()
    for b in ("on_crrt",):
        df[b] = df[b].astype("float")
    df = df.sort_values(["stay_id", "dose_time"]).reset_index(drop=True)
    g = df.groupby("stay_id")
    prev = g[LAGCOLS].shift(1)
    d = pd.DataFrame({"subject_id": df.subject_id, "stay_id": df.stay_id,
                      "hrs_since_prev": df.hrs_since_prev_glargine,
                      "egfr_now": df.egfr_now, "prev_night_nadir": prev.glu_nadir_night})
    d["d_dose"]     = df.dose_u_kg      - prev.dose_u_kg
    d["d_nadir"]    = df.glu_nadir_night - prev.glu_nadir_night
    d["d_glu_now"]  = df.glu_now        - prev.glu_now
    d["d_cho24"]    = df.cho_24h_cum    - prev.cho_24h_cum
    d["d_pred_eq"]  = df.pred_eq_active - prev.pred_eq_active
    d["d_glu_mean"] = df.glu_mean_24h   - prev.glu_mean_24h
    d["d_weight"]   = df.weight_kg      - prev.weight_kg
    d["d_egfr"]     = df.egfr_now       - prev.egfr_now
    d["d_on_crrt"]  = df.on_crrt        - prev.on_crrt
    d["d_iob"]      = df.iob_bolus_u    - prev.iob_bolus_u
    # gunluk doz araligi
    d = d[(d.hrs_since_prev >= 18) & (d.hrs_since_prev <= 30)]
    d["egfr_band"] = pd.cut(d.egfr_now, [-1, 15, 45, 90, 1e9],
                            labels=["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"])
    return df, d


def fit(dd, formula, y="d_nadir"):
    sub = dd.dropna(subset=[y] + [c for c in
          ["d_dose", "d_cho24", "d_pred_eq", "d_glu_mean", "d_weight", "d_egfr",
           "d_on_crrt", "d_iob", "prev_night_nadir"] if c in formula])
    m = smf.ols(f"{y} ~ {formula}", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub.subject_id})
    return m, len(sub)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    df, d = build_pairs(a.parquet)
    # hrs_since_prev toplam dagilim (dislanan araliklar)
    raw = df.hrs_since_prev_glargine.dropna()
    n_all = int((raw.notna()).sum())
    n_daily = len(d)
    dd = d.dropna(subset=["d_dose", "d_nadir"])
    n_var = int((dd.d_dose.abs() > 0.05).sum())
    print(f"[10] cift(gunluk 18-30h)={n_daily} / tum-ardisik~{n_all}; "
          f"nadir-cifti={len(dd)}, hasta={dd.subject_id.nunique()}")
    print(f"[10] |d_dose|>0.05 U/kg cift = {n_var}")

    res = dict(n_pairs_daily=n_daily, n_pairs_with_nadir=int(len(dd)),
               n_subjects=int(dd.subject_id.nunique()),
               n_pairs_dvar_gt005=n_var,
               d_dose_pctl=[round(float(x), 3) for x in np.percentile(dd.d_dose, [5, 25, 50, 75, 95])],
               hrs_since_prev_excluded={
                   "lt18": int((raw < 18).sum()), "in18_30": int(((raw >= 18) & (raw <= 30)).sum()),
                   "gt30": int((raw > 30).sum())})

    if n_var < 500:
        res["STOP"] = f"|d_dose|>0.05 cift={n_var} < 500 -- VARYASYON YOK"
        (out / "10_within_patient.json").write_text(json.dumps(res, indent=2))
        print("[10] STOP:", res["STOP"]); return

    # MODEL A
    mA, nA = fit(dd, RHS)
    # MODEL B
    mB, nB = fit(dd, RHS + " + prev_night_nadir")
    # MODEL C (etkilesim)
    mC, nC = fit(dd, "d_dose * C(egfr_band) + d_cho24 + d_pred_eq + d_glu_mean "
                     "+ d_weight + d_egfr + d_on_crrt + d_iob")

    def bd(m, term="d_dose"):
        ci = m.conf_int().loc[term]
        return dict(beta=round(float(m.params[term]), 2),
                    ci=[round(float(ci[0]), 2), round(float(ci[1]), 2)],
                    p=round(float(m.pvalues[term]), 4))
    res["modelA"] = dict(n=nA, d_dose=bd(mA))
    res["modelB"] = dict(n=nB, d_dose=bd(mB),
                         prev_night_nadir=bd(mB, "prev_night_nadir"))
    # C: per-band d_dose
    bands = ["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]
    ref = bands[0]; perband = {}
    base = float(mC.params["d_dose"])
    for b in bands:
        if b == ref:
            term = "d_dose"; beta = base
        else:
            it = f"d_dose:C(egfr_band)[T.{b}]"
            beta = base + float(mC.params.get(it, 0.0))
            term = it
        perband[b] = round(beta, 2)
    res["modelC_beta_by_band"] = perband

    # NEGATIF KONTROL (IS 3): d_glu_now (doz ANINDAKI = dozdan ONCE olculen glukoz)
    #   d_dose bunu ETKILEYEMEZ. beta anlamli -> rezidu karistirma, delta KIRLI.
    mNC, nNC = fit(dd, RHS, y="d_glu_now")
    res["negative_control_d_glu_now"] = dict(n=nNC, d_dose=bd(mNC),
        verdict=("KIRLI: d_dose->d_glu_now anlamli (rezidu karistirma; delta tasarimi da yanli)"
                 if mNC.pvalues["d_dose"] < 0.05 else
                 "TEMIZ: d_dose->d_glu_now sifira yakin (negatif kontrol gecti)"))
    print(f"[10] NEG-KONTROL d_glu_now: beta_d_dose={res['negative_control_d_glu_now']['d_dose']['beta']} "
          f"CI{res['negative_control_d_glu_now']['d_dose']['ci']} p={res['negative_control_d_glu_now']['d_dose']['p']} "
          f"-> {res['negative_control_d_glu_now']['verdict']}")

    res["interpretation"] = ("beta_d_dose NEGATIF beklenir (doz^ -> gece nadir v). "
                             "MODEL A birincil; B ile uyum zamanla-degisen karistirmayi test eder. "
                             "Negatif kontrol d_glu_now sifira yakin olmali.")
    (out / "10_within_patient.json").write_text(json.dumps(res, indent=2))
    print(f"[10] A beta_d_dose={res['modelA']['d_dose']['beta']} CI{res['modelA']['d_dose']['ci']} "
          f"| B={res['modelB']['d_dose']['beta']} | C bybands={perband}")
    print(f"Ciktilar -> {out.resolve()}")


if __name__ == "__main__":
    main()

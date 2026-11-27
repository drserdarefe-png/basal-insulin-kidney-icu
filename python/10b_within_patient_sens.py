#!/usr/bin/env python3
"""
10b — HASTA-ICI MODEL: DUYARLILIK KOSULARI (27.11.2026)                 v1.0
============================================================================
Ozgun 10_within_patient.py DEGISTIRILMEZ; modul olarak yuklenir (ayni RHS, ayni
cluster-robust SE). Uc kosul, her biri icin Model A ve Model B:
  orig       : ozgun tanimin yeniden uretimi (dogrulama; 10_within_patient.json ile AYNI olmali)
  trim995    : AIPW ile ayni aykiri kurali (dose_u_kg > p99.5 olan olayi iceren ciftler haric)
  fixed_wt   : doz, yatistaki ILK agirliga gore normalize (gunluk agirlik degisimi
               kaynakli d_dose devre disi; unite degismediyse d_dose = 0)
Yalniz TOPLU sonuc yazar: 10b_within_patient_sens.json

KULLANIM:
  python 10b_within_patient_sens.py <ozgun_10_within_patient.py> --parquet 'data/basal_events/*.parquet' --outdir <cikti>
"""
import argparse, importlib.util, json, sys
from pathlib import Path
import numpy as np, pandas as pd, duckdb

orig = sys.argv[1]
spec = importlib.util.spec_from_file_location("wp_orig", orig)
W = importlib.util.module_from_spec(spec); spec.loader.exec_module(W)


def load(parquet):
    df = duckdb.connect().execute(f"""
        SELECT * FROM read_parquet('{parquet}')
        WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL
    """).df()
    df["on_crrt"] = df["on_crrt"].astype("float")
    return df.sort_values(["stay_id", "dose_time"]).reset_index(drop=True)


def pairs(df, mode):
    df = df.copy()
    dose = df.dose_u_kg.copy()
    if mode == "fixed_wt":
        w0 = df.groupby("stay_id").weight_kg.transform("first")      # ilk (bos olmayan) agirlik
        dose = df.units_total / w0
    df["_dose"] = dose
    cols = [c for c in W.LAGCOLS if c != "dose_u_kg"] + ["_dose", "units_total"]
    prev = df.groupby("stay_id")[cols + ["dose_u_kg"]].shift(1)
    d = pd.DataFrame({"subject_id": df.subject_id, "stay_id": df.stay_id,
                      "hrs_since_prev": df.hrs_since_prev_glargine,
                      "egfr_now": df.egfr_now, "prev_night_nadir": prev.glu_nadir_night})
    d["d_dose"]     = df._dose           - prev._dose
    d["d_nadir"]    = df.glu_nadir_night - prev.glu_nadir_night
    d["d_glu_now"]  = df.glu_now         - prev.glu_now
    d["d_cho24"]    = df.cho_24h_cum     - prev.cho_24h_cum
    d["d_pred_eq"]  = df.pred_eq_active  - prev.pred_eq_active
    d["d_glu_mean"] = df.glu_mean_24h    - prev.glu_mean_24h
    d["d_weight"]   = df.weight_kg       - prev.weight_kg
    d["d_egfr"]     = df.egfr_now        - prev.egfr_now
    d["d_on_crrt"]  = df.on_crrt         - prev.on_crrt
    d["d_iob"]      = df.iob_bolus_u     - prev.iob_bolus_u
    d["units_changed"] = (df.units_total != prev.units_total) & prev.units_total.notna()
    keep = (d.hrs_since_prev >= 18) & (d.hrs_since_prev <= 30)
    n_excl_trim = 0
    if mode == "trim995":
        thr = df.dose_u_kg.quantile(0.995)
        bad = (df.dose_u_kg > thr) | (prev.dose_u_kg > thr)
        n_excl_trim = int((keep & bad).sum())
        keep = keep & ~bad
    return d[keep], n_excl_trim


def bd(m, term="d_dose"):
    ci = m.conf_int().loc[term]
    return dict(beta=round(float(m.params[term]), 2),
                ci=[round(float(ci[0]), 2), round(float(ci[1]), 2)],
                p=round(float(m.pvalues[term]), 4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args(sys.argv[2:])
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    df = load(a.parquet)
    res = dict(n_events=int(len(df)), n_stays=int(df.stay_id.nunique()),
               weight=dict(n_stays_weight_varies=int((df.groupby("stay_id").weight_kg.nunique() > 1).sum()),
                           min=round(float(df.weight_kg.min()), 1),
                           n_events_lt30kg=int((df.weight_kg < 30).sum())),
               runs={})
    for mode in ("orig", "trim995", "fixed_wt"):
        d, n_excl = pairs(df, mode)
        dd = d.dropna(subset=["d_dose", "d_nadir"])
        mA, nA = W.fit(dd, W.RHS)
        mB, nB = W.fit(dd, W.RHS + " + prev_night_nadir")
        used = dd.dropna(subset=["d_cho24", "d_pred_eq", "d_glu_mean", "d_weight", "d_egfr", "d_on_crrt", "d_iob"])
        res["runs"][mode] = dict(
            n_pairs_daily=int(len(d)), n_pairs_with_nadir=int(len(dd)),
            n_subjects=int(dd.subject_id.nunique()), n_pairs_excluded_by_trim=n_excl,
            n_pairs_units_changed=int(dd.units_changed.sum()),
            n_pairs_ddose_nonzero=int((dd.d_dose.abs() > 1e-9).sum()),
            d_dose_sd=round(float(dd.d_dose.std()), 4),
            d_dose_abs_max=round(float(dd.d_dose.abs().max()), 3),
            modelA=dict(n=int(nA), n_subjects=int(used.subject_id.nunique()), d_dose=bd(mA)),
            modelB=dict(n=int(nB), d_dose=bd(mB), prev_night_nadir=bd(mB, "prev_night_nadir")))
        r = res["runs"][mode]
        print(f"[10b] {mode:9} cift={r['n_pairs_with_nadir']:5} model_n={nA:5} "
              f"A={r['modelA']['d_dose']['beta']:8.2f} CI{r['modelA']['d_dose']['ci']} | "
              f"B={r['modelB']['d_dose']['beta']:8.2f} CI{r['modelB']['d_dose']['ci']} | "
              f"prev_nadir={r['modelB']['prev_night_nadir']['beta']}")
    (out / "10b_within_patient_sens.json").write_text(json.dumps(res, indent=2))
    print(f"Ciktilar -> {out.resolve()}")


if __name__ == "__main__":
    main()

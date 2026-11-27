#!/usr/bin/env python3
"""
16 — DUYARLILIK (Makale 2, FAZ 3): T3'un saglamligi.
is_first_dose ayrimi, NPH, ckd_stage, cr_baseline_method, 14-gun kirpma etkisi.
Cikti: qc/paper2_tables/T17_sensitivity.json
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy import stats

BANDS = ["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]
# BQ _tmp_sens_ckd.sql sonuclari (u_kg medyan):
CKD_STAGE = {"G1": 0.182, "G2": 0.198, "G3a": 0.185, "G3b": 0.198, "G4": 0.202, "G5": 0.159}
CR_METHOD = {"prior_median": 0.189, "stay_nadir": 0.188, "mdrd_backcalc_n4": 0.402}


def eband(e): return pd.cut(e, [-1, 15, 45, 90, 1e9], labels=BANDS).astype(str)


def dose_by_band(d):
    d = d[d.dose_u_kg.notna()]
    r = {}
    for b in BANDS:
        s = d[d.eb == b]
        r[b] = dict(n=int(len(s)), u_kg_med=round(float(np.median(s.dose_u_kg)), 3) if len(s) else None)
    pat = d.groupby("subject_id").agg(u=("dose_u_kg", "median"), b=("eb", "first"))
    try:
        p = float(stats.kruskal(*[pat[pat.b == b].u.dropna() for b in BANDS if (pat.b == b).any()])[1])
    except Exception:
        p = None
    r["kw_p"] = round(p, 4) if p is not None else None
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True); ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    T = Path(a.outdir) / "paper2_tables"; T.mkdir(parents=True, exist_ok=True)
    df = duckdb.connect().execute(f"SELECT * FROM read_parquet('{a.parquet}')").df()
    df["eb"] = eband(df.egfr_now)
    gl = df[df.basal_type == "glargine"]
    nph = df[df.basal_type == "nph"]

    sens = dict(
        T3_by_is_first_dose=dict(
            first_dose=dose_by_band(gl[gl.prev_glargine_dose_u_kg.isna()]),
            titration=dose_by_band(gl[gl.prev_glargine_dose_u_kg.notna()]),
            note="Ilk doz farkli karar; her ikisinde de eGFR bandlari arasi DUZ mu?"),
        T3_NPH=dict(**dose_by_band(nph),
                    note="NPH (BID, ~487 hasta). Desen glarjinle ayni mi?"),
        T3_ckd_stage_baseline=dict(u_kg_med_by_stage=CKD_STAGE,
                    note="Event eGFR yerine baseline CKD evresi: G1-G5 arasi DUZ (0.16-0.20) "
                         "-> T3 saglam, egfr_unsteady sorununa dayanikli."),
        T3_cr_baseline_method=dict(u_kg_med_by_method=CR_METHOD,
                    note="3 baseline-kreatinin proxy'si: prior_median/stay_nadir ~0.19 (stabil); "
                         "mdrd_backcalc yalniz 4 hasta (ihmal)."),
        crop_14day=dict(
            n_glargine_events=int(len(gl)),
            n_egfr_now_null=int(gl.egfr_now.isna().sum()),
            frac_egfr_null=round(float(gl.egfr_now.isna().mean()), 4),
            note="Makale 2 egfr_now OLAY-ANCHORLU (labevents<=dose_time 48h), GRID DEGIL. "
                 "14-gun grid kirpmasi UYGULANMAZ; egfr kaybi sadece kreatinin yoklugundan "
                 "(<1%). renal_hourly (grid) kullanilsaydi 336h otesi olaylar egfr kaybederdi."),
    )
    (T / "T17_sensitivity.json").write_text(json.dumps(sens, indent=2))

    def row(r): return "/".join(str(r[b]["u_kg_med"]) for b in BANDS) + f" (KW p={r['kw_p']})"
    print("[16] T3 first-dose U/kg by band:", row(sens["T3_by_is_first_dose"]["first_dose"]))
    print("[16] T3 titration:", row(sens["T3_by_is_first_dose"]["titration"]))
    print("[16] T3 NPH:", row(sens["T3_NPH"]))
    print(f"[16] egfr_now null (14-gun kirpma etkisi): {sens['crop_14day']['frac_egfr_null']*100:.1f}% "
          f"(olay-anchorlu -> kirpma yok)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
17 — T3d MEKANIZMA: acigin nedeni. eGFR<15 hastasi daha AGIR + daha DIYABETIK;
klinisyenin agirlik+DM sezgileri kilavuzun TERSI yone itiyor -> mutlak doz
eGFR<15'te DAHA YUKSEK (20 vs 17 U). Aracilik: agirlik+DM kontrol edilince
eGFR-doz isareti donuyor mu? Cikti: qc/paper2_tables/T3d_mechanism.json
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy import stats
import statsmodels.formula.api as smf

BANDS = ["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]


def boot_med_diff(a_s, a_v, b_s, b_v, B=2000, seed=0):
    a_v = np.asarray(a_v, float); b_v = np.asarray(b_v, float)
    a_s = np.asarray(a_s)[np.isfinite(a_v)]; a_v = a_v[np.isfinite(a_v)]
    b_s = np.asarray(b_s)[np.isfinite(b_v)]; b_v = b_v[np.isfinite(b_v)]
    def agg(s, v):
        g = pd.DataFrame({"s": s, "v": v}).groupby("s")["v"].apply(list)
        return g.values
    va, vb = agg(a_s, a_v), agg(b_s, b_v); rng = np.random.default_rng(seed)
    d = [np.median(np.concatenate([va[i] for i in rng.integers(0, len(va), len(va))]))
         - np.median(np.concatenate([vb[i] for i in rng.integers(0, len(vb), len(vb))]))
         for _ in range(B)]
    return float(np.median(a_v) - np.median(b_v)), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True); ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    T = Path(a.outdir) / "paper2_tables"; T.mkdir(parents=True, exist_ok=True)
    df = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{a.parquet}') WHERE basal_type='glargine' "
        "AND units_total IS NOT NULL AND weight_kg IS NOT NULL").df()
    for b in ("dm_t1", "dm_any"):
        df[b] = df[b].astype(float)
    df["eb"] = pd.cut(df.egfr_now, [-1, 15, 45, 90, 1e9], labels=BANDS).astype(str)
    df["dm_status"] = np.where(df.dm_t1 == 1, "T1DM", np.where(df.dm_any == 1, "T2DM", "stres"))

    # (a) mutlak doz farki eGFR<15 vs >=90
    lo = df[df.eb == "egfr<15"]; hi = df[df.eb == "egfr>=90"]
    d, cl, cu = boot_med_diff(lo.subject_id.values, lo.units_total.values,
                              hi.subject_id.values, hi.units_total.values)
    patu = df.groupby("subject_id").agg(u=("units_total", "median"), b=("eb", "first"))
    H, p = stats.kruskal(*[patu[patu.b == b].u.dropna() for b in BANDS if (patu.b == b).any()])

    # (b) aracilik: units_total ~ egfr_band (crude) vs + weight + dm (adjusted)
    dd = df[(df.eb != "nan") & df.units_total.notna() & df.weight_kg.notna()
            & df.dm_status.notna()].copy().reset_index(drop=True)
    dd = dd[np.isfinite(dd.weight_kg) & np.isfinite(dd.units_total)].reset_index(drop=True)
    dd["eb"] = pd.Categorical(dd.eb, categories=["egfr>=90", "egfr_45_90", "egfr_15_45", "egfr<15"])
    # NOT: cluster-SE statsmodels grouputils'te bu veride hata verdi -> HC1 robust SE.
    m1 = smf.ols("units_total ~ C(eb)", dd).fit(cov_type="HC1")
    m2 = smf.ols("units_total ~ C(eb) + weight_kg + C(dm_status)", dd).fit(cov_type="HC1")
    def coef(m, term="C(eb)[T.egfr<15]"):
        ci = m.conf_int().loc[term]
        return dict(beta=round(float(m.params[term]), 2),
                    ci=[round(float(ci[0]), 2), round(float(ci[1]), 2)],
                    p=round(float(m.pvalues[term]), 4))
    c1, c2 = coef(m1), coef(m2)

    res = dict(
        note="MEKANIZMA: eGFR<15 daha agir+diyabetik; klinisyenin agirlik+DM sezgileri "
             "kilavuzun TERSI yone itiyor. 'Bakmiyorlar' degil, 'baktiklari ters itiyor'.",
        absolute_units_egfr_lt15_med=round(float(np.median(lo.units_total)), 1),
        absolute_units_egfr_ge90_med=round(float(np.median(hi.units_total)), 1),
        median_diff_lt15_minus_ge90=round(d, 2), ci95=[round(cl, 2), round(cu, 2)],
        kw_units_by_band_p=round(float(p), 4),
        mediation=dict(
            crude_units_egfr_lt15_coef=c1,
            adjusted_units_egfr_lt15_coef_given_weight_dm=c2,
            interpretation=("Crude eGFR<15 katsayisi %s, agirlik+DM kontrolunde %s. %s"
                            % (("+" if c1["beta"] > 0 else ""), ("+" if c2["beta"] > 0 else ""),
                               "ISARET DONDU -> agirlik/DM gercek renal-azaltmayi MASKELIYOR"
                               if (c1["beta"] > 0 and c2["beta"] < 0) else
                               "isaret donmedi; katsayi degisimi rapor edildi"))),
    )
    (T / "T3d_mechanism.json").write_text(json.dumps(res, indent=2))
    print(f"[17] mutlak: eGFR<15 {res['absolute_units_egfr_lt15_med']}U vs >=90 "
          f"{res['absolute_units_egfr_ge90_med']}U; fark {d:.2f} [{cl:.2f},{cu:.2f}], KW p={res['kw_units_by_band_p']}")
    print(f"[17] aracilik units~eGFR<15: crude beta={c1['beta']} {c1['ci']} -> adjusted(+weight+dm) beta={c2['beta']} {c2['ci']}")
    print(f"[17] {res['mediation']['interpretation']}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
24 — eICU KAPANIS: eksik kapatma (yeni analiz yok, mevcut veriden).
3 Mundlak kontrast (between-within) coef+CI+p | 4 tier1 vs tier2+ AYRIK
5 secilim parse-fail vs success | 6 ICC @51 merkez | 2 UC-esik TOST tablosu.
Girdi: data/eicu_glargine_events.csv + data/eicu_supp.csv. Cikti: qc/eicu_closing.json.
"""
import json
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

TH = {"KDIGO25_0.0056": 0.0056, "KDIGO30_0.0068": 0.0068, "Aronoff50_0.01125": 0.01125}
OUT = {}


def load_main():
    d = pd.read_csv("data/eicu_glargine_events.csv")
    for c in ["dose_u_kg", "egfr", "prev_dose_u_kg", "hrs_since_prev", "age", "weight"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["dm"] = d.dm.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
    return d


def load_supp():
    d = pd.read_csv("data/eicu_supp.csv")
    for c in ["dose_u_kg", "egfr", "prev_dose_u_kg", "hrs_since_prev", "age", "weight", "dcr_frac"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d.age >= 18) & d.weight.between(30, 300) & d.dose_u_kg.between(0.01, 3.0)]
    d["is_first"] = d.is_first_dose.astype(str).str.lower().isin(["true", "1"])
    return d


def patient_level(d):
    d = d[d.dose_u_kg.notna() & d.egfr.notna() & (d.age >= 18) & d.weight.between(30, 300)]
    p = d.groupby(["hospitalid", "patientunitstayid"]).agg(
        u_kg=("dose_u_kg", "median"), egfr=("egfr", "median"),
        age=("age", "median"), weight=("weight", "median"), dm=("dm", "max")).reset_index()
    return p


# 3) MUNDLAK KONTRAST: dose ~ egfr + egfr_cmean -> egfr_cmean = between-within contextual
def item3(main):
    p = patient_level(main)
    p["egfr_cmean"] = p.groupby("hospitalid").egfr.transform("mean")
    m = smf.ols("u_kg ~ egfr + egfr_cmean + age + weight + dm", data=p).fit(
        cov_type="cluster", cov_kwds={"groups": p.hospitalid})
    ci = m.conf_int().loc["egfr_cmean"]
    OUT["item3_mundlak_contrast"] = dict(
        contextual_between_minus_within_per10=round(float(m.params["egfr_cmean"]) * 10, 5),
        ci95_per10=[round(float(ci[0]) * 10, 5), round(float(ci[1]) * 10, 5)],
        p=float(m.pvalues["egfr_cmean"]),
        within_egfr_per10=round(float(m.params["egfr"]) * 10, 5),
        note="egfr_cmean anlamli (+, p<0.05) ise vaka-karmasi DOGRULANIR.")


# 4) A7 AYRIK: tier1 vs tier2+ (disjoint)
def item4(supp):
    feats = ["egfr", "age", "weight", "dm", "hrs_since_prev"]

    def fit(dd, label):
        t = dd[(~dd.is_first) & dd.prev_dose_u_kg.notna() & dd[feats].notna().all(axis=1)].copy()
        if len(t) < 150:
            return dict(set=label, n=int(len(t)), note="n<150 -- guc dusuk")
        r2wo = float(sm.OLS(t.dose_u_kg, sm.add_constant(t[feats])).fit().rsquared)
        r2w = float(sm.OLS(t.dose_u_kg, sm.add_constant(t[feats + ["prev_dose_u_kg"]])).fit().rsquared)
        m = smf.ols("dose_u_kg ~ egfr + age + weight + dm + hrs_since_prev", data=t).fit(
            cov_type="cluster", cov_kwds={"groups": t.hospitalid})
        ci = m.conf_int().loc["egfr"]
        return dict(set=label, n=int(len(t)), r2_without=round(r2wo, 4), r2_with=round(r2w, 4),
                    egfr_per10=round(float(m.params["egfr"]) * 10, 5),
                    egfr_ci_per10=[round(float(ci[0]) * 10, 5), round(float(ci[1]) * 10, 5)])
    t1 = supp[supp.tier1_exact == 1]; t2 = supp[supp.tier1_exact == 0]
    OUT["item4_tier_disjoint"] = dict(tier1_exact=fit(t1, "tier1_exact"),
                                      tier2plus=fit(t2, "tier2plus_DISJOINT"))


# 5) SECILIM: parse-fail vs parse-success
def item5(main):
    adult = main[(main.age >= 18) & main.weight.between(30, 300)]
    fail = adult[adult.dose_u_kg.isna()]; ok = adult[adult.dose_u_kg.notna()]

    def desc(x):
        return dict(n=int(len(x)), n_hosp=int(x.hospitalid.nunique()),
                    med_egfr=round(float(x.egfr.median()), 1),
                    med_weight=round(float(x.weight.median()), 1),
                    pct_dm=round(100 * float(x.dm.mean()), 1))
    # merkez-basi parse-basari orani -> merkez dozuyla iliski (kumelenme yanliligi)
    cf = adult.assign(ok=adult.dose_u_kg.notna().astype(int)).groupby("hospitalid").agg(
        parse_rate=("ok", "mean"), n=("ok", "size")).reset_index()
    cf = cf[cf.n >= 20]
    cmed = ok.groupby("hospitalid").dose_u_kg.median().rename("cdose")
    j = cf.set_index("hospitalid").join(cmed).dropna()
    corr = float(np.corrcoef(j.parse_rate, j.cdose)[0, 1]) if len(j) > 5 else None
    OUT["item5_selection"] = dict(parse_fail=desc(fail), parse_success=desc(ok),
        corr_center_parserate_vs_centerdose=round(corr, 3) if corr is not None else None,
        note="parse-fail ile success benzer eGFR/agirlik/DM ise klinik secilim yok; "
             "corr~0 ise merkez-parse-orani merkez-dozuyla iliskisiz -> ICC/between yanli DEGIL.")


# 6) ICC @ >=20 analiz-edilebilir hasta merkezleri
def item6(main):
    p = patient_level(main)
    nper = p.groupby("hospitalid").patientunitstayid.nunique()
    keep = nper[nper >= 20].index
    p20 = p[p.hospitalid.isin(keep)]
    md = smf.mixedlm("u_kg ~ 1", data=p20, groups=p20.hospitalid).fit(reml=True)
    vh = float(md.cov_re.iloc[0, 0]); ve = float(md.scale)
    OUT["item6_icc"] = dict(n_centers=int(p20.hospitalid.nunique()),
                            icc=round(vh / (vh + ve), 4),
                            note="ICC yalniz >=20 analiz-edilebilir hastali merkezlerde (A8 ile tutarli).")


# 2) UC-ESIK TOST tablosu (mevcut egim tahminleri)
def item2_tost():
    # (analiz, egim_per10, ust_CI_per10)  -- kapanis raporundaki degerler
    rows = {
        "A1_first_dose": (0.00125, 0.00225),
        "A2_titration": (0.00304, 0.00444),
        "A4_within_center": (0.00200, 0.00263),
        "A6_steady_within": (0.00190, 0.00265),
        "A6_aki_within": (0.00350, 0.00463),
        "A7_tier1_egfr": (0.00317, 0.00487),
        "A9_within_pat_simple": (0.00254, 0.00510),
        "A9_within_pat_adj": (0.00272, 0.00678),
    }
    tab = {}
    for k, (slope, up) in rows.items():
        tab[k] = dict(slope_per10=slope, upper_ci_per10=up,
                      equiv={name: bool(up < thr) for name, thr in TH.items()},
                      margin_to_KDIGO25_pct=round(100 * (0.0056 - up) / 0.0056, 1))
    OUT["item2_tost_3thresholds"] = dict(thresholds_per10=TH, table=tab)


def main():
    m = load_main(); s = load_supp()
    item3(m); item4(s); item5(m); item6(m); item2_tost()
    json.dump(OUT, open("qc/eicu_closing.json", "w"), indent=2)
    print(json.dumps(OUT, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
23 — eICU A6 (steady-state/AKI) + A7 (parse gurultusu) + A8 (merkez sayisi tutarsizlik).
Girdi: data/eicu_supp.csv (gitignore). Cikti: qc/eicu_a6a7.json.
"""
import json
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

ARONOFF10 = 0.09 / 80.0 * 10   # 0.01125 U/kg per 10 mL/dk
OUT = {}


def load():
    d = pd.read_csv("data/eicu_supp.csv")
    for c in ["dose_u_kg", "egfr", "prev_dose_u_kg", "hrs_since_prev", "age", "weight", "dcr_frac"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d.age >= 18) & d.weight.between(30, 300) & d.dose_u_kg.between(0.01, 3.0)]
    d["is_first"] = d.is_first_dose.astype(str).str.lower().isin(["true", "1"])
    return d


def patient_level(d):
    p = d.groupby(["hospitalid", "patientunitstayid"]).agg(
        u_kg=("dose_u_kg", "median"), egfr=("egfr", "median"),
        age=("age", "median"), weight=("weight", "median"), dm=("dm", "max")).reset_index()
    p["band"] = pd.cut(p.egfr, [-1, 15, 45, 90, 1e9], labels=["lt15", "b15_45", "b45_90", "ge90"]).astype(str)
    return p


def band_means(p):
    r = {}
    for b in ["lt15", "b15_45", "b45_90", "ge90"]:
        s = p[p.band == b]
        if len(s) < 20:
            r[b] = dict(n=int(len(s)), note="n<20"); continue
        m = smf.ols("u_kg ~ 1", data=s).fit(cov_type="cluster", cov_kwds={"groups": s.hospitalid})
        ci = m.conf_int().loc["Intercept"]
        r[b] = dict(n=int(len(s)), median=round(float(s.u_kg.median()), 4),
                    ci=[round(float(ci[0]), 4), round(float(ci[1]), 4)])
    return r


def mundlak_within(p):
    cmean = p.groupby("hospitalid").egfr.transform("mean")
    p = p.assign(egfr_cmean=cmean, egfr_within=p.egfr - cmean)
    m = smf.ols("u_kg ~ egfr_within + egfr_cmean + age + weight + dm", data=p).fit(
        cov_type="cluster", cov_kwds={"groups": p.hospitalid})
    bw = float(m.params["egfr_within"]); se = float(m.bse["egfr_within"])
    ci90 = [(bw - 1.645 * se) * 10, (bw + 1.645 * se) * 10]
    return dict(within_per10=round(bw * 10, 5), within_90ci=[round(ci90[0], 5), round(ci90[1], 5)],
                between_per10=round(float(m.params["egfr_cmean"]) * 10, 5),
                tost_negligible=bool(ci90[1] < ARONOFF10 and ci90[0] > -ARONOFF10))


def a6(d):
    steady = d[d.dcr_frac.notna() & (d.dcr_frac < 0.20)]
    aki = d[d.dcr_frac.notna() & (d.dcr_frac >= 0.20)]
    ps, pa, pf = patient_level(steady), patient_level(aki), patient_level(d)
    OUT["A6_steady_state"] = dict(
        aronoff_per10=round(ARONOFF10, 5),
        n_steady_orders=int(len(steady)), n_steady_patients=int(steady.patientunitstayid.nunique()),
        n_aki_orders=int(len(aki)),
        steady_band_means=band_means(ps), steady_mundlak=mundlak_within(ps),
        aki_mundlak=(mundlak_within(pa) if pa.egfr.nunique() > 10 else "n_kucuk"),
        full_mundlak=mundlak_within(pf),
        note="AKI'da eGFR anlamsizlasir. Steady (|dCr|/Cr<%20) alt kumede within ~0 kaliyorsa bulgu SAGLAM.")


def a7(d):
    tier1 = d[d.tier1_exact == 1]
    feats = ["egfr", "age", "weight", "dm", "hrs_since_prev"]

    def a2fit(dd):
        t = dd[(~dd.is_first) & dd.prev_dose_u_kg.notna() & dd[feats].notna().all(axis=1)].copy()
        if len(t) < 200:
            return dict(n=int(len(t)), note="n<200")
        r2_wo = float(sm.OLS(t.dose_u_kg, sm.add_constant(t[feats])).fit().rsquared)
        r2_w = float(sm.OLS(t.dose_u_kg, sm.add_constant(t[feats + ["prev_dose_u_kg"]])).fit().rsquared)
        m = smf.ols("dose_u_kg ~ egfr + age + weight + dm + hrs_since_prev", data=t).fit(
            cov_type="cluster", cov_kwds={"groups": t.hospitalid})
        ci = m.conf_int().loc["egfr"]
        return dict(n=int(len(t)), r2_insample_without=round(r2_wo, 4), r2_insample_with=round(r2_w, 4),
                    egfr_per10=round(float(m.params["egfr"]) * 10, 5),
                    egfr_ci_per10=[round(float(ci[0]) * 10, 5), round(float(ci[1]) * 10, 5)])
    OUT["A7_parse_noise"] = dict(
        pct_tier1_exact=round(100 * float((d.tier1_exact == 1).mean()), 1),
        all_clean=a2fit(d), tier1_exact_only=a2fit(tier1),
        note="tier1=kesin 'N units'. tier1 ile tum-temiz benzerse -> egfr-yok bulgusu parse gurultusu DEGIL.")


def a8(d):
    pf = patient_level(d)
    nper = pf.groupby("hospitalid").patientunitstayid.nunique()
    OUT["A8_center_counts"] = dict(
        total_glargine_centers_IS2=91,
        centers_with_analyzable_patients=int(pf.hospitalid.nunique()),
        centers_ge20_analyzable=int((nper >= 20).sum()),
        aciklama=("Tutarsizlik cozuldu: 91 toplam glarjin merkezi (IS2); "
                  f"{int(pf.hospitalid.nunique())} merkezde analiz-edilebilir (doz+eGFR) hasta; "
                  f"{int((nper >= 20).sum())} merkezde >=20 analiz-edilebilir hasta. "
                  "v1'deki '51' (merkez-ici polyfit) ve '54' (>=20 dozlu) YERINE bu tek "
                  "cerceve; merkez-ici titrasyon artik A4-Mundlak (havuz) + A9 (within-patient) ile."))


def main():
    d = load()
    a6(d); a7(d); a8(d)
    json.dump(OUT, open("qc/eicu_a6a7.json", "w"), indent=2)
    print(json.dumps(OUT, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

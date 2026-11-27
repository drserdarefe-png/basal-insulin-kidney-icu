#!/usr/bin/env python3
"""
22 — eICU dozlama kolu REVIZYON (v2). Mevcut veri; yeni tablo YOK.
A1 ilk-doz (naive+cluster CI) | A2 4-sayi R2 in-sample+CV, ayni 5 ozellik, DUR(in-sample)
huni (attrition + dusen/kalan) | A3 kumelenmis band ort. | A4 TOST+Mundlak | A5 ICC+shrunken.
Birim: TUM egfr katsayilari U/kg per 10 mL/dk.
Girdi: data/eicu_glargine_events.csv, MIMIC data/basal_events/*.parquet.
Cikti: qc/eicu_replication_v2.{json,md}. Hasta-duzeyi veri COMMIT EDILMEZ.
"""
import json
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from sklearn.model_selection import GroupKFold
import lightgbm as lgb
import duckdb

OUT = {}
FEATS5 = ["egfr", "age", "weight", "dm", "hrs_since_prev"]   # iki kohortta ORTAK
ARONOFF_SLOPE_PER_ML = 0.09 / 80.0    # eGFR 90->10 icin %50 azalma ~ 0.00112 U/kg per mL/dk


def _b(s):
    return s.astype(str).str.lower().isin(["true", "1", "1.0"])


def load_eicu_raw():
    df = pd.read_csv("data/eicu_glargine_events.csv")
    for c in ["dose_u_kg", "prev_dose_u_kg", "egfr", "age", "weight", "scr", "hrs_since_prev", "order_seq"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["dm"] = _b(df.dm).astype(float)
    df["is_first"] = _b(df.is_first_dose)
    return df


def load_mimic():
    q = ("SELECT subject_id, egfr_now AS egfr, weight_kg AS weight, CAST(dm_any AS DOUBLE) AS dm, "
         "age, prev_glargine_dose_u_kg AS prev_dose_u_kg, hrs_since_prev_glargine AS hrs_since_prev, "
         "dose_u_kg FROM read_parquet('data/basal_events/*.parquet') "
         "WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL")
    df = duckdb.connect().execute(q).df()
    df = df[df.weight.between(30, 300) & df.dose_u_kg.between(0.01, 3.0)]
    df["is_first"] = df.prev_dose_u_kg.isna()
    return df


def egfr_ci(d, formula, group, per=10.0):
    """OLS: naive CI + cluster-robust CI, egfr katsayisi (per 10 mL/dk)."""
    m0 = smf.ols(formula, data=d).fit()
    mc = smf.ols(formula, data=d).fit(cov_type="cluster", cov_kwds={"groups": d[group]})
    c = float(m0.params["egfr"]) * per
    ci0 = [float(x) * per for x in m0.conf_int().loc["egfr"]]
    cic = [float(x) * per for x in mc.conf_int().loc["egfr"]]
    return dict(n=int(m0.nobs), egfr_per10=round(c, 5), ci_naive=[round(ci0[0], 5), round(ci0[1], 5)],
                ci_cluster=[round(cic[0], 5), round(cic[1], 5)], r2_insample=round(float(m0.rsquared), 4))


def ols_r2_insample(d, feats):
    X = sm.add_constant(d[feats]); return float(sm.OLS(d.dose_u_kg, X).fit().rsquared)


def cv_r2(d, feats, group):
    gkf = GroupKFold(5); sr = st = 0.0
    for tr, te in gkf.split(d, d.dose_u_kg, d[group]):
        m = lgb.LGBMRegressor(n_estimators=300, learning_rate=0.05, num_leaves=31,
                              min_child_samples=40, subsample=0.8, colsample_bytree=0.8, verbose=-1)
        m.fit(d[feats].iloc[tr], d.dose_u_kg.iloc[tr]); p = m.predict(d[feats].iloc[te])
        sr += ((d.dose_u_kg.iloc[te] - p) ** 2).sum(); st += ((d.dose_u_kg.iloc[te] - d.dose_u_kg.iloc[te].mean()) ** 2).sum()
    return 1.0 - sr / st


# ---------- HUNI ----------
def funnel(raw):
    steps = {}
    steps["1_order_events(tablo)"] = int(len(raw))
    a = raw[raw.age >= 18]; steps["2_adult>=18"] = int(len(a))
    w = a[a.weight.between(30, 300)]; steps["3_weight_30_300"] = int(len(w))
    dpar = w[w.dose_u_kg.notna()]; steps["4_dose_parseable"] = int(len(dpar))
    eg = dpar[dpar.egfr.notna()]; steps["5_creatinine<=48h(egfr)"] = int(len(eg))
    pv = eg[eg.prev_dose_u_kg.notna()]; steps["6_prev_dose(titrasyon)"] = int(len(pv))
    # dusen vs kalan: kreatinin adiminda (5) dropped=dpar without egfr
    dropped = dpar[dpar.egfr.isna()]; kept = eg
    cmp = dict(step="creatinine<=48h",
               dropped=dict(n=int(len(dropped)), n_hosp=int(dropped.hospitalid.nunique()),
                            med_dose_u_kg=round(float(dropped.dose_u_kg.median()), 3)),
               kept=dict(n=int(len(kept)), n_hosp=int(kept.hospitalid.nunique()),
                         med_dose_u_kg=round(float(kept.dose_u_kg.median()), 3),
                         med_egfr=round(float(kept.egfr.median()), 1)))
    OUT["A2c_funnel"] = dict(steps=steps, drop_vs_keep=cmp)


# ---------- A1 ----------
def a1(eicu, mimic):
    e = eicu[eicu.is_first & eicu.dose_u_kg.notna()].dropna(subset=FEATS5[:4]).copy()
    ee = egfr_ci(e, "dose_u_kg ~ egfr + age + weight + dm", "hospitalid")
    m = mimic[mimic.is_first].dropna(subset=["dose_u_kg"] + FEATS5[:4]).copy()
    m["subject_id"] = m.subject_id.astype(int)
    mm = egfr_ci(m, "dose_u_kg ~ egfr + age + weight + dm", "subject_id")
    OUT["A1_first_dose_per10mL"] = {"eICU": ee, "MIMIC": mm, "aronoff_slope_per10mL": round(ARONOFF_SLOPE_PER_ML * 10, 5)}


# ---------- A2 ----------
def a2(eicu, mimic):
    e = eicu[(~eicu.is_first) & eicu.dose_u_kg.notna() & eicu.prev_dose_u_kg.notna()
             & eicu[FEATS5].notna().all(axis=1)].copy()
    e["grp"] = e.patientunitstayid
    m = mimic[(~mimic.is_first) & mimic.prev_dose_u_kg.notna() & mimic[FEATS5].notna().all(axis=1)].copy()
    m["grp"] = m.subject_id.astype(int)

    def four(d):
        return dict(
            insample_without=round(ols_r2_insample(d, FEATS5), 4),
            insample_with=round(ols_r2_insample(d, FEATS5 + ["prev_dose_u_kg"]), 4),
            cv_without=round(cv_r2(d, FEATS5, "grp"), 4),
            cv_with=round(cv_r2(d, FEATS5 + ["prev_dose_u_kg"], "grp"), 4), n=int(len(d)))
    fe, fm = four(e), four(m)
    ce = egfr_ci(e, "dose_u_kg ~ egfr + age + weight + dm + hrs_since_prev", "hospitalid")
    dur = fe["insample_without"] > 0.3
    OUT["A2_same_rows"] = dict(features=FEATS5, eICU=fe, MIMIC=fm, eICU_egfr_per10=ce,
                               DUR_insample_without_gt0_3=bool(dur))
    return dur


# ---------- patient-level (T3/A4/A5) ----------
def patient_level(eicu):
    d = eicu[eicu.dose_u_kg.notna() & eicu.egfr.notna() & (eicu.age >= 18) & eicu.weight.between(30, 300)]
    pat = d.groupby(["hospitalid", "patientunitstayid"]).agg(
        u_kg=("dose_u_kg", "median"), egfr=("egfr", "median"),
        age=("age", "median"), weight=("weight", "median"), dm=("dm", "max")).reset_index()
    pat["band"] = pd.cut(pat.egfr, [-1, 15, 45, 90, 1e9], labels=["lt15", "b15_45", "b45_90", "ge90"]).astype(str)
    return pat


# ---------- A3 kumelenmis band ort. ----------
def a3(pat):
    d = pat.copy()
    m = smf.ols("u_kg ~ C(band, Treatment(reference='ge90'))", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": d.hospitalid})
    # band ortalamalari: naive (grup medyan/ort) + cluster CI on mean
    res = {}
    for b in ["lt15", "b15_45", "b45_90", "ge90"]:
        s = d[d.band == b]
        mm = smf.ols("u_kg ~ 1", data=s).fit(cov_type="cluster", cov_kwds={"groups": s.hospitalid})
        ci = mm.conf_int().loc["Intercept"]
        res[b] = dict(n=int(len(s)), mean_u_kg=round(float(s.u_kg.mean()), 4),
                      median_u_kg=round(float(s.u_kg.median()), 4),
                      ci_cluster=[round(float(ci[0]), 4), round(float(ci[1]), 4)])
    joint_p = float(m.f_pvalue)
    OUT["A3_clustered_bands"] = dict(by_band=res,
        joint_wald_p=("p<1e-16" if joint_p < 1e-16 else f"{joint_p:.2e}"),
        note="hasta-duzeyi u_kg, hospitalid-cluster CI. KW KALDIRILDI.")


# ---------- A4 TOST + Mundlak ----------
def a4(pat):
    d = pat.copy()
    cmean = d.groupby("hospitalid").egfr.transform("mean")
    d["egfr_cmean"] = cmean; d["egfr_within"] = d.egfr - cmean
    m = smf.ols("u_kg ~ egfr_within + egfr_cmean + age + weight + dm", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": d.hospitalid})
    bw = float(m.params["egfr_within"]); se = float(m.bse["egfr_within"])
    ci90 = [bw - 1.645 * se, bw + 1.645 * se]
    bound = ARONOFF_SLOPE_PER_ML
    tost_equiv = (ci90[1] < bound) and (ci90[0] > -bound)   # 90% CI guideline-esiginin ICINDE
    bc = float(m.params["egfr_cmean"]); bcci = m.conf_int().loc["egfr_cmean"]
    OUT["A4_tost_mundlak"] = dict(
        within_slope_per10mL=round(bw * 10, 5), within_90ci_per10mL=[round(ci90[0] * 10, 5), round(ci90[1] * 10, 5)],
        between_slope_per10mL=round(bc * 10, 5), between_ci95_per10mL=[round(float(bcci[0]) * 10, 5), round(float(bcci[1]) * 10, 5)],
        aronoff_bound_per10mL=round(bound * 10, 5), tost_equivalent_to_negligible=bool(tost_equiv),
        note="within=merkez-ici titrasyon; between=vaka-karmasi. TOST: within 90%CI Aronoff-esiginin icindeyse 'ihmal edilebilir'.")


# ---------- A5 ICC + shrunken ----------
def a5(pat):
    d = pat.copy()
    md = smf.mixedlm("u_kg ~ 1", data=d, groups=d.hospitalid).fit(reml=True)
    v_h = float(md.cov_re.iloc[0, 0]); v_e = float(md.scale)
    icc = v_h / (v_h + v_e)
    fe = float(md.fe_params["Intercept"])
    re = md.random_effects
    shrunk = np.array([fe + float(v.iloc[0]) for v in re.values()])
    # ham merkez ortalamalari (>=20 hasta)
    raw_c = d.groupby("hospitalid").u_kg.median()
    n_c = d.groupby("hospitalid").patientunitstayid.nunique()
    raw20 = raw_c[n_c >= 20]
    OUT["A5_icc_shrunken"] = dict(
        icc=round(icc, 3), n_centers=int(len(re)),
        shrunken_range=[round(float(shrunk.min()), 3), round(float(shrunk.max()), 3)],
        shrunken_p10_p90=[round(float(np.percentile(shrunk, 10)), 3), round(float(np.percentile(shrunk, 90)), 3)],
        raw_median_range_ge20=[round(float(raw20.min()), 3), round(float(raw20.max()), 3)],
        note="shrunken=EB/BLUP; ham aralik parantezde.")


# ---------- A9 WITHIN-PATIENT (sabit-etki) titrasyon ----------
def a9(eicu):
    d = eicu[eicu.dose_u_kg.notna() & eicu.egfr.notna() & (eicu.age >= 18) & eicu.weight.between(30, 300)].copy()
    # stay-ici eGFR araligi ve >=2 order
    g = d.groupby("patientunitstayid")
    rng = g.egfr.transform(lambda s: s.max() - s.min())
    cnt = g.dose_u_kg.transform("count")
    sub = d[(rng >= 15) & (cnt >= 2)].copy()
    # basit within: dose ~ egfr  (hasta-ici demean)
    for col in ["dose_u_kg", "egfr", "prev_dose_u_kg", "hrs_since_prev"]:
        sub[col + "_dm"] = sub[col] - sub.groupby("patientunitstayid")[col].transform("mean")
    m_simple = smf.ols("dose_u_kg_dm ~ egfr_dm - 1", data=sub).fit(
        cov_type="cluster", cov_kwds={"groups": sub.patientunitstayid})
    ci_s = m_simple.conf_int().loc["egfr_dm"]
    # ayarli within (prev + hrs mevcut olan order'lar): non-first
    adj = sub[sub.prev_dose_u_kg.notna() & sub.hrs_since_prev.notna()].copy()
    for col in ["dose_u_kg", "egfr", "prev_dose_u_kg", "hrs_since_prev"]:
        adj[col + "_dm"] = adj[col] - adj.groupby("patientunitstayid")[col].transform("mean")
    m_adj = smf.ols("dose_u_kg_dm ~ egfr_dm + prev_dose_u_kg_dm + hrs_since_prev_dm - 1", data=adj).fit(
        cov_type="cluster", cov_kwds={"groups": adj.patientunitstayid})
    ci_a = m_adj.conf_int().loc["egfr_dm"]
    OUT["A9_within_patient"] = dict(
        subset="stay-ici eGFR araligi>=15 mL/dk & >=2 order",
        n_patients=int(sub.patientunitstayid.nunique()), n_orders=int(len(sub)),
        simple_within_egfr_per10mL=round(float(m_simple.params["egfr_dm"]) * 10, 5),
        simple_ci95_per10mL=[round(float(ci_s[0]) * 10, 5), round(float(ci_s[1]) * 10, 5)],
        adj_n_orders=int(len(adj)),
        adj_within_egfr_per10mL=round(float(m_adj.params["egfr_dm"]) * 10, 5),
        adj_ci95_per10mL=[round(float(ci_a[0]) * 10, 5), round(float(ci_a[1]) * 10, 5)],
        aronoff_per10mL=round(ARONOFF_SLOPE_PER_ML * 10, 5),
        note="within-patient egim ~0 -> ayni hastada bobrek degisince bile doz degismiyor (EN GUCLU).")


def main():
    eicu = load_eicu_raw(); mimic = load_mimic()
    OUT["_n"] = dict(eicu_table_rows=int(len(eicu)), mimic_glargine_events=int(len(mimic)))
    funnel(eicu)
    a1(eicu, mimic)
    dur = a2(eicu, mimic)
    pat = patient_level(eicu)
    a3(pat); a4(pat); a5(pat); a9(eicu)
    json.dump(OUT, open("qc/eicu_replication_v2.json", "w"), indent=2)
    print(json.dumps(OUT, indent=2, ensure_ascii=False))
    print("\nDUR (A2 in-sample without>0.3):", "!!! TETIKLENDI" if dur else "yok")


if __name__ == "__main__":
    main()

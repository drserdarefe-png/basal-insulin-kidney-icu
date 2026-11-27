#!/usr/bin/env python3
"""
26 — F9 KAPANIS (yeni analiz yok, hukum netlestirme).
1 Mundlak kontrast (within-between) coef+CI+p
2 BASAMAK esik (lineerlestirme YOK): I(egfr<60) ve 3-kategori Aronoff
3 BANT-OZGU oran (politika x eGFR-bin) + kumele-bootstrap CI
4 BANT-N mutabakati (event vs patient, huniye bagli)
Girdi: data/eicu_glargine_events.csv, MIMIC parquet. Cikti: qc/fig4_f9.json.
p-degerleri yuvarlanir (p<1e-16). Hasta-duzeyi veri COMMIT EDILMEZ.
"""
import json
import numpy as np, pandas as pd, duckdb
import statsmodels.formula.api as smf

RNG = np.random.default_rng(20260717)
OUT = {}


def pround(p):
    return "p<1e-16" if p < 1e-16 else round(float(p), 4)


def mult_p1(e):
    return 1.0 if e > 50 else (0.75 if e >= 10 else 0.50)
def mult_p2(e):
    return 1.0 if e >= 60 else 0.70
def mult_p3(e):
    return 1.0 if e >= 60 else 0.75
POL = {"P1_Aronoff": mult_p1, "P2_KDIGO_glargine": mult_p2, "P3_KDIGO_TDD": mult_p3}


def load_eicu():
    e = pd.read_csv("data/eicu_glargine_events.csv")
    for c in ["dose_u_kg", "egfr", "weight", "age"]:
        e[c] = pd.to_numeric(e[c], errors="coerce")
    e = e[(e.age >= 18) & e.weight.between(30, 300) & e.dose_u_kg.between(0.01, 3.0) & e.egfr.notna()]
    e["dm"] = e.dm.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
    return e.rename(columns={"patientunitstayid": "pid", "hospitalid": "cid"})


def load_mimic():
    q = ("SELECT subject_id AS pid, subject_id AS cid, egfr_now AS egfr, dose_u_kg, "
         "weight_kg AS weight, age, CAST(dm_any AS DOUBLE) AS dm "
         "FROM read_parquet('data/basal_events/*.parquet') "
         "WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL AND egfr_now IS NOT NULL")
    d = duckdb.connect().execute(q).df()
    return d[d.weight.between(30, 300) & d.dose_u_kg.between(0.01, 3.0)]


def patient_level(order):
    agg = {"u_kg": ("dose_u_kg", "median"), "egfr": ("egfr", "median"),
           "age": ("age", "median"), "weight": ("weight", "median"), "dm": ("dm", "max")}
    for pol in POL:
        order[pol] = order.egfr.map(POL[pol])
        agg[pol] = (pol, "mean")
    return order.groupby(["cid", "pid"]).agg(**agg).reset_index()


# 1) MUNDLAK KONTRAST
def item1(pat, label):
    p = pat.copy()
    p["cmean"] = p.groupby("cid").egfr.transform("mean")
    p["ew"] = p.egfr - p.cmean
    m = smf.ols("u_kg ~ ew + cmean + age + weight + dm", data=p).fit(
        cov_type="cluster", cov_kwds={"groups": p.cid})
    # within - between = ew - cmean : lineer kombinasyon
    diff = float(m.params["ew"] - m.params["cmean"])
    # SE(diff) cluster kovaryanstan
    cov = m.cov_params()
    var = cov.loc["ew", "ew"] + cov.loc["cmean", "cmean"] - 2 * cov.loc["ew", "cmean"]
    se = float(np.sqrt(var))
    ci = [(diff - 1.96 * se) * 10, (diff + 1.96 * se) * 10]
    from scipy import stats
    z = diff / se; pval = 2 * (1 - stats.norm.cdf(abs(z)))
    OUT.setdefault("item1_mundlak_contrast", {})[label] = dict(
        within_per10=round(float(m.params["ew"]) * 10, 5),
        between_per10=round(float(m.params["cmean"]) * 10, 5),
        within_minus_between_per10=round(diff * 10, 5),
        ci95_per10=[round(ci[0], 5), round(ci[1], 5)], p=pround(pval),
        note="p>0.05 -> within!=between DOGRULANAMAZ -> 'vaka-karmasi' savunmasi ZAYIF.")


# 2) BASAMAK ESIK (lineerlestirme yok)
def item2(pat, label):
    p = pat.copy()
    anchor60 = float(p.loc[p.egfr >= 60, "u_kg"].median())
    p["lt60"] = (p.egfr < 60).astype(float)
    m = smf.ols("u_kg ~ lt60 + age + weight + dm", data=p).fit(
        cov_type="cluster", cov_kwds={"groups": p.cid})
    ci = m.conf_int().loc["lt60"]
    step = dict(anchor_ge60=round(anchor60, 4),
                observed_step_lt60=round(float(m.params["lt60"]), 4),
                ci95=[round(float(ci[0]), 4), round(float(ci[1]), 4)], p=pround(m.pvalues["lt60"]),
                guideline_P2_step_m30=round(-0.30 * anchor60, 4),
                guideline_P3_step_m25=round(-0.25 * anchor60, 4),
                excludes_P3=bool(float(ci[1]) < -0.25 * anchor60 or float(ci[0]) > -0.25 * anchor60))
    # 3-kategori Aronoff
    anchor50 = float(p.loc[p.egfr > 50, "u_kg"].median())
    p["cat"] = np.where(p.egfr < 10, "c_lt10", np.where(p.egfr <= 50, "b_10_50", "a_gt50"))
    m3 = smf.ols("u_kg ~ C(cat, Treatment(reference='a_gt50')) + age + weight + dm", data=p).fit(
        cov_type="cluster", cov_kwds={"groups": p.cid})
    cats = {}
    for k, gl in [("b_10_50", -0.25 * anchor50), ("c_lt10", -0.50 * anchor50)]:
        nm = [c for c in m3.params.index if k in c]
        if nm:
            co = float(m3.params[nm[0]]); ci3 = m3.conf_int().loc[nm[0]]
            cats[k] = dict(n=int((p.cat == k).sum()), observed=round(co, 4),
                           ci=[round(float(ci3[0]), 4), round(float(ci3[1]), 4)],
                           guideline_step=round(gl, 4), p=pround(m3.pvalues[nm[0]]))
        else:
            cats[k] = dict(n=int((p.cat == k).sum()), note="kategori bos/tekil")
    OUT.setdefault("item2_step_threshold", {})[label] = dict(
        anchor50=round(anchor50, 4), step_lt60=step, aronoff_3cat=cats,
        note="LINEERLESTIRME YOK; basamak katsayisi vs kilavuz-basamagi.")


# 3) BANT-OZGU ORAN + kumele bootstrap CI
def cluster_boot_median(vals, groups, n=1500):
    ug = np.unique(groups); idx = {g: np.where(groups == g)[0] for g in ug}
    out = []
    for _ in range(n):
        samp = RNG.choice(ug, size=len(ug), replace=True)
        take = np.concatenate([idx[g] for g in samp])
        out.append(np.median(vals[take]))
    return [round(float(np.percentile(out, 2.5)), 3), round(float(np.percentile(out, 97.5)), 3)]


def item3(pat, label):
    p = pat.copy()
    bins = [(-1, 10, "lt10"), (10, 30, "b10_30"), (30, 60, "b30_60"), (60, 90, "b60_90"), (90, 1e9, "ge90")]
    p["ebin"] = "na"
    for lo, hi, nm in bins:
        p.loc[(p.egfr > lo) & (p.egfr <= hi), "ebin"] = nm
    res = {}
    for pol in POL:
        anchor = float(p.loc[(p.egfr > 50) if pol == "P1_Aronoff" else (p.egfr >= 60), "u_kg"].median())
        p["_tgt"] = anchor * p[pol]
        p["_r"] = p.u_kg / p["_tgt"]
        per = {}
        for _, _, nm in bins:
            s = p[p.ebin == nm]
            if len(s) < 20:
                per[nm] = dict(n=int(len(s)), note="n<20"); continue
            ci = cluster_boot_median(s._r.to_numpy(), s.cid.to_numpy())
            per[nm] = dict(n=int(len(s)), median_ratio=round(float(s._r.median()), 3), ci=ci)
        res[pol] = dict(anchor=round(anchor, 4), by_egfr_bin=per)
    OUT.setdefault("item3_bandspecific_ratio", {})[label] = dict(
        by_policy=res, note="HAVUZ ORANI FIGURE GIRMEZ; capa bandinda oran~1 tanimla. eGFR<10 kritik.")


# 4) BANT-N mutabakati
def item4(eorder, morder):
    def counts(order, total_label):
        pl = patient_level(order.copy())
        ev = order.assign(band=pd.cut(order.egfr, [-1, 15, 45, 90, 1e9],
                          labels=["lt15", "b15_45", "b45_90", "ge90"]).astype(str))
        plb = pl.assign(band=pd.cut(pl.egfr, [-1, 15, 45, 90, 1e9],
                        labels=["lt15", "b15_45", "b45_90", "ge90"]).astype(str))
        return dict(total_events=int(len(order)), total_patients=int(len(pl)),
                    event_band={b: int((ev.band == b).sum()) for b in ["lt15", "b15_45", "b45_90", "ge90"]},
                    patient_band={b: int((plb.band == b).sum()) for b in ["lt15", "b15_45", "b45_90", "ge90"]})
    OUT["item4_band_n"] = dict(
        MIMIC=counts(morder, "MIMIC"), eICU=counts(eorder, "eICU"),
        note=("MIMIC <15: patient=85 vs event=774 -> Fig4 v9 HASTA-duzeyi (85); "
              "v8'in 774'u EVENT-duzeyi. Huni: MIMIC eGFR-seti 7.897 event / 1.108 hasta."))


def main():
    e = load_eicu(); m = load_mimic()
    pe = patient_level(e.copy()); pm = patient_level(m.copy())
    for pat, lab in [(pm, "MIMIC"), (pe, "eICU")]:
        item1(pat, lab); item2(pat, lab); item3(pat, lab)
    item4(e, m)
    json.dump(OUT, open("qc/fig4_f9.json", "w"), indent=2)
    print(json.dumps(OUT, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

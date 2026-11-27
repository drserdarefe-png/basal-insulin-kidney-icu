#!/usr/bin/env python3
"""
25 — FIGUR 4 v9 + s3.3/s3.4 yeniden hesap (iki kohort).
F0 KONVANSIYON: u_kg~egfr, coef>0 = yuksek eGFRde yuksek doz = KILAVUZ YONU.
F1 per-ORDER kilavuz carpani (3 politika) | F2 capa=azaltma-istenmeyen bant medyani
F3 oran iki yonlu | F4 kumelenmis Wald + slope vs 3 esik + kumele bootstrap
F5 iki-panel violin SVG | F6 N-huni | F7 s3.3 | F8 dil.
Girdi: data/eicu_glargine_events.csv, MIMIC parquet. Cikti: qc/fig4_v9.{json,md},
qc/paper2_figs/fig4_v9.svg. Hasta-duzeyi veri COMMIT EDILMEZ.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import statsmodels.formula.api as smf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RNG = np.random.default_rng(20260717)
BANDS = ["lt15", "b15_45", "b45_90", "ge90"]
BLBL = {"lt15": "<15", "b15_45": "15-45", "b45_90": "45-90", "ge90": ">=90"}
# 3 politika thresholdlari /10 mL (KDIGO-TDD %25 / KDIGO-glarjin %30 / Aronoff %50)
POL_SLOPE10 = {"P3_KDIGO_TDD25": 0.0056, "P2_KDIGO_glargine30": 0.0068, "P1_Aronoff50": 0.01125}
OUT = {}


def mult_p1(e):  # Aronoff
    return 1.0 if e > 50 else (0.75 if e >= 10 else 0.50)
def mult_p2(e):  # KDIGO glargine
    return 1.0 if e >= 60 else 0.70
def mult_p3(e):  # KDIGO TDD
    return 1.0 if e >= 60 else 0.75
POLICIES = {"P1_Aronoff": mult_p1, "P2_KDIGO_glargine": mult_p2, "P3_KDIGO_TDD": mult_p3}
ANCHOR_MASK = {"P1_Aronoff": lambda e: e > 50, "P2_KDIGO_glargine": lambda e: e >= 60,
               "P3_KDIGO_TDD": lambda e: e >= 60}


def band(egfr):
    return pd.cut(egfr, [-1, 15, 45, 90, 1e9], labels=BANDS).astype(str)


def load_eicu():
    e = pd.read_csv("data/eicu_glargine_events.csv")
    for c in ["dose_u_kg", "egfr", "weight", "age"]:
        e[c] = pd.to_numeric(e[c], errors="coerce")
    e = e[(e.age >= 18) & e.weight.between(30, 300) & e.dose_u_kg.between(0.01, 3.0) & e.egfr.notna()]
    return e.rename(columns={"patientunitstayid": "pid", "hospitalid": "cid"})[["pid", "cid", "dose_u_kg", "egfr"]]


def load_mimic():
    q = ("SELECT subject_id AS pid, subject_id AS cid, egfr_now AS egfr, dose_u_kg, weight_kg AS weight "
         "FROM read_parquet('data/basal_events/*.parquet') "
         "WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL AND egfr_now IS NOT NULL")
    d = duckdb.connect().execute(q).df()
    d = d[d.weight.between(30, 300) & d.dose_u_kg.between(0.01, 3.0)]
    return d[["pid", "cid", "dose_u_kg", "egfr"]]


def cluster_ci_mean(s, groups):
    m = smf.ols("y ~ 1", data=pd.DataFrame({"y": s.values, "g": groups.values})).fit(
        cov_type="cluster", cov_kwds={"groups": groups.values})
    ci = m.conf_int().loc["Intercept"]
    return float(s.mean()), [float(ci[0]), float(ci[1])]


def boot_banddiff(pat, n=2000):
    ug = pat.pid.unique()
    idx = {g: pat.index[pat.pid == g].to_numpy() for g in ug}
    a = pat.u_kg.to_numpy(); bnd = pat.band.to_numpy()
    diffs = []
    for _ in range(n):
        samp = RNG.choice(ug, size=len(ug), replace=True)
        take = np.concatenate([idx[g] for g in samp])
        av, bv = a[take], bnd[take]
        lo = np.median(av[bv == "lt15"]) if (bv == "lt15").any() else np.nan
        hi = np.median(av[bv == "ge90"]) if (bv == "ge90").any() else np.nan
        diffs.append(lo - hi)
    diffs = np.array(diffs)
    return round(float(np.nanmedian(diffs)), 4), [round(float(np.nanpercentile(diffs, 2.5)), 4),
                                                  round(float(np.nanpercentile(diffs, 97.5)), 4)]


def analyze(order, label):
    order = order.copy()
    order["band"] = band(order.egfr)
    # per-ORDER carpan
    for pol, f in POLICIES.items():
        order[pol] = order.egfr.map(f)
    # patient-duzeyi
    agg = {"u_kg": ("dose_u_kg", "median"), "egfr": ("egfr", "median")}
    for pol in POLICIES:
        agg[pol] = (pol, "mean")
    pat = order.groupby(["cid", "pid"]).agg(**agg).reset_index()
    pat["band"] = band(pat.egfr)
    res = {"n_orders": int(len(order)), "n_patients": int(len(pat)),
           "n_centers": int(pat.cid.nunique())}
    # band ort + cluster CI
    bres = {}
    for b in BANDS:
        s = pat[pat.band == b]
        if len(s) < 20:
            bres[b] = dict(n=int(len(s)), note="n<20"); continue
        mean, ci = cluster_ci_mean(s.u_kg, s.cid)
        bres[b] = dict(n=int(len(s)), mean=round(mean, 4), median=round(float(s.u_kg.median()), 4),
                       ci=[round(ci[0], 4), round(ci[1], 4)])
    res["bands"] = bres
    # slope + cluster CI /10mL + joint Wald
    m = smf.ols("u_kg ~ egfr", data=pat).fit(cov_type="cluster", cov_kwds={"groups": pat.cid})
    ci = m.conf_int().loc["egfr"]
    res["slope_per10"] = round(float(m.params["egfr"]) * 10, 5)
    res["slope_ci_per10"] = [round(float(ci[0]) * 10, 5), round(float(ci[1]) * 10, 5)]
    _p = float(m.pvalues["egfr"])
    res["slope_p"] = "p<1e-16" if _p < 1e-16 else round(_p, 4)
    mw = smf.ols("u_kg ~ C(band)", data=pat).fit(cov_type="cluster", cov_kwds={"groups": pat.cid})
    res["joint_wald_p"] = ("p<1e-16" if mw.f_pvalue < 1e-16 else f"{mw.f_pvalue:.2e}")
    res["slope_vs_policy_upperCI_below"] = {k: bool(res["slope_ci_per10"][1] < v) for k, v in POL_SLOPE10.items()}
    # band-diff bootstrap (kumele=patient)
    bd, bdci = boot_banddiff(pat)
    res["banddiff_lt15_minus_ge90"] = dict(est=bd, ci95=bdci,
        note="POZITIF=dusuk-eGFRde yuksek doz (anti-kilavuz); kumele bootstrap(hasta)")
    # F2/F3 capa + oran (per-order carpan, patient target)
    ratios = {}
    implied_band = {}   # figur icin: politika-ima dose per band
    for pol, f in POLICIES.items():
        anchor = float(pat.loc[pat.egfr.map(ANCHOR_MASK[pol]), "u_kg"].median())
        pat["_tgt"] = anchor * pat[pol]
        ra = float((pat.u_kg / pat["_tgt"]).median())              # (a) median(pratik/hedef) TERCIH
        rb = float(pat.u_kg.median() / pat["_tgt"].median())       # (b) medyan/medyan
        ratios[pol] = dict(anchor_u_kg=round(anchor, 4), ratio_a_median_of_ratio=round(ra, 3),
                           ratio_b_ratio_of_median=round(rb, 3))
        implied_band[pol] = {b: round(anchor * float(order.loc[order.band == b, pol].median()), 4) for b in BANDS}
    res["F2F3_anchor_ratio"] = ratios
    res["implied_dose_by_band"] = implied_band
    OUT[label] = res
    return pat, implied_band


def funnel_eicu():
    e = pd.read_csv("data/eicu_glargine_events.csv")
    for c in ["dose_u_kg", "egfr", "prev_dose_u_kg", "age", "weight"]:
        e[c] = pd.to_numeric(e[c], errors="coerce")
    st = {}
    st["1_order_events"] = int(len(e))
    a = e[e.age >= 18]; st["2_adult"] = int(len(a))
    w = a[a.weight.between(30, 300)]; st["3_weight"] = int(len(w))
    d = w[w.dose_u_kg.notna()]; st["4_dose_parse"] = int(len(d))
    g = d[d.egfr.notna()]; st["5_creatinine_48h"] = int(len(g))
    p = g[g.prev_dose_u_kg.notna()]; st["6_prev_dose"] = int(len(p))
    return dict(steps=st, uses=dict(F1_F4_s34="5_creatinine_48h (order); hasta-duzeyi ozetlenir",
                                    s33_titration="6_prev_dose", A1_first_dose="5 icinde is_first=1"))


def funnel_mimic():
    q = ("SELECT dose_u_kg, egfr_now, weight_kg, prev_glargine_dose_u_kg, age "
         "FROM read_parquet('data/basal_events/*.parquet') WHERE basal_type='glargine'")
    d = duckdb.connect().execute(q).df()
    st = {}
    st["1_glargine_events"] = int(len(d))
    a = d[d.dose_u_kg.notna()]; st["2_dose_nonnull"] = int(len(a))
    w = a[a.weight_kg.between(30, 300) & a.dose_u_kg.between(0.01, 3.0)]; st["3_weight_dose_sane"] = int(len(w))
    g = w[w.egfr_now.notna()]; st["4_egfr"] = int(len(g))
    p = g[g.prev_glargine_dose_u_kg.notna()]; st["5_prev_dose_titration"] = int(len(p))
    f = g[g.prev_glargine_dose_u_kg.isna()]; st["5b_first_dose"] = int(len(f))
    return dict(steps=st, note="MIMIC dolasimdaki 4 sayi -> tek huni; her analiz adima baglandi")


def figure(pats, implied, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    for ax, (label, pat) in zip(axes, pats.items()):
        data = [pat.loc[pat.band == b, "u_kg"].values for b in BANDS]
        xs = np.arange(1, 5)
        vp = ax.violinplot(data, positions=xs, showmedians=True, widths=0.8)
        for b in vp["bodies"]:
            b.set_facecolor("#9ecae1"); b.set_alpha(0.5)
        ax.boxplot(data, positions=xs, widths=0.25, showfliers=False,
                   medianprops=dict(color="black"))
        styles = {"P1_Aronoff": (":", "#d62728"), "P2_KDIGO_glargine": ("--", "#2ca02c"),
                  "P3_KDIGO_TDD": ("-.", "#9467bd")}
        for pol, (ls, col) in styles.items():
            ys = [implied[label][pol][b] for b in BANDS]
            ax.plot(xs, ys, ls=ls, color=col, lw=1.8, label=pol.replace("_", " "))
        for i, b in enumerate(BANDS):
            n = int((pat.band == b).sum())
            ax.text(xs[i], ax.get_ylim()[1] * 0.97, f"n={n}", ha="center", va="top", fontsize=8)
        r = OUT[label]
        ax.set_title(f"{label}\nslope/10mL={r['slope_per10']:+.4f} CI[{r['slope_ci_per10'][0]:+.4f},"
                     f"{r['slope_ci_per10'][1]:+.4f}] Wald {r['joint_wald_p']}", fontsize=9)
        ax.set_xticks(xs); ax.set_xticklabels([BLBL[b] for b in BANDS])
        ax.set_xlabel("eGFR band (mL/dk/1.73m2)")
    axes[0].set_ylabel("bazal glarjin U/kg (hasta-duzeyi medyan)")
    axes[0].set_ylim(0, 0.6)
    axes[0].legend(fontsize=7, loc="upper right")
    fig.suptitle("Figur 4 (v9): gozlenen bazal doz vs kilavuz-ima (per-order carpan). "
                 "Kesikli=P1/P2/P3 ima dose. KW YOK; kumelenmis Wald+slope CI.", fontsize=9)
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="svg"); plt.close(fig)


def main():
    e = load_eicu(); m = load_mimic()
    pe, ie = analyze(e, "eICU")
    pm, im = analyze(m, "MIMIC")
    OUT["F6_funnel"] = dict(eICU=funnel_eicu(), MIMIC=funnel_mimic())
    OUT["policy_thresholds_per10"] = POL_SLOPE10
    figure({"MIMIC": pm, "eICU": pe}, {"MIMIC": im, "eICU": ie}, "qc/paper2_figs/fig4_v9.svg")
    json.dump(OUT, open("qc/fig4_v9.json", "w"), indent=2)
    # ekrana ozet
    for lab in ["MIMIC", "eICU"]:
        r = OUT[lab]
        print(f"== {lab}: slope/10mL={r['slope_per10']:+.5f} CI{r['slope_ci_per10']} p={r['slope_p']} "
              f"Wald={r['joint_wald_p']} banddiff={r['banddiff_lt15_minus_ge90']}")
        print(f"   ratios: " + " | ".join(f"{k}: a={v['ratio_a_median_of_ratio']} b={v['ratio_b_ratio_of_median']} "
              f"anc={v['anchor_u_kg']}" for k, v in r["F2F3_anchor_ratio"].items()))
    print("F6 eICU funnel:", OUT["F6_funnel"]["eICU"]["steps"])
    print("F6 MIMIC funnel:", OUT["F6_funnel"]["MIMIC"]["steps"])


if __name__ == "__main__":
    main()

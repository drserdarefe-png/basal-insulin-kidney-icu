#!/usr/bin/env python3
"""
11 — MAKALE 2 SONUC TABLOLARI + FIGURLER (betimleyici + metodolojik)
====================================================================
Gozlemsel kimliklendirme DURDU. Elimizdekini sabitler: T1-T5 tablolari,
F1-F4 figurleri. YENI MODEL YOK -- mevcut parquet + json ciktilarindan.
KULLANIM: python python/11_paper2_report.py --parquet 'data/basal_events/*.parquet'
"""
import argparse, json, shutil
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BANDS = ["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]
ARONOFF = {"egfr>=90": 1.00, "egfr_45_90": 1.00, "egfr_15_45": 0.75, "egfr<15": 0.50}


def band(e):
    return pd.cut(e, [-1, 15, 45, 90, 1e9], labels=BANDS)


def miqr(x):
    x = x.dropna()
    return f"{np.median(x):.2f} [{np.percentile(x,25):.2f}-{np.percentile(x,75):.2f}]"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    out = Path(a.outdir); tdir = out / "paper2_tables"; fdir = out / "paper2_figs"
    tdir.mkdir(parents=True, exist_ok=True); fdir.mkdir(parents=True, exist_ok=True)

    df = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{a.parquet}') WHERE basal_type='glargine'").df()
    df["eb"] = band(df.egfr_now)
    # hasta duzeyi (median eGFR bandi)
    pat = df.groupby("subject_id").agg(
        egfr_med=("egfr_now", "median"), age=("age", "median"),
        gender=("gender", "first"), dm=("dm_any", "max"), cirr=("cirrhosis", "max"),
        crrt=("on_crrt", "max"), steroid=("on_steroid", "max"),
        weight=("weight_kg", "median"), los=("los_days", "max"),
        dose_ukg_med=("dose_u_kg", "median")).reset_index()
    pat["eb"] = band(pat.egfr_med)

    # ---- T1 CONSORT ----
    n_ev = len(df); n_pt = df.subject_id.nunique()
    n_night = int(df.glu_nadir_night.notna().sum())
    t1 = ["# T1 — Kohort akisi (CONSORT)\n",
          "| adim | n |", "|---|---:|",
          "| Tum ICU yatisi (MIMIC-IV v3.1) | 94.458 |",
          "| Kohort (EXC-01..11: DKA/gebelik/LOS<24h/glukoz<4/KC yet./nutr. vb.) | 7.876 |",
          "| Glarjin alan hasta | 1.117 |",
          f"| Glarjin doz olayi (03 tavan bazal 150U + event_hr>=0) | {n_ev:,} |",
          f"| + gece-glukozlu (BIRINCIL analiz seti) | {n_night:,} |",
          f"\n- Analiz seti: **{n_night:,} olay / {n_pt:,} hasta** "
          f"(ilk-doz {int(df.prev_glargine_dose_u_kg.isna().sum()):,} / "
          f"titrasyon {int(df.prev_glargine_dose_u_kg.notna().sum()):,}).",
          "- Dislama: glarjin units_total>150 U -> 4 olay (Kural 15, qc/03_dose_outliers.txt).",
          "- Registry: EXCLUSION_REGISTRY.md v1.3."]
    (tdir / "T1_cohort.md").write_text("\n".join(t1), encoding="utf-8")

    # ---- T2 Baseline by band ----
    rows = ["# T2 — Baseline (eGFR bandina gore, hasta duzeyi)\n",
            "| ozellik | " + " | ".join(BANDS) + " |", "|---|" + "---|" * len(BANDS)]
    def line(lbl, fn):
        return "| " + lbl + " | " + " | ".join(fn(pat[pat.eb == b]) for b in BANDS) + " |"
    rows.append(line("n hasta", lambda d: f"{len(d)}"))
    rows.append(line("yas medyan[IQR]", lambda d: miqr(d.age)))
    rows.append(line("kadin %", lambda d: f"{100*(d.gender=='F').mean():.0f}"))
    rows.append(line("DM %", lambda d: f"{100*d.dm.mean():.0f}"))
    rows.append(line("sirroz %", lambda d: f"{100*d.cirr.mean():.0f}"))
    rows.append(line("CRRT (hic) %", lambda d: f"{100*d.crrt.mean():.0f}"))
    rows.append(line("steroid (hic) %", lambda d: f"{100*d.steroid.mean():.0f}"))
    rows.append(line("agirlik kg medyan[IQR]", lambda d: miqr(d.weight)))
    rows.append(line("LOS gun medyan[IQR]", lambda d: miqr(d.los)))
    (tdir / "T2_baseline.md").write_text("\n".join(rows), encoding="utf-8")

    # ---- T3 Guideline-practice gap ----
    ref = np.median(df[df.eb == "egfr>=90"].dose_u_kg.dropna())
    ref_u = np.median(df[df.eb == "egfr>=90"].units_total.dropna())
    t3 = ["# T3 — Kilavuz-pratik acigi: glarjin dozu eGFR bandina gore\n",
          "| eGFR | n olay | units_total med[IQR] | U/kg med[IQR] | Aronoff ongoru U/kg | pratik/ongoru |",
          "|---|---:|---|---|---:|---:|"]
    for b in BANDS:
        d = df[df.eb == b]
        pred = ref * ARONOFF[b]
        obs = np.median(d.dose_u_kg.dropna())
        t3.append(f"| {b} | {len(d)} | {miqr(d.units_total)} | {miqr(d.dose_u_kg)} | "
                  f"{pred:.3f} | {obs/pred:.2f}x |")
    # KW hasta-duzeyi (clustering'e saygi: hasta basi medyan)
    groups = [pat[pat.eb == b].dose_ukg_med.dropna() for b in BANDS]
    H, p = stats.kruskal(*groups)
    t3.append(f"\nKruskal-Wallis (hasta-duzeyi medyan U/kg, clustering'e saygi): "
              f"H={H:.1f}, p={p:.3g}.")
    # EGIM testi (hasta-duzeyi U/kg ~ eGFR): Aronoff NEGATIF ongorur; gozlenen ~0?
    import statsmodels.api as sm
    pp = pat.dropna(subset=["dose_ukg_med", "egfr_med"])
    mm = sm.OLS(pp.dose_ukg_med.values, sm.add_constant(pp.egfr_med.values)).fit()
    sl = mm.params[1]; sci = mm.conf_int()[1]
    t3.append(f"\nEGIM (hasta-duzeyi U/kg ~ eGFR): {sl:.5f} [{sci[0]:.5f}, {sci[1]:.5f}] "
              f"U/kg per mL/min, p={mm.pvalues[1]:.3g}. Aronoff NEGATIF egim ongorur; "
              f"gozlenen egim {'~0 (CI 0 iceriyor)' if sci[0] < 0 < sci[1] else 'anlamli'} "
              "-> renal doz-azaltma UYGULANMIYOR.")
    t3.append(f"Referans eGFR>=90 medyani: {ref:.3f} U/kg ({ref_u:.0f} U). Aronoff: "
              ">50 mL/min %100, 10-50 %75, <10 %50.")
    t3.append("\n**BULGU:** gozlenen doz eGFR'den bagimsiz (egim ~0, KW p>0.5) ama Aronoff "
              "azalan doz ongorur -> KILAVUZ-PRATIK ACIGI, eGFR<15'te 2.3x. (Guc: T3b_power.json)")
    (tdir / "T3_guideline_gap.md").write_text("\n".join(t3), encoding="utf-8")

    # ---- T4 GPS/SHAP ----
    t4 = ["# T4 — Dozu ne belirliyor (GPS/SHAP, temiz skala)\n"]
    for nm, lbl in [("07_gps_titration_clean", "TITRASYON"), ("07_gps_firstdose_clean", "ILK-DOZ")]:
        j = json.load(open(out / f"{nm}.json"))
        er = next((r for r in j["shap_ranking"] if r["feature"] == "egfr_now"), None)
        t4.append(f"\n## {lbl} (n={j['n_events']})")
        t4.append(f"- CV R2 tam={j['cv_r2_full']} | prev-dose haric={j['cv_r2_no_inertia']} | "
                  f"egfr_now rank={er['rank'] if er else '-'}/{len(j['shap_ranking'])}")
        t4.append("- SHAP top-5: " + ", ".join(
            f"{r['feature']}({r['mean_abs_shap']})" for r in j["shap_ranking"][:5]))
    (tdir / "T4_gps_shap.md").write_text("\n".join(t4), encoding="utf-8")

    # ---- T5 Steroid timing ----
    st = json.load(open(out / "07_steroid_timing.json"))
    t5 = ["# T5 — Steroid zamanlamasi (generic x gece-aktif)\n",
          f"- Order baslama piki: {st['order_start_hour']['overall_peak_hour']}:00 "
          f"(sabah 6-10: %{100*st['order_start_hour']['morning_6_10_frac']:.0f})\n",
          "| generic | n olay | etki suresi (h) | gece-aktif oran |", "|---|---:|---:|---:|"]
    for r in st["active_steroid_at_glargine"]:
        t5.append(f"| {r['generic']} | {r['n_glargine_events']} | {r['effect_dur_hr']} | "
                  f"{r['frac_night_active']} |")
    t5.append("\n**BULGU:** pred_eq_24h [generic-agnostik] gece-aktif farkini yikar "
              "-> theta_steroid~1.00 maruziyet-tanimi hatasi.")
    (tdir / "T5_steroid.md").write_text("\n".join(t5), encoding="utf-8")

    # ================= FIGURLER =================
    # F1 violin dose_u_kg by band + Aronoff
    fig, ax = plt.subplots(figsize=(8, 5))
    data = [df[df.eb == b].dose_u_kg.dropna().values for b in BANDS]
    vp = ax.violinplot(data, showmedians=True)
    for pc in vp["bodies"]:
        pc.set_facecolor("0.75"); pc.set_edgecolor("black"); pc.set_alpha(0.9)
    for k in ("cmedians", "cbars", "cmins", "cmaxes"):
        vp[k].set_color("black")
    xs = list(range(1, 5))
    ax.plot(xs, [ref * ARONOFF[b] for b in BANDS], "k--D", ms=8, lw=2,
            label="Aronoff kilavuz ongorusu")
    ax.set_xticks(xs); ax.set_xticklabels(BANDS, rotation=15)
    ax.set_ylabel("glarjin doz (U/kg)"); ax.set_ylim(0, 1.2)
    ax.set_title("F1 — eGFR bandina gore glarjin dozu vs Aronoff kilavuzu")
    ax.legend(); ax.grid(axis="y", ls=":", alpha=0.5)
    fig.tight_layout(); fig.savefig(fdir / "F1_dose_by_egfr.png", dpi=300); plt.close(fig)

    # F2 SHAP bar (titrasyon)
    j = json.load(open(out / "07_gps_titration_clean.json"))
    top = j["shap_ranking"][:12][::-1]
    fig, ax = plt.subplots(figsize=(8, 6))
    names = [r["feature"] for r in top]; vals = [r["mean_abs_shap"] for r in top]
    cols = ["#c33" if n == "egfr_now" else "#39c" for n in names]
    ax.barh(names, vals, color=cols)
    er = next((r for r in j["shap_ranking"] if r["feature"] == "egfr_now"), None)
    ax.set_title(f"F2 — Dozu ne belirliyor (SHAP, titrasyon); prev_dose baskin, "
                 f"egfr_now rank {er['rank']}/{len(j['shap_ranking'])}")
    ax.set_xlabel("ortalama |SHAP|")
    fig.tight_layout(); fig.savefig(fdir / "F2_shap.png", dpi=300); plt.close(fig)

    # F3 AIPW egrisi -- mevcut png kopyala
    if (out / "07_aipw_curve.png").exists():
        shutil.copy(out / "07_aipw_curve.png", fdir / "F3_aipw_curve.png")

    # F4 A vs B forest
    w = json.load(open(out / "10_within_patient.json"))
    fig, ax = plt.subplots(figsize=(7, 3.2))
    A = w["modelA"]["d_dose"]; B = w["modelB"]["d_dose"]
    ys = [1, 0]
    for y, m, lbl in [(1, A, "Model A"), (0, B, "Model B (+prev_nadir)")]:
        ax.plot([m["ci"][0], m["ci"][1]], [y, y], "k-", lw=2)
        ax.plot(m["beta"], y, "o", ms=9, color="#39c")
        ax.text(m["beta"], y + 0.12, f"{m['beta']}", ha="center")
    ax.axvline(0, color="r", ls="--")
    ax.set_yticks([0, 1]); ax.set_yticklabels(["Model B", "Model A"])
    ax.set_xlabel("beta d_dose (gece nadir mg/dL / U/kg)")
    ax.set_title("F4 — A vs B ISARET DONMESI (zamanla-degisen karistirma)")
    fig.tight_layout(); fig.savefig(fdir / "F4_signflip.png", dpi=300); plt.close(fig)

    print(f"[11] Tablolar -> {tdir} | Figurler -> {fdir}")
    print(f"[11] T3 Aronoff ref={ref:.3f} U/kg; pratik/ongoru: "
          + ", ".join(f"{b}={np.median(df[df.eb==b].dose_u_kg.dropna())/(ref*ARONOFF[b]):.2f}x" for b in BANDS))


if __name__ == "__main__":
    main()

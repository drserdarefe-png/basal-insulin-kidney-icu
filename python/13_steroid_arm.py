#!/usr/bin/env python3
"""
13 — STEROID KOLU (Makale 2, FAZ 3). Betimleyici; KARISTIRMA AGIR (bkz. NOT).
Safak-hipo hipotezi: kisa etkili steroid (HC ~10h) sonunce, gunduze gore
ayarlanmis glarjin gece HIPO yapar. steroid_night_active=0 grubunda hipo daha sik?
Cikti: qc/paper2_tables/T10*.json, T11*.json
KULLANIM: python python/13_steroid_arm.py --parquet 'data/basal_events/*.parquet'
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy import stats

WARN = ("KARISTIRMA AGIR: steroid SECIMI endikasyona bagli (deksa->odem/onkoloji, "
        "HC->SEPTIK SOK). Septik sok hem yuksek hipo hem farkli glukoz dinamigi. "
        "'HC'de hipo sik' bulgusu steroid sonmesi mi sepsis mi AYIRT EDILEMEZ. "
        "NEDENSEL YORUM YOK.")


def load(parquet):
    df = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{parquet}') WHERE basal_type='glargine' "
        "AND dose_u_kg IS NOT NULL").df()
    df["on_steroid"] = df.on_steroid.astype(float)
    try:
        sep = pd.read_csv("/tmp/sepsis_stay.csv")
        df["sepsis3"] = df.stay_id.isin(set(sep.stay_id)).astype(int)
    except Exception:
        df["sepsis3"] = np.nan
    df["dtert"] = pd.qcut(df.dose_u_kg, 3, labels=["D1_dusuk", "D2", "D3_yuksek"])
    return df


def mi(x):
    x = pd.Series(x).dropna()
    return [round(float(np.median(x)), 3), round(float(np.percentile(x, 25)), 3),
            round(float(np.percentile(x, 75)), 3)] if len(x) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    T = Path(a.outdir) / "paper2_tables"; T.mkdir(parents=True, exist_ok=True)
    df = load(a.parquet)
    act = df[df.on_steroid == 1]

    # ---- T10a fizibilite ----
    gen = {}
    for g, d in act.groupby("steroid_generic"):
        gen[str(g)] = dict(n_events=int(len(d)), n_subjects=int(d.subject_id.nunique()))
    t10a = dict(steroid_active_events=int(len(act)), steroid_active_subjects=int(act.subject_id.nunique()),
                by_generic=gen,
                night_active_0=int((act.steroid_night_active == 0).sum()),
                night_active_1=int((act.steroid_night_active == 1).sum()),
                sepsis3_flag_n_stays=int(df.sepsis3.sum()) if df.sepsis3.notna().any() else None,
                sepsis3_note="Sepsis-3 flag muhtemelen EKSIK-yakalanmis (~100 stay); sinirli guc.")
    small = [g for g, v in gen.items() if v["n_subjects"] < 20]
    t10a["DUR_small_generics"] = small or None
    (T / "T10a_steroid_feasibility.json").write_text(json.dumps(t10a, indent=2))

    # ---- T10 dozlama: steroid alan doz artiriyor mu ----
    non = df[df.on_steroid == 0]
    t10 = dict(question="Klinisyen steroid alanda glarjin dozunu artiriyor mu?",
               steroid1_u_kg=mi(act.dose_u_kg), steroid0_u_kg=mi(non.dose_u_kg),
               steroid1_pred_eq_active_med=round(float(np.median(act.pred_eq_active.dropna())), 1),
               by_generic={g: mi(d.dose_u_kg) for g, d in act.groupby("steroid_generic")},
               kw_u_kg_steroid_p=round(float(stats.mannwhitneyu(
                   act.dose_u_kg.dropna(), non.dose_u_kg.dropna()).pvalue), 4))
    # pred_eq tertil
    a2 = act.dropna(subset=["pred_eq_active"]).copy()
    if len(a2) > 30:
        a2["pt"] = pd.qcut(a2.pred_eq_active, 3, labels=["P1", "P2", "P3"], duplicates="drop")
        t10["pred_eq_tertile_u_kg"] = {str(t): mi(a2[a2.pt == t].dose_u_kg) for t in a2.pt.dropna().unique()}
    (T / "T10_steroid_dosing.json").write_text(json.dumps(t10, indent=2))

    # ---- T11 safak hipo: night_active x hipo, DOZ TERTILINE gore ----
    t11 = {"NOT": WARN, "ongoru": "night_active=0 (steroid sondu) -> hipo DAHA SIK (doz gunduze gore ayarli)",
           "by_dose_tertile": {}}
    for dt in ["D1_dusuk", "D2", "D3_yuksek"]:
        sub = act[act.dtert == dt]
        row = {}
        for na in (0, 1):
            s = sub[sub.steroid_night_active == na]
            row[f"night_active={na}"] = dict(n=int(len(s)),
                hypo_night_l1=round(float(s.hypo_night_l1.mean()), 4) if len(s) else None)
        t11["by_dose_tertile"][dt] = row
    # sepsis-stratifiye [kismi kontrol]
    if df.sepsis3.notna().any():
        t11["by_sepsis"] = {}
        for sp in (0, 1):
            s = act[act.sepsis3 == sp]
            t11["by_sepsis"][f"sepsis={sp}"] = {
                f"night_active={na}": round(float(s[s.steroid_night_active == na].hypo_night_l1.mean()), 4)
                if len(s[s.steroid_night_active == na]) else None for na in (0, 1)}
    (T / "T11_dawn_hypo.json").write_text(json.dumps(t11, indent=2))

    # ---- T11b jenerik x hipo, doz tertiline gore ----
    t11b = {"NOT": WARN, "ongoru": "HC (%15.8 gece-aktif) > deksa (%56.7) hipo",
            "by_generic_overall": {}, "by_generic_x_tertile": {}}
    for g, d in act.groupby("steroid_generic"):
        t11b["by_generic_overall"][str(g)] = dict(n=int(len(d)),
            hypo_night_l1=round(float(d.hypo_night_l1.mean()), 4),
            night_active_frac=round(float(d.steroid_night_active.mean()), 3))
    for dt in ["D1_dusuk", "D2", "D3_yuksek"]:
        sub = act[act.dtert == dt]
        t11b["by_generic_x_tertile"][dt] = {str(g): round(float(d.hypo_night_l1.mean()), 4)
            for g, d in sub.groupby("steroid_generic") if len(d) >= 10}
    (T / "T11b_generic_hypo.json").write_text(json.dumps(t11b, indent=2))

    print(f"[13] T10a steroid aktif: {len(act)} olay / {act.subject_id.nunique()} hasta; "
          f"kucuk jenerik(<20 subj): {small}")
    print(f"[13] T10 steroid1 u_kg={t10['steroid1_u_kg'][0]} vs steroid0={t10['steroid0_u_kg'][0]} "
          f"MWU p={t10['kw_u_kg_steroid_p']}")
    for dt, r in t11["by_dose_tertile"].items():
        print(f"[13] T11 {dt}: na=0 hypo={r['night_active=0']['hypo_night_l1']} "
              f"(n{r['night_active=0']['n']}) vs na=1 hypo={r['night_active=1']['hypo_night_l1']} "
              f"(n{r['night_active=1']['n']})")


if __name__ == "__main__":
    main()

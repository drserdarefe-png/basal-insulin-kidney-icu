#!/usr/bin/env python3
"""
14 — DIYABET DURUMU KOLU (Makale 2, FAZ 3). Betimleyici.
Uc fenotip: T1DM (tampon yok) / T2DM (frenler) / dm_any=0 STRES HIPERGLISEMISI
(pankreas normal, stres gecince devralir -> glarjin devam ederse hipo).
SINIRLILIK: dm_t1 ICD-tabanli, T1/T2 ayriminda ZAYIF (C-peptid yok).
Cikti: qc/paper2_tables/T12*.json, T13.json, T14.json
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy import stats

# T12a BQ sonuclari (qc/_tmp_dmdef.sql):
CPEP_N_HADM = 0        # C-peptid kohortta yok
N_HADM_T1 = 54
N_T1_AND_T2 = 1


def mi(x):
    x = pd.Series(x).dropna()
    return [round(float(np.median(x)), 3), round(float(np.percentile(x, 25)), 3),
            round(float(np.percentile(x, 75)), 3)] if len(x) else None


def load(parquet):
    df = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{parquet}') WHERE basal_type='glargine' "
        "AND dose_u_kg IS NOT NULL").df()
    for b in ("dm_t1", "dm_any", "on_steroid"):
        df[b] = df[b].astype(float)
    df["dm_status"] = np.where(df.dm_t1 == 1, "T1DM",
                        np.where(df.dm_any == 1, "T2DM", "stres_hiperglisemi"))
    df["eb"] = pd.cut(df.egfr_now, [-1, 15, 45, 90, 1e9],
                      labels=["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"])
    df["dtert"] = pd.qcut(df.dose_u_kg, 3, labels=["D1_dusuk", "D2", "D3_yuksek"])
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True); ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    T = Path(a.outdir) / "paper2_tables"; T.mkdir(parents=True, exist_ok=True)
    df = load(a.parquet)
    t1 = df[df.dm_status == "T1DM"]

    # ---- T12a tanim dogrulama ----
    t12a = dict(dm_t1_definition="ICD REGEXP ^E10|^2500[13]",
                c_peptide_n_hadm_cohort=CPEP_N_HADM,
                c_peptide_verdict="C-peptid kohortta YOK (<100) -> tanim ICD'de kalir, SINIRLILIK.",
                n_hadm_t1_cohort=N_HADM_T1, n_hadm_t1_and_t2_contradiction=N_T1_AND_T2,
                dm_t1_glargine_subjects=int(t1.subject_id.nunique()),
                dm_t1_age_med_iqr=mi(t1.groupby("subject_id").age.first()),
                age_note="Gercek T1DM genc-baslangicli; medyan yas >60 ise ICD T1'i uzun-sureli "
                         "insulin-bagimli T2DM'yi yakaliyor olabilir.")
    (T / "T12a_dm_definition.json").write_text(json.dumps(t12a, indent=2))

    # ---- T12b fizibilite ----
    fb = {}
    for k in ["T1DM", "T2DM", "stres_hiperglisemi"]:
        d = df[df.dm_status == k]
        fb[k] = dict(n_events=int(len(d)), n_subjects=int(d.subject_id.nunique()))
    fb["DUR_dm_t1_small"] = ("dm_t1 subj<20" if fb["T1DM"]["n_subjects"] < 20 else None)
    fb["note_06"] = "06 positivity'de dm_t1=1 sadece 39 hasta ile GECEMEDI -- zayif kol."
    (T / "T12b_dm_feasibility.json").write_text(json.dumps(fb, indent=2))

    # ---- T12 dozlama ----
    t12 = dict(question="Klinisyen diyabet fenotipini ayirt ediyor mu? Ongoru: T1DM farkli.",
               by_dm_status={k: dict(n=int((df.dm_status == k).sum()),
                                     u_kg_med_iqr=mi(df[df.dm_status == k].dose_u_kg),
                                     units_total_med_iqr=mi(df[df.dm_status == k].units_total))
                             for k in ["T1DM", "T2DM", "stres_hiperglisemi"]})
    patm = df.groupby("subject_id").agg(u=("dose_u_kg", "median"), s=("dm_status", "first"))
    t12["kw_u_kg_p"] = round(float(stats.kruskal(
        *[patm[patm.s == k].u.dropna() for k in ["T1DM", "T2DM", "stres_hiperglisemi"]])[1]), 4)
    (T / "T12_dm_dosing.json").write_text(json.dumps(t12, indent=2))

    # ---- T13 hipo, doz tertiline gore stratifiye ----
    t13 = {"NOT": "BETIMLEYICI -- SEBEP DEGIL. Ongoru: ayni U/kg dozda T1DM>T2DM.",
           "by_dose_tertile": {}}
    for dt in ["D1_dusuk", "D2", "D3_yuksek"]:
        sub = df[df.dtert == dt]
        t13["by_dose_tertile"][dt] = {k: dict(n=int((sub.dm_status == k).sum()),
            hypo_night_l1=round(float(sub[sub.dm_status == k].hypo_night_l1.mean()), 4)
            if (sub.dm_status == k).any() else None) for k in ["T1DM", "T2DM", "stres_hiperglisemi"]}
    t13["stres_hiperglisemi_overall"] = dict(
        n=int((df.dm_status == "stres_hiperglisemi").sum()),
        hypo_night_l1=round(float(df[df.dm_status == "stres_hiperglisemi"].hypo_night_l1.mean()), 4),
        note="Glarjin baslanmasi EN TARTISMALI grup.")
    (T / "T13_hypo_dm.json").write_text(json.dumps(t13, indent=2))

    # ---- T14 etkilesim ----
    both = df[(df.dm_t1 == 1) & (df.eb == "egfr<15")]
    stress_ster = df[(df.dm_status == "stres_hiperglisemi") & (df.on_steroid == 1)]
    t14 = dict(dm_t1_x_egfr_lt15_subjects=int(both.subject_id.nunique()),
               dm_t1_x_egfr_lt15_events=int(len(both)),
               stress_x_steroid_events=int(len(stress_ster)),
               stress_x_steroid_subjects=int(stress_ster.subject_id.nunique()),
               note="Stres-hiperglisemi + steroid: iyatrojenik; glarjin haklilastirilir mi? "
                    "dm_t1 x eGFR<15 muhtemelen cok az (n<20 -> yorum yok).")
    (T / "T14_dm_interaction.json").write_text(json.dumps(t14, indent=2))

    print(f"[14] T12b: T1DM {fb['T1DM']['n_subjects']}h / T2DM {fb['T2DM']['n_subjects']}h / "
          f"stres {fb['stres_hiperglisemi']['n_subjects']}h")
    print(f"[14] T12 u_kg: T1={t12['by_dm_status']['T1DM']['u_kg_med_iqr'][0]} "
          f"T2={t12['by_dm_status']['T2DM']['u_kg_med_iqr'][0]} "
          f"stres={t12['by_dm_status']['stres_hiperglisemi']['u_kg_med_iqr'][0]} KW p={t12['kw_u_kg_p']}")
    print(f"[14] T13 stres hipo_l1={t13['stres_hiperglisemi_overall']['hypo_night_l1']}; "
          f"T12a dm_t1 yas={t12a['dm_t1_age_med_iqr']}")
    print(f"[14] T14 dm_t1xeGFR<15 subj={t14['dm_t1_x_egfr_lt15_subjects']}, stresxsteroid subj={t14['stress_x_steroid_subjects']}")


if __name__ == "__main__":
    main()

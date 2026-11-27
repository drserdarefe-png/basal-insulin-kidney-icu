#!/usr/bin/env python3
"""27 — F10 son sayilar: eGFR<10 n (event+patient), basamak-CI politika DISLAMA, capalar.
Girdi: data/eicu_glargine_events.csv, MIMIC parquet + qc/fig4_f9.json. Cikti: qc/fig4_f10.json."""
import json
import numpy as np, pandas as pd, duckdb

f9 = json.load(open("qc/fig4_f9.json"))
OUT = {}


def egfr_counts():
    e = pd.read_csv("data/eicu_glargine_events.csv")
    for c in ["dose_u_kg", "egfr", "weight", "age"]:
        e[c] = pd.to_numeric(e[c], errors="coerce")
    e = e[(e.age >= 18) & e.weight.between(30, 300) & e.dose_u_kg.between(0.01, 3.0) & e.egfr.notna()]
    q = ("SELECT subject_id, egfr_now AS egfr, dose_u_kg, weight_kg AS w "
         "FROM read_parquet('data/basal_events/*.parquet') "
         "WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL AND egfr_now IS NOT NULL")
    m = duckdb.connect().execute(q).df()
    m = m[m.w.between(30, 300) & m.dose_u_kg.between(0.01, 3.0)]
    OUT["egfr_lt10_n"] = dict(
        eICU=dict(events=int((e.egfr < 10).sum()),
                  patients=int(e.loc[e.egfr < 10, "patientunitstayid"].nunique())),
        MIMIC=dict(events=int((m.egfr < 10).sum()),
                   patients=int(m.loc[m.egfr < 10, "subject_id"].nunique())),
        note="item3 oran (2.01/1.61) HASTA-duzeyi; buradaki patients ile ayni.")


def exclusion():
    tab = {}
    for coh in ["eICU", "MIMIC"]:
        s = f9["item2_step_threshold"][coh]
        st = s["step_lt60"]
        lo, hi = st["ci95"]
        pol60 = {"P2_m30": st["guideline_P2_step_m30"], "P3_m25": st["guideline_P3_step_m25"]}
        excl60 = {k: bool(v < lo or v > hi) for k, v in pol60.items()}
        c10 = s["aronoff_3cat"]["c_lt10"]
        lo10, hi10 = c10["ci"]
        pol10 = c10["guideline_step"]  # P1 Aronoff -50%
        tab[coh] = dict(
            anchor_ge60=st["anchor_ge60"], anchor_gt50=s["anchor50"],
            step_lt60=dict(obs=st["observed_step_lt60"], ci=[lo, hi], p=st["p"],
                           policy=pol60, excludes=excl60,
                           anchor_note="-0.045 = P3(KDIGO-TDD -%25) x capa~0.178 -> TEYIT"),
            step_lt10_P1=dict(obs=c10["observed"], ci=[lo10, hi10], n_patients=c10["n"], p=c10["p"],
                              policy_P1_m50=pol10, excludes=bool(pol10 < lo10 or pol10 > hi10)))
    OUT["step_policy_exclusion"] = tab
    OUT["note_exclusion"] = ("Politika adimi CI'nin ALTINDA (daha negatif) -> DISLANIR = gozlenen "
                             "azaltma kilavuzdan anlamli KUCUK. p(vs0) degil, CI(vs politika) hukum verir.")


def ratio_overlap():
    e10 = f9["item3_bandspecific_ratio"]["eICU"]["by_policy"]["P1_Aronoff"]["by_egfr_bin"]["lt10"]
    m10 = f9["item3_bandspecific_ratio"]["MIMIC"]["by_policy"]["P1_Aronoff"]["by_egfr_bin"]["lt10"]
    overlap = not (e10["ci"][1] < m10["ci"][0] or m10["ci"][1] < e10["ci"][0])
    OUT["egfr_lt10_ratio_P1"] = dict(eICU=e10, MIMIC=m10, ci_overlap=bool(overlap),
        verdict=("CI'lar ORTUSUYOR -> GERCEK REPLIKASYON (iki merkez portfoyu ayni yonu/buyuklugu)"
                 if overlap else "CI ayrik"))


egfr_counts(); exclusion(); ratio_overlap()
json.dump(OUT, open("qc/fig4_f10.json", "w"), indent=2)
print(json.dumps(OUT, indent=2, ensure_ascii=False))

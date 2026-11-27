#!/usr/bin/env python3
"""
12 — MAKALE 2 EK ANALIZLER (T3b guc, T3c neg-kontrol, T6 hipo; hepatik T7-T9)
============================================================================
Betimleyici; nedensel iddia YOK. Cikti: qc/paper2_tables/*.json
KULLANIM: python python/12_paper2_supp.py --parquet 'data/basal_events/*.parquet'
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy import stats
from statsmodels.stats.power import TTestIndPower

BANDS = ["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]


def eband(e): return pd.cut(e, [-1, 15, 45, 90, 1e9], labels=BANDS)


def boot_med_diff(a_subj, a_val, b_subj, b_val, B=2000, seed=0):
    """iki grup medyan farki, subject-cluster bootstrap %95 CI."""
    a_subj, a_val = np.asarray(a_subj), np.asarray(a_val, float)
    b_subj, b_val = np.asarray(b_subj), np.asarray(b_val, float)
    ma, mb = np.isfinite(a_val), np.isfinite(b_val)
    a_subj, a_val, b_subj, b_val = a_subj[ma], a_val[ma], b_subj[mb], b_val[mb]
    rng = np.random.default_rng(seed)
    def one(subj, val):
        g = pd.DataFrame({"s": subj, "v": val}).groupby("s")["v"].apply(list)
        return g.index.values, g.values
    sa, va = one(a_subj, a_val); sb, vb = one(b_subj, b_val)
    diffs = []
    for _ in range(B):
        ra = np.concatenate([va[i] for i in rng.integers(0, len(va), len(va))])
        rb = np.concatenate([vb[i] for i in rng.integers(0, len(vb), len(vb))])
        diffs.append(np.median(ra) - np.median(rb))
    return float(np.median(a_val) - np.median(b_val)), \
        float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def load(parquet):
    df = duckdb.connect().execute(
        f"SELECT * FROM read_parquet('{parquet}') WHERE basal_type='glargine'").df()
    for b in ("cirrhosis", "dm_t1", "dm_any", "on_crrt"):
        df[b] = df[b].astype("float")
    df["eb"] = eband(df.egfr_now)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    T = Path(a.outdir) / "paper2_tables"; T.mkdir(parents=True, exist_ok=True)
    df = load(a.parquet)
    out = Path(a.outdir)

    # ===== T3b: etki buyuklugu + guc =====
    lo = df[df.eb == "egfr<15"]; hi = df[df.eb == "egfr>=90"]
    d, cl, cu = boot_med_diff(lo.subject_id.values, lo.dose_u_kg.values,
                              hi.subject_id.values, hi.dose_u_kg.values)
    refmed = float(np.nanmedian(hi.dose_u_kg))
    sd = float(np.nanstd(pd.concat([lo.dose_u_kg, hi.dose_u_kg])))
    eff50 = 0.5 * refmed            # tespit edilmek istenen %50 fark
    power = float(TTestIndPower().power(effect_size=eff50 / sd,
                  nobs1=lo.subject_id.nunique(), alpha=0.05,
                  ratio=hi.subject_id.nunique() / lo.subject_id.nunique()))
    t3b = dict(
        comparison="eGFR<15 vs eGFR>=90, dose_u_kg",
        median_diff=round(d, 4), ci95=[round(cl, 4), round(cu, 4)],
        n_subj_lo=int(lo.subject_id.nunique()), n_subj_hi=int(hi.subject_id.nunique()),
        target_50pct_effect_u_kg=round(eff50, 4), pooled_sd=round(sd, 4),
        posthoc_power_detect_50pct=round(power, 3),
        verdict=("GUC YETERLI (>=0.8): %50 fark tespit edilebilirdi -> 'fark YOK' denebilir"
                 if power >= 0.8 else
                 "GUC DUSUK: %50 fark tespit edilemeyebilirdi -> 'gosteremedik'"))
    (T / "T3b_power.json").write_text(json.dumps(t3b, indent=2))

    # ===== T3c: negatif kontrol -- agirlik ve DM tipi =====
    df["wband"] = pd.cut(df.weight_kg, [0, 70, 90, 1e9], labels=["<70", "70-90", ">90"])
    wtab = {str(b): dict(n=int((df.wband == b).sum()),
                         units_total_med=round(float(np.median(df[df.wband == b].units_total.dropna())), 1),
                         u_kg_med=round(float(np.median(df[df.wband == b].dose_u_kg.dropna())), 3))
            for b in ["<70", "70-90", ">90"]}
    df["dmtype"] = np.where(df.dm_t1 == 1, "T1DM",
                    np.where(df.dm_any == 1, "T2DM", "nonDM"))
    dtab = {k: dict(n=int((df.dmtype == k).sum()),
                    units_total_med=round(float(np.median(df[df.dmtype == k].units_total.dropna())), 1),
                    u_kg_med=round(float(np.median(df[df.dmtype == k].dose_u_kg.dropna())), 3))
            for k in ["T1DM", "T2DM", "nonDM"]}
    # KW hasta-duzeyi
    patw = df.groupby("subject_id").agg(u=("units_total", "median"), wb=("wband", "first")).dropna()
    Hw, pw = stats.kruskal(*[patw[patw.wb == b].u for b in ["<70", "70-90", ">90"] if (patw.wb == b).any()])
    t3c = dict(weight_bands=wtab, kw_units_by_weight_p=round(float(pw), 4),
               dm_type=dtab,
               interpretation=("Klinisyen AGIRLIGA gore ayarliyorsa units_total agirlikla artar; "
                               "T1DM daha yuksek doz almali. eGFR'ye gore ayarlamiyor (T3) ama "
                               "bunlara gore ayarliyorsa -> 'ayarlayabiliyorlar, eGFR icin degil'."))
    (T / "T3c_negcontrol.json").write_text(json.dumps(t3c, indent=2))

    # ===== T6: betimleyici hipoglisemi (eGFR) =====
    t6 = {"NOT": "BETIMLEYICI -- nedensel iddia YOK. 'daha sik goruluyor', 'doz yuzunden' DEGIL.",
          "by_egfr": {}}
    for b in BANDS:
        d_ = df[df.eb == b]
        t6["by_egfr"][b] = dict(n=int(len(d_)),
            hypo_night_l1=round(float(d_.hypo_night_l1.mean()), 4),
            hypo_night_l2=round(float(d_.hypo_night_l2.mean()), 4))
    (T / "T6_hypo_by_egfr.json").write_text(json.dumps(t6, indent=2))

    # ===== T7a: hepatik fizibilite =====
    cir = df[df.cirrhosis == 1]
    ncs = int(cir.subject_id.nunique())
    cp = cir.dropna(subset=["child_pugh_score"])
    def cpclass(s): return "A" if s <= 6 else ("B" if s <= 9 else "C")
    cp_cls = cp.assign(cls=cp.child_pugh_score.map(cpclass)) if len(cp) else cp
    t7a = dict(cirrhosis_n_events=int(len(cir)), cirrhosis_n_subjects=ncs,
               child_pugh_available_events=int(len(cp)),
               child_pugh_available_subjects=int(cp.subject_id.nunique()) if len(cp) else 0,
               meld_na_pctl=[round(float(x), 1) for x in np.percentile(cir.meld_na.dropna(), [25, 50, 75])]
                            if cir.meld_na.notna().any() else None,
               child_pugh_class_subjects=({c: int(cp_cls[cp_cls.cls == c].subject_id.nunique())
                                           for c in ["A", "B", "C"]} if len(cp) else {}))
    t7a["DUR"] = (f"cirrhosis n_subjects={ncs} < 20 -> tablo kurma" if ncs < 20 else None)
    (T / "T7a_hepatic_feasibility.json").write_text(json.dumps(t7a, indent=2))
    if ncs < 20:
        print("[12] T7a DUR:", t7a["DUR"]); return

    # ===== T7: hepatik dozlama =====
    def mi(x): x = x.dropna(); return [round(float(np.median(x)), 3),
        round(float(np.percentile(x, 25)), 3), round(float(np.percentile(x, 75)), 3)]
    non = df[df.cirrhosis == 0]
    dU, clU, cuU = boot_med_diff(cir.subject_id.values, cir.dose_u_kg.values,
                                 non.subject_id.values, non.dose_u_kg.values)
    t7 = dict(frame="Karacigerde SAYISAL kural YOK (renaldeki Aronoff gibi degil) -> "
                    "'hic azaltiyorlar mi' soruluyor; referanssiz, gozlenen fark+CI.",
              cirrhosis1=dict(n=int(len(cir)), units_total_med_iqr=mi(cir.units_total),
                              u_kg_med_iqr=mi(cir.dose_u_kg)),
              cirrhosis0=dict(n=int(len(non)), units_total_med_iqr=mi(non.units_total),
                              u_kg_med_iqr=mi(non.dose_u_kg)),
              u_kg_median_diff_cirr_minus_noncirr=round(dU, 4), ci95=[round(clU, 4), round(cuU, 4)])
    # Child-Pugh kirilimi
    t7["child_pugh"] = {c: dict(n=int((cp_cls.cls == c).sum()),
                                u_kg_med_iqr=mi(cp_cls[cp_cls.cls == c].dose_u_kg))
                        for c in ["A", "B", "C"] if (cp_cls.cls == c).any()}
    # MELD tertil
    cir2 = cir.dropna(subset=["meld_na"]).copy()
    if len(cir2) > 30:
        cir2["mt"] = pd.qcut(cir2.meld_na, 3, labels=["T1", "T2", "T3"], duplicates="drop")
        t7["meld_tertile"] = {str(t): dict(n=int((cir2.mt == t).sum()),
                                           u_kg_med_iqr=mi(cir2[cir2.mt == t].dose_u_kg))
                              for t in cir2.mt.dropna().unique()}
    (T / "T7_hepatic_dosing.json").write_text(json.dumps(t7, indent=2))

    # ===== T7b: cirrhosis SHAP rank =====
    g = json.load(open(out / "07_gps_titration_clean.json"))
    cr = next((r for r in g["shap_ranking"] if r["feature"] == "cirrhosis"), None)
    t7b = dict(source="07_gps_titration_clean.json",
               cirrhosis_shap_rank=(cr["rank"] if cr else None),
               cirrhosis_shap=(cr["mean_abs_shap"] if cr else None),
               n_features=len(g["shap_ranking"]),
               note="egfr_now rank ~11-17. Dusuk rank -> klinisyen sirroza gore de az ayarliyor.")
    (T / "T7b_cirrhosis_shap.json").write_text(json.dumps(t7b, indent=2))

    # ===== T8: betimleyici hipoglisemi (hepatik) =====
    t8 = {"NOT": "BETIMLEYICI -- SEBEP DEGIL.",
          "cirrhosis1": dict(n=int(len(cir)), hypo_night_l1=round(float(cir.hypo_night_l1.mean()), 4),
                             hypo_night_l2=round(float(cir.hypo_night_l2.mean()), 4)),
          "cirrhosis0": dict(n=int(len(non)), hypo_night_l1=round(float(non.hypo_night_l1.mean()), 4),
                             hypo_night_l2=round(float(non.hypo_night_l2.mean()), 4)),
          "child_pugh": {c: dict(n=int((cp_cls.cls == c).sum()),
                                 hypo_night_l1=round(float(cp_cls[cp_cls.cls == c].hypo_night_l1.mean()), 4))
                         for c in ["A", "B", "C"] if (cp_cls.cls == c).any()},
          "vs_egfr_lt15_hypo_l1": round(float(df[df.eb == "egfr<15"].hypo_night_l1.mean()), 4)}
    (T / "T8_hypo_hepatic.json").write_text(json.dumps(t8, indent=2))

    # ===== T9: cirrhosis x eGFR<15 (hepatorenal) =====
    both = df[(df.cirrhosis == 1) & (df.eb == "egfr<15")]
    t9 = dict(cirr_and_egfr_lt15_events=int(len(both)),
              cirr_and_egfr_lt15_subjects=int(both.subject_id.nunique()),
              note="Hepatorenal sinyali; muhtemelen az. n<20 ise yorum yapma.")
    (T / "T9_hepatorenal.json").write_text(json.dumps(t9, indent=2))

    print(f"[12] T3b power={t3b['posthoc_power_detect_50pct']} diff={t3b['median_diff']}{t3b['ci95']}")
    print(f"[12] T3c: units by weight <70/70-90/>90 = "
          + "/".join(str(wtab[b]['units_total_med']) for b in ['<70', '70-90', '>90'])
          + f" (KW p={t3c['kw_units_by_weight_p']}); DM T1/T2/non u_kg = "
          + "/".join(str(dtab[k]['u_kg_med']) for k in ['T1DM', 'T2DM', 'nonDM']))
    print(f"[12] T6 hypo_l1 eGFR<15={t6['by_egfr']['egfr<15']['hypo_night_l1']} vs >=90={t6['by_egfr']['egfr>=90']['hypo_night_l1']}")
    print(f"[12] T7a cirrhosis subj={ncs} CP-class={t7a['child_pugh_class_subjects']}")
    print(f"[12] T7 u_kg cirr-noncirr diff={t7['u_kg_median_diff_cirr_minus_noncirr']}{t7['ci95']}")
    print(f"[12] T7b cirrhosis SHAP rank={t7b['cirrhosis_shap_rank']}/{t7b['n_features']}")
    print(f"[12] T8 hypo_l1 cirr={t8['cirrhosis1']['hypo_night_l1']} vs noncirr={t8['cirrhosis0']['hypo_night_l1']} vs eGFR<15={t8['vs_egfr_lt15_hypo_l1']}")
    print(f"[12] T9 hepatorenal n_subj={t9['cirr_and_egfr_lt15_subjects']}")


if __name__ == "__main__":
    main()

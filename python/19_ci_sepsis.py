#!/usr/bin/env python3
"""
19 — Wilson CI + T11c sepsis 3-yonlu stratifikasyon (Makale 2 revizyon).
BETIMLEYICI. Kaynak = parquet (Kural 16). Sepsis3 = /tmp/sepsis_stay.csv uyeligi.
Cikti:
  qc/paper2_tables/T8_ci.json       (cirrhosis/non-cirr/eGFR<15 hipo + Wilson CI)
  qc/paper2_tables/T11_ci.json      (T11 hucreleri + Wilson CI)
  qc/paper2_tables/T11c_sepsis_strat.json  (sepsis x night_active x doz-tertili)
KULLANIM: .venv/bin/python python/19_ci_sepsis.py --parquet 'data/basal_events/*.parquet'
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb

Z = 1.959963985  # %95


def wilson(k, n):
    if n == 0:
        return None
    p = k / n
    z2 = Z * Z
    center = (k + z2 / 2) / (n + z2)
    half = Z * np.sqrt((k * (n - k) / n + z2 / 4)) / (n + z2)
    return dict(n=int(n), k=int(k), p=round(p, 4),
                lo=round(max(0.0, center - half), 4),
                hi=round(min(1.0, center + half), 4))


def cell(s):
    """s: pandas Series of 0/1 hypo flags."""
    s = s.dropna()
    return wilson(int(s.sum()), int(len(s)))


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
    # GLOBAL doz tertili (13_steroid_arm.py ile AYNI: tum glarjin olaylari uzerinde)
    df["dtert"] = pd.qcut(df.dose_u_kg, 3, labels=["D1_dusuk", "D2", "D3_yuksek"])
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    T = Path(a.outdir) / "paper2_tables"; T.mkdir(parents=True, exist_ok=True)
    df = load(a.parquet)

    # ---------- T8 CI (hepatik hipo) ----------
    t8 = {"NOT": "BETIMLEYICI. n kucuk, guc sinirli; 'sirotikte hipo AZ' DIYE YAZILMAZ.",
          "cirrhosis1_l1": cell(df[df.cirrhosis == 1].hypo_night_l1),
          "cirrhosis0_l1": cell(df[df.cirrhosis == 0].hypo_night_l1),
          "egfr_lt15_l1": cell(df[df.egfr_now < 15].hypo_night_l1)}
    (T / "T8_ci.json").write_text(json.dumps(t8, indent=2))

    # ---------- T11 CI (safak hipo, steroid-aktif) ----------
    act = df[df.on_steroid == 1]
    t11 = {"NOT": "Wilson %95 CI. CI'lar ortusuyorsa 'hipotez destekli' DEME.",
           "by_dose_tertile": {}}
    for dt in ["D1_dusuk", "D2", "D3_yuksek"]:
        sub = act[act.dtert == dt]
        t11["by_dose_tertile"][dt] = {
            f"night_active={na}": cell(sub[sub.steroid_night_active == na].hypo_night_l1)
            for na in (0, 1)}
    (T / "T11_ci.json").write_text(json.dumps(t11, indent=2))

    # ---------- T11c: sepsis x night_active x doz-tertili (KONTROL) ----------
    t11c = {"NOT": ("KISMI KONTROL. Sepsis-3 flag EKSIK-yakalanmis olabilir (uyeligi "
                    "/tmp/sepsis_stay.csv). Desen sepsis icinde de duruyorsa karistirma "
                    "uyarisi ZAYIFLAR; kayboluyorsa karistirma DOGRULANIR."),
            "sepsis3_available": bool(df.sepsis3.notna().any()),
            "sepsis3_n_stays": int(df.loc[df.sepsis3 == 1, "stay_id"].nunique()),
            "by_sepsis_tertile": {}}
    for sp in (0, 1):
        s_sp = act[act.sepsis3 == sp]
        block = {}
        for dt in ["D1_dusuk", "D2", "D3_yuksek"]:
            sub = s_sp[s_sp.dtert == dt]
            block[dt] = {f"night_active={na}": cell(sub[sub.steroid_night_active == na].hypo_night_l1)
                         for na in (0, 1)}
        t11c["by_sepsis_tertile"][f"sepsis={sp}"] = block
    # sepsis x night_active (tertilsiz, guc icin)
    t11c["by_sepsis_pooled"] = {}
    for sp in (0, 1):
        s_sp = act[act.sepsis3 == sp]
        t11c["by_sepsis_pooled"][f"sepsis={sp}"] = {
            f"night_active={na}": cell(s_sp[s_sp.steroid_night_active == na].hypo_night_l1)
            for na in (0, 1)}
    (T / "T11c_sepsis_strat.json").write_text(json.dumps(t11c, indent=2))

    # ---- ekrana ----
    def fmt(c):
        return f"n={c['n']} k={c['k']} p={c['p']} CI[{c['lo']},{c['hi']}]" if c else "n=0"
    print("[19] T8 cirr1:", fmt(t8["cirrhosis1_l1"]), "| cirr0:", fmt(t8["cirrhosis0_l1"]),
          "| eGFR<15:", fmt(t8["egfr_lt15_l1"]))
    for dt, r in t11["by_dose_tertile"].items():
        print(f"[19] T11 {dt}: na0 {fmt(r['night_active=0'])} || na1 {fmt(r['night_active=1'])}")
    for sp, block in t11c["by_sepsis_tertile"].items():
        print(f"[19] T11c {sp} pooled: na0 {fmt(t11c['by_sepsis_pooled'][sp]['night_active=0'])} "
              f"|| na1 {fmt(t11c['by_sepsis_pooled'][sp]['night_active=1'])}")
        for dt, r in block.items():
            print(f"        {dt}: na0 {fmt(r['night_active=0'])} || na1 {fmt(r['night_active=1'])}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
06 — POSITIVITY / COMMON-SUPPORT TANISI                        v1.0
===================================================================
DURUM: Bu script depoda YOKTU; gece-modu oturumunda YAZILDI (Claude).
  Positivity mantigini S. Efe DOGRULAMALI (esikler + stratum tanimlari).
  07/08 karar kuralinin 'in_support' girdisini ve model_card_support.md'yi
  bu script uretir (bkz. 08_dose_optimizer.py sat. 50-51, 81-84).

AMAC (Kural 4 cercevesinde nedensel dogrulama):
  Surekli-doz ISF modeli d_glu_t2h = b0 - ISF*doz + ... egimini ancak dozun
  DEGISKEN oldugu (deterministik olmayan) ve doz-gridinin AMPIRIK DESTEGI olan
  stratumlarda kimliklendirebilir. Positivity/overlap ihlali -> ekstrapolasyon.

  Her stratum icin:
    - doz>0 ve doz=0 birlikte var mi (deterministik degil mi)
    - hasta/satir sayisi yeterli mi
    - doz-gridinin (0..12 U, 0.5 adim; 08 ile ayni) kac noktasi destekli
  -> 'curve_identifiable' bayragi. HICBIR stratum kimliklendirilemezse
     08 icin kritik STOP sinyali (main ciktisinda POSITIVITY VERDICT).

REGISTRY NOTU (EXCLUSION_REGISTRY 'C_suspect_undocumented'):
  sens_suspect_undoc=TRUE satirlar dose=0 hucresini KIRLETEBILIR (pompa/dis
  sistem, belgelenmemis doz). 0_no_dose hucresi hem ham hem 'temiz'
  (sens_suspect_undoc haric) raporlanir.

KURULUM: pip install duckdb pandas numpy matplotlib pyarrow
KULLANIM: python python/06_positivity_check.py --parquet 'data/*.parquet' --outdir qc/
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 08_dose_optimizer.py ile AYNI doz gridi (np.arange(0,12.5,0.5))
DOSE_GRID = np.arange(0.0, 12.5, 0.5)
HALF_BIN  = 0.25                      # grid noktasi +-0.25 U penceresi

# Esikler -- CLAUDE.md Kural 10 (BAGLAYICI; degisiklik kullanici onayi ister):
#   "hucrede >= 50 satir VE >= 20 hasta; egri icin >= 3 destekli doz seviyesi"
MIN_CELL_ROWS  = 50                   # bir doz-grid hucresi destekli: >= 50 satir
MIN_CELL_SUBJ  = 20                   #                             VE >= 20 hasta
MIN_SUP_LEVELS = 3                    # egri kimliklendirilebilir: >= 3 destekli seviye
FRAC_DOSED_LO, FRAC_DOSED_HI = 0.0, 1.0     # yalniz raporlama (Kural 10 kapi degil)


def load(parquet, split, route_filter=None, extra_filter=None):
    conds = ["TRUE"]
    if split:
        conds.append(f"split = '{split}'")
    if route_filter:                      # kapsam duyarliligi (or. EN,PN,EN+PN)
        routes = ",".join("'" + r.replace("'", "''") + "'" for r in route_filter)
        conds.append(f"nutrition_route IN ({routes})")
    if extra_filter:                      # serbest SQL WHERE parcasi (parquet kolonlari)
        conds.append(f"({extra_filter})")
    q = f"SELECT * FROM read_parquet('{parquet}') WHERE {' AND '.join(conds)}"
    return duckdb.connect().execute(q).df()


def strata_defs(df):
    """(stratum_adi, {deger_etiketi: bool_maske}) listesi. Sutunu OLMAYAN
    stratumlar atlanir (model_set vs olay-duzeyi basal_events uyumu)."""
    C = df.columns
    out = [("overall", {"all": pd.Series(True, index=df.index)})]
    if "cirrhosis" in C:
        out.append(("cirrhosis", {
            "cirrhosis=1": df.cirrhosis == True, "cirrhosis=0": df.cirrhosis == False}))
    if "egfr_now" in C:
        egfr = df.egfr_now
        out.append(("egfr_band", {
            "egfr<15":    egfr < 15,
            "egfr_15_45": (egfr >= 15) & (egfr < 45),
            "egfr_45_90": (egfr >= 45) & (egfr < 90),
            "egfr>=90":   egfr >= 90}))
    if "on_crrt" in C:
        out.append(("on_crrt", {"crrt=1": df.on_crrt == 1, "crrt=0": df.on_crrt == 0}))
    if "on_steroid" in C:
        out.append(("on_steroid", {
            "steroid=1": df.on_steroid == 1, "steroid=0": df.on_steroid == 0}))
    if "dm_t1" in C:
        out.append(("dm_t1", {"dm_t1=1": df.dm_t1 == True, "dm_t1=0": df.dm_t1 == False}))
    if "basal_type" in C:                 # olay-duzeyi basal veri seti
        out.append(("basal_type", {f"bt={v}": df.basal_type == v
                                    for v in sorted(df.basal_type.dropna().unique())}))
    if "nutrition_route" in C:
        nr = df.nutrition_route.fillna("NA")
        out.append(("nutrition_route", {f"nr={v}": nr == v for v in sorted(nr.unique())}))
    if "cirrhosis" in C and "cho_interrupted_4h" in C:
        out.append(("int01_cirr_x_interrupt", {
            "cirr&interrupt": (df.cirrhosis == True) & (df.cho_interrupted_4h == 1),
            "cirr&cont":      (df.cirrhosis == True) & (df.cho_interrupted_4h == 0)}))
    return out


def grid_support(dose_vals, subj_vals):
    """Her grid noktasi: satir sayisi, HASTA sayisi, destekli mi.
    Kural 10: destekli = satir >= MIN_CELL_ROWS VE hasta >= MIN_CELL_SUBJ."""
    counts = np.zeros(len(DOSE_GRID), dtype=int)
    subj_counts = np.zeros(len(DOSE_GRID), dtype=int)
    for i, g in enumerate(DOSE_GRID):
        if i == 0:                       # 0 hucresi: [0, 0.25)
            mask = (dose_vals >= 0) & (dose_vals < HALF_BIN)
        else:
            mask = np.abs(dose_vals - g) < HALF_BIN
        counts[i] = int(mask.sum())
        subj_counts[i] = int(np.unique(subj_vals[mask]).size)
    supported = (counts >= MIN_CELL_ROWS) & (subj_counts >= MIN_CELL_SUBJ)
    return counts, subj_counts, supported


def eval_cell(sub):
    """Bir stratum-hucresi icin positivity metrikleri (Kural 10)."""
    n_rows = len(sub)
    n_subj = sub.subject_id.nunique()
    dose = sub.dose.values
    subj = sub.subject_id.values
    pos = dose[dose > 0]
    n_dosed = int((dose > 0).sum())
    frac_dosed = n_dosed / n_rows if n_rows else 0.0
    counts, subj_counts, supported = grid_support(dose, subj)
    n_sup_levels = int(supported.sum())                 # 0 dahil destekli seviye
    n_pos_supported = int(supported[1:].sum())          # 0 haric (bilgi)
    frac_above_grid = float((dose > DOSE_GRID[-1]).mean()) if n_rows else 0.0
    # Kural 10: egri kimliklendirilebilir <=> >= 3 destekli doz seviyesi
    identifiable = bool(n_sup_levels >= MIN_SUP_LEVELS)
    pq = (np.percentile(pos, [1, 25, 50, 75, 99]).round(2).tolist()
          if len(pos) else [None] * 5)
    return dict(
        n_rows=n_rows, n_subj=n_subj, n_dosed=n_dosed,
        frac_dosed=round(frac_dosed, 4),
        pos_dose_p01_25_50_75_99=pq,
        n_supported_levels=n_sup_levels,
        n_pos_grid_supported=n_pos_supported,
        frac_dose_above_grid=round(frac_above_grid, 4),
        in_support=supported.tolist(),
        grid_counts=counts.tolist(),
        grid_subj_counts=subj_counts.tolist(),
        curve_identifiable=identifiable)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--outdir", default="qc")
    ap.add_argument("--split", default=None,
                    help="belirtilirse yalniz o split; yoksa TUM satirlar")
    ap.add_argument("--route-filter", default=None,
                    help="virgullu nutrition_route listesi; or. 'EN,PN,EN+PN' "
                         "(kapsam duyarliligi). Yoksa TUM rotalar.")
    ap.add_argument("--extra-filter", default=None,
                    help="serbest SQL WHERE parcasi (parquet kolonlari); or. "
                         "\"nutrition_route IN ('EN','PN','EN+PN') OR npo_active = 1\"")
    ap.add_argument("--dose-col", default="dose_u_this_hr",
                    help="doz sutunu (or. model_set: dose_u_this_hr; basal: dose_u_kg)")
    ap.add_argument("--dose-grid", default="12,0.5",
                    help="'max,step' -- doz gridi (or. model_set '12,0.5'; U/kg '1.0,0.05')")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    route_filter = [r.strip() for r in a.route_filter.split(",")] if a.route_filter else None

    global DOSE_GRID, HALF_BIN
    dmax, dstep = [float(x) for x in a.dose_grid.split(",")]
    DOSE_GRID = np.round(np.arange(0.0, dmax + dstep / 2, dstep), 6)
    HALF_BIN = dstep / 2

    df = load(a.parquet, a.split, route_filter, a.extra_filter)
    df["dose"] = df[a.dose_col]
    df = df[df["dose"].notna()].reset_index(drop=True)
    print(f"[06] yuklendi: {len(df):,} satir, {df.subject_id.nunique():,} hasta"
          f"  (split={a.split or 'ALL'}, route_filter={route_filter or 'ALL'}, "
          f"extra_filter={a.extra_filter or 'None'}, dose_col={a.dose_col}, "
          f"grid=0..{dmax}/{dstep})")

    # ---- sens_suspect_undoc: 0_no_dose hucresi kirliligi (varsa) ----------
    zero = df[df.dose == 0]
    n_zero = len(zero)
    n_zero_suspect = (int((zero.sens_suspect_undoc == True).sum())
                      if "sens_suspect_undoc" in df.columns else 0)
    contam = dict(
        n_zero_dose=n_zero,
        n_zero_suspect_undoc=n_zero_suspect,
        frac_zero_contaminated=round(n_zero_suspect / n_zero, 4) if n_zero else 0.0)
    print(f"[06] 0_no_dose hucresi: {n_zero:,} satir, "
          f"{n_zero_suspect:,} ({contam['frac_zero_contaminated']:.1%}) sens_suspect_undoc")

    # ---- stratum bazli positivity ----------------------------------------
    results, summary_rows, identifiable_strata = {}, [], []
    for sname, cells in strata_defs(df):
        results[sname] = {}
        for clabel, mask in cells.items():
            sub = df[mask.values if hasattr(mask, "values") else mask]
            if len(sub) == 0:
                continue
            m = eval_cell(sub)
            results[sname][clabel] = m
            summary_rows.append(dict(stratum=sname, cell=clabel, **{
                k: m[k] for k in ("n_rows", "n_subj", "n_dosed", "frac_dosed",
                                   "n_supported_levels", "n_pos_grid_supported",
                                   "frac_dose_above_grid", "curve_identifiable")}))
            if m["curve_identifiable"] and sname != "overall":
                identifiable_strata.append(f"{sname}:{clabel}")

    overall_ok = results["overall"]["all"]["curve_identifiable"]
    any_curve = overall_ok or len(identifiable_strata) > 0

    # ---- ciktilar ---------------------------------------------------------
    support_json = dict(
        dose_grid=DOSE_GRID.tolist(),
        params=dict(MIN_CELL_ROWS=MIN_CELL_ROWS, MIN_CELL_SUBJ=MIN_CELL_SUBJ,
                    MIN_SUP_LEVELS=MIN_SUP_LEVELS, HALF_BIN=HALF_BIN,
                    dose_col=a.dose_col, dose_grid=a.dose_grid,
                    rule="CLAUDE.md Kural 10"),
        split=a.split or "ALL",
        route_filter=route_filter or "ALL",
        extra_filter=a.extra_filter or "None",
        n_rows=int(len(df)), n_subjects=int(df.subject_id.nunique()),
        zero_dose_contamination=contam,
        overall_in_support=results["overall"]["all"]["in_support"],
        overall_curve_identifiable=overall_ok,
        curve_identifiable_strata=identifiable_strata,
        any_curve=any_curve,
        strata=results)
    (out / "positivity_support.json").write_text(json.dumps(support_json, indent=2))

    summ = pd.DataFrame(summary_rows)
    summ.to_csv(out / "positivity_summary.csv", index=False)

    _write_model_card(out, support_json, summ)
    _plot(out, results)

    # ---- VERDICT (STOP sinyali) ------------------------------------------
    print("\n=== POSITIVITY VERDICT ===")
    print(f"  overall curve_identifiable : {overall_ok}")
    print(f"  identifiable strata (n={len(identifiable_strata)}): "
          f"{identifiable_strata[:8]}{' ...' if len(identifiable_strata)>8 else ''}")
    if any_curve:
        print("  VERDICT: OK -- en az bir stratumda dozyanit egrisi kurulabilir.")
    else:
        print("  VERDICT: STOP -- HICBIR STRATUMDA EGRI KURULAMIYOR (positivity ihlali).")
    print(f"\nCiktilar -> {out.resolve()}")
    return 0 if any_curve else 3


def _write_model_card(out, sj, summ):
    L = []
    L.append("# Model Card — Positivity / Common-Support Destegi\n")
    L.append("_06_positivity_check.py ciktisi. Script gece-modu oturumunda YAZILDI; "
             "esikler CLAUDE.md Kural 10'a hizali. Stratum tanimlari S. Efe onayi bekliyor._\n")
    L.append(f"- Veri: {sj['n_rows']:,} satir, {sj['n_subjects']:,} hasta "
             f"(split={sj['split']}, route_filter={sj.get('route_filter','ALL')})")
    p = sj["params"]
    L.append(f"- Doz sutunu: `{p.get('dose_col','dose_u_this_hr')}`; grid 0..{DOSE_GRID[-1]:.2f} "
             f"(max,step={p.get('dose_grid','12,0.5')}, {len(DOSE_GRID)} nokta, "
             f"pencere +-{p['HALF_BIN']})")
    L.append(f"- Esikler (**{p['rule']}**): bir doz-seviyesi hucresi destekli <=> "
             f">= {p['MIN_CELL_ROWS']} satir VE >= {p['MIN_CELL_SUBJ']} hasta; "
             f"egri kimliklendirilebilir <=> >= {p['MIN_SUP_LEVELS']} destekli seviye.\n")
    c = sj["zero_dose_contamination"]
    L.append("## 0_no_dose hucresi kirliligi (registry C_suspect_undocumented)\n")
    L.append(f"- 0-doz satir: {c['n_zero_dose']:,}; "
             f"sens_suspect_undoc: {c['n_zero_suspect_undoc']:,} "
             f"({c['frac_zero_contaminated']:.1%}). Duyarlilik analizinde cikarilir.\n")
    L.append("## Genel karar\n")
    L.append(f"- **overall curve_identifiable**: {sj['overall_curve_identifiable']}")
    L.append(f"- **any_curve** (STOP degil): {sj['any_curve']}")
    L.append(f"- kimliklendirilebilir stratum sayisi: "
             f"{len(sj['curve_identifiable_strata'])}\n")
    L.append("## Stratum ozeti\n")
    L.append("| stratum | hucre | n_satir | n_hasta | n_dozlu | doz_orani | "
             "destekli_seviye | (poz) | egri? |")
    L.append("|---|---|---:|---:|---:|---:|---:|---:|:--:|")
    for _, r in summ.iterrows():
        L.append(f"| {r.stratum} | {r.cell} | {r.n_rows:,} | {r.n_subj:,} | "
                 f"{r.n_dosed:,} | {r.frac_dosed:.2f} | {r.n_supported_levels} | "
                 f"{r.n_pos_grid_supported} | "
                 f"{'EVET' if r.curve_identifiable else 'hayir'} |")
    L.append("\n_Not: 'egri?' sutunu her stratum-hucresinde dozyanit egiminin "
             "AMPIRIK olarak kimliklendirilebilirligidir; nedensel gecerlilik "
             "degil, yalniz overlap/positivity kosuludur._\n")
    (out / "model_card_support.md").write_text("\n".join(L), encoding="utf-8")


def _plot(out, results):
    keys = [("overall", "all"), ("cirrhosis", "cirrhosis=1"),
            ("egfr_band", "egfr<15"), ("on_steroid", "steroid=1")]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, (s, c) in zip(axes.ravel(), keys):
        m = results.get(s, {}).get(c)
        if not m:
            ax.set_visible(False); continue
        ax.bar(DOSE_GRID, m["grid_counts"], width=0.4,
               color=["#c44" if not ok else "#39c" for ok in m["in_support"]])
        ax.axhline(MIN_CELL_ROWS, color="k", ls="--", lw=0.7)
        ax.set_title(f"{s}:{c}  (egri={'E' if m['curve_identifiable'] else 'H'})")
        ax.set_xlabel("doz (U)"); ax.set_ylabel("satir")
    fig.suptitle("Positivity: doz-grid destegi (kirmizi=destegsiz hucre)")
    fig.tight_layout()
    fig.savefig(out / "positivity_support.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

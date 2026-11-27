import duckdb, pandas as pd, numpy as np, json
from scipy import stats
df = duckdb.connect().execute(
    "SELECT * FROM read_parquet('data/basal_events/*.parquet') "
    "WHERE basal_type='glargine' AND units_total IS NOT NULL").df()
h = pd.read_csv("/tmp/ibw.csv")
df = df.merge(h, on="stay_id", how="left")
df["eb"] = pd.cut(df.egfr_now, [-1, 15, 45, 90, 1e9],
                  labels=["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]).astype(str)
inch = df.height_cm / 2.54
df["ibw"] = np.where(df.gender == "M", 50 + 2.3 * (inch - 60), 45.5 + 2.3 * (inch - 60))
df.loc[df.ibw < 30, "ibw"] = np.nan
df["u_ibw"] = df.units_total / df.ibw
df["wt_gain"] = df.weight_kg - df.admit_wt
BANDS = ["egfr<15", "egfr_15_45", "egfr_45_90", "egfr>=90"]


def med(x):
    x = pd.Series(x).dropna()
    return round(float(np.median(x)), 3) if len(x) else None


out = {"height_avail_frac": round(float(df.height_cm.notna().mean()), 3),
       "note": "IBW=Devine (50/45.5 + 2.3*(inch-60)). Odemli agirlik U/kg'yi bastirir; "
               "IBW gercek metabolik kutleye dozu gosterir.",
       "by_band": {}}
for b in BANDS:
    d = df[df.eb == b]
    out["by_band"][b] = dict(n=int(len(d)), u_kg_actual_med=med(d.dose_u_kg),
                             u_kg_IBW_med=med(d.u_ibw),
                             admit_to_dosetime_wt_gain_med=med(d.wt_gain))
lo, hi = df[df.eb == "egfr<15"], df[df.eb == "egfr>=90"]
out["gap_actual_ukg"] = round(med(lo.dose_u_kg) / med(hi.dose_u_kg), 2) if med(hi.dose_u_kg) else None
out["gap_IBW_ukg"] = round(med(lo.u_ibw) / med(hi.u_ibw), 2) if med(hi.u_ibw) else None
patm = df.dropna(subset=["u_ibw"]).groupby("subject_id").agg(u=("u_ibw", "median"), b=("eb", "first"))
grps = [patm[patm.b == b].u.dropna() for b in BANDS if (patm.b == b).any()]
out["kw_ukg_IBW_p"] = round(float(stats.kruskal(*grps)[1]), 4)
json.dump(out, open("qc/paper2_tables/T3e_ibw.json", "w"), indent=2)
print("height_avail", out["height_avail_frac"])
for b in BANDS:
    r = out["by_band"][b]
    print("  %-11s U/kg_actual=%s U/kg_IBW=%s wt_gain=%s" %
          (b, r["u_kg_actual_med"], r["u_kg_IBW_med"], r["admit_to_dosetime_wt_gain_med"]))
print("gap actual U/kg:", out["gap_actual_ukg"], "-> IBW U/kg:", out["gap_IBW_ukg"],
      "| KW_IBW p:", out["kw_ukg_IBW_p"])

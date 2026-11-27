import duckdb, pandas as pd, numpy as np, json
from scipy import stats
df = duckdb.connect().execute(
    "SELECT * FROM read_parquet('data/basal_events/*.parquet') "
    "WHERE basal_type='glargine' AND dose_u_kg IS NOT NULL").df()
fl = pd.read_csv("/tmp/dmflags.csv")
df = df.merge(fl, on="hadm_id", how="left")
df["oral"] = df.oral.fillna(0)


def ph(r):
    if r.oral == 1:
        return "T2DM_oral"
    if pd.notna(r.hba1c) and r.hba1c >= 6.5:
        return "kronikDM_HbA1c"
    return "stres_veya_olcumsuz"


df["alt_dm"] = df.apply(ph, axis=1)


def mi(x):
    x = pd.Series(x).dropna()
    if len(x) == 0:
        return None
    return [round(float(np.median(x)), 3), round(float(np.percentile(x, 25)), 3),
            round(float(np.percentile(x, 75)), 3)]


cats = ["T2DM_oral", "kronikDM_HbA1c", "stres_veya_olcumsuz"]
tab = {k: dict(n_events=int((df.alt_dm == k).sum()),
               n_subj=int(df[df.alt_dm == k].subject_id.nunique()),
               u_kg_med_iqr=mi(df[df.alt_dm == k].dose_u_kg),
               hypo_night_l1=round(float(df[df.alt_dm == k].hypo_night_l1.mean()), 4)
               if (df.alt_dm == k).any() else None) for k in cats}
pat = df.groupby("subject_id").agg(u=("dose_u_kg", "median"), s=("alt_dm", "first"))
grps = [pat[pat.s == k].u.dropna().values for k in cats]
print("grp sizes:", {k: len(g) for k, g in zip(cats, grps)})
grps = [g for g in grps if len(g) > 1]
H, p = stats.kruskal(*grps)
out = dict(
    proxy_availability=dict(hba1c_n=int(df.hba1c.notna().sum()),
                            hba1c_med=round(float(np.nanmedian(df.hba1c)), 1),
                            oral_n_events=int((df.oral == 1).sum())),
    note="OLCULU fenotip: oral antidiyabetik=T2DM; HbA1c>=6.5 (oralsiz)=kronikDM; "
         "else stres/olcumsuz. home-insulin MIMIC-IV'te temiz degil -> HbA1c+oral. "
         "stres kovasi olcumsuz DM'yi de icerir (SINIRLILIK).",
    by_alt_dm=tab, kw_u_kg_p=round(float(p), 4),
    verdict=("DM ayagi SAGLAM: alternatif OLCULU fenotiple de anlamli (p<0.05)"
             if p < 0.05 else
             "DM ayagi ZAYIF: alternatif fenotiple anlamli DEGIL -> T3c yalniz agirlik"))
json.dump(out, open("qc/paper2_tables/T12c_dm_proxy.json", "w"), indent=2)
print("alt_dm u_kg med:", {k: (tab[k]["u_kg_med_iqr"] or [None])[0] for k in cats})
print("n_subj:", {k: tab[k]["n_subj"] for k in cats})
print("KW p =", round(p, 4), "->", out["verdict"])

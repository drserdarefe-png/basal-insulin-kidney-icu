#!/usr/bin/env python3
"""
07 — IKI ASAMALI ISF MODELI (Bayesci hiyerarsik, numpyro)     v2.0
==================================================================
Uygulanmis: PATCH-02 (z_cp_high SIRROZ KAPISI -- kritik)

ASAMA 1: hasta-spesifik insulin duyarlilik faktoru (ISF)
    d_glu_t2h = b0_i - ISF_i * doz + gamma*X + eps
ASAMA 2 ayni modelde: theta = komorbiditenin ISF uzerindeki etkisi
    "eGFR<15 ISF'i %X artirir (%95 CrI ...)" -> yorumlanabilir cikti

KURAL 5 (bozulamaz): ISF_i = exp(mu + Z_i@theta + u_i) > 0
  -> doz artinca glukoz DUSER, monotonluk MIMARIDEN garanti.
  exp() KALDIRILAMAZ, lineer ISF ile DEGISTIRILEMEZ.

NEDEN GBM DEGIL: hasta-spesifik rastgele egim yok, monotonluk garantisi yok,
  belirsizlik kalibrasyonu yok (hipoglisemi kisiti icin SART).

KURULUM: pip install numpyro jax duckdb pandas arviz matplotlib
KULLANIM: python python/07_isf_two_stage.py --parquet 'data/*.parquet' --split train
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import jax, jax.numpy as jnp
import numpyro, numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS

numpyro.set_host_device_count(4)

Z_COLS = ["z_egfr_low", "z_egfr_mid", "z_crrt", "z_cirrhosis", "z_cp_high",
          "z_bun", "z_bmi", "z_dm_t1", "z_steroid"]
X_COLS_BASE = ["c_cho", "c_dcho", "c_iob", "c_basal", "c_glu", "c_vaso", "c_poc"]
# YAMA: glukoz-egilimi (c_dglu1/3) + rota kuklalari (is_none/is_maint).
#   is_none/is_maint o rotalardaki ORTALAMA olculemeyen CHO'yu sogurur -> b0_i
#   sismez, ISF korunur. c_dglu1 hem beslenme hem onceki insulini yansitir ->
#   c_iob ile kolinearite riski (bu yuzden --old-xcols duyarlilik kolu var).
X_COLS_FULL = X_COLS_BASE + ["c_dglu1", "c_dglu3", "is_none", "is_maint"]
X_COLS = X_COLS_FULL   # varsayilan; main() --old-xcols ile X_COLS_BASE secebilir


def load(parquet, split, extra_filter=None):
    where = (f"WHERE split = '{split}' AND d_glu_t2h IS NOT NULL "
             "AND glu_now IS NOT NULL AND dose_u_this_hr IS NOT NULL")
    if extra_filter:                      # kapsam duyarliligi (or. clean = fed OR npo)
        where += f" AND ({extra_filter})"
    return duckdb.connect().execute(f"""
        SELECT subject_id, stay_id, hr, d_glu_t2h, glu_now,
               dose_u_this_hr, iob_bolus_u, iob_basal_u,
               cho_g_hr, d_cho_2h, cho_interrupted_4h,
               d_glu_1h, d_glu_3h, nutrition_route,
               egfr_now, on_crrt, bun, cirrhosis, child_pugh_score,
               dm_t1, pred_eq_24h, weight_kg, on_vasopressor, glu_src
        FROM read_parquet('{parquet}')
        {where}
    """).df()


def prep(df):
    df = df.copy()
    g = df.groupby("subject_id").agg(
        egfr=("egfr_now", "median"), crrt=("on_crrt", "max"),
        bun=("bun", "median"),       cirr=("cirrhosis", "max"),
        cp=("child_pugh_score", "median"), dm_t1=("dm_t1", "max"),
        ster=("pred_eq_24h", "median"), wt=("weight_kg", "median"),
    ).reset_index()

    g["z_egfr_low"]  = (g.egfr < 15).astype(float)
    g["z_egfr_mid"]  = ((g.egfr >= 15) & (g.egfr < 45)).astype(float)
    g["z_crrt"]      = g.crrt.fillna(0).astype(float)
    g["z_cirrhosis"] = g.cirr.fillna(0).astype(float)

    # ---- PATCH-02: SIRROZ KAPISI. KALDIRILAMAZ. -------------------------
    # Child-Pugh = bilirubin + INR + albumin + asit + HE.
    # Non-sirotik SEPTIK hastada ucu de bozulur -> cp>=10 -> z_cp_high=1
    # -> Asama 2'de "ileri KC hastaligi ISF'i degistiriyor" katsayisi
    #    aslinda SEPSIS SIDDETINI olcer. theta_cirrhosis ile karisir,
    #    ikisi birlikte yorumlanamaz hale gelir.
    # 04 (v2.0) zaten cp'yi non-sirotikte NULL birakiyor; bu IKINCI katman.
    g["z_cp_high"] = ((g.cirr == 1) & (g.cp >= 10)).fillna(False).astype(float)

    g["z_bun"]     = ((g.bun - 25) / 20).fillna(0).clip(-2, 4)     # merkezlenmis
    g["z_bmi"]     = ((g.wt - 80) / 20).fillna(0).clip(-2, 4)      # agirlik proxy
    g["z_dm_t1"]   = g.dm_t1.fillna(0).astype(float)
    g["z_steroid"] = (g.ster.fillna(0) / 40).clip(0, 3)

    # ---- PATCH-02 DOGRULAMASI: kesisim BOS olmali --------------------
    assert (g.loc[g.cirr == 0, "z_cp_high"] == 0).all(), \
        "PATCH-02 IHLALI: non-sirotikte z_cp_high=1. Sirroz kapisi calismyor."
    print("[PATCH-02] z_cp_high prevalansi (cirr grubuna gore):")
    print(g.groupby("cirr")["z_cp_high"].mean().to_string())

    df = df.merge(g[["subject_id"] + Z_COLS], on="subject_id", how="left")

    df["c_cho"]   = df.cho_g_hr.fillna(0) / 10.0        # 10 g/h birimi
    df["c_dcho"]  = df.d_cho_2h.fillna(0) / 10.0
    df["c_iob"]   = df.iob_bolus_u.fillna(0)
    df["c_basal"] = df.iob_basal_u.fillna(0)
    df["c_glu"]   = (df.glu_now - 150) / 50.0           # regresyona-ortalama
    df["c_vaso"]  = df.on_vasopressor.fillna(0).astype(float)
    df["c_poc"]   = (df.glu_src == "POC").astype(float)
    # YAMA: glukoz-egilimi + rota kuklalari (X_COLS_FULL)
    df["c_dglu1"]  = df.d_glu_1h.fillna(0) / 10.0
    df["c_dglu3"]  = df.d_glu_3h.fillna(0) / 10.0
    df["is_none"]  = (df.nutrition_route == 'none').astype(float)
    df["is_maint"] = (df.nutrition_route == 'MAINT_ONLY').astype(float)
    df["dose"]    = df.dose_u_this_hr.fillna(0)
    df["y"]       = df.d_glu_t2h / 10.0                 # 10 mg/dL birimi

    df = df.dropna(subset=["y", "dose"] + Z_COLS)
    codes, uniq = pd.factorize(df.subject_id)
    df["pid"] = codes
    return df, uniq


def model(pid, dose, Z, X, n_pat, y=None):
    theta = numpyro.sample("theta", dist.Normal(0, 0.5).expand([Z.shape[1]]))
    mu_b  = numpyro.sample("mu_log_isf", dist.Normal(np.log(1.5), 0.5))
    sig_b = numpyro.sample("sigma_log_isf", dist.HalfNormal(0.5))
    with numpyro.plate("pat", n_pat):
        u  = numpyro.sample("u",  dist.Normal(0, 1))     # non-centered
        b0 = numpyro.sample("b0", dist.Normal(0, 1.0))
    # KURAL 5: exp() -> ISF > 0 -> doz artinca glukoz DUSER (monotonluk)
    isf = numpyro.deterministic("isf", jnp.exp(mu_b + Z @ theta + u * sig_b))
    gamma = numpyro.sample("gamma", dist.Normal(0, 1).expand([X.shape[1]]))
    mu = b0[pid] - isf[pid] * dose + X @ gamma
    sigma = numpyro.sample("sigma", dist.HalfNormal(3.0))
    nu = numpyro.sample("nu", dist.Gamma(2, 0.1))        # Student-t: outlier'a dayanikli
    numpyro.sample("obs", dist.StudentT(nu, mu, sigma), obs=y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--split", default="train")
    ap.add_argument("--outdir", default="qc/isf_out")
    ap.add_argument("--warmup", type=int, default=1000)
    ap.add_argument("--samples", type=int, default=1000)
    ap.add_argument("--max-rows", type=int, default=400_000)
    ap.add_argument("--old-xcols", action="store_true",
                    help="eski 7-terim X_COLS (glukoz-egilimi/rota kuklalari YOK) "
                         "-- kolinearite duyarlilik kolu")
    ap.add_argument("--drop-xcols", default=None,
                    help="X_COLS'tan cikarilacak virgullu liste (or. clean scope'ta "
                         "'is_none,is_maint' -> rank-eksik/huni onlenir)")
    ap.add_argument("--extra-filter", default=None,
                    help="serbest SQL WHERE parcasi (parquet kolonlari); or. "
                         "\"nutrition_route IN ('EN','PN','EN+PN') OR npo_active = 1\"")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    xcols = X_COLS_BASE if a.old_xcols else X_COLS_FULL
    if a.drop_xcols:                      # rank-eksik kuklalari cikar (clean scope huni fix)
        drop = {x.strip() for x in a.drop_xcols.split(",")}
        xcols = [c for c in xcols if c not in drop]

    df, uniq = prep(load(a.parquet, a.split, a.extra_filter))
    if len(df) > a.max_rows:
        # HASTA bazli orneklе -- satir bazli orneklem hiyerarsiyi BOZAR
        keep = np.random.default_rng(42).choice(
            df.pid.unique(), size=int(df.pid.nunique()*a.max_rows/len(df)), replace=False)
        df = df[df.pid.isin(keep)]
        codes, uniq = pd.factorize(df.subject_id); df["pid"] = codes
    print(f"\nn satir={len(df):,}  n hasta={df.pid.nunique():,}  X_COLS={xcols}")

    Zp = df.groupby("pid")[Z_COLS].first().sort_index().values
    X  = df[xcols].values

    mcmc = MCMC(NUTS(model, target_accept_prob=0.9),
                num_warmup=a.warmup, num_samples=a.samples, num_chains=4)
    mcmc.run(jax.random.PRNGKey(0),
             pid=jnp.array(df.pid.values), dose=jnp.array(df.dose.values),
             Z=jnp.array(Zp), X=jnp.array(X), n_pat=df.pid.nunique(),
             y=jnp.array(df.y.values))
    mcmc.print_summary(exclude_deterministic=True)
    s = mcmc.get_samples()

    rows = []
    for j, c in enumerate(Z_COLS):
        mult = np.exp(s["theta"][:, j])          # log-ISF -> carpan
        rows.append(dict(variable=c,
                         isf_multiplier=float(np.median(mult)),
                         cri_lo=float(np.percentile(mult, 2.5)),
                         cri_hi=float(np.percentile(mult, 97.5)),
                         p_gt_1=float((mult > 1).mean())))
    res = pd.DataFrame(rows); res.to_csv(out/"stage2_isf_effects.csv", index=False)

    print("\n=== ASAMA 2: komorbidite -> ISF carpani (MAKALENIN ANA TABLOSU) ===")
    print("(>1 = insulin daha ETKILI = doz azaltilmali)")
    for _, r in res.iterrows():
        print(f"  {r.variable:<14} x{r.isf_multiplier:5.2f} "
              f"[{r.cri_lo:.2f}-{r.cri_hi:.2f}]  P(>1)={r.p_gt_1:.3f}")

    base = float(np.median(np.exp(s["mu_log_isf"]))) * 10
    print(f"\nPopulasyon medyan ISF: {base:.1f} mg/dL dusus / unite")
    print(f"Hasta arasi SD (log)  : {float(np.median(s['sigma_log_isf'])):.3f}")

    # (n_ornek, n_hasta) tutarli sekil: mu[:,None]=(S,1), theta@Zp.T=(S,P),
    # u*sigma[:,None]=(S,P). (BUG DUZ: eski Zp@theta.T=(P,S) yayinlanamiyordu.)
    isf = np.exp(s["mu_log_isf"][:, None] + s["theta"] @ Zp.T
                 + s["u"] * s["sigma_log_isf"][:, None]) * 10
    pd.DataFrame({"subject_id": uniq,
                  "isf_median": np.median(isf, axis=0),
                  "isf_lo": np.percentile(isf, 2.5, axis=0),
                  "isf_hi": np.percentile(isf, 97.5, axis=0),
                  }).to_csv(out/"patient_isf.csv", index=False)
    # NOT: patient_isf.csv HASTA DUZEYI -> .gitignore kapsaminda, commit YOK

    (out/"fixed_effects.json").write_text(json.dumps({
        c: dict(median=float(np.median(s["gamma"][:, j])),
                lo=float(np.percentile(s["gamma"][:, j], 2.5)),
                hi=float(np.percentile(s["gamma"][:, j], 97.5)))
        for j, c in enumerate(xcols)}, indent=2))
    print(f"\nCiktilar -> {out.resolve()}")


if __name__ == "__main__":
    main()

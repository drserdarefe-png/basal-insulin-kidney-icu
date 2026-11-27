#!/usr/bin/env python3
"""
20 — UC DUNYA simulasyonu (Makale 2 metodoloji: kimliklendirilebilirlik bariyerleri).

Gercek nedensel etki beta = -20 mg/dL per U/kg (ayni uc dunyada).
Dozun VARYANSI uc dunyada AYNI (SD~0.10 U/kg); degisen tek sey varyansin KAYNAGI:
  A: gurultu            (geri besleme YOK,  U YOK)   -> naif regresyon kurtarir
  B: gurultu+geribesleme(geri besleme VAR,  U YOK)   -> OLCULU zaman-degisen karistirma
                                                        -> MSM/IPTW kurtarmali (BARIYER 1)
  C: U + gurultu        (geri besleme VAR,  U GIZLI) -> OLCULMEZ karistirma
                                                        -> hicbir sey kurtarmaz (BARIYER 2)

DUZELTMELER (kullanici):
  1) POZITIVITE: AR katsayisi phi=0.30 (eski ~0.8 -> R2=0.96 -> IPTW patladi).
     Doz hata terimi SD=0.04. MSM'de ESS raporla; ESS<%30 -> "uygulanamaz", UYDURMA.
  2) VARYANS TUTARLILIGI: her dunyada gain oyle ki SD(doz)~=0.10 (0.45*U ACILIMI YOK).
     Dunya C'de c (U katsayisi) GRID SEARCH:
       hedef R2_olculmus=0.08+-0.03, R2_gecmisli=0.67+-0.05,
       kisit SD(doz)=0.10+-0.02, doz araligi 0.02-0.60.
  3) 500 tekrar, n=5000 ve n=1117; iki panel. NE CIKARSA RAPORLA.

KULLANIM:
  .venv/bin/python python/20_three_worlds.py --calibrate-only
  .venv/bin/python python/20_three_worlds.py --reps 500
"""
import argparse, json
from pathlib import Path
import numpy as np

# ---- SABIT DGP parametreleri ----
BETA = -20.0          # gercek nedensel etki (mg/dL per U/kg)
M = 0.20              # ortalama doz (U/kg)
PHI = 0.30            # AR/geri-besleme katsayisi (DUZELTME 1)
ERR_SD = 0.04         # doz hata terimi SD (DUZELTME 1)
TARGET_SD = 0.10      # hedef SD(doz) (DUZELTME 2)
DOSE_LO, DOSE_HI = 0.02, 0.60
T = 6                 # hasta basina doz olayi (~6*1117 ~ 6700 satir, gercek 7560'a yakin)

# OLCULU zaman-degisen durum L (glukoz-benzeri, standart) + sonuc katsayilari
# TASARIM: dunya C'de dozu GIZLI U (persistan) surukler; L yalniz ZAYIF proxy.
#          analist prev_dose'u DEGIL, olculu L'yi ayarlar (R2_meas=measured-only).
RHO_L = 0.5           # L kalicilik
A1 = -2.0             # onceki doz -> L (tedavi-karistirici geri besleme)
KUL = 0.36            # U -> L (ZAYIF proxy; dunya C -> R2_meas ~ hedef 0.08)
SL = 1.0              # L gurultu SD
KFB_B = 0.022         # L -> doz geri besleme, dunya B (ORTA: pozitivite korunur)
KFB_C = 0.005         # L -> doz geri besleme, dunya C (tiny; ana surucu U)
KYL = 5.0             # L -> sonuc, dunya B (OLCULU karistirma)
KYU = 3.5             # U -> sonuc, dunya C (OLCULMEZ karistirma)
SY = 12.0             # sonuc gurultu SD


def simulate(world, c, gain, n, rng):
    """Uzunlamasina panel. Doner: A, Y, L(olculu), Aprev, U(gizli) (satir duzeyi, T*n)."""
    U = rng.normal(size=n)
    L_prev = rng.normal(0.0, 1.0, size=n)
    A_prev = np.full(n, M)
    A_all, Y_all, L_all, Ap_all, U_all = [], [], [], [], []
    for t in range(T):
        eL = rng.normal(0.0, SL, size=n)
        L = RHO_L * L_prev + A1 * (A_prev - M) + (KUL * U if world == "C" else 0.0) + eL
        if world == "A":
            syst = rng.normal(size=n)                      # KAYNAK: saf gurultu
        elif world == "B":
            syst = PHI * (A_prev - M) + KFB_B * L           # KAYNAK: OLCULU geri besleme
        else:  # C
            syst = PHI * (A_prev - M) + c * U + KFB_C * L    # KAYNAK: GIZLI U + gurultu
        err = rng.normal(0.0, ERR_SD, size=n)
        A = np.clip(M + gain * syst + err, DOSE_LO, DOSE_HI)
        eY = rng.normal(0.0, SY, size=n)
        conf = (KYL * L if world == "B" else (KYU * U if world == "C" else 0.0))
        Y = 140.0 + BETA * A + conf + eY                    # beta GERCEK
        A_all.append(A); Y_all.append(Y); L_all.append(L)
        Ap_all.append(A_prev); U_all.append(U)
        A_prev = A; L_prev = L
    return (np.concatenate(A_all), np.concatenate(Y_all), np.concatenate(L_all),
            np.concatenate(Ap_all), np.concatenate(U_all))


def calibrate_gain(world, c, rng, n=40000, target=TARGET_SD, err_sd=ERR_SD):
    """gain'i SD(doz)=target olacak sekilde bisection ile bul (AR ozyineli -> sim gerek)."""
    if world == "A":
        # A: syst=N(0,1) bagimsiz, err bagimsiz -> kapali form
        return float(np.sqrt(max(target**2 - err_sd**2, 1e-9)))
    lo, hi = 1e-4, 5.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        A = simulate(world, c, mid, n, np.random.default_rng(rng.integers(1 << 31)))[0]
        if A.std() < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    sst = ((y - y.mean()) ** 2).sum()
    r2 = 1.0 - (resid ** 2).sum() / sst if sst > 0 else 0.0
    return beta, resid, r2


def gps_r2(A, L, Aprev):
    """R2_olculmus = A~L (olculu durum) ; R2_gecmisli = A~L+Aprev (prev_dose eklenince)."""
    one = np.ones_like(A)
    _, _, r2m = ols(np.column_stack([one, L]), A)
    _, _, r2h = ols(np.column_stack([one, L, Aprev]), A)
    return r2m, r2h


def normlogpdf(x, mu, sd):
    return -0.5 * ((x - mu) / sd) ** 2 - np.log(sd) - 0.5 * np.log(2 * np.pi)


def estimators(A, Y, L, Aprev):
    """naive_crude, naive_adj, msm(IPTW stabilize) + ESS_frac.
    Analist yalniz OLCULU durum L'yi ayarlar (prev_dose'u DEGIL -> R2_meas seti)."""
    one = np.ones_like(A)
    # naif ham: Y ~ A
    b_crude, _, _ = ols(np.column_stack([one, A]), Y)
    # naif duzeltilmis: Y ~ A + L (olculu durum)
    b_adj, _, _ = ols(np.column_stack([one, A, L]), Y)
    # MSM stabilize agirlik: pay f(A) marjinal ; payda f(A|L olculu karistirici)
    muN = np.full_like(A, A.mean()); sN = A.std()
    bd, rd, _ = ols(np.column_stack([one, L]), A)
    muD = np.column_stack([one, L]) @ bd; sD = rd.std()
    logw = normlogpdf(A, muN, sN) - normlogpdf(A, muD, sD)
    logw -= logw.max()
    w = np.exp(logw)
    ess_frac = (w.sum() ** 2) / (np.sum(w ** 2) * len(w))
    # WLS Y ~ A (agirlikli)
    Xw = np.column_stack([one, A])
    W = w / w.mean()
    XtW = Xw.T * W
    b_msm = np.linalg.solve(XtW @ Xw, XtW @ Y)
    return dict(naive_crude=float(b_crude[1]), naive_adj=float(b_adj[1]),
                msm=float(b_msm[1]), ess_frac=float(ess_frac))


def summarize(vals):
    a = np.array(vals, float)
    return dict(mean=round(float(a.mean()), 2), sd=round(float(a.std()), 2),
                bias=round(float(a.mean() - BETA), 2),
                lo=round(float(np.percentile(a, 2.5)), 2),
                hi=round(float(np.percentile(a, 97.5)), 2))


def grid_search_c(rng):
    """Dunya C icin c: R2_gecmisli~0.67 & R2_olculmus~0.08, SD~0.10 kisiti."""
    grid = np.linspace(0.3, 5.0, 24)
    rows = []
    for c in grid:
        g = calibrate_gain("C", c, rng)
        A, Y, G, Ap, U = simulate("C", c, g, 60000, np.random.default_rng(rng.integers(1 << 31)))
        r2m, r2h = gps_r2(A, G, Ap)
        rows.append(dict(c=round(float(c), 3), gain=round(g, 4), sd=round(float(A.std()), 4),
                         r2_measured=round(r2m, 3), r2_history=round(r2h, 3),
                         frac_at_bound=round(float(((A <= DOSE_LO + 1e-9) | (A >= DOSE_HI - 1e-9)).mean()), 4)))
    # hedefe en yakin: R2_history 0.67'ye, SD 0.10 kisiti icinde
    ok = [r for r in rows if abs(r["sd"] - 0.10) <= 0.02]
    pool = ok or rows
    best = min(pool, key=lambda r: abs(r["r2_history"] - 0.67))
    return best, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--calibrate-only", action="store_true")
    ap.add_argument("--seed", type=int, default=20260717)
    ap.add_argument("--outdir", default="qc")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.outdir); out.mkdir(exist_ok=True)

    # --- 1) Dunya C icin c grid search ---
    best, grid = grid_search_c(rng)
    C_STAR = best["c"]
    print("== GRID SEARCH (dunya C) ==")
    print(f"{'c':>6} {'gain':>7} {'SD':>6} {'R2_meas':>8} {'R2_hist':>8} {'bound%':>7}")
    for r in grid:
        mark = " <-- SECILDI" if r["c"] == C_STAR else ""
        print(f"{r['c']:6.3f} {r['gain']:7.4f} {r['sd']:6.3f} {r['r2_measured']:8.3f} "
              f"{r['r2_history']:8.3f} {r['frac_at_bound']*100:6.1f}%{mark}")
    print(f"\nSECILEN c*={C_STAR}  (R2_meas={best['r2_measured']}, R2_hist={best['r2_history']}, "
          f"SD={best['sd']}, gain={best['gain']})")

    # --- gain kalibrasyonu (A,B world C secilen c ile) ---
    gains = {"A": calibrate_gain("A", 0.0, rng),
             "B": calibrate_gain("B", 0.0, rng),
             "C": calibrate_gain("C", C_STAR, rng)}
    cs = {"A": 0.0, "B": 0.0, "C": C_STAR}
    print("\n== KALIBRE gain (SD~0.10) + dunya-basi R2 ==")
    for w in "ABC":
        A, Y, L, Ap, U = simulate(w, cs[w], gains[w], 60000, np.random.default_rng(rng.integers(1 << 31)))
        r2m, r2h = gps_r2(A, L, Ap)
        print(f"  dunya {w}: gain={gains[w]:.4f}  SD={A.std():.4f}  "
              f"R2_meas={r2m:.3f}  R2_hist={r2h:.3f}")

    if a.calibrate_only:
        print("\n[calibrate-only] tam kosum atlandi.")
        return

    # --- 3) 500 tekrar x {5000, 1117} ---
    result = {"beta_true": BETA, "phi": PHI, "err_sd": ERR_SD, "c_star": C_STAR,
              "grid": grid, "calib": {"gains": gains, "R2_target": {"measured": 0.08, "history": 0.67}},
              "panels": {}}
    for n in (5000, 1117):
        acc = {w: {k: [] for k in ("naive_crude", "naive_adj", "msm", "ess_frac")} for w in "ABC"}
        for _ in range(a.reps):
            r = np.random.default_rng(rng.integers(1 << 31))
            for w in "ABC":
                A, Y, G, Ap, U = simulate(w, cs[w], gains[w], n, r)
                e = estimators(A, Y, G, Ap)
                for k in e:
                    acc[w][k].append(e[k])
        panel = {}
        for w in "ABC":
            ess = np.array(acc[w]["ess_frac"])
            msm_applicable = bool(np.median(ess) >= 0.30)
            panel[w] = dict(
                naive_crude=summarize(acc[w]["naive_crude"]),
                naive_adj=summarize(acc[w]["naive_adj"]),
                msm=summarize(acc[w]["msm"]),
                ess_frac_median=round(float(np.median(ess)), 3),
                ess_frac_min=round(float(ess.min()), 3),
                msm_applicable=msm_applicable,
                msm_verdict=("UYGULANABILIR" if msm_applicable else "UYGULANAMAZ (ESS<%30)"))
        result["panels"][f"n={n}"] = panel
        print(f"\n===== PANEL n={n} ({a.reps} tekrar) =====")
        print(f"{'dunya':>6} {'estimator':>12} {'mean':>8} {'bias':>7} {'sd':>6} {'95%band':>16} {'ESS':>7}")
        for w in "ABC":
            for k in ("naive_crude", "naive_adj", "msm"):
                s = panel[w][k]
                ess = f"{panel[w]['ess_frac_median']:.2f}" if k == "msm" else ""
                print(f"{w:>6} {k:>12} {s['mean']:8.2f} {s['bias']:7.2f} {s['sd']:6.2f} "
                      f"[{s['lo']:6.2f},{s['hi']:6.2f}] {ess:>7}")
            print(f"{'':>6} MSM: {panel[w]['msm_verdict']}")

    (out / "three_worlds_sim.json").write_text(json.dumps(result, indent=2))
    print(f"\n[20] -> {out/'three_worlds_sim.json'}")

    # --- figur: iki panel ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        ests = ["naive_crude", "naive_adj", "msm"]
        colors = {"naive_crude": "#c44", "naive_adj": "#48c", "msm": "#2a2"}
        for ax, n in zip(axes, (5000, 1117)):
            panel = result["panels"][f"n={n}"]
            x = np.arange(3)
            for j, k in enumerate(ests):
                means = [panel[w][k]["mean"] for w in "ABC"]
                los = [panel[w][k]["mean"] - panel[w][k]["lo"] for w in "ABC"]
                his = [panel[w][k]["hi"] - panel[w][k]["mean"] for w in "ABC"]
                ax.errorbar(x + (j - 1) * 0.22, means, yerr=[los, his], fmt="o",
                            color=colors[k], capsize=3, label=k)
            ax.axhline(BETA, ls="--", color="k", lw=1, label="gercek beta=-20")
            ax.set_xticks(x); ax.set_xticklabels(["A\n(gurultu)", "B\n(olculu fb)", "C\n(gizli U)"])
            ax.set_title(f"n={n} hasta")
            ax.set_ylabel("kestirilen doz-etki (mg/dL per U/kg)")
        axes[0].legend(fontsize=8, loc="upper left")
        fig.suptitle("Uc dunya: naif vs MSM kimliklendirilebilirlik (500 tekrar, 95% band)")
        fig.tight_layout()
        fig.savefig(out / "three_worlds_sim.png", dpi=110)
        print(f"[20] -> {out/'three_worlds_sim.png'}")
    except Exception as ex:
        print(f"[20] figur atlandi: {ex}")


if __name__ == "__main__":
    main()

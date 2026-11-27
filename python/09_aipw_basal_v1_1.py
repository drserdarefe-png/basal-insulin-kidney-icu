#!/usr/bin/env python3
"""
09 v1.1 — AIPW sarmalayici (27.11.2026)
=======================================
Ozgun 09_aipw_basal.py DEGISTIRILMEZ; modul olarak yuklenir ve yalniz kovaryat
listesi guncellenir:
  pred_eq_24h, on_steroid  ->  pred_eq_active, steroid_night_active
Gerekce: guncel basal_events tablosunda pred_eq_24h sutunu yok (steroid zamanlama
duzeltmesinden sonra pred_eq_active + steroid_night_active geldi). Boylece AIPW,
doz-belirleyici modelle (07_gps_full_clean) AYNI 24 kovaryati kullanir.
Ek cikti: aipw_set_summary.json (analiz setinin TOPLU ozeti; simulasyon capalari icin).

KULLANIM:
  python 09_aipw_basal_v1_1.py <ozgun_09_aipw_basal.py> --parquet 'data/basal_events/*.parquet' --outdir <cikti>
"""
import importlib.util, json, sys
from pathlib import Path

orig = sys.argv[1]
spec = importlib.util.spec_from_file_location("aipw_orig", orig)
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

SWAP = {"pred_eq_24h": "pred_eq_active", "on_steroid": "steroid_night_active"}   # ayni konumda degistir
m.X_CONF = [SWAP.get(c, c) for c in m.X_CONF]
print(f"[AIPW v1.1] X_CONF ({len(m.X_CONF)}): {m.X_CONF}")
print(f"[AIPW v1.1] X_TITR ({len(m.X_TITR)}): {m.X_TITR}")

sys.argv = [orig] + sys.argv[2:]
m.main()

# ---- analiz setinin toplu ozeti (yalniz sayilar) ----
import argparse
ap = argparse.ArgumentParser(); ap.add_argument("--parquet", required=True); ap.add_argument("--outdir", default="qc")
a = ap.parse_args(sys.argv[1:])
df = m.load(a.parquet)
summ = dict(
    n_events=int(len(df)), n_subjects=int(df.subject_id.nunique()),
    n_first=int(df.is_first.sum()), n_titration=int((~df.is_first).sum()),
    dose_u_kg_mean=round(float(df.dose_u_kg.mean()), 4), dose_u_kg_sd=round(float(df.dose_u_kg.std()), 4),
    dose_u_kg_max_after_trim=round(float(df.dose_u_kg.max()), 4),
    n_units_gt50=int((df.units_total > 50).sum()) if "units_total" in df else None,
    units_max=float(df.units_total.max()) if "units_total" in df else None,
    nadir_mean=round(float(df.glu_nadir_night.mean()), 1), nadir_sd=round(float(df.glu_nadir_night.std()), 1),
    n_by_egfr_band={str(k): int(v) for k, v in df.egfr_band.value_counts(dropna=False).items()},
    x_conf=m.X_CONF, x_titr=m.X_TITR)
Path(a.outdir).mkdir(parents=True, exist_ok=True)
(Path(a.outdir) / "aipw_set_summary.json").write_text(json.dumps(summ, indent=2))
print(f"[AIPW v1.1] analiz seti: n={summ['n_events']}, hasta={summ['n_subjects']}, "
      f"doz SS={summ['dose_u_kg_sd']}, >50U olay={summ['n_units_gt50']}")

#!/usr/bin/env bash
# =====================================================================================
# Insulin makalesi: birincil analizlerin GUNCEL veriyle yeniden calistirilmasi (27.11.2026)
#
# PowerShell'de tek komut:
#     cd "$env:USERPROFILE\Desktop\insülin"; wsl bash ./rerun_27.11.2026.sh
#
# Ne yapar:
#   1) WSL'de projeyi bulur (data/basal_events) ve parquet'in guncel oldugunu dogrular
#      (beklenen: 7954 glarjin olayi, 1117 hasta, 386 olay > 50 U).
#   2) Dogrulama tutarsa uc analizi calistirir: AIPW, hasta-ici model (+ duyarlilik), doz-belirleyici model.
#   3) Ciktilari bu klasordeki qc_rerun_27.11.2026 altina yazar.
# Ne YAPMAZ: proje dosyalarini degistirmez; hasta duzeyinde veri kopyalamaz (yalniz JSON + PNG, toplu sayilar).
# =====================================================================================
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="$HERE/qc_rerun_27.11.2026"
mkdir -p "$OUT"
exec > >(tee "$OUT/rerun_log.txt") 2>&1
EXP_N="${EXP_N:-7954}"; EXP_PAT="${EXP_PAT:-1117}"; EXP_GT50="${EXP_GT50:-386}"
PARQ="data/basal_events/*.parquet"
echo "== BASLADI $(date '+%d.%m.%Y %H:%M:%S') | WSL kullanici: $(whoami)"

# ---------- 1) proje klasoru ----------
if [ -n "${PROJ:-}" ]; then
  CANDS=("$PROJ/data/basal_events")
else
  mapfile -t CANDS < <(find "$HOME" -maxdepth 5 -type d -path '*/data/basal_events' -not -path '*/.*/*' 2>/dev/null | sort)
  if [ "${#CANDS[@]}" -eq 0 ]; then          # WSL ev dizininde yoksa Windows kullanici klasorune bak (yavas olabilir)
    WINHOME="$(echo "$HERE" | cut -d/ -f1-5)"
    echo "   WSL ev dizininde bulunamadi; $WINHOME altinda araniyor (1 dakika surebilir)..."
    mapfile -t CANDS < <(find "$WINHOME" -maxdepth 4 -type d -path '*/data/basal_events' -not -path '*/AppData/*' -not -path '*/.*/*' 2>/dev/null | sort)
  fi
fi
if [ "${#CANDS[@]}" -eq 0 ]; then
  echo "== DURDU: \$HOME ($HOME) altinda data/basal_events klasoru bulunamadi."
  echo "   Proje baska yerdeyse:  PROJ=/proje/yolu wsl bash ./rerun_27.11.2026.sh"
  exit 2
fi

P=""; PY=""
for D in "${CANDS[@]}"; do
  C="$(dirname "$(dirname "$D")")"
  CPY="$C/.venv/bin/python"; [ -x "$CPY" ] || CPY="$(command -v python3)"
  echo "== ADAY PROJE: $C"
  echo "   python: $CPY | parquet dosyasi: $(ls "$D"/*.parquet 2>/dev/null | wc -l) adet, en yeni: $(ls -t --time-style=+%d.%m.%Y_%H:%M -l "$D"/*.parquet 2>/dev/null | head -1 | awk '{print $6}')"
  if (cd "$C" && "$CPY" - "$EXP_N" "$EXP_PAT" "$EXP_GT50" <<'PYEOF'
import sys, importlib
exp_n, exp_pat, exp_gt50 = map(int, sys.argv[1:4])
try:
    import duckdb
except Exception as e:
    print("   PAKET EKSIK duckdb:", type(e).__name__); sys.exit(4)
con = duckdb.connect(); src = "read_parquet('data/basal_events/*.parquet')"
cols = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {src}").fetchall()]
need = ["units_total", "stay_id", "dose_time", "hrs_since_prev_glargine", "prev_glargine_dose_u_kg",
        "pred_eq_active", "steroid_night_active", "glu_nadir_night", "weight_kg", "egfr_now", "event_hr"]
miss = [c for c in need if c not in cols]
print(f"   sutun: {len(cols)} | pred_eq_24h var mi: {'pred_eq_24h' in cols} | eksik gerekli sutun: {miss or 'yok'}")
if miss: sys.exit(5)
n, pat, gt50, mx, nad = con.execute(f"""SELECT COUNT(*), COUNT(DISTINCT subject_id),
    SUM(CASE WHEN units_total > 50 THEN 1 ELSE 0 END), MAX(units_total), COUNT(glu_nadir_night)
    FROM {src} WHERE basal_type='glargine'""").fetchone()
print(f"   SAYIM: glarjin_olay={n} hasta={pat} olay_gt50U={gt50} max_U={mx} gece_glukozlu={nad}")
bad = []
for mod in ("lightgbm", "shap", "sklearn", "statsmodels", "pandas", "numpy", "matplotlib"):
    try:
        v = importlib.import_module(mod).__version__
    except Exception as e:
        v = "EKSIK"; bad.append(mod)
    print(f"   paket {mod}: {v}")
if bad: print("   PAKET EKSIK:", bad); sys.exit(4)
sys.exit(0 if (n, pat, gt50) == (exp_n, exp_pat, exp_gt50) else 3)
PYEOF
  ); then P="$C"; PY="$CPY"; echo "   -> DOGRULANDI"; break
  else echo "   -> bu aday uygun degil (cikis kodu $?: 3=sayim farkli, 4=paket eksik, 5=sutun eksik)"; fi
done
if [ -z "$P" ]; then
  echo "== DURDU: beklenen sayilara (glarjin_olay=$EXP_N hasta=$EXP_PAT olay_gt50U=$EXP_GT50) uyan guncel parquet bulunamadi."
  echo "   Analiz CALISTIRILMADI. Yukaridaki SAYIM satirini Claude'a iletin."
  exit 3
fi

# ---------- betikler: proje kopyasi varsa o, yoksa bu klasordeki kopya ----------
pick() { if [ -f "$P/python/$1" ]; then echo "$P/python/$1"; else echo "$HERE/python/$1"; fi; }
S07="$(pick 07_gps_diagnostic.py)"; S09="$(pick 09_aipw_basal.py)"; S10="$(pick 10_within_patient.py)"
for f in 07_gps_diagnostic.py 09_aipw_basal.py 10_within_patient.py; do
  a="$( [ -f "$P/python/$f" ] && md5sum "$P/python/$f" | cut -c1-8 || echo yok )"; b="$(md5sum "$HERE/python/$f" 2>/dev/null | cut -c1-8)"
  echo "   betik $f: proje=$a masaustu=$b $( [ "$a" = "$b" ] && echo AYNI || echo FARKLI )"
done

cd "$P" || exit 9
STATUS=""
run() { local name="$1"; shift; echo; echo "== CALISIYOR: $name"; local t0=$SECONDS
  if "$@"; then STATUS="$STATUS $name=TAMAM"; echo "== $name TAMAM ($((SECONDS - t0)) sn)"
  else STATUS="$STATUS $name=HATA"; echo "== $name HATA (cikis kodu $?)"; fi; }

CONF24="prev_glargine_dose_u_kg,hrs_since_prev_glargine,dose_24h_cum,weight_kg,cho_g_hr,glu_min_24h,iob_bolus_u,event_hr,egfr_now,dm_any,glu_mean_24h,glu_now,cho_24h_cum,age,bun,glu_sd_24h,kdigo_stage,pred_eq_active,on_vasopressor,dm_t1,child_pugh_score,on_crrt,steroid_night_active,cirrhosis"

run AIPW        "$PY" "$HERE/python/09_aipw_basal_v1_1.py" "$S09" --parquet "$PARQ" --outdir "$OUT"
run HASTA_ICI   "$PY" "$S10" --parquet "$PARQ" --outdir "$OUT"
run HASTA_ICI_DUYARLILIK "$PY" "$HERE/python/10b_within_patient_sens.py" "$S10" --parquet "$PARQ" --outdir "$OUT"
run DOZ_BELIRLEYICI "$PY" "$S07" --parquet "$PARQ" --outdir "$OUT" --out-name 07_gps_full_clean --clip-q 0.995 \
    --where "basal_type='glargine' AND prev_glargine_dose_u_kg IS NOT NULL" --confounders "$CONF24"

echo
echo "== CIKTI DOSYALARI ($OUT):"
ls -l --time-style=+%H:%M "$OUT" | awk 'NR>1{print "   ", $5, $6, $7}'
echo "== DURUM:$STATUS"
case "$STATUS" in *HATA*) echo "== BITTI (HATALI): yukaridaki HATA satirlarini Claude'a iletin.";; *) echo "== BITTI: 4 adim tamam. Claude'a 'bitti' yazmaniz yeterli.";; esac

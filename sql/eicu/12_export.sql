SELECT patientunitstayid, hospitalid, dose_u_kg, dose_u, weight, scr, egfr, egfr_band,
       dm, age, age_topcoded, order_seq, is_first_dose, prev_dose_u_kg, hrs_since_prev
FROM `uludag-mimic-insulin.insulin.eicu_glargine_events`

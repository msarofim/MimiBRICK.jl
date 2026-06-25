# BRICK-FM Observational Discrepancies

Known discrepancies between BRICK-FM posterior output and observational products,
with quantified magnitudes, root causes, and status. Intended for collaborators
reviewing PR2.

## Summary

| Component | Discrepancy | Magnitude | Root cause | Status |
|-----------|-------------|-----------|------------|--------|
| Thermal expansion | Overshoots steric obs | +0.51 cm vs NOAA steric @2025 | SNEASY MAP OHC pre-1971 ramp inflates `te_α` calibration | Known; calibration fix available |
| GSIC | Undershoots 1900 cumulative | −4 cm vs Frederikse (−3.25 vs −7.27 cm) | Mengel melt rate ∝ T-level, not warming rate | Structural emulator limit |
| GIS | Undershoots obs at 1900 | ~0.53 vs 2.1 cm (Frederikse), 1.45 cm (IGCC) | FaIR GMST ~0.1 °C cooler than IGCC historically | Forcing-driven; not a BRICK issue |
| Total GMSL | Slightly low vs Dangendorf 2024 | Consistent with GMST gap | Same as GIS | Forcing-driven |
| AIS high-forcing tail | Conservative at SSP5-8.5 | — | MICI not in BRICK AIS formulation | Known BRICK limitation |

---

## 1. Thermal Expansion (TE) Overshoot vs Steric Obs

**Magnitude:** At 2025, BRICK TE overshoots NOAA steric by approximately +0.51 cm.

**Root cause:** The FM posterior calibrates `te_α` to ~0.164, roughly 3× the
physics-based value in Tony's original calibration (~0.057). This traces to the
FaIR-mean OHC forcing used in calibration: FaIR's pre-1971 OHC trajectory
carries a large positive ramp (~+35 ZJ above 1850) relative to modern observational
products (~+14 ZJ, e.g. Cheng IAPv4.2). BRICK compensated by inflating `te_α`
to match the high OHC input. When driven by FaIR or modern obs (lower OHC), TE
is consequently overpredicted.

Note: the IGCC multi-product OHC average is closer to Gouretski 2007 (~+26 ZJ)
than to Cheng — so "overshoot vs steric" depends on which obs product is used.
The Cheng product is the low-side outlier; against Gouretski the overshoot
largely disappears.

**Fix:** Recalibrate with a modern OHC product (IGCC or Zanna+Cheng splice
anchored to FaIR mean). No model changes needed.

**Status:** Known and accepted for this PR. Flagged for PR3/recalibration.

---

## 2. GSIC Undershoot at 1900

**Magnitude:** GSIC cumulative at 1900: BRICK −3.25 cm vs Frederikse −7.27 cm
(undershoot ~4 cm).

**Root cause:** The Mengel emulator parameterizes melt rate as proportional to
the temperature *level* relative to equilibrium. It cannot simultaneously produce
the fast early-20th-century melt rate (response to a large LIA disequilibrium)
and the slower modern rate without the calibrated timescales (`tau_fast`, `tau_slow`)
hitting physical limits. The two-timescale structure is designed to address exactly
this tension, but the calibrated parameter values still leave a gap at 1900.

**Status:** Structural limitation of this class of emulator. GSIC projections
remain physically reasonable. The gap is largest in the early historical period
and does not represent a comparable bias in future projections.

---

## 3. GIS Undershoot vs Observations

**Magnitude:** FaIR-driven GIS at 1900: ~0.53 cm vs Frederikse 2.1 cm and IGCC obs 1.45 cm.

**Root cause:** GIS in BRICK scales with GMST only (OHC input is irrelevant for
this component). The gap is driven by FaIR GMST running ~0.1 °C cooler than IGCC
observational estimates over the historical period — not a BRICK calibration issue.

**Status:** Accepted. Would narrow if FaIR were calibrated to IGCC rather than
the CMIP ensemble mean. No BRICK recalibration needed.

---

## 4. Total GMSL vs Dangendorf 2024

**Magnitude:** FaIR-forced BRICK GMSL runs slightly low relative to Dangendorf 2024.

**Root cause:** Directly consistent with the GMST gap (item 3). GIS, GSIC, and TE
all scale with GMST, so a cooler FaIR forcing drives lower total SLR.

**Status:** Not a BRICK issue. Resolves if the FaIR forcing is updated.

---

## 5. AIS High-Forcing Tail (SSP5-8.5)

**Magnitude:** 2100 upper-tail AIS projections under SSP5-8.5 are conservative
relative to process-model studies that include marine ice cliff instability (MICI).

**Root cause:** BRICK's Antarctic Ice Sheet component does not include MICI. At
very high warming rates the model cannot produce the rapid retreat that some
process models project.

**Status:** Known limitation of BRICK's AIS formulation, not specific to the FM
calibration. Documented in the BRICK literature. Users needing MICI-inclusive
projections should combine with dedicated ice-sheet models.

---

## Calibration Targets

### Historical calibration (1900–2018, Frederikse 2020 + extensions)

| Target | Product | Period |
|--------|---------|--------|
| AIS | Frederikse 2020 component reconstruction | 1900–2018 |
| GIS | Frederikse 2020 component reconstruction | 1900–2018 |
| GSIC | Frederikse 2020 + Dyurgerov 2003 point term | 1900–2018 / 1961–2003 |
| Steric (TE) | Frederikse 2020 steric component | 1900–2018 |
| Total GMSL | Dangendorf 2024 tide gauge + altimetry | 1900–2018 |
| AIS modern rate | IMBIE 1992–2017 point term | 1992–2017 |

### Post-2018 extension

| Target | Product | Period |
|--------|---------|--------|
| AIS | GRACE-FO JPL mascon | 2018–2024 |
| GIS | GRACE-FO JPL mascon | 2018–2024 |
| GSIC | GlaMBIE satellite glacier mass balance | 2018–2024 |
| Steric | NOAA NCEI | 2018–2024 |
| Total GMSL | NOAA STAR altimetry | 2018–2024 |

# BRICK-FM Observational Discrepancies

Known discrepancies between BRICK-FM posterior output and observational products,
with quantified magnitudes, root causes, and status. Intended for collaborators
reviewing PR2.

## Summary

| Component | Discrepancy | Magnitude | Root cause | Status |
|-----------|-------------|-----------|------------|--------|
| Thermal expansion | Overshoots steric obs | +0.51 cm vs NOAA steric @2025 | Unconstrained 1900–1953 OHC window; persists after post-2018 re-fit | Known; calibration fix available |
| GSIC | Undershoots 1900 cumulative | −4 cm vs Frederikse (−3.25 vs −7.27 cm) | Mengel melt rate ∝ T-level, not warming rate | Structural emulator limit |
| GIS | Undershoots obs at 1900 | ~0.53 vs 2.1 cm (Frederikse), 1.45 cm (IGCC) | FaIR GMST ~0.1 °C cooler than IGCC historically | Forcing-driven; not a BRICK issue |
| Total GMSL | Slightly low vs Dangendorf 2024 | Consistent with GMST gap | Same as GIS | Forcing-driven |
| AIS high-forcing tail | Upper tail bounded | — | MICI not in BRICK AIS formulation | Known BRICK limitation |

---

## 1. Thermal Expansion (TE) Overshoot vs Steric Obs

**Magnitude:** At 2025, BRICK TE overshoots NOAA steric by approximately +0.51 cm.

**Root cause:** The FM `te_α` (~0.164) and Wong et al.'s original value (~0.057)
differ because they were calibrated against different OHC forcing trajectories.
Wong calibrated against SNEASY's internal OHC, which rises ~+35 ZJ over
1900–1971 — substantially larger than modern observation-anchored products over
the same period (~+14 ZJ). The FM calibrates against the FaIR mean, which is
close to the IGCC 2024 multi-product compilation (~+37–38 ZJ over 1971–2018).

The +0.51 cm overshoot vs NOAA steric at 2025 persists even after re-fitting
against post-2018 NOAA steric data (te_α shifted only 0.164→0.159). Switching
OHC forcing products would not resolve this: FaIR OHC and IGCC already agree
within ~3% over 1971–2018. The residual reflects a product-consistency gap
between FaIR OHC and NOAA thermosteric that `te_α` as a single scalar cannot
bridge.

**Fix:** A richer TE parameterization, or an OHC forcing trajectory constructed
to be internally consistent with the thermosteric calibration target.

**Status:** Known and accepted for this PR. Likely structural within the current
TE formulation.

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

**Magnitude:** BRICK-FM AIS at 2100 under SSP2-4.5 is ~43 cm (median) vs
MAGICC-Nauels ~11 cm — BRICK-FM runs higher in the scenario level. For pulse
marginals (SC-CO2), the ordering reverses: MAGICC AIS marginal T-sensitivity
is ~6× higher than BRICK's at 2100 (see `brick_fm_changes_and_rationale.md`
section C).

**Root cause:** BRICK does not include marine ice cliff instability (MICI),
which formally bounds the extreme upper tail. The Mengel committed-melt
mechanism drives a large time-integrated AIS contribution that produces a
high scenario level; MAGICC-Nauels has higher instantaneous T-sensitivity.
The two emulators represent different physical mechanisms, not a
conservative-vs-aggressive ordering.

**Status:** Known limitation of BRICK's AIS formulation. Users needing
MICI-inclusive projections should combine with dedicated ice-sheet models.

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

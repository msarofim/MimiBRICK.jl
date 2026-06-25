# BRICK-FM: Changes from the Original MimiBRICK Calibration and Their Rationale

This document describes each methodological departure from the original MimiBRICK
(Wong et al. 2022) calibration, the physical or statistical motivation for the
change, and known remaining limitations flagged as candidates for future work.

---

## Changes

### 1. Mengel-2016 two-timescale glacier emulator (replaces Wigley-Raper-Bakker)

**Change:** The default single-reservoir glaciers-and-small-ice-caps component
(Wigley-Raper-Bakker, WRB) is replaced by the Mengel et al. 2016 (PNAS 113:2597)
temperature-dependent-equilibrium emulator with two relaxation timescales.

**Rationale:** The WRB formulation leads to the behavior that sustained
temperature above the equilibrium threshold eventually melts the *entire* glacier
reservoir. Physical glaciers and ice caps retain a temperature-appropriate remnant
at any finite warming level (Marzeion et al. 2012). This matters for projections
because WRB over-depletes the reservoir under sustained future warming, suppressing
late-21st-century and post-2100 glacier contributions.

The Mengel emulator uses a saturating equilibrium S_eq(T) = a(1 − exp(−b(T − T_LIA))),
which asymptotes to a finite maximum `a`. A sustained warming T* commits only S_eq(T*)
< a, preserving a temperature-appropriate remnant. The two timescales (τ_fast ≈ 40 yr,
τ_slow ≈ 250 yr) reflect the observed glacier size distribution: small, low-elevation
glaciers respond quickly; large, high-altitude ice caps respond slowly. A single
timescale cannot simultaneously fit the rapid early-20th-century discharge and the
slower modern rate.

The LIA-disequilibrium offset (T_LIA < 0, relative to 1850-1900) gives glaciers a
non-zero equilibrium contribution at the start of the simulation, directly simulating
the "committed" melt from post-LIA warming that predates the forcing window — without
requiring an external natural-melt budget or an anthropogenic/natural forcing split.

---

### 2. Dangendorf 2024 total GMSL (replaces Church & White 2011)

**Change:** The total GMSL calibration target is updated from Church & White (2011)
to Dangendorf et al. (2024).

**Rationale:** Church & White (2011) relied on tide-gauge networks without GPS-based
vertical land motion (VLM) corrections, introducing spatially correlated errors from
subsidence and uplift at gauge sites. Dangendorf et al. (2024) applies GPS-derived
VLM corrections to the full tide-gauge record and uses an improved spatial
interpolation scheme, producing a more accurate reconstruction of the global-mean
signal. The Dangendorf reconstruction also has better coverage of the Southern
Hemisphere and updated uncertainty quantification.

---

### 3. Frederikse 2020 component data (replaces older component reconstructions)

**Change:** Per-component calibration targets (AIS, GIS, GSIC, steric) are taken
from Frederikse et al. (2020, Nature 584:70), which provides a closed sea-level
budget for 1900–2018 using consistent data sources and uncertainty propagation.

**Rationale:** Frederikse 2020 is the first reconstruction to close the global
sea-level budget from individual components over the full 20th century. The
consistent treatment across components — rather than assembling targets from
disparate sources — reduces the risk of double-counting or mismatched reference
periods. The closed-budget property provides an implicit consistency check: the
sum of components matches the total GMSL within uncertainty, so the individual
component likelihoods and the total GMSL likelihood (Dangendorf) are not
independent — they reinforce each other rather than conflicting.

---

### 4. Post-2018 observational extension (GRACE-FO, GlaMBIE, NOAA)

**Change:** The calibration window is extended from 1900–2018 to 1900–2026 using
modern satellite products for the post-2018 period.

**Rationale:** A longer calibration window reduces posterior uncertainty and
provides stronger constraints on parameters that govern future projections. The
post-2018 period is also the most dynamically active in the observational record —
accelerated Greenland and Antarctic mass loss, continued glacier retreat — and
including it anchors the model's response to the warming regime most relevant for
near-term projections.

Sources: GRACE-FO JPL mascon (AIS, GIS), GlaMBIE 2025 (GSIC), NOAA NCEI
thermosteric (steric), NOAA STAR altimetry (total GMSL, spliced to Dangendorf at
2018).

---

### 5. Calibration driven by FaIR v2.2 (replaces SNEASY internal forcing)

**Change:** The historical GMST and OHC trajectories used to drive BRICK during
calibration are taken from the FaIR v2.2 ensemble mean rather than SNEASY's own
internal output.

**Rationale:** The original MimiBRICK calibration was conducted with SNEASY
(Simple Nonlinear EArth SYstem model) as the coupled climate driver, so the BRICK
parameters were implicitly conditioned on SNEASY's historical temperature and OHC
trajectories. Using the same model for calibration and projection is internally
consistent, but it makes the calibrated BRICK posterior specific to SNEASY's
climate response — the posterior does not transfer cleanly to other climate drivers.

FaIR v2.2 is more widely used in climate economics and impact modeling (EPA SC-GHG
work, RFF-SP ensemble, IPI, GIVE) and provides OHC directly as an output, making
it the natural choice for studies that couple BRICK to FaIR-driven ensembles. By
calibrating BRICK against the FaIR mean trajectory, the posterior is conditioned on
a climate forcing that is consistent with the projection runs.

Note: the FaIR mean trajectory used for calibration runs approximately 0.1 °C
below IGCC observational estimates at 2024. This reflects the emissions scenario
(RFF-SP draws follow SSP2-4.5-like paths; real-world emissions post-2015 ran
warmer due to faster aerosol reductions and higher CH4), not an intrinsic FaIR
model bias. The effect slightly depresses GIS and GSIC contributions relative to
observation-forced runs; see `brick_fm_obs_discrepancies.md` for details.

---

### 6. Deterministic land-water storage for reproducible ensemble runs

**Change:** `create_brick_fair` (the recommended FM entry point) defaults to
`lws=:central` (deterministic 0.3 mm/yr mean rate). `get_model` retains
`lws=:random` for backward compatibility with upstream.

**Rationale:** The original `get_model` draws LWS from N(0.0003, 0.00018) m/yr
unseeded on every call, making results irreproducible build-to-build and
representing LWS uncertainty by a single arbitrary realization rather than
propagating it through the ensemble. For ensemble runs (the standard FM usage),
a single fixed LWS realization is appropriate; the LWS uncertainty (~0.16 cm by
2100) is small relative to the AIS/posterior spread. The `:central` option fixes
this at the distribution mean without requiring users to manage external seeds.

---

## Targets for Future Improvement

### A. OHC product for thermal expansion calibration

The FM posterior calibrates `te_α` to ~0.164 (≈3× Tony's original 0.057) because
the FaIR mean OHC used for calibration has a large pre-1971 positive ramp
(~+35 ZJ above 1850) relative to modern observational products (~+14 ZJ Cheng
IAPv4.2, ~+26 ZJ Gouretski/IGCC). BRICK compensates by inflating `te_α`.

**Fix:** Recalibrate using a modern OHC splice (IGCC multi-product mean, or
Zanna + Cheng anchored to the FaIR mean at the calibration start). This would
bring `te_α` closer to the physics-based value and reduce the ~+0.5 cm TE
overshoot at 2025.

### B. GSIC structural undershoot at 1900

The Mengel emulator still undershoots the Frederikse GSIC contribution at 1900
by ~4 cm (model −3.25 cm vs obs −7.27 cm). This is a structural limitation: the
melt rate scales with temperature *level* rather than warming *rate*, so the model
cannot simultaneously fit the rapid post-LIA discharge of the early 20th century
and the slower modern rate without pushing `tau_fast` to its physical lower bound.

**Potential fix:** A three-timescale emulator, or an explicit LIA committed-melt
budget added to the two-timescale structure. Alternatively, a Marzeion-style
natural melt baseline (separate from the temperature-forced term) could absorb
the early-century residual.

### C. MICI in Antarctic Ice Sheet projections

BRICK's AIS component does not include marine ice cliff instability (MICI). Under
high-forcing scenarios (SSP5-8.5, late 21st century), the upper tail of AIS
projections is consequently conservative relative to process-model studies that
include MICI.

**Path forward:** This is a known BRICK limitation (discussed in Wong et al. 2022)
and would require structural changes to the AIS component, likely meriting its own
PR.

### D. GIS undershoot from FaIR GMST bias

FaIR GMST runs ~0.1 °C cooler than IGCC observational estimates historically,
causing a ~0.5–1 cm GIS undershoot at 2024. This is a FaIR property, not a
BRICK calibration issue, but it propagates into BRICK-FM projections.

**Fix:** No BRICK change needed. Resolves if FaIR calibration is updated toward
IGCC or if a GMST bias correction is applied upstream.

### E. LWS uncertainty propagation

LWS uncertainty (~0.16 cm by 2100) is currently either ignored (`:central`) or
represented by a single unseeded draw (`:random`). Proper propagation would
require either a per-draw LWS sample (adding one dimension to the ensemble) or
an importance-weighted LWS treatment.

**Assessment:** LWS is a minor term relative to AIS/GIS uncertainty; this is low
priority unless the analysis specifically targets near-term SLR where LWS is a
larger fractional contributor.

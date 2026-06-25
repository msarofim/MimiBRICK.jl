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

### A. Residual thermal expansion overshoot

The FM posterior calibrates `te_α` to ~0.164, roughly 3× Wong et al.'s
original value of 0.057. The two values reflect different OHC forcing
assumptions during calibration. Wong calibrated against SNEASY's internal
OHC trajectory, which rises ~+35 ZJ over 1900–1971 — substantially larger
than modern observation-anchored products over the same period (~+14 ZJ for
FaIR mean, Zanna+Cheng). The FM calibrates against the FaIR mean, which lies
close to the IGCC 2024 multi-product compilation (Palmer & von Schuckmann;
+38 vs +37 ZJ over 1971–2018). Among modern products, Cheng et al. 2024 IAPv4.2 (ESSD 16:3517, +31 ZJ)
sits on the low side; IGCC and FaIR are in close agreement.

A residual ~+0.5 cm TE overshoot vs NOAA steric observations persists at 2025,
even after the posterior was re-fit against post-2018 NOAA steric data (which
changed te_α only slightly, 0.164→0.159). The overshoot appears to be driven
by the 1900–1953 period, where no direct ocean heat observations constrain
the calibration and the model accumulates more OHC than the observation-anchored
products suggest.

**Potential fix:** Extend the calibration OHC target back to 1850 using a
pre-ARGO reconstruction (e.g. Zanna 2019 spliced to IGCC at 1971), which
would better constrain `te_α` over the full historical window including the
early-century period.

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

### C. Antarctic Ice Sheet projections and MICI

BRICK's AIS component does not include marine ice cliff instability (MICI).
However, BRICK-Mengel does not produce conservative AIS projections overall:
in a direct comparison against the MAGICC-Nauels 2025 emulator (SSP2-4.5,
600-member AR6 drawnset), BRICK-Mengel AIS @2100 is ~43 cm vs MAGICC ~11 cm
at the median, and BRICK's p95 total SLR (108 cm) exceeds MAGICC's (87 cm).
The Mengel two-timescale committed-melt mechanism drives a large time-integrated
AIS contribution that more than offsets the absence of MICI.

There is a level-vs-marginal inversion: for pulse experiments (SC-CO2), MAGICC's
AIS *marginal* T-sensitivity is ~6× higher than BRICK-Mengel's at 2100, meaning
MAGICC assigns a higher SC-CO2 despite a lower scenario AIS level. The two
emulators represent genuinely different physical mechanisms (MAGICC: high
instantaneous T-sensitivity; BRICK-Mengel: high committed slow-timescale melt),
not simply a conservative-vs-aggressive ordering.

Formally, MICI absence does place a structural bound on BRICK's extreme upper tail
under very high forcing, and this remains a limitation relative to process models
that include it. Adding MICI would require structural changes to the AIS component
and would merit its own PR.

### D. GIS undershoot from emissions-scenario GMST gap

The FaIR mean trajectory used for calibration runs ~0.1 °C below recent IGCC
observational estimates at 2024, causing a ~0.5–1 cm GIS undershoot at that
date. This gap reflects the emissions scenario (RFF-SP draws follow SSP2-4.5-like
paths; real-world emissions post-2015 ran warmer due to faster aerosol reductions
and higher CH4) rather than an intrinsic FaIR model bias — all calibration
versions show the same gap when fed the same emissions.

**Fix:** No BRICK change needed. Resolves if the calibration forcing is updated
to use observed historical emissions (e.g. Smith 2024 splice) rather than an
SSP2-4.5-based trajectory.

### E. LWS uncertainty propagation

LWS uncertainty (~0.16 cm by 2100) is currently either ignored (`:central`) or
represented by a single unseeded draw (`:random`). Proper propagation would
require either a per-draw LWS sample (adding one dimension to the ensemble) or
an importance-weighted LWS treatment.

**Assessment:** LWS is a minor term relative to AIS/GIS uncertainty; this is low
priority unless the analysis specifically targets near-term SLR where LWS is a
larger fractional contributor.

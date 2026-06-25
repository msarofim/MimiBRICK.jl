# Calibration Guide: BRICK-Mengel

## When to re-run the calibration

The shipped posterior subsample (`data/MimiBRICK/parameters_subsample_brick_mengel.csv`)
is the output of the full MCMC calibration described below. You need to re-run it if:

- You change the **observational targets** (e.g., update to a new Frederikse/Dangendorf
  version, add post-2024 data).
- You change the **FaIR forcing** used to drive calibration (currently FaIR v2.2 SSP2-4.5
  mean trajectory; stored in `data/observations/fair_mean_{gmst,ohc}.csv`).
- You change the **prior structure** or free-parameter set in `calibrate_mcmc_mengel.jl`.
- You change the **Mengel component** equations (e.g., add a third timescale).

You do NOT need to re-run for:
- Applying the posterior to a different SSP scenario (use `create_brick_fair` with
  the relevant scenario forcing; the posterior is scenario-independent).
- Changing `start_year`/`end_year` for projections.
- Updating the Mengel `gic_sl0` initial condition or other non-free parameters.

## Calibration workflow

The calibration runs in three steps. Production runs need a multi-core machine
(NYU Torch `cs` partition recommended; see the `nyu-torch-hpc` skill for setup).

### Step 1: MAP point estimate

```
julia --project=. calibration/calibrate_full_joint_mengel.jl
```

Outputs `outputs/calib_full_joint_params.csv` — the maximum a posteriori (MAP)
parameter vector. Used as the MCMC starting point. Runtime: ~5 minutes local.

### Step 2: MCMC chains (production)

Run one chain per seed. Recommended: ≥4 seeds × ≥500 000 iterations on Torch.
The default targets are the 1900–2026 extended set (GRACE-FO/GlaMBIE/NOAA post-2018);
pass `--base` to use the 1900–2018 Frederikse-only targets instead.

```bash
# Local smoke test (2000 iter, ~2 min):
julia --project=. calibration/calibrate_mcmc_mengel.jl 2000 2026

# Base targets only (1900-2018, Frederikse 2020):
julia --project=. calibration/calibrate_mcmc_mengel.jl 2000 2026 --base

# Production (example SLURM array, seeds 1-4, 500k iter each):
sbatch calibration/run_mcmc_torch.sbatch
```

Each run writes `outputs/mcmc/chain_seed<SEED>_n<NITER>.csv`. If a previous
run produced `outputs/mcmc/adapted_cov.csv`, it is automatically used to seed
the proposal (far better mixing than the diagonal start).

Key acceptance-rate target: 0.10–0.40 (RAM targets 0.234). Chains below 0.05
or above 0.50 suggest a misspecified proposal or a likelihood issue.

### Step 3: Convergence diagnostics + posterior subsample

```
julia --project=. calibration/postprocess_mcmc_mengel.jl [n_subsample]
```

Reads all `outputs/mcmc/chain_seed*.csv`, burns the first half of each chain,
computes Gelman–Rubin R̂ and ESS per parameter, and writes a thinned subsample
to `data/MimiBRICK/parameters_subsample_brick_mengel.csv`.

Convergence targets: R̂ < 1.05 and ESS > 400 for all parameters. Parameters
flagged as not converged require longer chains.

Also writes `outputs/mcmc/adapted_cov.csv` — the empirical posterior covariance —
which seeds the next run's proposal for faster mixing.

## Free parameters (28 total)

18 physical parameters are free in the MCMC. The remaining BRICK shape parameters
(AIS geometry, precipitation) are fixed at the medoid of the prior to avoid
unidentifiable directions in parameter space.

| Group | Free params |
|-------|-------------|
| Antarctic Ice Sheet (key dynamics) | `ais_ocean_temperature₀`, `ais_α`, `ais_ν`, `temperature_threshold` |
| Antarctic Ocean | `anto_α`, `anto_β` |
| Greenland Ice Sheet | `greenland_a`, `greenland_b`, `greenland_α`, `greenland_β`, `greenland_v₀` |
| Thermal expansion | `te_α` |
| Mengel glacier (2-τ) | `gic_a`, `gic_b`, `gic_T_lia`, `gic_f`, `gic_tau_fast`, `gic_tau_slow` |
| AR(1) noise (per series) | `sd_*`, `rho_*` for AIS / GSIC / GIS / steric / total (10 params) |

## Notes on the `precip_log` convention

BRICK v2.0.0 reparameterized `ais_precipitation₀` to log-space (the model
computes `exp(ais_precipitation₀) * exp(κ·T)` internally). The shipped posterior
was calibrated WITH v2.0.0, so `ais_precipitation₀` is stored in log-space and
should be passed directly to v2.0.0 models — no `precip_log` shim needed.

If you apply the BRICK-Mengel posterior to a v1.x model (not recommended), you
would need to exponentiate `ais_precipitation₀` before passing it.

## Observation files

All calibration targets are in `data/observations/calibration_targets_brick_mengel.csv`.
The format is one row per year with columns: `year`, `ais`, `ais_lo`, `ais_hi`,
`gsic`, `gsic_lo`, `gsic_hi`, `gis`, `gis_lo`, `gis_hi`, `steric`, `steric_lo`,
`steric_hi`, `dang`, `dang_sig`, `lws`, `lws_lo`, `lws_hi`.

All values in cm relative to the 1995–2005 mean. Uncertainty bands are 5–95%.

## ============================================================================
## calibrate_mcmc_mengel.jl  —  Bayesian MCMC calibration of BRICK-Mengel
##
## Promotes the MAP point estimate to a full posterior using Robust Adaptive
## Metropolis (RAM) sampling with heteroscedastic AR(1) likelihoods per
## fitted series (Ruckert et al. 2017, following MimiBRICK's approach).
##
## Model: FaIR v2.2 SSP2-4.5 mean GMST + OHC → BRICK v2.0.0 with the
## Mengel-2016 two-timescale glacier emulator.
##
## Free parameters (28 total):
##   Physical (18): ais_ocean_temperature₀, ais_α, ais_ν, temperature_threshold,
##     anto_α, anto_β, greenland_{a,b,α,β,v₀}, te_α,
##     gic_{a,b,T_lia,f,tau_fast,tau_slow}
##   AR(1) noise (10): (sd, rho) × {ais, gsic, gis, steric, total}
## Weakly-constrained DAIS shape parameters are fixed at the prior medoid.
##
## Calibration targets (1900–2026 default, re-referenced to 1995–2005):
##   Frederikse 2020 AIS/GSIC/GIS/steric 1900–2018 + post-2018 extension:
##     GRACE-FO AIS/GIS, GlaMBIE GSIC, NOAA steric, NOAA STAR total (AR(1) likelihood)
##   Dangendorf 2024 total GMSL (AR(1) likelihood; LWS uncertainty folded in)
##   IMBIE 1992–2017 AIS rate point constraint (Gaussian)
##   Dyurgerov 1961–2003 GSIC rate point constraint (Gaussian)
##   Pass --base to use the 1900–2018 base targets only (Frederikse 2020 period).
##
## See docs/calibration_guide_mengel.md for when to re-run and runtime guidance.
##
## Usage:
##   julia --project=. calibration/calibrate_mcmc_mengel.jl [n_iter] [seed] [--base]
##   n_iter : MCMC iterations (default 2000 for smoke test; 500000 for production)
##   seed   : random seed (default 2026; use different seeds for parallel chains)
##   --base : use 1900–2018 base targets instead of the 1900–2026 extended default
## ============================================================================

using CSV, DataFrames, Mimi, MimiBRICK, Statistics, LinearAlgebra, Distributions, Random, Printf
using RobustAdaptiveMetropolisSampler

const REPO    = abspath(joinpath(@__DIR__, ".."))
const OBS_DIR = joinpath(REPO, "data", "observations")
const OUT_DIR = joinpath(REPO, "outputs", "mcmc")

N_ITER   = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
SEED     = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2026
USE_BASE = "--base" in ARGS      # pass --base to use 1900-2018 targets only

const Y0     = 1850
const Y1     = USE_BASE ? 2018 : 2026   # extended targets end at 2026
const B0, B1 = 1995, 2005              # re-reference window

years = collect(Y0:Y1)
idx(y) = findfirst(==(y), years)
ib     = [idx(y) for y in B0:B1]

# ---- AR(1) heteroscedastic log-likelihood (Ruckert et al. 2017) ----
function hetero_logl_ar1(res::Vector{Float64}, σ::Float64, ρ::Float64, ϵ::Vector{Float64})
    n  = length(res)
    σp = σ^2 / (1 - ρ^2)
    H  = abs.(collect(1:n)' .- collect(1:n))
    Σ  = σp .* ρ.^H .+ Diagonal(ϵ.^2)
    return logpdf(MvNormal(Symmetric(Σ)), res)
end

# ---- load forcing + calibration targets ----
function lc(path, col)
    d = CSV.read(path, DataFrame)
    Dict(Int(d[i, "year"]) => Float64(d[i, col]) for i in 1:nrow(d))
end

gmst = [lc(joinpath(OBS_DIR, "fair_mean_gmst.csv"), "gmst_C")[y]   for y in years]
ohc  = [lc(joinpath(OBS_DIR, "fair_mean_ohc.csv"),  "ohc_1e22J")[y] for y in years]

targets_file = USE_BASE ? "calibration_targets_brick_mengel.csv" :
                           "calibration_targets_brick_mengel_ext.csv"
targets_path = joinpath(OBS_DIR, targets_file)
isfile(targets_path) || error("Calibration targets not found: $targets_path\n" *
    "See data/observations/README.md for the expected format.")
tg = CSV.read(targets_path, DataFrame)
tgi(y) = findfirst(==(y), tg.year)
println("Calibration targets: $targets_file  ($(USE_BASE ? "1900-2018 base" : "1900-2026 extended"))")

FY  = collect(1900:Y1)    # fitted years match the targets file
fyi = [tgi(y) for y in FY]
myi = [idx(y)  for y in FY]

# Per-year obs uncertainty (floor at 0.05 cm to prevent likelihood collapse)
ϵband(lo, hi) = max.((hi .- lo) ./ (2 * 1.645), 0.05)

obs = (
    ais    = Float64.(tg.ais[fyi]),
    gsic   = Float64.(tg.gsic[fyi]),
    gis    = Float64.(tg.gis[fyi]),
    steric = Float64.(tg.steric[fyi]),
    dang   = Float64.(tg.dang[fyi]),
    lws    = Float64.(tg.lws[fyi]),
)
ϵ = (
    ais    = ϵband(tg.ais_lo[fyi],    tg.ais_hi[fyi]),
    gsic   = ϵband(tg.gsic_lo[fyi],   tg.gsic_hi[fyi]),
    gis    = ϵband(tg.gis_lo[fyi],    tg.gis_hi[fyi]),
    steric = ϵband(tg.steric_lo[fyi], tg.steric_hi[fyi]),
    # Total obs error = Dangendorf σ ⊕ LWS-budget σ (propagate LWS uncertainty)
    dang   = sqrt.(Float64.(tg.dang_sig[fyi]).^2 .+ ϵband(tg.lws_lo[fyi], tg.lws_hi[fyi]).^2),
)

# Modern-rate point constraints (Gaussian)
const IMBIE_MU, IMBIE_SIG = 0.72, 0.156    # AIS 1992–2017 cm (Bamber et al. 2018 Nature 558:219; TODO: verify σ against IMBIE 2023 Otosaka et al. ESSD which reports 2720±1390 Gt → σ≈0.38cm at 1σ)
const DYU_MU,   DYU_SIG   = 2.127, 0.148   # GSIC 1961–2003 cm

# ---- free physical parameters (name, component, symbol, prior) ----
G = :glaciers_small_icecaps
FREE = [
    # DAIS key dynamics (weak shape params fixed at prior medoid)
    (name="ais_ocean_temperature₀", comp=:antarctic_icesheet, sym=:ais_ocean_temperature₀, μ=0.72, σ=0.50, lo=0.50, hi=2.00),
    (name="antarctic_alpha",        comp=:antarctic_icesheet, sym=:ais_α,                  μ=0.23, σ=0.15, lo=0.05, hi=0.60),
    (name="antarctic_nu",           comp=:antarctic_icesheet, sym=:ais_ν,                  μ=0.009,σ=0.005,lo=0.001,hi=0.025),
    (name="antarctic_temp_threshold",comp=:antarctic_icesheet,sym=:temperature_threshold,  μ=-15., σ=5.,   lo=-25., hi=-5.),
    (name="anto_alpha",             comp=:antarctic_ocean,    sym=:anto_α,                 μ=0.28, σ=0.10, lo=0.05, hi=0.70),
    (name="anto_beta",              comp=:antarctic_ocean,    sym=:anto_β,                 μ=0.95, σ=0.10, lo=0.50, hi=1.50),
    # Greenland
    (name="greenland_a",            comp=:greenland_icesheet, sym=:greenland_a,            μ=-1.37,σ=0.5,  lo=-3.0, hi=0.0),
    (name="greenland_b",            comp=:greenland_icesheet, sym=:greenland_b,            μ=8.06, σ=2.0,  lo=2.0,  hi=15.0),
    (name="greenland_alpha",        comp=:greenland_icesheet, sym=:greenland_α,            μ=8e-4, σ=4e-4, lo=1e-4, hi=2e-3),
    (name="greenland_beta",         comp=:greenland_icesheet, sym=:greenland_β,            μ=9e-5, σ=4e-5, lo=1e-5, hi=2e-4),
    (name="greenland_v0",           comp=:greenland_icesheet, sym=:greenland_v₀,           μ=7.52, σ=0.50, lo=6.0,  hi=9.0),
    # Thermal expansion
    (name="thermal_alpha",          comp=:thermal_expansion,  sym=:te_α,                  μ=0.16, σ=0.05, lo=0.05, hi=0.35),
    # Mengel glacier (2-timescale)
    (name="gic_a",                  comp=G, sym=:gic_a,        μ=0.45, σ=0.08, lo=0.32,  hi=0.55),
    (name="gic_b",                  comp=G, sym=:gic_b,        μ=0.52, σ=0.25, lo=0.25,  hi=1.00),
    (name="gic_T_lia",              comp=G, sym=:gic_T_lia,    μ=-0.45,σ=0.30, lo=-1.00, hi=-0.10),
    (name="gic_f",                  comp=G, sym=:gic_f,        μ=0.50, σ=0.30, lo=0.02,  hi=0.98),
    (name="gic_tau_fast",           comp=G, sym=:gic_tau_fast, μ=40.,  σ=30.,  lo=5.,    hi=80.),
    (name="gic_tau_slow",           comp=G, sym=:gic_tau_slow, μ=300., σ=200., lo=80.,   hi=800.),
]

const NP     = length(FREE)
const SERIES = [:ais, :gsic, :gis, :steric, :dang]
const NN     = 2 * length(SERIES)   # (sd, rho) per series
const NK     = NP + NN
println("MCMC: $NP physical + $NN AR(1)-noise = $NK free params")

# ---- build model once; fix medoid params; forcing fixed throughout ----
m = get_model(glacier_model=:mengel, lws=:central, start_year=Y0, end_year=Y1)
# Glacier params at the prior means (medoid); recalibration updates these via FREE
_apply_mengel_defaults!(m)
set_external_forcing!(m, gmst, ohc)

reref(v) = 100 .* (v .- sum(v[ib]) / length(ib))

function logposterior(θ)
    # Hard bounds on physical and AR(1) noise params
    @inbounds for k in 1:NP
        (θ[k] < FREE[k].lo || θ[k] > FREE[k].hi) && return -Inf
    end
    σn = θ[NP+1:2:NK]
    ρn = θ[NP+2:2:NK]
    (any(σn .<= 0) || any(ρn .< 0) || any(ρn .>= 0.99)) && return -Inf

    # Update physical params and run
    @inbounds for k in 1:NP
        update_param!(m, FREE[k].comp, FREE[k].sym, θ[k])
    end
    run(m)

    ais   = reref(m[:antarctic_icesheet,    :ais_sea_level])[myi]
    gsic  = reref(m[:glaciers_small_icecaps,:gsic_sea_level])[myi]
    gis   = reref(m[:greenland_icesheet,    :greenland_sea_level])[myi]
    te    = reref(m[:thermal_expansion,     :te_sea_level])[myi]
    total = ais .+ gsic .+ gis .+ te .+ obs.lws   # LWS from obs budget

    ll = 0.0
    for (i, (mod, ob, ev)) in enumerate([
            (ais,  obs.ais,    ϵ.ais),
            (gsic, obs.gsic,   ϵ.gsic),
            (gis,  obs.gis,    ϵ.gis),
            (te,   obs.steric, ϵ.steric),
            (total,obs.dang,   ϵ.dang)])
        ll += hetero_logl_ar1(mod .- ob, σn[i], ρn[i], ev)
    end

    # Modern-rate point constraints
    full_ais  = reref(m[:antarctic_icesheet,    :ais_sea_level])
    full_gsic = reref(m[:glaciers_small_icecaps,:gsic_sea_level])
    ll += logpdf(Normal(IMBIE_MU, IMBIE_SIG), full_ais[idx(2017)]  - full_ais[idx(1992)])
    ll += logpdf(Normal(DYU_MU,   DYU_SIG),   full_gsic[idx(2003)] - full_gsic[idx(1961)])

    # Gaussian priors on physical params; weak half-normal on noise σ
    lp = 0.0
    @inbounds for k in 1:NP
        lp += logpdf(Normal(FREE[k].μ, FREE[k].σ), θ[k])
    end
    for i in 1:length(SERIES)
        lp += logpdf(truncated(Normal(0, 5), 0, Inf), σn[i])
    end

    return ll + lp
end

# ---- starting point: MAP (if available) else prior means ----
map_path = joinpath(REPO, "outputs", "calib_full_joint_params.csv")
θ0 = if isfile(map_path)
    println("Starting from MAP: $map_path")
    mapp = CSV.read(map_path, DataFrame)
    name_to_map = Dict(mapp.param[i] => Float64(mapp.MAP[i]) for i in 1:nrow(mapp))
    [get(name_to_map, k.name, Float64(k.μ)) for k in FREE]
else
    println("MAP file not found — starting from prior means (slower mixing).")
    [Float64(k.μ) for k in FREE]
end
append!(θ0, repeat([1.0, 0.5], length(SERIES)))   # initial σ=1cm, ρ=0.5

# Initial proposal: diagonal scaled to prior σ; adapted cov seeds the next run
prop  = vcat([0.1 * Float64(k.σ) for k in FREE], repeat([0.3, 0.1], length(SERIES)))
adcov = joinpath(OUT_DIR, "adapted_cov.csv")
cov0  = isfile(adcov) ? Matrix(CSV.read(adcov, DataFrame)) : Matrix(Diagonal(prop.^2))
isfile(adcov) && println("Seeding proposal from adapted covariance ($adcov).")

@printf("logpost(θ0) = %.2f  (start = %s)\n", logposterior(θ0), isfile(map_path) ? "MAP" : "prior means")

Random.seed!(SEED)
mkpath(OUT_DIR)
@time chain, accept, covout, lp = RAM_sample(logposterior, θ0, cov0, N_ITER; opt_α=0.234, output_log_probability_x=true)

CSV.write(joinpath(OUT_DIR, "adapted_cov_seed$(SEED).csv"), DataFrame(covout, :auto))
@printf("Finished: %d iter, acceptance = %.3f\n", N_ITER, accept)
@printf("Target acceptance: 0.10–0.40 (RAM targets 0.234).\n")

pnames = vcat([k.name for k in FREE], vcat([["sd_$s", "rho_$s"] for s in SERIES]...))
burn   = chain[(N_ITER÷2+1):end, :]
println("\nPosterior (2nd-half) medians for key params:")
for nm in ["ais_ocean_temperature₀", "anto_alpha", "gic_T_lia", "gic_f", "gic_tau_fast", "gic_tau_slow", "gic_a"]
    c = burn[:, findfirst(==(nm), pnames)]
    @printf("  %-26s %.3g ± %.2g\n", nm, median(c), std(c))
end

out = joinpath(OUT_DIR, "chain_seed$(SEED)_n$(N_ITER).csv")
df  = DataFrame(chain, pnames)
df.log_post    = lp
df.accept_rate = fill(accept, nrow(df))
CSV.write(out, df)
println("\nWrote $out")
println("Next: run postprocess_mcmc_mengel.jl to diagnose convergence and build the posterior subsample.")
println("Production: ≥4 chains × ≥500k iterations on Torch `cs` partition (see docs/calibration_guide_mengel.md).")

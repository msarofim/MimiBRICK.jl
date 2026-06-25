## ============================================================================
## postprocess_mcmc_mengel.jl  —  combine chains, diagnostics, posterior subsample
##
## Reads all outputs/mcmc/chain_seed*_n*.csv produced by calibrate_mcmc_mengel.jl,
## burns the first half of each chain, computes Gelman–Rubin R̂ and ESS per
## parameter (MCMCDiagnosticTools), and writes a thinned posterior subsample to
## data/MimiBRICK/parameters_subsample_brick_mengel.csv.
##
## Also writes outputs/mcmc/adapted_cov.csv — the empirical posterior covariance —
## which seeds the next calibration run's proposal for faster mixing.
##
## Convergence targets: R̂ < 1.05 and ESS > 400 for all parameters.
## Params that fail → run longer chains (re-run calibrate_mcmc_mengel.jl with
## the seeded adapted_cov.csv and a larger n_iter).
##
## Usage:
##   julia --project=. calibration/postprocess_mcmc_mengel.jl [n_subsample]
##   n_subsample : number of draws in the output subsample (default 10000)
## ============================================================================

using CSV, DataFrames, Statistics, Printf, LinearAlgebra
using MCMCDiagnosticTools

const REPO    = abspath(joinpath(@__DIR__, ".."))
const MCMCDIR = joinpath(REPO, "outputs", "mcmc")
const OUT_CSV = joinpath(REPO, "data", "MimiBRICK", "parameters_subsample_brick_mengel.csv")

N_SUB = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10000

files = [joinpath(MCMCDIR, f) for f in readdir(MCMCDIR)
         if startswith(f, "chain_seed") && endswith(f, ".csv")]
isempty(files) && error("No chain files found in $MCMCDIR.\n" *
    "Run calibrate_mcmc_mengel.jl first.")

println("Combining $(length(files)) chain(s):")
chains = DataFrame[]
for f in files
    d    = CSV.read(f, DataFrame)
    burn = d[(nrow(d) ÷ 2 + 1):end, :]   # discard first half (burn-in)
    push!(chains, burn)
    @printf("  %-50s  %d post-burn rows, acceptance %.3f\n",
            basename(f), nrow(burn), d.accept_rate[1])
end

# Parameter names: all columns except log_post and accept_rate
pnames = [n for n in names(chains[1]) if !(n in ["log_post", "accept_rate"])]
nmin   = minimum(nrow.(chains))
nc     = length(chains)
println("\n$nc chain(s) × $nmin draws. Convergence (target R̂ < 1.05, ESS > 400):")

bad = String[]
for p in pnames
    arr = Array{Float64}(undef, nmin, nc)
    for (ci, ch) in enumerate(chains)
        arr[:, ci] = Float64.(ch[1:nmin, p])
    end
    r = rhat(arr)
    e = ess(arr)
    if r > 1.05 || e < 400
        push!(bad, p)
        @printf("  %-28s R̂ = %.3f  ESS = %.0f  <-- needs longer chains\n", p, r, e)
    end
end
isempty(bad) ? println("  All params converged (R̂ < 1.05, ESS > 400).") :
               @printf("  %d param(s) not converged → run longer chains.\n", length(bad))

# Thinned posterior subsample
pool = vcat(chains...)
n    = nrow(pool)
step = max(1, n ÷ N_SUB)
sub  = pool[1:step:end, :][1:min(N_SUB, nrow(pool[1:step:end, :])), :]
mkpath(dirname(OUT_CSV))
CSV.write(OUT_CSV, sub[:, pnames])
@printf("\nWrote %s  (%d-member subsample of %d pooled draws)\n", OUT_CSV, nrow(sub), n)

# Empirical posterior covariance → seeds next run's proposal
M    = Matrix{Float64}(pool[:, pnames])
Σ    = cov(M) .+ 1e-10 * I(size(M, 2))   # small ridge for numerical stability
cov_path = joinpath(MCMCDIR, "adapted_cov.csv")
CSV.write(cov_path, DataFrame(Σ, :auto))
@printf("Wrote %s  (%d-D empirical posterior covariance)\n", cov_path, size(M, 2))
println("Pass this to the next calibrate_mcmc_mengel.jl run for faster mixing.")

if !isempty(bad)
    println("\n** NOT CONVERGED ** Re-run with longer chains, seeded by adapted_cov.csv:")
    println("   julia --project=. calibration/calibrate_mcmc_mengel.jl 500000 <seed>")
end

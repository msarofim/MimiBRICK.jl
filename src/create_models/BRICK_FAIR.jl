using CSVFiles, DataFrames

# ============================================================================
# BRICK_FAIR.jl
#
# Convenience constructor that builds a BRICK model driven by pre-computed
# FaIR v2.2 GMST + OHC trajectories shipped with the package. This is the
# recommended entry point for running BRICK with an external climate model:
#
#   m = create_brick_fair(ssprcp_scenario="ssp245", glacier_model=:mengel)
#   run(m)   # uses FaIR SSP2-4.5 GMST + OHC
#
# For custom forcing (e.g. per-draw FaIR ensembles), build with get_model and
# call set_external_forcing!(m, gmst, ohc) directly:
#
#   m = get_model(glacier_model=:mengel, lws=:central)
#   set_external_forcing!(m, my_gmst, my_ohc)
#   run(m)
#
# FaIR forcing file format (annual, must cover start_year:end_year):
#   GMST: CSV with columns `year` (Int), `gmst_C` (°C anomaly rel 1850-1900)
#   OHC:  CSV with columns `year` (Int), `ohc_1e22J` (cumulative stock, 1e22 J)
# ============================================================================

"""
    create_brick_fair(; ssprcp_scenario, start_year, end_year, glacier_model, lws)

Return a MimiBRICK model instance driven by FaIR v2.2 GMST + OHC forcing.

Wraps `get_model` + `set_external_forcing!`, loading the forcing CSVs from the
package data directory. Default glacier model is `:mengel` and LWS treatment is
`:central` (deterministic 0.3 mm/yr) — appropriate for reproducible ensemble runs.

Function Arguments:

      ssprcp_scenario = SSP scenario tag; one of ssp119, ssp126, ssp245 (default),
                        ssp370, ssp460, ssp585. Resolves to the corresponding
                        data/observations/fair_mean_{gmst,ohc}_<scenario>.csv.
      start_year      = initial year of the simulation period (default 1850)
      end_year        = ending year of the simulation period (default 2100)
      glacier_model   = :mengel (default) or :gsic
      lws             = :central (default), :zero, or :random
                        (see get_model docstring)
"""
function create_brick_fair(; ssprcp_scenario::String="ssp245", start_year::Int=1850, end_year::Int=2100,
                             glacier_model::Symbol=:mengel, lws::Symbol=:central)

    obs_dir   = joinpath(@__DIR__, "..", "..", "data", "observations")
    gmst_path = joinpath(obs_dir, "fair_mean_gmst_$(ssprcp_scenario).csv")
    ohc_path  = joinpath(obs_dir, "fair_mean_ohc_$(ssprcp_scenario).csv")

    isfile(gmst_path) || error("create_brick_fair: no FaIR GMST file for scenario \"$ssprcp_scenario\" " *
                               "(expected $gmst_path)")
    isfile(ohc_path)  || error("create_brick_fair: no FaIR OHC file for scenario \"$ssprcp_scenario\" " *
                               "(expected $ohc_path)")

    years = collect(start_year:end_year)

    # Load and index by year — robust to extra rows in the CSV.
    function load_by_year(path, vcol)
        df  = DataFrame(load(path))
        lut = Dict(Int(row.year) => Float64(row[vcol]) for row in eachrow(df))
        missing_years = filter(y -> !haskey(lut, y), years)
        isempty(missing_years) || error("create_brick_fair: $path missing years $(first(missing_years))…$(last(missing_years))")
        return [lut[y] for y in years]
    end

    gmst = load_by_year(gmst_path, :gmst_C)
    ohc  = load_by_year(ohc_path,  :ohc_1e22J)

    m = get_model(ssprcp_scenario=ssprcp_scenario, start_year=start_year, end_year=end_year,
                  glacier_model=glacier_model, lws=lws)
    set_external_forcing!(m, gmst, ohc)
    return m
end

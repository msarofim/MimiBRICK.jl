module MimiBRICK

# Load required packages
using CSVFiles, DataFrames, Distributions, Mimi

# Load BRICK Mimi components.
include(joinpath("components", "antarctic_icesheet_component.jl"))
include(joinpath("components", "antarctic_ocean_component.jl"))
include(joinpath("components", "glaciers_small_icecaps_component.jl"))
include(joinpath("components", "glaciers_mengel_component.jl"))
include(joinpath("components", "global_sea_level_component.jl"))
include(joinpath("components", "greenland_icesheet_component.jl"))
include(joinpath("components", "landwater_storage_component.jl"))
include(joinpath("components", "thermal_expansion_component.jl"))

# include creation functions for two other model variants
include(joinpath("create_models", "SNEASY_BRICK.jl"))
include(joinpath("create_models", "BRICK_DOECLIM.jl"))
include(joinpath("create_models", "BRICK_FAIR.jl"))

# include other helper functions
include(joinpath("downscale.jl"))
include(joinpath("brick_param_updates.jl"))

# include calibration functions
include(joinpath("calibration/main_calibration.jl"))

# function to create 'Building blocks for Relevant Ice and Climate Knowledge' (BRICK) model.
"""
    get_model(; ssprcp_scenario, start_year, end_year, glacier_model, lws)

Return a MimiBRICK model instance driven by pre-computed SNEASY forcing.

Function Arguments:

      ssprcp_scenario = (SSP-)RCP scenario for exogenous forcing;
                        valid values: RCP26, RCP45, RCP60, RCP85, ssp119, ssp126,
                        ssp245, ssp370, ssp460, ssp585, ssp534-over
      start_year      = initial year of the simulation period
      end_year        = ending year of the simulation period
      glacier_model   = :gsic (default, Wigley-Raper-Bakker) or :mengel
                        (Mengel-2016 two-timescale emulator; see PR1)
      lws             = land-water-storage treatment:
                          :random  (default) = unseeded draw from N(0.0003, 0.00018) m/yr
                          :central           = deterministic 0.3 mm/yr mean (recommended
                                               for reproducible ensemble runs)
                          :zero              = no LWS contribution
"""
function get_model(;ssprcp_scenario::String="ssp245", start_year::Int=1850, end_year::Int=2020,
                    glacier_model::Symbol=:gsic, lws::Symbol=:random)

    glacier_model in (:gsic, :mengel) || error("get_model: glacier_model must be :gsic or :mengel (got :$glacier_model)")
    lws in (:random, :central, :zero) || error("get_model: lws must be :random, :central, or :zero (got :$lws)")

    #-----------------------#
    # ----- Load Data ----- #
    #-----------------------#

    # Load exogenous time-series for global surface temperature and ocean heat content (output from SNEASY under RCP8.5).
    temperature_scenario = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "sneasy_temperature_"*ssprcp_scenario*"_1850_2300_07-06-2026.csv")))
    oceanheat_scenario   = DataFrame(load(joinpath(@__DIR__, "..", "data", "model_data", "sneasy_oceanheat_"*ssprcp_scenario*"_1850_2300_07-06-2026.csv")))

    #-------------------------#
    # ----- Build BRICK ----- #
    #-------------------------#

    # Create a Mimi-model instance.
    brick = Model()

    # Set model time horizon, defaults to 1850-2300.
    set_dimension!(brick, :time, start_year:end_year)

    # Add in BRICK components.
    add_comp!(brick, antarctic_ocean)
    add_comp!(brick, antarctic_icesheet)
    add_comp!(brick, glaciers_small_icecaps)
    add_comp!(brick, greenland_icesheet)
    add_comp!(brick, thermal_expansion)
    add_comp!(brick, landwater_storage)
    add_comp!(brick, global_sea_level)

    #-------------------------------------#
    # ----- Assign Model Parameters ----- #
    #-------------------------------------#

    # ----- Antarctic Ocean ----- #

    update_param!(brick, :antarctic_ocean, :anto_α, 0.28)
    update_param!(brick, :antarctic_ocean, :anto_β, 0.95)

    # ----- Antarctic Ice Sheet ----- #

    update_param!(brick, :antarctic_icesheet, :ais_ρ_ice, 917.0)
    update_param!(brick, :antarctic_icesheet, :ais_ρ_seawater, 1030.0)
    update_param!(brick, :antarctic_icesheet, :ais_ρ_rock, 4000.0)
    update_param!(brick, :antarctic_icesheet, :ais_sea_level₀, 0.0)
    update_param!(brick, :antarctic_icesheet, :ais_ocean_temperature₀, 0.72)
    update_param!(brick, :antarctic_icesheet, :ais_radius₀, 1.864e6)
    update_param!(brick, :antarctic_icesheet, :ais_bedheight₀, 781.0)
    update_param!(brick, :antarctic_icesheet, :ais_slope, 0.0006)
    update_param!(brick, :antarctic_icesheet, :ais_μ, 11.0)
    update_param!(brick, :antarctic_icesheet, :ais_runoffline_snowheight₀, 1400.0)
    update_param!(brick, :antarctic_icesheet, :ais_c, 100.0)
    update_param!(brick, :antarctic_icesheet, :ais_precipitation₀, log(0.37))
    update_param!(brick, :antarctic_icesheet, :ais_κ, 0.062)
    update_param!(brick, :antarctic_icesheet, :ais_ν, 0.0086)
    update_param!(brick, :antarctic_icesheet, :ais_iceflow₀, 1.2)
    update_param!(brick, :antarctic_icesheet, :ais_γ, 2.9)
    update_param!(brick, :antarctic_icesheet, :ais_α, 0.23)
    update_param!(brick, :antarctic_icesheet, :ais_temperature_coefficient, 0.8365)
    update_param!(brick, :antarctic_icesheet, :ais_temperature_intercept, 15.42)
    update_param!(brick, :antarctic_icesheet, :ais_local_fingerprint, -1.18)
    update_param!(brick, :antarctic_icesheet, :ocean_surface_area, 3.619e14)
    update_param!(brick, :antarctic_icesheet, :temperature_threshold, -15.0)
    update_param!(brick, :antarctic_icesheet, :λ, 0.0093)
    update_param!(brick, :antarctic_icesheet, :include_ais_DSL, true)

    # ----- Glaciers & Small Ice Caps ----- #

    update_param!(brick, :glaciers_small_icecaps, :gsic_β₀, 0.0013)
    update_param!(brick, :glaciers_small_icecaps, :gsic_v₀, 0.376)
    update_param!(brick, :glaciers_small_icecaps, :gsic_s₀, -0.0138)
    update_param!(brick, :glaciers_small_icecaps, :gsic_n, 0.847)
    update_param!(brick, :glaciers_small_icecaps, :gsic_teq, -0.15)

    # ----- Greenland Ice Sheet ----- #

    update_param!(brick, :greenland_icesheet, :greenland_a, -1.37)
    update_param!(brick, :greenland_icesheet, :greenland_b, 8.06)
    update_param!(brick, :greenland_icesheet, :greenland_α, 0.0008)
    update_param!(brick, :greenland_icesheet, :greenland_β, 0.00009)
    update_param!(brick, :greenland_icesheet, :greenland_v₀, 7.52)

    # ----- Thermal Expansion ----- #

    update_param!(brick, :thermal_expansion, :te_A, 3.619e14)
    update_param!(brick, :thermal_expansion, :te_C, 3991.86795711963)
    update_param!(brick, :thermal_expansion, :te_ρ, 1027.0)
    update_param!(brick, :thermal_expansion, :te_α, 0.16)
    update_param!(brick, :thermal_expansion, :te_s₀, 0.0)
    update_param!(brick, :thermal_expansion, :ocean_heat_mixed, zeros(length(start_year:end_year)))
    oceanheat_idx = findall((in)(start_year:end_year), oceanheat_scenario[!,:Year])
    update_param!(brick, :thermal_expansion, :ocean_heat_interior, oceanheat_scenario[oceanheat_idx, :"MAP Ocean Heat"])

    # ----- Landwater Storage ----- #

    update_param!(brick, :landwater_storage, :lws₀, 0.0)
    update_param!(brick, :landwater_storage, :first_projection_year, 2018)
    n = length(start_year:end_year)
    lws_rates = lws === :central ? fill(0.0003, n) :
                lws === :zero    ? zeros(n) :
                rand(Normal(0.0003, 0.00018), n)   # :random (default)
    update_param!(brick, :landwater_storage, :lws_random_sample, lws_rates)

    # ----- Set Parameters With Common Values Across Components ----- #

    add_shared_param!(brick, :model_seawater_freeze, -1.8)
    connect_param!(brick, :antarctic_icesheet, :seawater_freeze, :model_seawater_freeze)
    connect_param!(brick, :antarctic_ocean, :seawater_freeze, :model_seawater_freeze)

    temperature_idx = findall((in)(start_year:end_year), temperature_scenario[!,:Year])

    add_shared_param!(brick, :model_global_surface_temperature,  temperature_scenario[temperature_idx,:"MAP Temperature"], dims = [:time])
    connect_param!(brick, :antarctic_icesheet, :global_surface_temperature, :model_global_surface_temperature)
    connect_param!(brick, :antarctic_ocean, :global_surface_temperature, :model_global_surface_temperature)
    connect_param!(brick, :glaciers_small_icecaps, :global_surface_temperature, :model_global_surface_temperature)
    connect_param!(brick, :greenland_icesheet, :global_surface_temperature, :model_global_surface_temperature)

    #-----------------------------------------#
    #----- Create Component Connections ----- #
    #-----------------------------------------#

    # Set parameter connections (:component => :parameter).
    connect_param!(brick, :global_sea_level => :slr_glaciers_small_ice_caps, :glaciers_small_icecaps => :gsic_sea_level)
    connect_param!(brick, :global_sea_level => :slr_greeland_icesheet,       :greenland_icesheet     => :greenland_sea_level)
    connect_param!(brick, :global_sea_level => :slr_antartic_icesheet,       :antarctic_icesheet     => :ais_sea_level)
    connect_param!(brick, :global_sea_level => :slr_thermal_expansion,       :thermal_expansion      => :te_sea_level)
    connect_param!(brick, :global_sea_level => :slr_landwater_storage,       :landwater_storage      => :lws_sea_level)

    connect_param!(brick, :antarctic_icesheet => :antarctic_ocean_temperature, :antarctic_ocean  => :anto_temperature)
    connect_param!(brick, :antarctic_icesheet => :global_sea_level,            :global_sea_level => :sea_level_rise)

    # ----- Optional: swap in the Mengel-2016 glacier emulator -----
    # See create_brick_doeclim for details. Both glacier models cover the SAME
    # inventory (glaciers AND small ice caps); Mimi.replace! preserves the
    # :glaciers_small_icecaps slot name and the name-matched I/O wiring. The default
    # :gsic path is left untouched.
    if glacier_model == :mengel
        Mimi.replace!(brick, :glaciers_small_icecaps => glaciers_mengel)
        _apply_mengel_defaults!(brick)
    end

    # Return BRICK model.
    return brick
end

"""
    set_external_forcing!(m, gmst, ohc)

Override the exogenous GMST and ocean interior heat content on a `get_model`-built
BRICK instance with external (e.g. FaIR) forcing trajectories.

Arguments:
  gmst : global mean surface temperature anomaly relative to 1850–1900 (°C),
         length must equal the model time dimension
  ohc  : cumulative ocean heat content stock (1×10²² J), same length

Note: only valid for `get_model`-built models, which expose a shared
`:model_global_surface_temperature` parameter. Not applicable to
`create_brick_doeclim` or `create_sneasy_brick` models, where temperature
is computed internally by the coupled climate component.
"""
function set_external_forcing!(m, gmst::Vector{<:Real}, ohc::Vector{<:Real})
    update_param!(m, :model_global_surface_temperature, gmst)
    update_param!(m, :thermal_expansion, :ocean_heat_interior, ohc)
    return m
end

end # Module

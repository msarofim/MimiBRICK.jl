using Mimi

# ============================================================================
# brick_param_updates.jl
#
# Single source of truth for mapping a BRICK posterior CSV row onto a
# MimiBRICK model instance. Called by the FaIR-driven projection loop and
# the MCMC calibration; centralising here prevents the mapping from drifting
# across callers if the posterior column set ever changes.
#
# Greek-letter Mimi parameter names (anto_α, ais_μ, te_α, …) use Unicode to
# match the physics notation in the BRICK paper. The posterior CSV uses ASCII
# transliterations (anto_alpha, antarctic_mu, thermal_alpha, …) because CSV
# column names must be plain ASCII. The mapping below is the authoritative
# correspondence.
# ============================================================================

"""
    update_brick_params!(m, prow; precip_log=false, skip_glaciers=false)

Apply one posterior CSV row `prow` to a built MimiBRICK model `m` in place.

`precip_log` — only needed when applying an **old v1.x posterior** to a v2.0.0
model. MimiBRICK v2.0.0 stores `ais_precipitation₀` in log-space and computes
`exp(ais_precipitation₀) * exp(κ·T)`. The pre-v2.0.0 posterior stores
`antarctic_precip0` in linear m/yr (~0.80). Passing that directly into v2.0.0's
`exp()` overloads Antarctic snowfall by 2–4×, causing the +100 cm AIS@1900 blowup.
Set `precip_log=true` to apply `log()` before setting the parameter, restoring
v1.x behaviour (`exp(log(p)) = p`) without touching the model equations.
The BRICK-Mengel posterior shipped with this package was calibrated on v2.0.0 and
stores `antarctic_precip0` already in log-space — leave `precip_log=false`.

`skip_glaciers` — set `true` when the model's glacier slot has been replaced by
the Mengel emulator. The default WRB params (gsic_β₀, gsic_v₀, gsic_s₀, gsic_n)
do not exist on the Mengel component; set the Mengel params separately via
`update_param!` or the projection loop's PHYS mapping.
"""
function update_brick_params!(m, prow; precip_log::Bool=false, skip_glaciers::Bool=false)

    # Antarctic Ocean
    update_param!(m, :antarctic_ocean, :anto_α, prow.anto_alpha)
    update_param!(m, :antarctic_ocean, :anto_β, prow.anto_beta)

    # Antarctic Ice Sheet (15 parameters; ais_precipitation₀ log-shim if needed)
    update_param!(m, :antarctic_icesheet, :ais_sea_level₀,             prow.antarctic_s0)
    update_param!(m, :antarctic_icesheet, :ais_bedheight₀,             prow.antarctic_bed_height0)
    update_param!(m, :antarctic_icesheet, :ais_slope,                  prow.antarctic_slope)
    update_param!(m, :antarctic_icesheet, :ais_μ,                      prow.antarctic_mu)
    update_param!(m, :antarctic_icesheet, :ais_runoffline_snowheight₀, prow.antarctic_runoff_height0)
    update_param!(m, :antarctic_icesheet, :ais_c,                      prow.antarctic_c)
    update_param!(m, :antarctic_icesheet, :ais_precipitation₀,
                  precip_log ? log(prow.antarctic_precip0) : prow.antarctic_precip0)
    update_param!(m, :antarctic_icesheet, :ais_κ,                      prow.antarctic_kappa)
    update_param!(m, :antarctic_icesheet, :ais_ν,                      prow.antarctic_nu)
    update_param!(m, :antarctic_icesheet, :ais_iceflow₀,               prow.antarctic_flow0)
    update_param!(m, :antarctic_icesheet, :ais_γ,                      prow.antarctic_gamma)
    update_param!(m, :antarctic_icesheet, :ais_α,                      prow.antarctic_alpha)
    update_param!(m, :antarctic_icesheet, :temperature_threshold,      prow.antarctic_temp_threshold)
    update_param!(m, :antarctic_icesheet, :λ,                          prow.antarctic_lambda)

    # Glaciers & small ice caps (WRB single-reservoir; skip for Mengel-swapped models)
    if !skip_glaciers
        update_param!(m, :glaciers_small_icecaps, :gsic_β₀, prow.glaciers_beta0)
        update_param!(m, :glaciers_small_icecaps, :gsic_v₀, prow.glaciers_v0)
        update_param!(m, :glaciers_small_icecaps, :gsic_s₀, prow.glaciers_s0)
        update_param!(m, :glaciers_small_icecaps, :gsic_n,  prow.glaciers_n)
    end

    # Greenland Ice Sheet
    update_param!(m, :greenland_icesheet, :greenland_a,  prow.greenland_a)
    update_param!(m, :greenland_icesheet, :greenland_b,  prow.greenland_b)
    update_param!(m, :greenland_icesheet, :greenland_α,  prow.greenland_alpha)
    update_param!(m, :greenland_icesheet, :greenland_β,  prow.greenland_beta)
    update_param!(m, :greenland_icesheet, :greenland_v₀, prow.greenland_v0)

    # Thermal expansion
    update_param!(m, :thermal_expansion, :te_α,  prow.thermal_alpha)
    update_param!(m, :thermal_expansion, :te_s₀, prow.thermal_s0)
end

#!/usr/bin/env julia
#
# earth_model_comparison.jl
#
# Compare the solid-earth response at the end of two otherwise-identical
# lgm_ice_nh simulations that differ only in the geodynamic earth model:
#
#   igeo2  ->  VILMA1  (i_geo = 2)
#   igeo3  ->  VILMA2  (i_geo = 3)
#
# The comparison uses geo_restart.nc from each run, which is written by the
# CLIMBER-X geo module on a common 2880x1440 (0.125 deg) grid, so the fields
# are directly comparable between the two earth models. VILMA1's own rsl.nc is
# on a different grid and has no VILMA2 counterpart, so it is not used.
#
# Fields plotted (each as three panels: VILMA1, VILMA2, difference):
#   z_bed              bedrock elevation                 [m]
#   rsl                relative sea level                [m]
#   z_bed - z_bed_ref  bedrock deflection (GIA response) [m]
#
# Figures are written to analysis/figures/<run>/ (gitignored).
#
# Run:  julia --project=analysis analysis/earth_model_comparison.jl [run]
#       where [run] is the output/ subdirectory holding igeo2/ and igeo3/
#       (default: lgm_ice_nh).

using NCDatasets
using CairoMakie
using Printf
using Statistics

# --- Configuration ----------------------------------------------------------

const RUN = isempty(ARGS) ? "lgm_ice_nh" : ARGS[1]

const OUTDIR = joinpath(@__DIR__, "..", "output", RUN)

const RUNS = (
    vilma1  = joinpath(OUTDIR, "igeo2", "restart_out", "year_500", "geo_restart.nc"),
    vilma2  = joinpath(OUTDIR, "igeo3", "restart_out", "year_500", "geo_restart.nc"),
)

const LABELS = (vilma1 = "VILMA1 (i_geo=2)", vilma2 = "VILMA2 (i_geo=3)")

const FIGDIR = joinpath(@__DIR__, "figures", RUN)

# --- Data loading -----------------------------------------------------------

"""Load one geo_restart.nc; return NamedTuple of lon, lat and 2D fields."""
function load_geo(path)
    isfile(path) || error("missing file: $path")
    NCDataset(path) do ds
        lon = Array(ds["lon"][:])
        lat = Array(ds["lat"][:])
        # NCDatasets presents (lat,lon) CDL vars as (lon,lat) Julia arrays.
        z_bed     = Array(ds["z_bed"][:, :])
        z_bed_ref = Array(ds["z_bed_ref"][:, :])
        rsl       = Array(ds["rsl"][:, :])
        sea_level = Array(ds["sea_level"][:])[1]
        @assert size(z_bed) == (length(lon), length(lat)) "unexpected field orientation"
        (; lon, lat, z_bed, z_bed_ref, rsl, sea_level)
    end
end

# --- Plotting helpers -------------------------------------------------------

"""Symmetric limit for a difference field, ignoring NaNs."""
function sym_limit(a)
    m = maximum(abs, filter(isfinite, a); init = 0.0)
    m == 0 ? 1.0 : m
end

"""Shared limits across two comparable fields."""
function shared_limits(a, b)
    v = filter(isfinite, vcat(vec(a), vec(b)))
    isempty(v) ? (-1.0, 1.0) : (minimum(v), maximum(v))
end

"""Overlay present-day continent outlines (z_bed_ref = 0) on an axis."""
add_coast!(ax, a) =
    contour!(ax, a.lon, a.lat, a.z_bed_ref; levels = [0.0],
             color = :black, linewidth = 0.5)

"""
Three-panel comparison map for one field, stacked vertically:
row 1 = VILMA1, row 2 = VILMA2 (sharing one colorbar), row 3 = difference.

`get` extracts the 2D array (lon,lat) from a run NamedTuple. `coast` overlays
present-day continent outlines. `bsl`, if given as (vilma1=, vilma2=), prints
the barystatic sea level on each absolute panel.
"""
function panel_field(a, b, get; title, unit, cmap = :vik, coast = false, bsl = nothing)
    da, db = get(a), get(b)
    ddiff = db .- da
    lims = shared_limits(da, db)
    dlim = sym_limit(ddiff)

    fig = Figure(size = (1000, 1400), fontsize = 15)
    Label(fig[0, 1:2], title; fontsize = 20, font = :bold)

    abs_panels = (
        (1, da, LABELS.vilma1, bsl === nothing ? nothing : bsl.vilma1),
        (2, db, LABELS.vilma2, bsl === nothing ? nothing : bsl.vilma2),
    )

    local hm_abs
    for (row, data, ttl, bslval) in abs_panels
        ax = Axis(fig[row, 1]; title = ttl, ylabel = "Latitude",
                  xlabel = "", aspect = DataAspect())
        hm_abs = heatmap!(ax, a.lon, a.lat, data; colorrange = lims, colormap = cmap)
        coast && add_coast!(ax, a)
        if bslval !== nothing
            text!(ax, -160, -60; text = @sprintf("bsl = %.1f m", bslval),
                  align = (:left, :center), fontsize = 15, color = :black)
        end
    end
    # One shared colorbar for the two absolute panels, half their combined height.
    Colorbar(fig[1:2, 2], hm_abs; label = "$title [$unit]", height = Relative(1 / 2))

    axd = Axis(fig[3, 1]; title = "Difference (VILMA2 − VILMA1)",
               xlabel = "Longitude", ylabel = "Latitude", aspect = DataAspect())
    hmd = heatmap!(axd, a.lon, a.lat, ddiff; colorrange = (-dlim, dlim), colormap = :RdBu)
    coast && add_coast!(axd, a)
    Colorbar(fig[3, 2], hmd; label = "$title [$unit]")

    rowgap!(fig.layout, 6)
    fig
end

# --- Main -------------------------------------------------------------------

function main()
    mkpath(FIGDIR)

    println("Loading runs:")
    for (k, p) in pairs(RUNS)
        println("  $(LABELS[k]): $p")
    end
    d = map(load_geo, RUNS)

    bsl = (vilma1 = d.vilma1.sea_level, vilma2 = d.vilma2.sea_level)
    println("\nBarystatic sea level (sea_level scalar in restart):")
    @printf("  %-22s % .3f m\n", LABELS.vilma1, bsl.vilma1)
    @printf("  %-22s % .3f m\n", LABELS.vilma2, bsl.vilma2)
    if bsl.vilma2 == 0
        println("  (VILMA2 sea_level == 0 => this output predates the geo.f90 fix that " *
                "computes the scalar for i_geo==3; the rsl field is valid regardless.)")
    end

    figs = (
        ("z_bed",  panel_field(d.vilma1, d.vilma2, x -> x.z_bed;
                               title = "Bedrock elevation", unit = "m", cmap = :bukavu)),
        ("rsl",    panel_field(d.vilma1, d.vilma2, x -> x.rsl;
                               title = "Relative sea level", unit = "m", cmap = :vik,
                               coast = true, bsl = bsl)),
        ("dz_bed", panel_field(d.vilma1, d.vilma2, x -> x.z_bed .- x.z_bed_ref;
                               title = "Bedrock deflection (z_bed − z_bed_ref)", unit = "m",
                               cmap = :vik, coast = true)),
    )

    println("\nField statistics (global), mean ± std over (VILMA2 − VILMA1):")
    for (name, _) in figs
        get = name == "z_bed"  ? (x -> x.z_bed) :
              name == "rsl"    ? (x -> x.rsl)   :
                                 (x -> x.z_bed .- x.z_bed_ref)
        diff = filter(isfinite, vec(get(d.vilma2) .- get(d.vilma1)))
        @printf("  %-8s  mean=% .3f  std=% .3f  max|Δ|=% .3f m\n",
                name, mean(diff), std(diff), maximum(abs, diff))
    end

    println("\nWriting figures to $(FIGDIR):")
    for (name, fig) in figs
        path = joinpath(FIGDIR, "compare_$(name).png")
        save(path, fig; px_per_unit = 2)
        println("  $path")
    end
    println("\nDone.")
end

main()

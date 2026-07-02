#!/usr/bin/env julia
#
# earth_model_comparison.jl
#
# Compare the solid-earth response at the end of two otherwise-identical
# lgm_ice_nh simulations that differ only in the geodynamic earth model:
#
#   igeo2  ->  VILMA        (i_geo = 2)
#   igeo3  ->  FastEarth3D  (i_geo = 3)
#
# The comparison uses geo_restart.nc from each run, which is written by the
# CLIMBER-X geo module on a common 2880x1440 (0.125 deg) grid, so the fields
# are directly comparable between the two earth models. VILMA's own rsl.nc is
# on a different grid and has no FastEarth3D counterpart, so it is not used.
#
# Fields plotted (each as three panels: VILMA, FastEarth3D, difference):
#   z_bed              bedrock elevation                 [m]
#   rsl                relative sea level                [m]
#   z_bed - z_bed_ref  bedrock deflection (GIA response) [m]
#
# Figures are written to analysis/figures/ (gitignored).
#
# Run:  julia --project=analysis analysis/earth_model_comparison.jl

using NCDatasets
using CairoMakie
using Printf
using Statistics

# --- Configuration ----------------------------------------------------------

const OUTDIR = joinpath(@__DIR__, "..", "output", "lgm_ice_nh")

const RUNS = (
    vilma      = joinpath(OUTDIR, "igeo2", "restart_out", "year_500", "geo_restart.nc"),
    fastearth  = joinpath(OUTDIR, "igeo3", "restart_out", "year_500", "geo_restart.nc"),
)

const LABELS = (vilma = "VILMA (i_geo=2)", fastearth = "FastEarth3D (i_geo=3)")

const FIGDIR = joinpath(@__DIR__, "figures")

# Northern-hemisphere crop (deg N) — where the LGM NH ice sheets sit.
const LAT_MIN = 20.0

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

"""Restrict fields to lat >= LAT_MIN, returning the cropped lat and views."""
function crop_nh(d)
    jj = findall(>=(LAT_MIN), d.lat)
    (; d.lon, lat = d.lat[jj],
       z_bed = d.z_bed[:, jj], z_bed_ref = d.z_bed_ref[:, jj], rsl = d.rsl[:, jj])
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

"""
Three-panel comparison map for one field.
`get` extracts the 2D array (lon,lat) from a cropped-run NamedTuple.
"""
function panel_field(a, b, get; title, unit, cmap = :vik)
    da, db = get(a), get(b)
    ddiff = db .- da
    lims = shared_limits(da, db)
    dlim = sym_limit(ddiff)

    fig = Figure(size = (1600, 460), fontsize = 15)
    Label(fig[0, 1:3], title; fontsize = 20, font = :bold)

    specs = (
        (1, da, LABELS.vilma,               lims, cmap),
        (2, db, LABELS.fastearth,           lims, cmap),
        (3, ddiff, "Difference (FE3D − VILMA)", (-dlim, dlim), :RdBu),
    )

    for (col, data, ttl, clim, cm) in specs
        ax = Axis(fig[1, col]; title = ttl, xlabel = "Longitude",
                  ylabel = col == 1 ? "Latitude" : "", aspect = DataAspect())
        hm = heatmap!(ax, a.lon, a.lat, data; colorrange = clim, colormap = cm)
        Colorbar(fig[2, col], hm; vertical = false, flipaxis = false,
                 label = "$title [$unit]")
    end
    rowgap!(fig.layout, 8)
    rowsize!(fig.layout, 2, Relative(0.12))
    fig
end

# --- Main -------------------------------------------------------------------

function main()
    mkpath(FIGDIR)

    println("Loading runs:")
    for (k, p) in pairs(RUNS)
        println("  $(LABELS[k]): $p")
    end
    raw = map(load_geo, RUNS)
    d = map(crop_nh, raw)

    println("\nGlobal-mean relative sea level (sea_level scalar in restart):")
    @printf("  %-22s % .3f m\n", LABELS.vilma,     raw.vilma.sea_level)
    @printf("  %-22s % .3f m  (FE3D leaves the scalar at 0; use the rsl field)\n",
            LABELS.fastearth, raw.fastearth.sea_level)

    figs = (
        ("z_bed",  panel_field(d.vilma, d.fastearth, x -> x.z_bed;
                               title = "Bedrock elevation", unit = "m", cmap = :bukavu)),
        ("rsl",    panel_field(d.vilma, d.fastearth, x -> x.rsl;
                               title = "Relative sea level", unit = "m", cmap = :vik)),
        ("dz_bed", panel_field(d.vilma, d.fastearth, x -> x.z_bed .- x.z_bed_ref;
                               title = "Bedrock deflection (z_bed − z_bed_ref)", unit = "m",
                               cmap = :vik)),
    )

    println("\nField statistics (NH, lat ≥ $(LAT_MIN)°), mean ± std over |FE3D − VILMA|:")
    for (name, _) in figs
        get = name == "z_bed"  ? (x -> x.z_bed) :
              name == "rsl"    ? (x -> x.rsl)   :
                                 (x -> x.z_bed .- x.z_bed_ref)
        diff = filter(isfinite, vec(get(d.fastearth) .- get(d.vilma)))
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

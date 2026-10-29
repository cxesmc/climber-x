#!/usr/bin/env julia
#
# lgm_vs_glac1d.jl
#
# Map of the simulated Northern Hemisphere ice sheets next to the GLAC1D
# (Tarasov) LGM reconstruction, on a common north polar stereographic grid:
#
#   fill   RSL relative to today [m] (> 0: land lower relative to sea than today)
#   lines  ice margin (thick) and ice surface elevation every 500 m
#   grey   present-day coastline
#
# Model fields come from geo_hires.nc of the run (0.125 deg): h_ice, z_sur and
# rsl_ref (relative sea level relative to the reference PI bedrock). GLAC1D
# topography is relative to contemporaneous sea level, so its RSL relative to
# today is -(z_bed_lgm - z_bed_0ka), from input/geo_ice_tarasov_lgm.nc and
# input/geo_ice_tarasov_lgc_0ka.nc.
#
# Figures are written to analysis/figures/<run>/ (gitignored).
#
# Run:  julia --project=analysis analysis/lgm_vs_glac1d.jl [run] [year]
#       where [run] is the output/ subdirectory (default: lgm_ice_nh) and
#       [year] the output time to plot (default: the last one).

using NCDatasets
using CairoMakie
using Printf

# --- Configuration ----------------------------------------------------------

const RUN = length(ARGS) >= 1 ? ARGS[1] : "lgm_ice_nh"
const YEAR = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : nothing

const GEO_HIRES = joinpath(@__DIR__, "..", "output", RUN, "geo_hires.nc")
const INPUT     = joinpath(@__DIR__, "..", "input")
const GLAC_LGM  = joinpath(INPUT, "geo_ice_tarasov_lgm.nc")
const GLAC_0KA  = joinpath(INPUT, "geo_ice_tarasov_lgc_0ka.nc")

const FIGDIR = joinpath(@__DIR__, "figures", RUN)

# Target grid: spherical polar stereographic [km], LON0 pointing down
const R_EARTH = 6371.0
const LON0    = -40.0
const XS      = range(-4800, 4800; step = 20.0)
const YS      = XS

const RSL_LEVELS = -600:50:600
const ZS_LEVELS  = 500:500:4000
const H_MIN      = 10.0     # [m] ice thickness taken as ice covered

# --- Regridding -------------------------------------------------------------

"""Latitude and longitude [deg] of a polar stereographic point (x, y) [km]."""
function stereo_to_ll(x, y)
    rho = hypot(x, y)
    lat = 90 - 2 * atand(rho / (2R_EARTH))
    lon = mod(LON0 + atand(x, -y) + 180, 360) - 180
    return lat, lon
end

"""Bilinear value of f(lon,lat) on a regular, lon-periodic grid; NaN outside lat."""
function sample(f, lon, lat, lo, la)
    fi = (la - lat[1]) / (lat[2] - lat[1]) + 1
    (fi < 1 || fi > length(lat)) && return NaN
    i0 = clamp(floor(Int, fi), 1, length(lat) - 1); wi = fi - i0
    fj = mod(lo - lon[1], 360) / (lon[2] - lon[1]) + 1
    j0 = floor(Int, fj); wj = fj - j0
    j0 = mod1(j0, length(lon)); j1 = mod1(j0 + 1, length(lon))
    (1 - wi) * ((1 - wj) * f[j0, i0]     + wj * f[j1, i0]) +
         wi  * ((1 - wj) * f[j0, i0 + 1] + wj * f[j1, i0 + 1])
end

"""Regrid a (lon,lat) field onto the stereographic target grid."""
function regrid(f, lon, lat)
    [sample(f, lon, lat, reverse(stereo_to_ll(x, y))...) for x in XS, y in YS]
end

field(v) = Float64.(nomissing(v, NaN))

# --- Data loading -----------------------------------------------------------

"""Model fields at the requested (or last) output time of geo_hires.nc."""
function load_model(path, year)
    isfile(path) || error("missing file: $path")
    NCDataset(path) do ds
        lon, lat = Float64.(ds["lon"][:]), Float64.(ds["lat"][:])
        time = ds["time"][:]
        n = isnothing(year) ? length(time) : findfirst(==(year), time)
        isnothing(n) && error("year $year not in $path (times: $time)")
        (; year  = time[n],
           H     = regrid(field(ds["h_ice"][:, :, n]),   lon, lat),
           z_sur = regrid(field(ds["z_sur"][:, :, n]),   lon, lat),
           rsl   = regrid(field(ds["rsl_ref"][:, :, n]), lon, lat),
           z_bed_ref = regrid(field(ds["z_bed_ref"][:, :]), lon, lat))
    end
end

"""GLAC1D LGM ice thickness, surface elevation and RSL relative to today."""
function load_glac1d(path_lgm, path_0ka)
    zb0 = NCDataset(ds -> field(ds["bedrock_topography"][:, :, 1]), path_0ka)
    NCDataset(path_lgm) do ds
        lon, lat = Float64.(ds["lon"][:]), Float64.(ds["lat"][:])
        zb = field(ds["bedrock_topography"][:, :, 1])
        H  = field(ds["ice_thickness"][:, :, 1])
        (; H     = regrid(H, lon, lat),
           z_sur = regrid(zb .+ H, lon, lat),
           rsl   = regrid(-(zb .- zb0), lon, lat))
    end
end

# --- Plotting ---------------------------------------------------------------

"""RSL fill, present-day coast, ice surface contours and ice margin."""
function panel!(ax, d, z_bed_ref)
    lim = last(RSL_LEVELS) - 0.1
    hm = contourf!(ax, XS, YS, clamp.(d.rsl, -lim, lim);
                   levels = RSL_LEVELS, colormap = :BrBG)
    contour!(ax, XS, YS, z_bed_ref; levels = [0.0], color = (:gray35, 0.8), linewidth = 0.6)
    contour!(ax, XS, YS, ifelse.(d.H .> H_MIN, d.z_sur, NaN);
             levels = ZS_LEVELS, color = :black, linewidth = 0.7)
    contour!(ax, XS, YS, d.H; levels = [H_MIN], color = :black, linewidth = 1.8)
    hidedecorations!(ax); hidespines!(ax)
    return hm
end

function main()
    model = load_model(GEO_HIRES, YEAR)
    glac  = load_glac1d(GLAC_LGM, GLAC_0KA)

    fig = Figure(size = (1300, 720), fontsize = 15)
    ax1 = Axis(fig[1, 1]; aspect = DataAspect(),
               title = @sprintf("%s, year %d", RUN, model.year))
    ax2 = Axis(fig[1, 2]; aspect = DataAspect(), title = "GLAC1D (Tarasov), LGM")
    hm = panel!(ax1, model, model.z_bed_ref)
    panel!(ax2, glac, model.z_bed_ref)
    Colorbar(fig[1, 3], hm; label = "RSL relative to today [m]",
             ticks = first(RSL_LEVELS):200:last(RSL_LEVELS))
    Label(fig[2, 1:2], "Black lines: ice margin (thick) and surface elevation every 500 m. " *
          "Grey: present-day coastline. RSL > 0: land lower relative to sea than today.";
          fontsize = 13, color = :gray30)

    mkpath(FIGDIR)
    out = joinpath(FIGDIR, @sprintf("lgm_vs_glac1d_year%05d.png", model.year))
    save(out, fig; px_per_unit = 2)
    println(out)
end

main()

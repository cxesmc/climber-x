# Generating maps

**CLIMBER-X** couples several components that live on different grids
(atmosphere, ocean, land, ice sheets, solid Earth). To transfer fields between
them, the model uses precomputed **maps** — the interpolation weights that
transform a field from one coordinate grid to another. These maps are generated
internally by the model and stored in the `maps/` directory.

## First run of a configuration is slower

The maps a simulation needs depend on its grid configuration. The first time you
run a given configuration, the model checks `maps/` for the maps it requires and
**generates any that are missing**. This makes that first run noticeably slower
while the maps are built.

The generated maps are written to `maps/` and reused on every subsequent run.
Once they exist, the model finds them there and skips generation entirely, so all
later runs of the same configuration start at full speed. You do not need to
regenerate them unless the grid configuration changes.

## Ensembles: generate maps once, first

Every member of an ensemble writes into the **same shared `maps/` directory**. If
the maps do not yet exist and you launch the whole ensemble at once, all members
start generating the same map files simultaneously, write over each other, and
crash.

The rule is therefore simple: **make sure the maps exist before launching an
ensemble.** Do one map-generating run first, let it populate `maps/`, and only
then launch the ensemble — every member will then find the maps already in place
and run normally.

## `runme --gen-maps`

To make this easy, `runme` provides the `--gen-maps` option. It runs exactly one
map-generating simulation and then stops, without launching the real run or
ensemble.

- For a **single simulation**, it runs that one simulation.
- For an **ensemble**, it runs ensemble **member 0** only.

The map-generating run is placed in a `mapgen/` subdirectory of the output
directory (`<OUTDIR>/mapgen/`), so it stays out of the way of the real run
directories. Its only purpose is to populate the shared `maps/` directory.

The intended workflow is to issue the *same command twice* — first with
`--gen-maps` appended to build the maps, then again without it for the real run.
Keep the command exactly the same and only add or drop the flag:

```bash
# 1. Generate the maps (runs one simulation into OUTDIR/mapgen/)
runme -rs -o OUTDIR -p ctl.n_accel=1,5,10 --gen-maps

# 2. Run the real ensemble (maps now already exist and are reused)
runme -rs -o OUTDIR -p ctl.n_accel=1,5,10
```

Because it is your real command with the flag simply appended, the run flags
(`-r`, `-s`, etc.) must still be present — `--gen-maps` reuses them to actually
run the map-generating simulation. Once it has finished and the maps are in
`maps/`, drop the flag and launch the full ensemble.

For a single simulation the same two-step pattern works, though it is only needed
the first time a new configuration is run:

```bash
runme -rs -o RUNDIR --gen-maps    # build maps into RUNDIR/mapgen/
runme -rs -o RUNDIR               # real run, reusing the maps
```

See [Running with `runme`](runme-notes.md) for the full set of `runme` options.

## Map format and generating maps offline

**CLIMBER-X** maps use the **SCRIP** format — a widely used, self-describing
NetCDF convention for interpolation weights between two grids. Because it is a
standard format, the maps the model generates internally are interchangeable with
weights produced by other SCRIP-aware tools, most notably
[CDO](https://code.mpimet.mpg.de/projects/cdo/).

This means you can also **generate maps offline** and drop them into `maps/`. As
long as a weights file is correct (built between the right source and target
grids) and written in SCRIP format, **CLIMBER-X** will use it directly instead of
regenerating it. The model only builds a map when the corresponding file is
missing from `maps/`.

For example, CDO can produce first-order conservative remapping weights in SCRIP
format with `gencon`:

```bash
# Conservative weights from a source grid to a target grid, in SCRIP format
cdo gencon,target_grid.nc source_grid.nc weights_src_to_tgt.nc
```

Here `source_grid.nc` is a file on (or describing) the source grid,
`target_grid.nc` describes the target grid, and `weights_src_to_tgt.nc` is the
resulting SCRIP weights file. Place it in `maps/` under the filename the model
expects for that grid pair, and it will be picked up in place of an
internally-generated map.

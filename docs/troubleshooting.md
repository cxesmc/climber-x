# Troubleshooting CLIMBER-X

## Strange NetCDF writing error

Solution: 
## Segfault in geo initialisation on macOS (no backtrace)

On Apple silicon the main-thread stack is capped at 64 MB (`ulimit -s hard`), which `grid_init` on the 0.125° geo grid exceeds; the model dies right after the `RTopo ... Points summary` print with `Segmentation fault` and no backtrace even in a debug build. Raise the stack at link time (512 MB is the arm64 linker maximum), and set the OpenMP stack:

```bash
make climber-clim-ice LFLAGS_EXTRA='-Wl,-stack_size,0x20000000'
ulimit -s hard
export OMP_STACKSIZE=512M
```

Put `LFLAGS_EXTRA = -Wl,-stack_size,0x20000000` in the configme `macbook` machine fragment to make it permanent (that fragment also has to clear the default `-Wl,-zmuldefs`, which Apple's ld rejects).

Related macOS pitfalls: a run killed while generating a map leaves a truncated `maps/*.nc` that segfaults every later run on load (delete it); overwriting `climber.x` in place in a run directory can make macOS kill the new binary at launch with no output (`rm` it first, then copy); a debug build's `-ffpe-trap=underflow` trips inside fesm-utils map generation, so use a release build for runs.

# Running CLIMBER-X with `runme`

## General steps to prepare a simulation

After **CLIMBER-X** has been compiled, several steps must be completed to run the executable `climber.x`:

1. Create a run directory (`RUNDIR`).
2. Copy namelist parameter files to `RUNDIR`.
3. Make links to the `input`, `maps` and `restart` directories in `RUNDIR`.
4. Copy `VILMA` restart files to `RUNDIR` (since these are eventually modified by `climber.x`).
5. Copy the executable file `climber.x` to `RUNDIR`.
6. To run on the cluster, create a job submission script (e.g. `job.submit`) in `RUNDIR` to manage various computing options (number of processors, etc).

When these steps are completed, `climber.x` should be run directly from `RUNDIR` using the job submission script. E.g., using SLURM:

```bash
cd RUNDIR
sbatch job.submit
```

If not using the cluster, **CLIMBER-X** can also be run manually by entering the `RUNDIR` and running the following:

```bash
cd RUNDIR
./climber.x > out.out
```

Either of the above will run climber in `RUNDIR` with all simulation output stored in the same directory. Because the directory includes the executable and is self-contained, assuming the contents of the linked directories do not change, the simulation can be re-run at any time.

To perform all of the above steps manually for multiple simulations is very tedious and time-consuming, so a "run" script called `runme` is used to handle the process.

## Using `runme` for one simulation

Before using `runme`, the user should store a configuration file with some personal choices. To get started, copy the template config file to the main directory:

```bash
cp .runme/runme_config .runme_config
```

Next edit `.runme_config` so that the choices match your configuration. Mostly, this means setting `hpc` to the name of your current system and `account` to the default account to be used for your jobs submitted via SLURM. Also you can add your email address if you would like to receive notifications from SLURM about your jobs.

Now `runme` is ready for use. See `runme -h` for details on possible arguments.

Note that aside from possible optional arguments, `runme` is always called with the required argument `-o RUNDIR` that specifies the output directory. So, to run a `climber.x` simulation as a job on the cluster in `RUNDIR`, run the command:

```bash
runme -rs -o RUNDIR
```

where `RUNDIR` is the desired run directory. The option `-r` says that the job should actually be run (instead of just prepared) and `-s` specifies that the job should be run on the cluster. If `-r` is used alone, then the job is simply run as a background process. Other options include `-q, --queue` for the queue alias (short, priority, etc.), `-w, --wall` for the maximum wall clock time to allow in format HH:MM:SS, `--part` to name the processor partition (priority, standard, smp, etc), `--omp` to specify the number of processors, and others. See `runme -h` for all options.

When called as above, this script will run a simulation using the parameters as they are specified in the namelist parameter files in the `nml` directory. In addition, it is possible to modify the parameters of one simulation at the command line using the argument `-p KEY=VAL KEY=VAL ...`. So, for example, the following command:

```bash
runme -s -o RUNDIR -p ctl.n_accel=10
```

will run `climber.x` on the cluster in the output directory `RUNDIR` with the control parameter `control.n_accel` set to `10`. Note that `ctl` is a convenient alias for the namelist group `control`, as defined in `.runme/climberx_info.json`.

To run an ensemble of simulations with modified parameter values, `runme` handles the ensemble internally (see below) — no external tool is required.

## Running ensembles with `runme`

`runme` switches to ensemble mode automatically whenever a `-p` value is *ensemble-shaped* — a comma list (`a=1,2,3`), a range (`a=0:10:5`), or a distribution (`a=U?0,1`) — or whenever an ensemble parameter file is supplied with `-i FILE`. A single-valued `-p` entry is instead a *fixed override* applied to every member.

When running an ensemble, the `-o` argument no longer names a single `RUNDIR` but an encapsulating experiment directory `OUTDIR` that will contain one run directory per member.

The simplest ensemble varies one parameter over a comma list:

```bash
runme -rs -o OUTDIR -p ctl.n_accel=1,5,10
```

This ensemble of simulations will appear in `OUTDIR/0`, `OUTDIR/1` and `OUTDIR/2`, respectively.

A more informative directory naming can be obtained using the option `-a` (auto-dir) along with `-o`:

```bash
runme -rs -a -o OUTDIR -p ctl.n_accel=1,5,10
```

In this case, the run directories are named from the parameter values (group prefix dropped, vowels removed): `OUTDIR/nccl.1`, `OUTDIR/nccl.5` and `OUTDIR/nccl.10`, respectively.

General information about the ensemble can be found in the main experiment directory `OUTDIR`:

- `params.txt` : contains a table of the parameter combinations set on the command line (can be used to run a new ensemble with `-i`).
- `info.txt` : the same parameter table as `params.txt`, but also including an index of the `runid` (0,1,2, etc) and the `rundir`:

`info.txt`:

```python
  runid    ctl.n_accel  rundir
      0              1  nccl.1
      1              5  nccl.5
      2             10  nccl.10
```

It is of course possible to vary multiple parameters at once; `runme` takes the product of all ensemble-shaped dimensions:

```bash
runme -rs -o OUTDIR -p ctl.n_accel=1,5,10 smb.alb_ice=0.3,0.4
```

Single-valued entries can be mixed in as fixed overrides applied to every member:

```bash
runme -rs -o OUTDIR -p ctl.n_accel=1,5,10 ctl.flag_bgc=T
```

To generate a more complex ensemble, using e.g. Latin-Hypercube sampling, a two-step approach is often better. First use the `runme sample` subcommand to build the ensemble parameter file, then run it with `-i`:

```bash
# Generate ensemble parameters (LHS by default)
runme sample -o lhs.txt --seed 4 -N 100 atm.c_trop_2=U?0.8,1.2 smb.alb_ice=U?0.3,0.4

# Run ensemble
runme -rs -o OUTDIR -i lhs.txt
```

Here `U?0.8,1.2` denotes a uniform distribution between 0.8 and 1.2 (use `N?mean,std` for a normal distribution, or a plain comma list for discrete values). The companion `runme product` subcommand writes a full-factorial parameter file instead of sampling. This two-step method facilitates checking that the ensemble was generated properly and improves reproducibility, since the exact parameter values are available in the table.

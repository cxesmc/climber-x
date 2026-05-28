# Use case: MIS-3

(General description)

## Basic run

```bash
outdir=mis3
runme -rs -q 48h -w 24:00:00 --omp 32 -o output/$outdir/h4 -p ctl.nyears=5000 ctl.i_write_restart=1 ctl.iorbit=0 ctl.ifake_ice=2 ctl.ifake_geo=2 ctl.ico2=-1 ctl.year_ini=-39e3
```

## With Talento SMB

```bash
outdir=mis3-talento-smb
runme -rs -q 48h -w 24:00:00 --omp 32 -o output/$outdir/h4 -p ctl.nyears=5000 ctl.i_write_restart=1 ctl.iorbit=0 ctl.ifake_ice=2 ctl.ifake_geo=2 ctl.ico2=-1 ctl.year_ini=-39e3 smb.i_smb=3
```

## Basic run with ice sheets active

CRASHES!!

```bash
outdir=mis3-nh
runme -rs -q 48h -w 24:00:00 --omp 32 -o output/$outdir/h4 -p ctl.nyears=5000 ctl.i_write_restart=1 ctl.iorbit=0 ctl.ifake_ice=2 ctl.ifake_geo=2 ctl.ico2=-1 ctl.year_ini=-39e3 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM
```

## With Talento SMB

```bash
outdir=mis3-nh-talento-smb
runme -rs -q 48h -w 24:00:00 --omp 32 -o output/$outdir/h4 -p ctl.nyears=5000 ctl.i_write_restart=1 ctl.iorbit=0 ctl.ifake_ice=2 ctl.ifake_geo=2 ctl.ico2=-1 ctl.year_ini=-39e3 smb.i_smb=3  ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM
```

(Explain parameters choices, etc)

## Typical output

## Alternate configurations

...

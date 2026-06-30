module fastearth_model
  !! CLIMBER-X wrapper for the FastEarth3D solid-Earth / GIA model, a swappable
  !! alternative to VILMA selected at runtime with i_geo=3 (see src/geo/geo.f90).
  !!
  !! Mirrors the public surface of vilma_model (src/geo/vilma.F90): ice thickness
  !! goes in on the CLIMBER-X high-resolution geo grid, relative sea level and
  !! bedrock elevation come out on the same grid. Unlike the VILMA wrapper, the
  !! grid remapping (host lon-lat <-> the model's Gauss grid) is done INSIDE
  !! FastEarth3D — we just hand it geo_grid at init and read se%rsl / se%z_bed.
  !!
  !! All FastEarth3D-dependent code is guarded by #ifdef FASTEARTH so the module
  !! compiles to empty stubs in build variants that do not link libfastearth.a.

  use precision, only : wp
  use timer, only : n_year_geo
  use control, only : out_dir, geo_restart, restart_in_dir
  use coords, only : grid_class

#ifdef FASTEARTH
  ! FastEarth3D public API. The umbrella module re-exports the types, parameter
  ! record and restart I/O; the coupling procedures come from fe_coupling.
  use fastearth3d, only : solid_earth, fe_par_load, fe_restart_write, fe_restart_read
  use fe_coupling, only : solid_earth_init, solid_earth_update, solid_earth_finalize
#endif

  implicit none

  private
  public :: fastearth_init, fastearth_update, fastearth_end
  public :: fastearth_write_restart

#ifdef FASTEARTH
  ! single model instance held for the lifetime of the run
  type(solid_earth), save :: se
#endif

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f a s t e a r t h _ i n i t
  ! Purpose  :  initialize solid Earth (FastEarth3D)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fastearth_init(geo_grid, z_bed_eq_g, h_ice_eq_g, h_ice_g)

    implicit none

    type(grid_class), intent(in) :: geo_grid
    real(wp), intent(in) :: z_bed_eq_g(:,:)   ! equilibrium (relaxed) bedrock elevation [m]
    real(wp), intent(in) :: h_ice_eq_g(:,:)   ! ice thickness associated with equilibrium bedrock [m]
    real(wp), intent(in) :: h_ice_g(:,:)      ! current ice thickness [m]

#ifdef FASTEARTH

    ! Load configuration into se%par. The &fe3d group lives in geo_par.nml
    ! (compact overrides); the complete default set is input/fastearth.nml.
    call fe_par_load(se%par, trim(out_dir)//"/geo_par.nml", defaults_file="input/fastearth.nml", group="fe3d")

    ! Build the model and set the relaxed reference state. FastEarth3D owns its
    ! Gauss grid and builds the host<->Gauss remap from geo_grid internally.
    call solid_earth_init(se, z_bed_eq_g, h_ice_eq_g, grid=geo_grid, h_ice_init=h_ice_g)

    ! Restart: restore prognostic state on top of the freshly built model.
    if (geo_restart) then
      call fe_restart_read(se, trim(restart_in_dir)//"/fastearth/fe_restart.nc")
    endif

    print*
    print*,'======================================================='
    print*,' Initialisation of FastEarth3D complete'
    print*,'======================================================='
    print*

#else

    ! FastEarth3D backend not compiled in (built with fastearth=0). Fail fast
    ! rather than silently no-op (which would leave rsl/z_bed unset).
    print*,'======================================================='
    print*,' ERROR: i_geo=3 requires the FastEarth3D solid-earth backend,'
    print*,'        but this binary was built with fastearth=0.'
    print*,'        Rebuild with fastearth=1, or choose i_geo=0/1/2.'
    print*,'======================================================='
    stop 1

#endif

    return

  end subroutine fastearth_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f a s t e a r t h _ u p d a t e
  ! Purpose  :  update solid Earth (FastEarth3D)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fastearth_update(h_ice_g, rsl_g, z_bed_g)

    implicit none

    real(wp), intent(in) :: h_ice_g(:,:)   ! current ice thickness [m]
    real(wp), intent(out) :: rsl_g(:,:)    ! relative sea level [m]
    real(wp), intent(out) :: z_bed_g(:,:)  ! updated bedrock elevation [m]

#ifdef FASTEARTH

    ! advance by one geo coupling interval (FastEarth3D takes the interval in years)
    call solid_earth_update(se, h_ice_g, real(n_year_geo,wp))

    rsl_g = se%rsl
    z_bed_g = se%z_bed

#endif

    return

  end subroutine fastearth_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f a s t e a r t h _ e n d
  ! Purpose  :  end solid Earth (FastEarth3D)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fastearth_end

    implicit none

#ifdef FASTEARTH

    call solid_earth_finalize(se)

#endif

    return

  end subroutine fastearth_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f a s t e a r t h _ w r i t e _ r e s t a r t
  ! Purpose  :  write FastEarth3D restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fastearth_write_restart(dir)

    implicit none

    character (len=*) :: dir

#ifdef FASTEARTH

    integer :: cstat, estat
    character(len=256) :: cmsg

    ! ensure the restart subdirectory exists, then write the single-file restart
    call execute_command_line('mkdir -p '//trim(dir)//'/fastearth', &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
    call fe_restart_write(se, se%time, filename="fe_restart.nc", folder=trim(dir)//"/fastearth")

#endif

    return

  end subroutine fastearth_write_restart


end module fastearth_model

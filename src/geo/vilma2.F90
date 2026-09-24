module vilma2_model
  !! CLIMBER-X wrapper for the VILMA2 (fesmc/vilma) solid-Earth / GIA model, a swappable
  !! alternative to VILMA1 selected at runtime with i_geo=3 (see src/geo/geo.f90).
  !!
  !! Mirrors the public surface of vilma1_model (src/geo/vilma1.F90): ice thickness
  !! goes in on the CLIMBER-X high-resolution geo grid, relative sea level and
  !! bedrock elevation come out on the same grid. Unlike the VILMA1 wrapper, the
  !! grid remapping (host lon-lat <-> the model's Gauss grid) is done INSIDE
  !! VILMA2 — we just hand it geo_grid at init and read se%rsl / se%z_bed.
  !!
  !! All VILMA2-dependent code is guarded by #ifdef VILMA2 so the module
  !! compiles to empty stubs in build variants that do not link libvilma.a.

  use precision, only : wp
  use timer, only : n_year_geo
  use control, only : out_dir, geo_restart, restart_in_dir
  use coords, only : grid_class

#ifdef VILMA2
  ! VILMA2 public API (library modules keep the upstream vilma* names). The umbrella
  ! module re-exports the types, parameter record and restart I/O; the coupling
  ! procedures come from vilma_coupling.
  use vilma, only : solid_earth, vilma_par_load, vilma_restart_write, vilma_restart_read
  use vilma_coupling, only : solid_earth_init, solid_earth_update, solid_earth_finalize
#endif

  implicit none

  private
  public :: vilma2_init, vilma2_update, vilma2_end
  public :: vilma2_write_restart

#ifdef VILMA2
  ! single model instance held for the lifetime of the run
  type(solid_earth), save :: se
#endif

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a 2 _ i n i t
  ! Purpose  :  initialize solid Earth (VILMA2)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma2_init(geo_grid, z_bed_eq_g, h_ice_eq_g, h_ice_g)

    implicit none

    type(grid_class), intent(in) :: geo_grid
    real(wp), intent(in) :: z_bed_eq_g(:,:)   ! equilibrium (relaxed) bedrock elevation [m]
    real(wp), intent(in) :: h_ice_eq_g(:,:)   ! ice thickness associated with equilibrium bedrock [m]
    real(wp), intent(in) :: h_ice_g(:,:)      ! current ice thickness [m]

#ifdef VILMA2

    ! Load configuration into se%par. The &vilma group in geo_par.nml holds the
    ! minimal climber-x overrides; the complete default set is the canonical
    ! VILMA2 defaults file shipped in input/.
    call vilma_par_load(se%par, trim(out_dir)//"/geo_par.nml", &
      defaults_file="input/vilma_defaults.nml", group="vilma")

    ! Build the model and set the relaxed reference state. VILMA2 owns its
    ! Gauss grid and builds the host<->Gauss remap from geo_grid internally.
    call solid_earth_init(se, z_bed_eq_g, h_ice_eq_g, grid=geo_grid, h_ice_init=h_ice_g)

    ! Restart: restore prognostic state on top of the freshly built model.
    if (geo_restart) then
      call vilma_restart_read(se, trim(restart_in_dir)//"/vilma2/vilma_restart.nc")
    endif

    print*
    print*,'======================================================='
    print*,' Initialisation of VILMA2 complete'
    print*,'======================================================='
    print*

#else

    ! VILMA2 backend not compiled in (built with vilma2=0). Fail fast
    ! rather than silently no-op (which would leave rsl/z_bed unset).
    print*,'======================================================='
    print*,' ERROR: i_geo=3 requires the VILMA2 solid-earth backend,'
    print*,'        but this binary was built with vilma2=0.'
    print*,'        Rebuild with vilma2=1, or choose i_geo=0/1/2.'
    print*,'======================================================='
    stop 1

#endif

    return

  end subroutine vilma2_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a 2 _ u p d a t e
  ! Purpose  :  update solid Earth (VILMA2)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma2_update(h_ice_g, rsl_g, z_bed_g)

    implicit none

    real(wp), intent(in) :: h_ice_g(:,:)   ! current ice thickness [m]
    real(wp), intent(out) :: rsl_g(:,:)    ! relative sea level [m]
    real(wp), intent(out) :: z_bed_g(:,:)  ! updated bedrock elevation [m]

#ifdef VILMA2

    ! advance by one geo coupling interval (VILMA2 takes the interval in years)
    call solid_earth_update(se, h_ice_g, real(n_year_geo,wp))

    rsl_g = se%rsl
    z_bed_g = se%z_bed

#endif

    return

  end subroutine vilma2_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a 2 _ e n d
  ! Purpose  :  end solid Earth (VILMA2)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma2_end

    implicit none

#ifdef VILMA2

    call solid_earth_finalize(se)

#endif

    return

  end subroutine vilma2_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a 2 _ w r i t e _ r e s t a r t
  ! Purpose  :  write VILMA2 restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma2_write_restart(dir)

    implicit none

    character (len=*) :: dir

#ifdef VILMA2

    integer :: cstat, estat
    character(len=256) :: cmsg

    ! ensure the restart subdirectory exists, then write the single-file restart
    call execute_command_line('mkdir -p '//trim(dir)//'/vilma2', &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
    call vilma_restart_write(se, se%time, filename="vilma_restart.nc", folder=trim(dir)//"/vilma2")

#endif

    return

  end subroutine vilma2_write_restart


end module vilma2_model

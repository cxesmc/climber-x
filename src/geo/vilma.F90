module vilma_model

    use precision, only : wp
    use timer, only : dt_geo, sec_year

#ifdef VILMA
    ! modules from the VILMA library
    use mod_cons
    use mod_firstlevel
    use mod_struct_vg
    use mod_sle, only : rsl, dfgl
    use mod_io
#endif

    use timer, only : year_ini, nyears, n_year_geo
    use control, only : out_dir, geo_restart, restart_in_dir
    use constants, only : rho_i, rho_sw
    use geo_params, only : vilma_grid_file, l_visc_3d, visc_1d_file, visc_3d_file
    use coord, only : grid_class, grid_init
    use coord, only : map_scrip_class, map_scrip_init, map_scrip_field
    use ncio

    implicit none

    private
    public :: vilma_init, vilma_update, vilma_end
    public :: vilma_read_restart, vilma_write_restart

    integer :: ni, nj
    integer :: iepoch
    real(wp), dimension(:), allocatable :: epoch
    real(wp), dimension(:), allocatable :: lon, lat
    real(wp), dimension(:,:), allocatable :: tmp(:,:)
    real(wp), dimension(:,:), allocatable :: z_bed_ref(:,:)
    real(wp), dimension(:,:), allocatable :: z_bed_ref_g(:,:)
    real(wp), dimension(:,:), allocatable :: h_ice_ref(:,:)
    real(wp), dimension(:,:), allocatable :: h_ice(:,:)

    type(grid_class) :: vilma_grid
    type(map_scrip_class) :: maps_geo_to_vilma, maps_vilma_to_geo


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a _ u p d a t e
  ! Purpose  :  update solid Earth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma_update(geo_grid, h_ice_g, &
      rsl_g, z_bed_g)

    implicit none

    type(grid_class), intent(in) :: geo_grid

    real(wp), intent(in) :: h_ice_g(:,:)
    real(wp), intent(out) :: rsl_g(:,:)
    real(wp), intent(out) :: z_bed_g(:,:)

    real(wp), allocatable, dimension(:,:) :: mask_ice_g
    real(wp), allocatable, dimension(:,:) :: mask_ice

    integer :: stat, t
    character(len=256) :: fnm
    integer :: ncid


#ifdef VILMA
    ! delete old restart files
    open(unit=1, iostat=stat, file=io_rsl_rs%n, status='old')
    if (stat == 0) close(1, status='delete')
    open(unit=1, iostat=stat, file=io_dfgl_rs%n, status='old')
    if (stat == 0) close(1, status='delete')

    ! set begin and end time (in kyrs)
    vg%btime = vg%etime
    vg%etime = vg%btime + dble(n_year_geo*1.d-3)

!    fnm = trim(out_dir)//"/vilma_h_ice.nc"
!    call nc_create(fnm)
!    call nc_open(fnm,ncid)
!    call nc_write_dim(fnm,"epoch",x=epoch(iepoch-1:iepoch), units="ka BP", unlimited=.TRUE.,ncid=ncid)
!    call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
!    call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)
!    call nc_write(fnm,"Ice", h_ice, dims=["lon","lat","epoch"],start=[1,1,1],count=[ni,nj,1],long_name="Ice thickness",units="m",ncid=ncid)

    iepoch = iepoch + 1

    ! interpolate ice thickness to vilma grid
    call map_scrip_field(maps_geo_to_vilma,"hice",h_ice_g,h_ice,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5._dp*geo_grid%G%dx,geo_grid%G%dx])
    allocate(mask_ice_g(geo_grid%G%nx,geo_grid%G%ny))
    allocate(mask_ice(vilma_grid%G%nx,vilma_grid%G%ny))
    where (h_ice_g>0._wp) 
      mask_ice_g = 1.
    elsewhere
      mask_ice_g = 0.
    endwhere
    mask_ice = 1.
    call map_scrip_field(maps_geo_to_vilma,"mask",mask_ice_g,mask_ice,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5._dp*geo_grid%G%dx,geo_grid%G%dx])
    where (mask_ice<0.5) h_ice = 0._wp
    deallocate(mask_ice_g)
    deallocate(mask_ice)

    ! append h_ice to ice thickness netcdf file
    fnm = trim(out_dir)//"/vilma_h_ice.nc"
    call nc_open(fnm,ncid)
    call nc_write(fnm,"Ice", h_ice, dims=["lon","lat","epoch"],start=[1,1,iepoch],count=[ni,nj,1],long_name="Ice thickness",units="m",ncid=ncid)
    call nc_close(ncid)

!    call nc_write(fnm,"Ice", h_ice, dims=["lon","lat","epoch"],start=[1,1,2],count=[ni,nj,1],long_name="Ice thickness",units="m",ncid=ncid)
!    call nc_close(ncid)

    !print *,'ice vol vilma',sum(h_ice*vilma_grid%area, mask=vilma_grid%lat>0._wp)

    write (6,*) 'call time_evolution for ', vg%btime, 'to', vg%etime
    call time_evolution

    ! interpolate from Gauss-Legendre vilma_grid to geo_grid (regular lat-lon)
    ! relative sea level
    call map_scrip_field(maps_vilma_to_geo,"rsl",transpose(real(rsl,wp)),rsl_g,method="mean",missing_value=-9999._dp)

    ! update bedrock elevation
    z_bed_g = z_bed_ref_g - rsl_g

    !fnm = trim(out_dir)//"/test.nc"
    !call nc_create(fnm)
    !call nc_open(fnm,ncid)
    !call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
    !call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)
    !call nc_write(fnm,"dfgl", transpose(dfgl) ,dims=["lon","lat"],start=[1,1],count=[ni,nj],ncid=ncid)
    !call nc_write(fnm,"rsl", transpose(rsl) ,dims=["lon","lat"],start=[1,1],count=[ni,nj],ncid=ncid)
    !call nc_close(ncid)

    !fnm = trim(out_dir)//"/test_g.nc"
    !call nc_create(fnm)
    !call nc_open(fnm,ncid)
    !call nc_write_dim(fnm,"lon",x=geo_grid%G%x,axis="x",ncid=ncid)
    !call nc_write_dim(fnm,"lat",x=geo_grid%G%y,axis="y",ncid=ncid)
    !call nc_write(fnm,"rsl_g", rsl_g ,dims=["lon","lat"],start=[1,1],count=[geo_grid%G%nx,geo_grid%G%ny],ncid=ncid)
    !call nc_write(fnm,"h_ice_g", h_ice_g ,dims=["lon","lat"],start=[1,1],count=[geo_grid%G%nx,geo_grid%G%ny],ncid=ncid)
    !call nc_write(fnm,"z_bed_g", z_bed_g ,dims=["lon","lat"],start=[1,1],count=[geo_grid%G%nx,geo_grid%G%ny],ncid=ncid)
    !call nc_close(ncid)

#endif

   return

  end subroutine vilma_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a _ i n i t
  ! Purpose  :  initialize solid Earth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma_init(geo_grid, z_bed_ref_in, h_ice_ref_g, h_ice_g)

    implicit none

    type(grid_class), intent(in) :: geo_grid
    real(wp), intent(in) :: z_bed_ref_in(:,:)
    real(wp), intent(in) :: h_ice_ref_g(:,:)
    real(wp), intent(in) :: h_ice_g(:,:)

    real(wp), allocatable, dimension(:,:) :: mask_ice_g
    real(wp), allocatable, dimension(:,:) :: mask_ice

    integer :: t, nepoch, stat
    character(len=256) :: fnm
    integer :: ncid

    integer :: cstat, estat
    character(len=256) :: cmsg

#ifdef VILMA
    
    ni = nc_size(trim(vilma_grid_file),"lon")
    nj = nc_size(trim(vilma_grid_file),"lat")
    allocate( lon(ni) )
    allocate( lat(nj) )
    call nc_read(trim(vilma_grid_file),"lon",lon)
    call nc_read(trim(vilma_grid_file),"lat",lat)
   
    allocate(tmp(ni,nj))
    allocate(z_bed_ref(ni,nj))
    allocate(z_bed_ref_g(geo_grid%G%nx,geo_grid%G%ny))
    allocate(h_ice_ref(ni,nj))
    allocate(h_ice(ni,nj))

    ! generate grid object
    call grid_init(vilma_grid,name="vilma_grid",mtype="latlon",units="degrees", x=real(lon,dp),y=real(lat,dp), lon180=.true.)

    ! generate maps for mapping between geo and vilma
    call map_scrip_init(maps_geo_to_vilma,geo_grid,vilma_grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)
    call map_scrip_init(maps_vilma_to_geo,vilma_grid,geo_grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)

    ! initialize VILMA

    ! definition of general parameters (replaces call to r_cli subroutine where parameteres are read from stdin)
    ! 1.jmin and jmax // spectral resolution
    vg%jmin=0
    vg%jmax=170
    !2. time step, if 0, is read from tint // integration interval
    !vg%dt=0
    vg%dt=dt_geo  ! s
    !3. polynomial representation of PREM // earth structure is generated in program
    vg%l_prem=1
    !4. modified 3d viscosity (0/1/2) = (no, read, adjust) // character of lateral viscosity structure
    if (l_visc_3d) then
      vg%l_mod=1
    else
      vg%l_mod=0
    endif
    !5. toroidal loading // toroidal loading is not relevvant for GIA
    vg%l_toro=0
    !6. variation of rotation // rotational variations are considered in perturbation of potential
    vg%l_rot=31 
    !7. load grid type (3/2/1/0) // character of loading structure (load love numbers / grid / spectral grid / disc load
    vg%l_grid=2
    !8. =1, will stop after reading of environment // for testing, only input is read in
    vg%l_envonly=0
    !9. number of timesteps, if =0, whole history // number of integration steps
    !vg%ntime=5000
    vg%ntime=10000000
    !10. writing out is not each 1 kyr, but epochs are defined in wepochs.inp // defines for which epochs the output fields are written
    vg%l_wepoch=1

    ! adapt names of input/output files
    io_densi%n  = 'input/vilma/densi.inp' ! standard input of density structure, mainly radial discretisation
    io_visko%n  = 'input/vilma/'//trim(visc_1d_file) ! standard input of 1d viscosity structure
    io_nc3in%n  = 'input/vilma/'//trim(visc_3d_file) ! 3d viscosity input file
    io_tint%n   = 'input/vilma/tint.inp' ! if dt == 0
    io_sliin%n  = 'input/vilma/SLI_data.inp' ! units in m_sle (l_grid == 2)

    io_tmp%n     = trim(out_dir)//'/io.tmp' ! new temporary file for stdin
    io_lis%n     = trim(out_dir)//'/vega.lis' ! with lis_{param, lis_array} tracable log file
    io_wepoch%n  = trim(out_dir)//'wepochs.inp' ! if l_wepochs
    io_hist%n    = trim(out_dir)//'load_hist.inp'  ! units in m_sle (l_grid == 2)
    io_surf%n    = trim(out_dir)//'/vega1.nc' ! if l_wsurf output for spectral displacement field
    io_radii%n   = trim(out_dir)//'/radii.dat'  ! output of ascii of considered radial grid lines
    io_nc3out%n  = trim(out_dir)//'/visc3d.nc'  ! output of 3d viscosity file
    io_sliout%n  = trim(out_dir)//'/SLI_data.out' ! units in m_sle (l_grid == 2) if ikmax
    io_hsliout%n = trim(out_dir)//'/hSLI_data.out' ! units in m_sle (l_grid == 2) if ikmax
    io_rsl%n     = trim(out_dir)//'/rsl.nc' ! new io_rsl
    io_dfgl%n    = trim(out_dir)//'/dflag.nc' ! new io_dfgl
    io_rsl_rs%n  = trim(out_dir)//'rsl_rs.nc' ! for restart
    io_dfgl_rs%n = trim(out_dir)//'dflg_rs.nc' ! for restart
    io_deg1%n    = trim(out_dir)//'vega_deg1.dat' ! deg1 coefficients and cof_[xyz]
    io_oce%n     = trim(out_dir)//'vega_oce.dat' ! units in m_sle (l_grid == 2)
    io_rpt%n     = trim(out_dir)//'vega_rpt.dat' ! units in m_sle (l_grid == 2) if irot
    io_rslog%n   = trim(out_dir)//'restart.log'  ! new io_units for restart

    ! delete some netcdf files if already present in output directory (from old run)
    open(unit=1, iostat=stat, file=io_surf%n, status='old')
    if (stat == 0) close(1, status='delete')
    open(unit=1, iostat=stat, file=io_rsl%n, status='old')
    if (stat == 0) close(1, status='delete')
    open(unit=1, iostat=stat, file=io_dfgl%n, status='old')
    if (stat == 0) close(1, status='delete')

    io_mos_indx%n   = trim(restart_in_dir)//'/vilma//mos_indx.nc'  !matrix of Galerkin System
    io_mos_amtrx%n  = trim(restart_in_dir)//'/vilma//mos_amtrx.nc' !matrix of Galerkin System
    io_mos_acomp%n  = trim(restart_in_dir)//'/vilma//mos_acomp.nc' !matrix of Galerkin System
    io_mos_acompl%n = trim(restart_in_dir)//'/vilma//mos_acompl.nc'  !matrix of Galerkin System
    io_ve_struct%n  = trim(restart_in_dir)//'/vilma//ve_struct.nc' !restart array for ve structure
    io_pjj%n        = trim(restart_in_dir)//'/vilma//pjj.nc' !restart array for associated scalar sph
    io_pefgh%n      = trim(restart_in_dir)//'/vilma//pefgh.nc' !restart array for associated tensor sph
    io_nwl_struct%n = trim(restart_in_dir)//'/vilma//nwl_struct.nc'  !restart array for spatial discretisation in lon and lat
    io_visc3drs%n   = trim(restart_in_dir)//'/vilma//visc3d.nc'  !restart array 3d ve structure
    io_disp%n       = trim(restart_in_dir)//'/vilma//disp.nc'  !restart array for spectal displacement field
    io_stress%n     = trim(restart_in_dir)//'/vilma//stress.nc'  !restart array for spatial stress components (3d)
    io_1dstress%n   = trim(restart_in_dir)//'/vilma//ctc_stress.nc'  !restart array for sepctral stress components (1d)
    
    ! write io.tmp file with parameters, needed for restart capability
    open (1,file=io_tmp%n,form='formatted')
    write(1,fmt=*) vg%jmin, vg%jmax
    write(1,fmt=*) vg%dt
    write(1,fmt=*) vg%l_prem
    write(1,fmt=*) vg%l_mod
    write(1,fmt=*) vg%l_toro
    write(1,fmt=*) vg%l_rot
    write(1,fmt=*) vg%l_grid
    write(1,fmt=*) vg%l_envonly
    write(1,fmt=*) vg%ntime
    write(1,fmt=*) vg%l_wepoch
    close(1)

    ! specify load history file
    open (1,file=io_hist%n,form='formatted')
    write(1,fmt='(A)')trim(out_dir)//"loadh.inp"
    write(1,fmt='(A)')trim(out_dir)//"vilma_h_ice.nc"
    write(1,fmt='(A)',advance='no')trim(out_dir)//"vilma_z_bed_ref.nc"//" "//trim(out_dir)//"vilma_h_ice_ref.nc"
    close(1)

    ! epoch read from loadh.inp, then it knows which time slice in the netcdf in
    ! the ice history file it takes.
    !...has to be consistent with vilma_h_ice.nc

    ! check where epochs are read from!

    if (.not.geo_restart) then
      ! start from scratch
      vg%btime = 9999._dp
      vg%etime = dble(year_ini)*1.d-3
      vg%restart = .false. !(vg%btime /= 9999.0_dp) 
    else
      ! start from restart
      vg%btime = dble(year_ini)*1.d-3
      vg%etime = vg%btime !+ dble(n_year_geo)*1.d-3
      vg%restart = .true.
    endif
    print *,'btime, etime: ',vg%btime,vg%etime

    ! create specification file for ice load history input
    nepoch = nyears/n_year_geo+1
    allocate(epoch(1:nepoch))
    epoch(1) = dble(year_ini)*1.d-3
    do t=2,nepoch
      epoch(t) = epoch(t-1) + dble(n_year_geo)*1.d-3
    enddo

    open(1,file=trim(out_dir)//"loadh.inp",form='formatted')
    write(1,fmt=*) nj, ni, 910, 1020
    write(1,fmt=*) nepoch
    do t=1,nepoch
      write(1,fmt=*) epoch(t)
    enddo
    close(1)

    ! specify times to write vilma internal output
    open(1,file=io_wepoch%n,form='formatted')
    write(1,fmt=*) dble(year_ini*1.d-3)
    write(1,fmt=*) dble((year_ini+nyears)*1.d-3)
    close(1)

    iepoch = 1


    z_bed_ref_g = z_bed_ref_in
    ! interpolate z_bed_ref and h_ice to vilma grid
    call map_scrip_field(maps_geo_to_vilma,"zref",z_bed_ref_g,z_bed_ref,method="mean",missing_value=-9999._dp, &
      filt_method="none")
    call map_scrip_field(maps_geo_to_vilma,"hice",h_ice_ref_g,h_ice_ref,method="mean",missing_value=-9999._dp, &
      filt_method="none")
    call map_scrip_field(maps_geo_to_vilma,"hice",h_ice_g,h_ice,method="mean",missing_value=-9999._dp, &
      filt_method="none")
    allocate(mask_ice_g(geo_grid%G%nx,geo_grid%G%ny))
    allocate(mask_ice(vilma_grid%G%nx,vilma_grid%G%ny))
    where (h_ice_ref_g>0._wp) 
      mask_ice_g = 1.
    elsewhere
      mask_ice_g = 0.
    endwhere
    mask_ice = 1.
    call map_scrip_field(maps_geo_to_vilma,"mask",mask_ice_g,mask_ice,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5.0*geo_grid%G%dx,geo_grid%G%dx])
    where (mask_ice<0.5) h_ice_ref = 0._wp
    where (h_ice_g>0._wp) 
      mask_ice_g = 1.
    elsewhere
      mask_ice_g = 0.
    endwhere
    mask_ice = 1.
    call map_scrip_field(maps_geo_to_vilma,"mask",mask_ice_g,mask_ice,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5.0*geo_grid%G%dx,geo_grid%G%dx])
    where (mask_ice<0.5) h_ice = 0._wp
    deallocate(mask_ice_g)
    deallocate(mask_ice)

    ! write reference topography to netcdf file
    fnm = trim(out_dir)//"/vilma_z_bed_ref.nc"
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,"epoch",x=1._dp, units="ka BP", unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)
    call nc_write(fnm,"topo", z_bed_ref, dims=["lon","lat","epoch"],start=[1,1,1],count=[ni,nj,1],long_name="reference topography",units="m",ncid=ncid)
    call nc_close(ncid)

    ! write reference ice thickness to netcdf file
    fnm = trim(out_dir)//"/vilma_h_ice_ref.nc"
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,"epoch",x=1._dp, units="ka BP", unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)
    call nc_write(fnm,"Ice", h_ice_ref, dims=["lon","lat","epoch"],start=[1,1,1],count=[ni,nj,1],long_name="reference ice thickness",units="m",ncid=ncid)
    call nc_close(ncid)

    ! write h_ice to netcdf file
    fnm = trim(out_dir)//"/vilma_h_ice.nc"
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,"epoch",x=epoch(1:), units="ka BP", unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)
    !call nc_write(fnm,"Ice", h_ice, dims=["lon","lat","epoch"],start=[1,1,iepoch],count=[ni,nj,1],long_name="Ice thickness",units="m",ncid=ncid)
    call nc_write(fnm,"Ice", h_ice_ref, dims=["lon","lat","epoch"],start=[1,1,iepoch],count=[ni,nj,1],long_name="Ice thickness",units="m",ncid=ncid)
    call nc_close(ncid)


    if (vg%restart) then

      ! copy some files from vilma restart directory to output directory
      call execute_command_line('cp ' // adjustl(trim(restart_in_dir))//'/vilma/vega_deg1.dat '//trim(out_dir)//'/', &
        exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
      call execute_command_line('cp ' // adjustl(trim(restart_in_dir))//'/vilma/vega_oce.dat '//trim(out_dir)//'/', &
        exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
      call execute_command_line('cp ' // adjustl(trim(restart_in_dir))//'/vilma/vega_rpt.dat '//trim(out_dir)//'/', &
        exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
      call execute_command_line('cp ' // adjustl(trim(restart_in_dir))//'/vilma/restart.log '//trim(out_dir)//'/', &
        exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)

      ! read restart
      call r_restart        ! reads restart.log and checks for consistency between new btime and old etime

    else

      call setup

    endif
    
    print*
    print*,'======================================================='
    print*,' Initialisation of VILMA complete'
    print*,'======================================================='
    print*

#endif

    return

  end subroutine vilma_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a _ end
  ! Purpose  :  end solid Earth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma_end

    implicit none


#ifdef VILMA

    call close_evolution

#endif

    return

  end subroutine vilma_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma_write_restart(dir)

    implicit none

    character (len=*) :: dir

    integer :: cstat, estat
    character(len=256) :: cmsg

#ifdef VILMA
    
    ! rename restart files for restart output writing location
    io_mos_indx%n   = trim(dir)//'/vilma/mos_indx.nc'  !matrix of Galerkin System
    io_mos_amtrx%n  = trim(dir)//'/vilma/mos_amtrx.nc' !matrix of Galerkin System
    io_mos_acomp%n  = trim(dir)//'/vilma/mos_acomp.nc' !matrix of Galerkin System
    io_mos_acompl%n = trim(dir)//'/vilma/mos_acompl.nc'  !matrix of Galerkin System
    io_ve_struct%n  = trim(dir)//'/vilma/ve_struct.nc' !restart array for ve structure
    io_pjj%n        = trim(dir)//'/vilma/pjj.nc' !restart array for associated scalar sph
    io_pefgh%n      = trim(dir)//'/vilma/pefgh.nc' !restart array for associated tensor sph
    io_nwl_struct%n = trim(dir)//'/vilma/nwl_struct.nc'  !restart array for spatial discretisation in lon and lat
    io_visc3drs%n   = trim(dir)//'/vilma/visc3d.nc'  !restart array 3d ve structure
    io_disp%n       = trim(dir)//'/vilma/disp.nc'  !restart array for spectal displacement field
    io_stress%n     = trim(dir)//'/vilma/stress.nc'  !restart array for spatial stress components (3d)
    io_1dstress%n   = trim(dir)//'/vilma/ctc_stress.nc'  !restart array for sepctral stress components (1d)

    ! write restart
    call w_restart

    ! copy some files to vilma restart directory
    call execute_command_line('cp ' // adjustl(trim(out_dir))//'vega_deg1.dat '//trim(dir)//'/vilma/', &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
    call execute_command_line('cp ' // adjustl(trim(out_dir))//'vega_oce.dat '//trim(dir)//'/vilma/', &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
    call execute_command_line('cp ' // adjustl(trim(out_dir))//'vega_rpt.dat '//trim(dir)//'/vilma/', &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
    call execute_command_line('cp ' // adjustl(trim(out_dir))//'restart.log '//trim(dir)//'/vilma/', &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)

#endif

   return

  end subroutine vilma_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  v i l m a _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine vilma_read_restart(fnm)

    implicit none

    character (len=*) :: fnm


    !call nc_read(fnm,"uc",     sic%uc)


   return

  end subroutine vilma_read_restart


end module vilma_model


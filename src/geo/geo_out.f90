!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : g e o _ o u t
!
!  Purpose : geography model diagnostics and output
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! CLIMBER-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module geo_out

  use precision, only : wp
  use constants, only : rho_w
  use dim_name, only: dim_lon, dim_lat, dim_time
  use timer, only : year, year_geo, year_ini, year_now, nyears, sec_day, sec_year, nyout_geo, &
  n_year_geo, time_out_geo, ny_out_ts, y_out_ts_geo, time_out_ts_geo
  use control, only: flag_geo, ifake_geo, flag_lakes, out_dir
  use geo_params, only : l_output_hires, n_coast_cells, i_geo
  use geo_def, only : geo_class
  use ncio

  implicit none

  integer :: nout

  type ts_out
     real(wp) :: sea_level, A_bering, Aocn, Aveg, Aice, Alake
  end type

  type(ts_out), allocatable :: ann_ts(:)

  private
  public :: geo_diag_init, geo_diag

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g e o _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for geography
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_diag_init(geo)

    implicit none

    type(geo_class), intent(in) :: geo

    character (len=256) :: fnm
    integer :: ncid, i
    real(wp) :: empty_time(0)


    nout = 0

    ! allocate
    allocate(ann_ts(ny_out_ts))

    fnm = trim(out_dir)//"/geo_ts.nc"

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.)
    call nc_write_dim(fnm,dim_lat, x=1, axis="y", units="1")
    call nc_write_dim(fnm,dim_lon, x=1, axis="x", units="1")

    fnm = trim(out_dir)//"/geo.nc"

    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm, "nbr",x=[(i,i=1,n_coast_cells)], axis="e", ncid=ncid)
    call nc_write_dim(fnm, dim_lat, x=geo%grid%G%y, axis="y", units="degrees_north", ncid=ncid)
    call nc_write_dim(fnm, dim_lon, x=geo%grid%G%x, axis="x", units="degrees_east", ncid=ncid)
    call nc_close(ncid)

    if (l_output_hires) then
      fnm = trim(out_dir)//"/geo_hires.nc"
      call nc_create(fnm)
      call nc_open(fnm,ncid)
      call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
      unlimited=.TRUE., ncid=ncid)
      call nc_write_dim(fnm, dim_lat, x=geo%hires%grid%G%y, axis="y", units="degrees_north", ncid=ncid)
      call nc_write_dim(fnm, dim_lon, x=geo%hires%grid%G%x, axis="x", units="degrees_east", ncid=ncid)
      call nc_close(ncid)
    endif

   return

  end subroutine geo_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  g e o _ d i a g
  !   Purpose    :  geography diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_diag(geo)

    implicit none

    type(geo_class), intent(inout) :: geo

    real(wp), dimension(:,:), allocatable :: vol
    real(wp), dimension(:,:), allocatable :: vol_min
    real(wp), dimension(:,:), allocatable :: vol_pot
    real(wp), dimension(:,:), allocatable :: dh_p_e
    real(wp), dimension(:,:), allocatable :: dh_run
    integer, dimension(:,:), allocatable :: i_runoff
    integer, dimension(:,:), allocatable :: j_runoff
    real(wp), dimension(:,:), allocatable :: runoff
    integer :: i, j, n, y
    integer :: ni, nj
    integer :: n1out, ncid
    character (len=256) :: fnm
    real(wp) :: fac
    logical, save :: firstcall=.true.


    fac = 1.e-12_wp ! m2 -> mln km2

    if (mod(year,10).eq.1) then
      print '(a7,a9,6a7)','geo','year','seal','Aocn','Aveg','Aice','Alake','Abrng'
    endif
    print '(a7,i9,6F7.1)','geo',year_now,geo%sea_level,geo%ocn_area_tot*fac,geo%veg_area_tot*fac,geo%ice_area_tot*fac,geo%lake_area_tot*fac,geo%A_bering*1.e-6


    ! supress inital write to netCDF time series
    if(.not.firstcall) then
      ! current index
      y = y_out_ts_geo

      !print '(a8,2i9)','geotest ', ny_out_ts, y

      ann_ts(y)%sea_level = geo%sea_level
      ann_ts(y)%A_bering  = geo%A_bering
      ann_ts(y)%Aocn = geo%ocn_area_tot*fac
      ann_ts(y)%Aveg = geo%veg_area_tot*fac
      ann_ts(y)%Aice = geo%ice_area_tot*fac
      ann_ts(y)%Alake = geo%lake_area_tot*fac

      ! write to netcdf file 
      if (time_out_ts_geo) then
        n1out = year_geo-y+1    
        fnm = trim(out_dir)//"/geo_ts.nc"
        call nc_open(fnm,ncid)

        call nc_write(fnm,"time", dble([(i,i=(year_now-(y-1)*n_year_geo),(year_now),n_year_geo)]), &
        dim1=dim_time,start=[n1out],count=[y],ncid=ncid)    
        call nc_write(fnm,"sea_level", ann_ts(1:y)%sea_level, dim1=dim_time, start=[n1out], &
        count=[y],long_name="Global mean relative sea level",units="m",ncid=ncid)
        call nc_write(fnm,"A_bering", ann_ts(1:y)%A_bering, dim1=dim_time, start=[n1out], &
        count=[y],long_name="Bering Strait cross-sectional area",units="m2",ncid=ncid)
        call nc_write(fnm,"Aocn", ann_ts(1:y)%Aocn, dim1=dim_time, start=[n1out], &
        count=[y],long_name="Total ocean surface area",units="mln km2",ncid=ncid)
        call nc_write(fnm,"Aveg", ann_ts(1:y)%Aveg, dim1=dim_time, start=[n1out], &
        count=[y],long_name="Total ice-free lake-free land  surface area", &
        units="mln km2",ncid=ncid)
        call nc_write(fnm,"Aice", ann_ts(1:y)%Aice, dim1=dim_time, start=[n1out], &
        count=[y],long_name="Total ice sheet surface area",units="mln km2",ncid=ncid)
        call nc_write(fnm,"Alake", ann_ts(1:y)%Alake, dim1=dim_time, start=[n1out], &
        count=[y],long_name="Total lake surface area",units="mln km2",ncid=ncid)
        call nc_close(ncid)
      endif
    end if

    ! write 2D output
    if (time_out_geo .or. firstcall) then

      nout = nout +1 

      fnm = trim(out_dir)//"/geo.nc"
      call nc_open(fnm,ncid)
      if (firstcall) then
        call nc_write(fnm,dim_time,dble(year_ini), dim1=dim_time, start=[nout], count=[1],ncid=ncid)
      else
        call nc_write(fnm,dim_time,dble(year_now), dim1=dim_time, start=[nout], count=[1],ncid=ncid)
      endif
      call nc_write(fnm,"f_ocn", geo%f_ocn,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="ocean fraction, including floating ice",units="/",ncid=ncid)
      call nc_write(fnm,"f_ocn2", geo%f_ocn2,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="ocean fraction, excluding floating ice",units="/",ncid=ncid)
      call nc_write(fnm,"f_lnd", geo%f_lnd,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="land fraction",units="/",ncid=ncid)
      call nc_write(fnm,"f_ice", geo%f_ice,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="ice fraction",units="/",ncid=ncid)
      call nc_write(fnm,"f_ice_grd", geo%f_ice_grd,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="grounded ice fraction",units="/",ncid=ncid)
      call nc_write(fnm,"f_ice_flt", geo%f_ice_flt,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="floating ice fraction",units="/",ncid=ncid)
      call nc_write(fnm,"f_lake", geo%f_lake,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="lake fraction",units="/",ncid=ncid)
      call nc_write(fnm,"z_bed", geo%z_bed,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="mean bedrock elevation",units="m",ncid=ncid)
      call nc_write(fnm,"z_sur", geo%z_sur,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="grid cell mean surface elevation",units="m",ncid=ncid)
      call nc_write(fnm,"z_ocn", geo%z_ocn,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="mean elevation of ocean",units="m",ncid=ncid)
      call nc_write(fnm,"z_ocn_min", geo%z_ocn_min,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="min elevation of ocean",units="m",ncid=ncid)
      call nc_write(fnm,"z_ocn_max", geo%z_ocn_max,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="max elevation of ocean",units="m",ncid=ncid)
      call nc_write(fnm,"z_ocn_max_q", geo%z_ocn_max_q,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="quantile elevation of ocean",units="m",ncid=ncid)
      call nc_write(fnm,"z_ice", geo%z_ice,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="mean ice elevation",units="m",ncid=ncid)
      call nc_write(fnm,"z_lake", geo%z_lake,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="mean lake elevation",units="m",ncid=ncid)
      call nc_write(fnm,"z_veg", geo%z_veg,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="mean elevation of ice-free land",units="m",ncid=ncid)
      call nc_write(fnm,"z_veg_min", geo%z_veg_min,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="min elevation of ice-free land",units="m",ncid=ncid)
      call nc_write(fnm,"z_veg_max", geo%z_veg_max,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="max elevation of ice-free land",units="m",ncid=ncid)
      call nc_write(fnm,"z_sur_std", geo%z_sur_std,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="standard deviation of surface elevation",units="m",ncid=ncid)
      call nc_write(fnm,"z_sur_smooth_std", geo%z_sur_smooth_std,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="standard deviation of surface elevation",units="m",ncid=ncid)
      call nc_write(fnm,"z_veg_std", geo%z_veg_std,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="standard deviation of surface elevation, ice-free land only",units="m",ncid=ncid)
      call nc_write(fnm,"z_sur_lnd_std", geo%z_sur_lnd_std,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="standard deviation of surface elevation, ice-free land only",units="m",ncid=ncid)
      call nc_write(fnm,"i_runoff", geo%i_runoff,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="i index of runoff destination",units="/",ncid=ncid)
      call nc_write(fnm,"j_runoff", geo%j_runoff,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="j index of runoff destination",units="/",ncid=ncid)
      call nc_write(fnm,"i_runoff_veg", geo%i_runoff_veg,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="i index of runoff destination from ice-free land",units="/",ncid=ncid)
      call nc_write(fnm,"j_runoff_veg", geo%j_runoff_veg,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="j index of runoff destination from ice-free land",units="/",ncid=ncid)
      call nc_write(fnm,"i_runoff_ice", geo%i_runoff_ice,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="i index of runoff destination from ice-sheet",units="/",ncid=ncid)
      call nc_write(fnm,"j_runoff_ice", geo%j_runoff_ice,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="j index of runoff destination from ice-sheet",units="/",ncid=ncid)
      call nc_write(fnm,"f_drain_veg", geo%f_drain_veg(0,:,:),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="drainage fraction to ocean from ice-free land",units="/",ncid=ncid)
      call nc_write(fnm,"f_drain_ice", geo%f_drain_ice(0,:,:),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="drainage fration to ocean from ice sheets",units="/",ncid=ncid)
      call nc_write(fnm,"drain_basin", geo%drain_basins_ocn,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="drainage basins to ocean",units="/",ncid=ncid)
      call nc_write(fnm,"idivide_pac_atl", geo%idivide_pac_atl,dims=[dim_lat,dim_time],start=[1,nout],count=[geo%grid%G%ny,1],long_name="i-index of continental divide between Pacific and Atlantic",units="/",ncid=ncid)
      call nc_write(fnm,"idivide_atl_indpac", geo%idivide_atl_indpac,dims=[dim_lat,dim_time],start=[1,nout],count=[geo%grid%G%ny,1],long_name="i-index of continental divide between Atlantic and Pacific",units="/",ncid=ncid)
      call nc_write(fnm,"mask_coast", geo%mask_coast,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="coast mask",units="/",ncid=ncid)
      call nc_write(fnm,"mask_coast2", geo%mask_coast2,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="coast mask",units="/",ncid=ncid)
      call nc_write(fnm,"i_coast_nbr", geo%i_coast_nbr,dims=[dim_lon,dim_lat,"nbr",dim_time],start=[1,1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,n_coast_cells,1],long_name="i index of neighbors of coastal cells",units="/",ncid=ncid)
      call nc_write(fnm,"j_coast_nbr", geo%j_coast_nbr,dims=[dim_lon,dim_lat,"nbr",dim_time],start=[1,1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,n_coast_cells,1],long_name="j index of neighbors of coastal cells",units="/",ncid=ncid)
      call nc_write(fnm,"q_geo", geo%q_geo,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[geo%grid%G%nx,geo%grid%G%ny,1],long_name="geothermal heat flux",units="W/m2",ncid=ncid)
      call nc_close(ncid)

      if (l_output_hires) then

        ni = geo%hires%grid%G%nx
        nj = geo%hires%grid%G%ny

        ! write to file
        fnm = trim(out_dir)//"/geo_hires.nc"
        call nc_open(fnm,ncid)
        if (firstcall) then
          call nc_write(fnm,dim_time,dble(year_ini), dim1=dim_time, start=[nout], count=[1],ncid=ncid)
        else
          call nc_write(fnm,dim_time,dble(year_now), dim1=dim_time, start=[nout], count=[1],ncid=ncid)
        endif
        call nc_write(fnm,"z_topo", sngl(geo%hires%z_topo),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="topography",units="m",ncid=ncid)
        call nc_write(fnm,"z_topo_fill", sngl(geo%hires%z_topo_fill),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="filled topography",units="m",ncid=ncid)
        call nc_write(fnm,"z_topo_fil", sngl(geo%hires%z_topo_fil),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="filtered topography",units="m",ncid=ncid)
        call nc_write(fnm,"z_bed", sngl(geo%hires%z_bed),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="bedrock elevation",units="m",ncid=ncid)
        call nc_write(fnm,"dz_bed", sngl(geo%hires%z_bed-geo%hires%z_bed_rel),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="difference between bedrock elevation and relaxed bedrock elevation",units="m",ncid=ncid)
        if (firstcall) then
          call nc_write(fnm,"z_bed_ref", sngl(geo%hires%z_bed_ref),dims=[dim_lon,dim_lat],start=[1,1],count=[ni,nj],long_name="present day reference bedrock elevation",units="m",ncid=ncid)
          if (i_geo.eq.1) then
            call nc_write(fnm,"z_bed_rel", sngl(geo%hires%z_bed_rel),dims=[dim_lon,dim_lat],start=[1,1],count=[ni,nj],long_name="relaxed ice-free bedrock elevation",units="m",ncid=ncid)
          endif
        endif
        call nc_write(fnm,"z_sur", sngl(geo%hires%z_sur),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="surface elevation",units="m",ncid=ncid)
        call nc_write(fnm,"h_ice", sngl(geo%hires%h_ice),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="ice sheet thickness",units="m",ncid=ncid)
        call nc_write(fnm,"mask", geo%hires%mask,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="mask",units="0=ice, 1=land, 2=ocean, 3=floating ice, 4=lake",ncid=ncid)
        call nc_write(fnm,"rsl", sngl(geo%hires%rsl),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="relative sea level",units="m",ncid=ncid)
        call nc_write(fnm,"mask_lake", geo%hires%mask_lake,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="lakes mask",units="/",ncid=ncid)
        call nc_write(fnm,"mask_lake_pot", geo%hires%mask_lake_pot,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="potential lakes mask",units="/",ncid=ncid)
        call nc_write(fnm,"z_lake", sngl(geo%hires%z_lake),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="lake surface elevation",units="m",ncid=ncid)
        call nc_write(fnm,"z_sea_lake", sngl(geo%hires%z_sea_lake),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="sea/lake level",units="m",ncid=ncid)
        call nc_write(fnm,"i_runoff", geo%hires%i_runoff,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="i-index of runoff destination cell",units="/",ncid=ncid)
        call nc_write(fnm,"j_runoff", geo%hires%j_runoff,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="j-index of runoff destination cell",units="/",ncid=ncid)
        call nc_write(fnm,"map_runoff", geo%hires%map_runoff,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="runoff mapping to ocean/lakes",units="/",ncid=ncid)
        call nc_write(fnm,"flow_acc", sngl(geo%hires%flow_acc)*1e-12,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="river flow accumulation, upstream drainage area",units="mln km2",ncid=ncid)
        call nc_write(fnm,"drain_basin", geo%hires%drain_basins_ocn,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="drainage basins to ocean",units="/",ncid=ncid)

        if (flag_lakes) then

          allocate(vol(ni,nj))
          allocate(vol_min(ni,nj))
          allocate(vol_pot(ni,nj))
          allocate(dh_p_e(ni,nj))
          allocate(dh_run(ni,nj))
          allocate(i_runoff(ni,nj))
          allocate(j_runoff(ni,nj))
          allocate(runoff(ni,nj))

          do i=1,ni
            do j=1,nj
              n = geo%hires%mask_lake(i,j)
              if (n>0) then
                vol(i,j)    = geo%lake(n)%vol
                vol_min(i,j)= geo%lake(n)%vol_min
                vol_pot(i,j)= geo%lake(n)%vol_pot
                dh_p_e(i,j) = geo%lake(n)%dh_p_e
                dh_run(i,j) = geo%lake(n)%dh_run
                i_runoff(i,j) = geo%lake(n)%hires%i_runoff
                j_runoff(i,j) = geo%lake(n)%hires%j_runoff
                runoff(i,j) = geo%lake(n)%runoff_diag
              else
                vol(i,j)    = 0._wp
                vol_min(i,j)= 0._wp 
                vol_pot(i,j)= 0._wp 
                dh_p_e(i,j) = 0._wp
                dh_run(i,j) = 0._wp
                i_runoff(i,j) = 0._wp
                j_runoff(i,j) = 0._wp
                runoff(i,j) = 0._wp
              endif
            enddo
          enddo

          call nc_write(fnm,"lake_vol", sngl(vol),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="lake volume",units="m3",ncid=ncid)
          call nc_write(fnm,"lake_vol_min", sngl(vol_min),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="minimum lake volume",units="m3",ncid=ncid)
          call nc_write(fnm,"lake_vol_pot", sngl(vol_pot),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="potential lake volume",units="m3",ncid=ncid)
          call nc_write(fnm,"lake_dh_p_e", sngl(dh_p_e),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="lake depth tendency due to surface net water balance",units="m/a",ncid=ncid)
          call nc_write(fnm,"lake_dh_run", sngl(dh_run),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="lake depth tendency due to runoff input",units="m/a",ncid=ncid)
          call nc_write(fnm,"lake_dh", sngl(dh_run+dh_p_e),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="lake depth tendency",units="m/a",ncid=ncid)
          call nc_write(fnm,"lake_runoff", sngl(runoff),dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="runoff from lake to the ocean",units="Sv",ncid=ncid)
          call nc_write(fnm,"lake_i_runoff", i_runoff,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="i-index of lake runoff destination cell",units="/",ncid=ncid)
          call nc_write(fnm,"lake_j_runoff", j_runoff,dims=[dim_lon,dim_lat,dim_time],start=[1,1,nout],count=[ni,nj,1],long_name="j-index of lake runoff destination cell",units="/",ncid=ncid)

          deallocate(vol)
          deallocate(vol_min)
          deallocate(vol_pot)
          deallocate(dh_p_e)
          deallocate(dh_run)
          deallocate(i_runoff)
          deallocate(j_runoff)
          deallocate(runoff)

        endif

        call nc_close(ncid)

      endif

    endif

    firstcall = .false.

    return

  end subroutine geo_diag

end module geo_out

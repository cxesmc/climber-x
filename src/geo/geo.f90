!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : g e o
!
!  Purpose : main geography model 
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
module geo_mod

  use precision, only : wp, dp
  use control, only : in_dir, sea_level_init, flag_geo, flag_lakes, flag_ice, flag_smb, flag_bgc
  use control, only : geo_restart, restart_in_dir
  use constants, only : rho_i, rho_sw, pi, R_earth
  use dim_name, only: dim_lon, dim_lat, dim_depth, dim_lon1, dim_lat1, dim_time
  use coord, only : grid_class, grid_init
  use coord, only : map_scrip_init, map_scrip_class, map_scrip_field
  use ncio
  use geo_params, only : i_geo, geo_ref_file
  use geo_params, only : z_bed_rel_file, z_bed_1min_file, sigma_filter
  use geo_params, only : geo_params_init, l_connect_ocn, f_crit, n_coast_cells, i_fix_cell_grl
  use geo_params, only : h_ice_min
  use geo_params, only : l_close_hudson, l_close_baltic
  use geo_params, only : lon_ocn_origin, lat_ocn_origin
  use geo_params, only : l_write_timer
  use geo_grid, only : geo_grid_init
  use geo_grid, only : ni, nj, ni_topo, nj_topo, n_topo_sur, lon ,lat, lon_topo, lat_topo, area, area_dp
  use geo_grid, only : i0_topo, i1_topo, j0_topo, j1_topo, i_lowres, j_lowres
  use geo_def, only : geo_class
  use runoff_routing_mod, only : runoff_routing, runoff_routing_lakes
  use lake_mod, only : lake_pot, update_lake_mask, lake_type, ocean_lake_mask_level
  use fill_ocean_mod, only : fill_ocean
  use topo_fill_mod, only : topo_fill
  use topo_filter_mod, only : topo_filter
  use vilma_model, only : vilma_init, vilma_update, vilma_end, vilma_write_restart
  use gia_mod, only : gia_init, gia_update
  use q_geo_mod, only : geo_heat
  use sed_mod, only : sed
  use connect_ocn_mod, only : connect_ocn
  use fix_runoff_mod, only : fix_runoff, fix_runoff_lake
  use hires_to_lowres_mod, only : hires_to_lowres
  use coast_cells_mod, only : coast_cells
  use drainage_basins_mod, only : drainage_basins
  use corals_topo_mod, only : corals_topo
  !$ use omp_lib

  implicit none

   private
   public :: geo_init, geo_update, geo_end
   public :: geo_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g e o _ i n i t
  ! Purpose  :  initialize geo
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_init(geo, bnd_geo_grid, z_bed, z_bed_ref, bnd_ice_grid, h_ice)

    implicit none

    type(geo_class) :: geo 
    type(grid_class), intent(in) :: bnd_geo_grid
    real(wp), dimension(:,:), intent(in) :: z_bed
    real(wp), dimension(:,:), intent(in) :: z_bed_ref
    type(grid_class), intent(in) :: bnd_ice_grid
    real(wp), dimension(:,:), intent(in) :: h_ice

    integer :: i, j, nocn
    integer :: ni_rel, nj_rel
    integer :: ni_1min, nj_1min
    integer :: ppos
    real(wp) :: dlon_rel, dlat_rel
    character(len=256) :: fnm
    logical :: l_coral_exist
    real(wp), dimension(:), allocatable :: lon_rel
    real(wp), dimension(:), allocatable :: lat_rel
    integer, dimension(:), allocatable :: mask_cell
    real(wp), dimension(:,:), allocatable :: tmp
    real(wp), dimension(:,:), allocatable :: z_bed_anom
    real(wp), allocatable, dimension(:,:) :: mask_ice
    real(wp), allocatable, dimension(:,:) :: mask_ice_geo
    type(grid_class) :: grid_rel
    type(map_scrip_class) :: maps_hires_to_lowres
    type(map_scrip_class) :: maps_rel_to_geo
    type(map_scrip_class) :: maps_bndgeo_to_geo
    type(map_scrip_class) :: maps_bndice_to_geo


    ! read geo parameters
    call geo_params_init

    ! initialize grids
    call geo_grid_init(geo%grid, geo%hires%grid)

    allocate( geo%hires%z_bed_ref(ni_topo,nj_topo))
    allocate( geo%hires%h_ice_ref(ni_topo,nj_topo))
    allocate( geo%hires%z_topo_ref(ni_topo,nj_topo))
    allocate( geo%hires%mask_ref(ni_topo,nj_topo))
    allocate( geo%hires%z_bed_rel(ni_topo,nj_topo))
    allocate( geo%hires%z_bed(ni_topo,nj_topo))
    allocate( geo%hires%z_topo(ni_topo,nj_topo))
    allocate( geo%hires%z_topo_fil(ni_topo,nj_topo))
    allocate( geo%hires%z_topo_fill(ni_topo,nj_topo))
    allocate( geo%hires%z_sur(ni_topo,nj_topo))
    allocate( geo%hires%h_ice(ni_topo,nj_topo))
    allocate( geo%hires%rsl(ni_topo,nj_topo))
    allocate( geo%hires%mask(ni_topo,nj_topo))
    allocate( geo%hires%q_geo(ni_topo,nj_topo))
    allocate( geo%hires%h_sed(ni_topo,nj_topo))
    allocate( geo%hires%i_runoff(ni_topo,nj_topo))
    allocate( geo%hires%j_runoff(ni_topo,nj_topo))
    allocate( geo%hires%i_runoff_coarse(ni_topo,nj_topo))
    allocate( geo%hires%j_runoff_coarse(ni_topo,nj_topo))
    allocate( geo%hires%map_runoff(ni_topo,nj_topo))
    allocate( geo%hires%flow_acc(ni_topo,nj_topo))
    allocate( geo%hires%mask_lake_pot(ni_topo,nj_topo))
    geo%hires%mask_lake_pot = 0
    allocate( geo%hires%mask_lake(ni_topo,nj_topo))
    allocate( geo%hires%z_lake(ni_topo,nj_topo))
    allocate( geo%hires%z_sea_lake(ni_topo,nj_topo))
    allocate( geo%hires%mask_ocn_lake(ni_topo,nj_topo))
    allocate( geo%hires%drain_basins_ocn(ni_topo,nj_topo))

    allocate(geo%f_ocn0(ni,nj))
    allocate(geo%f_lnd0(ni,nj))
    allocate(geo%f_ocn(ni,nj))
    allocate(geo%f_ocn2(ni,nj))
    allocate(geo%f_lnd(ni,nj))
    allocate(geo%f_lake(ni,nj))
    allocate(geo%f_ice(ni,nj))
    allocate(geo%f_ice_grd(ni,nj))
    allocate(geo%f_ice_flt(ni,nj))
    allocate(geo%z_sur(ni,nj))
    allocate(geo%z_sur_std(ni,nj))
    allocate(geo%z_sur_smooth_std(ni,nj))
    allocate(geo%z_ocn(ni,nj))
    allocate(geo%z_ocn_min(ni,nj))
    allocate(geo%z_ocn_max(ni,nj))
    allocate(geo%z_ocn_max_q(ni,nj))
    allocate(geo%z_veg(ni,nj))
    allocate(geo%z_veg_min(ni,nj))
    allocate(geo%z_veg_max(ni,nj))
    allocate(geo%z_veg_std(ni,nj))
    allocate(geo%z_sur_lnd_std(ni,nj))
    allocate(geo%z_ice(ni,nj))
    allocate(geo%z_lake(ni,nj))
    allocate(geo%z_bed(ni,nj))
    allocate(geo%mask_coast(ni,nj))
    allocate(geo%mask_coast2(ni,nj))
    allocate(geo%i_coast_nbr(ni,nj,n_coast_cells))
    allocate(geo%j_coast_nbr(ni,nj,n_coast_cells))
    allocate(geo%coast_nbr(ni,nj))
    allocate(geo%i_runoff(ni,nj))
    allocate(geo%j_runoff(ni,nj))
    allocate(geo%i_runoff_veg(ni,nj))
    allocate(geo%j_runoff_veg(ni,nj))
    allocate(geo%i_runoff_ice(ni,nj))
    allocate(geo%j_runoff_ice(ni,nj))
    allocate(geo%drain_basins_ocn(ni,nj))
    allocate(geo%idivide_pac_atl(nj))
    allocate(geo%idivide_atl_indpac(nj))
    allocate(geo%coral_f_area(ni,nj,-250:50))
    allocate(geo%coral_f_topo(ni,nj,-250:50))
    allocate(geo%q_geo(ni,nj))

    !-------------------------------------------------------------------
    ! reference present day bedrock topography
    !-------------------------------------------------------------------

    call nc_read(trim(geo_ref_file),"bedrock_topography",geo%hires%z_bed_ref)
    call nc_read(trim(geo_ref_file),"ice_thickness",geo%hires%h_ice_ref)

    !------------------------------------------------------------------------
    ! restart 
    !------------------------------------------------------------------------
    if (geo_restart) then

      !------------------------------------------------------------------------
      ! read restart file
      !------------------------------------------------------------------------
      call geo_read_restart(trim(restart_in_dir)//"/geo_restart.nc",geo)

    else  

      !-------------------------------------------------------------------
      ! set initial topography and mask 
      !-------------------------------------------------------------------

      ! initialize bedrock topography, map from bnd to geo
      if (bnd_geo_grid%name==geo%hires%grid%name) then
        geo%hires%z_bed = z_bed
      else
        allocate(z_bed_anom(geo%hires%grid%G%nx,geo%hires%grid%G%ny))
        call map_scrip_init(maps_bndgeo_to_geo,bnd_geo_grid,geo%hires%grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
        call map_scrip_field(maps_bndgeo_to_geo,"z_bed",z_bed-z_bed_ref,z_bed_anom,method="mean",missing_value=-9999._dp)
        geo%hires%z_bed = geo%hires%z_bed_ref + z_bed_anom 
        deallocate(z_bed_anom)
      endif

      ! option to close Hudson Bay
      if (l_close_hudson) then
        where (geo%hires%grid%lon<-69._wp .and. geo%hires%grid%lon>-71._wp .and. geo%hires%grid%lat>60._wp .and. geo%hires%grid%lat<65._wp)
          geo%hires%z_bed = max(0.1_wp, geo%hires%z_bed)
        endwhere
      endif
      ! option to close Baltic Sea
      if (l_close_baltic) then
        where (geo%hires%grid%lon>11._wp .and. geo%hires%grid%lon<13._wp .and. geo%hires%grid%lat>53._wp .and. geo%hires%grid%lat<56._wp)
          geo%hires%z_bed = max(0.1_wp, geo%hires%z_bed)
        endwhere
      endif

      ! initialize ice sheet thickness with reference field
      geo%hires%h_ice = geo%hires%h_ice_ref
      ! map from bnd to geo
      if (bnd_ice_grid%name==geo%hires%grid%name) then
        geo%hires%h_ice = h_ice
      else
        call map_scrip_init(maps_bndice_to_geo,bnd_ice_grid,geo%hires%grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)
        call map_scrip_field(maps_bndice_to_geo,"h_ice",h_ice,geo%hires%h_ice,method="mean",missing_value=-9999._dp,reset=.false.)
        allocate(mask_ice(bnd_ice_grid%G%nx,bnd_ice_grid%G%ny))
        allocate(mask_ice_geo(geo%hires%grid%G%nx,geo%hires%grid%G%ny))
        where (h_ice>h_ice_min) 
          mask_ice = 1.
        elsewhere
          mask_ice = 0.
        endwhere
        mask_ice_geo = 1.
        call map_scrip_field(maps_bndice_to_geo,"mask",mask_ice,mask_ice_geo,method="mean",missing_value=-9999._dp,reset=.false., &
          filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
        where (mask_ice_geo<0.5) geo%hires%h_ice = 0._wp
        deallocate(mask_ice)
        deallocate(mask_ice_geo)
      endif

      ! derive mask and topography
      ! initially land everywhere
      geo%hires%mask = 1 
      ! preliminary set all ice to grounded
      where (geo%hires%h_ice.gt.0._wp) 
        geo%hires%mask = 0
      endwhere
      ! preliminary new topography
      geo%hires%z_topo = geo%hires%z_bed + geo%hires%h_ice
      ! floating ice
      where (geo%hires%h_ice.gt.0._wp .and. geo%hires%h_ice.lt.(-geo%hires%z_bed*rho_sw/rho_i) .and. geo%hires%z_bed.lt.0._wp)
        geo%hires%mask = 3
        geo%hires%z_topo = geo%hires%h_ice*(1._wp-rho_i/rho_sw)
      endwhere
      ! fill ocean domain starting from deep ocean 'origin' point
      call fill_ocean(geo%hires%z_topo,lon_topo,lat_topo,lon_ocn_origin,lat_ocn_origin, &
        geo%hires%mask) 

      ! set to lake where land but topography below sea level
      if (.not.flag_lakes) then
        where (geo%hires%mask.eq.1 .and. geo%hires%z_topo.le.0._wp)
          geo%hires%mask = 4  ! lake
        endwhere
      endif

      !-------------------------------------------------------------------
      ! fill topography depressions 
      !-------------------------------------------------------------------
      call topo_fill(geo%hires%z_topo, geo%hires%mask, &
        geo%hires%z_topo_fill)

      !-------------------------------------------------------------------
      ! initialize lakes
      !-------------------------------------------------------------------
      if (flag_lakes) then

        ! determine potential lakes
        call lake_pot(geo%hires%z_topo, geo%hires%z_topo_fill, geo%hires%mask, area, & 
          geo%hires%mask_lake_pot, geo%lake)
        geo%n_lakes = size(geo%lake)

        ! compute actual lake mask and elevation
        call update_lake_mask(geo%hires%z_topo, area, geo%lake, geo%hires%mask_lake_pot, &
          geo%hires%mask_lake, geo%hires%z_lake)

        ! exclude lake points over floating ice, assume ice is grounded instead
        where (geo%hires%mask_lake.ne.0 .and. geo%hires%mask.eq.3) 
          geo%hires%mask_lake = 0
          geo%hires%mask = 0
        endwhere
        ! add lakes to mask
        where (geo%hires%mask_lake.ne.0)
          geo%hires%mask = 4  ! lake
        endwhere

      else
        ! no lakes

        geo%n_lakes = 0                            
        geo%hires%z_lake = 0._wp

      endif

      ! derive ocean/lake mask
      call ocean_lake_mask_level(geo%hires%mask, geo%hires%z_lake, &
        geo%hires%mask_ocn_lake, geo%hires%z_sea_lake)

      !-------------------------------------------------------------------
      ! initialize surface elevation
      !-------------------------------------------------------------------
      geo%hires%z_sur = geo%hires%z_topo
      where (geo%hires%mask.eq.2) geo%hires%z_sur = 0._wp   ! over ocean
      where (geo%hires%mask.eq.4) geo%hires%z_sur = geo%hires%z_lake   ! over lake
      where (geo%hires%z_sur.lt.0._wp) geo%hires%z_sur = 0._wp  ! negative surface elevation not allowed, fixme? could be negative over lakes...

      !-------------------------------------------------------------------
      ! initialize runoff directions of filled topography
      !-------------------------------------------------------------------
      call runoff_routing(geo%hires%z_topo_fill, geo%hires%mask, real(geo%hires%grid%area,wp), &
        geo%hires%i_runoff, geo%hires%j_runoff, geo%hires%flow_acc)

      !-------------------------------------------------------------------
      ! map runoff to lakes and route lake runoff to ocean
      !-------------------------------------------------------------------
      if (flag_lakes) then
        call runoff_routing_lakes(geo%hires%z_topo, geo%hires%z_topo_fill, geo%hires%mask, geo%hires%mask_lake_pot, &
          geo%hires%i_runoff, geo%hires%j_runoff, &
          geo%hires%map_runoff, geo%lake)
      else
        geo%hires%map_runoff = 0    ! all runoff directly rooted to ocean
      endif
                   
      !-------------------------------------------------------------------
      ! initialize sea level
      !-------------------------------------------------------------------
      if (flag_geo) then
        geo%sea_level = sea_level_init
      endif

    endif

    geo%ocn_area_tot = sum(area_dp,geo%hires%mask.eq.2 .or. geo%hires%mask.eq.3)

    !-------------------------------------------------------------------
    ! reference present day topography and fractions
    !-------------------------------------------------------------------

    ! derive mask and topography
    ! initially land everywhere
    geo%hires%mask_ref = 1 
    ! preliminary set all ice to grounded
    where (geo%hires%h_ice_ref.gt.0._wp) 
      geo%hires%mask_ref = 0
    endwhere
    ! preliminary new topography
    geo%hires%z_topo_ref = geo%hires%z_bed_ref + geo%hires%h_ice_ref
    ! floating ice
    where (geo%hires%h_ice_ref.gt.0._wp .and. geo%hires%h_ice_ref.lt.(-geo%hires%z_bed_ref*rho_sw/rho_i) .and. geo%hires%z_bed_ref.lt.0._wp)
      geo%hires%mask_ref = 3
      geo%hires%z_topo_ref = geo%hires%h_ice_ref*(1._wp-rho_i/rho_sw)
    endwhere
    ! fill ocean domain starting from deep ocean 'origin' point
    call fill_ocean(geo%hires%z_topo_ref,lon_topo,lat_topo,lon_ocn_origin,lat_ocn_origin, &
      geo%hires%mask_ref) 

    ! set to lake where land but topography below sea level
    where (geo%hires%mask_ref.eq.1 .and. geo%hires%z_topo_ref.le.0._wp)
      geo%hires%mask_ref = 4  ! lake
    endwhere

    ! compute land and ocean fraction for present day sea level

    ! first set all ocean to be land
    where (geo%hires%mask_ref==2) geo%hires%mask_ref = 1

    ! fill ocean domain starting from deep ocean 'origin' point
    call fill_ocean(geo%hires%z_topo_ref,lon_topo,lat_topo,lon_ocn_origin,lat_ocn_origin, &
      geo%hires%mask_ref) 

    allocate( mask_cell(n_topo_sur) )
    do j=1,nj
      do i=1,ni
        mask_cell = pack(geo%hires%mask_ref(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
        ! ocean fraction
        ! count number of ocean points in grid cell
        nocn = count(mask_cell.eq.2) 
        geo%f_ocn0(i,j) = real(nocn,wp)/real(n_topo_sur,wp)
        ! apply critical threshold for ocean/land
        if (geo%f_ocn0(i,j).lt.f_crit) then
          geo%f_ocn0(i,j) = 0._wp
        endif
        if (geo%f_ocn0(i,j).ge.(1._wp-f_crit)) then
          geo%f_ocn0(i,j) = 1._wp
        endif
      enddo
    enddo
    deallocate(mask_cell)

    ! close Panama, ocean fraction=0
    geo%f_ocn0(16,23) = 0._wp
    geo%f_ocn0(17:18,22) = 0._wp
    geo%f_ocn0(19:20,21) = 0._wp
    geo%f_ocn0(21,20) = 0._wp

    ! land fraction
    geo%f_lnd0 = 1._wp - geo%f_ocn0

    !-------------------------------------------------------------------
    ! relaxed icefree bedrock topography 
    !-------------------------------------------------------------------

    if (flag_ice .or. i_geo.eq.1) then

      ni_rel = nc_size(trim(z_bed_rel_file),"lon")
      nj_rel = nc_size(trim(z_bed_rel_file),"lat")
      allocate( lon_rel(ni_rel) )
      allocate( lat_rel(nj_rel) )
      call nc_read(trim(z_bed_rel_file),"lon",lon_rel)
      call nc_read(trim(z_bed_rel_file),"lat",lat_rel)
      dlon_rel = lon_rel(2)-lon_rel(1)
      dlat_rel = lat_rel(2)-lat_rel(1)
      if (ni_rel.eq.ni_topo .and. nj_rel.eq.nj_topo) then
        call nc_read(trim(z_bed_rel_file),"z_bed_rel",geo%hires%z_bed_rel)
      else
        ! generate grid object
        ppos = scan(trim(z_bed_rel_file),".", BACK= .true.)-1
        call grid_init(grid_rel,name=trim(z_bed_rel_file(7:ppos)),mtype="latlon",units="degrees", &
          x0=real(lon_rel(1),dp),dx=real(dlon_rel,dp),nx=ni_rel,y0=real(lat_rel(1),dp),dy=real(dlat_rel,dp),ny=nj_rel)

        allocate( tmp(ni_rel,nj_rel) )
        call nc_read(trim(z_bed_rel_file),"z_bed_rel",tmp)
        call map_scrip_init(maps_rel_to_geo,grid_rel,geo%hires%grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)
        call map_scrip_field(maps_rel_to_geo,"z_bed_rel",tmp,geo%hires%z_bed_rel,method="mean",missing_value=-9999._dp, &
          filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])
        deallocate( tmp )

      endif
      deallocate( lon_rel )
      deallocate( lat_rel )

    endif

    !-------------------------------------------------------------------
    ! bedrock topography at 1 min resolution
    !-------------------------------------------------------------------

    if (flag_smb) then
      ni_1min = nc_size(trim(z_bed_1min_file),"lon")
      nj_1min = nc_size(trim(z_bed_1min_file),"lat")
      allocate( geo%hires%lon_1min(ni_1min) )
      allocate( geo%hires%lat_1min(nj_1min) )
      call nc_read(trim(z_bed_1min_file),"lon",geo%hires%lon_1min)
      call nc_read(trim(z_bed_1min_file),"lat",geo%hires%lat_1min)
      allocate( geo%hires%z_bed_1min(ni_1min,nj_1min))
      call nc_read(trim(z_bed_1min_file),"bedrock_topography",geo%hires%z_bed_1min)
    endif

    !-------------------------------------------------------------------
    ! topography factor for corals (Kleypas 1997)
    !-------------------------------------------------------------------

    if (flag_bgc) then

      fnm = trim(in_dir)//"corals_"//trim(geo_ref_file(7:))
      inquire( file=fnm, exist=l_coral_exist ) 
      if (.not.l_coral_exist) then  

        ! recompute only if coral file does not exist (has not been generated during previous runs)
        call corals_topo(geo%hires%z_topo_ref, &  ! in
          geo%coral_f_area, geo%coral_f_topo)    ! out

        ! write coral topographic variables to netcdf file
        print *,'create new coral file: ',trim(fnm)
        call nc_create(fnm)
        call nc_write_dim(fnm,dim_lon,x=lon,axis="x")
        call nc_write_dim(fnm,dim_lat,x=lat,axis="y")
        call nc_write_dim(fnm,dim_depth,x=-250._wp,dx=1._wp,nx=301,units="m",axis="z")
        call nc_write_dim(fnm,dim_lon1,x=lon_topo,axis="x")
        call nc_write_dim(fnm,dim_lat1,x=lat_topo,axis="y")
        call nc_write(fnm,"coral_f_area",geo%coral_f_area,  dims=[dim_lon,dim_lat,dim_depth],start=[1,1,1],count=[ni,nj,301])
        call nc_write(fnm,"coral_f_topo",geo%coral_f_topo,  dims=[dim_lon,dim_lat,dim_depth],start=[1,1,1],count=[ni,nj,301])

      else

        ! read from existing netcdf file
        print *,'read from existing coral file: ',trim(fnm)
        call nc_read(trim(fnm),"coral_f_area",geo%coral_f_area)
        call nc_read(trim(fnm),"coral_f_topo",geo%coral_f_topo)

      endif

    endif

    !-------------------------------------------------------------------
    ! geothermal heat flux
    !-------------------------------------------------------------------

    call geo_heat(geo%hires%grid, geo%hires%q_geo)
    call map_scrip_init(maps_hires_to_lowres,geo%hires%grid,geo%grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)
    call map_scrip_field(maps_hires_to_lowres,"q_geo",geo%hires%q_geo,geo%q_geo,method="mean",missing_value=-9999._dp, &
      filt_method="none",filt_par=[5._dp*geo%hires%grid%G%dx,geo%hires%grid%G%dx])

    !-------------------------------------------------------------------
    ! sediment thickness and mask
    !-------------------------------------------------------------------

    if (flag_ice) then
      call sed(geo%hires%grid, geo%hires%h_sed)
    endif

    !-------------------------------------------------------------------
    ! initialize VILMA
    !-------------------------------------------------------------------
    if (flag_geo .and. i_geo.eq.2) then
      call vilma_init(geo%hires%grid, geo%hires%z_bed_ref, geo%hires%h_ice_ref, geo%hires%h_ice) 
    endif

    print*
    print*,'======================================================='
    print*,' Initialisation of GEOGRAPHY complete'
    print*,'======================================================='
    print*

  return

  end subroutine geo_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g e o _ u p d a t e
  ! Purpose  :  update sea level, topography, geography and runoff directions
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_update(geo)

  implicit none

  type(geo_class), intent(inout) :: geo

  integer :: i, j, n
  integer :: nocn
  logical :: l_was_ice_grd
  real(wp) :: A_bering, dx
  logical, save :: firstcall = .true.

  !$ real(wp) :: time1,time2

  !-------------------------------------------------------------------
  ! update sea level and bedrock topography 
  !-------------------------------------------------------------------
  if (flag_geo) then

    !$ time1 = omp_get_wtime()
    if (i_geo==0) then
      !-------------------------------------------------------------------
      ! Uniform sea level from ice volume

      ! update sea level
      geo%sea_level = geo%sea_level + geo%d_sea_level
      geo%hires%rsl = geo%sea_level   ! just for output

    else if (i_geo==1) then
      !-------------------------------------------------------------------
      ! LLRA GIA model + uniform sea level from ice volume

      ! update sea level
      geo%sea_level = geo%sea_level + geo%d_sea_level
      geo%hires%rsl = geo%sea_level   ! just for output

     ! compute glacial isostatic adjustment and update bedrock elevation (relative to current sea level)
      call gia_update(geo%hires%z_bed_rel, geo%hires%h_ice, geo%hires%mask, geo%sea_level, geo%d_sea_level, &    ! in
        geo%hires%z_bed)  ! inout

    else if (i_geo==2) then
      !-------------------------------------------------------------------
      ! VILMA solid Earth model

      if (.not.firstcall) then
        ! takes current ice load (thickness) as input
        ! returns relative sea level and new bedrock elevation
        call vilma_update(geo%hires%grid, geo%hires%h_ice, & ! in
          geo%hires%rsl, geo%hires%z_bed)    ! out
      endif

    endif
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'vilma',time2-time1

  endif

  ! option to close Hudson Bay
  if (l_close_hudson) then
    where (geo%hires%grid%lon<-69._wp .and. geo%hires%grid%lon>-71._wp .and. geo%hires%grid%lat>60._wp .and. geo%hires%grid%lat<65._wp)
      geo%hires%z_bed = max(0.1_wp, geo%hires%z_bed)
    endwhere
  endif
  ! option to close Baltic Sea
  if (l_close_baltic) then
    where (geo%hires%grid%lon>11._wp .and. geo%hires%grid%lon<13._wp .and. geo%hires%grid%lat>53._wp .and. geo%hires%grid%lat<56._wp)
      geo%hires%z_bed = max(0.1_wp, geo%hires%z_bed)
    endwhere
  endif

  !-------------------------------------------------------------------
  ! derive ocean/land/ice mask and topography
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  !$omp parallel do private(i,j,l_was_ice_grd)
  do j=1,nj_topo
    do i=1,ni_topo
      l_was_ice_grd = geo%hires%mask(i,j)==0
      ! initially land everywhere
      geo%hires%mask(i,j) = 1 
      ! preliminary set all ice to grounded
      if (geo%hires%h_ice(i,j).gt.0._wp) then
        geo%hires%mask(i,j) = 0
      endif
      ! preliminary new topography
      geo%hires%z_topo(i,j) = geo%hires%z_bed(i,j) + geo%hires%h_ice(i,j)
      ! floating ice
      if (l_was_ice_grd) then
        ! ice was grounded, use lower threshold of ice thickness for floating
        ! conditions in order to avoid jumps between grounded and floating ice
        if (geo%hires%h_ice(i,j).gt.0._wp .and. geo%hires%h_ice(i,j).lt.(-geo%hires%z_bed(i,j)) .and. geo%hires%z_bed(i,j).lt.0._wp) then
          geo%hires%mask(i,j) = 3
          geo%hires%z_topo(i,j) = geo%hires%h_ice(i,j)*(1._wp-rho_i/rho_sw)
        endif
      else
        if (geo%hires%h_ice(i,j).gt.0._wp .and. geo%hires%h_ice(i,j).lt.(-geo%hires%z_bed(i,j)*rho_sw/rho_i) .and. geo%hires%z_bed(i,j).lt.0._wp) then
          geo%hires%mask(i,j) = 3
          geo%hires%z_topo(i,j) = geo%hires%h_ice(i,j)*(1._wp-rho_i/rho_sw)
        endif
      endif
    enddo
  enddo
  !$omp end parallel do
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'mask and topo update',time2-time1

  ! update land-ocean mask 
  ! fill ocean domain starting from deep ocean 'origin' point
  !$ time1 = omp_get_wtime()
  call fill_ocean(geo%hires%z_topo,lon_topo,lat_topo,lon_ocn_origin,lat_ocn_origin, &
    geo%hires%mask) 
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'fill_ocean',time2-time1

  ! set to lake where land but topography below sea level
  if (.not.flag_lakes) then
    where (geo%hires%mask.eq.1 .and. geo%hires%z_topo.le.0._wp)
      geo%hires%mask = 4  ! lake
    endwhere
  endif

  if (i_fix_cell_grl.eq.1) then
    ! always ocean
    geo%hires%mask(i0_topo(32):i1_topo(32),j0_topo(34):j1_topo(34)) = 2
  else if (i_fix_cell_grl.eq.2) then
    ! never ocean
    where (geo%hires%mask(i0_topo(32):i1_topo(32),j0_topo(34):j1_topo(34))==2)
      geo%hires%mask(i0_topo(32):i1_topo(32),j0_topo(34):j1_topo(34)) = 1
    endwhere
  endif

  ! filtered topography 
  !$ time1 = omp_get_wtime()
  call topo_filter(geo%hires%grid, geo%hires%z_topo, sigma_filter, &    ! in
    geo%hires%z_topo_fil)   ! out
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'topo_filter',time2-time1

  !-------------------------------------------------------------------
  ! compute ocean fraction of coarse grid and make sure all ocean cells are connected
  !-------------------------------------------------------------------
  if (l_connect_ocn) then
    !$ time1 = omp_get_wtime()
    call connect_ocn(geo%hires%mask,real(geo%grid%G%x,wp),real(geo%grid%G%y,wp),lon_ocn_origin,lat_ocn_origin, &
      geo%f_ocn) 
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'connect_ocn',time2-time1
  endif

  !-------------------------------------------------------------------
  ! fill topography depressions 
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  call topo_fill(geo%hires%z_topo, geo%hires%mask, &
    geo%hires%z_topo_fill)
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'topo_fill',time2-time1

  !-------------------------------------------------------------------
  ! update potential lakes
  !-------------------------------------------------------------------
  if (flag_lakes) then

    !$ time1 = omp_get_wtime()
    call lake_pot(geo%hires%z_topo, geo%hires%z_topo_fill, geo%hires%mask, area, & 
      geo%hires%mask_lake_pot, geo%lake)
    geo%n_lakes = size(geo%lake)
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'lake_pot',time2-time1

    ! compute actual lake mask and elevation
    !$ time1 = omp_get_wtime()
    call update_lake_mask(geo%hires%z_topo,area,geo%lake,geo%hires%mask_lake_pot, &
                          geo%hires%mask_lake, geo%hires%z_lake)
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'update_lake_mask',time2-time1

    ! exclude lake points over floating ice and assume ice is grounded
    where (geo%hires%mask_lake.ne.0 .and. geo%hires%mask.eq.3) 
      geo%hires%mask_lake = 0
      geo%hires%mask = 0
    endwhere
    ! add lakes to mask
    where (geo%hires%mask_lake.ne.0)
      geo%hires%mask = 4  ! lake
    endwhere

  endif

  ! update ocean/lake mask
  call ocean_lake_mask_level(geo%hires%mask, geo%hires%z_lake, &
    geo%hires%mask_ocn_lake, geo%hires%z_sea_lake)

  !-------------------------------------------------------------------
  ! update surface elevation
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  geo%hires%z_sur = geo%hires%z_topo
  where (geo%hires%mask.eq.2) geo%hires%z_sur = 0._wp   ! over ocean
  where (geo%hires%mask.eq.4) geo%hires%z_sur = geo%hires%z_lake   ! over lake
  where (geo%hires%z_sur.lt.0._wp) geo%hires%z_sur = 0._wp  ! negative surface elevation not allowed, fixme? could be negative over lakes...
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'update z_sur',time2-time1

  !-------------------------------------------------------------------
  ! update runoff directions of filled topography
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  call runoff_routing(geo%hires%z_topo_fill, geo%hires%mask, real(geo%hires%grid%area,wp), &
    geo%hires%i_runoff, geo%hires%j_runoff, geo%hires%flow_acc)
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'runoff_routing',time2-time1

  !-------------------------------------------------------------------
  ! map runoff to lakes and route lake runoff to ocean
  !-------------------------------------------------------------------
  if (flag_lakes) then
  !$ time1 = omp_get_wtime()
    call runoff_routing_lakes(geo%hires%z_topo, geo%hires%z_topo_fill, geo%hires%mask, geo%hires%mask_lake_pot, &
      geo%hires%i_runoff, geo%hires%j_runoff, &
      geo%hires%map_runoff, geo%lake)
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'runoff_routing_lakes',time2-time1
  endif
                   
  !-------------------------------------------------------------------
  ! runoff destinations to coarse grid indexes
  !-------------------------------------------------------------------
  !$omp parallel do private(i,j)
  do j=1,nj_topo
    do i=1,ni_topo
      if (geo%hires%mask(i,j).ne.2) then  ! not ocean point
        geo%hires%i_runoff_coarse(i,j) = i_lowres(geo%hires%i_runoff(i,j))
        geo%hires%j_runoff_coarse(i,j) = j_lowres(geo%hires%j_runoff(i,j))
      else
        geo%hires%i_runoff_coarse(i,j) = 0
        geo%hires%j_runoff_coarse(i,j) = 0
      endif
    enddo
  enddo
  !$omp end parallel do

  ! lakes routing
  if (flag_lakes) then
    do n=1,geo%n_lakes
      ! determine lake location on coarse grid
      geo%lake(n)%i = i_lowres(geo%lake(n)%hires%i)
      geo%lake(n)%j = j_lowres(geo%lake(n)%hires%j)
      ! determine lake runoff location on coarse grid
      geo%lake(n)%i_runoff = i_lowres(geo%lake(n)%hires%i_runoff)
      geo%lake(n)%j_runoff = j_lowres(geo%lake(n)%hires%j_runoff)
    enddo
  endif

  if (allocated(geo%f_lake_n)) deallocate(geo%f_lake_n)
  if (allocated(geo%f_drain_veg)) deallocate(geo%f_drain_veg)
  if (allocated(geo%f_drain_ice)) deallocate(geo%f_drain_ice)
  allocate(geo%f_lake_n(geo%n_lakes,ni,nj))
  allocate(geo%f_drain_veg(0:geo%n_lakes,ni,nj))
  allocate(geo%f_drain_ice(0:geo%n_lakes,ni,nj))

  !-------------------------------------------------------------------
  ! derive fractions, topography and runoff on coarse resolution
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  call hires_to_lowres(geo%n_lakes, geo%hires%mask, geo%hires%mask_lake, geo%hires%z_topo, geo%hires%z_topo_fil, geo%hires%z_bed, geo%hires%z_sur, & ! in
    geo%hires%map_runoff, geo%hires%i_runoff_coarse, geo%hires%j_runoff_coarse, &   ! in
    geo%f_ocn, geo%f_ocn2, geo%f_lnd, geo%f_ice, geo%f_ice_grd, geo%f_ice_flt, geo%f_lake, geo%f_lake_n, &
    geo%z_sur, geo%z_ocn, geo%z_ocn_min, geo%z_ocn_max, geo%z_ocn_max_q, geo%z_ice, geo%z_lake, geo%z_veg, geo%z_veg_min, geo%z_veg_max, geo%z_bed, & ! out
    geo%z_sur_std, geo%z_sur_smooth_std, geo%z_veg_std, geo%z_sur_lnd_std, &    ! out
    geo%f_drain_veg, geo%f_drain_ice, geo%i_runoff, geo%j_runoff, geo%i_runoff_veg, geo%j_runoff_veg, geo%i_runoff_ice, geo%j_runoff_ice)   ! out
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'hires_to_lowres',time2-time1

  !-------------------------------------------------------------------
  ! make sure runoff ends up in connected ocean cells
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  call fix_runoff(geo%f_ocn, &
    geo%i_runoff, geo%j_runoff)
  ! runoff from vegetated part
  call fix_runoff(geo%f_ocn, &
    geo%i_runoff_veg, geo%j_runoff_veg)
  ! runoff from ice sheets
  call fix_runoff(geo%f_ocn, &
    geo%i_runoff_ice, geo%j_runoff_ice)

  ! make sure lake runoff reaches ocean
  do n=1,geo%n_lakes
    call fix_runoff_lake(geo%f_ocn, &
      geo%lake(n)%i_runoff, geo%lake(n)%j_runoff)
  enddo
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'fix_runoff',time2-time1

  !-------------------------------------------------------------------
  ! derive ocean drainage basins
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  call drainage_basins(geo%hires%i_runoff_coarse, geo%hires%j_runoff_coarse, geo%i_runoff, geo%j_runoff, &
    geo%hires%drain_basins_ocn, geo%drain_basins_ocn, geo%idivide_pac_atl, geo%idivide_atl_indpac)
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'drainage_basins',time2-time1

  !-------------------------------------------------------------------
  ! find coastal cells 
  !-------------------------------------------------------------------
  !$ time1 = omp_get_wtime()
  call coast_cells(geo%f_lnd, geo%i_runoff, geo%j_runoff, &      ! in
    geo%mask_coast, geo%mask_coast2, geo%i_coast_nbr, geo%j_coast_nbr, geo%coast_nbr)   ! out
  !$ time2 = omp_get_wtime()
  !$ if(l_write_timer) print *,'coast_cells',time2-time1

  !-------------------------------------------------------------------
  ! update real ocean volume and area
  !-------------------------------------------------------------------
  geo%ocn_vol_tot = sum(-geo%hires%z_topo*area_dp, mask=geo%hires%mask.eq.2 .or. geo%hires%mask.eq.3)  
  geo%ocn_area_tot = sum(area_dp,geo%hires%mask.eq.2 .or. geo%hires%mask.eq.3)

  !-------------------------------------------------------------------
  ! other (diagnostic only) total areas
  !-------------------------------------------------------------------
  geo%veg_area_tot = sum(area_dp,geo%hires%mask.eq.1)
  geo%ice_area_tot = sum(area_dp,geo%hires%mask.eq.0 .or. geo%hires%mask.eq.3)
  geo%lake_area_tot = sum(area_dp,geo%hires%mask.eq.4)

  !-------------------------------------------------------------------
  ! diagnostic global mean relative sea level 
  !-------------------------------------------------------------------
  if (i_geo==2) then
    geo%sea_level = sum(geo%hires%rsl*area_dp, mask=geo%hires%mask.eq.2 .or. geo%hires%mask.eq.3)/geo%ocn_area_tot
  endif

  !-------------------------------------------------------------------
  ! compute Bering Strait cross-sectional area
  !-------------------------------------------------------------------
  geo%A_bering = 1.e15_wp
  do j=1,nj_topo
    if (geo%hires%grid%lat(1,j)>55._wp .and. geo%hires%grid%lat(1,j)<75._wp) then
      dx = 2._wp*pi*R_earth*cos(pi*geo%hires%grid%lat(1,j)/180._wp)*geo%hires%grid%G%dx/360._wp
      A_bering = 0._wp
      do i=1,ni_topo
        if (geo%hires%grid%lon(i,1)>-180._wp .and. geo%hires%grid%lon(i,1)<-160._wp) then
          A_bering = A_bering + max(0._wp,-geo%hires%z_bed(i,j))*dx
        endif
      enddo
      geo%A_bering = min(geo%A_bering,A_bering)
    endif
  enddo

  if (firstcall) firstcall = .false.


  return

  end subroutine geo_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g e o _ e n d
  ! Purpose  :  end geo model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_end(geo)

    implicit none 

    type(geo_class) :: geo


    if (flag_geo .and. i_geo==2) then
      ! end VILMA
      call vilma_end
    endif

    call geo_dealloc(geo)


   return

  end subroutine geo_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g e o _ d e a l l o c
  ! Purpose  :  deallocate all state variables 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_dealloc(geo)

    implicit none 

    type(geo_class) :: geo


    deallocate(geo%hires%z_bed_ref)
    deallocate(geo%hires%h_ice_ref)
    deallocate(geo%hires%z_topo_ref)
    deallocate(geo%hires%mask_ref)
    deallocate(geo%hires%z_bed_rel)
    deallocate(geo%hires%z_bed)
    deallocate(geo%hires%z_topo)
    deallocate(geo%hires%z_topo_fil)
    deallocate(geo%hires%z_topo_fill)
    deallocate(geo%hires%z_sur)
    deallocate(geo%hires%h_ice)
    deallocate(geo%hires%rsl)
    deallocate(geo%hires%mask)
    deallocate(geo%hires%q_geo)
    deallocate(geo%hires%h_sed)
    deallocate(geo%hires%i_runoff)
    deallocate(geo%hires%j_runoff)
    deallocate(geo%hires%i_runoff_coarse)
    deallocate(geo%hires%j_runoff_coarse)
    deallocate(geo%hires%map_runoff)
    deallocate(geo%hires%flow_acc)
    deallocate(geo%hires%mask_lake_pot)
    deallocate(geo%hires%mask_lake)
    deallocate(geo%hires%z_lake)
    deallocate(geo%hires%z_sea_lake)
    deallocate(geo%hires%mask_ocn_lake)
    deallocate(geo%hires%drain_basins_ocn)

    deallocate(geo%f_ocn0)
    deallocate(geo%f_lnd0)
    deallocate(geo%f_ocn)
    deallocate(geo%f_ocn2)
    deallocate(geo%f_lnd)
    deallocate(geo%f_lake)
    deallocate(geo%f_ice)
    deallocate(geo%f_ice_grd)
    deallocate(geo%f_ice_flt)
    deallocate(geo%z_sur)
    deallocate(geo%z_sur_std)
    deallocate(geo%z_sur_smooth_std)
    deallocate(geo%z_ocn)
    deallocate(geo%z_ocn_min)
    deallocate(geo%z_ocn_max)
    deallocate(geo%z_ocn_max_q)
    deallocate(geo%z_veg)
    deallocate(geo%z_veg_min)
    deallocate(geo%z_veg_max)
    deallocate(geo%z_veg_std)
    deallocate(geo%z_sur_lnd_std)
    deallocate(geo%z_ice)
    deallocate(geo%z_lake)
    deallocate(geo%z_bed)
    deallocate(geo%mask_coast)
    deallocate(geo%mask_coast2)
    deallocate(geo%i_coast_nbr)
    deallocate(geo%j_coast_nbr)
    deallocate(geo%coast_nbr)
    deallocate(geo%i_runoff)
    deallocate(geo%j_runoff)
    deallocate(geo%i_runoff_veg)
    deallocate(geo%j_runoff_veg)
    deallocate(geo%i_runoff_ice)
    deallocate(geo%j_runoff_ice)
    deallocate(geo%drain_basins_ocn)
    deallocate(geo%idivide_pac_atl)
    deallocate(geo%idivide_atl_indpac)
    deallocate(geo%coral_f_area)
    deallocate(geo%coral_f_topo)
    deallocate(geo%q_geo)


   return

  end subroutine geo_dealloc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g e o _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_write_restart(dir,fnm,geo)

    implicit none

    character (len=*) :: dir
    character (len=*) :: fnm
    type(geo_class) :: geo


    call nc_create(fnm)
    call nc_write_dim(fnm,"g",x=[1])
    call nc_write_dim(fnm,"lon",x=geo%hires%grid%lon(:,1),axis="x")
    call nc_write_dim(fnm,"lat",x=geo%hires%grid%lat(1,:),axis="y")

    call nc_write(fnm,"sea_level", geo%sea_level, dim1="g", long_name="sea_level",units="m")
    call nc_write(fnm,"V_ice_af", geo%V_ice_af, dim1="g", long_name="volume of ice above floatation",units="m3")

    call nc_write(fnm,"z_bed ", geo%hires%z_bed , dims=["lon","lat"],long_name="z_bed ",units="m")
    call nc_write(fnm,"z_topo", geo%hires%z_topo, dims=["lon","lat"],long_name="z_topo",units="m") 
    call nc_write(fnm,"h_ice ", geo%hires%h_ice , dims=["lon","lat"],long_name="h_ice ",units="m")
    call nc_write(fnm,"rsl   ", geo%hires%rsl   , dims=["lon","lat"],long_name="rsl   ",units="m")
    call nc_write(fnm,"mask  ", geo%hires%mask  , dims=["lon","lat"],long_name="mask  ",units="1")

    if (flag_geo .and. i_geo.eq.2) then
      call vilma_write_restart(dir)
    endif

   return

  end subroutine geo_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g e o _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine geo_read_restart(fnm,geo)

    implicit none

    character (len=*) :: fnm
    type(geo_class) :: geo

    call nc_read(fnm,"sea_level", geo%sea_level      )
    call nc_read(fnm,"V_ice_af",  geo%V_ice_af       )

    call nc_read(fnm,"z_bed    ", geo%hires%z_bed    )
    call nc_read(fnm,"z_topo   ", geo%hires%z_topo   ) 
    call nc_read(fnm,"h_ice    ", geo%hires%h_ice    )
    call nc_read(fnm,"rsl      ", geo%hires%rsl      )
    call nc_read(fnm,"mask     ", geo%hires%mask     )

    print *,'read restart file ',fnm


   return

  end subroutine geo_read_restart

end module geo_mod

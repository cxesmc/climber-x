!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b m b _ o u t
!
!  Purpose : diagnostics and output of BMB model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module bmb_out

  use precision, only : wp
  use dim_name, only: dim_x, dim_y, dim_time, dim_month
  use constants, only : rho_i
  use timer, only : time_soy_bmb, time_eoy_bmb, time_eom_bmb, time_out_bmb
  use timer, only : year_clim, year_now, mon, time_out_ts, ny_out_ts, y_out_ts_clim, n_accel, nmon_year, nstep_mon_bmb, nstep_year_bmb
  use control, only : out_dir
  use bmb_def, only : bmb_class, ts_out, s_out 
  use bmb_params, only : l_monthly_output
  use ncio
  use coord, only : grid_class

  implicit none

  private
  public :: bmb_diag, bmb_diag_init


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for bmb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_diag_init(bmb)

    implicit none

    type(bmb_class) :: bmb

    integer :: k
    character (len=256) :: fnm

    bmb%nout = 0

    ! allocate
    allocate(bmb%ann_ts(ny_out_ts))

    ! initialize netcdf output
    fnm = trim(out_dir)//"/bmb_"//trim(bmb%grid%name)//"_ts.nc"
    call ts_nc(fnm)

    fnm = trim(out_dir)//"/bmb_"//trim(bmb%grid%name)//".nc"
    call bmb_nc(fnm,bmb%grid)

    allocate(bmb%ann_s%mask_ocn_lake   (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%mask_ice_shelf   (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%zb         (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%t_bmb    (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%s_bmb    (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%t_disc    (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%s_disc    (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%t_freeze   (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%bmb        (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%bmb_ann    (bmb%grid%G%nx,bmb%grid%G%ny))
    allocate(bmb%ann_s%bmb_ann_mask(bmb%grid%G%nx,bmb%grid%G%ny))

    do k=1,nmon_year
      allocate(bmb%mon_s(k)%mask_ice_shelf   (bmb%grid%G%nx,bmb%grid%G%ny))
      allocate(bmb%mon_s(k)%zb         (bmb%grid%G%nx,bmb%grid%G%ny))
      allocate(bmb%mon_s(k)%t_bmb      (bmb%grid%G%nx,bmb%grid%G%ny))
      allocate(bmb%mon_s(k)%s_bmb      (bmb%grid%G%nx,bmb%grid%G%ny))
      allocate(bmb%mon_s(k)%t_freeze   (bmb%grid%G%nx,bmb%grid%G%ny))
      allocate(bmb%mon_s(k)%bmb        (bmb%grid%G%nx,bmb%grid%G%ny))
    enddo


   return

  end subroutine bmb_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b m b _ d i a g
  !   Purpose    :  sea ice diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_diag(bmb)

    implicit none

    type(bmb_class) :: bmb

    integer :: m, y
    integer :: ppos
    character (len=256) :: fnm
    character (len=256) :: dom
    real(wp) :: mon_avg, ann_avg


    ! current index
    y = y_out_ts_clim

    mon_avg = 1._wp/nstep_mon_bmb
    ann_avg = 1._wp/nstep_year_bmb

    if( time_eoy_bmb ) then

     ! integrated basal mass balance of ice shelf
     bmb%ann_ts(y)%bmb = sum(bmb%bmb_ann*bmb%grid%area,mask=bmb%mask_ice_shelf==1) * 1.e-6_wp ! kg/m2/a * km2 * m2/km2 * Gt/kg = Gt/a

     ! write to netcdf file 
     if (time_out_ts) then
       fnm = trim(out_dir)//"/bmb_"//trim(bmb%grid%name)//"_ts.nc"
       call ts_nc_write(fnm,bmb%ann_ts(1:y),year_clim-y+1,y)
     endif
     ppos = scan(trim(bmb%grid%name),"-")-1
     dom = trim(bmb%grid%name(1:ppos)) 
     ! print header
     if (mod(year_clim,10).eq.1) then
       print '(a7,a9,a7)','bmb_'//dom,'year','bmb'
     endif

     ! print values
     print '(a7,i9,1F7.0)', &
       'bmb_'//dom,year_now,bmb%ann_ts(y)%bmb

    endif


    ! spatially explicit output
    if ( time_out_bmb ) then

      if( time_soy_bmb ) then
        do m=1,nmon_year
          bmb%mon_s(m)%mask_ice_shelf   = 0._wp 
          bmb%mon_s(m)%zb         = 0._wp 
          bmb%mon_s(m)%t_bmb      = 0._wp 
          bmb%mon_s(m)%s_bmb      = 0._wp 
          bmb%mon_s(m)%t_freeze   = 0._wp 
          bmb%mon_s(m)%bmb        = 0._wp 
        enddo
      endif

      bmb%mon_s(mon)%mask_ice_shelf   = bmb%mon_s(mon)%mask_ice_shelf   + bmb%mask_ice_shelf   * mon_avg
      bmb%mon_s(mon)%zb         = bmb%mon_s(mon)%zb         + bmb%zb         * mon_avg
      bmb%mon_s(mon)%t_bmb      = bmb%mon_s(mon)%t_bmb      + bmb%t_bmb      * mon_avg
      bmb%mon_s(mon)%s_bmb      = bmb%mon_s(mon)%s_bmb      + bmb%s_bmb      * mon_avg
      bmb%mon_s(mon)%t_freeze   = bmb%mon_s(mon)%t_freeze   + bmb%t_freeze   * mon_avg
      bmb%mon_s(mon)%bmb        = bmb%mon_s(mon)%bmb        + bmb%bmb        * mon_avg

    endif

    if (time_out_bmb .and. time_eoy_bmb) then

      bmb%ann_s%bmb_ann = bmb%bmb_ann / rho_i   ! m(ice)/a

      bmb%ann_s%t_disc = bmb%t_disc
      bmb%ann_s%s_disc = bmb%s_disc

      bmb%nout = bmb%nout+1
      call bmb_diag_out(bmb,bmb%grid)
    endif


   return

  end subroutine bmb_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ d i a g _ o u t
  ! Purpose  :  write sea ice netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_diag_out(bmb,grid)

    implicit none

    type(bmb_class) :: bmb
    type(grid_class) :: grid
    integer :: k, ncid
    character (len=256) :: fnm


    ! Get annual values
    call bmb_ave( bmb%mon_s,bmb%ann_s )
    
    bmb%ann_s%mask_ocn_lake = bmb%mask_ocn_lake
    bmb%ann_s%bmb_ann_mask = bmb%ann_s%bmb_ann * bmb%ann_s%mask_ice_shelf  ! m(ice)/a

    ! write to file
    fnm = trim(out_dir)//"/bmb_"//trim(bmb%grid%name)//".nc"
    call nc_open(fnm,ncid)
    call nc_write(fnm,dim_time,dble(year_now), dim1=dim_time, start=[bmb%nout], count=[1],ncid=ncid)    
    do k = 1, nmon_year
       call bmb_nc_write(fnm,ncid,bmb%mon_s(k),grid,k,bmb%nout)
    end do
    call bmb_nc_write(fnm,ncid,bmb%ann_s,grid,nmon_year+1,bmb%nout)
    call nc_close(ncid)


   return

  end subroutine bmb_diag_out
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  initialize netcdf file for time series output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", unlimited=.TRUE.)

    return

  end subroutine ts_nc
 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  write time series to netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,nout,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: nout, ncid, y, i

    call nc_open(fnm,ncid)
    call nc_write(fnm,"time",  real([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))],wp), dim1=dim_time,start=[nout],count=[y],ncid=ncid)    
    call nc_write(fnm,"bmb", vars%bmb, dims=[dim_time],start=[nout],count=[y],&
    long_name="integrated floating ice basal mass balance of ice sheets",units="Gt/a",ncid=ncid)
    call nc_close(ncid)


   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ n c
  ! Purpose  :  Initialize bmb netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_nc(fnm,grid)

    implicit none

    character (len=*) :: fnm
    type(grid_class) :: grid
    integer :: ncid
    real(wp) :: empty_time(0)


    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm, dim_month, x=1._wp, dx=1._wp, nx=13, axis="e", &
    units="months", ncid=ncid)
    call nc_write_dim(fnm, dim_y, x=grid%G%y0, dx=grid%G%dy, nx=grid%G%ny,&
    axis="y", units="km", ncid=ncid)
    call nc_write_dim(fnm, dim_x, x=grid%G%x0, dx=grid%G%dx, nx=grid%G%nx,&
    axis="x", units="km", ncid=ncid)

    call nc_write(fnm,"phi",sngl(grid%lat), dims=[dim_x,dim_y],start=[1,1],count=[grid%G%nx,grid%G%ny],ncid=ncid)
    call nc_write(fnm,"lambda",sngl(grid%lon), dims=[dim_x,dim_y],start=[1,1],count=[grid%G%nx,grid%G%ny],ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine bmb_nc


  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ n c _ w r i t e
  ! Purpose  :  Output of bmb netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_nc_write(fnm,ncid,vars,grid,ndat,nout)

    implicit none

    type(s_out) :: vars
    type(grid_class) :: grid

    character (len=*) :: fnm
    integer :: ndat, nout, ncid


    if (l_monthly_output) then
      call nc_write(fnm,"t_bmb", sngl(vars%t_bmb   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water temperature used for basal mass balance",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"s_bmb", sngl(vars%s_bmb   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water salinity used for basal mass balance",grid_mapping="polar_stereographic",units="psu",ncid=ncid) 
      call nc_write(fnm,"t_freeze", sngl(vars%t_freeze   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="freezing temperature at the ice shelf base",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"bmb", sngl(vars%bmb   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="basal mass balance of floating ice",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid) 
      call nc_write(fnm,"t_bmb_mask", sngl(vars%t_bmb*vars%mask_ice_shelf), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water temperature used for basal mass balance, masked with ice shelf mask",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"s_bmb_mask", sngl(vars%s_bmb*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water salinity used for basal mass balance, masked with ice shelf mask",grid_mapping="polar_stereographic",units="psu",ncid=ncid) 
      call nc_write(fnm,"t_freeze_mask", sngl(vars%t_freeze*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="freezing temperature at the ice shelf base, masked with ice shelf mask",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"bmb_mask", sngl(vars%bmb*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_month,dim_time],start=[1,1,ndat,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="basal mass balance of floating ice, masked with ice shelf mask",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid) 
    else
    if (ndat.eq.13) then
      call nc_write(fnm,"t_bmb", sngl(vars%t_bmb   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water temperature used for basal mass balance",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"s_bmb", sngl(vars%s_bmb   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water salinity used for basal mass balance",grid_mapping="polar_stereographic",units="psu",ncid=ncid) 
      call nc_write(fnm,"t_freeze", sngl(vars%t_freeze   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="freezing temperature at the ice shelf base",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"bmb", sngl(vars%bmb   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="basal mass balance of floating ice",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid) 
      call nc_write(fnm,"t_bmb_mask", sngl(vars%t_bmb*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water temperature used for basal mass balance, masked with ice shelf mask",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"s_bmb_mask", sngl(vars%s_bmb*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water salinity used for basal mass balance, masked with ice shelf mask",grid_mapping="polar_stereographic",units="psu",ncid=ncid) 
      call nc_write(fnm,"t_freeze_mask", sngl(vars%t_freeze*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="freezing temperature at the ice shelf base, masked with ice shelf mask",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"bmb_mask", sngl(vars%bmb*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="basal mass balance of floating ice, masked with ice shelf mask",grid_mapping="polar_stereographic",units="kg/m2/s",ncid=ncid) 
    endif
    endif
    if (ndat.eq.13) then
      call nc_write(fnm,"t_disc", sngl(vars%t_disc   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water temperature used for small-scale basal melt in ice sheet_model",grid_mapping="polar_stereographic",units="degC",ncid=ncid) 
      call nc_write(fnm,"s_disc", sngl(vars%s_disc   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="water salinity used for small-scale basal melt in ice sheet model",grid_mapping="polar_stereographic",units="psu",ncid=ncid) 
      call nc_write(fnm,"mask_ocn_lake", vars%mask_ocn_lake, dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="ocean/lake mask",grid_mapping="polar_stereographic",units="",ncid=ncid) 
      call nc_write(fnm,"mask_ice_shelf", sngl(vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="ice mask",grid_mapping="polar_stereographic",units="",ncid=ncid) 
      call nc_write(fnm,"zb", sngl(vars%zb   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="elevation of ice base relative to sea level",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
      call nc_write(fnm,"zb_mask", sngl(vars%zb*vars%mask_ice_shelf   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1,1], &
        long_name="elevation of ice base relative to sea level, masked with ice shelf mask",grid_mapping="polar_stereographic",units="m",ncid=ncid) 
      call nc_write(fnm,"bmb_ann", sngl(vars%bmb_ann   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1], &
        long_name="annual basal mass balance of floating ice",grid_mapping="polar_stereographic",units="m/a",ncid=ncid) 
      call nc_write(fnm,"bmb_ann_mask", sngl(vars%bmb_ann_mask   ), dims=[dim_x,dim_y,dim_time],start=[1,1,nout],count=[grid%G%nx,grid%G%ny,1], &
        long_name="annual basal mass balance of floating ice, masked with ice shelf mask",grid_mapping="polar_stereographic",units="m/a",ncid=ncid) 
    endif

   return

  end subroutine bmb_nc_write



  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ a v e
  ! Purpose  :  Average (or sum) the sea ice fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_ave(d,ave)

    implicit none

    type(s_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div

    n = size(d)
    div = dble(n)

    ! Set all values to zero
    ave%mask_ice_shelf = 0._wp 
    ave%zb       = 0._wp 
    ave%t_bmb    = 0._wp 
    ave%s_bmb    = 0._wp 
    ave%t_freeze = 0._wp 
    ave%bmb      = 0._wp 

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
      ave%mask_ice_shelf  = ave%mask_ice_shelf  + d(k)%mask_ice_shelf  / div       
      ave%zb        = ave%zb        + d(k)%zb        / div       
      ave%t_bmb     = ave%t_bmb     + d(k)%t_bmb     / div       
      ave%s_bmb     = ave%s_bmb     + d(k)%s_bmb     / div       
      ave%t_freeze  = ave%t_freeze  + d(k)%t_freeze  / div       
      ave%bmb       = ave%bmb       + d(k)%bmb       / div       
    end do

   return

  end subroutine bmb_ave


end module bmb_out

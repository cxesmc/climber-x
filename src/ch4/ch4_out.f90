!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c h 4 _ o u t
!
!  Purpose : atmospheric CH4 model diagnostic and output
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module ch4_out

  use precision, only : wp
  use dim_name, only: dim_time, dim_lon, dim_lat
  use timer, only : n_accel, year, year_clim, year_now, ny_out_ts, y_out_ts_clim, &
  time_out_ts_clim
  use control, only : out_dir
  use ch4_def, only : ch4_class
  use ncio

  implicit none

  type ts_out
      real(wp) :: ch4       !! atmospheric ch4 concentration [ppm]
      real(wp) :: dch4ocn_dt
      real(wp) :: dch4lnd_dt
      real(wp) :: dch4emis_dt
      real(wp) :: dch4ox_dt
      real(wp) :: demis_other_dt
      real(wp) :: tau
  end type

  type(ts_out), allocatable :: ann_ts(:)

  private
  public :: ch4_diag_init, ch4_diag

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c h 4 _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_diag_init

    implicit none


    ! allocate
    allocate(ann_ts(ny_out_ts))

    ! initialize time series output
    call ts_nc(trim(out_dir)//"/ch4_ts.nc")

   return

  end subroutine ch4_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c h 4 _ d i a g
  !   Purpose    :  ch4 diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ch4_diag(ch4)

    implicit none

    type(ch4_class), intent(in) :: ch4

    integer :: y


    ! current index
    y = y_out_ts_clim

    ! sum up and average over the year
    ann_ts(y)%ch4        = ch4%ch4           
    ann_ts(y)%dch4ocn_dt   = ch4%dch4ocn_dt * 1.e-9_wp    ! TgCH4/yr, positive into the atmosphere
    ann_ts(y)%dch4lnd_dt   = ch4%dch4lnd_dt * 1.e-9_wp    ! TgCH4/yr, positive into the atmosphere
    ann_ts(y)%dch4emis_dt  = ch4%dch4emis_dt * 1.e-9_wp   ! TgCH4/yr, positive into the atmosphere
    ann_ts(y)%dch4ox_dt    = ch4%dch4ox_dt * 1.e-9_wp     ! TgCH4/yr, positive into the atmosphere
    ann_ts(y)%demis_other_dt   = ch4%demis_other_dt * 1.e-9_wp     ! TgCH4/yr, positive into the atmosphere
    ann_ts(y)%tau          = ch4%tau    ! years 

    ! write to standard output
    if (mod(year,10).eq.1) then
      print '(a7,a9,6a7)','ch4','year','CH4','CH4ocn','CH4lnd','CH4emis','CH4ox','CH4otherEmis','CH4tau'
    endif

    print '(a7,i9,F7.1,5F7.1)', &
      'ch4',year_now,ann_ts(y)%ch4,ann_ts(y)%dch4ocn_dt,ann_ts(y)%dch4lnd_dt,ann_ts(y)%dch4emis_dt,ann_ts(y)%dch4ox_dt,ann_ts(y)%demis_other_dt,ann_ts(y)%tau

    ! write to netcdf file 
    if (time_out_ts_clim) then
      call ts_nc_write(trim(out_dir)//"/ch4_ts.nc",ann_ts(1:y),year_clim-y+1,y)
    endif


   return

  end subroutine ch4_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  initialize netcdf file for time series output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(8) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm, dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.)
    call nc_write_dim(fnm, dim_lat, x=1, axis="y", units="1")
    call nc_write_dim(fnm, dim_lon, x=1, axis="x", units="1")

    return

  end subroutine ts_nc
 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  write time series to netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,ndat,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: ndat, y, ncid, i

    call nc_open(fnm,ncid)
    call nc_write(fnm,"time",  dble([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))]), dim1=dim_time,start=[ndat],count=[y],ncid=ncid)    
    call nc_write(fnm,"ch4       ",  vars%ch4       , dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric ch4 concentration",units="ppb",ncid=ncid) 
    call nc_write(fnm,"CH4ocn_dt  ",  vars%dch4ocn_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="ocean methane flux to atmosphere",units="TgCH4/yr",ncid=ncid)  
    call nc_write(fnm,"CH4lnd_dt  ",  vars%dch4lnd_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="natural land methane flux to atmosphere",units="TgCH4/yr",ncid=ncid)  
    call nc_write(fnm,"CH4emis_dt ",  vars%dch4emis_dt , dim1=dim_time,start=[ndat],count=[y],long_name="anthropogenic methane emissions to atmosphere",units="TgCH4/yr",ncid=ncid)  
    call nc_write(fnm,"CH4ox_dt   ",  vars%dch4ox_dt   , dim1=dim_time,start=[ndat],count=[y],long_name="methane oxidation in the atmosphere",units="TgCH4/yr",ncid=ncid)  
    call nc_write(fnm,"tau   ",  vars%tau   , dim1=dim_time,start=[ndat],count=[y],long_name="methane lifetime in the atmosphere",units="years",ncid=ncid)  
    call nc_close(ncid)

   return

  end subroutine ts_nc_write

end module ch4_out

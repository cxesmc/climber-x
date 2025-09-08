!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : n 2 o _ o u t
!
!  Purpose : atmospheric N2O model diagnostic and output
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
module n2o_out

  use precision, only : wp
  use dim_name, only: dim_time, dim_lon, dim_lat
  use timer, only : n_accel, year, year_clim, year_now, ny_out_ts, y_out_ts_clim, &
  time_out_ts_clim
  use control, only : out_dir
  use n2o_def, only : n2o_class
  use ncio

  implicit none

  type ts_out
      real(wp) :: n2o       !! atmospheric n2o concentration [ppb]
      real(wp) :: dn2oocn_dt
      real(wp) :: dn2olnd_dt
      real(wp) :: dn2oemis_dt
      real(wp) :: dn2oox_dt
      real(wp) :: tau
  end type

  type(ts_out), allocatable :: ann_ts(:)

  private
  public :: n2o_diag_init, n2o_diag

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  n 2 o _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_diag_init

    implicit none


    ! allocate
    allocate(ann_ts(ny_out_ts))

    ! initialize time series output
    call ts_nc(trim(out_dir)//"/n2o_ts.nc")

   return

  end subroutine n2o_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  n 2 o _ d i a g
  !   Purpose    :  n2o diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_diag(n2o)

    implicit none

    type(n2o_class), intent(in) :: n2o

    integer :: y


    ! current index
    y = y_out_ts_clim

    ! sum up and average over the year
    ann_ts(y)%n2o        = n2o%n2o           
    ann_ts(y)%dn2oocn_dt   = n2o%dn2oocn_dt * 1.e-9_wp    ! TgN2O-N/yr, positive into the atmosphere
    ann_ts(y)%dn2olnd_dt   = n2o%dn2olnd_dt * 1.e-9_wp    ! TgN2O-N/yr, positive into the atmosphere
    ann_ts(y)%dn2oemis_dt  = n2o%dn2oemis_dt * 1.e-9_wp   ! TgN2O-N/yr, positive into the atmosphere
    ann_ts(y)%dn2oox_dt    = n2o%dn2oox_dt * 1.e-9_wp     ! TgN2O-N/yr
    ann_ts(y)%tau          = n2o%tau    ! years 

    ! write to standard output
    if (mod(year,10).eq.1) then
      print '(a7,a9,6a7)','n2o','year','N2O','N2Oocn','N2Olnd','N2Oe','N2Oox','N2Otau'
    endif

    print '(a7,i9,6F7.1)', &
      'n2o',year_now,ann_ts(y)%n2o,ann_ts(y)%dn2oocn_dt,ann_ts(y)%dn2olnd_dt,ann_ts(y)%dn2oemis_dt,ann_ts(y)%dn2oox_dt,ann_ts(y)%tau

    ! write to netcdf file 
    if (time_out_ts_clim) then
      call ts_nc_write(trim(out_dir)//"/n2o_ts.nc",ann_ts(1:y),year_clim-y+1,y)
    endif


   return

  end subroutine n2o_diag


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
    call nc_write(fnm,"n2o       ",  vars%n2o       , dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric n2o concentration",units="ppb",ncid=ncid) 
    call nc_write(fnm,"N2Oocn_dt  ",  vars%dn2oocn_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="ocean n2o flux to atmosphere",units="TgN2O-N/yr",ncid=ncid)  
    call nc_write(fnm,"N2Olnd_dt  ",  vars%dn2olnd_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="natural land n2o flux to atmosphere",units="TgN2O-N/yr",ncid=ncid)  
    call nc_write(fnm,"N2Oemis_dt ",  vars%dn2oemis_dt , dim1=dim_time,start=[ndat],count=[y],long_name="anthropogenic n2o emissions to atmosphere",units="TgN2O-N/yr",ncid=ncid)  
    call nc_write(fnm,"N2Oox_dt   ",  vars%dn2oox_dt   , dim1=dim_time,start=[ndat],count=[y],long_name="n2o removal rate in the atmosphere",units="TgN2O-N/yr",ncid=ncid)  
    call nc_write(fnm,"tau   ",  vars%tau   , dim1=dim_time,start=[ndat],count=[y],long_name="n2o lifetime in the atmosphere",units="years",ncid=ncid)  
    call nc_close(ncid)

   return

  end subroutine ts_nc_write

end module n2o_out

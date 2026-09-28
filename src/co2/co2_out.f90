!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o 2 _ o u t
!
!  Purpose : atmospheric CO2 model diagnostic and output
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
module co2_out

  use precision, only : wp
  use dim_name, only: dim_time, dim_lon, dim_lat
  use timer, only : n_accel, year, year_clim, year_now, ny_out_ts, y_out_ts_clim, time_out_ts_clim
  use control, only : out_dir
  use constants, only : c13_c12_std, c14_c_std
  use co2_def, only : co2_class
  use ncio

  implicit none

  type ts_out
      real(wp) :: co2       !! atmospheric CO2 concentration [ppm]
      real(wp) :: Catm      !! atmospheric carbon content [GtC]
      real(wp) :: C13atm      !! atmospheric carbon 13 content [GtC]
      real(wp) :: C14atm      !! atmospheric carbon 14 content [kgC]
      real(wp) :: d13C
      real(wp) :: D14C
      real(wp) :: dCocn_dt
      real(wp) :: dCocn_cum
      real(wp) :: dC13ocn_dt
      real(wp) :: dC14ocn_dt
      real(wp) :: d13C_ocn
      real(wp) :: dClnd_dt
      real(wp) :: dClnd_cum
      real(wp) :: dC13lnd_dt
      real(wp) :: dC14lnd_dt
      real(wp) :: d13C_lnd
      real(wp) :: dCemis_dt
      real(wp) :: dCemis_cum
      real(wp) :: dCemis_extra_dt
      real(wp) :: dCemis_extra_cum
      real(wp) :: dCweath_dt
      real(wp) :: dCweath_cum
      real(wp) :: dC13weath_dt
      real(wp) :: dC14weath_dt
      real(wp) :: dCvolc_dt
      real(wp) :: dCvolc_cum
      real(wp) :: dC13volc_dt
      real(wp) :: dC14dec_dt
      real(wp) :: dC14prod_dt
      real(wp) :: dCH4_dt
  end type

  type(ts_out), allocatable :: ann_ts(:)

  real(wp) :: dCocn_cum, dClnd_cum, dCemis_cum, dCemis_extra_cum, dCweath_cum, dCvolc_cum

  private
  public :: co2_diag_init, co2_diag

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o 2 _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_diag_init

    implicit none

    ! allocate
    allocate(ann_ts(ny_out_ts))

    ! initialize time series output
    call ts_nc(trim(out_dir)//"/co2_ts.nc")

    ! initialize cumulated values to zero
    dCocn_cum = 0._wp
    dClnd_cum = 0._wp
    dCweath_cum = 0._wp
    dCvolc_cum = 0._wp
    dCemis_cum = 0._wp
    dCemis_extra_cum = 0._wp

   return

  end subroutine co2_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c o 2 _ d i a g
  !   Purpose    :  CO2 diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_diag(co2)

    implicit none

    type(co2_class), intent(in) :: co2

    integer :: y


    ! current index
    y = y_out_ts_clim

    ! sum up and average over the year
    ann_ts(y)%co2        = co2%co2           
    ann_ts(y)%Catm       = co2%Catm    * 1.e-12_wp     ! GtC
    ann_ts(y)%C13atm     = co2%C13atm  * 1.e-12_wp     ! GtC
    ann_ts(y)%C14atm     = co2%C14atm       ! kgC
    ann_ts(y)%d13C       = (co2%c13_c12/c13_c12_std - 1._wp) * 1000._wp
    ann_ts(y)%D14C       = (co2%c14_c*(0.975_wp/(1._wp+max(-999._wp,ann_ts(y)%d13C)/1000._wp))**2 / c14_c_std - 1._wp) * 1000._wp
    ann_ts(y)%dCocn_dt   = -co2%dCocn_dt * 1.e-12_wp    ! GtC/yr, positive into the atmosphere
    dCocn_cum  = dCocn_cum - co2%dCocn_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dCocn_cum  = dCocn_cum 
    ann_ts(y)%dC13ocn_dt = -co2%dC13ocn_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dC14ocn_dt = -co2%dC14ocn_dt ! kgC
    ann_ts(y)%d13C_ocn   = ((co2%dC13ocn_dt/co2%dCocn_dt)/c13_c12_std - 1._wp) * 1000._wp
    ann_ts(y)%dClnd_dt   = -co2%dClnd_dt * 1.e-12_wp    ! GtC/yr, positive into the atmosphere
    dClnd_cum  = dClnd_cum - co2%dClnd_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dClnd_cum  = dClnd_cum 
    ann_ts(y)%dC13lnd_dt = -co2%dC13lnd_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dC14lnd_dt = -co2%dC14lnd_dt  ! kgC 
    ann_ts(y)%d13C_lnd   = ((co2%dC13lnd_dt/co2%dClnd_dt)/c13_c12_std - 1._wp) * 1000._wp
    ann_ts(y)%dCemis_dt  = co2%dCemis_dt * 1.e-12_wp    ! GtC/yr, positive into the atmosphere
    dCemis_cum  = dCemis_cum + co2%dCemis_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dCemis_cum  = dCemis_cum 
    ann_ts(y)%dCemis_extra_dt  = co2%dCemis_extra_dt * 1.e-12_wp    ! GtC/yr, positive into the atmosphere
    dCemis_extra_cum  = dCemis_extra_cum + co2%dCemis_extra_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dCemis_extra_cum  = dCemis_extra_cum 
    ann_ts(y)%dCweath_dt   = -co2%dCweath_dt * 1.e-12_wp    ! GtC/yr, positive into the atmosphere
    dCweath_cum  = dCweath_cum - co2%dCweath_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dCweath_cum  = dCweath_cum 
    ann_ts(y)%dC13weath_dt = -co2%dC13weath_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dC14weath_dt = -co2%dC14weath_dt  ! kgC
    ann_ts(y)%dCvolc_dt   = co2%dCvolc_dt * 1.e-12_wp    ! GtC/yr, positive into the atmosphere
    dCvolc_cum  = dCvolc_cum + co2%dCvolc_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dCvolc_cum  = dCvolc_cum 
    ann_ts(y)%dC13volc_dt = co2%dC13volc_dt * 1.e-12_wp ! GtC
    ann_ts(y)%dC14dec_dt = co2%dC14dec_dt ! kgC 
    ann_ts(y)%dC14prod_dt = co2%dC14prod_dt ! kgC 
    ann_ts(y)%dCH4_dt   = co2%dCH4_dt * 1.e-12_wp    ! GtC/yr, positive into the atmosphere


    ! write to standard output
    if (mod(year,10).eq.1) then
      print '(a7,a9,10a7)','co2','year','co2','dCocn','dClnd','Cemis','CH4ox','Cfb','dCweath','dCvolc','d13C','D14C'
    endif

    print '(a7,i9,F7.1,7F7.3,2F7.1)', &
      'co2',year_now,ann_ts(y)%co2,ann_ts(y)%dCocn_dt,ann_ts(y)%dClnd_dt,ann_ts(y)%dCemis_dt,ann_ts(y)%dCH4_dt,ann_ts(y)%dCemis_extra_dt, &
      ann_ts(y)%dCweath_dt,ann_ts(y)%dCvolc_dt,ann_ts(y)%d13C,ann_ts(y)%D14C

    ! write to netcdf file 
    if (time_out_ts_clim) then
      call ts_nc_write(trim(out_dir)//"/co2_ts.nc",ann_ts(1:y),year_clim-y+1,y)
    endif


   return

  end subroutine co2_diag


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
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
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

    call nc_open(fnm,ncid,writable=.TRUE.)
    call nc_write(fnm,"time", dble([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))]), &
    dim1=dim_time,start=[ndat],count=[y],ncid=ncid)    
    call nc_write(fnm,"co2", vars%co2, dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric co2 concentration",units="ppm",ncid=ncid) 
    call nc_write(fnm,"Catm      ",  vars%Catm      , dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric carbon",units="GtC",ncid=ncid)  
    call nc_write(fnm,"C13atm    ",  vars%C13atm    , dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric carbon 13",units="GtC",ncid=ncid)  
    call nc_write(fnm,"C14atm    ",  vars%C14atm    , dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric carbon 14",units="kgC",ncid=ncid)  
    call nc_write(fnm,"d13C      ",  vars%d13C      , dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric delta 13 C",units="permil",ncid=ncid)  
    call nc_write(fnm,"D14C      ",  vars%D14C      , dim1=dim_time,start=[ndat],count=[y],long_name="atmospheric Delta 14 C",units="permil",ncid=ncid)  
    call nc_write(fnm,"dCocn_dt  ",  vars%dCocn_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="net ocean carbon flux to atmosphere",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dCocn_cum ",  vars%dCocn_cum , dim1=dim_time,start=[ndat],count=[y],long_name="cumulative net ocean carbon flux to atmosphere",units="GtC",ncid=ncid)  
    call nc_write(fnm,"dC13ocn_dt",  vars%dC13ocn_dt, dim1=dim_time,start=[ndat],count=[y],long_name="net ocean carbon 13 flux to atmosphere",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dC14ocn_dt",  vars%dC14ocn_dt, dim1=dim_time,start=[ndat],count=[y],long_name="net ocean carbon 14 flux to atmosphere",units="kgC/yr",ncid=ncid)  
    call nc_write(fnm,"d13C_ocn  ",  vars%d13C_ocn  , dim1=dim_time,start=[ndat],count=[y],long_name="delta 13 C of ocean-atmosphere flux",units="permil",ncid=ncid)  
    call nc_write(fnm,"dClnd_dt  ",  vars%dClnd_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="net land carbon flux to atmosphere",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dClnd_cum ",  vars%dClnd_cum , dim1=dim_time,start=[ndat],count=[y],long_name="cumulative net land carbon flux to atmosphere ",units="GtC",ncid=ncid)  
    call nc_write(fnm,"dC13lnd_dt",  vars%dC13lnd_dt, dim1=dim_time,start=[ndat],count=[y],long_name="net land carbon 13 flux to atmosphere",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dC14lnd_dt",  vars%dC14lnd_dt, dim1=dim_time,start=[ndat],count=[y],long_name="net land carbon 14 flux to atmosphere",units="kgC/yr",ncid=ncid)  
    call nc_write(fnm,"d13C_lnd  ",  vars%d13C_lnd  , dim1=dim_time,start=[ndat],count=[y],long_name="delta 13 C of land-atmosphere flux",units="permil",ncid=ncid)  
    call nc_write(fnm,"dCemis_dt ",  vars%dCemis_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="carbon emissions to atmosphere",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dCemis_cum",  vars%dCemis_cum , dim1=dim_time,start=[ndat],count=[y],long_name="cumulative carbon emissions to atmosphere",units="GtC",ncid=ncid)  
    call nc_write(fnm,"dCemis_extra_dt ",  vars%dCemis_extra_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="carbon emissions to atmosphere from additional feedbacks",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dCemis_extra_cum",  vars%dCemis_extra_cum , dim1=dim_time,start=[ndat],count=[y],long_name="cumulative carbon emissions to atmosphere from additional feedbacks",units="GtC",ncid=ncid)  
    call nc_write(fnm,"dCweath_dt  ",  vars%dCweath_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="CO2 consumption by weathering",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dCweath_cum ",  vars%dCweath_cum , dim1=dim_time,start=[ndat],count=[y],long_name="cumulative CO2 consumption by weathering",units="GtC",ncid=ncid)  
    call nc_write(fnm,"dC13weath_dt",  vars%dC13weath_dt, dim1=dim_time,start=[ndat],count=[y],long_name="CO2 13 consumption by weathering",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dC14weath_dt",  vars%dC14weath_dt, dim1=dim_time,start=[ndat],count=[y],long_name="CO2 14 consumption by weathering",units="kgC/yr",ncid=ncid)  
    call nc_write(fnm,"dCvolc_dt  ",  vars%dCvolc_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="volcanic CO2 degassing",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dCvolc_cum ",  vars%dCvolc_cum , dim1=dim_time,start=[ndat],count=[y],long_name="cumulative volcanic CO2 degassing",units="GtC",ncid=ncid)  
    call nc_write(fnm,"dC13volc_dt",  vars%dC13volc_dt, dim1=dim_time,start=[ndat],count=[y],long_name="volcaninc CO2 13 degassing",units="GtC/yr",ncid=ncid)  
    call nc_write(fnm,"dC14dec_dt",  vars%dC14dec_dt, dim1=dim_time,start=[ndat],count=[y],long_name="C14 decay",units="kgC/yr",ncid=ncid)  
    call nc_write(fnm,"dC14prod_dt",  vars%dC14prod_dt, dim1=dim_time,start=[ndat],count=[y],long_name="C14 production rate",units="kgC/yr",ncid=ncid)  
    call nc_write(fnm,"dCH4_dt  ",  vars%dCH4_dt  , dim1=dim_time,start=[ndat],count=[y],long_name="carbon flux to atmosphere due to oxidation of anthropogenic CH4",units="GtC/yr",ncid=ncid)  
    call nc_close(ncid)

   return

  end subroutine ts_nc_write

end module co2_out

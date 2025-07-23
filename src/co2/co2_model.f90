!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o 2 _ m o d e l
!
!  Purpose : atmospheric CO2 model
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
module co2_model

    use precision, only : wp
    use control, only: co2_restart, restart_in_dir, co2_const, d13c_atm_const, D14C_atm_const, l_c13, l_c14, l_ocn_co2, l_lnd_co2
    use control, only : iC14_production, C14_production_file, C14_production_const
    use control, only : ico2_degas, co2_degas_file, co2_degas_const, d13c_degas
    use control, only : ico2_emis, co2_emis_file, co2_emis_const, co2_pulse, k_emis_fb, C_emis_fb, co2_emis_min
    use control, only : id13C_emis, d13C_emis_file, d13C_emis_const
    use timer, only : sec_year, year
    use constants, only : c13_c12_std, c14_c_std, c14_tdec
    use co2_def, only : co2_class

    use ncio

    implicit none

    ! conversion factor from kgC to ppm CO2
    real(wp), parameter :: kgC_to_ppm = 1._wp/2.12_wp*1e-12_wp  ! ppm/PgC * PgC/kgC 
    ! http://how-it-looks.blogspot.com/2011/07/petagrams-of-carbon.html
    ! The mean mass of the atmosphere is 5.1480e18 kg. The molar mass of air is 28.966 g/mol. 
    ! So the atmosphere contains 5.1480e21 g / 28.966 g/mol = 1.7773e20 moles of air.
    ! Mole fractions of carbon dioxide are expressed in ppm and directly convertible from parts-per-million by volume.
    ! So 1 ppm CO2 = 1.7773e20 / 1e6 = 1.7773e14 moles CO2
    ! For every mole of carbon dioxide, there is one mole of carbon.  The molar mass of carbon is 12.01 g/mol.
    ! So 1 ppm carbon dioxide = 1.7773e14 mol * 12.01 g/mol = 2.134 Pg of carbon.
    ! 2.12 PgC per ppm used in IPCC (Prather et al., 2012)

    ! conversion factor from atoms/cm2/s -> kgC14/yr, see also Roth and Joos 2013 (1 atoms/cm2/s = 267 mol/yr)
    real(wp), parameter :: conv_prod = sec_year*510.1e12_wp*1.e4_wp / 6.022e23_wp * 12._wp/1000._wp   ! atoms/cm2/s -> kgC14/yr
                                       ! atoms/cm2/s * s/yr * m2 * cm2/m2 * mol/atoms * g/mol * kg/g = kg/yr

    real(wp), dimension(:), allocatable :: C14_production_time, C14_production_data
    real(wp), dimension(:), allocatable :: co2_degas_time, co2_degas_data
    real(wp), dimension(:), allocatable :: co2_emis_time, co2_emis_data
    real(wp), dimension(:), allocatable :: d13C_emis_time, d13C_emis_data

    private
    public :: co2_init, co2_update, co2_end
    public :: co2_read_restart, co2_write_restart


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o 2 _ u p d a t e
  ! Purpose  :  update atmospheric CO2
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_update(co2,time_now)

    implicit none

    type(co2_class) :: co2
    real(wp), intent(in) :: time_now

    integer :: imin, i0, i1
    real(wp) :: w0, w1
    real(wp) :: dCatm_dt


    ! update C14 production rate if needed
    if (iC14_production.eq.1) then
      ! interpolate from C14_production_data to current year
      if (time_now.le.C14_production_time(lbound(C14_production_time,1))) then
        co2%dC14prod_dt = C14_production_data(lbound(C14_production_data,1)) * conv_prod ! kgC
      elseif (time_now.ge.C14_production_time(ubound(C14_production_time,1))) then
        co2%dC14prod_dt = C14_production_data(ubound(C14_production_data,1)) * conv_prod ! kgC
      else
        imin = minloc(abs(C14_production_time-time_now),1) 
        if (C14_production_time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(C14_production_time(i0)-time_now)/(C14_production_time(i1)-C14_production_time(i0))
        w1 = 1._wp - w0
        co2%dC14prod_dt = (w0*C14_production_data(i0) + w1*C14_production_data(i1)) * conv_prod ! kgC
      endif
    endif

    ! update CO2 degassing rate if needed    
    if (ico2_degas.eq.1 .and. year.eq.1) then
      ! half of silicate weathering rate
      co2%dCvolc_dt = 0.5_wp*co2%weath_sil_avg
    endif
    if (ico2_degas.eq.2 .and. year.eq.1) then
      ! explicitely compute volcanic degassing needed to balance carbon fluxes
      ! to atm, from restart
      co2%dCvolc_dt = co2%dCvolc_dt_eq 
    endif
    if (ico2_degas.eq.3) then
      ! interpolate from co2_degas_data to current year
      if (time_now.le.co2_degas_time(lbound(co2_degas_time,1))) then
        co2%dCvolc_dt = co2_degas_data(lbound(co2_degas_data,1))  ! kgC/yr
      elseif (time_now.ge.co2_degas_time(ubound(co2_degas_time,1))) then
        co2%dCvolc_dt = co2_degas_data(ubound(co2_degas_data,1))  ! kgC/yr
      else
        imin = minloc(abs(co2_degas_time-time_now),1) 
        if (co2_degas_time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(co2_degas_time(i0)-time_now)/(co2_degas_time(i1)-co2_degas_time(i0))
        w1 = 1._wp - w0
        co2%dCvolc_dt = (w0*co2_degas_data(i0) + w1*co2_degas_data(i1)) ! kgC
      endif
    endif
    
    co2%dC13volc_dt = (d13c_degas*1.e-3_wp+1._wp)*c13_c12_std * co2%dCvolc_dt ! kgC

    ! update CO2 emission rate if needed    
    if (ico2_emis.eq.1) then
      ! interpolate from co2_emis_data to current year
      if (time_now.le.co2_emis_time(lbound(co2_emis_time,1))) then
        co2%dCemis_dt = co2_emis_data(lbound(co2_emis_data,1)) ! kgC
      elseif (time_now.ge.co2_emis_time(ubound(co2_emis_time,1))) then
        co2%dCemis_dt = co2_emis_data(ubound(co2_emis_data,1)) ! kgC
      else
        imin = minloc(abs(co2_emis_time-time_now),1) 
        if (co2_emis_time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(co2_emis_time(i0)-time_now)/(co2_emis_time(i1)-co2_emis_time(i0))
        w1 = 1._wp - w0
        co2%dCemis_dt = (w0*co2_emis_data(i0) + w1*co2_emis_data(i1)) ! kgC
      endif
      if (time_now.gt.0._wp) then
        ! apply minimum emission rate for future projections
        co2%dCemis_dt = max(co2_emis_min*1.e12_wp,co2%dCemis_dt)
      endif
    else if (ico2_emis.eq.2) then
      ! CO2 pulse, constant emission rate
      if (time_now.le.0._wp) then
        co2%dCemis_dt = 0._wp
      else
        if (time_now.lt.co2_pulse/10.) then
          co2%dCemis_dt = 10._wp*1.e12_wp     ! kgC/yr
        else
          co2%dCemis_dt = 0._wp
        endif
      endif
    else if (ico2_emis.eq.3) then
      ! CO2 pulse, triangle shape
      if (time_now.le.0._wp) then
        co2%dCemis_dt = 0._wp
      else
        if (time_now.lt.0.5_wp*co2_pulse/10._wp) then
          ! up 
          co2%dCemis_dt = 20._wp*time_now/(0.5_wp*co2_pulse/10._wp)*1.e12_wp     ! kgC/yr
        else if (time_now.lt.co2_pulse/10.) then
          ! down
          co2%dCemis_dt = 20._wp*(co2_pulse/10._wp-time_now)/(0.5_wp*co2_pulse/10._wp)*1.e12_wp     ! kgC/yr
        else
          co2%dCemis_dt = 0._wp
        endif
      endif
    else if (ico2_emis.eq.4) then
      ! ZECMIP, fit emissions to follow 1pctCO2 scenario up to a given cumulative emission
      if (time_now.le.0._wp) then
        co2%dCemis_dt = 0._wp
      else
        ! 1 % increase 
        dCatm_dt = co2%Catm*0.01_wp
        if (co2%Cemis_cum*1e-12_wp.lt.co2_pulse) then
          co2%dCemis_dt = dCatm_dt + co2%dClnd_dt + co2%dCocn_dt 
        else
          co2%dCemis_dt = 0._wp
        endif
      endif
    else if (ico2_emis.eq.5) then
      ! negative co2 emissions followed by positive emissions 
      if (time_now.lt.1000._wp) then
        co2%dCemis_dt = -0.2_wp*1.e12_wp     ! kgC/yr
      else if (time_now.ge.1000._wp .and. time_now.le.2000._wp) then
        co2%dCemis_dt = 0.2_wp*1.e12_wp     ! kgC/yr
      endif
    else if (ico2_emis.eq.6) then
      ! negative co2 emissions followed by positive emissions 
      if (time_now.lt.2000._wp) then
        co2%dCemis_dt = -0.5_wp*1.e12_wp     ! kgC/yr
      else if (time_now.ge.2000._wp .and. time_now.le.4000._wp) then
        co2%dCemis_dt = 0.5_wp*1.e12_wp     ! kgC/yr
      endif
    endif

    ! cumulated emissions
    co2%Cemis_cum = co2%Cemis_cum + co2%dCemis_dt   ! GtC

    ! additional emission term based on global temperature change
    co2%dCemis_extra_dt = (C_emis_fb-co2%Cemis_extra_cum)/C_emis_fb * k_emis_fb*co2%dT_glob_cum * 1.e12_wp  ! PgC/yr/K * K * kgC/PgC = kgC/yr
    co2%Cemis_extra_cum = co2%Cemis_extra_cum + co2%dCemis_extra_dt*1.e-12_wp  ! PgC

    ! update d13C of CO2 emissions if needed    
    if (id13C_emis.eq.1) then
      ! interpolate from d13c_emis_data to current year
      if (time_now.le.d13c_emis_time(lbound(d13c_emis_time,1))) then
        co2%d13C_emis = d13c_emis_data(lbound(d13c_emis_data,1)) ! kgC
      elseif (time_now.ge.d13c_emis_time(ubound(d13c_emis_time,1))) then
        co2%d13C_emis = d13c_emis_data(ubound(d13c_emis_data,1)) ! kgC
      else
        imin = minloc(abs(d13c_emis_time-time_now),1) 
        if (d13c_emis_time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(d13c_emis_time(i0)-time_now)/(d13c_emis_time(i1)-d13c_emis_time(i0))
        w1 = 1._wp - w0
        co2%d13C_emis = (w0*d13c_emis_data(i0) + w1*d13c_emis_data(i1)) ! kgC
      endif
    endif

    co2%dC13emis_dt = (co2%d13c_emis*1.e-3_wp+1._wp)*c13_c12_std * co2%dCemis_dt ! kgC

    ! option to exclude carbon flux from ocn
    if (.not.l_ocn_co2) then
      co2%dCocn_dt = 0._wp
      co2%dC13ocn_dt = 0._wp
      co2%dC14ocn_dt = 0._wp
    endif
    ! option to exclude carbon flux from lnd
    if (.not.l_lnd_co2) then
      co2%dClnd_dt = 0._wp
      co2%dC13lnd_dt = 0._wp
      co2%dC14lnd_dt = 0._wp
    endif

    ! radiocarbon decay
    if (l_c14) then
      co2%dC14dec_dt = c14_tdec*sec_year*co2%C14atm
    endif

    ! add fluxes from ocean, land, weathering, volcanic degassing, emissions, CH4 oxidation, production and decay
    co2%Catm   = co2%Catm - co2%dCocn_dt - co2%dClnd_dt + co2%dCemis_dt + co2%dCH4_dt + co2%dCvolc_dt - co2%dCweath_dt + co2%dCemis_extra_dt ! kgC
    if (l_c13) co2%C13atm = co2%C13atm - co2%dC13ocn_dt - co2%dC13lnd_dt + co2%dC13emis_dt + co2%dC13volc_dt - co2%dC13weath_dt    ! kgC13
    if (l_c14) co2%C14atm = co2%C14atm - co2%dC14ocn_dt - co2%dC14lnd_dt - co2%dC14weath_dt - co2%dC14dec_dt + co2%dC14prod_dt

    ! isotopic ratios
    co2%c13_c12 = co2%C13atm / co2%Catm
    co2%c14_c   = co2%C14atm / co2%Catm

    ! CO2 concentration (ppm)
    co2%co2 = co2%Catm * kgC_to_ppm

   return

  end subroutine co2_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o 2 _ i n i t
  ! Purpose  :  initialize atmospheric CO2
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_init(co2)

    implicit none

    type(co2_class) :: co2

    integer :: ntime


    if (co2_restart) then

      ! read restart file
      call co2_read_restart(trim(restart_in_dir)//"/co2_restart.nc",co2)

      ! derive CO2 concentration
      co2%co2 = co2%Catm * kgC_to_ppm
      ! isotopic ratios
      co2%c13_c12 = co2%C13atm / co2%Catm
      co2%c14_c   = co2%C14atm / co2%Catm

    else

      ! initialise
      co2%co2 = co2_const  ! ppm
      co2%Catm = co2%co2 / kgC_to_ppm  ! kgC
      ! isotopes
      co2%c13_c12 = (d13c_atm_const/1000._wp+1._wp)*c13_c12_std
      co2%c14_c   = (D14c_atm_const/1000._wp+1._wp)*c14_c_std / ((0.975_wp/(1._wp+d13c_atm_const/1000._wp))**2)
      co2%C13atm  = co2%Catm*co2%c13_c12
      co2%C14atm  = co2%Catm*co2%c14_c

    endif

    ! C14 production
    if (iC14_production.eq.0) then

      co2%dC14prod_dt = C14_production_const * conv_prod  ! kgC/yr

    else if (iC14_production.eq.1) then

      ! read C14 production rate file
      ntime = nc_size(trim(C14_production_file),"time")
      allocate( C14_production_time(ntime) )
      allocate( C14_production_data(ntime) )
      call nc_read(trim(C14_production_file),"time",C14_production_time)    ! time in years BP
      call nc_read(trim(C14_production_file),"C14_production",C14_production_data) ! atoms/cm2/s

    endif

    ! volcanic degassing
    if (ico2_degas.eq.0) then

      ! constant volcanic CO2 degassing
      co2%dCvolc_dt = co2_degas_const * 1.e12_wp ! kgC/yr

    else if (ico2_degas.eq.1) then

      co2%dCvolc_dt = 0._wp

    else if (ico2_degas.eq.2) then

      co2%dCvolc_dt = 0._wp

    else if (ico2_degas.eq.3) then

      stop 'ico2_degas==3 not supported yet'

      ! read CO2 degassing rate file
      ntime = nc_size(trim(co2_degas_file),"time")
      allocate( co2_degas_time(ntime) )
      allocate( co2_degas_data(ntime) )
      call nc_read(trim(co2_degas_file),"time",co2_degas_time)    ! time in years BP
      call nc_read(trim(co2_degas_file),"co2_degas",co2_degas_data) ! PgC/yr
      co2_degas_data = co2_degas_data * 1.e12 ! kgC/yr

    endif

    ! emissions
    if (ico2_emis.eq.0) then

      ! constant CO2 emissions
      co2%dCemis_dt = co2_emis_const * 1.e12_wp ! kgC/yr

    else if (ico2_emis.eq.1) then

      ! read CO2 emissions file
      ntime = nc_size(trim(co2_emis_file),"time")
      allocate( co2_emis_time(ntime) )
      allocate( co2_emis_data(ntime) )
      call nc_read(trim(co2_emis_file),"time",co2_emis_time)    ! time in years BP
      call nc_read(trim(co2_emis_file),"co2_emis",co2_emis_data) ! PgC/yr
      co2_emis_data = co2_emis_data * 1.e12_wp  ! kgC/yr

    endif

    co2%Cemis_cum = 0._wp

    ! d13C of emissions
    if (id13C_emis.eq.0) then

      ! constant d13C of CO2 emissions
      co2%d13C_emis = d13C_emis_const   ! permil

    else if (id13C_emis.eq.1) then

      ! read d13C of CO2 emissions file
      ntime = nc_size(trim(d13C_emis_file),"time")
      allocate( d13c_emis_time(ntime) )
      allocate( d13c_emis_data(ntime) )
      call nc_read(trim(d13c_emis_file),"time",d13c_emis_time)    ! time in years BP
      call nc_read(trim(d13c_emis_file),"d13C_emis",d13c_emis_data) ! permil

    endif

    ! CH4 oxidation
    co2%dCH4_dt = 0._wp

    co2%T_glob = 0._wp
    co2%dT_glob_cum = 0._wp
    co2%Cemis_extra_cum = 0._wp

  return

  end subroutine co2_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o 2 _ e n d 
  ! Purpose  :  end atmospheric CO2
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_end(co2)

    implicit none

    type(co2_class) :: co2



   return

  end subroutine co2_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o 2 _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_write_restart(fnm,co2)

    implicit none

    character (len=*) :: fnm
    type(co2_class) :: co2


    call nc_create(fnm)
    call nc_write_dim(fnm,"x",x=[1])
    call nc_write(fnm,"Catm",   co2%Catm,     dim1="x",  long_name="atmospheric carbon",units="kgC")
    call nc_write(fnm,"C13atm", co2%C13atm,   dim1="x",  long_name="atmospheric carbon 13",units="kgC")
    call nc_write(fnm,"C14atm", co2%C14atm,   dim1="x",  long_name="atmospheric carbon 14",units="kgC")

   return

  end subroutine co2_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c o 2 _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine co2_read_restart(fnm,co2)

    implicit none

    character (len=*) :: fnm
    type(co2_class) :: co2


    call nc_read(fnm,"Catm",     co2%Catm)
    call nc_read(fnm,"C13atm",   co2%C13atm)
    call nc_read(fnm,"C14atm",   co2%C14atm)

   return

  end subroutine co2_read_restart


end module co2_model


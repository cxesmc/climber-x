!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : n 2 o _ m o d e l
!
!  Purpose : atmospheric N2O model
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
module n2o_model

    use precision, only : wp
    use control, only: n2o_restart, restart_in_dir, n2o_const, i_n2o_tau, n2o_tau_const, n2o_ref
    use control, only : in2o_emis, n2o_emis_file, n2o_emis_const
    use timer, only : sec_year
    use n2o_def, only : n2o_class

    use ncio

    implicit none

    ! conversion factor from kgN2O-N to ppb N2O
    real(wp), parameter :: kgN_to_ppb = 1._wp/4.976_wp*1e-9_wp  ! ppb/TgN2O-N * TgN2O-N/kgN2O-N
    ! The mean mass of the atmosphere is 5.1480e18 kg. The molar mass of air is 28.966 g/mol. 
    ! So the atmosphere contains 5.1480e21 g / 28.966 g/mol = 1.7773e20 moles of air.
    ! Mole fractions of carbon dioxide are expressed in ppm and directly convertible from parts-per-million by volume.
    ! So 1 ppb N2O = 1.7773e20 / 1e9 = 1.7773e11 moles N2O.
    ! The molar mass of N is 14.0 g/mol.
    ! So 1 ppb N2O = 1.7773e11 mol * 2*14.0 g/mol = 4.976 Tg N2O-N.
    
    real(wp), dimension(:), allocatable :: n2o_emis_time, n2o_emis_data, f_n2o_emis_agro_data, n2o_tau_time, n2o_tau_data

    private
    public :: n2o_init, n2o_update, n2o_end
    public :: n2o_read_restart, n2o_write_restart


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  n 2 o _ u p d a t e
  ! Purpose  :  update atmospheric n2o
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_update(n2o,time_now)

    implicit none

    type(n2o_class) :: n2o
    real(wp), intent(in) :: time_now

    integer :: imin, i0, i1
    real(wp) :: w0, w1
    real(wp) :: n2o_tau
    real(wp) :: emi_CO, emi_NOx, emi_VOC, tau_OH_0, d_OH, d_OH_rel, tau_OH


    ! update n2o emission rate if needed    
    if (in2o_emis.eq.1) then
      ! interpolate from n2o_emis_data to current year
      if (time_now.le.n2o_emis_time(lbound(n2o_emis_time,1))) then
        n2o%dn2oemis_dt = n2o_emis_data(lbound(n2o_emis_data,1)) ! kgN2O
      elseif (time_now.ge.n2o_emis_time(ubound(n2o_emis_time,1))) then
        n2o%dn2oemis_dt = n2o_emis_data(ubound(n2o_emis_data,1)) ! kgN2O
      else
        imin = minloc(abs(n2o_emis_time-time_now),1) 
        if (n2o_emis_time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        w0 = 1._wp - abs(n2o_emis_time(i0)-time_now)/(n2o_emis_time(i1)-n2o_emis_time(i0))
        w1 = 1._wp - w0
        n2o%dn2oemis_dt = (w0*n2o_emis_data(i0) + w1*n2o_emis_data(i1)) !  kgN2O-N/yr
      endif
    endif

    if (i_n2o_tau.eq.1) then
      ! constant lifetime

      n2o%tau = n2o_tau_const

    endif

    n2o%dn2oox_dt = n2o%n2om/n2o%tau   ! kgN2O-N/yr

    ! add fluxes from ocean, land, anthropogenic emissions and decay (oxidation)
    n2o%n2om   = n2o%n2om   + n2o%dn2oocn_dt   + n2o%dn2olnd_dt   + n2o%dn2oemis_dt - n2o%dn2oox_dt ! kgN2O-N

    ! n2o concentration
    n2o%n2o = n2o%n2om * kgN_to_ppb   ! ppb

   return

  end subroutine n2o_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  n 2 o _ i n i t
  ! Purpose  :  initialize atmospheric n2o
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_init(n2o)

    implicit none

    type(n2o_class) :: n2o

    integer :: ntime


    if (n2o_restart) then

      ! read restart file
      call n2o_read_restart(trim(restart_in_dir)//"/n2o_restart.nc",n2o)

    else

      ! initialise
      n2o%n2o = n2o_const  ! ppb

    endif

    n2o%n2om = n2o%n2o / kgN_to_ppb   ! kgN2O-N

    ! emissions
    if (in2o_emis.eq.0) then

      ! constant n2o emissions
      n2o%dn2oemis_dt = n2o_emis_const * 1.e9_wp ! kgN2O-N/yr

    else if (in2o_emis.eq.1) then

      ! read n2o emissions file
      ntime = nc_size(trim(n2o_emis_file),"time")
      allocate( n2o_emis_time(ntime) )
      allocate( n2o_emis_data(ntime) )
      call nc_read(trim(n2o_emis_file),"time",n2o_emis_time)    ! time in years BP
      call nc_read(trim(n2o_emis_file),"n2o_emis",n2o_emis_data) ! TgN2O-N/yr
      n2o_emis_data = n2o_emis_data * 1.e9_wp  ! kgN2O-N/yr

    endif

  return

  end subroutine n2o_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  n 2 o _ e n d 
  ! Purpose  :  end atmospheric n2o
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_end(n2o)

    implicit none

    type(n2o_class) :: n2o



   return

  end subroutine n2o_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  n 2 o _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_write_restart(fnm,n2o)

    implicit none

    character (len=*) :: fnm
    type(n2o_class) :: n2o


    call nc_create(fnm)
    call nc_write_dim(fnm,"x",x=[1])
    call nc_write(fnm,"n2o",   n2o%n2o,     dim1="x",  long_name="atmospheric N2O concentration",units="ppb")

   return

  end subroutine n2o_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  n 2 o _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine n2o_read_restart(fnm,n2o)

    implicit none

    character (len=*) :: fnm
    type(n2o_class) :: n2o


    call nc_read(fnm,"n2o",     n2o%n2o)

   return

  end subroutine n2o_read_restart


end module n2o_model


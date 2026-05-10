!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o 2 _ m o d
!
!  Purpose : atmospheric CO2 concentration
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
module co2_mod

   use precision, only : wp
   use timer, only : year
   use control, only : co2_const, dco2_dt, co2_max
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time, co2_data

contains

   subroutine co2_init(co2_file)

   implicit none

   character (len=*), intent(in), optional :: co2_file

   integer :: ntime


   ntime = nc_size(trim(co2_file),"time")
   allocate( time(ntime) )
   allocate( co2_data(ntime) )
   call nc_read(trim(co2_file),"time",time)    ! time in years BP
   call nc_read(trim(co2_file),"co2",co2_data) 


   end subroutine co2_init


    subroutine co2_update(ico2, time_now, co2)

    implicit none

    integer, intent(in) :: ico2
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: co2  ! current atmospheric CO2 concentration

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (ico2.eq.1 .or. ico2.eq.-1) then

      ! interpolate from co2_data to current year
      if (time_now.le.minval(time)) then
        co2 = co2_data(lbound(co2_data,1))
      elseif (time_now.ge.maxval(time)) then
        co2 = co2_data(ubound(co2_data,1))
      else
        imin = minloc(abs(time-time_now),1) 
        if (time(imin).lt.time_now) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        if (time(i1)-time(i0).eq.0._wp) then
          w0 = 1._wp
        else
          w0 = 1._wp - abs(time(i0)-time_now)/(time(i1)-time(i0))
        endif
        w1 = 1._wp - w0
        co2 = w0*co2_data(i0) + w1*co2_data(i1)
      endif

    else if (ico2.eq.2) then

      ! exponential increase with rate dco2_dt (in %/yr) up to co2_max
      co2 = co2*(1._wp+dco2_dt/100._wp)
      !co2 = min(1120._wp+10._wp,co2)
      co2 = min(co2_max,co2)

    else if (ico2.eq.3) then

      ! linear CO2 increase
      co2 = co2 + dco2_dt

    else if (ico2.eq.4) then

      ! linear CO2 decrease
      co2 = co2 - dco2_dt

    else if (ico2.eq.5) then

      ! linear CO2 increase
      co2 = co2 + 0.003

    else if (ico2.eq.6) then

      ! linear CO2 decrease
      co2 = co2 - 0.003

    else if (ico2.eq.7) then

      ! linear CO2 increase
      co2 = co2 + 0.005

    else if (ico2.eq.8) then

      ! linear CO2 decrease
      co2 = co2 - 0.005

    else if (ico2.eq.9) then

      ! linear CO2 ramp up
      if (year.le.5000) then
        co2 = co2 + (co2_const-280._wp)/5000._wp
      endif

    endif


   return

  end subroutine co2_update

end module co2_mod

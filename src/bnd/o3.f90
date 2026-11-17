!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o 3 _ m o d
!
!  Purpose : atmospheric O3 concentration
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
module o3_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time
   real(wp), dimension(:,:,:,:), allocatable :: o3_data
   real(wp), dimension(:), allocatable :: o3_pl    !! relative pressure (p/1000 hPa) of the O3 levels, surface first

contains

   subroutine o3_init(io3, o3_file)

   implicit none

   integer, intent(in) :: io3
   character (len=*), intent(in), optional :: o3_file

   integer :: ntime, nlon, nlat, nlev
   real(wp), dimension(:), allocatable :: lev
   character (len=64) :: lev_units


   if (io3.eq.1 .or. io3.eq.2) then

     nlon  = nc_size(trim(o3_file),"lon")
     nlat  = nc_size(trim(o3_file),"lat")
     nlev  = nc_size(trim(o3_file),"lev")
     if (io3.eq.1) then
       ntime = 1
     else
       ntime = nc_size(trim(o3_file),"time")
       allocate( time(ntime) )
       call nc_read(trim(o3_file),"time",time)    ! time in years BP
     endif
     allocate( o3_data(nlon,nlat,nlev,ntime) )
     call nc_read(trim(o3_file),"o3",o3_data) 

     ! vertical levels of the O3 field from the file, as relative pressure p/p0 (p0 = 1000 hPa)
     allocate( lev(nlev) )
     allocate( o3_pl(nlev) )
     call nc_read(trim(o3_file),"lev",lev)
     call nc_read_attr(trim(o3_file),"lev","units",lev_units)
     select case (trim(lev_units))
     case ("hPa", "mbar", "millibar")
       o3_pl = lev/1000._wp
     case ("Pa")
       o3_pl = lev/1.e5_wp
     case default
       print *,'ERROR: unknown units of lev in O3 file ',trim(o3_file),': ',trim(lev_units)
       stop
     end select

     ! surface first, as expected by the interpolation to the atmospheric levels in the coupler
     if (nlev.gt.1 .and. o3_pl(1).lt.o3_pl(nlev)) then
       o3_pl = o3_pl(nlev:1:-1)
       o3_data = o3_data(:,:,nlev:1:-1,:)
     endif
     if (any(o3_pl(2:nlev).ge.o3_pl(1:nlev-1))) then
       print *,'ERROR: O3 levels in ',trim(o3_file),' are not strictly monotonic'
       stop
     endif

   endif


   end subroutine o3_init


    subroutine o3_update(io3, time_now, o3)

    implicit none

    integer, intent(in) :: io3
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: o3(:,:,:)  ! current atmospheric o3 concentration

    integer :: imin, i0, i1
    real(wp) :: w0, w1


    if (io3.eq.1) then

        o3(:,:,:) = o3_data(:,:,:,1)

    else if (io3.eq.2) then

      ! interpolate from o3_data to current year
      if (time_now.le.minval(time)) then
        o3(:,:,:) = o3_data(:,:,:,lbound(o3_data,4))
      elseif (time_now.ge.maxval(time)) then
        o3(:,:,:) = o3_data(:,:,:,ubound(o3_data,4))
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
        o3(:,:,:) = w0*o3_data(:,:,:,i0) + w1*o3_data(:,:,:,i1)
      endif

    endif


   return

  end subroutine o3_update

end module o3_mod

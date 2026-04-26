!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l u c _ m o d
!
!  Purpose : land use change
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
module luc_mod

   use precision, only : wp
   use ncio

   implicit none

   real(wp), dimension(:), allocatable :: time
   real(wp), dimension(:,:,:), allocatable :: f_crop_data
   real(wp), dimension(:,:,:), allocatable :: f_pasture_data

contains

   subroutine luc_init(iluc, luc_file)

   implicit none

   integer, intent(in) :: iluc
   character (len=*), intent(in) :: luc_file

   integer :: ntime, nlon, nlat


   if (iluc.eq.1) then

     nlon  = nc_size(trim(luc_file),"lon")
     nlat  = nc_size(trim(luc_file),"lat")
     allocate( f_crop_data(nlon,nlat,1) )
     allocate( f_pasture_data(nlon,nlat,1) )
     call nc_read(trim(luc_file),"crops",f_crop_data) 
     call nc_read(trim(luc_file),"pasture",f_pasture_data) 

   else if (iluc.eq.2) then

     ntime = nc_size(trim(luc_file),"time")
     nlon  = nc_size(trim(luc_file),"lon")
     nlat  = nc_size(trim(luc_file),"lat")
     allocate( time(ntime) )
     allocate( f_crop_data(nlon,nlat,ntime) )
     allocate( f_pasture_data(nlon,nlat,ntime) )
     call nc_read(trim(luc_file),"time",time)    ! time in years BP
     call nc_read(trim(luc_file),"crops",f_crop_data) 
     call nc_read(trim(luc_file),"pasture",f_pasture_data) 

   endif

   end subroutine luc_init


    subroutine luc_update(iluc, time_now, f_crop, f_pasture, df_crop, df_pasture)

    implicit none

    integer, intent(in) :: iluc
    real(wp), intent(in) :: time_now       ! current year BP
    real(wp), intent(inout) :: f_crop(:,:)  ! current crop fraction
    real(wp), intent(inout) :: f_pasture(:,:)  ! current pasture fraction
    real(wp), intent(inout) :: df_crop(:,:)  ! change in crop fraction
    real(wp), intent(inout) :: df_pasture(:,:)  ! change in pasture fraction

    integer :: imin, i0, i1
    real(wp) :: w0, w1
    real(wp), allocatable :: f_crop_old(:,:)
    real(wp), allocatable :: f_pasture_old(:,:)


    ! save fields from previous step
    f_crop_old(:,:)    = f_crop(:,:)
    f_pasture_old(:,:) = f_pasture(:,:)

    if (iluc.eq.1) then

      f_crop(:,:)    = f_crop_data(:,:,1)
      f_pasture(:,:) = f_pasture_data(:,:,1)

    else if (iluc.eq.2) then

      ! interpolate to current year
      if (time_now.le.minval(time)) then
        f_crop(:,:)    = f_crop_data(:,:,lbound(f_crop_data,3))
        f_pasture(:,:) = f_pasture_data(:,:,lbound(f_pasture_data,3))
      elseif (time_now.ge.maxval(time)) then
        f_crop(:,:)    = f_crop_data(:,:,ubound(f_crop_data,3))
        f_pasture(:,:) = f_pasture_data(:,:,ubound(f_pasture_data,3))
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
        f_crop(:,:)    = w0*f_crop_data(:,:,i0)    + w1*f_crop_data(:,:,i1)
        f_pasture(:,:) = w0*f_pasture_data(:,:,i0) + w1*f_pasture_data(:,:,i1)
      endif

    endif

    ! compute changes in crop and pasture fraction
    df_crop(:,:) = f_crop(:,:) - f_crop_old(:,:)
    df_pasture(:,:) = f_pasture(:,:) - f_pasture_old(:,:)

   return

  end subroutine luc_update

end module luc_mod

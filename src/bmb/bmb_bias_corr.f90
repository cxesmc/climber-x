!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b i a s _ c o r r _ m o d 
!
!  Purpose : climate bias correction
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
module bmb_bias_corr_mod

   use precision, only : wp
   use bmb_params, only : bias_corr_file, Tocn_bias_scale_fac
   use bmb_params, only : Tocn_offset, Tocn_offset_file  
   use ncio

   implicit none

   real(wp), allocatable, dimension(:) :: time_bias
   real(wp), allocatable, dimension(:,:) :: Tocn_bias_data_const
   real(wp), allocatable, dimension(:,:,:) :: Tocn_bias_data_var

   real(wp), allocatable, dimension(:) :: time_offset
   real(wp), allocatable, dimension(:) :: Tocn_offset_data

   private
   public :: Tocn_bias_corr_init, Tocn_bias_corr_update
   public :: Tocn_offset_init, Tocn_offset_update

contains


   subroutine Tocn_bias_corr_init(i_Tocn_bias_corr)

   implicit none

   integer, intent(in) :: i_Tocn_bias_corr

   integer :: nx, ny, ntime


   nx = nc_size(trim(bias_corr_file),"lon")
   ny = nc_size(trim(bias_corr_file),"lat")

   if (i_Tocn_bias_corr.eq.1) then

     allocate(Tocn_bias_data_const(nx,ny))
     ! read temperature bias file
     call nc_read(trim(bias_corr_file),"Tocn_bias",Tocn_bias_data_const) 

   else if (i_Tocn_bias_corr.eq.2) then

     ntime = nc_size(trim(bias_corr_file),"time")
     allocate(Tocn_bias_data_var(nx,ny,ntime))
     ! read temperature bias file
     call nc_read(trim(bias_corr_file),"Tocn_bias",Tocn_bias_data_var) 

   endif

   return

   end subroutine Tocn_bias_corr_init

   
   subroutine Tocn_bias_corr_update(i_Tocn_bias_corr, time, Tocn_bias)

   implicit none

   integer, intent(in) :: i_Tocn_bias_corr
   real(wp), intent(in) :: time
   real(wp), intent(out) :: Tocn_bias(:,:)

   integer :: imin, i0, i1
   real(wp) :: w0, w1


   if (i_Tocn_bias_corr.eq.1) then

     Tocn_bias = Tocn_bias_data_const

   else if (i_Tocn_bias_corr.eq.2) then

      ! interpolate to current year
      if (time.le.minval(time_bias)) then
        Tocn_bias(:,:) = Tocn_bias_data_var(:,:,lbound(time_bias,1))
      elseif (time.ge.maxval(time_bias)) then
        Tocn_bias(:,:) = Tocn_bias_data_var(:,:,ubound(time_bias,1))
      else
        imin = minloc(abs(time_bias-time),1) 
        if (time_bias(imin).lt.time) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        if (time_bias(i1)-time_bias(i0).eq.0._wp) then
          w0 = 1._wp
        else
          w0 = 1._wp - abs(time_bias(i0)-time)/(time_bias(i1)-time_bias(i0))
        endif
        w1 = 1._wp - w0
        Tocn_bias(:,:) = w0*Tocn_bias_data_var(:,:,i0) + w1*Tocn_bias_data_var(:,:,i1)
      endif

   endif

   ! apply scaling factor 
   Tocn_bias = Tocn_bias_scale_fac*Tocn_bias 

   return

   end subroutine Tocn_bias_corr_update


   subroutine Tocn_offset_init(i_Tocn_offset)

   implicit none

   integer, intent(in) :: i_Tocn_offset

   integer :: ntime


   if (i_Tocn_offset.eq.2) then

     ntime = nc_size(trim(Tocn_offset_file),"time")
     allocate( time_offset(ntime) )
     allocate( Tocn_offset_data(ntime) )
     call nc_read(trim(Tocn_offset_file),"time",time_offset)    ! time in years BP
     call nc_read(trim(Tocn_offset_file),"Tocn_offset",Tocn_offset_data) 

   endif

   end subroutine Tocn_offset_init

   
   subroutine Tocn_offset_update(i_Tocn_offset, time, Tocn_offset_now)

   implicit none

   integer, intent(in) :: i_Tocn_offset
   real(wp), intent(in) :: time
   real(wp), intent(out) :: Tocn_offset_now

   integer :: imin, i0, i1
   real(wp) :: w0, w1


   if (i_Tocn_offset.eq.1) then

     Tocn_offset_now = Tocn_offset

   else if (i_Tocn_offset.eq.2) then

      ! interpolate to current year
      if (time.le.minval(time_offset)) then
        Tocn_offset_now = Tocn_offset_data(lbound(time_offset,1))
      elseif (time.ge.maxval(time_offset)) then
        Tocn_offset_now = Tocn_offset_data(ubound(time_offset,1))
      else
        imin = minloc(abs(time_offset-time),1) 
        if (time_offset(imin).lt.time) then
          i0 = imin
          i1 = imin+1
        else
          i0 = imin-1
          i1 = imin
        endif
        if (time_offset(i1)-time_offset(i0).eq.0._wp) then
          w0 = 1._wp
        else
          w0 = 1._wp - abs(time_offset(i0)-time)/(time_offset(i1)-time_offset(i0))
        endif
        w1 = 1._wp - w0
        Tocn_offset_now = w0*Tocn_offset_data(i0) + w1*Tocn_offset_data(i1)
      endif

   endif


   return

   end subroutine Tocn_offset_update

end module bmb_bias_corr_mod

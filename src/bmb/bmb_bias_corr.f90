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
   use ncio

   implicit none

   private
   public :: Tocn_bias_corr

contains


   subroutine Tocn_bias_corr(i_Tocn_bias_corr, bias_corr_file, Tocn_bias)

   implicit none

   integer, intent(in) :: i_Tocn_bias_corr
   character (len=*), intent(in) :: bias_corr_file

   real(wp), intent(out) :: Tocn_bias(:,:)


   if (i_Tocn_bias_corr.eq.1) then

     ! read temperature bias file
     call nc_read(trim(bias_corr_file),"Tocn_bias",Tocn_bias) 

   endif

   return

   end subroutine Tocn_bias_corr

end module bmb_bias_corr_mod

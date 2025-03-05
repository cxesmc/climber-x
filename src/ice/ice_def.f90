!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i c e _ d e f
!
!  Purpose : definition of ice sheet model class 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module ice_def

    use precision, only : wp
    use coord, only : grid_class

    implicit none
    
    type grid_ice_to_cmn_type
      integer, dimension(:,:), allocatable :: i_lowres 
      integer, dimension(:,:), allocatable :: j_lowres
      integer, dimension(:,:), allocatable :: ncells
    end type

    type ice_class
        ! Variables on ice sheet grid that interact with the coupler 
        ! All information will be on the ice-sheet native grid 

        ! Ice grid information 
        type(grid_class)      :: grid 

        ! ice <-> cmn grid correspondence
        type(grid_ice_to_cmn_type) :: grid_ice_to_cmn

        ! out
        logical :: error
        integer,  allocatable :: mask_extent(:,:) ! should be external in bnd?  fixme
        real(wp), allocatable :: z_sur(:,:) 
        real(wp), allocatable :: z_sur_std(:,:) 
        real(wp), allocatable :: z_bed_std(:,:) 
        real(wp), allocatable :: z_base(:,:) 
        real(wp), allocatable :: H_ice(:,:) 
        real(wp), allocatable :: calv(:,:) 
        real(wp), allocatable :: Q_b(:,:) 

        ! in
        integer, allocatable :: mask_ocn_lake(:,:)
        real(wp), allocatable :: z_sl(:,:)
        real(wp), allocatable :: z_bed(:,:) 
        real(wp), allocatable :: z_bed_fil(:,:) 
        real(wp), allocatable :: smb(:,:) 
        real(wp), allocatable :: accum(:,:) 
        real(wp), allocatable :: runoff(:,:) 
        real(wp), allocatable :: bmb_shlf(:,:) 
        real(wp), allocatable :: temp_s(:,:) 
        real(wp), allocatable :: temp_g(:,:) 
        real(wp), allocatable :: q_geo(:,:) 
        real(wp), allocatable :: H_sed(:,:) 
        real(wp), allocatable :: t_ocn(:,:) 
        real(wp), allocatable :: s_ocn(:,:) 
        
    end type 

    private
    public :: ice_class

end module

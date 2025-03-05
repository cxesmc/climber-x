!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b m b _ d e f
!
!  Purpose : definition of basal mass balance model class 
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
module bmb_def

    use precision, only : wp
    use coord, only : grid_class
    use coord, only : map_class, map_scrip_class
    use timer, only : nmon_year, nday_year

    implicit none
    
    type ts_out
      real(wp) :: bmb
    end type

    type s_out
      ! input variables on coupler (low) resolution
      real(wp), dimension(:,:,:), allocatable :: t_ocn_in  !! ocean temperature on coupler grid [K]
      real(wp), dimension(:,:,:), allocatable :: s_ocn_in  !! ocean salinity on coupler grid [psu]
      ! input variables from ice sheet (matches bmb grid)
      integer, dimension(:,:), allocatable :: mask_ocn_lake   !! lake mask on ice sheet grid []
      real(wp), dimension(:,:), allocatable :: mask_ice_shelf   !! ice mask on ice sheet grid []
      real(wp), dimension(:,:), allocatable :: zb   !! ice base elevation on ice sheet grid []
      ! interpolated variables on ice sheet grid
      real(wp), dimension(:,:,:), allocatable :: t_ocn  !! ocean temperature interpolated to ice sheet grid [K]
      real(wp), dimension(:,:,:), allocatable :: s_ocn  !! ocean salinity interpolated to ice sheet grid [psu]
      ! bmb variables
      real(wp), dimension(:,:), allocatable :: t_bmb   !! temperature used for basal mass balance [K]
      real(wp), dimension(:,:), allocatable :: s_bmb   !! salinity used for basal mass balance [psu]
      real(wp), dimension(:,:), allocatable :: t_disc   !! temperature used for small-scale basal melt [K]
      real(wp), dimension(:,:), allocatable :: s_disc   !! salinity used for small-scale basal melt [psu]
      real(wp), dimension(:,:), allocatable :: t_freeze   !! freezing point temperature [K]
      real(wp), dimension(:,:), allocatable :: bmb   !! basal mass balance of floating ice [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: bmb_ann   !! annual integral of basal mass balance of floating ice [kg/m2/a]
      real(wp), dimension(:,:), allocatable :: bmb_ann_mask   !! annual integral of basal mass balance of floating ice only for shelf mask [kg/m2/a]
    end type

    type bmb_class

      ! grid definition
      type(grid_class) :: grid
      type(grid_class) :: grid_in
      type(map_class) :: map_cmn_to_ice
      type(map_scrip_class) :: maps_cmn_to_ice

      integer :: ncells                                        !! total number of active bmb columns
      integer, dimension(:), allocatable :: idx_cell_active    !! index of active bmb cells
      integer, dimension(:, :), allocatable :: ij_1d           !! n --> i,j  (2 x ncol)
      integer, dimension(:, :), allocatable :: id_map          !! i,j --> n  (ni, nj)

      ! input variables on coupler (low) resolution
      integer, dimension(:,:,:), allocatable :: mask_ocn_in  !! ocean mask on coupler grid [/]
      integer, dimension(:,:), allocatable :: mask_ocn_in_tmp  !! ocean mask on coupler grid [/]
      real(wp), dimension(:), allocatable :: z_ocn_in  !! depth of ocean layers for temperature and salinity input [m]
      real(wp), dimension(:,:,:), allocatable :: t_ocn_in  !! ocean temperature on coupler grid [K]
      real(wp), dimension(:,:,:), allocatable :: s_ocn_in  !! ocean salinity on coupler grid [psu]
      integer, dimension(:,:), allocatable :: mask_lake_in  !! lake mask on coupler grid [/]
      integer, dimension(:,:), allocatable :: mask_lake_in_tmp  !! lake mask on coupler grid [/]
      real(wp), dimension(:), allocatable :: z_lake_in  !! depth of lake layers for temperature and salinity input [m]
      real(wp), dimension(:,:,:), allocatable :: t_lake_in  !! lake temperature on coupler grid [K]
      real(wp), dimension(:,:,:), allocatable :: s_lake_in  !! lake salinity on coupler grid [psu]
      ! input variables from ice sheet (matches bmb grid)
      integer, dimension(:,:), allocatable :: mask_ocn_lake  !! lake mask on ice sheet grid [/]
      integer, dimension(:,:), allocatable :: mask_ice_shelf   !! ice mask on ice sheet grid []
      real(wp), dimension(:,:), allocatable :: zb   !! ice base elevation on ice sheet grid []
      real(wp), dimension(:,:), allocatable :: zl_fil   !! filtered lithosphere elevation on ice sheet grid []
      ! interpolated variables on ice sheet grid
      real(wp), dimension(:,:,:), allocatable :: t_ocn  !! ocean temperature interpolated to ice sheet grid [K]
      real(wp), dimension(:,:,:), allocatable :: s_ocn  !! ocean salinity interpolated to ice sheet grid [psu]
      real(wp), dimension(:,:,:), allocatable :: t_lake  !! lake temperature interpolated to ice sheet grid [K]
      real(wp), dimension(:,:,:), allocatable :: s_lake  !! lake salinity interpolated to ice sheet grid [psu]
      ! bmb variables
      real(wp), dimension(:,:), allocatable :: t_bmb   !! temperature used for basal mass balance [K]
      real(wp), dimension(:,:), allocatable :: s_bmb   !! salinity used for basal mass balance [psu]
      real(wp), dimension(:,:), allocatable :: t_disc   !! temperature used for small-scale basal melt [K]
      real(wp), dimension(:,:), allocatable :: s_disc   !! salinity used for small-scale basal melt [psu]
      real(wp), dimension(:,:), allocatable :: t_freeze   !! freezing point temperature [K]
      real(wp), dimension(:,:), allocatable :: bmb   !! basal mass balance of floating ice [kg/m2/s]
      real(wp), dimension(:,:), allocatable :: bmb_ann   !! annual integral of basal mass balance of floating ice [kg/m2/a]

      ! output types
      integer :: nout
      type(ts_out), allocatable :: ann_ts(:)
      type(ts_out) :: mon_ts(nmon_year)
      type(s_out) :: mon_s(nmon_year), ann_s

    end type

    private
    public :: bmb_class, ts_out, s_out

end module

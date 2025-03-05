!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : g e o _ d e f
!
!  Purpose : definition of geo model class 
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
module geo_def

    use precision, only : wp, dp
    use coord, only : grid_class

    implicit none
    
    type hires_lake_type
      integer :: i, j !! index of lake deepest point position on high resolution
      integer :: i_runoff, j_runoff   !! index of lake runoff destination cells on high resolution
    end type

    integer, parameter :: n_lev = 50

    type lake_type
      type(hires_lake_type) :: hires
      integer :: i, j !! index of lake deepest point position on coarse resolution
      integer :: i_runoff, j_runoff   !! index of lake runoff destination cells on coarse resolution
      real(wp) :: area_min, vol_min
      real(wp) :: area_pot, vol_pot
      real(wp) :: area, vol, vol_old
      real(wp), dimension(0:n_lev) :: z_lev
      real(wp), dimension(0:n_lev) :: vol_lev
      real(wp), dimension(0:n_lev) :: area_lev
      real(wp) :: z
      real(wp) :: depth
      real(wp) :: depth_max
      real(wp) :: z_bot
      real(wp) :: dh_p_e
      real(wp) :: dh_run
      real(wp) :: runoff
      real(wp) :: runoff_geo
      real(wp) :: runoff_diag
      real(wp) :: calving
    end type

    type hires_type
      type(grid_class) :: grid
      real(wp), dimension(:), allocatable :: lon_1min
      real(wp), dimension(:), allocatable :: lat_1min
      real(wp), dimension(:,:), allocatable :: z_bed_ref     !! present day reference bedrock topography [m]
      real(wp), dimension(:,:), allocatable :: h_ice_ref     !! present day ice thickness [m]
      real(wp), dimension(:,:), allocatable :: z_topo_ref     !! present day reference topography [m]
      integer,  dimension(:,:), allocatable :: mask_ref     !! present day reference mask [1]
      real(wp), dimension(:,:), allocatable :: z_topo     !! topography [m]
      real(wp), dimension(:,:), allocatable :: z_topo_fil     !! filtered topography [m]
      real(wp), dimension(:,:), allocatable :: z_topo_fill     !! filled topography [m]
      real(wp), dimension(:,:), allocatable :: z_sur     !! surface elevation [m]
      real(wp), dimension(:,:), allocatable :: h_ice     !! ice thickness [m]
      real(wp), dimension(:,:), allocatable :: rsl     !! relative sea level [m]
      integer, dimension(:,:), allocatable :: mask        !! ice/land/ocean/floating ice mask (following sicopolis convention)
      real(wp), dimension(:,:), allocatable :: z_bed  !! bedrock elevation [m]
      real(wp), dimension(:,:), allocatable :: z_bed_rel  !! isostatically relaxed bedrock elevation [m]
      real(wp), dimension(:,:), allocatable :: z_bed_1min  !! bedrock elevation on 1 min resolution [m]
      real(wp), dimension(:,:), allocatable :: q_geo     !! geothermal heat flux [W/m2]
      real(wp), dimension(:,:), allocatable :: h_sed     !! sediment thickness [m]
      integer, dimension(:,:), allocatable :: i_runoff, j_runoff  !! (i,j) indexes of runoff destination cell
      integer, dimension(:,:), allocatable :: i_runoff_coarse, j_runoff_coarse !! (i,j) indexes of runoff destination cell on the coearse grid
      integer, dimension(:,:), allocatable :: map_runoff      !! map of runoff to ocean/lakes
      real(wp), dimension(:,:), allocatable :: flow_acc       !! flow accumulation [m2]
      integer, dimension(:,:), allocatable :: mask_lake_pot   !! potential lake mask
      integer, dimension(:,:), allocatable :: mask_lake       !! lake mask
      real(wp), dimension(:,:), allocatable :: z_lake         !! lake surface elevation
      real(wp), dimension(:,:), allocatable :: z_sea_lake         !! sea/lake surface elevation
      integer, dimension(:,:), allocatable :: mask_ocn_lake       !! ocean/lake mask
      integer, dimension(:,:), allocatable :: drain_basins_ocn 
    end type hires_type

    type geo_class
      real(wp) :: sea_level  !! sea level [m]
      real(wp) :: d_sea_level      !! change in sea level [m]
      real(wp) :: V_ice_af       !! global ice sheet volume above floatation [m3]
      real(wp) :: A_bering       !! Bering Strait cross-sectional area [m2]      
      type(grid_class) :: grid
      real(wp), dimension(:,:), allocatable :: f_ocn0     !! ocean fraction for present day sea level
      real(wp), dimension(:,:), allocatable :: f_lnd0     !! land fraction for present day sea level (1-ocean fraction)
      real(wp), dimension(:,:), allocatable :: f_ocn      !! ocean fraction including floating ice
      real(wp), dimension(:,:), allocatable :: f_ocn2     !! ocean fraction excluding floating ice
      real(wp), dimension(:,:), allocatable :: f_lnd     !! land fraction (1-ocean fraction)
      real(wp), dimension(:,:), allocatable :: f_lake    !! total lake fraction
      real(wp), dimension(:,:,:), allocatable :: f_lake_n  !! lakes fraction
      real(wp), dimension(:,:), allocatable :: f_ice     !! ice sheet fraction 
      real(wp), dimension(:,:), allocatable :: f_ice_grd     !! grounded ice sheet fraction 
      real(wp), dimension(:,:), allocatable :: f_ice_flt     !! floating ice sheet fraction 
      real(wp), dimension(:,:), allocatable :: z_sur     !! grid cell mean surface elevation [m]
      real(wp), dimension(:,:), allocatable :: z_sur_std     !! standard deviation of grid cell surface elevation [m]
      real(wp), dimension(:,:), allocatable :: z_sur_smooth_std     !! standard deviation of smoothed grid cell surface elevation [m]
      real(wp), dimension(:,:), allocatable :: z_ocn         !! grid cell mean elevation of ocean part [m]
      real(wp), dimension(:,:), allocatable :: z_ocn_min     !! grid cell min elevation of ocean part [m]
      real(wp), dimension(:,:), allocatable :: z_ocn_max     !! grid cell max elevation of ocean part [m]
      real(wp), dimension(:,:), allocatable :: z_ocn_max_q   !! grid cell quantile elevation of ocean part [m]
      real(wp), dimension(:,:), allocatable :: z_veg     !! grid cell mean elevation of ice-free land part [m]
      real(wp), dimension(:,:), allocatable :: z_veg_min     !! grid cell max elevation of ice-free land part [m]
      real(wp), dimension(:,:), allocatable :: z_veg_max     !! grid cell min elevation of ice-free land part [m]
      real(wp), dimension(:,:), allocatable :: z_veg_std     !! standard deviation of grid cell surface elevation, ice-free land part [m]
      real(wp), dimension(:,:), allocatable :: z_sur_lnd_std     !! standard deviation of grid cell surface elevation over land [m]
      real(wp), dimension(:,:), allocatable :: z_ice     !! grid gell mean elevation of ice sheet [m]
      real(wp), dimension(:,:), allocatable :: z_lake    !! grid gell mean elevation of lakes [m]
      real(wp), dimension(:,:), allocatable :: z_bed     !! grid gell mean bedrock elevation [m]
      integer, dimension(:,:), allocatable :: mask_coast     !! coast mask []
      integer, dimension(:,:), allocatable :: mask_coast2     !! coast mask []
      integer, dimension(:,:,:), allocatable :: i_coast_nbr     !! i index of neighbors of coastal cells []
      integer, dimension(:,:,:), allocatable :: j_coast_nbr     !! j index of neighbors of coastal cells []
      integer, dimension(:,:), allocatable :: coast_nbr     !! number of coastal cell neighbors []
      real(wp), dimension(:,:,:), allocatable :: coral_f_area
      real(wp), dimension(:,:,:), allocatable :: coral_f_topo
      real(wp), dimension(:,:), allocatable :: q_geo     !! geothermal_heat flux [W/m2]
      type(hires_type) :: hires
      integer :: n_lakes !! number of lakes
      type(lake_type), dimension(:), allocatable :: lake
      real(wp), dimension(:,:,:), allocatable :: f_drain_veg     !! drainage fractions for land 
      real(wp), dimension(:,:,:), allocatable :: f_drain_ice     !! drainage fractions for ice
      real(dp) :: ocn_area_tot
      real(dp) :: ocn_vol_tot
      real(wp) :: veg_area_tot
      real(wp) :: ice_area_tot
      real(wp) :: lake_area_tot
      integer, dimension(:,:), allocatable :: i_runoff, j_runoff
      integer, dimension(:,:), allocatable :: i_runoff_veg, j_runoff_veg 
      integer, dimension(:,:), allocatable :: i_runoff_ice, j_runoff_ice
      integer, dimension(:,:), allocatable :: drain_basins_ocn 
      integer, dimension(:), allocatable :: idivide_pac_atl 
      integer, dimension(:), allocatable :: idivide_atl_indpac
    end type geo_class

    private
    public :: geo_class, lake_type, n_lev

end module

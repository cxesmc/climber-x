!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ d e f
!
!  Purpose : definition of ocean model class 
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
module ocn_def

    use precision, only : wp
    use ocn_grid, only : grid_class

    implicit none
    
    ! define class with daily input fields
    type daily_input_save_class
      real(wp), dimension(:,:,:), allocatable :: stressxu
      real(wp), dimension(:,:,:), allocatable :: stressyv
      real(wp), dimension(:,:,:), allocatable :: stressxv
      real(wp), dimension(:,:,:), allocatable :: stressyu
      real(wp), dimension(:,:,:), allocatable :: wind
      real(wp), dimension(:,:,:), allocatable :: f_sic
      real(wp), dimension(:,:,:), allocatable :: slp
      real(wp), dimension(:,:,:), allocatable :: p_e_sic
      real(wp), dimension(:,:,:), allocatable :: fw_brines
      real(wp), dimension(:,:,:), allocatable :: runoff
      real(wp), dimension(:,:,:), allocatable :: runoff_veg
      real(wp), dimension(:,:,:), allocatable :: runoff_ice
      real(wp), dimension(:,:,:), allocatable :: runoff_lake
      real(wp), dimension(:,:,:), allocatable :: calving
      real(wp), dimension(:,:,:), allocatable :: bmelt_grd
      real(wp), dimension(:,:,:), allocatable :: bmelt_flt
      real(wp), dimension(:,:,:), allocatable :: bmelt
      real(wp), dimension(:,:,:), allocatable :: flx
    end type


    type ocn_class

     logical :: error                   !! error flag

     type(daily_input_save_class) :: daily_input_save

     type(grid_class) :: grid
     real(wp) :: saln0                    !! reference salinity [psu]
     real(wp) :: fw_glob                  !! annual net freshwater flux to ocean [kg/m2/s]
     real(wp) :: fw_glob_tmp                  !! annual global integral net freshwater flux to ocean [kg]
     real(wp) :: dvsf
     real(wp) :: cfc11_atm                  !! atmospheric cfc11 [ppt]
     real(wp) :: cfc12_atm                  !! atmospheric cfc12 [ppt]
     real(wp) :: hosing                 !! freshwater hosing flux [Sv]
     real(wp) :: noise_fw                 !! freshwater flux noise [kg/m2/s]
     real(wp) :: noise_flx                !! heat flux noise [W/m2]
     real(wp) :: amoc                   !! AMOC strength [Sv]
     real(wp) :: bering_tf              !! northward Bering Strait throughflow from parameterisation [Sv]
     real(wp) :: bering_fw              !! northwater freshwater transport through the Bering Strait from parameterisation [Sv]
     real(wp) :: A_bering              !! Bering Strait cross-sectional area [m2]     
     real(wp) :: buoy_sic_NA(6)   !! buoyancy flux from sea ice export over North Atlantic [N]
     logical, allocatable :: l_tracers_trans(:)    !! flags for bgc tracer transport
     logical, allocatable :: l_tracer_dic(:)    !! flags for bgc DIC tracers
     logical, allocatable :: l_tracers_isodiff(:)    !! flags for isopycnal diffusion of tracers 
     real(wp), allocatable :: z_ocn_max(:,:)    !! ocean floor elevation [m]
     real(wp), allocatable :: f_ocn(:,:)    !! ocean fraction in grid cell, including floating ice [1]
     real(wp), allocatable :: f_ocn_old(:,:)    !! old ocean fraction in grid cell [1]
     real(wp), allocatable :: f_ocn2(:,:)    !! ocean fraction in grid cell, excluding floating ice [1]
     integer, allocatable :: mask_coast(:,:)    !! mask of coastal cells []

     real(wp), allocatable :: f_sic(:,:)    !! sea ice fraction (of ocean fraction) [1]
     real(wp), allocatable :: wind(:,:)     !! 10 m wind speed [m/s]
     real(wp), allocatable :: slp(:,:)      !! sea level pressure [Pa]
     real(wp), allocatable :: stressxu(:,:), stressxv(:,:), stressyu(:,:), stressyv(:,:)
     real(wp), allocatable :: tau(:,:,:)    !! zonal and meridional wind stress [N/m2]
     real(wp), allocatable :: dtau_dz2(:,:,:)  !! second z-derivative of zonal wind stress [N/m4]
     real(wp), allocatable :: dtav_dz2(:,:,:)  !! second z-derivative of meridional wind stress [N/m4]
     real(wp), allocatable :: u(:,:,:,:)    !! 3D velocity field [m/s] 
     real(wp), allocatable :: ub(:,:,:)     !! barotropic velocity [m/s]
     real(wp), allocatable :: ub_isl(:,:,:,:)     !! barotropic velocity around islands[m/s]
     real(wp), allocatable :: psi(:,:)      !! barotropic stream function [kg/s]
  
     real(wp), allocatable :: flx(:,:)  !! heat flux into the ocean [W/m2]
     real(wp), allocatable :: fw(:,:)   !! freshwater flux into the ocean [kg/m2/s] 
     real(wp), allocatable :: fw_corr(:,:)   !! corrected (zero-mean annual) freshwater flux into the ocean [kg/m2/s] 
     real(wp), allocatable :: p_e_sic(:,:)   !! freshwater flux, precipitation-Evaporation+sea ice fluxes [kg/m2/s] 
     real(wp), allocatable :: fw_brines(:,:)   !! freshwater flux from brine rejection [kg/m2/s] 
     real(wp), allocatable :: runoff(:,:)   !! runoff into the ocean [kg/m2/s] 
     real(wp), allocatable :: runoff_veg(:,:)   !! runoff into the ocean [kg/m2/s] 
     real(wp), allocatable :: runoff_ice(:,:)   !! runoff into the ocean [kg/m2/s] 
     real(wp), allocatable :: melt_ice(:,:)   !! freshwater flux into the ocean from ice sheet melt [kg/m2/s] 
     real(wp), allocatable :: runoff_lake(:,:)   !! runoff into the ocean [kg/m2/s] 
     real(wp), allocatable :: calving(:,:)   !! ice calving into the ocean [kg/m2/s] 
     real(wp), allocatable :: bmelt(:,:)   !! basal melt of ice into the ocean [kg/m2/s] 
     real(wp), allocatable :: bmelt_grd(:,:)   !! basal melt of grounded ice sheet into the ocean [kg/m2/s] 
     real(wp), allocatable :: bmelt_flt(:,:)   !! basal melt of floating shelf ice into the ocean [kg/m2/s] 
     real(wp), allocatable :: fw_hosing(:,:)    !! freshwater hosing forcing [kg/m2/s]
     real(wp), allocatable :: fw_flux_adj(:,:)    !! freshwater flux adjustment [kg/m2/s]
     real(wp), allocatable :: flx_sur(:,:,:)  !! surface tracer input flux to the ocean [m/s * tracer concentration]
     real(wp), allocatable :: flx_bot(:,:,:)  !! bottom tracer input flux to the ocean [m/s * tracer concentration]
     real(wp), allocatable :: fw_noise(:,:)    !! freshwater noise [kg/m2/s]
     real(wp), allocatable :: flx_noise(:,:)    !! heat flux noise [W/m2]
     real(wp), allocatable :: ts(:,:,:,:)   !! 3D tracer fields [tracer concentration]
     real(wp), allocatable :: fdx(:,:,:,:)   !! zonal diffusive tracer volume flux [m3 * tracer conc]
     real(wp), allocatable :: fdy(:,:,:,:)   !! meridional diffusive tracer volume flux [m3 * tracer conc]
     real(wp), allocatable :: fdz(:,:,:,:)   !! vertical diffusive tracer volume flux [m3 * tracer conc]
     real(wp), allocatable :: fax(:,:,:,:)   !! zonal advective tracer volume flux [m3 * tracer conc]
     real(wp), allocatable :: fay(:,:,:,:)   !! meridional advective tracer volume flux [m3 * tracer conc]
     real(wp), allocatable :: faz(:,:,:,:)   !! vertical advective tracer volume flux [m3 * tracer conc]
     real(wp), allocatable :: dts_dt_adv(:,:,:,:)   !! tracer tendency due to advection [tracer conc/s]
     real(wp), allocatable :: dts_dt_diff(:,:,:,:)   !! tracer tendency due to diffusion [tracer conc/s]
     real(wp), allocatable :: sst_min(:,:)  !! minimum annual surface layer temperature for corals [degC]
     real(wp), allocatable :: sst_max(:,:)  !! maximum annual surface layer temperature for corals [degC]
     real(wp), allocatable :: rho(:,:,:)    !! density [kg/m3]
     real(wp), allocatable :: ke_tau(:,:) !! kinetic energy by wind input into the mixed layer []
     real(wp), allocatable :: mld(:,:)      !! mixed layer depth [m]
     real(wp), allocatable :: dconv(:,:)    !! depth of convection [m] 
     real(wp), allocatable :: dven(:,:)     !! depth of ventilation [m]
     real(wp), allocatable :: conv_pe(:,:)     !! potential energy released by convection [J/m2]
     integer, allocatable :: nconv(:,:)             !! number of layers mixed by convection []
     integer, allocatable :: kven(:,:)              !! number of layers mixed my ventialtion []
     real(wp), allocatable :: ssh(:,:)      !! elevation of the free surface [m]
     real(wp), allocatable :: q_geo(:,:)      !! geothermal heat flux [W/m2]
    end type

    private
    public :: ocn_class

end module

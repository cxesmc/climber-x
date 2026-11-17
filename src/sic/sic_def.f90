!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s i c _ d e f
!
!  Purpose : definition of sea ice model class 
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
module sic_def

    use precision, only : wp

    implicit none
    
    type sic_class   

      logical :: error                   !! error flag

      integer,  allocatable, dimension(:,:) :: mask_new_cell   !! mask of newly formed cells []
      real(wp), allocatable, dimension(:,:) :: f_ocn    !! ocean fraction in grid cell []
      real(wp), allocatable, dimension(:,:) :: u        !! zonal sea ice drift velocity on u-grid [m/s]
      real(wp), allocatable, dimension(:,:) :: v        !! meridional sea ice drift velocity on v-grid [m/s]
      real(wp), allocatable, dimension(:,:) :: uo      !! zonal top ocean layer velocity on u-grid [m/s]
      real(wp), allocatable, dimension(:,:) :: vo      !! meridional top ocean layer velocity on v-grid [m/s]
      real(wp), allocatable, dimension(:,:) :: tauxa    !! zonal wind stress over sea ice [N/m2]
      real(wp), allocatable, dimension(:,:) :: tauya    !! meridional wind stress over sea ice [N/m2]
      real(wp), allocatable, dimension(:,:) :: tauxo    !! zonal sea ice stress on ocean on u-grid [N/m2]
      real(wp), allocatable, dimension(:,:) :: tauyo    !! meridional sea ice stress on ocean on v-grid [N/m2]
      real(wp), allocatable, dimension(:,:) :: ssh    !! elevation of the ocean free surface [m]
      real(wp), allocatable, dimension(:,:) :: str_d  ! The divergence stress tensor component [Pa m].
      real(wp), allocatable, dimension(:,:) :: str_t  ! The tension stress tensor component [Pa m].
      real(wp), allocatable, dimension(:,:) :: str_s  ! The shearing stress tensor component [Pa m].
      real(wp), allocatable, dimension(:,:) :: fxic   ! Zonal force due to internal stresses [Pa].
      real(wp), allocatable, dimension(:,:) :: fxic_d ! Zonal force due to divergence internal stress [Pa].
      real(wp), allocatable, dimension(:,:) :: fxic_t ! Zonal force due to tension internal stress [Pa].
      real(wp), allocatable, dimension(:,:) :: fxic_s ! Zonal force due to shearing internal stress [Pa].
      real(wp), allocatable, dimension(:,:) :: Cor_u  ! Zonal Coriolis acceleration [m s-2].
      real(wp), allocatable, dimension(:,:) :: PFu    ! Zonal hydrostatic pressure driven acceleration [m s-2].
      real(wp), allocatable, dimension(:,:) :: fyic   ! Meridional force due to internal stresses [Pa].
      real(wp), allocatable, dimension(:,:) :: fyic_d ! Meridional force due to divergence internal stress [Pa].
      real(wp), allocatable, dimension(:,:) :: fyic_t ! Meridional force due to tension internal stress [Pa].
      real(wp), allocatable, dimension(:,:) :: fyic_s ! Meridional force due to shearing internal stress [Pa].
      real(wp), allocatable, dimension(:,:) :: Cor_v  ! Meridional Coriolis acceleration [m s-2].
      real(wp), allocatable, dimension(:,:) :: PFv    ! Meridional hydrostatic pressure driven acceleration [m s-2].
      real(wp), allocatable, dimension(:,:) :: sst  !! sea surface temperature [degree]
      real(wp), allocatable, dimension(:,:) :: sss  !! sea surface salinity [psu]
      real(wp), allocatable, dimension(:,:) :: wind     !! surface wind speed [m/s]
      real(wp), allocatable, dimension(:,:) :: pressure !! surface pressure [Pa]
      real(wp), allocatable, dimension(:,:) :: snow     !! snowfall rate [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: rain     !! rainfall rate [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: flx_lwd_sic  !! longwave downwelling radiation at the surface over sea ice [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_lwd_ocn  !! longwave downwelling radiation at the surface over ocean water [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_swnet_sic    !! net surface shortwave radiation absorbed over sea ice [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_swnet_ocn    !! net surface shortwave radiation absorbed over ocean water [W/m2]
      real(wp), allocatable, dimension(:,:) :: t_air_sic    !! surface air temperature over sea ice [K]
      real(wp), allocatable, dimension(:,:) :: t_air_ocn    !! surface air temperature over ocean water [K]
      real(wp), allocatable, dimension(:,:) :: q_air_sic    !! surface air specific humidity over sea ice [kg/kg]
      real(wp), allocatable, dimension(:,:) :: q_air_ocn    !! surface air specific humidity over ocean water [kg/kg]
      real(wp), allocatable, dimension(:,:) :: coszm    !! cosine of the radiation weighted daily solar zenith angle []
      real(wp), allocatable, dimension(:,:) :: h_sic_mean   !! grid cell mean sea ice thickness [m]
      real(wp), allocatable, dimension(:,:) :: h_sic   !! sea ice thickness over ice covered area [m]
      real(wp), allocatable, dimension(:,:) :: f_sic   !! sea ice fraction of available ocean fraction (f_ocn) []
      real(wp), allocatable, dimension(:,:) :: dh_sic_dt_dyn   !! grid cell mean sea ice thickness change from dynamics (transport) [m/s]
      real(wp), allocatable, dimension(:,:) :: r_open_shear   !! opening rate by shear deformation 0.5*(del-|div|) from the EVP dynamics [1/s]
      real(wp), allocatable, dimension(:,:) :: df_sic_shear   !! sea ice fraction change from shear opening (ridging) [1/s]
      real(wp), allocatable, dimension(:,:) :: dh_sic_dt_therm   !! grid cell mean sea ice thickness change from thermodynamics [m/s]
      real(wp), allocatable, dimension(:,:) :: dh_sic_sic  !! change in sea ice thickness from thermodynamics over ice-covered areas [m]
      real(wp), allocatable, dimension(:,:) :: dh_sic_ocn  !! change in sea ice thickness from thermodynamics over ice-free areas [m]
      real(wp), allocatable, dimension(:,:) :: fax_sic   !! zonal sea ice volume flux from advection [m3]
      real(wp), allocatable, dimension(:,:) :: fay_sic   !! meridional sea ice volume flux from advection [m3]
      real(wp), allocatable, dimension(:,:) :: fdx_sic   !! zonal sea ice volume flux from diffusion [m3]
      real(wp), allocatable, dimension(:,:) :: fdy_sic   !! meridional sea ice volume flux from diffusion [m3]
      real(wp), allocatable, dimension(:,:) :: fax_snow   !! zonal snow volume flux from advection [m3]
      real(wp), allocatable, dimension(:,:) :: fay_snow   !! meridional snow volume flux from advection [m3]
      real(wp), allocatable, dimension(:,:) :: fdx_snow   !! zonal snow volume flux from diffusion [m3]
      real(wp), allocatable, dimension(:,:) :: fdy_snow   !! meridional snow volume flux from diffusion [m3]
      real(wp), allocatable, dimension(:,:) :: h_snow_mean   !! grid cell mean snow thickness [m]
      real(wp), allocatable, dimension(:,:) :: h_snow   !! snow thickness over ice covered area [m]
      real(wp), allocatable, dimension(:,:) :: w_snow   !! snow water equivalent over ice covered area [m]
      real(wp), allocatable, dimension(:,:) :: w_snow_max   !! maximum seasonal snow water equivalent over ice covered area [m]
      real(wp), allocatable, dimension(:,:) :: dh_snow  !! change in snow thickness [m]
      real(wp), allocatable, dimension(:,:) :: dh_snow_dt_dyn   !! grid cell mean snow thickness change from dynamics (transport) [m/s]
      real(wp), allocatable, dimension(:,:) :: dh_snow_dt_therm   !! grid cell mean snow thickness change from thermodynamics [m/s]     
      real(wp), allocatable, dimension(:,:) :: flx_melt_top !! snow/ice melt flux from the top [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_melt_bot !! snow/ice melt flux from the bottom [W/m2]
      real(wp), allocatable, dimension(:,:) :: t_skin_sic   !! sea ice skin temperature [K]
      real(wp), allocatable, dimension(:,:) :: t_skin_ocn   !! ocean skin temperature [K]
      real(wp), allocatable, dimension(:,:) :: flx_sh_sic   !! sensible heat flux over sea ice (positive upwards) [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_lh_sic   !! latent heat flux over sea ice (positive upwards) [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_lwu_sic  !! upwelling longwave radiation flux over sea ice [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_ocn_sic  !! heat flux into the ocean from sea ice [W/m2]
      real(wp), allocatable, dimension(:,:) :: fw_ocn_sic   !! freshwater flux into the ocean from sea ice [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: evp_sic      !! evaporation from sea ice [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: flx_sh_ocn   !! sensible heat flux over ice-free ocean (positive upwards) [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_lh_ocn   !! latent heat flux over ice-free ocean (positive upwards) [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_lwu_ocn  !! upwelling longwave radiation flux over ice-free ocean [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_ocn_ocn  !! heat flux into the ice-free ocean [W/m2]
      real(wp), allocatable, dimension(:,:) :: fw_ocn_ocn   !! freshwater flux into the ice-free ocean [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: evp_ocn      !! evaporation from ice-free ocean [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: flx_sh       !! sensible heat flux (positive upwards) [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_lh       !! latent heat flux (positive upwards) [W/m2]
      real(wp), allocatable, dimension(:,:) :: flx_lwu      !! upwelling longwave radiation flux[W/m2]
      real(wp), allocatable, dimension(:,:) :: evp          !! evaporation [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: fw_ocn       !! grid averaged freshwater flux to the ocean (excluding runoff) [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: fw_brines    !! grid averaged freshwater flux to the ocean from brine rejection [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: flx_ocn      !! grid averaged heat flux into the ocean [W/m2]
      real(wp), allocatable, dimension(:,:) :: snow_grain   !! snow grain size [micro m]
      real(wp), allocatable, dimension(:,:) :: dust_dep     !! dust deposition [kg/m2/s]
      real(wp), allocatable, dimension(:,:) :: dust_con     !! dust concentration in snow [kg/kg]
      real(wp), allocatable, dimension(:,:) :: alb_snow_vis_dir  !! visible, direct beam snow albedo []
      real(wp), allocatable, dimension(:,:) :: alb_snow_vis_dif  !! visible, diffuse snow albedo []
      real(wp), allocatable, dimension(:,:) :: alb_snow_nir_dir  !! near infrared, direct beam snow albedo []
      real(wp), allocatable, dimension(:,:) :: alb_snow_nir_dif  !! near infrared, diffuse snow albedo []
      real(wp), allocatable, dimension(:,:) :: alb_ocn_vis_dir  !! visible, direct beam ocean albedo []
      real(wp), allocatable, dimension(:,:) :: alb_ocn_vis_dif  !! visible, diffuse ocean albedo []
      real(wp), allocatable, dimension(:,:) :: alb_ocn_nir_dir  !! near infrared, direct beam ocean albedo []
      real(wp), allocatable, dimension(:,:) :: alb_ocn_nir_dif  !! near infrared, diffuse ocean albedo []
      real(wp), allocatable, dimension(:,:) :: alb_sic_vis_dir  !! visible, direct beam sea ice albedo []
      real(wp), allocatable, dimension(:,:) :: alb_sic_vis_dif  !! visible, diffuse sea ice albedo []
      real(wp), allocatable, dimension(:,:) :: alb_sic_nir_dir  !! near infrared, direct beam sea ice albedo []
      real(wp), allocatable, dimension(:,:) :: alb_sic_nir_dif  !! near infrared, diffuse sea ice albedo []
      real(wp), allocatable, dimension(:,:) :: albedo_sic   !! broadband sea ice albedo, assuming half cloud cover []
      real(wp), allocatable, dimension(:,:) :: albedo_ocn   !! broadband ocean albedo, assuming half cloud cover []
      real(wp), allocatable, dimension(:,:) :: rough_m_ocn  !! surface roughness for momentum over ice-free ocean []
      real(wp), allocatable, dimension(:,:) :: rough_h_ocn  !! surface roughness for heat/moisture over ice-free ocean []
      real(wp), allocatable, dimension(:,:) :: rough_m_sic  !! surface roughness for momentum over sea ice []
      real(wp), allocatable, dimension(:,:) :: rough_h_sic  !! surface roughness for heat/moisture over sea ice []
      real(wp), allocatable, dimension(:,:) :: Cde_ocn  !! drag coefficient for evaporation over ice-free ocean []
      real(wp), allocatable, dimension(:,:) :: Cdh_ocn  !! drag coefficient for sensible heat flux over ice-free ocean []
      real(wp), allocatable, dimension(:,:) :: Cde_sic  !! drag coefficient for evaporation over sea ice []
      real(wp), allocatable, dimension(:,:) :: Cdh_sic  !! drag coefficient for sensible heat flux over sea ice []

      real(wp) :: buoy_sic_NA(6)   !! buoyancy flux from sea ice export over North Atlantic [N]

    end type

    private
    public :: sic_class

end module

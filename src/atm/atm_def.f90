!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : a t m _ d e f
!
!  Purpose : definition of atmosphere model class 
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
module atm_def

    use atm_params, only : wp
    use atm_grid, only : nm
    use atm_grid, only : i_ocn, i_sic, i_lnd, i_ice, i_lake

    implicit none
    
    type atm_class

      logical :: error

      ! grid information
      integer :: nm = nm
      integer :: i_ocn = i_ocn
      integer :: i_sic = i_sic
      integer :: i_lnd = i_lnd
      integer :: i_ice = i_ice
      integer :: i_lake = i_lake

      real(wp), allocatable :: pl(:)
      real(wp), allocatable :: zl(:)

      real(wp) :: co2      !! atmospheric CO2 (ppmv)
      real(wp) :: co2e     !! equivalent atmospheric CO2 (ppmv)
      real(wp) :: ch4      !! atmospheric CH4 (ppb)
      real(wp) :: n2o      !! atmospheric N2O (ppb)
      real(wp) :: cfc11    !! atmospheric CFC11 (ppt)
      real(wp) :: cfc12    !! atmospheric CFC12 (ppt)

      real(wp) :: had_fi   !! ITCZ position latitude (radians)
      real(wp) :: had_width   !! width of Hadley cell (radians)

      real(wp) :: eccentricity
      real(wp) :: precession
      real(wp) :: obliquity

      real(wp) :: t2m_glob_ann
      real(wp) :: dt2m_glob_ann_cum

      real(wp), allocatable, dimension(:,:,:) :: solar
      real(wp), allocatable, dimension(:,:) :: solarm
      real(wp), allocatable, dimension(:,:,:) :: cosz
      real(wp), allocatable, dimension(:,:) :: coszm
      
      real(wp), allocatable, dimension(:,:) :: cam 
      real(wp), allocatable, dimension(:,:) :: co2flx
      real(wp), allocatable, dimension(:,:) :: C13flx
      real(wp), allocatable, dimension(:,:) :: C14flx

      integer, allocatable, dimension(:) :: idivide_pac_atl
      integer, allocatable, dimension(:) :: idivide_atl_indpac

      real(wp), allocatable, dimension(:,:,:) :: zs
      real(wp), allocatable, dimension(:,:) :: zsa
      real(wp), allocatable, dimension(:,:) :: zsa_smooth
      real(wp), allocatable, dimension(:,:) :: slope
      real(wp), allocatable, dimension(:,:) :: slope_x
      real(wp), allocatable, dimension(:,:) :: slope_y
      real(wp), allocatable, dimension(:,:) :: pzsa0
      real(wp), allocatable, dimension(:,:) :: pzsa
      real(wp), allocatable, dimension(:,:) :: frlnd
      real(wp), allocatable, dimension(:,:) :: frocn
      real(wp), allocatable, dimension(:,:) :: f_ice_lake
      real(wp), allocatable, dimension(:,:) :: sigoro

      real(wp), allocatable, dimension(:,:) :: tam        !! extrapolated surface temperature (K)
      real(wp), allocatable, dimension(:,:) :: qam        !! extrapolated surface specific humidity (kg/kg)
      real(wp), allocatable, dimension(:,:) :: ram        !! extrapolated surface relative humidity (/)
      real(wp), allocatable, dimension(:,:) :: gams       !! lapse rate in the boundary layer (K/m)
      real(wp), allocatable, dimension(:,:) :: gamb       !! lapse rate at the lower troposphere (K/m)
      real(wp), allocatable, dimension(:,:) :: gamt       !! lapse rate in the upper troposphere (K/m)
      real(wp), allocatable, dimension(:,:) :: dam        !! surface dust mass mixing ratio (kg/kg)
      real(wp), allocatable, dimension(:,:) :: hrm        !! vertical scale for relative humidity(m)
      real(wp), allocatable, dimension(:,:) :: hqeff      !! effective vertical scale for specific humidity (m)
      real(wp), allocatable, dimension(:,:) :: wcon       !! atmospheric water content (kg m-2)
      real(wp), allocatable, dimension(:,:) :: cld_rh        !! cloud fraction (.)
      real(wp), allocatable, dimension(:,:) :: cld_low        !! cloud fraction (.)
      real(wp), allocatable, dimension(:,:) :: cld        !! cloud fraction (.)
      real(wp), allocatable, dimension(:,:) :: cld_dat        !! cloud fraction (.)
      real(wp), allocatable, dimension(:,:) :: cld_day_dat        !! daytime cloud fraction (.)
      real(wp), allocatable, dimension(:,:) :: prc        !! total precipitation (kg m-2 s-1)
      real(wp), allocatable, dimension(:,:,:) :: prcw       !! rain(kg m-2 s-1)
      real(wp), allocatable, dimension(:,:,:) :: prcs       !! snowfall(kg m-2 s-1)
      real(wp), allocatable, dimension(:,:) :: prc_conv   !! precipitation from supersaturation (kg m-2 s-1)
      real(wp), allocatable, dimension(:,:) :: prc_wcon  
      real(wp), allocatable, dimension(:,:) :: prc_over
      real(wp), allocatable, dimension(:,:) :: hcld       !! cloud height (m)
      real(wp), allocatable, dimension(:,:) :: clot       !! cloud optical thickness (.)
      real(wp), allocatable, dimension(:,:) :: alb_cld    !! cloud albedo (.)
      real(wp), allocatable, dimension(:,:) :: htrop      !! tropopause height (m)
      real(wp), allocatable, dimension(:) :: ptrop        !! pressure level of tropopause, zonal mean (.)
      real(wp), allocatable, dimension(:,:) :: ttrop      !! tropopause temperature (K)
      real(wp), allocatable, dimension(:,:) :: winda      !! average surface wind magnitude (m s-1)
      real(wp), allocatable, dimension(:,:,:) :: wind     !! surface wind magnitude (m s-1)
      real(wp), allocatable, dimension(:,:) :: aerosol_ot !! aerosol optical thickness (.)
      real(wp), allocatable, dimension(:,:) :: aerosol_im !! aerosol imaginary refractive index (.)
      real(wp), allocatable, dimension(:,:) :: so4        !! atmospheric SO4 load (kg/m2)
      real(wp), allocatable, dimension(:,:,:) :: o3       !! atmospheric O3 concentration (mol/mol)

      real(wp), allocatable, dimension(:,:) :: hdust      !! dust height scale (m)
      real(wp), allocatable, dimension(:,:) :: dust_load  !! dust load (kg/m2)
      real(wp), allocatable, dimension(:,:) :: dust_emis  !! dust emissions (kg/m2/s) 
      real(wp), allocatable, dimension(:,:) :: dust_dep   !! dust deposition (kg/m2/s)
      real(wp), allocatable, dimension(:,:) :: dust_dep_dry !! dust dry deposition (kg/m2/s)
      real(wp), allocatable, dimension(:,:) :: dust_dep_wet !! dust wet deposition (kg/m2/s)
      real(wp), allocatable, dimension(:,:) :: dust_ot    !! dust optical thickness (.)

      real(wp), allocatable, dimension(:,:,:) :: frst
      real(wp), allocatable, dimension(:,:,:) :: tskin     !! skin temperature (K)
      real(wp), allocatable, dimension(:,:,:) :: t2        !! 2m surface air temperature (K)
      real(wp), allocatable, dimension(:,:,:) :: ra2
      real(wp), allocatable, dimension(:,:,:) :: q2        !! 2m surface specific humidity (kg/kg)
      real(wp), allocatable, dimension(:,:,:) :: r2        !! 2m surface relative humidity
      real(wp), allocatable, dimension(:,:,:) :: alb_vu_s  !! visible+UV clear sky surface albedo (.)
      real(wp), allocatable, dimension(:,:,:) :: alb_vu_c  !! visible+UV cloudy sky surface albedo (.)
      real(wp), allocatable, dimension(:,:,:) :: alb_ir_s  !! infrared clear sky surface albedo (.)
      real(wp), allocatable, dimension(:,:,:) :: alb_ir_c  !! infrared cloudy sky surface albedo (.)
      real(wp), allocatable, dimension(:,:,:) :: cd        !! drag coefficient
      real(wp), allocatable, dimension(:,:,:) :: cd0       !! drag coefficient without orographic component
      real(wp), allocatable, dimension(:,:,:) :: z0m       !! surface roughness for momentum
      real(wp), allocatable, dimension(:,:) :: zoro        !! orographic roughness

      real(wp), allocatable, dimension(:,:) :: ra2a
      real(wp), allocatable, dimension(:,:) :: cda
      real(wp), allocatable, dimension(:,:) :: cd0a       !! drag coefficient without orographic component
      real(wp), allocatable, dimension(:,:) :: sha
      real(wp), allocatable, dimension(:,:) :: lha
      real(wp), allocatable, dimension(:,:) :: evpa
      real(wp), allocatable, dimension(:,:) :: tskina
      real(wp), allocatable, dimension(:,:) :: t2a
      real(wp), allocatable, dimension(:,:) :: q2a
      real(wp), allocatable, dimension(:,:) :: r2a
      real(wp), allocatable, dimension(:,:) :: rskina

      ! 3D longitude-latitude-height
      real(wp), allocatable, dimension(:,:,:) :: t3   !! atmospheric temperature (K)
      real(wp), allocatable, dimension(:,:,:) :: q3   !! atmospheric specific humidity (kg/kg)
      real(wp), allocatable, dimension(:,:,:) :: tp   !! potential temperature
      real(wp), allocatable, dimension(:,:,:) :: d3   !! atmospheric dust mass mixing ratio (kg/kg)

      real(wp), allocatable, dimension(:,:) :: acbar
      real(wp), allocatable, dimension(:,:) :: sin_cos_acbar
      real(wp), allocatable, dimension(:,:,:) :: cos_acbar
      real(wp), allocatable, dimension(:,:,:) :: sin_acbar
      real(wp), allocatable, dimension(:,:,:) :: epsa
      real(wp), allocatable, dimension(:,:) :: atsl
      real(wp), allocatable, dimension(:,:) :: aslp
      real(wp), allocatable, dimension(:,:) :: aslp_temp
      real(wp), allocatable, dimension(:,:) :: aslp_topo
      real(wp), allocatable, dimension(:,:) :: dz500
      real(wp), allocatable, dimension(:,:) :: slp
      real(wp), allocatable, dimension(:,:) :: slp_dat
      real(wp), allocatable, dimension(:,:) :: tsl_dat
      real(wp), allocatable, dimension(:,:,:) :: ps
      real(wp), allocatable, dimension(:,:) :: psa
      real(wp), allocatable, dimension(:,:,:) :: us
      real(wp), allocatable, dimension(:,:,:) :: vs
      real(wp), allocatable, dimension(:,:) :: usk
      real(wp), allocatable, dimension(:,:) :: vsk
      real(wp), allocatable, dimension(:,:) :: ugb
      real(wp), allocatable, dimension(:,:) :: vgb
      real(wp), allocatable, dimension(:,:) :: ugbf
      real(wp), allocatable, dimension(:,:) :: vgbf
      real(wp), allocatable, dimension(:,:) :: uab
      real(wp), allocatable, dimension(:,:) :: vab
      real(wp), allocatable, dimension(:,:,:) :: taux
      real(wp), allocatable, dimension(:,:,:) :: tauy
      real(wp), allocatable, dimension(:) :: uz500
      real(wp), allocatable, dimension(:,:) :: wcld
      real(wp), allocatable, dimension(:,:) :: woro
      real(wp), allocatable, dimension(:,:) :: wsyn       !! vertical synoptic wind component at ~700 hPa (m/s)
      real(wp), allocatable, dimension(:,:) :: weff 
      real(wp), allocatable, dimension(:,:) :: fweff 

      real(wp), allocatable, dimension(:,:,:) :: ua
      real(wp), allocatable, dimension(:,:,:) :: va
      real(wp), allocatable, dimension(:,:,:) :: u3
      real(wp), allocatable, dimension(:,:,:) :: v3
      real(wp), allocatable, dimension(:,:,:) :: uter
      real(wp), allocatable, dimension(:,:,:) :: vter
      real(wp), allocatable, dimension(:,:,:) :: uterf
      real(wp), allocatable, dimension(:,:,:) :: vterf
      real(wp), allocatable, dimension(:,:,:) :: fax
      real(wp), allocatable, dimension(:,:,:) :: faxo
      real(wp), allocatable, dimension(:,:,:) :: fay
      real(wp), allocatable, dimension(:,:,:) :: fayo
      real(wp), allocatable, dimension(:,:) :: fac
      real(wp), allocatable, dimension(:,:,:) :: w3

      real(wp), allocatable, dimension(:,:) :: convdse
      real(wp), allocatable, dimension(:,:) :: convwtr
      real(wp), allocatable, dimension(:,:) :: convdst
      real(wp), allocatable, dimension(:,:) :: convco2
      real(wp), allocatable, dimension(:,:) :: faxdse
      real(wp), allocatable, dimension(:,:) :: faxwtr
      real(wp), allocatable, dimension(:,:) :: faxdst
      real(wp), allocatable, dimension(:,:) :: faxco2
      real(wp), allocatable, dimension(:,:) :: faydse
      real(wp), allocatable, dimension(:,:) :: faywtr
      real(wp), allocatable, dimension(:,:) :: faydst
      real(wp), allocatable, dimension(:,:) :: fayco2
      real(wp), allocatable, dimension(:,:) :: fdxdse
      real(wp), allocatable, dimension(:,:) :: fdxwtr
      real(wp), allocatable, dimension(:,:) :: fdxdst
      real(wp), allocatable, dimension(:,:) :: fdxco2
      real(wp), allocatable, dimension(:,:) :: fdydse
      real(wp), allocatable, dimension(:,:) :: fdywtr
      real(wp), allocatable, dimension(:,:) :: fdydst
      real(wp), allocatable, dimension(:,:) :: fdyco2

      real(wp), allocatable, dimension(:,:,:) :: fswr_sur
      real(wp), allocatable, dimension(:,:,:) :: fswr_sur_cs
      real(wp), allocatable, dimension(:,:,:) :: fswr_sur_cld
      real(wp), allocatable, dimension(:,:,:) :: flwr_dw_sur
      real(wp), allocatable, dimension(:,:,:) :: flwr_dw_sur_cs
      real(wp), allocatable, dimension(:,:,:) :: flwr_dw_sur_cld
      real(wp), allocatable, dimension(:,:,:) :: flwr_up_sur

      real(wp), allocatable, dimension(:,:) :: dswd_dalb_vu_cs
      real(wp), allocatable, dimension(:,:) :: dswd_dalb_ir_cs
      real(wp), allocatable, dimension(:,:) :: dswd_dalb_vu_cld
      real(wp), allocatable, dimension(:,:) :: dswd_dalb_ir_cld
      real(wp), allocatable, dimension(:,:) :: dswd_dz_ir_cs
      real(wp), allocatable, dimension(:,:) :: dswd_dz_ir_cld
      real(wp), allocatable, dimension(:,:) :: swr_dw_sur_vis_cs
      real(wp), allocatable, dimension(:,:) :: swr_dw_sur_nir_cs
      real(wp), allocatable, dimension(:,:) :: swr_dw_sur_vis_cld
      real(wp), allocatable, dimension(:,:) :: swr_dw_sur_nir_cld

      real(wp), allocatable, dimension(:,:) :: rb_top
      real(wp), allocatable, dimension(:,:) :: rb_sur
      real(wp), allocatable, dimension(:,:) :: rb_atm
      real(wp), allocatable, dimension(:,:) :: rb_str     
      real(wp), allocatable, dimension(:,:) :: swr_dw_top
      real(wp), allocatable, dimension(:,:) :: swr_top
      real(wp), allocatable, dimension(:,:) :: swr_top_cs
      real(wp), allocatable, dimension(:,:) :: swr_top_cld
      real(wp), allocatable, dimension(:,:) :: swr_sur
      real(wp), allocatable, dimension(:,:) :: lwr_top
      real(wp), allocatable, dimension(:,:) :: lwr_top_cs
      real(wp), allocatable, dimension(:,:) :: lwr_top_cld
      real(wp), allocatable, dimension(:,:) :: lwr_sur
      real(wp), allocatable, dimension(:,:) :: lwr_tro
      real(wp), allocatable, dimension(:,:) :: lwr_cld

      real(wp), allocatable, dimension(:,:) :: tsl
      real(wp), allocatable, dimension(:,:) :: tsksl

      real(wp), allocatable, dimension(:,:) :: eke
      real(wp), allocatable, dimension(:,:) :: sam
      real(wp), allocatable, dimension(:,:) :: sam2
      real(wp), allocatable, dimension(:,:) :: synprod
      real(wp), allocatable, dimension(:,:) :: syndiss
      real(wp), allocatable, dimension(:,:) :: synadv
      real(wp), allocatable, dimension(:,:) :: syndif
      real(wp), allocatable, dimension(:,:,:) :: synsur
      real(wp), allocatable, dimension(:,:) :: cdif     
      real(wp), allocatable, dimension(:,:) :: diffxdse
      real(wp), allocatable, dimension(:,:) :: diffydse
      real(wp), allocatable, dimension(:,:) :: diffxwtr
      real(wp), allocatable, dimension(:,:) :: diffywtr
      real(wp), allocatable, dimension(:,:) :: diffxdst
      real(wp), allocatable, dimension(:,:) :: diffydst

    end type

    private
    public :: atm_class

end module

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : a t m _ o u t
!
!  Purpose : diagnostics and output of atmospheric model
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
module atm_out

  use atm_params, only : wp
  use constants, only : pi, T0, g, Rd, fqsat, frac_vu, z_sfl
  use dim_name, only: dim_month, dim_time, dim_x, dim_y
  use timer, only : n_accel, year, year_clim, year_now, sec_day, mon, doy, nday_year, &
  nmon_year, nstep_year_atm, nstep_mon_atm
  use timer, only : nyout_atm, ny_out_ts, y_out_ts_clim, time_out_ts_clim
  use timer, only : time_soy_atm, time_eoy_atm, time_eom_atm, time_out_atm
  use control, only : out_dir
  use climber_grid, only : lon, lat, lon0, lat0, dlon, dlat
  use atm_grid, only : im, jm, km, imc, jmc, kmc, nm, aim, k850, k500, dpl, sqr, esqr, zl, zc, pl, &
  dxt, dy, cost, sint, i_ice
  use atm_params, only : l_daily_output
  use atm_params, only : amas, p0, ra, gad, cle, cls, cp, cv, hatm, sigma_so4, c_syn_6, l_co2d
  use atm_params, only : l_output_flx3d, l_output_extended
  use vesta_mod, only : t_prof, rh_prof
  use atm_def, only : atm_class
  use ncio

  implicit none

  private
  public :: atm_diag, atm_diag_init

  integer :: nout

  integer :: i_grl, j_grl
  integer :: i_grl_1, i_grl_2, j_grl_1, j_grl_2
  integer :: i_natl_1, i_natl_2, j_natl_1, j_natl_2
  integer :: i_ant, j_ant
  integer :: i_wais, j_wais

  type ts_out
    real(wp) :: eccentricity
    real(wp) :: precession
    real(wp) :: obliquity
    real(wp) :: Smax65N
    real(wp) :: Smax65S
    real(wp) :: Smax30N
    real(wp) :: Smax30S
    real(wp) :: co2
    real(wp) :: co2e
    real(wp) :: co2flx
    real(wp) :: co2d_avg
    real(wp) :: co2d_wais
    real(wp) :: ch4
    real(wp) :: n2o
    real(wp) :: so4
    real(wp) :: t2m    
    real(wp) :: tnh
    real(wp) :: tsh    
    real(wp) :: dtnsh    
    real(wp) :: tn6090
    real(wp) :: ts6090
    real(wp) :: tgrl
    real(wp) :: tgrl1
    real(wp) :: tgrl1_jja
    real(wp) :: tnatl
    real(wp) :: tant    
    real(wp) :: dust_grl
    real(wp) :: dust_wais
    real(wp) :: q2m    
    real(wp) :: r2m    
    real(wp) :: tam
    real(wp) :: tskin
    real(wp) :: prc       
    real(wp) :: evp       
    real(wp) :: wcon      
    real(wp) :: cld       
    real(wp) :: sha 
    real(wp) :: lha
    real(wp) :: rbtop     
    real(wp) :: rbatm     
    real(wp) :: rbsur     
    real(wp) :: snettop   
    real(wp) :: snetsur   
    real(wp) :: lnettop   
    real(wp) :: lnetsur 
    real(wp) :: ldwnsur   
    real(wp) :: rbtopcs     
    real(wp) :: rbatmcs     
    real(wp) :: rbsurcs     
    real(wp) :: snettopcs   
    real(wp) :: snetsurcs   
    real(wp) :: lnettopcs   
    real(wp) :: lnetsurcs 
    real(wp) :: ldwnsurcs   
    real(wp) :: lupsur   
    real(wp) :: sol       
    real(wp) :: plan
    real(wp) :: sreftop   
    real(wp) :: slp       
    real(wp) :: swrcrf    
    real(wp) :: lwrcrf    
    real(wp) :: htrop     
    real(wp) :: ttrop     
    real(wp) :: hcld      
    real(wp) :: dustload
    real(wp) :: dustdep
    real(wp) :: dustdrydep
    real(wp) :: dustwetdep
    real(wp) :: fw_pac_atl 
    real(wp) :: fw_atl_indpac
    real(wp) :: t2m_l
    real(wp) :: q2m_l    
    real(wp) :: r2m_l    
    real(wp) :: tam_l
    real(wp) :: tskin_l
    real(wp) :: prc_l       
    real(wp) :: evp_l       
    real(wp) :: wcon_l      
    real(wp) :: cld_l       
    real(wp) :: sha_l 
    real(wp) :: lha_l
    real(wp) :: rbtop_l     
    real(wp) :: rbatm_l     
    real(wp) :: rbsur_l     
    real(wp) :: snettop_l   
    real(wp) :: snetsur_l   
    real(wp) :: lnettop_l   
    real(wp) :: lnetsur_l 
    real(wp) :: itcz
    real(wp) :: itcz_min
    real(wp) :: itcz_max
    real(wp) :: had_width
  end type

  type a_out
    real(wp), allocatable, dimension(:,:,:) :: solar
    real(wp), allocatable, dimension(:,:,:) :: cosz
    real(wp), allocatable, dimension(:,:) :: solarm
    real(wp), allocatable, dimension(:,:) :: coszm

    real(wp) :: had_fi
    real(wp) :: had_width

    real(wp), allocatable, dimension(:) :: ptrop
    real(wp), allocatable, dimension(:) :: tslz
    real(wp), allocatable, dimension(:) :: tskslz
    real(wp), allocatable, dimension(:) :: ekez
    real(wp), allocatable, dimension(:) :: slpz
    real(wp), allocatable, dimension(:) :: vabz
    real(wp), allocatable, dimension(:) :: acbarz
    real(wp), allocatable, dimension(:) :: uz850
    real(wp), allocatable, dimension(:) :: uz500

    real(wp), allocatable, dimension(:,:,:) :: zs
    real(wp), allocatable, dimension(:,:) :: zsa
    real(wp), allocatable, dimension(:,:) :: zsa_smooth
    real(wp), allocatable, dimension(:,:) :: slope
    real(wp), allocatable, dimension(:,:) :: slope_x
    real(wp), allocatable, dimension(:,:) :: slope_y
    real(wp), allocatable, dimension(:,:) :: frlnd
    real(wp), allocatable, dimension(:,:) :: frocn

    real(wp), allocatable, dimension(:,:) :: tam        !! extrapolated surface temperature (K)
    real(wp), allocatable, dimension(:,:) :: dtamdt     !! rate of change of tam (K/day)
    real(wp), allocatable, dimension(:,:) :: qam        !! extrapolated surface specific humidity (kg/kg)
    real(wp), allocatable, dimension(:,:) :: ram        !! extrapolated surface relative humidity (/)
    real(wp), allocatable, dimension(:,:) :: gams       !! lapse rate in the boundary layer
    real(wp), allocatable, dimension(:,:) :: gamb       !! lapse rate at the lower troposphere
    real(wp), allocatable, dimension(:,:) :: gamt       !! lapse rate in the upper troposphere
    real(wp), allocatable, dimension(:,:) :: dam        !! surface dust mass mixing ratio (kg/kg)
    real(wp), allocatable, dimension(:,:) :: hqeff      !! vertical moisture scale (m)
    real(wp), allocatable, dimension(:,:) :: hrm        !! vertical scale for relative humidity(m)
    real(wp), allocatable, dimension(:,:) :: wcon       !! atmospheric water content (kg m-2)
    real(wp), allocatable, dimension(:,:) :: cld_rh        !! cloud fraction (.)
    real(wp), allocatable, dimension(:,:) :: cld_low        !! cloud fraction (.)
    real(wp), allocatable, dimension(:,:) :: cld        !! cloud fraction (.)
    real(wp), allocatable, dimension(:,:) :: prc        !! total precipitation (kg m-2 s-1)
    real(wp), allocatable, dimension(:,:) :: prcw       !! rain(kg m-2 s-1)
    real(wp), allocatable, dimension(:,:) :: prcs       !! snowfall(kg m-2 s-1)
    real(wp), allocatable, dimension(:,:) :: prc_conv 
    real(wp), allocatable, dimension(:,:) :: prc_wcon  
    real(wp), allocatable, dimension(:,:) :: prc_over
    real(wp), allocatable, dimension(:,:) :: hcld       !! cloud height (m)
    real(wp), allocatable, dimension(:,:) :: ctt        !! top cloud temperature (K)
    real(wp), allocatable, dimension(:,:) :: clot       !! cloud optical thickness (.)
    real(wp), allocatable, dimension(:,:) :: alb_cld    !! cloud albedo (.)
    real(wp), allocatable, dimension(:,:) :: alb_sur_cs    !! surface albedo for clear sky(.)
    real(wp), allocatable, dimension(:,:) :: alb_sur_cld   !! surface albedo for cloudy sky (.)
    real(wp), allocatable, dimension(:,:) :: htrop      !! tropopause height (m)
    real(wp), allocatable, dimension(:,:) :: ttrop      !! tropopause temperature (K)
    real(wp), allocatable, dimension(:,:) :: wind       !! average surface wind magnitude (m s-1)

    real(wp), allocatable, dimension(:,:,:) :: frst
    real(wp), allocatable, dimension(:,:,:) :: tskin     !! skin temperature (K)
    real(wp), allocatable, dimension(:,:,:) :: t2        !! 2m surface air temperature (K)
    real(wp), allocatable, dimension(:,:,:) :: q2        !! 2m surface specific humidity (kg/kg)
    real(wp), allocatable, dimension(:,:,:) :: r2        !! 2m surface relative humidity
    real(wp), allocatable, dimension(:,:,:) :: ra2
    real(wp), allocatable, dimension(:,:,:) :: alb_vu_s  !! visible+UV clear sky surface albedo (.)
    real(wp), allocatable, dimension(:,:,:) :: alb_vu_c  !! visible+UV cloudy sky surface albedo (.)
    real(wp), allocatable, dimension(:,:,:) :: alb_ir_s  !! infrared clear sky surface albedo (.)
    real(wp), allocatable, dimension(:,:,:) :: alb_ir_c  !! infrared cloudy sky surface albedo (.)
    real(wp), allocatable, dimension(:,:,:) :: cd

    real(wp), allocatable, dimension(:,:) :: tskina     !! skin temperature (K)
    real(wp), allocatable, dimension(:,:) :: t2a        !! 2m surface air temperature (K)
    real(wp), allocatable, dimension(:,:) :: t2a_dat    !! observed present-day 2m surface air temperature (K)
    real(wp), allocatable, dimension(:,:) :: q2a        !! 2m surface specific humidity (kg/kg)
    real(wp), allocatable, dimension(:,:) :: thetae
    real(wp), allocatable, dimension(:,:) :: dq 
    real(wp), allocatable, dimension(:,:) :: dr
    real(wp), allocatable, dimension(:,:) :: r2a        !! 2m surface relative humidity
    real(wp), allocatable, dimension(:,:) :: rskina 
    real(wp), allocatable, dimension(:,:) :: ra2a
    real(wp), allocatable, dimension(:,:) :: cda
    real(wp), allocatable, dimension(:,:) :: cd0a
    real(wp), allocatable, dimension(:,:) :: sha
    real(wp), allocatable, dimension(:,:) :: lha
    real(wp), allocatable, dimension(:,:) :: evpa
    real(wp), allocatable, dimension(:,:) :: Ri

    real(wp), allocatable, dimension(:,:) :: hdust      !! dust height scale (m)
    real(wp), allocatable, dimension(:,:) :: dust_load  !! dust load (kg/m2)
    real(wp), allocatable, dimension(:,:) :: dust_emis  !! dust emissions (kg/m2/s) 
    real(wp), allocatable, dimension(:,:) :: dust_dep   !! dust deposition (kg/m2/s)
    real(wp), allocatable, dimension(:,:) :: dust_dep_dry !! dust dry deposition (kg/m2/s)
    real(wp), allocatable, dimension(:,:) :: dust_dep_wet !! dust wet deposition (kg/m2/s)
    real(wp), allocatable, dimension(:,:) :: dust_ot    !! dust optical thickness (.)

    real(wp), allocatable, dimension(:,:) :: cam    
    real(wp), allocatable, dimension(:,:) :: co2d
    real(wp), allocatable, dimension(:,:) :: co2flx    

    real(wp), allocatable, dimension(:,:) :: so4_ot    !! SO4 optical thickness (.)

    ! 3D longitude-latitude-height
    real(wp), allocatable, dimension(:,:,:) :: t3   !! atmospheric temperature (K)
    real(wp), allocatable, dimension(:,:,:) :: q3   !! atmospheric specific humidity (kg/kg)
    real(wp), allocatable, dimension(:,:,:) :: r3   !! atmospheric relative humidity
    real(wp), allocatable, dimension(:,:,:) :: tp   !! potential temperature
    real(wp), allocatable, dimension(:,:,:) :: gamma   !! atmospheric temperature lapse rate (K/m)
    real(wp), allocatable, dimension(:,:,:) :: rho 

    real(wp), allocatable, dimension(:,:) :: acbar
    real(wp), allocatable, dimension(:,:) :: epsa
    real(wp), allocatable, dimension(:,:) :: slp
    real(wp), allocatable, dimension(:,:) :: slp1
    real(wp), allocatable, dimension(:,:) :: tsksl
    real(wp), allocatable, dimension(:,:) :: atsksl
    real(wp), allocatable, dimension(:,:) :: atsl
    real(wp), allocatable, dimension(:,:) :: atsli
    real(wp), allocatable, dimension(:,:) :: aslp
    real(wp), allocatable, dimension(:,:) :: aslp_temp
    real(wp), allocatable, dimension(:,:) :: aslp_topo
    real(wp), allocatable, dimension(:,:) :: dz500
    real(wp), allocatable, dimension(:,:) :: us
    real(wp), allocatable, dimension(:,:) :: vs
    real(wp), allocatable, dimension(:,:) :: usk
    real(wp), allocatable, dimension(:,:) :: vsk
    real(wp), allocatable, dimension(:,:) :: ugb
    real(wp), allocatable, dimension(:,:) :: vgb
    real(wp), allocatable, dimension(:,:) :: uab
    real(wp), allocatable, dimension(:,:) :: vab
    real(wp), allocatable, dimension(:,:,:) :: taux
    real(wp), allocatable, dimension(:,:,:) :: tauy
    real(wp), allocatable, dimension(:,:) :: wcld
    real(wp), allocatable, dimension(:,:) :: woro
    real(wp), allocatable, dimension(:,:) :: wsyn       !! vertical synoptic wind component (m/s)
    real(wp), allocatable, dimension(:,:) :: weff    
    real(wp), allocatable, dimension(:,:) :: fweff    

    real(wp), allocatable, dimension(:,:,:) :: ua
    real(wp), allocatable, dimension(:,:,:) :: va
    real(wp), allocatable, dimension(:,:,:) :: u3
    real(wp), allocatable, dimension(:,:,:) :: v3
    real(wp), allocatable, dimension(:,:,:) :: w3
    real(wp), allocatable, dimension(:,:,:) :: uter
    real(wp), allocatable, dimension(:,:,:) :: vter

    real(wp), allocatable, dimension(:,:,:) :: fax
    real(wp), allocatable, dimension(:,:,:) :: faxo
    real(wp), allocatable, dimension(:,:,:) :: fay
    real(wp), allocatable, dimension(:,:,:) :: fayo
    real(wp), allocatable, dimension(:,:) :: fac

    real(wp), allocatable, dimension(:,:) :: xz

    real(wp), allocatable, dimension(:,:) :: convdse
    real(wp), allocatable, dimension(:,:) :: convadse
    real(wp), allocatable, dimension(:,:) :: convddse
    real(wp), allocatable, dimension(:,:) :: convwtr
    real(wp), allocatable, dimension(:,:) :: convawtr
    real(wp), allocatable, dimension(:,:) :: convdwtr
    real(wp), allocatable, dimension(:,:) :: faxmas
    real(wp), allocatable, dimension(:,:) :: faymas
    real(wp), allocatable, dimension(:,:) :: faxdse
    real(wp), allocatable, dimension(:,:) :: faxcpt
    real(wp), allocatable, dimension(:,:) :: faxwtr
    real(wp), allocatable, dimension(:,:) :: faydse
    real(wp), allocatable, dimension(:,:) :: faycpt
    real(wp), allocatable, dimension(:,:) :: faywtr
    real(wp), allocatable, dimension(:,:) :: fdxdse
    real(wp), allocatable, dimension(:,:) :: fdxwtr
    real(wp), allocatable, dimension(:,:) :: fdydse
    real(wp), allocatable, dimension(:,:) :: fdywtr
    real(wp), allocatable, dimension(:) :: fayg         ! zonal integral of meridional mass flux
    real(wp), allocatable, dimension(:) :: faydseg      ! zonal integral of meridional dry static energy flux by mean circulation (W)
    real(wp), allocatable, dimension(:) :: faycptg      ! zonal integral of meridional thermal energy flux by mean circulation (W)
    real(wp), allocatable, dimension(:) :: faygzg       ! zonal integral of meridional geopotential energy flux by mean circulation (W)
    real(wp), allocatable, dimension(:) :: fayleg       ! zonal integral of meridional latent energy flux by mean circulation (W)
    real(wp), allocatable, dimension(:) :: faywtrg      ! zonal integral of meridional water flux by mean circulation (kg/s)
    real(wp), allocatable, dimension(:) :: fdydseg      ! zonal integral of meridional dry static energy flux by eddies (W)
    real(wp), allocatable, dimension(:) :: fdyleg       ! zonal integral of meridional latent energy flux by eddies (W)
    real(wp), allocatable, dimension(:) :: fdywtrg      ! zonal integral of meridional water flux by eddies (kg/s)
    real(wp), allocatable, dimension(:) :: fydseg       ! zonal integral of meridional dry static energy flux (W)
    real(wp), allocatable, dimension(:) :: fyheatg      ! zonal integral of meridional heat flux (W)
    real(wp), allocatable, dimension(:) :: fyleg        ! zonal integral of meridional latent energy flux (W)
    real(wp), allocatable, dimension(:) :: fywtrg       ! zonal integral of meridional moisture flux (kg/s) 

    real(wp), allocatable, dimension(:,:,:) :: fswr_sur
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

    real(wp), allocatable, dimension(:,:) :: rb_top
    real(wp), allocatable, dimension(:,:) :: rb_sur
    real(wp), allocatable, dimension(:,:) :: rb_atm
    real(wp), allocatable, dimension(:,:) :: rb_str     
    real(wp), allocatable, dimension(:,:) :: swr_atm
    real(wp), allocatable, dimension(:,:) :: swr_top
    real(wp), allocatable, dimension(:,:) :: swr_top_cs
    real(wp), allocatable, dimension(:,:) :: swr_sur
    real(wp), allocatable, dimension(:,:) :: swr_sur_cs
    real(wp), allocatable, dimension(:,:) :: swr_dw_sur
    real(wp), allocatable, dimension(:,:) :: swr_dw_top
    real(wp), allocatable, dimension(:,:) :: lwr_atm
    real(wp), allocatable, dimension(:,:) :: lwr_atm_cld
    real(wp), allocatable, dimension(:,:) :: lwr_top
    real(wp), allocatable, dimension(:,:) :: lwr_top_cs
    real(wp), allocatable, dimension(:,:) :: lwr_sur
    real(wp), allocatable, dimension(:,:) :: lwr_sur_cs
    real(wp), allocatable, dimension(:,:) :: lwr_tro
    real(wp), allocatable, dimension(:,:) :: lwr_cld
    real(wp), allocatable, dimension(:,:) :: cre_top
    real(wp), allocatable, dimension(:,:) :: swr_cre_top
    real(wp), allocatable, dimension(:,:) :: lwr_cre_top
    real(wp), allocatable, dimension(:,:) :: cre_sur
    real(wp), allocatable, dimension(:,:) :: swr_cre_sur
    real(wp), allocatable, dimension(:,:) :: lwr_cre_sur
    real(wp), allocatable, dimension(:,:) :: aplan      !! planetary albedo

    real(wp), allocatable, dimension(:,:) :: tsl

    real(wp), allocatable, dimension(:,:) :: eke
    real(wp), allocatable, dimension(:,:) :: sam
    real(wp), allocatable, dimension(:,:) :: sam2
    real(wp), allocatable, dimension(:,:) :: synprod
    real(wp), allocatable, dimension(:,:) :: syndiss
    real(wp), allocatable, dimension(:,:) :: synadv
    real(wp), allocatable, dimension(:,:) :: syndif
    real(wp), allocatable, dimension(:,:) :: synsur
    real(wp), allocatable, dimension(:,:) :: cdif     
    real(wp), allocatable, dimension(:,:) :: diffxdse
    real(wp), allocatable, dimension(:,:) :: diffydse
    real(wp), allocatable, dimension(:,:) :: diffxwtr
    real(wp), allocatable, dimension(:,:) :: diffywtr

    real(wp), allocatable, dimension(:) :: fw_pac_atl 
    real(wp), allocatable, dimension(:) :: fw_atl_indpac

  end type

  type(ts_out), allocatable :: ann_ts(:)
  type(a_out) :: day_a(nday_year), mon_a(nmon_year), ann_a


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for atmosphere
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_diag_init

    implicit none

    integer :: k
    real(wp), parameter :: lon_grl = -42.5_wp   ! NGRIP
    real(wp), parameter :: lat_grl = 75.1_wp    ! NGRIP
    real(wp), parameter :: lon_grl_1=-61_wp, lon_grl_2=-21._wp, lat_grl_1=61._wp , lat_grl_2=79._wp
    real(wp), parameter :: lon_natl_1=-51_wp, lon_natl_2=9._wp, lat_natl_1=61._wp , lat_natl_2=79._wp
    real(wp), parameter :: lon_ant = 0._wp 
    real(wp), parameter :: lat_ant = -82.5_wp 
    real(wp), parameter :: lon_wais = -112.5_wp 
    real(wp), parameter :: lat_wais = -77.5_wp 


    nout = 0

    ! initialise netcdf output
    call ts_nc(trim(out_dir)//"/atm_ts.nc")
    call atm_nc(trim(out_dir)//"/atm.nc")
    if (l_daily_output) call atm_daily_nc(trim(out_dir)//"/atm_daily.nc")

    ! allocate
    allocate(ann_ts(ny_out_ts))

     allocate(ann_a%solar(nday_year,24,jm))
     allocate(ann_a%cosz(nday_year,24,jm))
     allocate(ann_a%solarm(nday_year,jm))
     allocate(ann_a%coszm(nday_year,jm))

     allocate(ann_a%ptrop(jm))
     allocate(ann_a%tslz(jm))
     allocate(ann_a%tskslz(jm))
     allocate(ann_a%ekez(jm))
     allocate(ann_a%slpz(jm))
     allocate(ann_a%vabz(jmc))
     allocate(ann_a%acbarz(jm))
     allocate(ann_a%uz850(jm))
     allocate(ann_a%uz500(jm))

     allocate(ann_a%zs(im,jm,nm))
     allocate(ann_a%zsa(im,jm))
     allocate(ann_a%zsa_smooth(im,jm))
     allocate(ann_a%slope(im,jm))
     allocate(ann_a%slope_x(im,jm))
     allocate(ann_a%slope_y(im,jm))
     allocate(ann_a%frlnd(im,jm))
     allocate(ann_a%frocn(im,jm))

     allocate(ann_a%tam(im,jm))  
     allocate(ann_a%dtamdt(im,jm))  
     allocate(ann_a%qam(im,jm))  
     allocate(ann_a%ram(im,jm))  
     allocate(ann_a%gams(im,jm)) 
     allocate(ann_a%gamb(im,jm)) 
     allocate(ann_a%gamt(im,jm)) 
     allocate(ann_a%dam(im,jm))  
     allocate(ann_a%hqeff(im,jm))  
     allocate(ann_a%hrm(im,jm))  
     allocate(ann_a%wcon(im,jm)) 
     allocate(ann_a%cld_rh(im,jm)) 
     allocate(ann_a%cld_low(im,jm)) 
     allocate(ann_a%cld(im,jm)) 
     allocate(ann_a%prc(im,jm)) 
     allocate(ann_a%prcw(im,jm))
     allocate(ann_a%prcs(im,jm))
     allocate(ann_a%prc_conv(im,jm))
     allocate(ann_a%prc_wcon(im,jm))
     allocate(ann_a%prc_over(im,jm))
     allocate(ann_a%hcld(im,jm))
     allocate(ann_a%ctt(im,jm)) 
     allocate(ann_a%clot(im,jm)) 
     allocate(ann_a%alb_cld(im,jm)) 
     allocate(ann_a%alb_sur_cs(im,jm)) 
     allocate(ann_a%alb_sur_cld(im,jm)) 
     allocate(ann_a%htrop(im,jm)) 
     allocate(ann_a%ttrop(im,jm))  
     allocate(ann_a%wind(im,jm))   

     allocate(ann_a%frst(im,jm,nm))
     allocate(ann_a%tskin(im,jm,nm)) 
     allocate(ann_a%t2(im,jm,nm)) 
     allocate(ann_a%q2(im,jm,nm))  
     allocate(ann_a%r2(im,jm,nm))   
     allocate(ann_a%ra2(im,jm,nm))
     allocate(ann_a%alb_vu_s(im,jm,nm))
     allocate(ann_a%alb_vu_c(im,jm,nm))
     allocate(ann_a%alb_ir_s(im,jm,nm))
     allocate(ann_a%alb_ir_c(im,jm,nm))
     allocate(ann_a%cd(im,jm,nm))

     allocate(ann_a%tskina(im,jm)) 
     allocate(ann_a%t2a(im,jm)) 
     allocate(ann_a%t2a_dat(im,jm)) 
     allocate(ann_a%q2a(im,jm))  
     allocate(ann_a%thetae(im,jm)) 
     allocate(ann_a%dq(im,jm))  
     allocate(ann_a%dr(im,jm))  
     allocate(ann_a%r2a(im,jm))   
     allocate(ann_a%rskina(im,jm))   
     allocate(ann_a%ra2a(im,jm))
     allocate(ann_a%cda(im,jm))
     allocate(ann_a%cd0a(im,jm))
     allocate(ann_a%sha(im,jm))
     allocate(ann_a%lha(im,jm))
     allocate(ann_a%evpa(im,jm))
     allocate(ann_a%Ri(im,jm))
 
     allocate(ann_a%hdust(im,jm)) 
     allocate(ann_a%dust_load(im,jm)) 
     allocate(ann_a%dust_emis(im,jm)) 
     allocate(ann_a%dust_dep(im,jm)) 
     allocate(ann_a%dust_dep_dry(im,jm)) 
     allocate(ann_a%dust_dep_wet(im,jm)) 
     allocate(ann_a%dust_ot(im,jm)) 

     allocate(ann_a%cam(im,jm)) 
     allocate(ann_a%co2d(im,jm)) 
     allocate(ann_a%co2flx(im,jm)) 

     allocate(ann_a%so4_ot(im,jm)) 

     allocate(ann_a%t3(im,jm,kmc))
     allocate(ann_a%q3(im,jm,kmc))
     allocate(ann_a%r3(im,jm,kmc))
     allocate(ann_a%tp(im,jm,kmc)) 
     allocate(ann_a%gamma(im,jm,km))
     allocate(ann_a%rho(im,jm,kmc))
 
     allocate(ann_a%acbar(im,jm))
     allocate(ann_a%epsa(im,jm))
     allocate(ann_a%tsksl(im,jm))
     allocate(ann_a%atsksl(im,jm))
     allocate(ann_a%atsl(im,jm))
     allocate(ann_a%atsli(im,jm))
     allocate(ann_a%aslp(im,jm))
     allocate(ann_a%aslp_temp(im,jm))
     allocate(ann_a%aslp_topo(im,jm))
     allocate(ann_a%dz500(im,jm))
     allocate(ann_a%slp(im,jm))
     allocate(ann_a%slp1(im,jm))
     allocate(ann_a%us(im,jm))
     allocate(ann_a%vs(im,jm))
     allocate(ann_a%usk(im,jm))
     allocate(ann_a%vsk(im,jm))
     allocate(ann_a%ugb(im,jm))
     allocate(ann_a%vgb(im,jm))
     allocate(ann_a%uab(imc,jm))
     allocate(ann_a%vab(im,jmc))
     allocate(ann_a%taux(im,jm,nm))
     allocate(ann_a%tauy(im,jm,nm))
     allocate(ann_a%wcld(im,jm))
     allocate(ann_a%woro(im,jm))
     allocate(ann_a%wsyn(im,jm))
     allocate(ann_a%weff(im,jm))   
     allocate(ann_a%fweff(im,jm))   

     allocate(ann_a%ua(im,jm,km))
     allocate(ann_a%va(im,jm,km))
     allocate(ann_a%u3(im,jm,km))
     allocate(ann_a%v3(im,jm,km))
     allocate(ann_a%w3(im,jm,kmc))
     allocate(ann_a%uter(im,jm,km))
     allocate(ann_a%vter(im,jm,km))

     if (l_output_flx3d) then
       allocate(ann_a%fax(imc,jm,km))
       allocate(ann_a%faxo(imc,jm,km))
       allocate(ann_a%fay(im,jmc,km))
       allocate(ann_a%fayo(im,jmc,km))
       allocate(ann_a%fac(im,jm))
     endif

     allocate(ann_a%xz(jmc,kmc))

     allocate(ann_a%convdse(im,jm))
     allocate(ann_a%convadse(im,jm))
     allocate(ann_a%convddse(im,jm))
     allocate(ann_a%convwtr(im,jm))
     allocate(ann_a%convawtr(im,jm))
     allocate(ann_a%convdwtr(im,jm))
     allocate(ann_a%faxmas(im,jm))
     allocate(ann_a%faymas(im,jm))
     allocate(ann_a%faxdse(im,jm))
     allocate(ann_a%faxcpt(im,jm))
     allocate(ann_a%faxwtr(im,jm))
     allocate(ann_a%faydse(im,jm))
     allocate(ann_a%faycpt(im,jm))
     allocate(ann_a%faywtr(im,jm))
     allocate(ann_a%fdxdse(im,jm))
     allocate(ann_a%fdxwtr(im,jm))
     allocate(ann_a%fdydse(im,jm))
     allocate(ann_a%fdywtr(im,jm))
     allocate(ann_a%fayg(jmc))
     allocate(ann_a%faydseg(jmc))
     allocate(ann_a%faycptg(jmc))
     allocate(ann_a%faygzg(jmc))
     allocate(ann_a%fayleg(jmc))
     allocate(ann_a%faywtrg(jmc))
     allocate(ann_a%fdydseg(jmc))
     allocate(ann_a%fdyleg(jmc))
     allocate(ann_a%fdywtrg(jmc))
     allocate(ann_a%fydseg(jmc))
     allocate(ann_a%fyheatg(jmc))
     allocate(ann_a%fyleg(jmc))
     allocate(ann_a%fywtrg(jmc))

     allocate(ann_a%fswr_sur(im,jm,nm))
     allocate(ann_a%flwr_dw_sur(im,jm,nm))
     allocate(ann_a%flwr_dw_sur_cs(im,jm,nm))
     allocate(ann_a%flwr_dw_sur_cld(im,jm,nm))
     allocate(ann_a%flwr_up_sur(im,jm,nm))

     allocate(ann_a%dswd_dalb_vu_cs(im,jm))
     allocate(ann_a%dswd_dalb_ir_cs(im,jm))
     allocate(ann_a%dswd_dalb_vu_cld(im,jm))
     allocate(ann_a%dswd_dalb_ir_cld(im,jm))
     allocate(ann_a%dswd_dz_ir_cs(im,jm))
     allocate(ann_a%dswd_dz_ir_cld(im,jm))

     allocate(ann_a%rb_top(im,jm))
     allocate(ann_a%rb_sur(im,jm))
     allocate(ann_a%rb_atm(im,jm))
     allocate(ann_a%swr_atm(im,jm))
     allocate(ann_a%swr_top(im,jm))
     allocate(ann_a%swr_top_cs(im,jm))
     allocate(ann_a%swr_sur(im,jm))
     allocate(ann_a%swr_sur_cs(im,jm))
     allocate(ann_a%swr_dw_sur(im,jm))
     allocate(ann_a%swr_dw_top(im,jm))
     allocate(ann_a%lwr_atm(im,jm))
     allocate(ann_a%lwr_atm_cld(im,jm))
     allocate(ann_a%lwr_top(im,jm))
     allocate(ann_a%lwr_top_cs(im,jm))
     allocate(ann_a%lwr_sur(im,jm))
     allocate(ann_a%lwr_sur_cs(im,jm))
     allocate(ann_a%lwr_tro(im,jm))
     allocate(ann_a%lwr_cld(im,jm))
     allocate(ann_a%rb_str(im,jm))     
     allocate(ann_a%cre_top(im,jm))
     allocate(ann_a%swr_cre_top(im,jm))
     allocate(ann_a%lwr_cre_top(im,jm))
     allocate(ann_a%cre_sur(im,jm))
     allocate(ann_a%swr_cre_sur(im,jm))
     allocate(ann_a%lwr_cre_sur(im,jm))
     allocate(ann_a%aplan(im,jm))

     allocate(ann_a%tsl(im,jm))

     allocate(ann_a%eke(im,jm))
     allocate(ann_a%sam(im,jm))
     allocate(ann_a%sam2(im,jm))
     allocate(ann_a%synprod(im,jm))
     allocate(ann_a%syndiss(im,jm))
     allocate(ann_a%synadv(im,jm))
     allocate(ann_a%syndif(im,jm))
     allocate(ann_a%synsur(im,jm))
     allocate(ann_a%cdif(im,jm) )    
     allocate(ann_a%diffxdse(imc,jm))
     allocate(ann_a%diffydse(im,jmc))
     allocate(ann_a%diffxwtr(imc,jm))
     allocate(ann_a%diffywtr(im,jmc))

     allocate(ann_a%fw_pac_atl(jm))    
     allocate(ann_a%fw_atl_indpac(jm))    

    do k=1,nmon_year
     allocate(mon_a(k)%ptrop(jm))
     allocate(mon_a(k)%tslz(jm))
     allocate(mon_a(k)%tskslz(jm))
     allocate(mon_a(k)%ekez(jm))
     allocate(mon_a(k)%slpz(jm))
     allocate(mon_a(k)%vabz(jmc))
     allocate(mon_a(k)%acbarz(jm))
     allocate(mon_a(k)%uz850(jm))
     allocate(mon_a(k)%uz500(jm))
     allocate(mon_a(k)%tam(im,jm))  
     allocate(mon_a(k)%dtamdt(im,jm))  
     allocate(mon_a(k)%qam(im,jm))  
     allocate(mon_a(k)%ram(im,jm))  
     allocate(mon_a(k)%gams(im,jm)) 
     allocate(mon_a(k)%gamb(im,jm)) 
     allocate(mon_a(k)%gamt(im,jm)) 
     allocate(mon_a(k)%dam(im,jm))  
     allocate(mon_a(k)%hqeff(im,jm))  
     allocate(mon_a(k)%hrm(im,jm))  
     allocate(mon_a(k)%wcon(im,jm)) 
     allocate(mon_a(k)%cld_rh(im,jm)) 
     allocate(mon_a(k)%cld_low(im,jm)) 
     allocate(mon_a(k)%cld(im,jm)) 
     allocate(mon_a(k)%prc(im,jm)) 
     allocate(mon_a(k)%prcw(im,jm))
     allocate(mon_a(k)%prcs(im,jm))
     allocate(mon_a(k)%prc_conv(im,jm))
     allocate(mon_a(k)%prc_wcon(im,jm))
     allocate(mon_a(k)%prc_over(im,jm))
     allocate(mon_a(k)%hcld(im,jm))
     allocate(mon_a(k)%ctt(im,jm)) 
     allocate(mon_a(k)%clot(im,jm)) 
     allocate(mon_a(k)%alb_cld(im,jm)) 
     allocate(mon_a(k)%alb_sur_cs(im,jm)) 
     allocate(mon_a(k)%alb_sur_cld(im,jm)) 
     allocate(mon_a(k)%htrop(im,jm)) 
     allocate(mon_a(k)%ttrop(im,jm))  
     allocate(mon_a(k)%wind(im,jm))   

     allocate(mon_a(k)%frst(im,jm,nm)) 
     allocate(mon_a(k)%tskin(im,jm,nm)) 
     allocate(mon_a(k)%t2(im,jm,nm)) 
     allocate(mon_a(k)%q2(im,jm,nm))  
     allocate(mon_a(k)%r2(im,jm,nm))   
     allocate(mon_a(k)%ra2(im,jm,nm))
     allocate(mon_a(k)%alb_vu_s(im,jm,nm))
     allocate(mon_a(k)%alb_vu_c(im,jm,nm))
     allocate(mon_a(k)%alb_ir_s(im,jm,nm))
     allocate(mon_a(k)%alb_ir_c(im,jm,nm))
     allocate(mon_a(k)%cd(im,jm,nm))

     allocate(mon_a(k)%tskina(im,jm)) 
     allocate(mon_a(k)%t2a(im,jm)) 
     allocate(mon_a(k)%t2a_dat(im,jm)) 
     allocate(mon_a(k)%thetae(im,jm)) 
     allocate(mon_a(k)%q2a(im,jm))  
     allocate(mon_a(k)%dq(im,jm))  
     allocate(mon_a(k)%dr(im,jm))  
     allocate(mon_a(k)%r2a(im,jm))   
     allocate(mon_a(k)%rskina(im,jm))   
     allocate(mon_a(k)%ra2a(im,jm))
     allocate(mon_a(k)%cd0a(im,jm))
     allocate(mon_a(k)%cda(im,jm))
     allocate(mon_a(k)%sha(im,jm))
     allocate(mon_a(k)%lha(im,jm))
     allocate(mon_a(k)%evpa(im,jm))
     allocate(mon_a(k)%Ri(im,jm))
 
     allocate(mon_a(k)%hdust(im,jm)) 
     allocate(mon_a(k)%dust_load(im,jm)) 
     allocate(mon_a(k)%dust_emis(im,jm)) 
     allocate(mon_a(k)%dust_dep(im,jm)) 
     allocate(mon_a(k)%dust_dep_dry(im,jm)) 
     allocate(mon_a(k)%dust_dep_wet(im,jm)) 
     allocate(mon_a(k)%dust_ot(im,jm)) 

     allocate(mon_a(k)%cam(im,jm)) 
     allocate(mon_a(k)%co2d(im,jm)) 
     allocate(mon_a(k)%co2flx(im,jm)) 

     allocate(mon_a(k)%so4_ot(im,jm)) 

     allocate(mon_a(k)%t3(im,jm,kmc))
     allocate(mon_a(k)%q3(im,jm,kmc))
     allocate(mon_a(k)%r3(im,jm,kmc))
     allocate(mon_a(k)%tp(im,jm,kmc)) 
     allocate(mon_a(k)%gamma(im,jm,km))
     allocate(mon_a(k)%rho(im,jm,kmc))
 
     allocate(mon_a(k)%acbar(im,jm))
     allocate(mon_a(k)%epsa(im,jm))
     allocate(mon_a(k)%tsksl(im,jm))
     allocate(mon_a(k)%atsksl(im,jm))
     allocate(mon_a(k)%atsl(im,jm))
     allocate(mon_a(k)%atsli(im,jm))
     allocate(mon_a(k)%aslp(im,jm))
     allocate(mon_a(k)%aslp_temp(im,jm))
     allocate(mon_a(k)%aslp_topo(im,jm))
     allocate(mon_a(k)%dz500(im,jm))
     allocate(mon_a(k)%slp(im,jm))
     allocate(mon_a(k)%slp1(im,jm))
     allocate(mon_a(k)%us(im,jm))
     allocate(mon_a(k)%vs(im,jm))
     allocate(mon_a(k)%usk(im,jm))
     allocate(mon_a(k)%vsk(im,jm))
     allocate(mon_a(k)%ugb(im,jm))
     allocate(mon_a(k)%vgb(im,jm))
     allocate(mon_a(k)%uab(imc,jm))
     allocate(mon_a(k)%vab(im,jmc))
     allocate(mon_a(k)%taux(im,jm,nm))
     allocate(mon_a(k)%tauy(im,jm,nm))
     allocate(mon_a(k)%wcld(im,jm))
     allocate(mon_a(k)%woro(im,jm))
     allocate(mon_a(k)%wsyn(im,jm))
     allocate(mon_a(k)%weff(im,jm))   
     allocate(mon_a(k)%fweff(im,jm))   

     allocate(mon_a(k)%ua(im,jm,km))
     allocate(mon_a(k)%va(im,jm,km))
     allocate(mon_a(k)%u3(im,jm,km))
     allocate(mon_a(k)%v3(im,jm,km))
     allocate(mon_a(k)%w3(im,jm,kmc))
     allocate(mon_a(k)%uter(im,jm,km))
     allocate(mon_a(k)%vter(im,jm,km))

     if (l_output_flx3d) then
       allocate(mon_a(k)%fax(imc,jm,km))
       allocate(mon_a(k)%faxo(imc,jm,km))
       allocate(mon_a(k)%fay(im,jmc,km))
       allocate(mon_a(k)%fayo(im,jmc,km))
       allocate(mon_a(k)%fac(im,jm))
     endif

     allocate(mon_a(k)%xz(jmc,kmc))

     allocate(mon_a(k)%convdse(im,jm))
     allocate(mon_a(k)%convadse(im,jm))
     allocate(mon_a(k)%convddse(im,jm))
     allocate(mon_a(k)%convwtr(im,jm))
     allocate(mon_a(k)%convawtr(im,jm))
     allocate(mon_a(k)%convdwtr(im,jm))
     allocate(mon_a(k)%faxmas(im,jm))
     allocate(mon_a(k)%faymas(im,jm))
     allocate(mon_a(k)%faxdse(im,jm))
     allocate(mon_a(k)%faxcpt(im,jm))
     allocate(mon_a(k)%faxwtr(im,jm))
     allocate(mon_a(k)%faydse(im,jm))
     allocate(mon_a(k)%faycpt(im,jm))
     allocate(mon_a(k)%faywtr(im,jm))
     allocate(mon_a(k)%fdxdse(im,jm))
     allocate(mon_a(k)%fdxwtr(im,jm))
     allocate(mon_a(k)%fdydse(im,jm))
     allocate(mon_a(k)%fdywtr(im,jm))
     allocate(mon_a(k)%fayg(jmc))
     allocate(mon_a(k)%faydseg(jmc))
     allocate(mon_a(k)%faycptg(jmc))
     allocate(mon_a(k)%faygzg(jmc))
     allocate(mon_a(k)%fayleg(jmc))
     allocate(mon_a(k)%faywtrg(jmc))
     allocate(mon_a(k)%fdydseg(jmc))
     allocate(mon_a(k)%fdyleg(jmc))
     allocate(mon_a(k)%fdywtrg(jmc))
     allocate(mon_a(k)%fydseg(jmc))
     allocate(mon_a(k)%fyheatg(jmc))
     allocate(mon_a(k)%fyleg(jmc))
     allocate(mon_a(k)%fywtrg(jmc))

     allocate(mon_a(k)%fswr_sur(im,jm,nm))
     allocate(mon_a(k)%flwr_dw_sur(im,jm,nm))
     allocate(mon_a(k)%flwr_dw_sur_cs(im,jm,nm))
     allocate(mon_a(k)%flwr_dw_sur_cld(im,jm,nm))
     allocate(mon_a(k)%flwr_up_sur(im,jm,nm))

     allocate(mon_a(k)%dswd_dalb_vu_cs(im,jm))
     allocate(mon_a(k)%dswd_dalb_ir_cs(im,jm))
     allocate(mon_a(k)%dswd_dalb_vu_cld(im,jm))
     allocate(mon_a(k)%dswd_dalb_ir_cld(im,jm))
     allocate(mon_a(k)%dswd_dz_ir_cs(im,jm))
     allocate(mon_a(k)%dswd_dz_ir_cld(im,jm))

     allocate(mon_a(k)%rb_top(im,jm))
     allocate(mon_a(k)%rb_sur(im,jm))
     allocate(mon_a(k)%rb_atm(im,jm))
     allocate(mon_a(k)%swr_atm(im,jm))
     allocate(mon_a(k)%swr_top(im,jm))
     allocate(mon_a(k)%swr_top_cs(im,jm))
     allocate(mon_a(k)%swr_sur(im,jm))
     allocate(mon_a(k)%swr_sur_cs(im,jm))
     allocate(mon_a(k)%swr_dw_sur(im,jm))
     allocate(mon_a(k)%swr_dw_top(im,jm))
     allocate(mon_a(k)%lwr_atm(im,jm))
     allocate(mon_a(k)%lwr_atm_cld(im,jm))
     allocate(mon_a(k)%lwr_top(im,jm))
     allocate(mon_a(k)%lwr_top_cs(im,jm))
     allocate(mon_a(k)%lwr_sur(im,jm))
     allocate(mon_a(k)%lwr_sur_cs(im,jm))
     allocate(mon_a(k)%lwr_tro(im,jm))
     allocate(mon_a(k)%lwr_cld(im,jm))
     allocate(mon_a(k)%rb_str(im,jm))     
     allocate(mon_a(k)%cre_top(im,jm))
     allocate(mon_a(k)%swr_cre_top(im,jm))
     allocate(mon_a(k)%lwr_cre_top(im,jm))
     allocate(mon_a(k)%cre_sur(im,jm))
     allocate(mon_a(k)%swr_cre_sur(im,jm))
     allocate(mon_a(k)%lwr_cre_sur(im,jm))
     allocate(mon_a(k)%aplan(im,jm))

     allocate(mon_a(k)%tsl(im,jm))

     allocate(mon_a(k)%eke(im,jm))
     allocate(mon_a(k)%sam(im,jm))
     allocate(mon_a(k)%sam2(im,jm))
     allocate(mon_a(k)%synprod(im,jm))
     allocate(mon_a(k)%syndiss(im,jm))
     allocate(mon_a(k)%synadv(im,jm))
     allocate(mon_a(k)%syndif(im,jm))
     allocate(mon_a(k)%synsur(im,jm))
     allocate(mon_a(k)%cdif(im,jm) )    
     allocate(mon_a(k)%diffxdse(imc,jm))
     allocate(mon_a(k)%diffydse(im,jmc))
     allocate(mon_a(k)%diffxwtr(imc,jm))
     allocate(mon_a(k)%diffywtr(im,jmc))

     allocate(mon_a(k)%fw_pac_atl(jm))    
     allocate(mon_a(k)%fw_atl_indpac(jm))    

    enddo


    ! find i,j indexes for Greenland and Antarctica
    i_grl = minloc(abs(lon-lon_grl),1)
    j_grl = minloc(abs(-lat-lat_grl),1)
    i_ant = minloc(abs(lon-lon_ant),1)
    j_ant = minloc(abs(-lat-lat_ant),1)
    i_grl_1 = minloc(abs(lon-lon_grl_1),1)
    i_grl_2 = minloc(abs(lon-lon_grl_2),1)
    j_grl_1 = minloc(abs(-lat-lat_grl_1),1)
    j_grl_2 = minloc(abs(-lat-lat_grl_2),1)
    i_natl_1 = minloc(abs(lon-lon_natl_1),1)
    i_natl_2 = minloc(abs(lon-lon_natl_2),1)
    j_natl_1 = minloc(abs(-lat-lat_natl_1),1)
    j_natl_2 = minloc(abs(-lat-lat_natl_2),1)
    i_wais = minloc(abs(lon-lon_wais),1)
    j_wais = minloc(abs(-lat-lat_wais),1)


   return

  end subroutine atm_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a t m _ d i a g
  !   Purpose    :  atmosphere diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_diag(atm)

    implicit none

    type(atm_class), intent(in) :: atm

    integer :: i, j, m, k, y, ind, imi, jmi, ipl, jpl
    character (len=256) :: fnm
    real(wp) :: mon_avg, ann_avg, jja_avg, favg, favgl, favg6090, fcum, xsum, lsqr, esqr6090
    real(wp) :: t, rh
    real(wp) :: faxmasi, faymasi
    real(wp) :: tup
    real(wp) :: ps, slp
    real(wp) :: ctt

    real(wp) :: fayg(jmc)
    real(wp) :: faydseg(jmc)
    real(wp) :: faycptg(jmc)
    real(wp) :: faygzg(jmc)
    real(wp) :: fayleg(jmc)
    real(wp) :: faywtrg(jmc)
    real(wp) :: fdydseg(jmc)
    real(wp) :: fdyleg(jmc)
    real(wp) :: fdywtrg(jmc)
    real(wp) :: fydseg(jmc)
    real(wp) :: fyheatg(jmc)
    real(wp) :: fyleg(jmc)
    real(wp) :: fywtrg(jmc)
    real(wp) :: fw_pac_atl(jm)
    real(wp) :: fw_atl_indpac(jm)

    real(wp) :: ann_co2d_avg
    real(wp) :: ann_co2flx
    real(wp) :: ann_co2d_wais
    real(wp) :: ann_t2m
    real(wp) :: ann_t2mnh
    real(wp) :: ann_t2msh
    real(wp) :: ann_t2mn6090
    real(wp) :: ann_t2ms6090
    real(wp) :: ann_q2m
    real(wp) :: ann_r2m
    real(wp) :: ann_tam
    real(wp) :: ann_tskin
    real(wp) :: ann_prc
    real(wp) :: ann_evp
    real(wp) :: ann_wcon
    real(wp) :: ann_cld
    real(wp) :: ann_sha
    real(wp) :: ann_lha
    real(wp) :: ann_rbtop
    real(wp) :: ann_rbatm
    real(wp) :: ann_rbsur
    real(wp) :: ann_snettop
    real(wp) :: ann_snetsur
    real(wp) :: ann_lnettop
    real(wp) :: ann_lnetsur
    real(wp) :: ann_ldwnsur
    real(wp) :: ann_rbtopcs
    real(wp) :: ann_rbsurcs
    real(wp) :: ann_snettopcs
    real(wp) :: ann_snetsurcs
    real(wp) :: ann_lnettopcs
    real(wp) :: ann_lnetsurcs
    real(wp) :: ann_ldwnsurcs
    real(wp) :: ann_lupsur 
    real(wp) :: ann_slp
    real(wp) :: ann_htrop
    real(wp) :: ann_ttrop
    real(wp) :: ann_hcld
    real(wp) :: ann_plan
    real(wp) :: ann_sol
    real(wp) :: ann_sreftop
    real(wp) :: ann_swrcrf
    real(wp) :: ann_lwrcrf
    real(wp) :: ann_dustload
    real(wp) :: ann_dustdep
    real(wp) :: ann_dustdrydep
    real(wp) :: ann_dustwetdep
    real(wp) :: ann_fw_pac_atl
    real(wp) :: ann_fw_atl_indpac
    real(wp) :: ann_t2m_l
    real(wp) :: ann_q2m_l
    real(wp) :: ann_r2m_l
    real(wp) :: ann_tam_l
    real(wp) :: ann_tskin_l
    real(wp) :: ann_prc_l
    real(wp) :: ann_evp_l
    real(wp) :: ann_wcon_l
    real(wp) :: ann_cld_l
    real(wp) :: ann_sha_l
    real(wp) :: ann_lha_l
    real(wp) :: ann_rbtop_l
    real(wp) :: ann_rbatm_l
    real(wp) :: ann_rbsur_l
    real(wp) :: ann_snettop_l 
    real(wp) :: ann_snetsur_l 
    real(wp) :: ann_lnettop_l
    real(wp) :: ann_lnetsur_l

    real(wp), allocatable, dimension(:,:) :: xz
    real(wp), allocatable, dimension(:,:) :: faxcpt
    real(wp), allocatable, dimension(:,:) :: faycpt

    real(wp) :: tsl(im,jm)

    ! current index
    y = y_out_ts_clim

    mon_avg = 1._wp/nstep_mon_atm
    ann_avg = 1._wp/nstep_year_atm
    jja_avg = 1._wp/(nstep_mon_atm*3._wp)

    if (time_soy_atm) then
      ann_ts(y)%co2flx=0
      ann_ts(y)%co2d_avg=0
      ann_ts(y)%co2d_wais=0
      ann_ts(y)%t2m=0
      ann_ts(y)%tnh=0
      ann_ts(y)%tsh=0
      ann_ts(y)%dtnsh=0
      ann_ts(y)%tn6090=0
      ann_ts(y)%ts6090=0
      ann_ts(y)%tgrl=0
      ann_ts(y)%tgrl1=0
      ann_ts(y)%tgrl1_jja=0
      ann_ts(y)%tnatl=0
      ann_ts(y)%tant=0
      ann_ts(y)%dust_grl=0
      ann_ts(y)%dust_wais=0
      ann_ts(y)%q2m=0
      ann_ts(y)%r2m=0
      ann_ts(y)%tam=0
      ann_ts(y)%tskin=0
      ann_ts(y)%prc=0
      ann_ts(y)%evp=0
      ann_ts(y)%wcon=0
      ann_ts(y)%cld=0
      ann_ts(y)%sha = 0
      ann_ts(y)%lha = 0
      ann_ts(y)%rbtop = 0
      ann_ts(y)%rbatm = 0
      ann_ts(y)%rbsur = 0
      ann_ts(y)%snettop =0
      ann_ts(y)%snetsur =0
      ann_ts(y)%lnettop =0
      ann_ts(y)%lnetsur =0
      ann_ts(y)%ldwnsur =0
      ann_ts(y)%rbtopcs = 0
      ann_ts(y)%rbatmcs = 0
      ann_ts(y)%rbsurcs = 0
      ann_ts(y)%snettopcs =0
      ann_ts(y)%snetsurcs =0
      ann_ts(y)%lnettopcs =0
      ann_ts(y)%lnetsurcs =0
      ann_ts(y)%ldwnsurcs =0
      ann_ts(y)%lupsur  =0
      ann_ts(y)%slp=0
      ann_ts(y)%htrop=0
      ann_ts(y)%ttrop=0
      ann_ts(y)%hcld=0      
      ann_ts(y)%plan=0
      ann_ts(y)%sol=0
      ann_ts(y)%sreftop =0
      ann_ts(y)%swrcrf =0    
      ann_ts(y)%lwrcrf =0     
      ann_ts(y)%dustload=0      
      ann_ts(y)%dustdep=0      
      ann_ts(y)%dustdrydep=0      
      ann_ts(y)%dustwetdep=0      
      ann_ts(y)%fw_pac_atl =0     
      ann_ts(y)%fw_atl_indpac =0     
      ann_ts(y)%t2m_l=0
      ann_ts(y)%q2m_l=0
      ann_ts(y)%r2m_l=0
      ann_ts(y)%tam_l=0
      ann_ts(y)%tskin_l=0
      ann_ts(y)%prc_l=0
      ann_ts(y)%evp_l=0
      ann_ts(y)%wcon_l=0
      ann_ts(y)%cld_l=0
      ann_ts(y)%sha_l = 0
      ann_ts(y)%lha_l = 0
      ann_ts(y)%rbtop_l = 0
      ann_ts(y)%rbatm_l = 0
      ann_ts(y)%rbsur_l = 0
      ann_ts(y)%snettop_l =0
      ann_ts(y)%snetsur_l =0
      ann_ts(y)%lnettop_l =0
      ann_ts(y)%lnetsur_l =0
      ann_ts(y)%itcz =0
      ann_ts(y)%itcz_min =0
      ann_ts(y)%itcz_max =0
      ann_ts(y)%had_width =0
    endif

    ann_co2d_avg=0
    ann_co2flx=0
    ann_co2d_wais=0
    ann_t2m=0
    ann_t2mnh=0
    ann_t2msh=0
    ann_t2mn6090=0
    ann_t2ms6090=0
    ann_q2m=0
    ann_r2m=0
    ann_tam=0
    ann_tskin=0
    ann_prc=0
    ann_evp=0
    ann_wcon=0
    ann_cld=0
    ann_sha = 0
    ann_lha = 0
    ann_rbtop = 0
    ann_rbatm = 0
    ann_rbsur = 0
    ann_snettop =0
    ann_snetsur =0
    ann_lnettop =0
    ann_lnetsur =0
    ann_ldwnsur =0
    ann_rbtopcs = 0
    ann_rbsurcs = 0
    ann_snettopcs =0
    ann_snetsurcs =0
    ann_lnettopcs =0
    ann_lnetsurcs =0
    ann_ldwnsurcs =0
    ann_lupsur  =0
    ann_slp=0
    ann_htrop=0
    ann_ttrop=0
    ann_hcld=0      
    ann_plan=0
    ann_sol=0
    ann_sreftop =0
    ann_swrcrf =0    
    ann_lwrcrf =0     
    ann_dustload=0      
    ann_dustdep=0      
    ann_dustdrydep=0      
    ann_dustwetdep=0      
    ann_fw_pac_atl =0     
    ann_fw_atl_indpac =0     
    ann_t2m_l=0
    ann_q2m_l=0
    ann_r2m_l=0
    ann_tam_l=0
    ann_tskin_l=0
    ann_prc_l=0
    ann_evp_l=0
    ann_wcon_l=0
    ann_cld_l=0
    ann_sha_l = 0
    ann_lha_l = 0
    ann_rbtop_l = 0
    ann_rbatm_l = 0
    ann_rbsur_l = 0
    ann_snettop_l =0
    ann_snetsur_l =0
    ann_lnettop_l =0
    ann_lnetsur_l =0

    lsqr = sum(sqr(:,:)*atm%frlnd(:,:))
    esqr6090 = sum(sqr(1,:),mask=-lat.ge.60._wp)*im

    !$omp parallel do private(i,j,ind,favg,favgl,fcum) reduction(+:ann_co2d_avg, ann_co2flx, ann_co2d_wais, ann_t2m ,ann_t2mnh) &
    !$omp reduction(+:ann_t2msh ,ann_t2mn6090 ,ann_t2ms6090 ,ann_q2m ,ann_r2m ,ann_tam ,ann_tskin ,ann_prc ,ann_evp ,ann_wcon ,ann_cld ,ann_sha) &
    !$omp reduction(+:ann_lha ,ann_rbtop ,ann_rbatm ,ann_rbsur ,ann_snettop ,ann_snetsur ,ann_lnettop ,ann_lnetsur ,ann_ldwnsur, ann_rbtopcs) &
    !$omp reduction(+:ann_rbsurcs ,ann_snettopcs ,ann_snetsurcs ,ann_lnettopcs ,ann_lnetsurcs ,ann_ldwnsurcs ,ann_lupsur ,ann_slp ,ann_htrop) &
    !$omp reduction(+:ann_ttrop, ann_hcld ,ann_plan ,ann_sol ,ann_sreftop ,ann_swrcrf ,ann_lwrcrf ,ann_dustload ,ann_dustdep ,ann_dustdrydep) &
    !$omp reduction(+:ann_dustwetdep ,ann_fw_pac_atl ,ann_fw_atl_indpac ,ann_t2m_l ,ann_q2m_l ,ann_r2m_l ,ann_tam_l ,ann_tskin_l ,ann_prc_l) &
    !$omp reduction(+:ann_evp_l ,ann_wcon_l ,ann_cld_l ,ann_sha_l ,ann_lha_l ,ann_rbtop_l ,ann_rbatm_l ,ann_rbsur_l ,ann_snettop_l ,ann_snetsur_l) &
    !$omp reduction(+:ann_lnettop_l ,ann_lnetsur_l)
    do j=1,jm
      do i=1,im

        favg = ann_avg * sqr(i,j)/esqr
        if (lsqr.gt.0._wp) then
          favgl = ann_avg * sqr(i,j)*atm%frlnd(i,j)/lsqr
        else
          favgl = 0._wp
        endif
        favg6090 = ann_avg * sqr(i,j)/esqr6090
        fcum = sqr(i,j)

        ann_co2d_avg   = ann_co2d_avg   + atm%cam(i,j)*1.e6*28.97/44.0095*favg ! ppm
        ann_co2flx     = ann_co2flx     + atm%co2flx(i,j)*12._wp/44.0095*sqr(i,j)*sec_day   ! kgCO2/m2/s * kgC/kgCO2 * m2 * s = kgC
        ann_t2m        = ann_t2m        + sum(atm%t2(i,j,:)*atm%frst(i,j,:))*favg
        if (-lat(j).ge.0._wp) then
          ann_t2mnh      = ann_t2mnh      + sum(atm%t2(i,j,:)*atm%frst(i,j,:))*favg*2._wp
        else
          ann_t2msh      = ann_t2msh      + sum(atm%t2(i,j,:)*atm%frst(i,j,:))*favg*2._wp
        endif
        if (-lat(j).ge.60._wp) then
          ann_t2mn6090 = ann_t2mn6090 + sum(atm%t2(i,j,:)*atm%frst(i,j,:))*favg6090
        endif
        if (lat(j).ge.60._wp) then
          ann_t2ms6090 = ann_t2ms6090 + sum(atm%t2(i,j,:)*atm%frst(i,j,:))*favg6090
        endif
        ann_q2m        = ann_q2m        + sum(atm%q2(i,j,:)*atm%frst(i,j,:))*favg
        ann_r2m        = ann_r2m        + sum(atm%r2(i,j,:)*atm%frst(i,j,:))*favg
        ann_tskin      = ann_tskin      + sum(atm%tskin(i,j,:)*atm%frst(i,j,:))*favg
        ann_tam        = ann_tam        + atm%tam(i,j)*favg
        ann_prc        = ann_prc        + atm%prc(i,j)*sec_day*fcum
        ann_evp        = ann_evp        + atm%evpa(i,j)*sec_day*fcum
        ann_wcon       = ann_wcon       + atm%wcon(i,j)*favg
        ann_cld        = ann_cld        + atm%cld(i,j)*favg
        ann_sha        = ann_sha        + atm%sha(i,j)*favg
        ann_lha        = ann_lha        + atm%lha(i,j)*favg
        ann_rbtop      = ann_rbtop      + atm%rb_top(i,j)*favg
        ann_rbatm      = ann_rbatm      + atm%rb_atm(i,j)*favg
        ann_rbsur      = ann_rbsur      + atm%rb_sur(i,j)*favg
        ann_snettop    = ann_snettop    + atm%swr_top(i,j)*favg
        ann_snetsur    = ann_snetsur    + sum(atm%fswr_sur(i,j,:)*atm%frst(i,j,:))*favg
        ann_lnettop    = ann_lnettop    + atm%lwr_top(i,j)*favg
        ann_lnetsur    = ann_lnetsur    + atm%lwr_sur(i,j)*favg
        ann_ldwnsur    = ann_ldwnsur    + sum(atm%flwr_dw_sur(i,j,:)*atm%frst(i,j,:))*favg
        ann_rbtopcs    = ann_rbtopcs    + (atm%swr_top_cs(i,j)+atm%lwr_top_cs(i,j))*favg
        ann_rbsurcs    = ann_rbsurcs    + (sum(atm%fswr_sur_cs(i,j,:)*atm%frst(i,j,:)) &
          +sum(atm%flwr_dw_sur_cs(i,j,:)*atm%frst(i,j,:))-sum(atm%flwr_up_sur(i,j,:)*atm%frst(i,j,:)))*favg
        ann_snettopcs  = ann_snettopcs  + atm%swr_top_cs(i,j)*favg
        ann_snetsurcs  = ann_snetsurcs  + sum(atm%fswr_sur_cs(i,j,:)*atm%frst(i,j,:))*favg
        ann_lnettopcs  = ann_lnettopcs  + atm%lwr_top_cs(i,j)*favg
        ann_lnetsurcs  = ann_lnetsurcs  + (sum(atm%flwr_dw_sur_cs(i,j,:)*atm%frst(i,j,:))-sum(atm%flwr_up_sur(i,j,:)*atm%frst(i,j,:)))*favg
        ann_ldwnsurcs  = ann_ldwnsurcs  + sum(atm%flwr_dw_sur_cs(i,j,:)*atm%frst(i,j,:))*favg
        ann_lupsur     = ann_lupsur     + sum(atm%flwr_up_sur(i,j,:)*atm%frst(i,j,:))*favg
        ann_slp        = ann_slp        + atm%slp(i,j)*1.e-2*favg
        ann_sol        = ann_sol        + atm%solarm(doy,j)*favg
        ann_sreftop    = ann_sreftop    + (atm%solarm(doy,j)-atm%swr_top(i,j))*favg
        ann_swrcrf     = ann_swrcrf     + (atm%swr_top_cld(i,j)-atm%swr_top_cs(i,j))*atm%cld(i,j)*favg      
        ann_lwrcrf     = ann_lwrcrf     + (atm%lwr_top_cld(i,j)-atm%lwr_top_cs(i,j))*atm%cld(i,j)*favg      
        ann_htrop      = ann_htrop      + atm%htrop(i,j)*favg
        ann_ttrop      = ann_ttrop      + atm%ttrop(i,j)*favg
        ann_hcld       = ann_hcld       + atm%hcld(i,j)*favg
        ann_dustload   = ann_dustload   + atm%dust_load(i,j)*fcum*1.e-12  ! Tg
        ann_dustdep    = ann_dustdep    + atm%dust_dep(i,j)*fcum*sec_day*1.e-12 ! Tg/yr
        ann_dustdrydep = ann_dustdrydep + atm%dust_dep_dry(i,j)*fcum*sec_day*1.e-12 ! Tg/yr
        ann_dustwetdep = ann_dustwetdep + atm%dust_dep_wet(i,j)*fcum*sec_day*1.e-12 ! Tg/yr

        ann_t2m_l      = ann_t2m_l     + sum(atm%t2(i,j,:)*atm%frst(i,j,:))*favgl
        ann_q2m_l      = ann_q2m_l     + sum(atm%q2(i,j,:)*atm%frst(i,j,:))*favgl
        ann_r2m_l      = ann_r2m_l     + sum(atm%r2(i,j,:)*atm%frst(i,j,:))*favgl
        ann_tam_l      = ann_tam_l     + atm%tam(i,j)*favgl
        ann_tskin_l    = ann_tskin_l   + sum(atm%tskin(i,j,:)*atm%frst(i,j,:))*favgl
        ann_prc_l      = ann_prc_l     + atm%prc(i,j)*sec_day*sqr(i,j)*atm%frlnd(i,j)
        ann_evp_l      = ann_evp_l     + atm%evpa(i,j)*sec_day*sqr(i,j)*atm%frlnd(i,j)
        ann_wcon_l     = ann_wcon_l    + atm%wcon(i,j)*favgl
        ann_cld_l      = ann_cld_l     + atm%cld(i,j)*favgl
        ann_sha_l      = ann_sha_l     + atm%sha(i,j)*favgl
        ann_lha_l      = ann_lha_l     + atm%lha(i,j)*favgl
        ann_rbtop_l    = ann_rbtop_l   + atm%rb_top(i,j)*favgl
        ann_rbatm_l    = ann_rbatm_l   + atm%rb_atm(i,j)*favgl
        ann_rbsur_l    = ann_rbsur_l   + atm%rb_sur(i,j)*favgl
        ann_snettop_l  = ann_snettop_l + atm%swr_top(i,j)*favgl
        ann_snetsur_l  = ann_snetsur_l + sum(atm%fswr_sur(i,j,:)*atm%frst(i,j,:))*favgl
        ann_lnettop_l  = ann_lnettop_l + atm%lwr_top(i,j)*favgl
        ann_lnetsur_l  = ann_lnetsur_l + atm%lwr_sur(i,j)*favgl
      enddo
      if (atm%idivide_pac_atl(j).gt.0) then
        ind = atm%idivide_pac_atl(j)+1
        if (ind.eq.im+1) ind=1
        ann_fw_pac_atl = ann_fw_pac_atl + (atm%faxwtr(ind,j)+atm%fdxwtr(ind,j))*ann_avg ! kg/s
      endif
      if (atm%idivide_atl_indpac(j).gt.0) then
        ind = atm%idivide_atl_indpac(j)+1
        if (ind.eq.im+1) ind=1
        ann_fw_atl_indpac = ann_fw_atl_indpac + (atm%faxwtr(ind,j)+atm%fdxwtr(ind,j))*ann_avg   ! kg/s
      endif
    enddo
    !$omp end parallel do
        
    ann_ts(y)%co2d_avg      = ann_ts(y)%co2d_avg      + ann_co2d_avg        
    ann_ts(y)%co2flx        = ann_ts(y)%co2flx        + ann_co2flx*1e-12    ! kgC -> PgC
    ann_ts(y)%co2d_wais     = ann_ts(y)%co2d_wais     + atm%cam(i_wais,j_wais)*1.e6*28.97/44.0095 * ann_avg ! ppm
    ann_ts(y)%t2m           = ann_ts(y)%t2m           + ann_t2m        
    ann_ts(y)%tnh           = ann_ts(y)%tnh           + ann_t2mnh        
    ann_ts(y)%tsh           = ann_ts(y)%tsh           + ann_t2msh        
    ann_ts(y)%dtnsh         = ann_ts(y)%dtnsh         + (ann_t2mnh-ann_t2msh)
    ann_ts(y)%tn6090        = ann_ts(y)%tn6090        + ann_t2mn6090        
    ann_ts(y)%ts6090        = ann_ts(y)%ts6090        + ann_t2ms6090        
    ann_ts(y)%tgrl          = ann_ts(y)%tgrl          + sum(atm%t2(i_grl,j_grl,:)*atm%frst(i_grl,j_grl,:)) * ann_avg
    ann_ts(y)%tgrl1         = ann_ts(y)%tgrl1         + sum(atm%t2(i_grl_1:i_grl_2,j_grl_2:j_grl_1,i_ice)*atm%frst(i_grl_1:i_grl_2,j_grl_2:j_grl_1,i_ice)) &
                                                      / (sum(atm%frst(i_grl_1:i_grl_2,j_grl_2:j_grl_1,i_ice))+1.e-20_wp)* ann_avg
    if (mon.ge.6 .and. mon.le.8) then
      ann_ts(y)%tgrl1_jja   = ann_ts(y)%tgrl1_jja     + sum(atm%t2(i_grl_1:i_grl_2,j_grl_2:j_grl_1,i_ice)*atm%frst(i_grl_1:i_grl_2,j_grl_2:j_grl_1,i_ice)) &
                                                      / (sum(atm%frst(i_grl_1:i_grl_2,j_grl_2:j_grl_1,i_ice))+1.e-20_wp)* jja_avg
    endif
    ann_ts(y)%tnatl         = ann_ts(y)%tnatl         + sum(atm%t2(i_natl_1:i_natl_2,j_natl_2:j_natl_1,:)*atm%frst(i_natl_1:i_natl_2,j_natl_2:j_natl_1,:)) &
                                                      / size(atm%frst(i_natl_1:i_natl_2,j_natl_2:j_natl_1,1))* ann_avg
    ann_ts(y)%tant          = ann_ts(y)%tant          + sum(atm%t2(i_ant,j_ant,:)*atm%frst(i_ant,j_ant,:)) * ann_avg 
    ann_ts(y)%dust_grl      = ann_ts(y)%dust_grl      + atm%dust_dep(i_grl,j_grl) * ann_avg 
    ann_ts(y)%dust_wais     = ann_ts(y)%dust_wais     + atm%dust_dep(i_wais,j_wais) * ann_avg 
    ann_ts(y)%q2m           = ann_ts(y)%q2m           + ann_q2m        
    ann_ts(y)%r2m           = ann_ts(y)%r2m           + ann_r2m        
    ann_ts(y)%tskin         = ann_ts(y)%tskin         + ann_tskin      
    ann_ts(y)%tam           = ann_ts(y)%tam           + ann_tam        
    ann_ts(y)%prc           = ann_ts(y)%prc           + ann_prc        
    ann_ts(y)%evp           = ann_ts(y)%evp           + ann_evp        
    ann_ts(y)%wcon          = ann_ts(y)%wcon          + ann_wcon       
    ann_ts(y)%cld           = ann_ts(y)%cld           + ann_cld        
    ann_ts(y)%sha           = ann_ts(y)%sha           + ann_sha        
    ann_ts(y)%lha           = ann_ts(y)%lha           + ann_lha        
    ann_ts(y)%rbtop         = ann_ts(y)%rbtop         + ann_rbtop      
    ann_ts(y)%rbatm         = ann_ts(y)%rbatm         + ann_rbatm      
    ann_ts(y)%rbsur         = ann_ts(y)%rbsur         + ann_rbsur      
    ann_ts(y)%snettop       = ann_ts(y)%snettop       + ann_snettop    
    ann_ts(y)%snetsur       = ann_ts(y)%snetsur       + ann_snetsur    
    ann_ts(y)%lnettop       = ann_ts(y)%lnettop       + ann_lnettop    
    ann_ts(y)%lnetsur       = ann_ts(y)%lnetsur       + ann_lnetsur    
    ann_ts(y)%ldwnsur       = ann_ts(y)%ldwnsur       + ann_ldwnsur    
    ann_ts(y)%rbtopcs       = ann_ts(y)%rbtopcs       + ann_rbtopcs      
    ann_ts(y)%rbatmcs       = ann_ts(y)%rbatmcs       + (ann_rbtopcs-ann_rbsurcs)      
    ann_ts(y)%rbsurcs       = ann_ts(y)%rbsurcs       + ann_rbsurcs      
    ann_ts(y)%snettopcs     = ann_ts(y)%snettopcs     + ann_snettopcs    
    ann_ts(y)%snetsurcs     = ann_ts(y)%snetsurcs     + ann_snetsurcs    
    ann_ts(y)%lnettopcs     = ann_ts(y)%lnettopcs     + ann_lnettopcs    
    ann_ts(y)%lnetsurcs     = ann_ts(y)%lnetsurcs     + ann_lnetsurcs    
    ann_ts(y)%ldwnsurcs     = ann_ts(y)%ldwnsurcs     + ann_ldwnsurcs    
    ann_ts(y)%lupsur        = ann_ts(y)%lupsur        + ann_lupsur     
    ann_ts(y)%slp           = ann_ts(y)%slp           + ann_slp        
    ann_ts(y)%sol           = ann_ts(y)%sol           + ann_sol        
    ann_ts(y)%sreftop       = ann_ts(y)%sreftop       + ann_sreftop    
    ann_ts(y)%swrcrf        = ann_ts(y)%swrcrf        + ann_swrcrf     
    ann_ts(y)%lwrcrf        = ann_ts(y)%lwrcrf        + ann_lwrcrf     
    ann_ts(y)%htrop         = ann_ts(y)%htrop         + ann_htrop      
    ann_ts(y)%ttrop         = ann_ts(y)%ttrop         + ann_ttrop      
    ann_ts(y)%hcld          = ann_ts(y)%hcld          + ann_hcld       
    ann_ts(y)%dustload      = ann_ts(y)%dustload      + ann_dustload   
    ann_ts(y)%dustdep       = ann_ts(y)%dustdep       + ann_dustdep    
    ann_ts(y)%dustdrydep    = ann_ts(y)%dustdrydep    + ann_dustdrydep 
    ann_ts(y)%dustwetdep    = ann_ts(y)%dustwetdep    + ann_dustwetdep 
    ann_ts(y)%fw_pac_atl    = ann_ts(y)%fw_pac_atl    + ann_fw_pac_atl 
    ann_ts(y)%fw_atl_indpac = ann_ts(y)%fw_atl_indpac + ann_fw_atl_indpac 

    ann_ts(y)%t2m_l         = ann_ts(y)%t2m_l         + ann_t2m_l     
    ann_ts(y)%q2m_l         = ann_ts(y)%q2m_l         + ann_q2m_l     
    ann_ts(y)%r2m_l         = ann_ts(y)%r2m_l         + ann_r2m_l     
    ann_ts(y)%tam_l         = ann_ts(y)%tam_l         + ann_tam_l     
    ann_ts(y)%tskin_l       = ann_ts(y)%tskin_l       + ann_tskin_l  
    ann_ts(y)%prc_l         = ann_ts(y)%prc_l         + ann_prc_l     
    ann_ts(y)%evp_l         = ann_ts(y)%evp_l         + ann_evp_l     
    ann_ts(y)%wcon_l        = ann_ts(y)%wcon_l        + ann_wcon_l    
    ann_ts(y)%cld_l         = ann_ts(y)%cld_l         + ann_cld_l     
    ann_ts(y)%sha_l         = ann_ts(y)%sha_l         + ann_sha_l     
    ann_ts(y)%lha_l         = ann_ts(y)%lha_l         + ann_lha_l     
    ann_ts(y)%rbtop_l       = ann_ts(y)%rbtop_l       + ann_rbtop_l   
    ann_ts(y)%rbatm_l       = ann_ts(y)%rbatm_l       + ann_rbatm_l   
    ann_ts(y)%rbsur_l       = ann_ts(y)%rbsur_l       + ann_rbsur_l   
    ann_ts(y)%snettop_l     = ann_ts(y)%snettop_l     + ann_snettop_l 
    ann_ts(y)%snetsur_l     = ann_ts(y)%snetsur_l     + ann_snetsur_l 
    ann_ts(y)%lnettop_l     = ann_ts(y)%lnettop_l     + ann_lnettop_l 
    ann_ts(y)%lnetsur_l     = ann_ts(y)%lnetsur_l     + ann_lnetsur_l 
    ann_ts(y)%itcz          = ann_ts(y)%itcz          + atm%had_fi*180./pi * ann_avg
    ann_ts(y)%itcz_min      = min(ann_ts(y)%itcz_min, atm%had_fi*180./pi)
    ann_ts(y)%itcz_max      = max(ann_ts(y)%itcz_max, atm%had_fi*180./pi)
    ann_ts(y)%had_width     = ann_ts(y)%had_width     + atm%had_width*180./pi * ann_avg

    if( time_eoy_atm ) then

      ann_ts(y)%eccentricity = atm%eccentricity
      ann_ts(y)%precession   = atm%eccentricity*atm%precession
      ann_ts(y)%obliquity    = atm%obliquity

      ann_ts(y)%Smax65N = maxval(0.5_wp*(atm%solarm(:,minloc(abs(-lat-65._wp)))+atm%solarm(:,minloc(abs(-lat-65._wp),back=.true.))))
      ann_ts(y)%Smax65S = maxval(0.5_wp*(atm%solarm(:,minloc(abs(-lat+65._wp)))+atm%solarm(:,minloc(abs(-lat+65._wp),back=.true.))))
      ann_ts(y)%Smax30N = maxval(0.5_wp*(atm%solarm(:,minloc(abs(-lat-30._wp)))+atm%solarm(:,minloc(abs(-lat-30._wp),back=.true.))))
      ann_ts(y)%Smax30S = maxval(0.5_wp*(atm%solarm(:,minloc(abs(-lat+30._wp)))+atm%solarm(:,minloc(abs(-lat+30._wp),back=.true.))))

      ann_ts(y)%co2 = atm%co2
      ann_ts(y)%co2e= atm%co2e
      ann_ts(y)%ch4 = atm%ch4
      ann_ts(y)%n2o = atm%n2o
      ann_ts(y)%so4 = sum(atm%so4(:,:)*sqr(:,:))    ! kg SO4

      ann_ts(y)%t2m = ann_ts(y)%t2m - T0    ! degC
      ann_ts(y)%tnh = ann_ts(y)%tnh - T0    ! degC
      ann_ts(y)%tsh = ann_ts(y)%tsh - T0    ! degC
      ann_ts(y)%dtnsh = ann_ts(y)%dtnsh - T0    ! degC
      ann_ts(y)%tn6090 = ann_ts(y)%tn6090 - T0    ! degC
      ann_ts(y)%ts6090 = ann_ts(y)%ts6090 - T0    ! degC
      ann_ts(y)%tgrl = ann_ts(y)%tgrl - T0    ! degC
      ann_ts(y)%tgrl1 = ann_ts(y)%tgrl1 - T0    ! degC
      ann_ts(y)%tgrl1_jja = ann_ts(y)%tgrl1_jja - T0    ! degC
      ann_ts(y)%tnatl = ann_ts(y)%tnatl - T0    ! degC
      ann_ts(y)%tant = ann_ts(y)%tant - T0    ! degC
      ann_ts(y)%tam = ann_ts(y)%tam - T0    ! degC
      ann_ts(y)%tskin = ann_ts(y)%tskin - T0    ! degC
      ann_ts(y)%prc = ann_ts(y)%prc*1.e-15  ! 10^15 kg
      ann_ts(y)%evp = ann_ts(y)%evp*1.e-15  ! 10^15 kg

      ann_ts(y)%t2m_l = ann_ts(y)%t2m_l - T0    ! degC
      ann_ts(y)%tam_l = ann_ts(y)%tam_l - T0    ! degC
      ann_ts(y)%tskin_l = ann_ts(y)%tskin_l - T0    ! degC
      ann_ts(y)%prc_l = ann_ts(y)%prc_l*1.e-15  ! 10^15 kg
      ann_ts(y)%evp_l = ann_ts(y)%evp_l*1.e-15  ! 10^15 kg

      ann_ts(y)%plan = ann_ts(y)%sreftop/ann_ts(y)%sol  

      ann_ts(y)%fw_pac_atl = ann_ts(y)%fw_pac_atl * 1.e-9 ! kg/s -> Sv 
      ann_ts(y)%fw_atl_indpac = ann_ts(y)%fw_atl_indpac * 1.e-9 ! kg/s -> Sv 

     ! write to file
     if (time_out_ts_clim) then
       fnm = trim(out_dir)//"/atm_ts.nc"
       call ts_nc_write(fnm,ann_ts(1:y),year_clim-y+1,y)
     endif

     ! print header
     if (mod(year_clim,10).eq.1) then
        print '(a7,a9,17a7)','atm','year','S65N','CO2','CH4','N2O','t2m','prc','evp','alb','sol','swnett','swnets','lwnett','lwnets','lwdwns','rbtop','rbatm','rbsur'
     endif

     ! print values
     print '(a7,i9,7F7.1,1F7.3,6F7.1,F7.2,2F7.1)', &
     'atm',year_now,ann_ts(y)%Smax65N,ann_ts(y)%co2,ann_ts(y)%ch4,ann_ts(y)%n2o,ann_ts(y)%t2m,ann_ts(y)%prc,ann_ts(y)%evp, &
     ann_ts(y)%plan,ann_ts(y)%sol,ann_ts(y)%snettop,ann_ts(y)%snetsur,ann_ts(y)%lnettop,ann_ts(y)%lnetsur,ann_ts(y)%ldwnsur, &
     ann_ts(y)%rbtop,ann_ts(y)%rbatm,ann_ts(y)%rbsur

    endif


    ! spatially explicit output
    if ( time_out_atm ) then

      if (time_soy_atm) then
        do m=1,nmon_year
          mon_a(m)%had_fi      = 0. 
          mon_a(m)%had_width   = 0. 
          mon_a(m)%ptrop       = 0. 
          mon_a(m)%tslz        = 0. 
          mon_a(m)%tskslz       = 0. 
          mon_a(m)%ekez        = 0. 
          mon_a(m)%slpz        = 0. 
          mon_a(m)%vabz        = 0. 
          mon_a(m)%acbarz      = 0. 
          mon_a(m)%uz850       = 0. 
          mon_a(m)%uz500       = 0. 
          mon_a(m)%tam         = 0. 
          mon_a(m)%dtamdt      = 0. 
          mon_a(m)%qam         = 0. 
          mon_a(m)%ram         = 0. 
          mon_a(m)%gams        = 0. 
          mon_a(m)%gamb        = 0. 
          mon_a(m)%gamt        = 0. 
          mon_a(m)%dam         = 0. 
          mon_a(m)%hqeff       = 0. 
          mon_a(m)%hrm         = 0. 
          mon_a(m)%wcon        = 0. 
          mon_a(m)%cld_rh         = 0. 
          mon_a(m)%cld_low         = 0. 
          mon_a(m)%cld         = 0. 
          mon_a(m)%prc         = 0. 
          mon_a(m)%prcw        = 0. 
          mon_a(m)%prcs        = 0. 
          mon_a(m)%prc_conv    = 0. 
          mon_a(m)%prc_wcon    = 0. 
          mon_a(m)%prc_over    = 0. 
          mon_a(m)%hcld        = 0. 
          mon_a(m)%ctt         = 0. 
          mon_a(m)%clot        = 0. 
          mon_a(m)%alb_cld     = 0. 
          mon_a(m)%alb_sur_cs  = 0. 
          mon_a(m)%alb_sur_cld = 0. 
          mon_a(m)%htrop       = 0. 
          mon_a(m)%ttrop       = 0. 
          mon_a(m)%wind        = 0. 
          mon_a(m)%frst        = 0. 
          mon_a(m)%tskin       = 0. 
          mon_a(m)%t2          = 0. 
          mon_a(m)%q2          = 0. 
          mon_a(m)%r2          = 0. 
          mon_a(m)%ra2         = 0. 
          mon_a(m)%alb_vu_s    = 0. 
          mon_a(m)%alb_vu_c    = 0. 
          mon_a(m)%alb_ir_s    = 0. 
          mon_a(m)%alb_ir_c    = 0. 
          mon_a(m)%tskina      = 0. 
          mon_a(m)%t2a         = 0. 
          mon_a(m)%thetae      = 0. 
          mon_a(m)%q2a         = 0. 
          mon_a(m)%dq          = 0. 
          mon_a(m)%dr          = 0. 
          mon_a(m)%r2a         = 0. 
          mon_a(m)%rskina      = 0. 
          mon_a(m)%ra2a        = 0. 
          mon_a(m)%cd0a        = 0. 
          mon_a(m)%cda         = 0. 
          mon_a(m)%cd          = 0. 
          mon_a(m)%sha         = 0. 
          mon_a(m)%lha         = 0. 
          mon_a(m)%evpa        = 0. 
          mon_a(m)%Ri          = 0. 
          mon_a(m)%t3          = 0. 
          mon_a(m)%q3          = 0. 
          mon_a(m)%r3          = 0. 
          mon_a(m)%tp          = 0. 
          mon_a(m)%gamma       = 0. 
          mon_a(m)%rho         = 0. 
          mon_a(m)%acbar       = 0. 
          mon_a(m)%epsa        = 0. 
          mon_a(m)%tsksl      = 0. 
          mon_a(m)%atsksl      = 0. 
          mon_a(m)%atsl      = 0. 
          mon_a(m)%atsli       = 0. 
          mon_a(m)%aslp        = 0. 
          mon_a(m)%aslp_temp   = 0. 
          mon_a(m)%aslp_topo   = 0. 
          mon_a(m)%dz500       = 0. 
          mon_a(m)%slp         = 0. 
          mon_a(m)%slp1        = 0. 
          mon_a(m)%us         = 0. 
          mon_a(m)%vs         = 0. 
          mon_a(m)%usk        = 0. 
          mon_a(m)%vsk        = 0. 
          mon_a(m)%ugb        = 0. 
          mon_a(m)%vgb        = 0. 
          mon_a(m)%uab        = 0. 
          mon_a(m)%vab        = 0. 
          mon_a(m)%taux        = 0. 
          mon_a(m)%tauy        = 0. 
          mon_a(m)%wcld        = 0. 
          mon_a(m)%woro        = 0. 
          mon_a(m)%weff        = 0. 
          mon_a(m)%fweff       = 0. 
          mon_a(m)%ua         = 0. 
          mon_a(m)%va         = 0. 
          mon_a(m)%u3         = 0. 
          mon_a(m)%v3         = 0. 
          mon_a(m)%w3         = 0. 
          mon_a(m)%uter        = 0. 
          mon_a(m)%vter        = 0. 
          if (l_output_flx3d) then
            mon_a(m)%fax        = 0. 
            mon_a(m)%faxo       = 0. 
            mon_a(m)%fay        = 0. 
            mon_a(m)%fayo       = 0. 
            mon_a(m)%fac        = 0. 
          endif
          mon_a(m)%diffxdse    = 0. 
          mon_a(m)%diffydse    = 0. 
          mon_a(m)%diffxwtr    = 0. 
          mon_a(m)%diffywtr    = 0. 
          mon_a(m)%wsyn        = 0. 
          mon_a(m)%convdse     = 0. 
          mon_a(m)%convadse    = 0. 
          mon_a(m)%convddse    = 0. 
          mon_a(m)%convwtr     = 0. 
          mon_a(m)%convawtr    = 0. 
          mon_a(m)%convdwtr    = 0. 
          mon_a(m)%faxmas     = 0.   
          mon_a(m)%faymas     = 0.   
          mon_a(m)%faxdse     = 0.   
          mon_a(m)%faxcpt     = 0.   
          mon_a(m)%faxwtr     = 0.
          mon_a(m)%faydse     = 0.
          mon_a(m)%faycpt     = 0.
          mon_a(m)%faywtr     = 0.
          mon_a(m)%fdxdse     = 0.
          mon_a(m)%fdxwtr     = 0.
          mon_a(m)%fdydse     = 0.
          mon_a(m)%fdywtr     = 0.
          mon_a(m)%fayg       = 0. 
          mon_a(m)%faydseg    = 0. 
          mon_a(m)%faycptg    = 0. 
          mon_a(m)%faygzg     = 0. 
          mon_a(m)%fayleg     = 0. 
          mon_a(m)%faywtrg    = 0. 
          mon_a(m)%fdydseg    = 0. 
          mon_a(m)%fdyleg     = 0. 
          mon_a(m)%fdywtrg    = 0. 
          mon_a(m)%fydseg     = 0. 
          mon_a(m)%fyheatg    = 0. 
          mon_a(m)%fyleg      = 0. 
          mon_a(m)%fywtrg     = 0. 
          mon_a(m)%fswr_sur    = 0. 
          mon_a(m)%flwr_dw_sur = 0. 
          mon_a(m)%flwr_dw_sur_cs = 0. 
          mon_a(m)%flwr_dw_sur_cld = 0. 
          mon_a(m)%flwr_up_sur = 0. 
          mon_a(m)%dswd_dalb_vu_cs= 0.
          mon_a(m)%dswd_dalb_ir_cs= 0.
          mon_a(m)%dswd_dalb_vu_cld= 0.
          mon_a(m)%dswd_dalb_ir_cld= 0.
          mon_a(m)%dswd_dz_ir_cs= 0.
          mon_a(m)%dswd_dz_ir_cld= 0.
          mon_a(m)%rb_top      = 0. 
          mon_a(m)%rb_sur      = 0. 
          mon_a(m)%rb_atm      = 0. 
          mon_a(m)%swr_atm    = 0. 
          mon_a(m)%swr_top    = 0. 
          mon_a(m)%swr_top_cs    = 0. 
          mon_a(m)%swr_sur    = 0. 
          mon_a(m)%swr_sur_cs    = 0. 
          mon_a(m)%swr_dw_sur = 0. 
          mon_a(m)%swr_dw_top = 0. 
          mon_a(m)%lwr_atm    = 0. 
          mon_a(m)%lwr_atm_cld    = 0. 
          mon_a(m)%lwr_top    = 0. 
          mon_a(m)%lwr_top_cs    = 0. 
          mon_a(m)%lwr_sur    = 0. 
          mon_a(m)%lwr_sur_cs = 0. 
          mon_a(m)%lwr_tro = 0. 
          mon_a(m)%lwr_cld = 0. 
          mon_a(m)%rb_str      = 0. 
          mon_a(m)%cre_top     = 0. 
          mon_a(m)%swr_cre_top = 0. 
          mon_a(m)%lwr_cre_top = 0. 
          mon_a(m)%cre_sur     = 0. 
          mon_a(m)%swr_cre_sur = 0. 
          mon_a(m)%lwr_cre_sur = 0. 
          mon_a(m)%aplan       = 0. 
          mon_a(m)%tsl         = 0. 
          mon_a(m)%eke         = 0. 
          mon_a(m)%sam         = 0. 
          mon_a(m)%sam2        = 0. 
          mon_a(m)%synprod     = 0. 
          mon_a(m)%syndiss     = 0. 
          mon_a(m)%synadv      = 0. 
          mon_a(m)%syndif      = 0. 
          mon_a(m)%synsur      = 0. 
          mon_a(m)%cdif        = 0. 
          mon_a(m)%xz          = 0. 
          mon_a(m)%fw_pac_atl  = 0. 
          mon_a(m)%fw_atl_indpac  = 0. 
          mon_a(m)%hdust       = 0. 
          mon_a(m)%dust_load   = 0. 
          mon_a(m)%dust_emis   = 0. 
          mon_a(m)%dust_dep    = 0. 
          mon_a(m)%dust_dep_dry= 0. 
          mon_a(m)%dust_dep_wet= 0. 
          mon_a(m)%dust_ot     = 0. 
          mon_a(m)%cam         = 0. 
          mon_a(m)%co2d        = 0. 
          mon_a(m)%co2flx      = 0. 
          mon_a(m)%so4_ot      = 0. 
        enddo
      endif

      mon_a(mon)%had_fi      = mon_a(mon)%had_fi      + mon_avg * atm%had_fi*180./pi    ! deg
      mon_a(mon)%had_width   = mon_a(mon)%had_width   + mon_avg * atm%had_width*180./pi ! deg
      mon_a(mon)%ptrop       = mon_a(mon)%ptrop       + mon_avg * atm%ptrop
      mon_a(mon)%tslz        = mon_a(mon)%tslz        + mon_avg * sum(atm%tsl,1)/im
      mon_a(mon)%tskslz       = mon_a(mon)%tskslz       + mon_avg * sum(atm%tsksl,1)/im
      mon_a(mon)%ekez        = mon_a(mon)%ekez        + mon_avg * sum(atm%sam,1)/im
      mon_a(mon)%slpz        = mon_a(mon)%slpz        + mon_avg * sum(atm%slp,1)/im
      mon_a(mon)%vabz        = mon_a(mon)%vabz        + mon_avg * sum(atm%vab,1)/im
      mon_a(mon)%acbarz      = mon_a(mon)%acbarz      + mon_avg * sum(atm%acbar,1)/im
      mon_a(mon)%uz850       = mon_a(mon)%uz850       + mon_avg * sum(atm%u3(:,:,k850),1)/im
      mon_a(mon)%uz500       = mon_a(mon)%uz500       + mon_avg * sum(atm%u3(:,:,k500),1)/im
      mon_a(mon)%tam         = mon_a(mon)%tam         + mon_avg * atm%tam
      mon_a(mon)%qam         = mon_a(mon)%qam         + mon_avg * atm%qam
      mon_a(mon)%ram         = mon_a(mon)%ram         + mon_avg * atm%ram
      mon_a(mon)%gams        = mon_a(mon)%gams        + mon_avg * atm%gams
      mon_a(mon)%gamb        = mon_a(mon)%gamb        + mon_avg * atm%gamb
      mon_a(mon)%gamt        = mon_a(mon)%gamt        + mon_avg * atm%gamt
      mon_a(mon)%dam         = mon_a(mon)%dam         + mon_avg * atm%dam
      mon_a(mon)%hqeff       = mon_a(mon)%hqeff       + mon_avg * atm%hqeff
      mon_a(mon)%hrm         = mon_a(mon)%hrm         + mon_avg * atm%hrm
      mon_a(mon)%wcon        = mon_a(mon)%wcon        + mon_avg * atm%wcon
      mon_a(mon)%cld_rh         = mon_a(mon)%cld_rh         + mon_avg * atm%cld_rh
      mon_a(mon)%cld_low        = mon_a(mon)%cld_low        + mon_avg * atm%cld_low
      mon_a(mon)%cld         = mon_a(mon)%cld         + mon_avg * atm%cld
      mon_a(mon)%prc         = mon_a(mon)%prc         + mon_avg * atm%prc*sec_day
      mon_a(mon)%prcw        = mon_a(mon)%prcw        + mon_avg * sum(atm%prcw*atm%frst,3)*sec_day
      mon_a(mon)%prcs        = mon_a(mon)%prcs        + mon_avg * sum(atm%prcs*atm%frst,3)*sec_day
      mon_a(mon)%prc_conv    = mon_a(mon)%prc_conv    + mon_avg * atm%prc_conv*sec_day
      mon_a(mon)%prc_wcon    = mon_a(mon)%prc_wcon    + mon_avg * atm%prc_wcon*sec_day
      mon_a(mon)%prc_over    = mon_a(mon)%prc_over    + mon_avg * atm%prc_over*sec_day
      mon_a(mon)%hcld        = mon_a(mon)%hcld        + mon_avg * atm%hcld
      mon_a(mon)%clot        = mon_a(mon)%clot        + mon_avg * atm%clot
      mon_a(mon)%alb_cld     = mon_a(mon)%alb_cld     + mon_avg * atm%alb_cld
      mon_a(mon)%htrop       = mon_a(mon)%htrop       + mon_avg * atm%htrop
      mon_a(mon)%ttrop       = mon_a(mon)%ttrop       + mon_avg * atm%ttrop
      mon_a(mon)%wind        = mon_a(mon)%wind        + mon_avg * atm%winda
      mon_a(mon)%frst        = mon_a(mon)%frst        + mon_avg * atm%frst
      mon_a(mon)%tskin       = mon_a(mon)%tskin       + mon_avg * atm%tskin
      mon_a(mon)%t2          = mon_a(mon)%t2          + mon_avg * atm%t2
      mon_a(mon)%q2          = mon_a(mon)%q2          + mon_avg * atm%q2
      mon_a(mon)%r2          = mon_a(mon)%r2          + mon_avg * atm%r2
      mon_a(mon)%ra2         = mon_a(mon)%ra2         + mon_avg * atm%ra2
      mon_a(mon)%alb_vu_s    = mon_a(mon)%alb_vu_s    + mon_avg * atm%alb_vu_s
      mon_a(mon)%alb_vu_c    = mon_a(mon)%alb_vu_c    + mon_avg * atm%alb_vu_c
      mon_a(mon)%alb_ir_s    = mon_a(mon)%alb_ir_s    + mon_avg * atm%alb_ir_s
      mon_a(mon)%alb_ir_c    = mon_a(mon)%alb_ir_c    + mon_avg * atm%alb_ir_c
      mon_a(mon)%alb_sur_cs  = mon_a(mon)%alb_sur_cs  + mon_avg * (frac_vu*sum(atm%alb_vu_s*atm%frst,3) + (1.-frac_vu)*sum(atm%alb_ir_s*atm%frst,3))
      mon_a(mon)%alb_sur_cld = mon_a(mon)%alb_sur_cld + mon_avg * (frac_vu*sum(atm%alb_vu_c*atm%frst,3) + (1.-frac_vu)*sum(atm%alb_ir_c*atm%frst,3))
      mon_a(mon)%tskina      = mon_a(mon)%tskina      + mon_avg * sum(atm%tskin*atm%frst,3)
      mon_a(mon)%t2a         = mon_a(mon)%t2a         + mon_avg * sum(atm%t2*atm%frst,3)
      do i=1,im
        do j=1,jm
          mon_a(mon)%thetae(i,j) = mon_a(mon)%thetae(i,j) + mon_avg * (atm%t2a(i,j)*exp(cle*atm%ram(i,j)*fqsat(atm%t2a(i,j),101300._wp)/(cp*atm%t2a(i,j))))
        enddo
      enddo
      mon_a(mon)%q2a         = mon_a(mon)%q2a         + mon_avg * sum(atm%q2*atm%frst,3)
      mon_a(mon)%r2a         = mon_a(mon)%r2a         + mon_avg * sum(atm%r2*atm%frst,3)
      mon_a(mon)%rskina      = mon_a(mon)%rskina      + mon_avg * atm%rskina
      mon_a(mon)%ra2a        = mon_a(mon)%ra2a        + mon_avg * atm%ra2a
      mon_a(mon)%cd0a        = mon_a(mon)%cd0a        + mon_avg * atm%cd0a
      mon_a(mon)%cda         = mon_a(mon)%cda         + mon_avg * atm%cda
      mon_a(mon)%cd          = mon_a(mon)%cd          + mon_avg * atm%cd 
      mon_a(mon)%sha         = mon_a(mon)%sha         + mon_avg * atm%sha
      mon_a(mon)%lha         = mon_a(mon)%lha         + mon_avg * atm%lha
      mon_a(mon)%evpa        = mon_a(mon)%evpa        + mon_avg * atm%evpa*sec_day
      mon_a(mon)%Ri          = mon_a(mon)%Ri          + mon_avg * (g*z_sfl*(1._wp-atm%tskina/atm%tam)/(sqrt(atm%ugb**2+atm%vgb**2+(c_syn_6*sqrt(atm%sam))**2)))
      mon_a(mon)%acbar       = mon_a(mon)%acbar       + mon_avg * atm%acbar
      mon_a(mon)%epsa        = mon_a(mon)%epsa        + mon_avg * sum(atm%epsa*atm%frst,3)
      mon_a(mon)%tsksl       = mon_a(mon)%tsksl       + mon_avg * atm%tsksl
      do j=1,jm
        mon_a(mon)%atsksl(:,j)        = mon_a(mon)%atsksl(:,j)        + mon_avg * (atm%tsksl(:,j)-sum(atm%tsksl(:,j))*aim)
        mon_a(mon)%atsl(:,j)        = mon_a(mon)%atsl(:,j)        + mon_avg * (atm%tsl(:,j)-sum(atm%tsl(:,j))*aim)
      enddo
      do i=1,im
      do j=1,jm
        mon_a(mon)%atsli(i,j)        = mon_a(mon)%atsli(i,j)        + mon_avg * (sum(atm%t3(i,j,1:8))-sum(sum(atm%t3(:,j,1:8),2),1)*aim)
      enddo
      enddo
      mon_a(mon)%aslp        = mon_a(mon)%aslp        + mon_avg * atm%aslp
      mon_a(mon)%aslp_temp   = mon_a(mon)%aslp_temp   + mon_avg * atm%aslp_temp
      mon_a(mon)%aslp_topo   = mon_a(mon)%aslp_topo   + mon_avg * atm%aslp_topo
      mon_a(mon)%dz500       = mon_a(mon)%dz500       + mon_avg * atm%dz500 
      mon_a(mon)%slp         = mon_a(mon)%slp         + mon_avg * atm%slp
      mon_a(mon)%us         = mon_a(mon)%us         + mon_avg * sum(atm%us*atm%frst,3)
      mon_a(mon)%vs         = mon_a(mon)%vs         + mon_avg * sum(atm%vs*atm%frst,3)
      mon_a(mon)%usk        = mon_a(mon)%usk        + mon_avg * atm%usk
      mon_a(mon)%vsk        = mon_a(mon)%vsk        + mon_avg * atm%vsk
      mon_a(mon)%ugb        = mon_a(mon)%ugb        + mon_avg * atm%ugb
      mon_a(mon)%vgb        = mon_a(mon)%vgb        + mon_avg * atm%vgb
      mon_a(mon)%uab        = mon_a(mon)%uab        + mon_avg * atm%uab
      mon_a(mon)%vab        = mon_a(mon)%vab        + mon_avg * atm%vab
      mon_a(mon)%taux        = mon_a(mon)%taux        + mon_avg * atm%taux
      mon_a(mon)%tauy        = mon_a(mon)%tauy        + mon_avg * atm%tauy
      mon_a(mon)%wcld        = mon_a(mon)%wcld        + mon_avg * atm%wcld
      mon_a(mon)%woro        = mon_a(mon)%woro        + mon_avg * atm%woro
      mon_a(mon)%weff        = mon_a(mon)%weff        + mon_avg * atm%weff
      mon_a(mon)%fweff       = mon_a(mon)%fweff       + mon_avg * atm%fweff
      mon_a(mon)%ua         = mon_a(mon)%ua         + mon_avg * atm%ua
      mon_a(mon)%va         = mon_a(mon)%va         + mon_avg * atm%va
      mon_a(mon)%u3         = mon_a(mon)%u3         + mon_avg * atm%u3
      mon_a(mon)%v3         = mon_a(mon)%v3         + mon_avg * atm%v3
      mon_a(mon)%w3         = mon_a(mon)%w3         + mon_avg * atm%w3
      mon_a(mon)%uter        = mon_a(mon)%uter        + mon_avg * atm%uter
      mon_a(mon)%vter        = mon_a(mon)%vter        + mon_avg * atm%vter
      if (l_output_flx3d) then
        mon_a(mon)%fax        = mon_a(mon)%fax        + mon_avg * atm%fax
        mon_a(mon)%faxo       = mon_a(mon)%faxo       + mon_avg * atm%faxo
        mon_a(mon)%fay        = mon_a(mon)%fay        + mon_avg * atm%fay
        mon_a(mon)%fayo       = mon_a(mon)%fayo       + mon_avg * atm%fayo
        mon_a(mon)%fac        = mon_a(mon)%fac        + mon_avg * atm%fac
      endif
      mon_a(mon)%diffxdse    = mon_a(mon)%diffxdse        + mon_avg * atm%diffxdse
      mon_a(mon)%diffydse    = mon_a(mon)%diffydse        + mon_avg * atm%diffydse
      mon_a(mon)%diffxwtr    = mon_a(mon)%diffxwtr        + mon_avg * atm%diffxwtr
      mon_a(mon)%diffywtr    = mon_a(mon)%diffywtr        + mon_avg * atm%diffywtr
      mon_a(mon)%wsyn        = mon_a(mon)%wsyn        + mon_avg * atm%wsyn
      mon_a(mon)%convdse     = mon_a(mon)%convdse     + mon_avg * atm%convdse
      mon_a(mon)%convwtr     = mon_a(mon)%convwtr     + mon_avg * atm%convwtr*sec_day
      mon_a(mon)%hdust        = mon_a(mon)%hdust        + mon_avg * atm%hdust
      mon_a(mon)%dust_load    = mon_a(mon)%dust_load    + mon_avg * atm%dust_load   
      mon_a(mon)%dust_emis    = mon_a(mon)%dust_emis    + mon_avg * atm%dust_emis   
      mon_a(mon)%dust_dep     = mon_a(mon)%dust_dep     + mon_avg * atm%dust_dep    
      mon_a(mon)%dust_dep_dry = mon_a(mon)%dust_dep_dry + mon_avg * atm%dust_dep_dry
      mon_a(mon)%dust_dep_wet = mon_a(mon)%dust_dep_wet + mon_avg * atm%dust_dep_wet
      mon_a(mon)%dust_ot      = mon_a(mon)%dust_ot      + mon_avg * atm%dust_ot     
      mon_a(mon)%cam          = mon_a(mon)%cam  + mon_avg * atm%cam    ! kgCO2/kg 
      mon_a(mon)%co2d         = mon_a(mon)%co2d + mon_avg * atm%cam*1.e6*28.97/44.0095    ! kgCO2/kg -> ppm
      mon_a(mon)%co2flx       = mon_a(mon)%co2flx + mon_avg * atm%co2flx*1.e3*sec_day   ! gCO2/m2/day
      mon_a(mon)%so4_ot       = mon_a(mon)%so4_ot       + mon_avg * sigma_so4*atm%so4    

      do j=1,jm
        do i=1,im

          ipl=i+1
          if (ipl.gt.im) ipl=1
          imi=i-1
          if (imi.lt.1) imi=im  
          jmi=max(1,j-1)
          jpl=min(jm,j+1)

          ! derive sea level pressure to be compared to other models, using constant lapse rate of 6.5 K/km
          ! compute surface pressure from sea level pressure using actual lapse rate
          ps = atm%slp(i,j)*(1.+atm%gamb(i,j)*atm%zsa(i,j)/(atm%tam(i,j)-atm%gamb(i,j)*atm%zsa(i,j)))**(-g/(Rd*atm%gamb(i,j)))
          ! reduce pressure to sea level using constant lapse rate of 6.5 K/km
          ! https://www.wmo.int/pages/prog/www/IMOP/meetings/SI/ET-Stand-1/Doc-10_Pressure-red.pdf
          slp = ps*(1.+6.5e-3*atm%zsa(i,j)/(atm%tam(i,j)-6.5e-3*atm%zsa(i,j)))**(g/(Rd*6.5e-3))
          mon_a(mon)%slp1(i,j) = mon_a(mon)%slp1(i,j) + mon_avg * slp 

          ! surface humidity gradient
          mon_a(mon)%dq(i,j) = mon_a(mon)%dq(i,j) + mon_avg * (fqsat(atm%tskina(i,j),atm%psa(i,j))-atm%q2a(i,j))
          mon_a(mon)%dr(i,j) = mon_a(mon)%dr(i,j) + mon_avg * (atm%rskina(i,j)-atm%ram(i,j))
        enddo
      enddo

      do j=1,jm
        do i=1,im
          ! heat+geopotential energy convergence from advection only
          mon_a(mon)%convadse(i,j) = mon_a(mon)%convadse(i,j) &
            + (atm%faxdse(i,j)-atm%faxdse(i+1,j) + atm%faydse(i,j+1)-atm%faydse(i,j))/sqr(i,j)*cp*mon_avg  ! K * kg/s / m2 * J/kg/K = W/m2
          ! heat + geopotential energy convergence from diffusion only
          mon_a(mon)%convddse(i,j) = mon_a(mon)%convddse(i,j) &
            + (atm%fdxdse(i,j)-atm%fdxdse(i+1,j) + atm%fdydse(i,j+1)-atm%fdydse(i,j))/sqr(i,j)*cp*mon_avg
          ! moisture convergence from advection only
          mon_a(mon)%convawtr(i,j) = mon_a(mon)%convawtr(i,j) &
            + (atm%faxwtr(i,j)-atm%faxwtr(i+1,j) + atm%faywtr(i,j+1)-atm%faywtr(i,j))/sqr(i,j)*sec_day*mon_avg  ! kg/s / m2 *s/day = kg/m2/day
          ! moisture convergence from diffusion only
          mon_a(mon)%convdwtr(i,j) = mon_a(mon)%convdwtr(i,j) &
            + (atm%fdxwtr(i,j)-atm%fdxwtr(i+1,j) + atm%fdywtr(i,j+1)-atm%fdywtr(i,j))/sqr(i,j)*sec_day*mon_avg
        enddo
      enddo

      ! compute cp*T fluxes for output
      allocate(faxcpt(imc,jm))
      allocate(faycpt(im,jmc))
      do j=1,jm
        jmi=max(1,j-1)
        do i=1,im
          imi=i-1
          if (imi.eq.0) imi=im
          faxcpt(i,j)=0._wp      
          faycpt(i,j)=0._wp      
          do k=1,km
            if (atm%fax(i,j,k).gt.0.) then
              tup=atm%t3(imi,j,k)
            else    
              tup=atm%t3(i,j,k)
            endif
            faxcpt(i,j)=faxcpt(i,j)+atm%fax(i,j,k)*tup ! kg/s * K
            if (j.gt.1) then
              if (atm%fay(i,j,k).gt.0._wp) then
                tup=atm%t3(i,j,k)
              else    
                tup=atm%t3(i,jmi,k)
              endif
              faycpt(i,j)=faycpt(i,j)+atm%fay(i,j,k)*tup    ! kg/s * K
            endif
          enddo
          if (j.eq.jm) then
            faycpt(i,jmc)=0.              
          endif
        enddo
        faxcpt(imc,j)=faxcpt(1,j)              
      enddo

      do j=1,jm
        do i=1,im
          ! vertically integrated mass flux
          faxmasi = 0.
          faymasi = 0.
          do k=1,km-1
            faxmasi = faxmasi + 0.5*(atm%fax(i,j,k) + atm%fax(i+1,j,k))  ! kg/m2/s
            faymasi = faymasi + 0.5*(atm%fay(i,j,k) + atm%fay(i,j+1,k))  ! kg/m2/s
          enddo
          mon_a(mon)%faxmas(i,j)     = mon_a(mon)%faxmas(i,j)     + mon_avg * faxmasi/dy      ! kg/m2/s /m  = kg/s/m
          mon_a(mon)%faymas(i,j)     = mon_a(mon)%faymas(i,j)     + mon_avg * faymasi/dy      ! kg/m2/s /m  = kg/s/m
          ! vertically integrated potential temperature and water vapor flux
          mon_a(mon)%faxdse(i,j)     = mon_a(mon)%faxdse(i,j)     + mon_avg * 0.5*(atm%faxdse(i,j) + atm%faxdse(i+1,j))*cp/dy      ! kg/s*K *J/K/kg /m  = W/m
          mon_a(mon)%faxcpt(i,j)     = mon_a(mon)%faxcpt(i,j)     + mon_avg * 0.5*(faxcpt(i,j) + faxcpt(i+1,j))*cp/dy      ! kg/s*K *J/K/kg /m  = W/m
          mon_a(mon)%faxwtr(i,j)     = mon_a(mon)%faxwtr(i,j)     + mon_avg * 0.5*(atm%faxwtr(i,j) + atm%faxwtr(i+1,j))/dy 
          mon_a(mon)%faydse(i,j)     = mon_a(mon)%faydse(i,j)     + mon_avg * 0.5*(atm%faydse(i,j) + atm%faydse(i,j+1))*cp/dxt(j) 
          mon_a(mon)%faycpt(i,j)     = mon_a(mon)%faycpt(i,j)     + mon_avg * 0.5*(faycpt(i,j) + faycpt(i,j+1))*cp/dxt(j) 
          mon_a(mon)%faywtr(i,j)     = mon_a(mon)%faywtr(i,j)     + mon_avg * 0.5*(atm%faywtr(i,j) + atm%faywtr(i,j+1))/dxt(j) 
          mon_a(mon)%fdxdse(i,j)     = mon_a(mon)%fdxdse(i,j)     + mon_avg * 0.5*(atm%fdxdse(i,j) + atm%fdxdse(i+1,j))*cp/dy 
          mon_a(mon)%fdxwtr(i,j)     = mon_a(mon)%fdxwtr(i,j)     + mon_avg * 0.5*(atm%fdxwtr(i,j) + atm%fdxwtr(i+1,j))/dy 
          mon_a(mon)%fdydse(i,j)     = mon_a(mon)%fdydse(i,j)     + mon_avg * 0.5*(atm%fdydse(i,j) + atm%fdydse(i,j+1))*cp/dxt(j) 
          mon_a(mon)%fdywtr(i,j)     = mon_a(mon)%fdywtr(i,j)     + mon_avg * 0.5*(atm%fdywtr(i,j) + atm%fdywtr(i,j+1))/dxt(j) 
        enddo
      enddo

      fayg(:)=0
      faydseg(:)=0
      faycptg(:)=0
      faygzg(:) =0
      fayleg(:)=0
      faywtrg(:)=0
      fdydseg(:)=0
      fdyleg(:)=0
      fdywtrg(:)=0
      fydseg(:)=0
      fyheatg(:)=0
      fyleg(:)=0
      fywtrg(:)=0

      do j=2,jm
        do i=1,im

          do k=1,km
            fayg(j) =fayg(j)+atm%fay(i,j,k)
          enddo 

          faydseg(j)=faydseg(j)+atm%faydse(i,j)*cp  ! kg/s*K *J/K/kg = W
          faycptg(j)=faycptg(j)+faycpt(i,j)*cp      ! kg/s*K *J/K/kg = W
          fayleg(j) =fayleg(j) +atm%faywtr(i,j)*cle ! kg/s * J/kg = W
          faywtrg(j)=faywtrg(j)+atm%faywtr(i,j)     ! kg/s
          fdydseg(j)=fdydseg(j)+atm%fdydse(i,j)*cp
          fdyleg(j) =fdyleg(j) +atm%fdywtr(i,j)*cle ! kg/s * J/kg = W
          fdywtrg(j)=fdywtrg(j)+atm%fdywtr(i,j)

        enddo

        faygzg(j)=faydseg(j)-faycptg(j)
        fydseg(j)=faydseg(j)+fdydseg(j)
        fyleg(j) =fayleg(j)+fdyleg(j)
        fyheatg(j)=fydseg(j)+fyleg(j) 
        fywtrg(j)=faywtrg(j)+fdywtrg(j)

      enddo

      deallocate(faxcpt)
      deallocate(faycpt)

      mon_a(mon)%fayg       = mon_a(mon)%fayg       + mon_avg * fayg
      mon_a(mon)%faydseg    = mon_a(mon)%faydseg    + mon_avg * faydseg*1.e-15    ! PW
      mon_a(mon)%faycptg    = mon_a(mon)%faycptg    + mon_avg * faycptg*1.e-15    ! PW
      mon_a(mon)%faygzg     = mon_a(mon)%faygzg     + mon_avg * faygzg *1.e-15    ! PW
      mon_a(mon)%fayleg     = mon_a(mon)%fayleg     + mon_avg * fayleg *1.e-15    ! PW
      mon_a(mon)%faywtrg    = mon_a(mon)%faywtrg    + mon_avg * faywtrg
      mon_a(mon)%fdydseg    = mon_a(mon)%fdydseg    + mon_avg * fdydseg*1.e-15    ! PW
      mon_a(mon)%fdyleg     = mon_a(mon)%fdyleg     + mon_avg * fdyleg*1.e-15     ! PW
      mon_a(mon)%fdywtrg    = mon_a(mon)%fdywtrg    + mon_avg * fdywtrg
      mon_a(mon)%fydseg     = mon_a(mon)%fydseg     + mon_avg * fydseg*1.e-15    ! PW
      mon_a(mon)%fyheatg    = mon_a(mon)%fyheatg    + mon_avg * fyheatg*1.e-15   ! PW
      mon_a(mon)%fyleg      = mon_a(mon)%fyleg      + mon_avg * fyleg*1.e-15      ! PW
      mon_a(mon)%fywtrg     = mon_a(mon)%fywtrg     + mon_avg * fywtrg

      ! water vapor transport in and out of the Atlantic catchment 
      fw_pac_atl = 0.
      fw_atl_indpac = 0.
      do j=1,jm
        if (atm%idivide_pac_atl(j).gt.0) then
          ind = atm%idivide_pac_atl(j)+1
          if (ind.eq.im+1) ind=1
          fw_pac_atl(j) = (atm%faxwtr(ind,j)+atm%fdxwtr(ind,j))*1.e-9 ! kg/s -> Sv
        endif
        if (atm%idivide_atl_indpac(j).gt.0) then
          ind = atm%idivide_atl_indpac(j)+1
          if (ind.eq.im+1) ind=1
          fw_atl_indpac(j) = (atm%faxwtr(ind,j)+atm%fdxwtr(ind,j))*1.e-9 ! kg/s -> Sv
        endif
      enddo

      mon_a(mon)%fw_pac_atl    = mon_a(mon)%fw_pac_atl    + mon_avg * fw_pac_atl
      mon_a(mon)%fw_atl_indpac = mon_a(mon)%fw_atl_indpac + mon_avg * fw_atl_indpac

      mon_a(mon)%fswr_sur    = mon_a(mon)%fswr_sur    + mon_avg * atm%fswr_sur
      mon_a(mon)%flwr_dw_sur = mon_a(mon)%flwr_dw_sur + mon_avg * atm%flwr_dw_sur
      mon_a(mon)%flwr_dw_sur_cs  = mon_a(mon)%flwr_dw_sur_cs  + mon_avg * atm%flwr_dw_sur_cs
      mon_a(mon)%flwr_dw_sur_cld = mon_a(mon)%flwr_dw_sur_cld + mon_avg * atm%flwr_dw_sur_cld
      mon_a(mon)%flwr_up_sur = mon_a(mon)%flwr_up_sur + mon_avg * atm%flwr_up_sur
      mon_a(mon)%dswd_dalb_vu_cs  = mon_a(mon)%dswd_dalb_vu_cs  + mon_avg * atm%dswd_dalb_vu_cs   
      mon_a(mon)%dswd_dalb_ir_cs  = mon_a(mon)%dswd_dalb_ir_cs  + mon_avg * atm%dswd_dalb_ir_cs  
      mon_a(mon)%dswd_dalb_vu_cld = mon_a(mon)%dswd_dalb_vu_cld + mon_avg * atm%dswd_dalb_vu_cld 
      mon_a(mon)%dswd_dalb_ir_cld = mon_a(mon)%dswd_dalb_ir_cld + mon_avg * atm%dswd_dalb_ir_cld 
      mon_a(mon)%dswd_dz_ir_cs  = mon_a(mon)%dswd_dz_ir_cs  + mon_avg * atm%dswd_dz_ir_cs 
      mon_a(mon)%dswd_dz_ir_cld = mon_a(mon)%dswd_dz_ir_cld + mon_avg * atm%dswd_dz_ir_cld 
      mon_a(mon)%swr_atm    = mon_a(mon)%swr_atm    + mon_avg * (atm%swr_top-atm%swr_sur)
      mon_a(mon)%swr_top    = mon_a(mon)%swr_top    + mon_avg * atm%swr_top
      mon_a(mon)%swr_top_cs    = mon_a(mon)%swr_top_cs    + mon_avg * atm%swr_top_cs
      mon_a(mon)%swr_sur    = mon_a(mon)%swr_sur    + mon_avg * sum(atm%fswr_sur*atm%frst,3)
      mon_a(mon)%swr_sur_cs    = mon_a(mon)%swr_sur_cs    + mon_avg * sum(atm%fswr_sur_cs*atm%frst,3)
      mon_a(mon)%lwr_atm    = mon_a(mon)%lwr_atm    + mon_avg * (atm%lwr_top-atm%lwr_sur)
      mon_a(mon)%lwr_atm_cld    = mon_a(mon)%lwr_atm_cld    + mon_avg * (atm%lwr_cld-atm%lwr_sur)
      mon_a(mon)%lwr_top    = mon_a(mon)%lwr_top    + mon_avg * atm%lwr_top
      mon_a(mon)%lwr_top_cs    = mon_a(mon)%lwr_top_cs    + mon_avg * atm%lwr_top_cs
      mon_a(mon)%lwr_sur    = mon_a(mon)%lwr_sur    + mon_avg * atm%lwr_sur
      mon_a(mon)%lwr_sur_cs = mon_a(mon)%lwr_sur_cs + mon_avg * sum((atm%flwr_dw_sur_cs-atm%flwr_up_sur)*atm%frst,3)
      mon_a(mon)%lwr_tro = mon_a(mon)%lwr_tro + mon_avg * atm%lwr_tro
      mon_a(mon)%lwr_cld = mon_a(mon)%lwr_cld + mon_avg * atm%lwr_cld
      mon_a(mon)%rb_top      = mon_a(mon)%rb_top      + mon_avg * atm%rb_top
      mon_a(mon)%rb_sur      = mon_a(mon)%rb_sur      + mon_avg * atm%rb_sur
      mon_a(mon)%rb_atm      = mon_a(mon)%rb_atm      + mon_avg * atm%rb_atm
      mon_a(mon)%rb_str      = mon_a(mon)%rb_str      + mon_avg * atm%rb_str
      mon_a(mon)%cre_top     = mon_a(mon)%cre_top     + mon_avg * (atm%swr_top_cld+atm%lwr_top_cld-atm%swr_top_cs-atm%lwr_top_cs)
      mon_a(mon)%swr_cre_top = mon_a(mon)%swr_cre_top + mon_avg * (atm%swr_top_cld-atm%swr_top_cs)
      mon_a(mon)%lwr_cre_top = mon_a(mon)%lwr_cre_top + mon_avg * (atm%lwr_top_cld-atm%lwr_top_cs)
      mon_a(mon)%cre_sur     = mon_a(mon)%cre_sur     + mon_avg * (sum((atm%fswr_sur_cld-atm%fswr_sur_cs)*atm%frst,3)   &
        +sum(((atm%flwr_dw_sur_cld-atm%flwr_up_sur)-(atm%flwr_dw_sur_cs-atm%flwr_up_sur))*atm%frst,3))
      mon_a(mon)%swr_cre_sur = mon_a(mon)%swr_cre_sur + mon_avg * sum((atm%fswr_sur_cld-atm%fswr_sur_cs)*atm%frst,3)
      mon_a(mon)%lwr_cre_sur = mon_a(mon)%lwr_cre_sur + mon_avg * sum(((atm%flwr_dw_sur_cld-atm%flwr_up_sur)-(atm%flwr_dw_sur_cs-atm%flwr_up_sur))*atm%frst,3) 
      mon_a(mon)%tsl         = mon_a(mon)%tsl         + mon_avg * atm%tsl
      mon_a(mon)%eke         = mon_a(mon)%eke         + mon_avg * atm%eke
      mon_a(mon)%sam         = mon_a(mon)%sam         + mon_avg * atm%sam
      mon_a(mon)%sam2        = mon_a(mon)%sam2        + mon_avg * atm%sam2
      mon_a(mon)%synprod     = mon_a(mon)%synprod     + mon_avg * atm%synprod
      mon_a(mon)%syndiss     = mon_a(mon)%syndiss     + mon_avg * atm%syndiss
      mon_a(mon)%synadv      = mon_a(mon)%synadv      + mon_avg * atm%synadv            
      mon_a(mon)%syndif      = mon_a(mon)%syndif      + mon_avg * atm%syndif
      mon_a(mon)%synsur      = mon_a(mon)%synsur      + mon_avg * sum(atm%synsur*atm%frst,3)
      mon_a(mon)%cdif        = mon_a(mon)%cdif        + mon_avg * atm%cdif           

      do j=1,jm
        do i=1,im

          imi=i-1
          if (imi.lt.1) imi=im
          ipl=i+1
          if (ipl.gt.im) ipl=1
          jmi=j-1
          if (jmi.lt.1) jmi=1
          jpl=j+1
          if (jpl.gt.jm) jpl=jm        

          tsl(i,j) = atm%tam(i,j) - atm%gams(i,j)*1500._wp + atm%gamb(i,j)*atm%zsa(i,j)

        enddo
      enddo

      do j=1,jm
        do i=1,im

          ! incoming radiation at TOA
          mon_a(mon)%swr_dw_top(i,j) = mon_a(mon)%swr_dw_top(i,j) + mon_avg * atm%swr_dw_top(i,j)

          ! downward shortwave radiation at the surface
          mon_a(mon)%swr_dw_sur(i,j) = mon_a(mon)%swr_dw_sur(i,j) + mon_avg * &
            (atm%cld(i,j) * (frac_vu*atm%swr_dw_sur_vis_cld(i,j) + (1._wp-frac_vu)*atm%swr_dw_sur_nir_cld(i,j))  &
            + (1._wp-atm%cld(i,j)) * (frac_vu*atm%swr_dw_sur_vis_cs(i,j) + (1._wp-frac_vu)*atm%swr_dw_sur_nir_cs(i,j)))

          ! planetary albedo
          if (atm%solarm(doy,j).gt.1._wp) then
            mon_a(mon)%aplan(i,j) = mon_a(mon)%aplan(i,j) + mon_avg * (1._wp-atm%swr_top(i,j)/atm%solarm(doy,j))
          endif

          ! top cloud temperature
          ctt = t_prof(atm%zsa(i,j), atm%hcld(i,j), atm%tam(i,j), atm%gams(i,j), atm%gamb(i,j), atm%gamt(i,j), atm%htrop(i,j), 1)
          mon_a(mon)%ctt(i,j) = mon_a(mon)%ctt(i,j) + mon_avg * ctt

          ! recalculation of 3d fields to levels
          do k=1,kmc

            ! temperature
            t = t_prof(atm%zsa(i,j), zl(k), atm%tam(i,j), atm%gams(i,j), atm%gamb(i,j), atm%gamt(i,j), atm%htrop(i,j), 1)
            mon_a(mon)%t3(i,j,k) = mon_a(mon)%t3(i,j,k) + mon_avg * t

            ! potential temperature
            mon_a(mon)%tp(i,j,k) = mon_a(mon)%tp(i,j,k) + mon_avg * (t+gad*zl(k))

            ! relative humidity
            rh = rh_prof(atm%zsa(i,j), zl(k), atm%ram(i,j), atm%hrm(i,j), atm%htrop(i,j))
            mon_a(mon)%r3(i,j,k) = mon_a(mon)%r3(i,j,k) + mon_avg * rh 

            ! specific humidity 
            mon_a(mon)%q3(i,j,k) = mon_a(mon)%q3(i,j,k) + mon_avg * rh*fqsat(t,p0*exp(-zl(k)/hatm))
          enddo

          if (time_eom_atm) then
            do k=1,km
              ! lapse rate
              mon_a(mon)%gamma(i,j,k) = - (mon_a(mon)%t3(i,j,k+1)-mon_a(mon)%t3(i,j,k))/(zl(k+1)-zl(k))
            enddo
            do k=1,kmc
              ! density, eq. 4.70 in PIK-Report
              mon_a(mon)%rho(i,j,k) = ra*exp(-zl(k)/hatm) &
                * (1.-1./T0*(mon_a(mon)%t3(i,j,1)-T0)) &
                * (1.-1./mon_a(mon)%t3(i,j,1)*(mon_a(mon)%t3(i,j,k)-mon_a(mon)%t3(i,j,1)))
            enddo
          endif

        enddo
      enddo 

      ! atmospheric mean meridional circulation 

      allocate(xz(jmc,kmc))
      xz(:,:) = 0.
      do j=2,jm
        do k=1,km
          xsum=0.
          do i=1,im
            xsum=xsum+atm%fay(i,j,k)*1.E-10
          enddo
          xz(j,k+1)=xz(j,k)+xsum       
        enddo
      enddo       
      mon_a(mon)%xz = mon_a(mon)%xz + xz * mon_avg
      deallocate(xz)

      ! rate of change of atmospheric temperature
      mon_a(mon)%dtamdt = mon_a(mon)%dtamdt &   ! K/day
        + (atm%convdse + atm%rb_atm + atm%sha + cle*sum(atm%prcw*atm%frst,3)+cls*sum(atm%prcs*atm%frst,3))/(amas*cv) * sec_day * mon_avg   

!      ! todo meridional heat transport from heat fluxes
!      do i=1,im
!        do j=1,jm
!      ! net heat at the surface
!      mon_a(mon)%qnets(i,j) = mon_a(mon)%qnets(i,j) + (atm%aflwr_sur(i,j)+atm%afswr_sur(i,j)-atm%sha(i,j)-atm%lha(i,j)) * mon_avg
!      ! net heat at TOA
!      mon_a(mon)%qnett(i,j) = mon_a(mon)%qnett(i,j) + atm%rb_top(i,j) * mon_avg
!      enddo
!      enddo
!      
!      if (nday.eq.360) then 
!      
!      do i=1,im
!       do j=1,jm
!        qtopgan = qtopgan +sqr(i,j)*tqneta(i,j)/esqr
!        qsurgan = qsurgan +sqr(i,j)*sqneta(i,j)/esqr    
!       enddo
!      enddo 
!        
!      do j=1,jm
!        sumtqnet=0.
!        sumsqnet=0.        
!       do i=1,im
!        sumtqnet=sumtqnet-(tqneta(i,j)-qtopgan)*sqr(i,j)
!        sumsqnet=sumsqnet-(sqneta(i,j)-qsurgan)*sqr(i,j)
!        enddo      
!        aoimf(j+1)=aoimf(j)+sumtqnet
!        oimf(j+1) =oimf(j) +sumsqnet       
!        oimfa(j+1)=aoimf(j+1)-amf(j+1)
!      enddo
!      
!      if (nyr.eq.nyrq.and.ksetup.eq.2) then 
!
!       write (23,rec=1) (1.e-15*aoimf(j),j=jmc,1,-1),
!     >                  (1.e-15*amf(j),j=jmc,1,-1),
!     >                  (1.e-15*oimf(j),j=jmc,1,-1)


      if (l_daily_output) then

        day_a(doy)%tam     = atm%tam
        day_a(doy)%qam     = atm%qam
        day_a(doy)%cam     = atm%cam
        day_a(doy)%tskina  = atm%tskina
        day_a(doy)%lha     = atm%lha
        day_a(doy)%sha     = atm%sha
        day_a(doy)%weff    = atm%weff
        day_a(doy)%wcld    = atm%wcld
        day_a(doy)%cld     = atm%cld
        day_a(doy)%hcld    = atm%hcld
        day_a(doy)%clot    = atm%clot
        day_a(doy)%prc     = atm%prc*sec_day
        day_a(doy)%wind    = atm%winda
        day_a(doy)%slp     = atm%slp
        day_a(doy)%us     = sum(atm%us*atm%frst,3)
        day_a(doy)%vs     = sum(atm%vs*atm%frst,3)
        day_a(doy)%ugb    = atm%ugb
        day_a(doy)%vgb    = atm%vgb
        day_a(doy)%uab    = atm%uab
        day_a(doy)%vab    = atm%vab
        day_a(doy)%ua     = atm%ua
        day_a(doy)%va     = atm%va
        day_a(doy)%u3     = atm%u3
        day_a(doy)%v3     = atm%v3
        day_a(doy)%w3     = atm%w3
        day_a(doy)%convdse = atm%convdse
        day_a(doy)%convwtr = atm%convwtr*sec_day
        day_a(doy)%sam     = atm%sam
        day_a(doy)%gams    = atm%gams
        day_a(doy)%gamb    = atm%gamb
        day_a(doy)%gamt    = atm%gamt
        day_a(doy)%hrm     = atm%hrm
        day_a(doy)%flwr_dw_sur         = atm%flwr_dw_sur
        day_a(doy)%flwr_dw_sur_cs      = atm%flwr_dw_sur_cs
        day_a(doy)%flwr_dw_sur_cld     = atm%flwr_dw_sur_cld

      endif

    endif

    if (time_out_atm .and. time_eoy_atm) then

      ! Get annual values
      call atm_ave( mon_a,ann_a )

      ann_a%solar      = atm%solar               
      ann_a%cosz       = atm%cosz
      ann_a%solarm      = atm%solarm               
      ann_a%coszm       = atm%coszm
      ann_a%zs          = atm%zs
      ann_a%zsa         = atm%zsa
      ann_a%zsa_smooth  = atm%zsa_smooth
      ann_a%slope       = atm%slope
      ann_a%slope_x     = atm%slope_X
      ann_a%slope_y     = atm%slope_y
      ann_a%frlnd       = atm%frlnd
      ann_a%frocn       = atm%frocn

      do m=1,12
        call nc_read('input/ERA5_t2m_monclim_5x5.nc',"t2m",mon_a(m)%t2a_dat(:,:),start=[1,1,m],count=[im,jm,1]) 
      enddo

      call atm_diag_out
    endif


   return

  end subroutine atm_diag


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ d i a g _ o u t
  ! Purpose  :  write atmosphere netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_diag_out

    implicit none

    integer :: k, ncid
    character (len=256) :: fnm


    nout = nout +1 

    ! write to file
    fnm = trim(out_dir)//"/atm.nc"
    call nc_open(fnm,ncid)
    call nc_write(fnm,dim_time,real(year_now,wp), dim1=dim_time, start=[nout], count=[1],ncid=ncid)    
    do k = 1, nmon_year
       call atm_nc_write(fnm,ncid,mon_a(k),k,nout)
    end do
    call atm_nc_write(fnm,ncid,ann_a,nmon_year+1,nout)
    call nc_close(ncid)

    if (l_daily_output) then
      ! write to file
      fnm = trim(out_dir)//"/atm_daily.nc"
      call nc_open(fnm,ncid)
      call nc_write(fnm,dim_time,real(year_now,wp), dim1=dim_time, start=[nout], count=[1],ncid=ncid)    
      do k = 1, nday_year
        call atm_daily_nc_write(fnm,ncid,day_a(k),k,nout)
      end do
      call nc_close(ncid)
    endif


   return

  end subroutine atm_diag_out
  

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c
  ! Purpose  :  initialize netcdf file for time series output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc(fnm)

    implicit none

    character (len=*) :: fnm
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_write_dim(fnm,dim_time, x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.)
    call nc_write_dim(fnm,dim_y, x=1, axis="y", units="1")
    call nc_write_dim(fnm,dim_x, x=1, axis="x", units="1")
!    call nc_write_dim(fnm,dim_month,x=1._wp,dx=1._wp,nx=13,units="months")

    return

  end subroutine ts_nc
 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  t s _ n c _ w r i t e
  ! Purpose  :  write time series to netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ts_nc_write(fnm,vars,nout,y)

    implicit none

    type(ts_out) :: vars(:)

    character (len=*) :: fnm
    integer :: nout, y, ncid, i


    call nc_open(fnm,ncid)
    call nc_write(fnm,"time", real([(i,i=(year_now-(y-1)*n_accel),(year_now),(n_accel))],wp), &
    dim1=dim_time,start=[nout],count=[y],ncid=ncid)    
    call nc_write(fnm,"eccentricity", vars%eccentricity, dims=[dim_time],start=[nout],count=[y],long_name="eccentricity",units="1",ncid=ncid)
    call nc_write(fnm,"precession", vars%precession, dims=[dim_time],start=[nout],count=[y],long_name="climatic precession",units="1",ncid=ncid)
    call nc_write(fnm,"obliquity", vars%obliquity, dims=[dim_time],start=[nout],count=[y],long_name="obliquity",units="degrees",ncid=ncid)
    call nc_write(fnm,"Smax65N ", vars%Smax65N , dims=[dim_time],start=[nout],count=[y],long_name="maximum insolation at 65N",units="W/m2",ncid=ncid)
    call nc_write(fnm,"Smax65S ", vars%Smax65S , dims=[dim_time],start=[nout],count=[y],long_name="maximum insolation at 65S",units="W/m2",ncid=ncid)
    call nc_write(fnm,"Smax30N ", vars%Smax30N , dims=[dim_time],start=[nout],count=[y],long_name="maximum insolation at 30N",units="W/m2",ncid=ncid)
    call nc_write(fnm,"Smax30S ", vars%Smax30S , dims=[dim_time],start=[nout],count=[y],long_name="maximum insolation at 30S",units="W/m2",ncid=ncid)
    call nc_write(fnm,"co2     ", vars%co2     , dims=[dim_time],start=[nout],count=[y],long_name="atmospheric CO2 concentration for radiation",units="ppm",ncid=ncid)
    call nc_write(fnm,"co2e    ", vars%co2e    , dims=[dim_time],start=[nout],count=[y],long_name="equivalent atmospheric CO2 concentration for radiation",units="ppm",ncid=ncid)
    if (l_co2d) then
    call nc_write(fnm,"co2flx  ", vars%co2flx  , dims=[dim_time],start=[nout],count=[y],long_name="global carbon flux to atmosphere",units="PgC/yr",ncid=ncid)
    call nc_write(fnm,"co2d_avg", vars%co2d_avg, dims=[dim_time],start=[nout],count=[y],long_name="global mean atmospheric CO2 concentration from 2D field",units="ppm",ncid=ncid)
    call nc_write(fnm,"co2d_wais", vars%co2d_wais, dims=[dim_time],start=[nout],count=[y],long_name="atmospheric CO2 concentration at WAIS ice divide",units="ppm",ncid=ncid)
    endif
    call nc_write(fnm,"ch4     ", vars%ch4     , dims=[dim_time],start=[nout],count=[y],long_name="atmospheric CH4 concentration",units="ppb",ncid=ncid)
    call nc_write(fnm,"n2o     ", vars%n2o     , dims=[dim_time],start=[nout],count=[y],long_name="atmospheric N2O concentration",units="ppb",ncid=ncid)
    call nc_write(fnm,"so4     ", vars%so4     , dims=[dim_time],start=[nout],count=[y],long_name="atmospheric SO4 content",units="kg",ncid=ncid)
    call nc_write(fnm,"tg      ", vars%t2m     , dims=[dim_time],start=[nout],count=[y],long_name="global surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tnh     ", vars%tnh     , dims=[dim_time],start=[nout],count=[y],long_name="NH surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tsh     ", vars%tsh     , dims=[dim_time],start=[nout],count=[y],long_name="SH surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"dtnsh   ", vars%dtnsh   , dims=[dim_time],start=[nout],count=[y],long_name="difference between mean NH and SH surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tn6090  ", vars%tn6090  , dims=[dim_time],start=[nout],count=[y],long_name="surface air temperature between 60-90N",units="degC",ncid=ncid)
    call nc_write(fnm,"ts6090  ", vars%ts6090  , dims=[dim_time],start=[nout],count=[y],long_name="surface air temperature between 60-90S",units="degC",ncid=ncid)
    call nc_write(fnm,"tgrl    ", vars%tgrl    , dims=[dim_time],start=[nout],count=[y],long_name="Greenland surface air temperature NGRIP",units="degC",ncid=ncid)
    call nc_write(fnm,"tgrl1   ", vars%tgrl1   , dims=[dim_time],start=[nout],count=[y],long_name="Greenland surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tgrl1_jja", vars%tgrl1_jja, dims=[dim_time],start=[nout],count=[y],long_name="JJA Greenland surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tnatl   ", vars%tnatl   , dims=[dim_time],start=[nout],count=[y],long_name="North Atlantic surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tant    ", vars%tant    , dims=[dim_time],start=[nout],count=[y],long_name="Antarctica surface air temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"dust_grl", vars%dust_grl, dims=[dim_time],start=[nout],count=[y],long_name="Greenland dust deposition",units="kg/m2/s",ncid=ncid)
    call nc_write(fnm,"dust_wais", vars%dust_wais , dims=[dim_time],start=[nout],count=[y],long_name="WAIS dust deposition",units="kg/m2/s",ncid=ncid)
    call nc_write(fnm,"tam     ", vars%tam     , dims=[dim_time],start=[nout],count=[y],long_name="temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tskin   ", vars%tskin   , dims=[dim_time],start=[nout],count=[y],long_name="skin temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"q2      ", vars%q2m     , dims=[dim_time],start=[nout],count=[y],long_name="surface air humidity",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"r2      ", vars%r2m     , dims=[dim_time],start=[nout],count=[y],long_name="surface air relative humidity",units="1",ncid=ncid)
    call nc_write(fnm,"prc     ", vars%prc     , dims=[dim_time],start=[nout],count=[y],long_name="precipitation",units="10^15 kg",ncid=ncid)
    call nc_write(fnm,"evp     ", vars%evp     , dims=[dim_time],start=[nout],count=[y],long_name="evaporation",units="10^15 kg",ncid=ncid)
    call nc_write(fnm,"wcon    ", vars%wcon    , dims=[dim_time],start=[nout],count=[y],long_name="water content",units="",ncid=ncid)
    call nc_write(fnm,"cld     ", vars%cld     , dims=[dim_time],start=[nout],count=[y],long_name="cloud fraction",units="",ncid=ncid)
    call nc_write(fnm,"sol     ", vars%sol     , dims=[dim_time],start=[nout],count=[y],long_name="shortwave in at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"snettop ", vars%snettop , dims=[dim_time],start=[nout],count=[y],long_name="shortwave net at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"snetsur ", vars%snetsur , dims=[dim_time],start=[nout],count=[y],long_name="shortwave net at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lnettop ", vars%lnettop , dims=[dim_time],start=[nout],count=[y],long_name="longwave net at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lnetsur ", vars%lnetsur , dims=[dim_time],start=[nout],count=[y],long_name="longwave net at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"ldwnsur ", vars%ldwnsur , dims=[dim_time],start=[nout],count=[y],long_name="longwave down at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lupsur  ", vars%lupsur  , dims=[dim_time],start=[nout],count=[y],long_name="longwave up at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbtop   ", vars%rbtop   , dims=[dim_time],start=[nout],count=[y],long_name="radiative balance at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbatm   ", vars%rbatm   , dims=[dim_time],start=[nout],count=[y],long_name="radiative balance of the atmosphere",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbsur   ", vars%rbsur   , dims=[dim_time],start=[nout],count=[y],long_name="radiative balance at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"snettopcs ", vars%snettopcs , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky shortwave net at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"snetsurcs ", vars%snetsurcs , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky shortwave net at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lnettopcs ", vars%lnettopcs , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky longwave net at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lnetsurcs ", vars%lnetsurcs , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky longwave net at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"ldwnsurcs ", vars%ldwnsurcs , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky longwave down at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbtopcs   ", vars%rbtopcs   , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky radiative balance at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbatmcs   ", vars%rbatmcs   , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky radiative balance of the atmosphere",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbsurcs   ", vars%rbsurcs   , dims=[dim_time],start=[nout],count=[y],long_name="clear-sky radiative balance at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"sh      ", vars%sha     , dims=[dim_time],start=[nout],count=[y],long_name="surface sensible heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lh      ", vars%lha     , dims=[dim_time],start=[nout],count=[y],long_name="surface latent heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"slp     ", vars%slp     , dims=[dim_time],start=[nout],count=[y],long_name="sea level pressure",units="hPa",ncid=ncid)
    call nc_write(fnm,"swrcrf  ", vars%swrcrf  , dims=[dim_time],start=[nout],count=[y],long_name="shortwave cloud radiative forcing",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwrcrf  ", vars%lwrcrf  , dims=[dim_time],start=[nout],count=[y],long_name="longwave cloud radiative forcing",units="W/m2",ncid=ncid)
    call nc_write(fnm,"htrop   ", vars%htrop   , dims=[dim_time],start=[nout],count=[y],long_name="tropopause height",units="m",ncid=ncid)
    call nc_write(fnm,"ttrop   ", vars%ttrop   , dims=[dim_time],start=[nout],count=[y],long_name="tropopause temperature",units="K",ncid=ncid)
    call nc_write(fnm,"hcld    ", vars%hcld    , dims=[dim_time],start=[nout],count=[y],long_name="cloud height",units="m",ncid=ncid)
    call nc_write(fnm,"plan    ", vars%plan    , dims=[dim_time],start=[nout],count=[y],long_name="planetary albedo",units="",ncid=ncid)
    call nc_write(fnm,"dustload", vars%dustload, dims=[dim_time],start=[nout],count=[y],long_name="dust load",units="Tg",ncid=ncid)
    call nc_write(fnm,"dustdep ", vars%dustdep , dims=[dim_time],start=[nout],count=[y],long_name="dust deposition",units="Tg/yr",ncid=ncid)
    call nc_write(fnm,"dustdrydep", vars%dustdrydep, dims=[dim_time],start=[nout],count=[y],long_name="dust dry deposition",units="Tg/yr",ncid=ncid)
    call nc_write(fnm,"dustwetdep", vars%dustwetdep, dims=[dim_time],start=[nout],count=[y],long_name="dust wet deposition",units="Tg/yr",ncid=ncid)
    call nc_write(fnm,"fw_pac_atl", vars%fw_pac_atl , dims=[dim_time],start=[nout],count=[y],long_name="water vapor transport from Pacific to Atlantic",units="Sv",ncid=ncid)
    call nc_write(fnm,"fw_atl_indpac", vars%fw_atl_indpac , dims=[dim_time],start=[nout],count=[y],long_name="water vapor transport from Atlantic to Indo-Pacific",units="Sv",ncid=ncid)
    call nc_write(fnm,"itcz", vars%itcz , dims=[dim_time],start=[nout],count=[y],long_name="Latitude of the ITCZ (annual average)",units="deg",ncid=ncid)
    call nc_write(fnm,"itcz_min", vars%itcz_min , dims=[dim_time],start=[nout],count=[y],long_name="Min latitude of the ITCZ",units="deg",ncid=ncid)
    call nc_write(fnm,"itcz_max", vars%itcz_max , dims=[dim_time],start=[nout],count=[y],long_name="Max latitude of the ITCZ",units="deg",ncid=ncid)
    call nc_write(fnm,"had_width", vars%had_width , dims=[dim_time],start=[nout],count=[y],long_name="Width of the Hadley cells (annual average)",units="deg",ncid=ncid)

    call nc_write(fnm,"tg_l      ", vars%t2m_l     , dims=[dim_time],start=[nout],count=[y],long_name="temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tam_l     ", vars%tam_l     , dims=[dim_time],start=[nout],count=[y],long_name="temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"tskin_l   ", vars%tskin_l   , dims=[dim_time],start=[nout],count=[y],long_name="temperature",units="degC",ncid=ncid)
    call nc_write(fnm,"q2_l      ", vars%q2m_l     , dims=[dim_time],start=[nout],count=[y],long_name="surface air humidity",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"r2_l      ", vars%r2m_l     , dims=[dim_time],start=[nout],count=[y],long_name="surface air relative humidity",units="1",ncid=ncid)
    call nc_write(fnm,"prc_l     ", vars%prc_l     , dims=[dim_time],start=[nout],count=[y],long_name="precipitation",units="10^15 kg",ncid=ncid)
    call nc_write(fnm,"evp_l     ", vars%evp_l     , dims=[dim_time],start=[nout],count=[y],long_name="evaporation",units="10^15 kg",ncid=ncid)
    call nc_write(fnm,"wcon_l    ", vars%wcon_l    , dims=[dim_time],start=[nout],count=[y],long_name="water content",units="",ncid=ncid)
    call nc_write(fnm,"cld_l     ", vars%cld_l     , dims=[dim_time],start=[nout],count=[y],long_name="cloud fraction",units="",ncid=ncid)
    call nc_write(fnm,"snettop_l ", vars%snettop_l , dims=[dim_time],start=[nout],count=[y],long_name="shortwave net at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"snetsur_l ", vars%snetsur_l , dims=[dim_time],start=[nout],count=[y],long_name="shortwave net at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lnettop_l ", vars%lnettop_l , dims=[dim_time],start=[nout],count=[y],long_name="longwave net at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lnetsur_l ", vars%lnetsur_l , dims=[dim_time],start=[nout],count=[y],long_name="longwave net at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbtop_l   ", vars%rbtop_l   , dims=[dim_time],start=[nout],count=[y],long_name="radiative balance at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbatm_l   ", vars%rbatm_l   , dims=[dim_time],start=[nout],count=[y],long_name="radiative balance of the atmosphere",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rbsur_l   ", vars%rbsur_l   , dims=[dim_time],start=[nout],count=[y],long_name="radiative balance at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"sh_l      ", vars%sha_l     , dims=[dim_time],start=[nout],count=[y],long_name="surface sensible heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lh_l      ", vars%lha_l     , dims=[dim_time],start=[nout],count=[y],long_name="surface latent heat flux",units="W/m2",ncid=ncid)


    call nc_close(ncid)


   return

  end subroutine ts_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ n c
  ! Author   :  Matteo Willeit
  ! Purpose  :  Initialize atmosphere netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_nc(fnm)

    implicit none

    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,"time", x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE., ncid=ncid)
    call nc_write_dim(fnm, "mon", x=1._wp, dx=1._wp, nx=13, units="months", &
    axis="e", ncid=ncid)
    call nc_write_dim(fnm,"lat", x=lat, axis="y", ncid=ncid)
    call nc_write_dim(fnm,"lon", x=lon, axis="x", ncid=ncid)
    call nc_write_dim(fnm,"doy",x=1._wp,dx=1._wp,nx=nday_year,ncid=ncid)
    call nc_write_dim(fnm,"hour",x=1._wp,dx=1._wp,nx=24,ncid=ncid)
    call nc_write_dim(fnm,"zlay",x=zc(1:km),units="m",ncid=ncid)
    call nc_write_dim(fnm,"zlev",x=zl(1:km),units="m",ncid=ncid)
    call nc_write_dim(fnm,"zlevw",x=zl(1:kmc),units="m",ncid=ncid)
    call nc_write_dim(fnm,"plev",x=pl(1:km)*1000._wp,units="hPa",ncid=ncid)
    call nc_write_dim(fnm,"plevw",x=pl(1:kmc)*1000._wp,units="hPa",ncid=ncid)
    call nc_write_dim(fnm,"lonu",x=lon0,dx=dlon,nx=imc,ncid=ncid)
    call nc_write_dim(fnm,"latv",x=lat0,dx=dlat,nx=jmc,ncid=ncid)
    call nc_write_dim(fnm,"st",x=1._wp,dx=1._wp,nx=nm,ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine atm_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ n c _ w r i t e
  ! Author   :  Matteo Willeit
  ! Purpose  :  Output of atmosphere netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_nc_write(fnm,ncid,vars,ndat,nout)

    implicit none

    type(a_out) :: vars

    character (len=*) :: fnm
    integer :: ndat, nout, ncid


    if (ndat.eq.13) then

      call nc_write(fnm,"solar      ", sngl(vars%solar     (:,:,jm:1:-1)), dims=["doy ","hour","lat ","time"], &
        start=[1,1,1,nout],count=[nday_year,24,jm,1],long_name="TOA incoming solar radiation",units="W/m2",ncid=ncid)
      call nc_write(fnm,"cosz       ", sngl(vars%cosz      (:,:,jm:1:-1)), dims=["doy ","hour","lat ","time"], &
        start=[1,1,1,nout],count=[nday_year,24,jm,1],long_name="cosine of solar zenith angle",units="",ncid=ncid)
      call nc_write(fnm,"solarm     ", sngl(vars%solarm    (:,jm:1:-1)), dims=["doy ","lat ","time"], &
        start=[1,1,nout],count=[nday_year,jm,1],long_name="daily mean TOA incoming solar radiation",units="W/m2",ncid=ncid)
      call nc_write(fnm,"coszm      ", sngl(vars%coszm     (:,jm:1:-1)), dims=["doy ","lat ","time"], &
        start=[1,1,nout],count=[nday_year,jm,1],long_name="radiation weighted daily mean cosine of solar zenith angle",units="",ncid=ncid)

      call nc_write(fnm,"zsa        ", sngl(vars%zsa       (:,jm:1:-1) ), dims=["lon ","lat ","time"], &
        start=[1,1,nout],count=[im,jm,1],long_name="grid cell average surface elevation",units="m",ncid=ncid)
      call nc_write(fnm,"zsa_smooth ", sngl(vars%zsa_smooth(:,jm:1:-1) ), dims=["lon ","lat ","time"], &
        start=[1,1,nout],count=[im,jm,1],long_name="grid cell average smoothed surface elevation",units="m",ncid=ncid)
      call nc_write(fnm,"slope      ", sngl(vars%slope     (:,jm:1:-1) ), dims=["lon ","lat ","time"], &
        start=[1,1,nout],count=[im,jm,1],long_name="topography slope",units="m/m",ncid=ncid)
      call nc_write(fnm,"slope_x    ", sngl(vars%slope_x   (:,jm:1:-1) ), dims=["lon ","lat ","time"], &
        start=[1,1,nout],count=[im,jm,1],long_name="topography slope in zonal direction",units="m/m",ncid=ncid)
      call nc_write(fnm,"slope_y    ", sngl(vars%slope_y   (:,jm:1:-1) ), dims=["lon ","lat ","time"], &
        start=[1,1,nout],count=[im,jm,1],long_name="topography slope in meridional direction",units="m/m",ncid=ncid)
      call nc_write(fnm,"frlnd      ", sngl(vars%frlnd     (:,jm:1:-1) ), dims=["lon ","lat ","time"], &
        start=[1,1,nout],count=[im,jm,1],long_name="land fraction in grid cell",units="1",ncid=ncid)
      call nc_write(fnm,"frocn      ", sngl(vars%frocn     (:,jm:1:-1) ), dims=["lon ","lat ","time"], &
        start=[1,1,nout],count=[im,jm,1],long_name="ocean fraction in grid cell",units="1",ncid=ncid)
      call nc_write(fnm,"zs         ", sngl(vars%zs        (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","time"], &
        start=[1,1,1,nout],count=[im,jm,nm,1],long_name="surface elevation for each macro surface type",units="m",ncid=ncid)

      call nc_write(fnm,"fw_pac_atl    ", sngl(vars%fw_pac_atl(jm:1:-1)),    dims=["lat ","time"], &
        start=[1,nout],count=[jm,1],long_name="freshwater transport from Pacific to Atlantic catchment basin",units="Sv",ncid=ncid)
      call nc_write(fnm,"fw_atl_indpac ", sngl(vars%fw_atl_indpac(jm:1:-1)), dims=["lat ","time"], &
        start=[1,nout],count=[jm,1],long_name="freshwater transport from Atlantic to Indo-Pacific catchment basin",units="Sv",ncid=ncid)

    endif

    call nc_write(fnm,"had_fi     ", sngl(vars%had_fi     ), dims=["mon ","time"], &
      start=[ndat,nout],count=[1,1],long_name="",units="deg",ncid=ncid)
    call nc_write(fnm,"had_width  ", sngl(vars%had_width  ), dims=["mon ","time"], &
      start=[ndat,nout],count=[1,1],long_name="",units="deg",ncid=ncid)

    call nc_write(fnm,"ptrop      ", sngl(vars%ptrop      (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="relative pressure at tropopause",units="/",ncid=ncid)
    call nc_write(fnm,"tslz       ", sngl(vars%tslz       (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="zonal mean sea level temperature",units="K",ncid=ncid)
    call nc_write(fnm,"tskslz      ", sngl(vars%tskslz      (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="zonal mean skin temperature reduced to sea level",units="K",ncid=ncid)
    call nc_write(fnm,"ekez       ", sngl(vars%ekez       (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="zonal mean eddy kinetic energy",units="m2/s2",ncid=ncid)
    call nc_write(fnm,"slpz       ", sngl(vars%slpz       (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="zonal mean sea level pressure",units="Pa",ncid=ncid)
    call nc_write(fnm,"vabz       ", sngl(vars%vabz       (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="zonal mean meridional ageostrophic wind in PBL",units="m/s",ncid=ncid)
    call nc_write(fnm,"acbarz     ", sngl(vars%acbarz     (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="zonal mean cross-isobar angle",units="rad",ncid=ncid)
    call nc_write(fnm,"uz850      ", sngl(vars%uz850      (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="zonal mean 850 hPa zonal wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"uz500      ", sngl(vars%uz500      (jm:1:-1)), dims=["lat ","mon ","time"], &
      start=[1,ndat,nout],count=[jm,1,1],long_name="zonal mean 500 hPa zonal wind",units="m/s",ncid=ncid)

    call nc_write(fnm,"tam        ", sngl(vars%tam       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="prognostic atmospheric temperature",units="K",ncid=ncid)
    call nc_write(fnm,"qam        ", sngl(vars%qam       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="prognostic atmospheric specific humidity",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"ram        ", sngl(vars%ram       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="prognostic atmospheric relative humidity",units="1",ncid=ncid)
    call nc_write(fnm,"gams       ", sngl(vars%gams      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="temperature lapse rate in the boundary layer",units="K/m",ncid=ncid)
    call nc_write(fnm,"gamb       ", sngl(vars%gamb      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="temperature lapse rate at the surface",units="K/m",ncid=ncid)
    call nc_write(fnm,"gamt       ", sngl(vars%gamt      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="emperature lapse rate at 15 km",units="K/m",ncid=ncid)
    call nc_write(fnm,"dam        ", sngl(vars%dam       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="prognostic atmospheric dust mass mixing ratio",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"hqeff      ", sngl(vars%hqeff     (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="effective specific humidity height scale",units="m",ncid=ncid)
    call nc_write(fnm,"hrm        ", sngl(vars%hrm       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="relative humidity height scale",units="m",ncid=ncid)
    call nc_write(fnm,"wcon       ", sngl(vars%wcon      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated water content",units="kg/m2",ncid=ncid)
    call nc_write(fnm,"cld_rh        ", sngl(vars%cld_rh       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="large-scale cloud fraction in grid cell",units="1",ncid=ncid)
    call nc_write(fnm,"cld_low       ", sngl(vars%cld_low      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="low cloud fraction in grid cell",units="1",ncid=ncid)
    call nc_write(fnm,"cld        ", sngl(vars%cld       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="total cloud fraction in grid cell",units="1",ncid=ncid)
    call nc_write(fnm,"prc        ", sngl(vars%prc       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="total precipitation",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"prcw       ", sngl(vars%prcw      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="rainfall",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"prcs       ", sngl(vars%prcs      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="snowfall",units="kg/m2/day (water equivalent)",ncid=ncid)
    call nc_write(fnm,"hcld       ", sngl(vars%hcld      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="cloud top height",units="m",ncid=ncid)
    call nc_write(fnm,"ctt        ", sngl(vars%ctt       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="clout top temperature",units="K",ncid=ncid)
    call nc_write(fnm,"clot       ", sngl(vars%clot      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="cloud optical thickness",units="",ncid=ncid)
    call nc_write(fnm,"alb_cld    ", sngl(vars%alb_cld   (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="cloud albedo",units="1",ncid=ncid)
    call nc_write(fnm,"htrop      ", sngl(vars%htrop     (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="height of tropopause",units="m",ncid=ncid)
    call nc_write(fnm,"ttrop      ", sngl(vars%ttrop     (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="temperature of stratosphere",units="K",ncid=ncid)
    call nc_write(fnm,"wind       ", sngl(vars%wind      (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="10 m wind speed",units="m/s",ncid=ncid)

    call nc_write(fnm,"frst       ", sngl(vars%frst      (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="fraction of macro surface types",units="1",ncid=ncid)
    call nc_write(fnm,"tskin      ", sngl(vars%tskin     (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="skin temperature for each macro surface type",units="K",ncid=ncid)
    call nc_write(fnm,"t2         ", sngl(vars%t2        (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="2m temperature for each macro surface type",units="K",ncid=ncid)
    call nc_write(fnm,"q2         ", sngl(vars%q2        (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="2m specific humidity for each macro surface type",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"r2         ", sngl(vars%r2        (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="2m relative humidity for each macro surface type",units="1",ncid=ncid)
    call nc_write(fnm,"ra2        ", sngl(vars%ra2       (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="surface air density",units="kg/m2",ncid=ncid)
    call nc_write(fnm,"alb_vu_s   ", sngl(vars%alb_vu_s  (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="visible clear-sky surface albedo for each macro surface type",units="1",ncid=ncid)
    call nc_write(fnm,"alb_vu_c   ", sngl(vars%alb_vu_c  (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="visible diffuse surface albedo for each macro surface type",units="1",ncid=ncid)
    call nc_write(fnm,"alb_ir_s   ", sngl(vars%alb_ir_s  (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="near-infrared clear sky surface albedo for each macro surface type",units="1",ncid=ncid)
    call nc_write(fnm,"alb_ir_c   ", sngl(vars%alb_ir_c  (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="near-infrared diffuse surface albedo for each macro surface type",units="1",ncid=ncid)
    call nc_write(fnm,"cd         ", sngl(vars%cd        (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="drag coefficient for each macro surface type",units="",ncid=ncid)
    call nc_write(fnm,"taux       ", sngl(vars%taux      (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="zonal surface wind stress",units="N/m2",ncid=ncid)
    call nc_write(fnm,"tauy       ", sngl(vars%tauy      (:,jm:1:-1,:) ), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="meridional surface wind stress",units="N/m2",ncid=ncid)

    call nc_write(fnm,"tsl        ", sngl(vars%tsl       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="sea level temperature",units="K",ncid=ncid)
    call nc_write(fnm,"tskina     ", sngl(vars%tskina    (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell mean skin temperature",units="K",ncid=ncid)
    call nc_write(fnm,"t2a        ", sngl(vars%t2a       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell mean 2m temperature",units="K",ncid=ncid)
    call nc_write(fnm,"thetae     ", sngl(vars%thetae    (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell mean 2m temperature",units="K",ncid=ncid)
    call nc_write(fnm,"q2a        ", sngl(vars%q2a       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell mean 2m specific humidity",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"dq         ", sngl(vars%dq        (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="specific humidity gradient between saturated skin and 2m",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"dr         ", sngl(vars%dr        (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="specific humidity gradient between saturated skin and 2m",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"r2a        ", sngl(vars%r2a       (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell mean 2m relative humidity",units="1",ncid=ncid)
    call nc_write(fnm,"rskina     ", sngl(vars%rskina    (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell mean skin relative humidity",units="1",ncid=ncid)
    call nc_write(fnm,"rskin_ram  ", sngl(vars%rskina(:,jm:1:-1)-vars%ram(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="skin - atmospheric relative humidity",units="1",ncid=ncid)

    if (ndat.ne.13) then
      call nc_write(fnm,"t2m_bias   ", sngl(vars%t2a(:,jm:1:-1)-vars%t2a_dat(:,1:jm) ), dims=["lon ","lat ","mon ","time"], &
        start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell mean 2m temperature",units="K",ncid=ncid)
    endif

    call nc_write(fnm,"tskin_t2   ", sngl(vars%tskina(:,jm:1:-1)-vars%t2a(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="skin temperature - 2m temperature",units="K",ncid=ncid)
    call nc_write(fnm,"tskin_tam  ", sngl(vars%tskina(:,jm:1:-1)-vars%tam(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="skin temperature - atmospheric temperature",units="K",ncid=ncid)

    call nc_write(fnm,"t3         ", sngl(vars%t3        (:,jm:1:-1,:) ), dims=["lon  ","lat  ","zlevw","mon  ","time "], &
      start=[1,1,1,ndat,nout],count=[im,jm,kmc,1,1],long_name="3D atmospheric temperature",units="K",ncid=ncid)
    call nc_write(fnm,"q3         ", sngl(vars%q3        (:,jm:1:-1,:) ), dims=["lon  ","lat  ","zlevw","mon  ","time "], &
      start=[1,1,1,ndat,nout],count=[im,jm,kmc,1,1],long_name="3D specific humidity",units="kg/kg",ncid=ncid)
    call nc_write(fnm,"r3         ", sngl(vars%r3        (:,jm:1:-1,:) ), dims=["lon  ","lat  ","zlevw","mon  ","time "], &
      start=[1,1,1,ndat,nout],count=[im,jm,kmc,1,1],long_name="3D relative humidity",units="1",ncid=ncid)
    call nc_write(fnm,"tp         ", sngl(vars%tp        (:,jm:1:-1,:) ), dims=["lon  ","lat  ","zlevw","mon  ","time "], &
      start=[1,1,1,ndat,nout],count=[im,jm,kmc,1,1],long_name="3D potential temperature",units="K",ncid=ncid)
    call nc_write(fnm,"gamma      ", sngl(vars%gamma     (:,jm:1:-1,:) ), dims=["lon ","lat ","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="3D atmospheric temperature lapse rate",units="K/m",ncid=ncid)
    call nc_write(fnm,"rho        ", sngl(vars%rho       (:,jm:1:-1,:) ), dims=["lon  ","lat  ","zlevw","mon  ","time "], &
      start=[1,1,1,ndat,nout],count=[im,jm,kmc,1,1],long_name="3D atmospheric density",units="kg/m3",ncid=ncid)

    call nc_write(fnm,"convdse    ", sngl(vars%convdse    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated dry static energy convergence",units="W/m2",ncid=ncid)
    call nc_write(fnm,"convadse   ", sngl(vars%convadse   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated dry static energy convergence by mean circulation",units="W/m2",ncid=ncid)
    call nc_write(fnm,"convddse   ", sngl(vars%convddse   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated dry static energy convergence by synoptic eddies",units="",ncid=ncid)
    call nc_write(fnm,"convwtr    ", sngl(vars%convwtr    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated moisture convergence",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"convawtr   ", sngl(vars%convawtr   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated moisture convergence by mean circulation",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"convdwtr   ", sngl(vars%convdwtr   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated moisture convergence by synotic eddies",units="kg/m2/s",ncid=ncid)

    call nc_write(fnm,"ra2a       ", sngl(vars%ra2a       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average surface air density",units="kg/m3",ncid=ncid)
    call nc_write(fnm,"cda        ", sngl(vars%cda        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average drag coefficient",units="",ncid=ncid)
    call nc_write(fnm,"cd0a       ", sngl(vars%cd0a       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average drag coefficient without orographic component",units="",ncid=ncid)
    call nc_write(fnm,"sha        ", sngl(vars%sha        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average surface sensible heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lha        ", sngl(vars%lha        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average surface latent heat flux",units="W/m2",ncid=ncid)
    call nc_write(fnm,"evpa       ", sngl(vars%evpa       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average evaporation",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"Ri         ", sngl(vars%Ri         (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="bulk Richardson number",units="1",ncid=ncid)
    call nc_write(fnm,"acbar      ", sngl(vars%acbar      (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="cross-isobar angle",units="rad",ncid=ncid)
    call nc_write(fnm,"epsa       ", sngl(vars%epsa       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="epsilon",units="",ncid=ncid)
    call nc_write(fnm,"tsksl       ", sngl(vars%tsksl       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal sea level temperature",units="K",ncid=ncid)
    call nc_write(fnm,"atsksl       ", sngl(vars%atsksl       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal sea level temperature",units="K",ncid=ncid)
    call nc_write(fnm,"atsl       ", sngl(vars%atsl       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal sea level temperature",units="K",ncid=ncid)
    call nc_write(fnm,"atsli        ", sngl(vars%atsli       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal sea level temperature",units="K",ncid=ncid)
    call nc_write(fnm,"aslp       ", sngl(vars%aslp       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal sea level pressure",units="Pa",ncid=ncid)
    call nc_write(fnm,"aslp_temp  ", sngl(vars%aslp_temp  (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal sea level pressure from azonal temperature",units="Pa",ncid=ncid)
    call nc_write(fnm,"aslp_topo  ", sngl(vars%aslp_topo  (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal sea level pressure from topographic stationary waves",units="Pa",ncid=ncid)
    call nc_write(fnm,"slp        ", sngl(vars%slp        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="sea level pressure",units="Pa",ncid=ncid)
    call nc_write(fnm,"slp1       ", sngl(vars%slp1       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="sea level pressure derived using constant 6.5 K/km lapse rate",units="Pa",ncid=ncid)
    call nc_write(fnm,"dz500      ", sngl(vars%dz500      (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="azonal geopotential height at 500 hPa",units="m",ncid=ncid)
    call nc_write(fnm,"us        ", sngl(vars%us        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="zonal 10m wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"vs        ", sngl(vars%vs        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="meridional 10m wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"usk       ", sngl(vars%usk       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="zonal katabatic surface wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"vsk       ", sngl(vars%vsk       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="meridional katabatic surface wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"ugb       ", sngl(vars%ugb       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="zonal geostrophic wind in the boundary layer",units="m/s",ncid=ncid)
    call nc_write(fnm,"vgb       ", sngl(vars%vgb       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="meridional geostrophic wind in the boundary layer",units="m/s",ncid=ncid)
    call nc_write(fnm,"uab       ", sngl(vars%uab       (:,jm:1:-1)), dims=["lonu","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[imc,jm,1,1],long_name="zonal ageostrophic wind in the boundary layer",units="m/s",ncid=ncid)
    call nc_write(fnm,"vab       ", sngl(vars%vab       (:,jmc:1:-1)), dims=["lon ","latv","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jmc,1,1],long_name="meridional ageostrophic wind in the boundary layer",units="m/s",ncid=ncid)
    call nc_write(fnm,"wcld       ", sngl(vars%wcld       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="vertical velocity at cloud level",units="m/s",ncid=ncid)
    call nc_write(fnm,"woro       ", sngl(vars%woro       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="vertical velocity induced by sub-grid orography",units="m/s",ncid=ncid)
    call nc_write(fnm,"wsyn       ", sngl(vars%wsyn       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="synoptic vertical velocity at 700 hPa",units="m/s",ncid=ncid)
    call nc_write(fnm,"weff       ", sngl(vars%weff       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="effective vertical velocity",units="m/s",ncid=ncid)
    call nc_write(fnm,"fweff      ", sngl(vars%fweff      (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="effective vertical velocity factor for cloud parameterisation",units="m/s",ncid=ncid)

    call nc_write(fnm,"rb_top     ", sngl(vars%rb_top     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="radiative balance at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rb_sur     ", sngl(vars%rb_sur     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="radiative balance at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rb_atm     ", sngl(vars%rb_atm     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="radiative balance of the atmosphere",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_atm   ", sngl(vars%swr_atm   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net shortwave radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_top   ", sngl(vars%swr_top   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net shortwave radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_top_cs   ", sngl(vars%swr_top_cs   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net shortwave clear-sky radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_sur   ", sngl(vars%swr_sur   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net shortwave radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_sur_cs", sngl(vars%swr_sur_cs(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net shortwave clear-sky radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_dw_sur", sngl(vars%swr_dw_sur(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average downward shortwave radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_dw_top", sngl(vars%swr_dw_top(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="incoming shortwave radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_atm   ", sngl(vars%lwr_atm   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_atm_cld   ", sngl(vars%lwr_atm_cld   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_top   ", sngl(vars%lwr_top   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_top_cs   ", sngl(vars%lwr_top_cs   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave clear-sky radiation at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_sur   ", sngl(vars%lwr_sur   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_sur_cs   ", sngl(vars%lwr_sur_cs   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave clear-sky radiation at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_tro", sngl(vars%lwr_tro(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave radiation at the tropopause",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_cld", sngl(vars%lwr_cld(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="grid-cell average net longwave radiation at the tropopause",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rb_str     ", sngl(vars%rb_str     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="radiative balance of the stratosphere",units="W/m2",ncid=ncid)
    call nc_write(fnm,"cre_top    ", sngl(vars%cre_top    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="cloud radiative effect at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_cre_top", sngl(vars%swr_cre_top(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="shortwave cloud radiative effect at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_cre_top", sngl(vars%lwr_cre_top(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="longwave cloud radiative effect at TOA",units="W/m2",ncid=ncid)
    call nc_write(fnm,"cre_sur    ", sngl(vars%cre_sur    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="cloud radiative effect at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"swr_cre_sur", sngl(vars%swr_cre_sur(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="shortwave cloud radiative effect at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"lwr_cre_sur", sngl(vars%lwr_cre_sur(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="longwave cloud radiative effect at the surface",units="W/m2",ncid=ncid)
    call nc_write(fnm,"alb_plan   ", sngl(vars%aplan      (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="planetary albedo",units="1",ncid=ncid)
    call nc_write(fnm,"alb_sur_cs ", sngl(vars%alb_sur_cs (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="clear-sky surface albedo",units="1",ncid=ncid)
    call nc_write(fnm,"alb_sur_cld", sngl(vars%alb_sur_cld(:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="cloudy-sky surface albedo",units="1",ncid=ncid)

    call nc_write(fnm,"eke2_6        ", sngl(vars%eke        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="eddy kinetic energy",units="m2/s2",ncid=ncid)
    call nc_write(fnm,"eke        ", sngl(vars%sam        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="eddy kinetic energy",units="m2/s2",ncid=ncid)
    call nc_write(fnm,"eke2       ", sngl(vars%sam2       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="eddy kinetic energy in 2-6 day band",units="m2/s2",ncid=ncid)
    call nc_write(fnm,"ekeprod    ", sngl(vars%synprod    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="production of eddy kinetic energy",units="",ncid=ncid)
    call nc_write(fnm,"ekediss    ", sngl(vars%syndiss    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="dissipation of eddy kinetic energy",units="",ncid=ncid)
    call nc_write(fnm,"ekeadv     ", sngl(vars%synadv     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="advection of eddy kinetic energy",units="",ncid=ncid)
    call nc_write(fnm,"ekedif     ", sngl(vars%syndif     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="diffusion of eddy kinetic energy",units="",ncid=ncid)
    call nc_write(fnm,"wind_syn   ", sngl(vars%synsur     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="synoptic 10m wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"cdif       ", sngl(vars%cdif       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="effective macrodiffusive coefficient",units="m2/s",ncid=ncid)
    call nc_write(fnm,"diffxdse   ", sngl(vars%diffxdse   (:,jm:1:-1)), dims=["lonu","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[imc,jm,1,1],long_name="zonal effective macrodiffusivity for dry static energy",units="m2/s",ncid=ncid)
    call nc_write(fnm,"diffydse   ", sngl(vars%diffydse   (:,jmc:1:-1)), dims=["lon ","latv","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jmc,1,1],long_name="meridional effective macrodiffusivity for dry static energy",units="m2/s",ncid=ncid)
    call nc_write(fnm,"diffxwtr   ", sngl(vars%diffxwtr   (:,jm:1:-1)), dims=["lonu","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[imc,jm,1,1],long_name="zonal effective macrodiffusivity for dry water vapor",units="m2/s",ncid=ncid)
    call nc_write(fnm,"diffywtr   ", sngl(vars%diffywtr   (:,jmc:1:-1)), dims=["lon ","latv","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jmc,1,1],long_name="meridional effective macrodiffusivity for water vapor",units="m2/s",ncid=ncid)

    call nc_write(fnm,"u3         ", sngl(vars%u3        (:,jm:1:-1,:)), dims=["lon ","lat ","zlev","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="3D zonal wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"v3         ", sngl(vars%v3        (:,jm:1:-1,:)), dims=["lon ","lat ","zlev","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="3D meridional wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"w3         ", sngl(vars%w3        (:,jm:1:-1,:)), dims=["lon  ","lat  ","zlevw","mon  ","time "], &
      start=[1,1,1,ndat,nout],count=[im,jm,kmc,1,1],long_name="3D vertical velocity",units="m/s",ncid=ncid)

    call nc_write(fnm,"xz         ", sngl(vars%xz         (jmc:1:-1,:)), dims=["latv ","zlevw","mon  ","time "], &
      start=[1,1,ndat,nout],count=[jmc,kmc,1,1],long_name="mean meridional mass streamfunction",units="10**10 kg/s",ncid=ncid)

    call nc_write(fnm,"fayg      ", sngl(vars%fayg      (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"faydseg   ", sngl(vars%faydseg   (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional dry static energy transport by mean circulation",units="PW",ncid=ncid)
    call nc_write(fnm,"faycptg   ", sngl(vars%faycptg   (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional cp*T transport by mean circulation",units="PW",ncid=ncid)
    call nc_write(fnm,"faygzg    ", sngl(vars%faygzg    (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional geopotential energy transport by mean circulation",units="PW",ncid=ncid)
    call nc_write(fnm,"fayleg    ", sngl(vars%fayleg    (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional latent heat transport by mean circulation",units="PW",ncid=ncid)
    call nc_write(fnm,"faywtrg   ", sngl(vars%faywtrg   (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional moisture transport by mean circulation",units="kg/s",ncid=ncid)
    call nc_write(fnm,"fdydseg   ", sngl(vars%fdydseg   (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional dry static energy transport by synoptic eddies",units="PW",ncid=ncid)
    call nc_write(fnm,"fdyleg    ", sngl(vars%fdyleg    (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional latent heat transport by synoptic eddies",units="PW",ncid=ncid)
    call nc_write(fnm,"fdywtrg   ", sngl(vars%fdywtrg   (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="meridional moisture transport by synoptic eddies",units="kg/s",ncid=ncid)
    call nc_write(fnm,"fydseg    ", sngl(vars%fydseg    (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="total meridional dry static energy transport",units="PW",ncid=ncid)
    call nc_write(fnm,"fyheatg   ", sngl(vars%fyheatg   (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="total meridional heat transport",units="PW",ncid=ncid)
    call nc_write(fnm,"fyleg     ", sngl(vars%fyleg     (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="total meridional latent heat transport",units="PW",ncid=ncid)
    call nc_write(fnm,"fywtrg    ", sngl(vars%fywtrg    (jmc:1:-1)), dims=["latv","mon ","time"], &
      start=[1,ndat,nout],count=[jmc,1,1],long_name="total meridional moisture transport",units="kg/s",ncid=ncid)

    call nc_write(fnm,"fswr_sur       ", sngl(vars%fswr_sur       (:,jm:1:-1,:)), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="net surface shortwave radiation for each surface type",units="W/m2",ncid=ncid)
    call nc_write(fnm,"flwr_dw_sur    ", sngl(vars%flwr_dw_sur    (:,jm:1:-1,:)), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="downward surface longwave radiation for each surface type",units="W/m2",ncid=ncid)
    call nc_write(fnm,"flwr_dw_sur_cs ", sngl(vars%flwr_dw_sur_cs (:,jm:1:-1,:)), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="downward surface longwave clear-sky radiation for each surface type",units="W/m2",ncid=ncid)
    call nc_write(fnm,"flwr_dw_sur_cld", sngl(vars%flwr_dw_sur_cld(:,jm:1:-1,:)), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="downward surface longwave cloudy radiation for each surface type",units="W/m2",ncid=ncid)
    call nc_write(fnm,"flwr_up_sur    ", sngl(vars%flwr_up_sur    (:,jm:1:-1,:)), dims=["lon ","lat ","st  ","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="upward surface longwave radiation for each surface type",units="W/m2",ncid=ncid)

    call nc_write(fnm,"hdust       ", sngl(vars%hdust       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="dust height scale",units="m",ncid=ncid)
    call nc_write(fnm,"dust_load   ", sngl(vars%dust_load   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="atmospheric dust load",units="kg/m2",ncid=ncid)
    call nc_write(fnm,"dust_emis   ", sngl(vars%dust_emis   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="dust emissions",units="kg/m2/s",ncid=ncid)
    call nc_write(fnm,"dust_dep    ", sngl(vars%dust_dep    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="dust deposition",units="kg/m2/s",ncid=ncid)
    call nc_write(fnm,"dust_dep_dry", sngl(vars%dust_dep_dry(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="dry dust deposition",units="kg/m2/s",ncid=ncid)
    call nc_write(fnm,"dust_dep_wet", sngl(vars%dust_dep_wet(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="wet dust deposition",units="kg/m2/s",ncid=ncid)
    call nc_write(fnm,"dust_ot     ", sngl(vars%dust_ot     (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="dust optical thickness",units="1",ncid=ncid)

    if (l_co2d) then
    call nc_write(fnm,"cam         ", sngl(vars%cam         (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="CO2 mass mixing ratio",units="kg(CO2)/kg(air)",ncid=ncid)
    call nc_write(fnm,"co2d        ", sngl(vars%co2d        (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="CO2 concentration",units="ppmv",ncid=ncid)
    call nc_write(fnm,"co2flx      ", sngl(vars%co2flx      (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="surface CO2 flux",units="gCO2/m2/day",ncid=ncid)
    endif

    call nc_write(fnm,"so4_ot      ", sngl(vars%so4_ot      (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="SO4 optical thickness",units="1",ncid=ncid)

    if (l_output_extended) then
    call nc_write(fnm,"dtamdt     ", sngl(vars%dtamdt    (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="tendency of atmospheric temperature",units="K/day",ncid=ncid)
    call nc_write(fnm,"prc_conv   ", sngl(vars%prc_conv  (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="precipitation generated by moisture convergence",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"prc_wcon   ", sngl(vars%prc_wcon  (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="precipitation generated from residence time of water vapor",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"prc_over   ", sngl(vars%prc_over  (:,jm:1:-1) ), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="precipitation generated by crossing saturation level",units="kg/m2/day",ncid=ncid)
    call nc_write(fnm,"dswd_dalb_vu_cs ", sngl(vars%dswd_dalb_vu_cs (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="partial derivative of shortwave surface down radiation wrt surface albedo",units="W/m2",ncid=ncid)
    call nc_write(fnm,"dswd_dalb_ir_cs ", sngl(vars%dswd_dalb_ir_cs (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="partial derivative of shortwave surface down radiation wrt surface albedo",units="W/m2",ncid=ncid)
    call nc_write(fnm,"dswd_dalb_vu_cld", sngl(vars%dswd_dalb_vu_cld(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="partial derivative of shortwave surface down radiation wrt surface albedo",units="W/m2",ncid=ncid)
    call nc_write(fnm,"dswd_dalb_ir_cld", sngl(vars%dswd_dalb_ir_cld(:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="partial derivative of shortwave surface down radiation wrt surface albedo",units="W/m2",ncid=ncid)
    call nc_write(fnm,"dswd_dz_ir_cs   ", sngl(vars%dswd_dz_ir_cs   (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="partial derivative of shortwave surface down radiation wrt surface albedo",units="W/m2",ncid=ncid)
    call nc_write(fnm,"dswd_dz_ir_cld  ", sngl(vars%dswd_dz_ir_cld  (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="partial derivative of shortwave surface down radiation wrt surface albedo",units="W/m2",ncid=ncid)
    if (l_output_flx3d) then
    call nc_write(fnm,"fax       ", sngl(vars%fax       (:,jm:1:-1,:)), dims=["lonu","lat ","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[imc,jm,km,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"faxo      ", sngl(vars%faxo      (:,jm:1:-1,:)), dims=["lonu","lat ","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[imc,jm,km,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"fay       ", sngl(vars%fay       (:,jmc:1:-1,:)), dims=["lon ","latv","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jmc,km,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"fayo      ", sngl(vars%fayo      (:,jmc:1:-1,:)), dims=["lon ","latv","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jmc,km,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"fac       ", sngl(vars%fac       (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
      endif
    call nc_write(fnm,"u3a        ", sngl(vars%ua        (:,jm:1:-1,:)), dims=["lon ","lat ","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="3D ageostrophic zonal wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"v3a        ", sngl(vars%va        (:,jm:1:-1,:)), dims=["lon ","lat ","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="3D ageostrophic meridional wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"uter       ", sngl(vars%uter       (:,jm:1:-1,:)), dims=["lon ","lat ","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="zonal thermal wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"vter       ", sngl(vars%vter       (:,jm:1:-1,:)), dims=["lon ","lat ","zlay","mon ","time"], &
      start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="meridional thermal wind",units="m/s",ncid=ncid)
    call nc_write(fnm,"faxmas    ", sngl(vars%faxmas    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated zonal mass flux",units="",ncid=ncid)
    call nc_write(fnm,"faymas    ", sngl(vars%faymas    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated meridional mass flux",units="",ncid=ncid)
    call nc_write(fnm,"faxdse    ", sngl(vars%faxdse    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated zonal dry static energy flux by mean circulation",units="kg/s*K",ncid=ncid)
    call nc_write(fnm,"faxcpt    ", sngl(vars%faxcpt    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated zonal cp*T flux by mean circulation",units="kg/s*K",ncid=ncid)
    call nc_write(fnm,"faxwtr    ", sngl(vars%faxwtr    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated zonal moisture flux by mean circulation",units="kg/s",ncid=ncid)
    call nc_write(fnm,"faydse    ", sngl(vars%faydse    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated meridional dry static energy flux by mean circulation",units="kg/s*K",ncid=ncid)
    call nc_write(fnm,"faycpt    ", sngl(vars%faycpt    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated meridional cp*T flux by mean circulation",units="kg/s*K",ncid=ncid)
    call nc_write(fnm,"faywtr    ", sngl(vars%faywtr    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated meridional moisture flux by mean circulation",units="kg/s",ncid=ncid)
    call nc_write(fnm,"fdxdse    ", sngl(vars%fdxdse    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated zonal dry static energy flux by synoptic eddies",units="kg/s*K",ncid=ncid)
    call nc_write(fnm,"fdxwtr    ", sngl(vars%fdxwtr    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated zonal moisture flux by synoptic eddies",units="kg/s",ncid=ncid)
    call nc_write(fnm,"fdydse    ", sngl(vars%fdydse    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated meridional dry static energy flux by synoptic eddies",units="kg/s*K",ncid=ncid)
    call nc_write(fnm,"fdywtr    ", sngl(vars%fdywtr    (:,jm:1:-1)), dims=["lon ","lat ","mon ","time"], &
      start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="column integrated meridional moisture flux by synoptic eddies",units="kg/s",ncid=ncid)
    endif

   return

  end subroutine atm_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ d a i l y _ n c
  ! Author   :  Matteo Willeit
  ! Purpose  :  Initialize atmosphere netcdf output
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_daily_nc(fnm)

    implicit none

    character (len=*) :: fnm
    integer :: ncid
    real(wp) :: empty_time(0)

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,"time", x=empty_time, axis="t", units="years BP", &
    unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,"doy", x=1._wp, dx=1._wp, nx=nday_year, axis="e", &
    units="days", ncid=ncid)
    call nc_write_dim(fnm, "lon", x=lon, axis="x", ncid=ncid)
    call nc_write_dim(fnm, "lat", x=lat, axis="y", ncid=ncid)
    call nc_write_dim(fnm, "lev", x=zl(1:km), units="m", ncid=ncid)
    call nc_write_dim(fnm, "lonu",x=-180._wp,dx=5._wp,nx=imc,ncid=ncid)
    call nc_write_dim(fnm, "latv",x=-90._wp,dx=5._wp,nx=jmc,ncid=ncid)
    call nc_write_dim(fnm, "levw",x=zl,units="m",ncid=ncid)
    call nc_write_dim(fnm, "st",x=1._wp,dx=1._wp,nx=nm,ncid=ncid)
    call nc_close(ncid)

   return

  end subroutine atm_daily_nc


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ d a i l y _ n c _ w r i t e
  ! Author   :  Matteo Willeit
  ! Purpose  :  Output of atmosphere netcdf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_daily_nc_write(fnm,ncid,vars,ndat,nout)

    implicit none

    type(a_out) :: vars

    character (len=*) :: fnm
    integer :: ndat, nout, ncid


    call nc_write(fnm,"tam        ", sngl(vars%tam       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"tskina     ", sngl(vars%tskina    (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"qam        ", sngl(vars%qam       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"cam        ", sngl(vars%cam       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"lha        ", sngl(vars%lha       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"sha        ", sngl(vars%sha       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"weff       ", sngl(vars%weff      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"wcld       ", sngl(vars%wcld      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"cld        ", sngl(vars%cld       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"hcld       ", sngl(vars%hcld      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"clot       ", sngl(vars%clot      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"gams       ", sngl(vars%gams      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"gamb       ", sngl(vars%gamb      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"gamt       ", sngl(vars%gamt      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"hrm        ", sngl(vars%hrm       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"prc        ", sngl(vars%prc       (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"wind       ", sngl(vars%wind      (:,jm:1:-1) ), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"convdse    ", sngl(vars%convdse   (:,jm:1:-1)),  dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"convwtr    ", sngl(vars%convwtr   (:,jm:1:-1)),  dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"slp        ", sngl(vars%slp       (:,jm:1:-1)),  dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"us         ", sngl(vars%us        (:,jm:1:-1)),  dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"vs         ", sngl(vars%vs        (:,jm:1:-1)),  dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"ugb        ", sngl(vars%ugb       (:,jm:1:-1)),  dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"vgb        ", sngl(vars%vgb       (:,jm:1:-1)),  dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"uab        ", sngl(vars%uab       (:,jm:1:-1)),  dims=["lonu","lat ","doy ","time"],start=[1,1,ndat,nout],count=[imc,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"vab        ", sngl(vars%vab       (:,jmc:1:-1)), dims=["lon ","latv","doy ","time"],start=[1,1,ndat,nout],count=[im,jmc,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"sam        ", sngl(vars%sam        (:,jm:1:-1)), dims=["lon ","lat ","doy ","time"],start=[1,1,ndat,nout],count=[im,jm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"flwr_dw_sur            ", sngl(vars%flwr_dw_sur            (:,jm:1:-1,:)), dims=["lon ","lat ","st  ","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"flwr_dw_sur_cs         ", sngl(vars%flwr_dw_sur_cs         (:,jm:1:-1,:)), dims=["lon ","lat ","st  ","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="",units="",ncid=ncid)
    call nc_write(fnm,"flwr_dw_sur_cld        ", sngl(vars%flwr_dw_sur_cld        (:,jm:1:-1,:)), dims=["lon ","lat ","st  ","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,nm,1,1],long_name="",units="",ncid=ncid)

!    call nc_write(fnm,"ua        ", sngl(vars%ua        (:,jm:1:-1,:)), dims=["lon ","lat ","lev ","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="",units="",ncid=ncid)
!    call nc_write(fnm,"va        ", sngl(vars%va        (:,jmc:1:-1,:)),dims=["lon ","lat ","lev ","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="",units="",ncid=ncid)
!    call nc_write(fnm,"u3        ", sngl(vars%u3        (:,jm:1:-1,:)), dims=["lon ","lat ","lev ","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="",units="",ncid=ncid)
!    call nc_write(fnm,"v3        ", sngl(vars%v3        (:,jm:1:-1,:)), dims=["lon ","lat ","lev ","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,km,1,1],long_name="",units="",ncid=ncid)
!    call nc_write(fnm,"w3        ", sngl(vars%w3        (:,jm:1:-1,:)), dims=["lon ","lat ","levw","doy ","time"],start=[1,1,1,ndat,nout],count=[im,jm,kmc,1,1],long_name="",units="",ncid=ncid)

   return

  end subroutine atm_daily_nc_write


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ a v e
  ! Purpose  :  Average (or sum) the atmosphere fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_ave(d,ave)

    implicit none

    type(a_out) :: d(:), ave

    integer :: k, n
    real(wp) :: div

    n = size(d)
    div = real(n,wp)

    ! Set all values to zero
    ave%had_fi      = 0._wp
    ave%had_width   = 0._wp
    ave%ptrop       = 0._wp
    ave%tslz        = 0._wp
    ave%tskslz      = 0._wp
    ave%ekez        = 0._wp
    ave%slpz        = 0._wp
    ave%vabz        = 0._wp
    ave%acbarz      = 0._wp
    ave%uz850       = 0._wp
    ave%uz500       = 0._wp
    ave%tam         = 0._wp
    ave%dtamdt      = 0._wp
    ave%qam         = 0._wp
    ave%ram         = 0._wp
    ave%gams        = 0._wp
    ave%gamb        = 0._wp
    ave%gamt        = 0._wp
    ave%dam         = 0._wp
    ave%hqeff       = 0._wp
    ave%hrm         = 0._wp
    ave%wcon        = 0._wp
    ave%cld_rh        = 0._wp
    ave%cld_low        = 0._wp
    ave%cld         = 0._wp
    ave%prc         = 0._wp
    ave%prcw        = 0._wp
    ave%prcs        = 0._wp
    ave%prc_conv    = 0._wp
    ave%prc_wcon    = 0._wp
    ave%prc_over    = 0._wp
    ave%hcld        = 0._wp
    ave%ctt         = 0._wp
    ave%clot        = 0._wp
    ave%alb_cld     = 0._wp
    ave%alb_sur_cs  = 0._wp
    ave%alb_sur_cld = 0._wp
    ave%htrop       = 0._wp
    ave%ttrop       = 0._wp
    ave%wind        = 0._wp
    ave%frst        = 0._wp
    ave%tskin       = 0._wp
    ave%t2          = 0._wp
    ave%q2          = 0._wp
    ave%r2          = 0._wp
    ave%ra2         = 0._wp
    ave%alb_vu_s    = 0._wp
    ave%alb_vu_c    = 0._wp
    ave%alb_ir_s    = 0._wp
    ave%alb_ir_c    = 0._wp
    ave%cd          = 0._wp
    ave%tskina      = 0._wp
    ave%t2a         = 0._wp
    ave%thetae      = 0._wp
    ave%q2a         = 0._wp
    ave%dq          = 0._wp
    ave%dr          = 0._wp
    ave%r2a         = 0._wp
    ave%rskina      = 0._wp
    ave%ra2a        = 0._wp
    ave%cda         = 0._wp
    ave%cd0a        = 0._wp
    ave%sha         = 0._wp
    ave%lha         = 0._wp
    ave%evpa        = 0._wp
    ave%Ri          = 0._wp
    ave%t3          = 0._wp
    ave%q3          = 0._wp
    ave%r3          = 0._wp
    ave%tp          = 0._wp
    ave%gamma       = 0._wp
    ave%rho         = 0._wp
    ave%acbar       = 0._wp
    ave%epsa        = 0._wp
    ave%tsksl        = 0._wp
    ave%atsksl        = 0._wp
    ave%atsl        = 0._wp
    ave%atsli         = 0._wp
    ave%aslp        = 0._wp
    ave%aslp_temp   = 0._wp
    ave%aslp_topo   = 0._wp
    ave%dz500       = 0._wp
    ave%slp         = 0._wp
    ave%slp1        = 0._wp
    ave%us         = 0._wp
    ave%vs         = 0._wp
    ave%usk        = 0._wp
    ave%vsk        = 0._wp
    ave%ugb        = 0._wp
    ave%vgb        = 0._wp
    ave%uab        = 0._wp
    ave%vab        = 0._wp
    ave%taux        = 0._wp
    ave%tauy        = 0._wp
    ave%wcld        = 0._wp
    ave%woro        = 0._wp
    ave%weff        = 0._wp
    ave%fweff       = 0._wp
    ave%ua         = 0._wp
    ave%va         = 0._wp
    ave%u3         = 0._wp
    ave%v3         = 0._wp
    ave%w3         = 0._wp
    ave%uter        = 0._wp
    ave%vter        = 0._wp
    if (l_output_flx3d) then
      ave%fax        = 0._wp
      ave%faxo       = 0._wp
      ave%fay        = 0._wp
      ave%fayo       = 0._wp
      ave%fac        = 0._wp
    endif
    ave%xz          = 0._wp
    ave%diffxdse    = 0._wp
    ave%diffydse    = 0._wp
    ave%diffxwtr    = 0._wp
    ave%diffywtr    = 0._wp
    ave%wsyn        = 0._wp
    ave%convdse     = 0._wp
    ave%convadse    = 0._wp
    ave%convddse    = 0._wp
    ave%convwtr     = 0._wp
    ave%convawtr    = 0._wp
    ave%convdwtr    = 0._wp
    ave%faxmas     = 0._wp
    ave%faymas     = 0._wp
    ave%faxdse     = 0._wp
    ave%faxcpt     = 0._wp
    ave%faxwtr     = 0._wp
    ave%faydse     = 0._wp
    ave%faycpt     = 0._wp
    ave%faywtr     = 0._wp
    ave%fdxdse     = 0._wp
    ave%fdxwtr     = 0._wp
    ave%fdydse     = 0._wp
    ave%fdywtr     = 0._wp
    ave%fayg       = 0._wp
    ave%faydseg    = 0._wp
    ave%faycptg    = 0._wp
    ave%faygzg     = 0._wp
    ave%fayleg     = 0._wp
    ave%faywtrg    = 0._wp
    ave%fdydseg    = 0._wp
    ave%fdyleg     = 0._wp
    ave%fdywtrg    = 0._wp
    ave%fydseg     = 0._wp
    ave%fyheatg    = 0._wp
    ave%fyleg      = 0._wp
    ave%fywtrg     = 0._wp
    ave%fswr_sur    = 0._wp
    ave%flwr_dw_sur = 0._wp
    ave%flwr_dw_sur_cs = 0._wp
    ave%flwr_dw_sur_cld = 0._wp
    ave%flwr_up_sur = 0._wp
    ave%dswd_dalb_vu_cs  = 0._wp
    ave%dswd_dalb_ir_cs  = 0._wp
    ave%dswd_dalb_vu_cld = 0._wp
    ave%dswd_dalb_ir_cld = 0._wp
    ave%dswd_dz_ir_cs  = 0._wp
    ave%dswd_dz_ir_cld = 0._wp
    ave%rb_top      = 0._wp
    ave%rb_sur      = 0._wp
    ave%rb_atm      = 0._wp
    ave%swr_top    = 0._wp
    ave%swr_atm    = 0._wp
    ave%swr_top_cs    = 0._wp
    ave%swr_sur    = 0._wp
    ave%swr_sur_cs    = 0._wp
    ave%swr_dw_sur = 0._wp
    ave%swr_dw_top = 0._wp
    ave%lwr_atm    = 0._wp
    ave%lwr_atm_cld    = 0._wp
    ave%lwr_top    = 0._wp
    ave%lwr_top_cs    = 0._wp
    ave%lwr_sur    = 0._wp
    ave%lwr_sur_cs = 0._wp
    ave%lwr_tro = 0._wp
    ave%lwr_cld = 0._wp
    ave%rb_str      = 0._wp
    ave%cre_top     = 0._wp
    ave%swr_cre_top = 0._wp
    ave%lwr_cre_top = 0._wp
    ave%cre_sur     = 0._wp
    ave%swr_cre_sur = 0._wp
    ave%lwr_cre_sur = 0._wp
    ave%aplan       = 0._wp
    ave%tsl         = 0._wp
    ave%sam         = 0._wp
    ave%sam2        = 0._wp
    ave%eke         = 0._wp
    ave%synprod     = 0._wp
    ave%syndiss     = 0._wp
    ave%synadv      = 0._wp
    ave%syndif      = 0._wp
    ave%synsur      = 0._wp
    ave%cdif        = 0._wp
    ave%fw_pac_atl  = 0.
    ave%fw_atl_indpac = 0.
    ave%hdust         = 0.
    ave%dust_load     = 0.
    ave%dust_emis     = 0.
    ave%dust_dep      = 0.
    ave%dust_dep_dry  = 0.
    ave%dust_dep_wet  = 0.
    ave%dust_ot       = 0.
    ave%cam           = 0.
    ave%co2d          = 0.
    ave%co2flx        = 0.
    ave%so4_ot        = 0.

    ! Loop over the time indices to sum up and average (if necessary)
    do k = 1, n
       ave%had_fi      = ave%had_fi      + d(k)%had_fi      / div
       ave%had_width   = ave%had_width   + d(k)%had_width   / div
       ave%ptrop       = ave%ptrop       + d(k)%ptrop       / div
       ave%tslz        = ave%tslz        + d(k)%tslz        / div
       ave%tskslz       = ave%tskslz       + d(k)%tskslz       / div
       ave%ekez        = ave%ekez        + d(k)%ekez        / div
       ave%slpz        = ave%slpz        + d(k)%slpz        / div
       ave%vabz        = ave%vabz        + d(k)%vabz        / div
       ave%acbarz      = ave%acbarz      + d(k)%acbarz      / div
       ave%uz850       = ave%uz850       + d(k)%uz850       / div
       ave%uz500       = ave%uz500       + d(k)%uz500       / div
       ave%tam         = ave%tam         + d(k)%tam         / div
       ave%dtamdt      = ave%dtamdt      + d(k)%dtamdt      / div
       ave%qam         = ave%qam         + d(k)%qam         / div
       ave%ram         = ave%ram         + d(k)%ram         / div
       ave%gams        = ave%gams        + d(k)%gams        / div
       ave%gamb        = ave%gamb        + d(k)%gamb        / div
       ave%gamt        = ave%gamt        + d(k)%gamt        / div
       ave%dam         = ave%dam         + d(k)%dam         / div
       ave%hqeff       = ave%hqeff       + d(k)%hqeff       / div
       ave%hrm         = ave%hrm         + d(k)%hrm         / div
       ave%wcon        = ave%wcon        + d(k)%wcon        / div
       ave%cld_rh         = ave%cld_rh         + d(k)%cld_rh         / div
       ave%cld_low         = ave%cld_low         + d(k)%cld_low         / div
       ave%cld         = ave%cld         + d(k)%cld         / div
       ave%prc         = ave%prc         + d(k)%prc         / div
       ave%prcw        = ave%prcw        + d(k)%prcw        / div
       ave%prcs        = ave%prcs        + d(k)%prcs        / div
       ave%prc_conv    = ave%prc_conv    + d(k)%prc_conv    / div
       ave%prc_wcon    = ave%prc_wcon    + d(k)%prc_wcon    / div
       ave%prc_over    = ave%prc_over    + d(k)%prc_over    / div
       ave%hcld        = ave%hcld        + d(k)%hcld        / div
       ave%ctt         = ave%ctt         + d(k)%ctt         / div
       ave%clot        = ave%clot        + d(k)%clot        / div
       ave%alb_cld     = ave%alb_cld     + d(k)%alb_cld     / div
       ave%alb_sur_cs  = ave%alb_sur_cs  + d(k)%alb_sur_cs  / div
       ave%alb_sur_cld = ave%alb_sur_cld + d(k)%alb_sur_cld / div
       ave%htrop       = ave%htrop       + d(k)%htrop       / div
       ave%ttrop       = ave%ttrop       + d(k)%ttrop       / div
       ave%wind        = ave%wind        + d(k)%wind        / div
       ave%frst        = ave%frst        + d(k)%frst        / div
       ave%tskin       = ave%tskin       + d(k)%tskin       / div
       ave%t2          = ave%t2          + d(k)%t2          / div
       ave%q2          = ave%q2          + d(k)%q2          / div
       ave%r2          = ave%r2          + d(k)%r2          / div
       ave%ra2         = ave%ra2         + d(k)%ra2         / div
       ave%alb_vu_s    = ave%alb_vu_s    + d(k)%alb_vu_s    / div
       ave%alb_vu_c    = ave%alb_vu_c    + d(k)%alb_vu_c    / div
       ave%alb_ir_s    = ave%alb_ir_s    + d(k)%alb_ir_s    / div
       ave%alb_ir_c    = ave%alb_ir_c    + d(k)%alb_ir_c    / div
       ave%tskina      = ave%tskina      + d(k)%tskina      / div
       ave%t2a         = ave%t2a         + d(k)%t2a         / div
       ave%thetae      = ave%thetae      + d(k)%thetae      / div
       ave%q2a         = ave%q2a         + d(k)%q2a         / div
       ave%dq          = ave%dq          + d(k)%dq          / div
       ave%dr          = ave%dr          + d(k)%dr          / div
       ave%r2a         = ave%r2a         + d(k)%r2a         / div
       ave%rskina      = ave%rskina      + d(k)%rskina      / div
       ave%ra2a        = ave%ra2a        + d(k)%ra2a        / div
       ave%cda         = ave%cda         + d(k)%cda         / div
       ave%cd0a        = ave%cd0a        + d(k)%cd0a        / div
       ave%cd          = ave%cd          + d(k)%cd          / div
       ave%sha         = ave%sha         + d(k)%sha         / div
       ave%lha         = ave%lha         + d(k)%lha         / div
       ave%evpa        = ave%evpa        + d(k)%evpa        / div
       ave%Ri          = ave%Ri          + d(k)%Ri          / div
       ave%t3          = ave%t3          + d(k)%t3          / div
       ave%q3          = ave%q3          + d(k)%q3          / div
       ave%r3          = ave%r3          + d(k)%r3          / div
       ave%tp          = ave%tp          + d(k)%tp          / div
       ave%gamma       = ave%gamma       + d(k)%gamma       / div
       ave%rho         = ave%rho         + d(k)%rho         / div
       ave%acbar       = ave%acbar       + d(k)%acbar       / div
       ave%epsa        = ave%epsa        + d(k)%epsa        / div
       ave%slp         = ave%slp         + d(k)%slp         / div
       ave%slp1        = ave%slp1        + d(k)%slp1        / div
       ave%tsksl        = ave%tsksl        + d(k)%tsksl        / div
       ave%atsksl        = ave%atsksl        + d(k)%atsksl        / div
       ave%atsl        = ave%atsl        + d(k)%atsl        / div
       ave%atsli        = ave%atsli        + d(k)%atsli        / div
       ave%aslp        = ave%aslp        + d(k)%aslp        / div
       ave%aslp_temp   = ave%aslp_temp   + d(k)%aslp_temp   / div
       ave%aslp_topo   = ave%aslp_topo   + d(k)%aslp_topo   / div
       ave%dz500       = ave%dz500       + d(k)%dz500       / div
       ave%us         = ave%us         + d(k)%us         / div
       ave%vs         = ave%vs         + d(k)%vs         / div
       ave%usk        = ave%usk        + d(k)%usk        / div
       ave%vsk        = ave%vsk        + d(k)%vsk        / div
       ave%ugb        = ave%ugb        + d(k)%ugb        / div
       ave%vgb        = ave%vgb        + d(k)%vgb        / div
       ave%uab        = ave%uab        + d(k)%uab        / div
       ave%vab        = ave%vab        + d(k)%vab        / div
       ave%taux        = ave%taux        + d(k)%taux        / div
       ave%tauy        = ave%tauy        + d(k)%tauy        / div
       ave%wcld        = ave%wcld        + d(k)%wcld        / div
       ave%woro        = ave%woro        + d(k)%woro        / div
       ave%weff        = ave%weff        + d(k)%weff        / div
       ave%fweff       = ave%fweff       + d(k)%fweff       / div
       ave%ua         = ave%ua         + d(k)%ua         / div
       ave%va         = ave%va         + d(k)%va         / div
       ave%u3         = ave%u3         + d(k)%u3         / div
       ave%v3         = ave%v3         + d(k)%v3         / div
       ave%w3         = ave%w3         + d(k)%w3         / div
       ave%uter        = ave%uter        + d(k)%uter        / div
       ave%vter        = ave%vter        + d(k)%vter        / div
       if (l_output_flx3d) then
         ave%fax        = ave%fax        + d(k)%fax        / div
         ave%faxo       = ave%faxo       + d(k)%faxo       / div
         ave%fay        = ave%fay        + d(k)%fay        / div
         ave%fayo       = ave%fayo       + d(k)%fayo       / div
         ave%fac        = ave%fac        + d(k)%fac        / div
       endif
       ave%xz          = ave%xz          + d(k)%xz          / div
       ave%diffxdse        = ave%diffxdse        + d(k)%diffxdse        / div
       ave%diffydse        = ave%diffydse        + d(k)%diffydse        / div
       ave%diffxwtr        = ave%diffxwtr        + d(k)%diffxwtr        / div
       ave%diffywtr        = ave%diffywtr        + d(k)%diffywtr        / div
       ave%wsyn        = ave%wsyn        + d(k)%wsyn        / div
       ave%convdse     = ave%convdse     + d(k)%convdse     / div
       ave%convadse    = ave%convadse    + d(k)%convadse    / div
       ave%convddse    = ave%convddse    + d(k)%convddse    / div
       ave%convwtr     = ave%convwtr     + d(k)%convwtr     / div
       ave%convawtr    = ave%convawtr    + d(k)%convawtr    / div
       ave%convdwtr    = ave%convdwtr    + d(k)%convdwtr    / div
       ave%faxmas     = ave%faxmas     + d(k)%faxmas     / div 
       ave%faymas     = ave%faymas     + d(k)%faymas     / div 
       ave%faxdse     = ave%faxdse     + d(k)%faxdse     / div 
       ave%faxcpt     = ave%faxcpt     + d(k)%faxcpt     / div 
       ave%faxwtr     = ave%faxwtr     + d(k)%faxwtr     / div
       ave%faydse     = ave%faydse     + d(k)%faydse     / div
       ave%faycpt     = ave%faycpt     + d(k)%faycpt     / div
       ave%faywtr     = ave%faywtr     + d(k)%faywtr     / div
       ave%fdxdse     = ave%fdxdse     + d(k)%fdxdse     / div
       ave%fdxwtr     = ave%fdxwtr     + d(k)%fdxwtr     / div
       ave%fdydse     = ave%fdydse     + d(k)%fdydse     / div
       ave%fdywtr     = ave%fdywtr     + d(k)%fdywtr     / div
       ave%fayg       = ave%fayg       + d(k)%fayg       / div
       ave%faydseg    = ave%faydseg    + d(k)%faydseg    / div
       ave%faycptg    = ave%faycptg    + d(k)%faycptg    / div
       ave%faygzg     = ave%faygzg     + d(k)%faygzg     / div
       ave%fayleg     = ave%fayleg     + d(k)%fayleg     / div
       ave%faywtrg    = ave%faywtrg    + d(k)%faywtrg    / div
       ave%fdydseg    = ave%fdydseg    + d(k)%fdydseg    / div
       ave%fdyleg     = ave%fdyleg     + d(k)%fdyleg     / div
       ave%fdywtrg    = ave%fdywtrg    + d(k)%fdywtrg    / div
       ave%fydseg     = ave%fydseg     + d(k)%fydseg     / div
       ave%fyheatg    = ave%fyheatg    + d(k)%fyheatg    / div
       ave%fyleg      = ave%fyleg      + d(k)%fyleg      / div
       ave%fywtrg     = ave%fywtrg     + d(k)%fywtrg     / div
       ave%fswr_sur    = ave%fswr_sur    + d(k)%fswr_sur    / div
       ave%flwr_dw_sur = ave%flwr_dw_sur + d(k)%flwr_dw_sur / div
       ave%flwr_dw_sur_cs = ave%flwr_dw_sur_cs + d(k)%flwr_dw_sur_cs / div
       ave%flwr_dw_sur_cld = ave%flwr_dw_sur_cld + d(k)%flwr_dw_sur_cld / div
       ave%flwr_up_sur = ave%flwr_up_sur + d(k)%flwr_up_sur / div
       ave%dswd_dalb_vu_cs  = ave%dswd_dalb_vu_cs  + d(k)%dswd_dalb_vu_cs    / div 
       ave%dswd_dalb_ir_cs  = ave%dswd_dalb_ir_cs  + d(k)%dswd_dalb_ir_cs    / div 
       ave%dswd_dalb_vu_cld = ave%dswd_dalb_vu_cld + d(k)%dswd_dalb_vu_cld   / div 
       ave%dswd_dalb_ir_cld = ave%dswd_dalb_ir_cld + d(k)%dswd_dalb_ir_cld   / div 
       ave%dswd_dz_ir_cs  = ave%dswd_dz_ir_cs  + d(k)%dswd_dz_ir_cs    / div 
       ave%dswd_dz_ir_cld = ave%dswd_dz_ir_cld + d(k)%dswd_dz_ir_cld   / div 
       ave%rb_top      = ave%rb_top      + d(k)%rb_top      / div
       ave%rb_sur      = ave%rb_sur      + d(k)%rb_sur      / div
       ave%rb_atm      = ave%rb_atm      + d(k)%rb_atm      / div
       ave%swr_atm    = ave%swr_atm    + d(k)%swr_atm    / div
       ave%swr_top    = ave%swr_top    + d(k)%swr_top    / div
       ave%swr_top_cs    = ave%swr_top_cs    + d(k)%swr_top_cs    / div
       ave%swr_sur    = ave%swr_sur    + d(k)%swr_sur    / div
       ave%swr_sur_cs    = ave%swr_sur_cs    + d(k)%swr_sur_cs    / div
       ave%swr_dw_sur = ave%swr_dw_sur + d(k)%swr_dw_sur / div
       ave%swr_dw_top = ave%swr_dw_top + d(k)%swr_dw_top / div
       ave%lwr_atm    = ave%lwr_atm    + d(k)%lwr_atm    / div
       ave%lwr_atm_cld    = ave%lwr_atm_cld    + d(k)%lwr_atm_cld    / div
       ave%lwr_top    = ave%lwr_top    + d(k)%lwr_top    / div
       ave%lwr_top_cs    = ave%lwr_top_cs    + d(k)%lwr_top_cs    / div
       ave%lwr_sur    = ave%lwr_sur    + d(k)%lwr_sur    / div
       ave%lwr_sur_cs = ave%lwr_sur_cs + d(k)%lwr_sur_cs / div
       ave%lwr_tro = ave%lwr_tro + d(k)%lwr_tro / div
       ave%lwr_cld = ave%lwr_cld + d(k)%lwr_cld / div
       ave%rb_str      = ave%rb_str      + d(k)%rb_str      / div
       ave%cre_top     = ave%cre_top     + d(k)%cre_top     / div
       ave%swr_cre_top = ave%swr_cre_top + d(k)%swr_cre_top / div
       ave%lwr_cre_top = ave%lwr_cre_top + d(k)%lwr_cre_top / div
       ave%cre_sur     = ave%cre_sur     + d(k)%cre_sur     / div
       ave%swr_cre_sur = ave%swr_cre_sur + d(k)%swr_cre_sur / div
       ave%lwr_cre_sur = ave%lwr_cre_sur + d(k)%lwr_cre_sur / div
       ave%aplan       = ave%aplan       + d(k)%aplan       / div
       ave%tsl         = ave%tsl         + d(k)%tsl         / div
       ave%eke         = ave%eke         + d(k)%eke         / div
       ave%sam         = ave%sam         + d(k)%sam         / div
       ave%sam2        = ave%sam2        + d(k)%sam2        / div
       ave%synprod     = ave%synprod     + d(k)%synprod     / div
       ave%syndiss     = ave%syndiss     + d(k)%syndiss     / div
       ave%synadv      = ave%synadv      + d(k)%synadv      / div
       ave%syndif      = ave%syndif      + d(k)%syndif      / div
       ave%synsur      = ave%synsur      + d(k)%synsur      / div
       ave%cdif        = ave%cdif        + d(k)%cdif        / div
       ave%fw_pac_atl  = ave%fw_pac_atl  + d(k)%fw_pac_atl  / div
       ave%fw_atl_indpac = ave%fw_atl_indpac + d(k)%fw_atl_indpac / div
       ave%hdust        = ave%hdust        + d(k)%hdust        / div
       ave%dust_load    = ave%dust_load    + d(k)%dust_load    / div
       ave%dust_emis    = ave%dust_emis    + d(k)%dust_emis    / div
       ave%dust_dep     = ave%dust_dep     + d(k)%dust_dep     / div
       ave%dust_dep_dry = ave%dust_dep_dry + d(k)%dust_dep_dry / div
       ave%dust_dep_wet = ave%dust_dep_wet + d(k)%dust_dep_wet / div
       ave%dust_ot      = ave%dust_ot      + d(k)%dust_ot      / div
       ave%cam          = ave%cam     + d(k)%cam    / div
       ave%co2d         = ave%co2d    + d(k)%co2d   / div
       ave%co2flx       = ave%co2flx  + d(k)%co2flx / div
       ave%so4_ot       = ave%so4_ot       + d(k)%so4_ot       / div
    enddo

   return

  end subroutine atm_ave

end module atm_out

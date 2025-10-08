!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : p h o t o s y n t h e s i s _ m o d
!
!  Purpose : photosynthesis
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
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
module photosynthesis_mod

  use precision, only : wp
  use timer, only : sec_day, day_year, time_soy_lnd, time_eoy_lnd
  use constants, only : e_sat_w, q_to_e, T0
  use lnd_grid, only : nsurf, npft, nveg, nl, flag_tree, flag_veg
  use lnd_params, only : i_ci, i_beta, i_vcmax, i_disc, veg_par, pft_par, sw_par_frac

  real(wp), parameter :: alpha_c3 = 0.08_wp
  real(wp), parameter :: alpha_c4 = 0.053_wp
  real(wp), parameter :: a_c3 = 0.015_wp      ! leaf respiration as fraction of Vmax for C3 plants 
  real(wp), parameter :: a_c4 = 0.02_wp       ! leaf respiration as fraction of vmax for C4 plants 
  real(wp), parameter :: ne_c3 = 8.e-4_wp     ! molCO2 m−2 s−1 kg C (kg N)−1
  real(wp), parameter :: ne_c4 = 4.e-4_wp     ! molCO2 m−2 s−1 kg C (kg N)−1

  real(wp), parameter :: O2 = 20.9_wp      ! O2 concentration in air (%)

  real(wp), parameter :: thetar = 0.7_wp 

  real(wp), parameter :: q10ko  = 1.2_wp    ! q10 for temperature-sensitive parameter ko 
  real(wp), parameter :: q10kc  = 2.1_wp    ! q10 for temperature-sensitive parameter kc 
  real(wp), parameter :: q10tau = 0.57_wp   ! q10 for temperature-sensitive parameter tau 
  real(wp), parameter :: ko25   = 3.0e4_wp  ! value of ko at 25 deg C 
  real(wp), parameter :: kc25   = 30.0_wp   ! value of kc at 25 deg C 
  real(wp), parameter :: tau25  = 2600.0_wp ! value of tau at 25 deg C 
  real(wp), parameter :: eta25  = 9.1e-4_wp ! value of eta at 25 deg C

  real(wp), parameter :: cmass = 12.0_wp    ! atomic mass of carbon 
  real(wp), parameter :: cq    = 4.6e-6_wp  ! mol/J, conversion factor for solar radiation at 550 nm,  J/m2 --> (mol photons)/m2 

  real(wp), parameter :: alphaa = 0.4_wp  ! fraction of PAR assimilated at ecosystem level, relative to leaf level

  real(wp), parameter :: cn_sapwood = 330._wp
  real(wp), parameter :: cn_root = 29._wp

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p h o t o s y n t h e s i s
  !   Purpose    :  computes daily net carbon assimilation for each PFT
  !                 after Haxeltine & Prentice, 1996. Used also in LPJ
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine photosynthesis(co2,c13_c12_atm,c14_c_atm, &
                           frac_surf,t2m,t2m_min_mon,gdd5,t_soil,q2m,pressure,swnet,albedo,daylength, &
                           theta_w,theta_field,theta_wilt,wilt,root_frac, &
                           lai,phen,leaf_c,stem_c,root_c, &
                           discrimination,ci,g_can,gpp,npp,npp13,npp14,npp_cum,npp13_cum,npp14_cum,npp_ann,npp13_ann,npp14_ann, &
                           aresp)

    implicit none

    real(wp), intent(in) :: co2, c13_c12_atm, c14_c_atm
    real(wp), dimension(:), intent(in) :: frac_surf, pressure, t2m, q2m, swnet, albedo
    real(wp), intent(in) :: t2m_min_mon, gdd5, daylength
    real(wp), dimension(:), intent(in) :: theta_w, theta_field, theta_wilt, t_soil
    real(wp), dimension(:), intent(in) :: lai, phen
    real(wp), dimension(:), intent(in) :: leaf_c
    real(wp), dimension(:), intent(in) :: stem_c
    real(wp), dimension(:), intent(in) :: root_c
    real(wp), dimension(:,:), intent(out) :: wilt
    real(wp), dimension(:,:), intent(in) :: root_frac
    real(wp), dimension(:), intent(out) :: discrimination, ci, g_can, gpp, npp, npp13, npp14
    real(wp), dimension(:), intent(out) :: aresp
    real(wp), dimension(:), intent(inout) :: npp_cum, npp13_cum, npp14_cum
    real(wp), dimension(:), intent(inout) :: npp_ann, npp13_ann, npp14_ann

    integer :: n, k
    logical :: lim_bio
    real(wp) :: tleaf, vpd, sqrt_vpd, par, apar, temp_air, hum_air, swdown, f_veg, p
    real(wp) :: pO2, tau, kc, ko, gamma_star, eta_star, fac, fexp
    real(wp) :: beta_theta, dtheta, zeta, xi 
    real(wp) :: c_1, c_2, sigma, ca, pa, pi, b, s
    real(wp) :: vm, vm25, je, jc, tstress, beta
    real(wp) :: agd, rd, rm, and, adt, conv_fac
    real(wp) :: low, high
    real(wp) :: ftemp_air, ftemp_soil, resp10, wT
    real(wp) :: resp_leaf, resp_stem, resp_root
    real(wp) :: Ratm13
    real(wp) :: Ratm14
    real(wp) :: alpha13_gpp, alpha14_gpp
    real(wp), parameter :: eps = 1.e-30_wp


    ! compute mean climate over vegetated grid cell part
    f_veg = sum(frac_surf,mask=flag_veg.eq.1)
    p  = 0._wp
    temp_air = 0._wp
    hum_air = 0._wp
    swdown = 0._wp
    do n=1,nveg
      p = p + pressure(n)*frac_surf(n)/f_veg
      temp_air = temp_air + t2m(n)*frac_surf(n)/f_veg
      hum_air = hum_air + q2m(n)*frac_surf(n)/f_veg
      swdown = swdown + swnet(n)/(1._wp-albedo(n))*frac_surf(n)/f_veg
    enddo

    if (daylength.gt.0._wp .and. (temp_air-T0).gt.-3._wp) then

      ! air temperature (degC)
      tleaf = temp_air - T0
      ! vapor pressure deficit (Pa)
      vpd = max(e_sat_w(temp_air)-q_to_e(hum_air,p), 0.1_wp)
      sqrt_vpd = sqrt(vpd)

      tau = tau25 * q10tau**((tleaf-25._wp)*0.1_wp)
      kc  = kc25  * q10kc**((tleaf-25._wp)*0.1_wp)
      ko  = ko25  * q10ko**((tleaf-25._wp)*0.1_wp)
      eta_star = 1e-3_wp*exp(-3.719_wp+580._wp/(-138._wp+tleaf+T0))/eta25   ! water viscosity relative to reference at 25 degC, eq. 13 in Wang 2017

      ! elevation (pressure) dependence of O2 partial pressure
      pO2 = O2*1.e-2_wp * p  ! Pa

      gamma_star = pO2 / ( 2._wp * tau )    ! Pa
      fac = 1._wp + pO2 / ko

      ! temperature factors for autotrophic respiration
      if( t_soil(2) .gt. 240._wp ) then
        ftemp_soil = exp(308.56_wp * (1._wp/56.02_wp - 1._wp/(46.02_wp+t_soil(2)-T0)) )
      else
        ftemp_soil = 0._wp
      endif
      ftemp_air = exp(308.56_wp * (1._wp/56.02_wp - 1._wp/(46.02_wp+temp_air-T0)) )

      ! NET photosynthetically active radiation, in mol/m2/day
      par = swdown*sw_par_frac * sec_day * (1._wp-veg_par%leaf_alb) * cq 

      do n=1,npft

        ! bioclimatic limit, 10 degC tolerance
        lim_bio = (t2m_min_mon-T0).ge.(pft_par%t_cmon_min(n)-10._wp) .and. (t2m_min_mon-T0).le.(pft_par%t_cmon_max(n)+10._wp) &
          .and. gdd5.ge.(pft_par%gdd5_min(n)-200._wp)

        if(lai(n).gt.0._wp .and. lim_bio) then

          fexp = (1._wp-exp(-veg_par%ext_coef*lai(n)))
          ! absorbed PAR limited by the fraction of par assimilated at ecosystem level, the leaf scattering
          apar = alphaa * fexp * par

          ! temperature stress factor
          if( tleaf.lt.pft_par%t_co2_high(n) .and.  tleaf.gt.pft_par%t_co2_low(n)) then
            low = 1._wp / (1._wp + exp(pft_par%photos_k1(n) * (pft_par%photos_k2(n) - tleaf)))
            high = 1._wp - 0.01_wp * exp(pft_par%photos_k3(n) * (tleaf - pft_par%t_photos_high(n)))
            tstress = low * high
          else
            tstress = 0._wp
          endif

          ! air CO2 concentration
          ca = co2 * 1.e-6_wp  ! mol/mol
          ! air CO2 partial pressure
          pa = ca*p  ! Pa

          ! C3 path
          if( pft_par%i_c4(n) .eq. 0 ) then

            ! ratio of leaf-internal and air CO2 concentration
            if (i_ci.eq.1) then
              ! Cowan–Farquhar optimality hypothesis (Medlyn 2011)
              xi = 1._wp - 1.6_wp / (1._wp + pft_par%g1(n) / (sqrt_vpd*sqrt(1.e-3_wp)))
            else if (i_ci.eq.2) then
              ! least cost optimality model (Prentice 2014)
              if (i_beta.eq.0) then
                ! constant beta
                beta_theta = pft_par%beta(n)
              else if (i_beta.eq.1) then
                ! soil moisture dependent beta following Lavergne 2020, Fig. 6
                dtheta = sum((theta_w(:)-theta_field(:))*root_frac(:,n))
                beta_theta = pft_par%beta(n)*exp(1.81_wp*dtheta)   ! pft_par%beta(n) is beta at field capacity
              endif
              ! zeta parameter for ci/ca ratio, Lavergne 2019 eq. 4
              zeta = sqrt(beta_theta*(kc*fac+gamma_star)/(1.6_wp*eta_star))
              xi = gamma_star/pa + (1._wp-gamma_star/pa)*zeta/(zeta+sqrt_vpd)
            endif

            ! leaf internal CO2 concentration
            ci(n) = xi * ca  ! mol/mol
            ! leaf internal CO2 partial pressure
            pi = ci(n)*p  ! Pa

            c_1 = alpha_c3 * tstress * cmass * (pi - gamma_star) / (pi + 2._wp * gamma_star) 
            c_2 = (pi - gamma_star) / (pi + kc * fac)

            if (i_vcmax.eq.1) then
              ! 'strong optimality hypothesis', Haxeltine and Prentice 1996
              s = 24._wp / daylength * a_c3
              sigma = 1._wp - (c_2 - s) / (c_2 - thetar * s)
              if( sigma .gt. 0._wp ) then 
                sigma = sqrt( sigma )
              else
                sigma = 0._wp
              endif
              vm = (1.0_wp/a_c3)*(c_1/c_2)*((2.0_wp*thetar-1.0_wp)*s-(2.0_wp*thetar*s-c_2)*sigma)*apar  ! gC/m2/d
            else if (i_vcmax.eq.2) then
              ! CLM 4.5, eq. 8.17
              ! Vc,max at 25°C, molCO2/m2/s
              vm25 = 1._wp/(25._wp*pft_par%sla(n)*1.e-3_wp)*pft_par%flnr(n)*7.16_wp*60._wp *1.e-6_wp
              vm = tstress*vm25*2._wp**(0.1_wp*(tleaf-25._wp))
              ! convert to gC/m2/d
              vm = vm * sec_day*cmass
            else if (i_vcmax.eq.3) then
              ! coordination hypothesis (acclimation), see e.g. Harrison 2021 Box 2
              !vm = 1.02_wp*apar*(pi+kc*fac)/(pi+2._wp*gamma_star)  ! gC/mol * mol/m2/day = gC/m2/day
              vm = c_1*apar*(pi+kc*fac)/(pi-gamma_star)  ! gC/mol * mol/m2/day = gC/m2/day
            endif

            b = a_c3

            ! C4 path
          else

            ! ratio of leaf-internal and air CO2 concentration
            ! Cowan–Farquhar optimality hypothesis (Medlyn 2011)
            xi = 1._wp - 1.6_wp / (1._wp + pft_par%g1(n) / (sqrt_vpd*sqrt(1.e-3_wp)))

            ! leaf internal CO2 concentration
            ci(n) = xi * ca  ! mol/mol
            ! leaf internal CO2 partial pressure
            pi = ci(n)*p  ! Pa

            c_1 = alpha_c4 * tstress * cmass 
            c_2 = 1.0_wp

            if (i_vcmax.eq.1) then
              ! 'strong optimality hypothesis', Haxeltine and Prentice 1996
              s = 24._wp / daylength * a_c4
              sigma = 1._wp - (c_2 - s) / (c_2 - thetar * s)
              if( sigma .gt. 0._wp ) then 
                sigma = sqrt( sigma )
              else
                sigma = 0._wp
              endif
              vm = (1.0_wp/a_c4)*c_1/c_2*((2.0_wp*thetar-1.0_wp)*s-(2.0_wp*thetar*s-c_2)*sigma)*apar
            else if (i_vcmax.eq.2) then
              ! CLM 4.5, eq. 8.17
              ! Vc,max at 25°C, molCO2/m2/s
              vm25 = 1._wp/(25._wp*pft_par%sla(n)*1.e-3_wp)*pft_par%flnr(n)*7.16_wp*60._wp *1.e-6_wp
              vm = tstress*vm25*2._wp**(0.1_wp*(tleaf-25._wp))
              ! convert to gC/m2/d
              vm = vm * sec_day*cmass
            else if (i_vcmax.eq.3) then
              ! coordination hypothesis (acclimation), see e.g. Harrison 2021 Box 2
              !vm = 1.02_wp*apar*(pi+kc*fac)/(pi+2._wp*gamma_star)  ! gC/mol * mol/m2/day = gC/m2/day
              vm = c_1*apar  ! gC/mol * mol/m2/day = gC/m2/day
            endif

            b = a_c4

          endif

          ! Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
          ! Eqn 3, Haxeltine & Prentice 1996

          je = c_1*apar/daylength

          ! Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
          ! Eqn 5, Haxeltine & Prentice 1996

          jc = c_2*vm/24._wp

          ! limit by soil moisture, beta factor, weighted with root depth
          beta = 0._wp
          do k=1,nl
            if (theta_w(k).ge.theta_field(k)) then
              wilt(k,n) = 1._wp
            else if (theta_w(k).le.theta_wilt(k)) then
              wilt(k,n) = 0._wp
            else
              wilt(k,n) = (theta_w(k) - theta_wilt(k))/(theta_field(k) - theta_wilt(k))
            endif
            beta = beta + wilt(k,n) * root_frac(k,n)
          enddo

          ! Calculation of daily gross photosynthesis, Agd, gC/m2/day
          ! Eqn 2, Haxeltine & Prentice 1996
          agd = (je+jc-sqrt((je+jc)*(je+jc)-4.0_wp*thetar*je*jc))/(2.0_wp*thetar)*daylength &  ! leaf level, pot
            * beta &         ! soil moisture limitation
            * 1._wp/(1._wp+veg_par%gamma_down*log(co2/280._wp))    ! photosyntheis downregulation following Arora & Scinocca 2016, eq. 11/12

          agd = max(0._wp, agd)

          ! Daily leaf respiration, Rd, gC/m2/day
          ! Eqn 10, Haxeltine & Prentice 1996    
          rd = b*vm         &  ! leaf level, pot
            * beta           ! soil moisture limitation

          ! Daily net photosynthesis, And, gC/m2/day
          and = agd - rd

          if( and .lt. 0._wp ) then
            and = 0._wp
            rd = 0._wp
          endif

          ! Total DAYTIME net photosynthesis, Adt, gC/m2/day
          ! Eqn 19, Haxeltine & Prentice 1996
          adt = and + (1._wp-daylength/24._wp) * rd

          ! daily CANOPY GPP, kgC/m2/day
          gpp(n) = agd * 1.e-3_wp  ! convert from gC to kgC

          ! canopy conductance conversion factor, from mol/m2/s to m/s using ideal gas equation
          conv_fac = 8.314_wp * (tleaf+T0) / p  ! m3/mol

          ! average DAYTIME CANOPY conductance for CO2, m/s
          g_can(n) = pft_par%g_min(n) * fexp * beta * 1.e-3_wp & ! minimum canopy conductance, m/s, depends on LAI, limited by soil moisture
            + adt/cmass / ((ca - ci(n)) * 3600._wp*daylength) * conv_fac ! gC/m2/day * mol/gC / s/h / h/day * m3/mol

          if( flag_tree(n) .eq. 1) then
            wT = min(max(0._wp,(t2m_min_mon-T0)),20._wp)/20._wp
            resp10 = wT*0.011_wp + (1._wp-wT)*0.066_wp ! tropical trees
          else
            resp10 = 0.066_wp 
          endif

          ! leaf respiration
          resp_leaf = rd/1000._wp   ! converted from gC to kgC
          ! stem respiration
          resp_stem = resp10*stem_c(n)/pft_par%aws(n)/cn_sapwood*ftemp_air
          ! root respiration
          resp_root = resp10*root_c(n)/cn_root*ftemp_soil*phen(n)

          rm = resp_leaf + resp_stem + resp_root

          ! growth respiration, accounted for implicitely in the computation of NPP
          ! rg = 0.25 * NPP

          ! net primary productivity, kgC/m2/day, LPJ
          npp(n) = 0.75_wp * (gpp(n) - rm)
          npp(n) = max(0._wp, npp(n))  ! limit NPP to be positive, check


          !---------------------------
          ! isotopes

          Ratm13  = c13_c12_atm
          Ratm14  = c14_c_atm

          ! fractionation during assimilation, Farquar 1994, Lloyd 2004
          ! C3 path
          if( pft_par%i_c4(n) .eq. 0 ) then
            ! C3 path
            if (i_disc.eq.1) then
              discrimination(n) = 4.4_wp*(ca-ci(n))/ca + 27._wp*ci(n)/ca ! permil
            else if (i_disc.eq.2) then
              !including an explicit fractionation term for photorespiration as
              !recommended by several studies (e.g. Ubierna&Farquhar, 2014;
              !Schubert & Jahren, 2018; Lavergne et al., 2019)
              discrimination(n) = 4.4_wp*(ca-ci(n))/ca + 27._wp*ci(n)/ca - 12._wp*gamma_star/pa ! permil
            endif
          else
            ! C4 path
            discrimination(n) = 4.4_wp*(ca-ci(n))/ca + (-5.7_wp+20._wp*0.35_wp)*ci(n)/ca ! permil
          endif

          alpha13_gpp = 1._wp-discrimination(n)*1.e-3_wp
          alpha14_gpp = 1._wp-2._wp*discrimination(n)*1.e-3_wp

          npp13(n) = alpha13_gpp*Ratm13*npp(n)
          npp14(n) = alpha14_gpp*Ratm14*npp(n)

        else ! no photosynthesis

          wilt(:,n)      = 0._wp
          discrimination(n) = 0._wp ! check
          g_can(n)       = 0._wp
          gpp(n)         = 0._wp
          npp(n)         = 0._wp
          npp13(n)       = 0._wp
          npp14(n)       = 0._wp

        endif

      enddo


    else ! no light

      wilt      = 0._wp
      discrimination = 0._wp ! check
      g_can       = 0._wp
      gpp         = 0._wp
      npp         = 0._wp
      npp13       = 0._wp
      npp14       = 0._wp

    endif

    ! autotrophic respiration, kgC/m2/day
    aresp = gpp - npp

    ! compute annual mean NPP for vegetation carbon and dynamics
    if (time_soy_lnd) then
      npp_cum   = 0._wp
      npp13_cum = 0._wp
      npp14_cum = 0._wp
    endif
    npp_cum   = npp_cum   + npp   / sec_day / day_year ! kgC/m2/s
    npp13_cum = npp13_cum + npp13 / sec_day / day_year ! kgC/m2/s
    npp14_cum = npp14_cum + npp14 / sec_day / day_year ! kgC/m2/s
    if (time_eoy_lnd) then
      npp_ann   = npp_cum
      npp13_ann = npp13_cum
      npp14_ann = npp14_cum
    endif


    return

  end subroutine photosynthesis

end module photosynthesis_mod



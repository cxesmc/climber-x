!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m b _ s u r f a c e _ p a r 
!
!  Purpose : surface parameters for SEMIX model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit, Andrey Ganopolski and Reinhard Calov
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
module smb_surface_par

  use precision, only : wp
  use constants, only : karman, g, pi, T0, rho_i, z_sfl
  use smb_params, only : dt, l_neutral
  use smb_params, only : snow_par, surf_par

  implicit none

  private
  public :: surface_albedo, resistance

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u r f a c e _ a l b e d o
  !   Purpose    :  compute surface albedo 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surface_albedo(z_sur_std, f_ice, alb_ice, h_snow, w_snow, w_snow_max, snow, &
                        t_skin, dust_dep, coszm, &
                        snow_grain, dust_con, f_snow, dt_snowfree, alb_bg, &
                        alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif, &
                        alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif)

    implicit none

    real(wp), intent(in) :: z_sur_std
    real(wp), intent(in) :: f_ice
    real(wp), intent(in) :: alb_ice
    real(wp), intent(in) :: h_snow, w_snow, w_snow_max
    real(wp), intent(in) :: snow
    real(wp), intent(in) :: t_skin, dust_dep, coszm
    real(wp), intent(inout) :: snow_grain, dust_con
    real(wp), intent(inout) :: f_snow
    real(wp), intent(inout) :: dt_snowfree
    real(wp), intent(inout) :: alb_bg
    real(wp), intent(inout) :: alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif
    real(wp), intent(inout) :: alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif

    real(wp) :: z0m, f_snow_orog
    real(wp), parameter :: eps = 1.e-10_wp


    if (h_snow.gt.0._wp) then

      ! snow grain size
      if (snow_par%lsnow_aging) then
        call snow_grain_size(t_skin,snow, &
                             snow_grain)
      else
        snow_grain = snow_par%snow_grain_fresh
      endif

      ! dust effect on snow albedo
      if (snow_par%lsnow_dust) then
        ! dust concentration in top snow layer
        call dust_in_snow(dust_dep,snow,w_snow,w_snow_max, &
                          dust_con)
      else
        dust_con = 0._wp
      endif

      ! compute snow albedo
      if (snow_par%isnow_albedo.eq.1) then
        ! climber-2 snow albedo parameterisation, following Warren & Wiscombe 1980
        call snow_albedo_ww(z_sur_std,snow_grain, dust_con, coszm, &
                            alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)
      else if (snow_par%isnow_albedo.eq.2) then
        ! snow albedo parameterisation following Dang et al 2015
        call snow_albedo_dang(z_sur_std, snow_grain, dust_con, coszm, &
                              alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)
      endif

    else

      snow_grain = snow_par%snow_grain_old

    endif

    ! subgrid snow fraction
    ! orography factor for snow cover fraction 
    if (snow_par%l_fsnow_orog) then
      ! reduce snow cover fraction over rough topography following Roesch 2001, eq. 7, but without sqrt
      f_snow_orog = h_snow/(h_snow+snow_par%c_fsnow_orog*z_sur_std+eps)
    else
      f_snow_orog = 1._wp
    endif
    ! snow fraction after Niu and Yang 2007, Roesch 2001
    z0m = surf_par%z0m_ice
    f_snow = tanh(h_snow/(snow_par%c_fsnow*z0m)) * f_snow_orog

    ! compute time over which snow-free
    if (w_snow .lt. snow_par%w_snow_crit) then
      dt_snowfree = dt_snowfree + dt
    endif

    ! background albedo
    ! weight ice and bare soil according to sub-grid ice fraction
    alb_bg = f_ice*alb_ice + (1._wp-f_ice)*surf_par%alb_soil

    ! allband surface albedo for direct and diffuse radiation 
    alb_vis_dir = f_snow * alb_snow_vis_dir + (1._wp-f_snow) * alb_bg
    alb_nir_dir = f_snow * alb_snow_nir_dir + (1._wp-f_snow) * alb_bg
    alb_vis_dif = f_snow * alb_snow_vis_dif + (1._wp-f_snow) * alb_bg
    alb_nir_dif = f_snow * alb_snow_nir_dif + (1._wp-f_snow) * alb_bg


    return

    end subroutine surface_albedo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ g r a i n _ s i z e
  !   Purpose    :  compute snow grain size
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_grain_size(t_skin,snow, &
                             snow_grain)

    implicit none

    real(wp), intent(in) :: t_skin, snow
    real(wp), intent(inout) :: snow_grain

    real(wp) :: f_age, f_tage1, f_tage2, f_tage, f_p


    ! CLIMBER-2 parameterisation of snow grain size (age), tuned to MARv3.6 (using
    ! CROCUS snow model) for Greenland

    ! dry snow temperature dependence for snow grain size
    f_tage1 = exp( snow_par%f_age_t * min(0._wp, (t_skin-(T0-snow_par%dT_age))) )
    ! melting snow temperature dependence for snow grain size
    f_tage2 = exp( min(0._wp, t_skin-(T0-snow_par%dT_age)) )
    ! snow grain size temperature factor
    f_tage = f_tage1 + f_tage2
    ! averaged 'snow age'       
    f_p = f_tage * (snow_par%snow_0/max(1.e-20_wp,snow))**snow_par%snow_1
    f_age = 1._wp - log(1._wp+f_p)/f_p
    ! snow grain size
    snow_grain = snow_par%snow_grain_fresh + (snow_par%snow_grain_old-snow_par%snow_grain_fresh)*f_age 


    return

  end subroutine snow_grain_size


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d u s t _ i n _ s n o w
  !   Purpose    :  compute dust concentration in top snow layer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dust_in_snow(dust_dep,snow,w_snow,w_snow_max, &
                          dust_con)

    implicit none

    real(wp), intent(in) :: dust_dep
    real(wp), intent(in) :: snow
    real(wp), intent(in) :: w_snow, w_snow_max
    real(wp), intent(inout) :: dust_con

    real(wp) :: dust_con_melt_fac

    real(wp), parameter :: dust_con_max = 1000._wp*1.e-6_wp ! kg/kg


    ! compute dust concentration in snowfall
    dust_con = dust_dep/max(1.e-7_wp,snow) ! kg(dust)/m2/s * kg(snow)/m2/s = kg(dust)/kg(snow)
    ! increase dust concentration when snow melts, assuming dust is not removed with melted water
    ! scavenging efficiency of dust with meltwater is 10-30% (Doherty 2013)
    if (w_snow.gt.1._wp) then 
      dust_con_melt_fac = 1._wp + (w_snow_max-w_snow)/snow_par%w_snow_dust  ! melt of w_snow_dust kg/m2 swe causes a doubling of the dust concentration
      dust_con_melt_fac = max(dust_con_melt_fac,1._wp)
      dust_con_melt_fac = min(dust_con_melt_fac,5._wp)  ! limit to a five-fold increase
    else
      dust_con_melt_fac = 1._wp
    endif
    dust_con = dust_con_melt_fac*dust_con

    ! scaling of dust concentration to mimic different imaginary refractive indexes of dust
    dust_con = snow_par%dust_con_scale*dust_con

    ! avoid underflow and negative values
    if (dust_con.lt.1.e-15_wp) dust_con = 0._wp
    ! limit dust concentration in snow to dust_con_max
    dust_con = min(dust_con,dust_con_max)


    return

  end subroutine dust_in_snow


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o _ w w
  !   Purpose    :  compute snow albedo following Warren & Wiscombe 1980
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo_ww(z_sur_std, snow_grain, dust_con, coszm, &
                            alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)

    implicit none

    real(wp), intent(in) :: z_sur_std
    real(wp), intent(in) :: snow_grain, dust_con, coszm
    real(wp), intent(out) :: alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif

    integer :: k
    real(wp) :: f_cosz, f_age
    real(wp) :: d1, rint, c_dust_new_vis, c_dust_age_vis, c_dust_new_nir, c_dust_age_nir
    real(wp) :: c_age_vis, c_age_nir

    real(wp), dimension(4) :: tab0 = (/1.001_wp,10._wp,100._wp,1000._wp/)
    real(wp), dimension(4) :: tab1 = (/0.00_wp, 0.02_wp,0.10_wp,0.30_wp/)
    real(wp), dimension(4) :: tab2 = (/0.01_wp, 0.05_wp,0.15_wp,0.30_wp/)


    ! Clear sky snow albedo, zenit angle dependence
    ! zenith angle factor slightly modified from BATS
    ! 2 in the denominator instead of 4 and applied also to angles < 60
    ! => better agreement with Gardner & Sharp 2010
    f_cosz = 0.5_wp*(3._wp/(1._wp+2._wp*coszm)-1._wp)
    f_cosz = max(0._wp,f_cosz)

    ! compute effect of dust concentration on snow albedo after Warren and Wiscombe 1980, Figure 5
    d1 = min(dust_con*1.d6,999._wp)
    if (d1.gt.1.0001_wp) then
      if (d1.gt.1000._wp) d1=1000._wp
      if (d1.ge.1.0001_wp .and. d1.lt.10._wp) k=1
      if (d1.ge.10._wp .and. d1.lt.100._wp)   k=2
      if (d1.ge.100._wp .and. d1.le.1000._wp) k=3
      rint = (log(d1)-log(tab0(k)))/(log(tab0(k+1))-log(tab0(k)))
      c_dust_new_vis = (1._wp-rint)*tab1(k)+rint*(tab1(k+1))
      c_dust_age_vis = (1._wp-rint)*tab2(k)+rint*(tab2(k+1))
    else
      c_dust_new_vis = 0._wp
      c_dust_age_vis = 0._wp
    endif
    ! NIR reduction of snow albedo by dust largely overestimated by this scheme!
    c_dust_new_nir = 0.5_wp*c_dust_new_vis
    c_dust_age_nir = 0.5_wp*c_dust_age_vis
    c_age_vis = surf_par%d_alb_age_vis + c_dust_age_vis
    c_age_nir = surf_par%d_alb_age_nir + c_dust_age_nir

    !f_age = (snow_grain-snow_par%snow_grain_fresh)/(snow_par%snow_grain_old-snow_par%snow_grain_fresh)
    f_age = log10(1._wp+(snow_grain-snow_par%snow_grain_fresh)/200._wp)/log10(1._wp+(snow_par%snow_grain_old-snow_par%snow_grain_fresh)/200._wp)
    alb_snow_vis_dif = surf_par%alb_vis_dif_snow_new - f_age*c_age_vis - c_dust_new_vis
    alb_snow_nir_dif = surf_par%alb_nir_dif_snow_new - f_age*c_age_nir - c_dust_new_nir
    alb_snow_vis_dir = alb_snow_vis_dif + 0.4_wp*f_cosz*(1._wp-alb_snow_vis_dif)
    alb_snow_nir_dir = alb_snow_nir_dif + 0.4_wp*f_cosz*(1._wp-alb_snow_nir_dif)


    return

  end subroutine snow_albedo_ww


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o _ d a n g
  !   Purpose    :  compute snow albedo following Dang et al 2015
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo_dang(z_sur_std,snow_grain,dust_con,coszm, &
                              alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)

    implicit none

    real(wp), intent(in) :: z_sur_std
    real(wp), intent(in) :: snow_grain, dust_con, coszm
    real(wp), intent(out) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif

    real(wp) :: r, rn, c
    real(wp) :: x, f, H, p
    real(wp) :: dalb_sigma_orog
    real(wp) :: alpha_age_vis_dif, alpha_age_nir_dif, alpha_age_vis_dir, alpha_age_nir_dir
    real(wp) :: dalpha_vis_dif, dalpha_nir_dif, dalpha_vis_dir, dalpha_nir_dir
    real(wp), parameter :: r0 = 100._wp ! um
    real(wp), parameter :: c0 = 1.e-6_wp ! kg/kg
    real(wp), parameter :: dust_min = 1.e-8_wp  ! kg/kg

    ! snow grain radius
    r = snow_grain
    rn = log10(r/r0)
    if (dust_con.gt.dust_min) then
      x = log10(dust_con*1.e6_wp)
    endif

    dalb_sigma_orog = snow_par%k_sigma_orog*tanh(z_sur_std/snow_par%sigma_orog_crit)

    !----------------------------
    ! diffuse visible snow albedo

    ! diffuse vis albedo including snow aging effect, eq. 2
    alpha_age_vis_dif = 0.9856_wp + snow_par%dalb_snow_vis - 0.0202_wp*rn - 0.0125_wp*rn**2 
    ! effect of black carbon (equivalent)
    if (dust_con.gt.dust_min) then
      ! black carbon equivalent, diffuse vis radiation, eq. 9
      f = 152._wp + 15.92_wp*x - 0.39_wp*x**2
      c = dust_con/f 
      H = c/c0*(r/r0)**0.73_wp
      p = log10(H)
      dalpha_vis_dif = 10._wp**(-0.050_wp*p**2+0.514_wp*p-0.890_wp)
    else
      dalpha_vis_dif = 0._wp
    endif
    ! visible diffuse snow albedo
    alb_snow_vis_dif = min(1._wp, alpha_age_vis_dif - dalpha_vis_dif)

    !----------------------------
    ! diffuse near-infrared snow albedo

    ! diffuse nir albedo including snow aging effect
    alpha_age_nir_dif = 0.7493_wp + snow_par%dalb_snow_nir - 0.1820_wp*rn - 0.0388_wp*rn**2
    ! effect of dust in the infrared is small
    dalpha_nir_dif = 0._wp
    ! near-infrared diffuse snow albedo
    alb_snow_nir_dif = min(1._wp, alpha_age_nir_dif - dalpha_nir_dif)

    !----------------------------
    ! clear-sky visible snow albedo

    ! solar zenith angle effect
    r = snow_grain*(1._wp+0.781_wp*(coszm-0.65_wp)**2) ! effective grain size corrected for zenith angle
    rn = log10(r/r0)
    ! direct albedo including snow aging effect
    alpha_age_vis_dir = 0.9849_wp + snow_par%dalb_snow_vis - 0.0215_wp*rn - 0.0132_wp*rn**2
    ! effect of black carbon (equivalent)
    if (dust_con.gt.dust_min) then
      ! black carbon equivalent, direct radiation, equation (9)
      f = 155._wp + 17.15_wp*x + 0.27_wp*x**2
      c = dust_con/f 
      H = c/c0*(r/r0)**0.73_wp
      p = log10(H)
      dalpha_vis_dir = 10._wp**(-0.049_wp*p**2+0.525_wp*p-0.893_wp)
    else
      dalpha_vis_dir = 0._wp
    endif
    ! allband direct snow albedo
    alb_snow_vis_dir = min(1._wp, alpha_age_vis_dir - dalpha_vis_dir)

    !----------------------------
    ! clear-sky near-infrared snow albedo

    ! solar zenith angle effect
    r = snow_grain*(1._wp+0.791_wp*(coszm-0.65_wp)**2) ! effective grain size corrected for zenith angle
    rn = log10(r/r0)
    ! direct albedo including snow aging effect
    alpha_age_nir_dir = 0.6596_wp + snow_par%dalb_snow_nir - 0.1927_wp*rn - 0.0229_wp*rn**2
    ! effect of dust in the infrared is small
    dalpha_nir_dir = 0._wp
    ! allband direct snow albedo
    alb_snow_nir_dir = min(1._wp, alpha_age_nir_dir - dalpha_nir_dir)

    !-----------------------------
    ! reduce albedo over rough orography

    alb_snow_vis_dif = alb_snow_vis_dif - dalb_sigma_orog
    alb_snow_vis_dir = alb_snow_vis_dir - dalb_sigma_orog
    alb_snow_nir_dif = alb_snow_nir_dif - dalb_sigma_orog
    alb_snow_nir_dir = alb_snow_nir_dir - dalb_sigma_orog

    return

  end subroutine snow_albedo_dang

  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r e s i s t a n c e 
  !   Purpose    :  compute exchange coefficients and aerodynamic resistance given snow depth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine resistance(h_snow, t2m, t_skin, wind, &
                       r_a)

    implicit none

    real(wp), intent(in) :: t2m, t_skin
    real(wp), intent(in) :: h_snow
    real(wp), intent(in) :: wind
    real(wp), intent(out) :: r_a

    real(wp) :: rough_m, rough_h, fsnow
    real(wp) :: Ch, Ri
    real(wp) :: Ch_neutral


    fsnow = h_snow/(h_snow+10._wp*surf_par%z0m_ice)
    rough_m = fsnow * surf_par%z0m_snow + (1._wp-fsnow) * surf_par%z0m_ice
    rough_h = rough_m * surf_par%zm_to_zh

    ! neutral exchange coefficient
    Ch_neutral = (karman/log(z_sfl/rough_m)) * (karman/log(z_sfl/rough_h))

    ! Richardson stability number 
    Ri = g * z_sfl * (1._wp-t_skin/t2m) / wind**2

    if (l_neutral) then
      ! neutral stratification
      Ch = Ch_neutral
    else
      ! account for atmospheric stability through a Ri number dependence
      if (Ri.lt.0._wp) then ! "unstable" stratification
        Ch = Ch_neutral * (1._wp - 2._wp*Ri)
      else ! "stable" stratification
        Ch = Ch_neutral
      endif
    endif

    ! aerodynamic resistance
    if (Ch*wind.gt.0._wp) then
      r_a = 1._wp / (Ch * wind)
    else
      r_a = 1.e20_wp
    endif


    return

  end subroutine resistance


end module smb_surface_par

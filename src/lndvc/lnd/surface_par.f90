!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s u r f a c e _ p a r _ l n d
!
!  Purpose : land surface parameters
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
module lndvc_surface_par_lnd

  use precision, only : wp
  use constants, only : karman, g, T0, frac_vu
  use lnd_grid, only : i_ice, i_lake, i_bare
  use lnd_grid, only : npft
  use lnd_params, only : l_neutral, z_sfl
  use lnd_params, only : snow_par, surf_par, veg_par

  implicit none

  private
  public :: surface_frac_up, snow_albedo_lake, surface_albedo_lake, resist_aer_lake, resist_sur_lake

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u r f a c e _ f r a c _ u p
  !   Purpose    :  compute surface type fractions
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surface_frac_up(f_ice,f_ice_grd,f_shelf,f_lake,f_veg,pft_frac, &
                            frac_surf)

    implicit none

    real(wp), intent(in) :: f_ice, f_ice_grd, f_shelf, f_lake, f_veg
    real(wp), dimension(:), intent(in) :: pft_frac
    real(wp), dimension(:), intent(inout) :: frac_surf

    integer :: n
    real(wp), dimension(npft) :: pft_frac_tmp


    where (pft_frac.gt.veg_par%seed_fraction) 
      pft_frac_tmp = pft_frac
    elsewhere
      pft_frac_tmp = 0._wp
    endwhere

     ! ice, shelf and lake are in absolute grid fractions, PFT are in fraction of the vegetated part 
     frac_surf(i_ice) = f_ice
     frac_surf(i_lake) = f_lake
     do n=1,npft
      frac_surf(n) = pft_frac_tmp(n) * f_veg
     enddo
     frac_surf(i_bare) = (1._wp-sum(pft_frac_tmp)) * f_veg
     if (frac_surf(i_bare).lt.0._wp) frac_surf(i_bare) = 0._wp

     if(minval(frac_surf).lt.0._wp) then
       print *,'negative surface fraction!'
       print *,'frac_surf',frac_surf
       print *,'f_veg',f_veg
       if (frac_surf(i_bare).lt.0._wp) frac_surf(i_bare) = 0._wp
       if (minval(frac_surf).lt.-1.e-10_wp) stop
     endif
     if(abs(sum(frac_surf)+f_ice_grd-f_ice+f_shelf) .gt. (1._wp+1.e-5_wp))  then
      print *,'sum surface frac',sum(frac_surf)
      print *,'frac_surf',frac_surf
      print *,'f_veg',f_veg
      print *,'f_shelf',f_shelf
      print *,'f_ice,f_ice_grd',f_ice,f_ice_grd
      stop
     endif
     if(sum(frac_surf)+f_ice_grd-f_ice+f_shelf .lt. (1._wp-1.e-5_wp))  then
      print *,'sum surface frac',sum(frac_surf)+f_ice_grd-f_ice+f_shelf
      print *,'frac_surf',frac_surf
      print *,'f_veg',f_veg
      print *,'f_shelf',f_shelf
      print *,'f_ice,f_ice_grd',f_ice,f_ice_grd
      stop
     endif

     return

  end subroutine surface_frac_up


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o
  !   Purpose    :  albedo of snow for each surface type
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo_lake(t_skin, snow, w_snow, w_snow_max, dust_dep, coszm, &
                        alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif, &
                        snow_grain, dust_con)

    implicit none

    real(wp), intent(in) :: t_skin, snow, w_snow, w_snow_max
    real(wp), intent(in) :: dust_dep, coszm
    real(wp), intent(inout) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
    real(wp), intent(inout) :: snow_grain
    real(wp), intent(inout) :: dust_con


    ! snow grain size
    if (snow_par%l_snow_aging) then
      call snow_grain_size(t_skin, snow, snow_grain)
    else
      snow_grain = snow_par%snow_grain_fresh
    endif

    ! dust effect on snow albedo
    if (snow_par%l_snow_dust) then
      ! dust concentration in top snow layer
      call dust_in_snow(dust_dep, snow, w_snow, w_snow_max, dust_con)
    else
      dust_con = 0._wp
    endif

    ! compute snow albedo
    if (snow_par%i_snow_albedo.eq.1) then
      ! climber-2 snow albedo parameterisation, following Warren & Wiscombe 1980
      call snow_albedo_ww(snow_grain, dust_con, coszm, &
                          alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)
    else if (snow_par%i_snow_albedo.eq.2) then
      ! snow albedo parameterisation following Dang et al 2015
      call snow_albedo_dang(snow_grain, dust_con, coszm, &
                            alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)
    endif


    return

    end subroutine snow_albedo_lake


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o _ w w
  !   Purpose    :  compute snow albedo following Warren & Wiscombe 1980
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo_ww(snow_grain, dust_con, coszm, &
                            alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)

    implicit none

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
    alb_snow_vis_dif = snow_par%alb_snow_vis_dif_new - f_age*c_age_vis - c_dust_new_vis
    alb_snow_nir_dif = snow_par%alb_snow_nir_dif_new - f_age*c_age_nir - c_dust_new_nir
    alb_snow_vis_dir = alb_snow_vis_dif + 0.4_wp*f_cosz*(1._wp-alb_snow_vis_dif)
    alb_snow_nir_dir = alb_snow_nir_dif + 0.4_wp*f_cosz*(1._wp-alb_snow_nir_dif)


    return

  end subroutine snow_albedo_ww


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s n o w _ a l b e d o _ d a n g
  !   Purpose    :  compute snow albedo following Dang et al 2015
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine snow_albedo_dang(snow_grain,dust_con,coszm, &
                              alb_snow_vis_dir, alb_snow_nir_dir, alb_snow_vis_dif, alb_snow_nir_dif)

    implicit none

    real(wp), intent(in) :: snow_grain, dust_con, coszm
    real(wp), intent(out) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif

    real(wp) :: r, rn, c
    real(wp) :: x, f, H, p
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

    !----------------------------
    ! diffuse visible snow albedo

    ! diffuse vis albedo including snow aging effect, eq. 2
    alpha_age_vis_dif = 0.9856_wp - 0.0202_wp*rn - 0.0125_wp*rn**2 
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
    alb_snow_vis_dif = alpha_age_vis_dif - dalpha_vis_dif

    !----------------------------
    ! diffuse near-infrared snow albedo

    ! diffuse nir albedo including snow aging effect
    alpha_age_nir_dif = 0.7493_wp - 0.1820_wp*rn - 0.0388_wp*rn**2
    ! effect of dust in the infrared is small
    dalpha_nir_dif = 0._wp
    ! near-infrared diffuse snow albedo
    alb_snow_nir_dif = alpha_age_nir_dif - dalpha_nir_dif

    !----------------------------
    ! clear-sky visible snow albedo

    ! solar zenith angle effect
    r = snow_grain*(1._wp+0.781_wp*(coszm-0.65_wp)**2) ! effective grain size corrected for zenith angle
    rn = log10(r/r0)
    ! direct albedo including snow aging effect
    alpha_age_vis_dir = 0.9849_wp - 0.0215_wp*rn - 0.0132_wp*rn**2
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
    alb_snow_vis_dir = alpha_age_vis_dir - dalpha_vis_dir

    !----------------------------
    ! clear-sky near-infrared snow albedo

    ! solar zenith angle effect
    r = snow_grain*(1._wp+0.791_wp*(coszm-0.65_wp)**2) ! effective grain size corrected for zenith angle
    rn = log10(r/r0)
    ! direct albedo including snow aging effect
    alpha_age_nir_dir = 0.6596_wp - 0.1927_wp*rn - 0.0229_wp*rn**2
    ! effect of dust in the infrared is small
    dalpha_nir_dir = 0._wp
    ! allband direct snow albedo
    alb_snow_nir_dir = alpha_age_nir_dir - dalpha_nir_dir


    return

  end subroutine snow_albedo_dang


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
    f_tage1 = exp( snow_par%f_age_t * min(0._wp, (t_skin-T0)) )
    ! melting snow temperature dependence for snow grain size
    f_tage2 = exp( min(0._wp, t_skin-T0) )
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
      dust_con_melt_fac = 1._wp + (w_snow_max-w_snow)/10._wp  ! melt of 10 kg/m2 swe causes a doubling of the dust concentration
      dust_con_melt_fac = max(dust_con_melt_fac,1._wp)
      dust_con_melt_fac = min(dust_con_melt_fac,5._wp)  ! limit to a five-fold increase
    else
      dust_con_melt_fac = 1._wp
    endif
    dust_con = dust_con_melt_fac*dust_con

    ! avoid underflow and negative values
    if (dust_con.lt.1.e-15_wp) dust_con = 0._wp
    ! limit dust concentration in snow to dust_con_max
    dust_con = min(dust_con,dust_con_max)


    return

  end subroutine dust_in_snow


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u r f a c e _ a l b e d o
  !   Purpose    :  albedo for each surface type
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine surface_albedo_lake(h_snow,coszm,f_lake_ice, &
                           alb_snow_vis_dir,alb_snow_vis_dif,alb_snow_nir_dir,alb_snow_nir_dif, &
                           f_snow, &
                           alb_vis_dir,alb_vis_dif,alb_nir_dir,alb_nir_dif,albedo)

    implicit none

    real(wp), intent(in) :: h_snow
    real(wp), intent(in) :: coszm
    real(wp), intent(in) :: f_lake_ice
    real(wp), intent(in) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
    real(wp), intent(out) :: f_snow
    real(wp), intent(out) :: alb_vis_dir, alb_vis_dif, alb_nir_dir, alb_nir_dif, albedo

    real(wp) :: f_ice, alb_dir_water


    ! snow fraction after Niu and Yang 2007, Roesch 2001
    f_snow = tanh(h_snow/(snow_par%c_fsnow*surf_par%z0m_lake_ice))
    f_ice  = f_lake_ice
    alb_dir_water = 0.05_wp/(max(0.01,coszm)+0.15_wp)     ! CLM4.5, eq. 9.1, Pivoravov 1972
    alb_vis_dir = (1._wp-f_ice) * alb_dir_water &
      + f_ice * (f_snow * alb_snow_vis_dir + (1._wp-f_snow) * surf_par%alb_vis_dir_ice)
    alb_vis_dif = (1._wp-f_ice) * surf_par%alb_vis_dif_water &
      + f_ice * (f_snow * alb_snow_vis_dif + (1._wp-f_snow) * surf_par%alb_vis_dif_ice)
    alb_nir_dir = (1._wp-f_ice) * alb_dir_water &
      + f_ice * (f_snow * alb_snow_nir_dir + (1._wp-f_snow) * surf_par%alb_nir_dir_ice)
    alb_nir_dif = (1._wp-f_ice) * surf_par%alb_nir_dif_water &
      + f_ice * (f_snow * alb_snow_nir_dif + (1._wp-f_snow) * surf_par%alb_nir_dif_ice)

    ! composite albedo, assuming half cloud cover, diagnostic only
    albedo = frac_vu * 0.5_wp*(alb_vis_dir + alb_vis_dif) &
      + (1._wp-frac_vu) * 0.5_wp*(alb_nir_dir + alb_nir_dif)

    return

  end subroutine surface_albedo_lake


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r e s i s t _ a e r 
  !   Purpose    :  compute drag coefficients and aerodynamic resistance given snow depth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine resist_aer_lake(h_snow,tatm,t_skin,wind, &
                       z0m,rough_m,rough_h,Ch,r_a,Ri)

    implicit none

    real(wp), intent(in) :: h_snow
    real(wp), intent(in) :: tatm, t_skin
    real(wp), intent(in) :: wind
    real(wp), intent(in) :: z0m
    real(wp), intent(out) :: rough_m, rough_h
    real(wp), intent(out) :: Ch, r_a, Ri

    real(wp) :: fsnow
    real(wp) :: u_star, Re
    real(wp) :: log_m, log_h, Ch_neutral

    real(wp), parameter :: nu = 1.461e-5    ! kinematic molecular viscosity (m2/s)


    ! account for snow cover (lake is non-vegetated: flag_pft=0, z0m fixed)
    fsnow = h_snow/(h_snow+10._wp*z0m)
    rough_m = fsnow * surf_par%z0m_snow + (1._wp-fsnow) * z0m ! roughness including snow

    log_m = karman/log(z_sfl/rough_m)

    ! roughness for heat and water
    if (surf_par%i_z0h.eq.1) then
      rough_h = surf_par%zm_to_zh_const * rough_m
    else if (surf_par%i_z0h.eq.2) then
      ! formulation following Brutsaert 1982, Kanda 2007
      u_star = log_m*wind
      Re = u_star*rough_m/nu
      rough_h = rough_m*exp(-(1.29*Re**0.25_wp-2._wp))
    else if (surf_par%i_z0h.eq.3) then
      ! Zilitinkevich 1995
      u_star = log_m*wind
      Re = u_star*rough_m/nu
      rough_h = rough_m*exp(-karman*0.1_wp*sqrt(Re))
    else if (surf_par%i_z0h.eq.4) then
      ! Yang 2008, kB^-1~2 for bare soil and other non-vegetated surfaces
      rough_h = rough_m*exp(-2._wp)
    endif

    log_h = karman/log(z_sfl/rough_h)

    ! neutral heat exchange coefficient
    Ch_neutral = log_m * log_h

    ! Richardson number
    Ri = g * 100._wp * (1._wp - t_skin / tatm) / wind**2

    if( l_neutral ) then
      ! neutral stratification
      Ch = Ch_neutral
    else
      ! account for atmospheric stability through a Ri number dependence following BATS
      if( Ri .lt. 0._wp ) then ! "unstable" stratification
        Ch = Ch_neutral * (1._wp - surf_par%f_Ri_unstab * Ri)
      else ! "stable" stratification
        Ch = Ch_neutral / (1._wp + surf_par%f_Ri_stab * Ri)
      endif
    endif

    ! aerodynamic resistance
    r_a = 1._wp / (Ch * wind)

    return

  end subroutine resist_aer_lake


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r e s i s t _ s u r _ l a k e
  !   Purpose    :  surface resistance to evapotranspiration, lake tile
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine resist_sur_lake(beta_s,r_s)

    implicit none

    real(wp), intent(out) :: beta_s, r_s

    r_s = 0._wp
    beta_s = 1._wp

    return

  end subroutine resist_sur_lake

end module lndvc_surface_par_lnd

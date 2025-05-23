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
module surface_par_lnd

  use precision, only : wp
  !use constants, only : karman, g, pi, T0, z_sfl, frac_vu
  use constants, only : karman, g, pi, T0, frac_vu
  use lnd_grid, only : i_ice, i_lake, i_bare, i_trees, i_grass, i_shrub, is_veg, is_ice, is_lake, flag_pft, flag_veg
  use lnd_grid, only : npft, nsurf, nsoil, ngrass, ntrees, nshrub, nl, dz
  use lnd_params, only : l_neutral, i_racan, p_cdense, z_sfl
  use lnd_params, only : pft_par, snow_par, hydro_par, surf_par, veg_par 

  implicit none

  private
  public :: surface_frac_up, snow_albedo, surface_albedo, resist_aer, resist_sur

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
  subroutine snow_albedo(f_veg, f_ice, f_lake, w_snow, w_snow_max, &
                        t_skin_veg, t_skin_ice, t_skin_lake, snow, dust_dep, coszm, &
                        alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif, &
                        snow_grain, dust_con)

    implicit none

    real(wp), dimension(:), intent(in) :: w_snow, w_snow_max
    real(wp), intent(in) :: f_veg, f_ice, f_lake
    real(wp), intent(in) :: t_skin_veg, t_skin_ice, t_skin_lake
    real(wp), dimension(:), intent(in) :: snow
    real(wp), intent(in) :: dust_dep, coszm
    real(wp), dimension(:), intent(inout) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
    real(wp), dimension(:), intent(inout) :: snow_grain
    real(wp), dimension(:), intent(inout) :: dust_con


    ! for vegetated grid cell portion
    if (f_veg.gt.0._wp) then

      ! snow grain size
      if (snow_par%l_snow_aging) then
        call snow_grain_size(t_skin_veg, snow(i_bare), &
                             snow_grain(is_veg))
      else
        snow_grain(is_veg) = snow_par%snow_grain_fresh
      endif

      ! dust effect on snow albedo
      if (snow_par%l_snow_dust) then
        ! dust concentration in top snow layer
        call dust_in_snow(dust_dep, snow(i_bare), w_snow(is_veg), w_snow_max(is_veg), &
                          dust_con(is_veg))
      else
        dust_con(is_veg) = 0._wp
      endif

      ! compute snow albedo
      if (snow_par%i_snow_albedo.eq.1) then
        ! climber-2 snow albedo parameterisation, following Warren & Wiscombe 1980
        call snow_albedo_ww(snow_grain(is_veg), dust_con(is_veg), coszm, &
                            alb_snow_vis_dir(is_veg), alb_snow_nir_dir(is_veg), alb_snow_vis_dif(is_veg), alb_snow_nir_dif(is_veg))
      else if (snow_par%i_snow_albedo.eq.2) then
        ! snow albedo parameterisation following Dang et al 2015
        call snow_albedo_dang(snow_grain(is_veg), dust_con(is_veg), coszm, &
                              alb_snow_vis_dir(is_veg), alb_snow_nir_dir(is_veg), alb_snow_vis_dif(is_veg), alb_snow_nir_dif(is_veg))
      endif

    else

      alb_snow_vis_dif(is_veg) = 0._wp 
      alb_snow_nir_dif(is_veg) = 0._wp 
      alb_snow_vis_dir(is_veg) = 0._wp 
      alb_snow_nir_dir(is_veg) = 0._wp 

    endif

    ! for ice grid cell portion
    if (f_ice.gt.0._wp) then

      ! snow grain size
      if (snow_par%l_snow_aging) then
        call snow_grain_size(t_skin_ice, snow(i_ice), &
                             snow_grain(is_ice))
      else
        snow_grain(is_ice) = snow_par%snow_grain_fresh
      endif

      ! dust effect on snow albedo
      if (snow_par%l_snow_dust) then
        ! dust concentration in top snow layer
        call dust_in_snow(dust_dep, snow(i_ice), w_snow(is_ice), w_snow_max(is_ice), &
                          dust_con(is_ice))
      else
        dust_con(is_ice) = 0._wp
      endif

      ! compute snow albedo
      if (snow_par%i_snow_albedo.eq.1) then
        ! climber-2 snow albedo parameterisation, following Warren & Wiscombe 1980
        call snow_albedo_ww(snow_grain(is_ice), dust_con(is_ice), coszm, &
                            alb_snow_vis_dir(is_ice), alb_snow_nir_dir(is_ice), alb_snow_vis_dif(is_ice), alb_snow_nir_dif(is_ice))
      else if (snow_par%i_snow_albedo.eq.2) then
        ! snow albedo parameterisation following Dang et al 2015
        call snow_albedo_dang(snow_grain(is_ice), dust_con(is_ice), coszm, &
                              alb_snow_vis_dir(is_ice), alb_snow_nir_dir(is_ice), alb_snow_vis_dif(is_ice), alb_snow_nir_dif(is_ice))
      endif

    else

      alb_snow_vis_dif(is_ice) = 0._wp 
      alb_snow_nir_dif(is_ice) = 0._wp 
      alb_snow_vis_dir(is_ice) = 0._wp 
      alb_snow_nir_dir(is_ice) = 0._wp 

    endif


    ! lake
    if (f_lake.gt.0._wp) then

      ! snow grain size
      if (snow_par%l_snow_aging) then
        call snow_grain_size(t_skin_lake, snow(i_lake), &
                             snow_grain(is_lake))
      else
        snow_grain(is_lake) = snow_par%snow_grain_fresh
      endif

      ! dust effect on snow albedo
      if (snow_par%l_snow_dust) then
        ! dust concentration in top snow layer
        call dust_in_snow(dust_dep, snow(i_lake), w_snow(is_lake), w_snow_max(is_lake), &
                          dust_con(is_lake))
      else
        dust_con(is_lake) = 0._wp
      endif

      ! compute snow albedo
      if (snow_par%i_snow_albedo.eq.1) then
        ! climber-2 snow albedo parameterisation, following Warren & Wiscombe 1980
        call snow_albedo_ww(snow_grain(is_lake), dust_con(is_lake), coszm, &
                            alb_snow_vis_dir(is_lake), alb_snow_nir_dir(is_lake), alb_snow_vis_dif(is_lake), alb_snow_nir_dif(is_lake))
      else if (snow_par%i_snow_albedo.eq.2) then
        ! snow albedo parameterisation following Dang et al 2015
        call snow_albedo_dang(snow_grain(is_lake), dust_con(is_lake), coszm, &
                              alb_snow_vis_dir(is_lake), alb_snow_nir_dir(is_lake), alb_snow_vis_dif(is_lake), alb_snow_nir_dif(is_lake))
      endif

    else

      alb_snow_vis_dif(is_lake) = 0._wp 
      alb_snow_nir_dif(is_lake) = 0._wp 
      alb_snow_vis_dir(is_lake) = 0._wp 
      alb_snow_nir_dir(is_lake) = 0._wp 

    endif


    return

    end subroutine snow_albedo


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
    real(wp) :: c_age_vis, c_age_nir, c_soot_vis, c_soot_nir

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

    real(wp) :: dust_con_melt_fac, snowmelt_eff, snow_eff

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
  subroutine surface_albedo(frac_surf,z_veg_std,h_snow,coszm,lai,sai,z0m,f_snow_can,f_lake_ice, &
                           alb_snow_vis_dir,alb_snow_vis_dif,alb_snow_nir_dir,alb_snow_nir_dif, &
                           alb_bare_vis,alb_bare_nir, &
                           f_snow, &
                           alb_vis_dir,alb_vis_dif,alb_nir_dir,alb_nir_dif,albedo)

    implicit none

    real(wp), dimension(:), intent(in) :: h_snow
    real(wp), intent(in) :: z_veg_std
    real(wp), intent(in) :: coszm
    real(wp), dimension(:), intent(in) :: z0m, frac_surf, f_snow_can
    real(wp), intent(in) :: f_lake_ice 
    real(wp), dimension(:), intent(in) :: lai, sai
    real(wp), dimension(:), intent(in) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
    real(wp), intent(in) :: alb_bare_vis, alb_bare_nir
    real(wp), dimension(:), intent(out) :: f_snow
    real(wp), dimension(:), intent(out) :: alb_vis_dir, alb_vis_dif, alb_nir_dir, alb_nir_dif, albedo

    integer :: n, n_s
    real(wp) :: lsai_eff, svf, svf_dir, svf_dif, f_ice, f_snow_b, f_snow_t
    real(wp) :: alb_dir_water
    real(wp) :: alb_vis_bcan, alb_nir_bcan
    real(wp) :: alb_vis_dir_can, alb_nir_dir_can, alb_vis_dif_can, alb_nir_dif_can
    real(wp) :: f_snow_fac_orog
    real(wp), parameter :: sai_svf_scale = 2._wp  ! scale factor of stem area index for albedo 
    real(wp), parameter :: z0m_bcan = 0.02_wp   ! m, roughness length below canopy
    real(wp), parameter :: eps = 1.e-10_wp
    


    ! ice
    if( frac_surf(i_ice) .gt. 0._wp ) then
      ! snow fraction after Niu and Yang 2007, Roesch 2001
      f_snow(i_ice) = tanh(h_snow(is_ice)/(snow_par%c_fsnow*z0m(i_ice))) 
      alb_vis_dir(i_ice) = f_snow(i_ice) * alb_snow_vis_dir(is_ice) + (1._wp-f_snow(i_ice)) * surf_par%alb_vis_dir_ice
      alb_vis_dif(i_ice) = f_snow(i_ice) * alb_snow_vis_dif(is_ice) + (1._wp-f_snow(i_ice)) * surf_par%alb_vis_dif_ice
      alb_nir_dir(i_ice) = f_snow(i_ice) * alb_snow_nir_dir(is_ice) + (1._wp-f_snow(i_ice)) * surf_par%alb_nir_dir_ice
      alb_nir_dif(i_ice) = f_snow(i_ice) * alb_snow_nir_dif(is_ice) + (1._wp-f_snow(i_ice)) * surf_par%alb_nir_dif_ice
    else
      f_snow(i_ice) = 0._wp
    endif

    ! lake
    if( frac_surf(i_lake) .gt. 0._wp ) then
       ! snow fraction after Niu and Yang 2007, Roesch 2001
       f_snow(i_lake) = tanh(h_snow(is_lake)/(snow_par%c_fsnow*surf_par%z0m_lake_ice)) 
       f_ice  = f_lake_ice
       alb_dir_water = 0.05_wp/(max(0.01,coszm)+0.15_wp)     ! CLM4.5, eq. 9.1, Pivoravov 1972
       alb_vis_dir(i_lake) = (1._wp-f_ice) * alb_dir_water &
         + f_ice * (f_snow(i_lake) * alb_snow_vis_dir(is_lake) + (1._wp-f_snow(i_lake)) * surf_par%alb_vis_dir_ice)
       alb_vis_dif(i_lake) = (1._wp-f_ice) * surf_par%alb_vis_dif_water &
         + f_ice * (f_snow(i_lake) * alb_snow_vis_dif(is_lake) + (1._wp-f_snow(i_lake)) * surf_par%alb_vis_dif_ice)
       alb_nir_dir(i_lake) = (1._wp-f_ice) * alb_dir_water &
         + f_ice * (f_snow(i_lake) * alb_snow_nir_dir(is_lake) + (1._wp-f_snow(i_lake)) * surf_par%alb_nir_dir_ice)
       alb_nir_dif(i_lake) = (1._wp-f_ice) * surf_par%alb_nir_dif_water &
         + f_ice * (f_snow(i_lake) * alb_snow_nir_dif(is_lake) + (1._wp-f_snow(i_lake)) * surf_par%alb_nir_dif_ice)
    else
      f_snow(i_lake) = 0._wp
    endif

    ! bare soil
    if( frac_surf(i_bare) .gt. 0._wp ) then
      ! orography factor for snow cover fraction 
      if (snow_par%l_fsnow_orog) then
        ! reduce snow cover fraction over rough topography following Roesch 2001, eq. 7
        f_snow_fac_orog = h_snow(is_veg)/(h_snow(is_veg)+snow_par%c_fsnow_orog*max(snow_par%z_veg_std_min,z_veg_std)+eps)
      else
        f_snow_fac_orog = 1._wp
      endif
      ! snow fraction after Niu and Yang 2007, Roesch 2001
      f_snow(i_bare) = tanh(h_snow(is_veg)/(snow_par%c_fsnow*z0m(i_bare))) * f_snow_fac_orog
      alb_vis_dir(i_bare) = f_snow(i_bare) * alb_snow_vis_dir(is_veg) + (1._wp-f_snow(i_bare)) * alb_bare_vis
      alb_vis_dif(i_bare) = f_snow(i_bare) * alb_snow_vis_dif(is_veg) + (1._wp-f_snow(i_bare)) * alb_bare_vis
      alb_nir_dir(i_bare) = f_snow(i_bare) * alb_snow_nir_dir(is_veg) + (1._wp-f_snow(i_bare)) * alb_bare_nir
      alb_nir_dif(i_bare) = f_snow(i_bare) * alb_snow_nir_dif(is_veg) + (1._wp-f_snow(i_bare)) * alb_bare_nir
    endif

    ! vegetation, Otto 2011
    ! trees
    do n=1,ntrees
     n_s = i_trees(n) ! surface type index (also equal to PFT index)
     if( frac_surf(n_s) .gt. 0._wp ) then

      ! sky view factor, including zenith angle dependence
      ! for diffuse radiation, fitted relation from Verseghy 1993 based on 3 zenith angles
      ! use zenith angle=45° -> cos(45._wp*pi/180._wp) = 0.7071
      svf_dif = exp( - (lai(n_s)+sai_svf_scale*sai(n_s)) * veg_par%ext_coef/0.7071_wp )
      ! for direct radiation
      svf_dir = exp( - (lai(n_s)+sai_svf_scale*sai(n_s)) * veg_par%ext_coef/max(0.01_wp,coszm) )

      ! albedo of snow covered ground below the canopy, use always diffuse snow albedo!
      ! orography factor for snow cover fraction 
      if (snow_par%l_fsnow_orog) then
        ! reduce snow cover fraction over rough topography following Roesch 2001, eq. 7
        f_snow_fac_orog = h_snow(is_veg)/(h_snow(is_veg)+snow_par%c_fsnow_orog*max(snow_par%z_veg_std_min,z_veg_std)+eps)
      else
        f_snow_fac_orog = 1._wp
      endif
      ! snow fraction after Niu and Yang 2007
      f_snow_b = tanh(h_snow(is_veg)/(snow_par%c_fsnow*z0m_bcan)) * f_snow_fac_orog
      alb_vis_bcan = f_snow_b * alb_snow_vis_dif(is_veg) + (1._wp-f_snow_b) * pft_par%alb_bg_vis(n_s) 
      alb_nir_bcan = f_snow_b * alb_snow_nir_dif(is_veg) + (1._wp-f_snow_b) * pft_par%alb_bg_nir(n_s)

      ! albedo of snow covered canopy
      f_snow_t = f_snow_can(n_s)
      alb_vis_dir_can = (1._wp-f_snow_t) * pft_par%alb_can_vis_dir(n_s) + f_snow_t * pft_par%alb_can_vis_dir_snow(n_s) 
      alb_vis_dif_can = (1._wp-f_snow_t) * pft_par%alb_can_vis_dif(n_s) + f_snow_t * pft_par%alb_can_vis_dif_snow(n_s)
      alb_nir_dir_can = (1._wp-f_snow_t) * pft_par%alb_can_nir_dir(n_s) + f_snow_t * pft_par%alb_can_nir_dir_snow(n_s)
      alb_nir_dif_can = (1._wp-f_snow_t) * pft_par%alb_can_nir_dif(n_s) + f_snow_t * pft_par%alb_can_nir_dif_snow(n_s)

      f_snow(n_s) = svf_dir*f_snow_b + (1._wp-svf_dir)*f_snow_t

      ! tile albedo, weighted mean
      alb_vis_dir(n_s) = &
                   svf_dir        * alb_vis_bcan          &  ! ground below canopy
                 + (1._wp-svf_dir) * alb_vis_dir_can         ! canopy
      alb_vis_dif(n_s) = &
                   svf_dif        * alb_vis_bcan          &  ! ground below canopy
                 + (1._wp-svf_dif) * alb_vis_dif_can         ! canopy
      alb_nir_dir(n_s) = &
                   svf_dir        * alb_nir_bcan          &  ! ground below canopy
                 + (1._wp-svf_dir) * alb_nir_dir_can         ! canopy
      alb_nir_dif(n_s) = &
                   svf_dif        * alb_nir_bcan          &  ! ground below canopy
                 + (1._wp-svf_dif) * alb_nir_dif_can         ! canopy

     else
       f_snow(n_s) = 0._wp 
     endif
    enddo

    ! grass 
    do n=1,ngrass
     n_s = i_grass(n)
     if( frac_surf(n_s) .gt. 0._wp ) then
      ! sky view factor, no zenith angle dependence
      svf = exp( - (lai(n_s)+sai(n_s)) * veg_par%ext_coef)
      ! orography factor for snow cover fraction 
      if (snow_par%l_fsnow_orog) then
        ! reduce snow cover fraction over rough topography following Roesch 2001, eq. 7
        f_snow_fac_orog = h_snow(is_veg)/(h_snow(is_veg)+snow_par%c_fsnow_orog*max(snow_par%z_veg_std_min,z_veg_std)+eps)
      else
        f_snow_fac_orog = 1._wp
      endif
      ! snow fraction after Niu and Yang 2007
      f_snow(n_s) = tanh(h_snow(is_veg)/(snow_par%c_fsnow*z0m(n_s))) * f_snow_fac_orog
      alb_vis_dir(n_s) = &
                   f_snow(n_s)                       * alb_snow_vis_dir(is_veg)      &  ! snow
                 + (1._wp-f_snow(n_s)) * (1._wp-svf) * pft_par%alb_can_vis_dir(n_s)  &  ! snowfree canopy
                 + (1._wp-f_snow(n_s)) * svf         * pft_par%alb_bg_vis(n_s)
      alb_vis_dif(n_s) = &
                   f_snow(n_s)                       * alb_snow_vis_dif(is_veg)      &  ! snow
                 + (1._wp-f_snow(n_s)) * (1._wp-svf) * pft_par%alb_can_vis_dif(n_s)  &  ! snowfree canopy
                 + (1._wp-f_snow(n_s)) * svf         * pft_par%alb_bg_vis(n_s)
      alb_nir_dir(n_s) = &
                   f_snow(n_s)                       * alb_snow_nir_dir(is_veg)      &  ! snow
                 + (1._wp-f_snow(n_s)) * (1._wp-svf) * pft_par%alb_can_nir_dir(n_s)  &  ! snowfree canopy
                 + (1._wp-f_snow(n_s)) * svf         * pft_par%alb_bg_nir(n_s)
      alb_nir_dif(n_s) = &
                   f_snow(n_s)                       * alb_snow_nir_dif(is_veg)      &  ! snow
                 + (1._wp-f_snow(n_s)) * (1._wp-svf) * pft_par%alb_can_nir_dif(n_s)  &  ! snowfree canopy
                 + (1._wp-f_snow(n_s)) * svf         * pft_par%alb_bg_nir(n_s)

     endif
    enddo

    ! shrubs
    do n=1,nshrub
     n_s = i_shrub(n) ! surface type index (also equal to PFT index)
     if( frac_surf(n_s) .gt. 0._wp ) then

      ! sky view factor, including zenith angle dependence
      ! effective LAI+SAI index, accounting for snow thickness
      lsai_eff = (lai(n_s)+sai_svf_scale*sai(n_s)) * (1._wp-tanh(h_snow(is_veg)/(snow_par%c_fsnow*z0m(n_s))))
      ! for diffuse radiation, fitted relation from Verseghy 1993 based on 3 zenith angles
      ! use zenith angle=45° -> cos(45._wp*pi/180._wp) = 0.7071
      svf_dif = exp( -lsai_eff * veg_par%ext_coef/0.7071_wp )
      ! for direct radiation
      svf_dir = exp( -lsai_eff * veg_par%ext_coef/max(0.01_wp,coszm) )

      ! albedo of snow covered ground below the canopy, use always diffuse snow albedo!
      ! orography factor for snow cover fraction 
      if (snow_par%l_fsnow_orog) then
        ! reduce snow cover fraction over rough topography following Roesch 2001, eq. 7
        f_snow_fac_orog = h_snow(is_veg)/(h_snow(is_veg)+snow_par%c_fsnow_orog*max(snow_par%z_veg_std_min,z_veg_std)+eps)
      else
        f_snow_fac_orog = 1._wp
      endif
      ! snow fraction after Niu and Yang 2007
      f_snow_b = tanh(h_snow(is_veg)/(snow_par%c_fsnow*z0m_bcan)) * f_snow_fac_orog
      alb_vis_bcan = f_snow_b * alb_snow_vis_dif(is_veg) + (1._wp-f_snow_b) * pft_par%alb_bg_vis(n_s) 
      alb_nir_bcan = f_snow_b * alb_snow_nir_dif(is_veg) + (1._wp-f_snow_b) * pft_par%alb_bg_nir(n_s)

      ! albedo of snow covered canopy
      f_snow_t = f_snow_can(n_s)
      alb_vis_dir_can = (1._wp-f_snow_t) * pft_par%alb_can_vis_dir(n_s) + f_snow_t * pft_par%alb_can_vis_dir_snow(n_s) 
      alb_vis_dif_can = (1._wp-f_snow_t) * pft_par%alb_can_vis_dif(n_s) + f_snow_t * pft_par%alb_can_vis_dif_snow(n_s)
      alb_nir_dir_can = (1._wp-f_snow_t) * pft_par%alb_can_nir_dir(n_s) + f_snow_t * pft_par%alb_can_nir_dir_snow(n_s)
      alb_nir_dif_can = (1._wp-f_snow_t) * pft_par%alb_can_nir_dif(n_s) + f_snow_t * pft_par%alb_can_nir_dif_snow(n_s)

      f_snow(n_s) = svf_dir*f_snow_b + (1._wp-svf_dir)*f_snow_t

      ! tile albedo, weighted mean
      alb_vis_dir(n_s) = &
                   svf_dir        * alb_vis_bcan          &  ! ground below canopy
                 + (1._wp-svf_dir) * alb_vis_dir_can         ! canopy
      alb_vis_dif(n_s) = &
                   svf_dif        * alb_vis_bcan      &  ! ground below canopy
                 + (1._wp-svf_dif) * alb_vis_dif_can         ! canopy
      alb_nir_dir(n_s) = &
                   svf_dir        * alb_nir_bcan      &  ! ground below canopy
                 + (1._wp-svf_dir) * alb_nir_dir_can         ! canopy
      alb_nir_dif(n_s) = &
                   svf_dif        * alb_nir_bcan      &  ! ground below canopy
                 + (1._wp-svf_dif) * alb_nir_dif_can         ! canopy

     else
       f_snow(n_s) = 0._wp 
     endif
    enddo


    ! composite albedo, assuming half cloud cover, diagnostic only (except for offline simulations?)
    do n=1,nsurf
      if( frac_surf(n) .gt. 0._wp ) then
        albedo(n) = frac_vu * 0.5_wp*(alb_vis_dir(n) + alb_vis_dif(n)) &
          + (1._wp-frac_vu) * 0.5_wp*(alb_nir_dir(n) + alb_nir_dif(n))
      else
        albedo(n) = 0._wp
      endif
    enddo

    return

  end subroutine surface_albedo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r e s i s t _ a e r 
  !   Purpose    :  compute drag coefficients and aerodynamic resistance given snow depth
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine resist_aer(frac_surf,veg_h,lai,sai,h_snow,tatm,t_skin,wind, &
                       z0m,rough_m,rough_h,Ch,r_a,r_a_can,Ri)

    implicit none

    real(wp), dimension(:), intent(in) :: veg_h, lai, sai
    real(wp), dimension(:), intent(in) :: frac_surf, tatm, t_skin
    real(wp), dimension(:), intent(in) :: h_snow
    real(wp), dimension(:), intent(in) :: wind
    real(wp), dimension(:), intent(inout) :: z0m, rough_m, rough_h
    real(wp), dimension(:), intent(out) :: Ch, r_a, Ri
    real(wp), dimension(:), intent(out) :: r_a_can

    integer :: n
    real(wp) :: fsnow, hsnow
    real(wp) :: u_star, Re
    real(wp) :: log_m, log_h, Ch_neutral
    real(wp) :: flsai
    
    real(wp), parameter :: nu = 1.461e-5    ! kinematic molecular viscosity (m2/s)


    do n=1,nsurf

      if( frac_surf(n) .gt. 0._wp ) then

        if( n .eq. i_ice ) then
          hsnow = h_snow(is_ice)
        elseif( n .eq. i_lake ) then
          hsnow = h_snow(is_lake)
        else
          hsnow = h_snow(is_veg)
        endif

        ! roughness for momentum
        if( flag_pft(n) .eq. 1 ) then
          z0m(n) = pft_par%hveg_z0_scale(n) * veg_h(n)  ! roughness for snow free
        endif
        ! account for snow cover
        fsnow = hsnow/(hsnow+10._wp*z0m(n))
        rough_m(n) = fsnow * surf_par%z0m_snow + (1._wp-fsnow) * z0m(n) ! roughness including snow

        log_m = karman/log(z_sfl/rough_m(n))

        ! roughness_for heat and water
        if (surf_par%i_z0h.eq.1) then
          rough_h(n) = surf_par%zm_to_zh_const * rough_m(n)
        else if (surf_par%i_z0h.eq.2) then
          ! formulation following Brutsaert 1982, Kanda 2007
          ! surface friction velocity
          u_star = log_m*wind(n)
          ! roughness Reynolds number
          Re = u_star*rough_m(n)/nu
          rough_h(n) = rough_m(n)*exp(-(1.29*Re**0.25_wp-2._wp))
        else if (surf_par%i_z0h.eq.3) then
          ! Zilitinkevich 1995
          ! surface friction velocity
          u_star = log_m*wind(n)
          ! roughness Reynolds number
          Re = u_star*rough_m(n)/nu
          rough_h(n) = rough_m(n)*exp(-karman*0.1_wp*sqrt(Re))
        else if (surf_par%i_z0h.eq.4) then
          if( flag_pft(n) .eq. 1 ) then
            ! Zilitinkevich 1995
            ! surface friction velocity
            u_star = log_m*wind(n)
            ! roughness Reynolds number
            Re = u_star*rough_m(n)/nu
            rough_h(n) = rough_m(n)*exp(-karman*0.1_wp*sqrt(Re))
          else
            ! Yang 2008, kB^-1~2 for bare soil and other non-vegetated surfaces
            rough_h(n) = rough_m(n)*exp(-2._wp)
          endif
        endif

        log_h = karman/log(z_sfl/rough_h(n))

        ! neutral heat exchange coefficient
        Ch_neutral = log_m * log_h 

        ! Richardson number
        Ri(n) = g * 100._wp * (1._wp - t_skin(n) / tatm(n)) / wind(n)**2 

        if( l_neutral ) then
          ! neutral stratification
          Ch(n) = Ch_neutral
        else
          ! account for atmospheric stability through a Ri number dependence following BATS
          if( Ri(n) .lt. 0._wp ) then ! "unstable" stratification
            Ch(n) = Ch_neutral * (1._wp - surf_par%f_Ri_unstab * Ri(n))
          else ! "stable" stratification
            Ch(n) = Ch_neutral / (1._wp + surf_par%f_Ri_stab * Ri(n))
          endif
        endif

        ! aerodynamic resistance
        r_a(n) = 1._wp / (Ch(n) * wind(n))

        ! aerodynamic resistance for ground below canopy
        if( flag_pft(n) .eq. 1 ) then
          if (i_racan.eq.1) then
            r_a_can(n) = 1._wp/(p_cdense*wind(n)) 
          else if (i_racan.eq.2) then
            r_a_can(n) = (1._wp-exp(-(lai(n)+10._wp*sai(n))))/(p_cdense*wind(n)) 
          endif
        endif

      endif
    enddo

    if (i_racan.eq.3) then
      do n=1,npft
        r_a_can(n) = r_a(i_bare)
      enddo
    endif

    return

  end subroutine resist_aer


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r e s i s t _ s u r 
  !   Purpose    :  compute surface resistance to evapotranspiration
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine resist_sur(frac_surf,mask_snow,w_snow,theta,theta_w,theta_field,theta_sat,psi_sat,psi_exp,k_exp,g_can,beta_s,r_s,beta_s_can,r_s_can)

    implicit none

    integer,  dimension(:), intent(in) :: mask_snow
    real(wp), dimension(:), intent(in) :: w_snow
    real(wp), dimension(:), intent(in) :: frac_surf
    real(wp), dimension(:), intent(in) :: theta, theta_w, theta_field, theta_sat, psi_sat
    integer,  dimension(:), intent(in) :: psi_exp, k_exp
    real(wp), dimension(:), intent(in) :: g_can
    real(wp), dimension(:), intent(out) :: beta_s, r_s
    real(wp), dimension(:), intent(out) :: beta_s_can, r_s_can

    integer :: n
    real(wp) :: f_snow, beta_soil, beta_snow


    ! bare soil, use resistance OR beta factor
    r_s(i_bare) = 0._wp
    beta_snow = 1._wp
    if (hydro_par%i_evp_soil.eq.1) then
      ! CLM, Lee and Pielke 1992 
      if( theta_w(1) .lt. hydro_par%theta_crit_evp ) then
        beta_soil = 0.25_wp * (1._wp - cos(pi * theta_w(1) / hydro_par%theta_crit_evp))**2
      else
        beta_soil = 1._wp
      endif
    else if (hydro_par%i_evp_soil.eq.2) then
      ! CLIMBER-2
      if( theta_w(1) .lt. hydro_par%theta_crit_evp ) then
        beta_soil = (theta_w(1) / hydro_par%theta_crit_evp)**2
      else
        beta_soil = 1._wp
      endif
    endif
!    if (w_snow(is_veg).gt.snow_par%w_snow_crit) then
!      beta_s(i_bare) = beta_snow
!    else
!      f_snow = w_snow(is_veg)/snow_par%w_snow_crit
!      beta_s(i_bare) = f_snow*beta_snow + (1._wp-f_snow)*beta_soil
!    endif
    if (mask_snow(is_veg).eq.0) then
      beta_s(i_bare) = beta_soil
    else
      beta_s(i_bare) = beta_snow
    endif
    ! scale to account for fraction of top layer from which evaporation can occur
    beta_s(i_bare) = hydro_par%dz_evp/dz(1) * beta_s(i_bare)

    ! lake
    if( frac_surf(i_lake) .gt. 0._wp ) then
      r_s(i_lake) = 0._wp
      beta_s(i_lake) = 1._wp
    endif

    ! ice
    if( frac_surf(i_ice) .gt. 0._wp ) then
      r_s(i_ice) = 0._wp
      beta_s(i_ice) = 1._wp
    endif

    ! vegetation
    do n=1,npft
      if( frac_surf(n) .gt. 0._wp ) then
        r_s(n) = 1._wp / max(1.e-5_wp, 1.6_wp*g_can(n)) ! s/m   , factor 1.6 is to convert from CO2 to H2O conductance
        beta_s(n) = 1._wp
        ! for evaporation from soil below canopy use bare soil resistance
        r_s_can(n) = r_s(i_bare)
        beta_s_can(n) = beta_s(i_bare)
      endif
    enddo

    return

  end subroutine resist_sur

end module surface_par_lnd

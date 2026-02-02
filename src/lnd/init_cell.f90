!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i n i t _ c e l l
!
!  Purpose : initialize cells
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
module init_cell

  use precision, only : wp
  use constants, only : rho_w, rho_i, T0
  use lnd_grid, only : dz, nl, nlc, nl_l, nsurf, npft, i_surf, i_ice, i_lake, i_bare, flag_veg, dz_l, z_l, z_int_l
  use lnd_params, only : i_init_veg, veg_par, pft_par, snow_par, surf_par, peat_par, dz_lake_ice


  implicit none

  private
  public :: init_cell_veg, init_cell_ice, init_cell_ice_grd, init_cell_shelf, init_cell_lake

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n i t _ c e l l _ v e g
  !   Purpose    :  Initialize prognostic variables of vegetated grid cell 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine init_cell_veg(c13_c12_atm,c14_c_atm, &
                          f_veg,f_shelf_old,frac_surf,pft_frac,t_shelf,theta_i_shelf,theta_w_shelf,theta_sat, &
                          t_skin,t_skin_veg,t_soil,w_can,s_can, &
                          w_snow,w_snow_max,h_snow,mask_snow, &
                          theta,theta_w,theta_i,w_w,w_i, &
                          alt,gdd5,gdd,phen,phen_acc,lai_bal,lai,sai,root_frac,litter_in_frac,gamma_dist, gamma_fire, &
                          npp_ann, npp13_ann, npp14_ann, &
                          leaf_c,root_c, &
                          stem_c,veg_h,veg_c,veg_c13,veg_c14, &
                          litter_c,litter_c13,litter_c14,fast_c,fast_c13,fast_c14,slow_c,slow_c13,slow_c14, &
                          cato_c,cato_c13,cato_c14,litter_c_peat,litter_c13_peat,litter_c14_peat, &
                          acro_c,acro_c13,acro_c14,f_peat,f_peat_pot,w_table_min,w_table_peat,dCpeat_dt, &
                          z0m) 

    real(wp), intent(in) :: c13_c12_atm, c14_c_atm
    real(wp), intent(in) :: f_veg, f_shelf_old
    real(wp), dimension(:), intent(inout) :: frac_surf 
    real(wp), dimension(:), intent(inout) :: pft_frac
    real(wp), dimension(0:), intent(in) :: t_shelf
    real(wp), dimension(:), intent(in) :: theta_i_shelf, theta_w_shelf, theta_sat
    real(wp), dimension(0:), intent(inout) :: t_soil
    real(wp), dimension(:), intent(inout) :: theta, theta_w, theta_i, w_w, w_i
    real(wp), dimension(:), intent(inout) :: w_can, s_can, t_skin
    real(wp), intent(inout) :: w_snow, w_snow_max, h_snow
    integer, intent(inout) :: mask_snow
    real(wp), intent(inout) :: alt, gdd5, t_skin_veg
    real(wp), dimension(:), intent(inout) :: gdd, phen, phen_acc, lai_bal, lai, sai, gamma_dist, gamma_fire
    real(wp), dimension(:,:), intent(inout) :: root_frac
    real(wp), dimension(:), intent(inout) :: litter_in_frac
    real(wp), dimension(:), intent(inout) :: npp_ann, npp13_ann, npp14_ann
    real(wp), dimension(:), intent(inout) :: leaf_c, root_c
    real(wp), dimension(:), intent(inout) :: stem_c, veg_h, veg_c, veg_c13, veg_c14
    real(wp), dimension(:), intent(inout) :: litter_c,litter_c13,litter_c14,fast_c,fast_c13,fast_c14,slow_c,slow_c13,slow_c14
    real(wp), dimension(:), intent(inout) :: cato_c, cato_c13, cato_c14
    real(wp), intent(inout) :: litter_c_peat, litter_c13_peat, litter_c14_peat, acro_c, acro_c13, acro_c14
    real(wp), intent(inout) :: f_peat, f_peat_pot, w_table_min, w_table_peat, dCpeat_dt
    real(wp), dimension(:), intent(inout) :: z0m

    integer :: n
    real(wp) :: t_skin_mean, sum_frac


    ! initialise prognostic variables

    alt   = -1._wp

    ! vegetation properties
    gdd5       = 5000._wp
    gdd        = 0._wp
    phen       = 0._wp
    phen_acc   = 0._wp
    if (i_init_veg.eq.1) then
      ! desert
      pft_frac   = veg_par%seed_fraction
      lai_bal    = pft_par%lai_min
    else if (i_init_veg.eq.2) then
      ! forest
      pft_frac   = veg_par%seed_fraction
      lai_bal        = pft_par%lai_min
      pft_frac(1:2)  = 0.5_wp
      lai_bal(1:2)   = pft_par%lai_max(1:2)  
    endif
    lai        = lai_bal
    sai        = lai_bal * veg_par%sai_scale
    npp_ann    = 0._wp
    npp13_ann  = 0._wp
    npp14_ann  = 0._wp
    leaf_c     = lai_bal/pft_par%sla
    root_c     = leaf_c
    stem_c     = pft_par%awl*lai_bal**pft_par%bwl
    veg_h      = pft_par%awh * lai_bal
    veg_c      = leaf_c + root_c + stem_c
    veg_c13    = veg_c * c13_c12_atm
    veg_c14    = veg_c * c14_c_atm 
    gamma_dist = 0.001_wp
    gamma_fire = 0._wp
    root_frac  = pft_par%root_frac
    litter_in_frac  = pft_par%litter_in_frac

    w_can   = 0._wp
    s_can   = 0._wp

    ! compute grid cell mean skin temperature and use to initialize
    sum_frac = sum(frac_surf,mask=flag_veg.eq.0)
    if( sum_frac .gt. 0._wp ) then
      t_skin_mean = 0._wp
      do n=1,nsurf
        if( flag_veg(n) .eq. 0 ) t_skin_mean = t_skin_mean + t_skin(n)*frac_surf(n)/sum_frac
      enddo
    else if (f_shelf_old.eq.1._wp) then
      t_skin_mean = t_shelf(0)
    else
      ! initialize
      t_skin_mean = T0 
    endif
    do n=1,nsurf
      if( flag_veg(n) .eq. 1 ) t_skin(n) = t_skin_mean
    enddo
    t_skin_veg = t_skin_mean

    w_snow         = 0._wp
    w_snow_max     = w_snow
    h_snow         = 0._wp
    mask_snow      = 0

    if (f_shelf_old.gt.0._wp) then
      ! use old shelf soil column state to initialize soil column
      t_soil(0)    = T0
      t_soil(1:nl) = t_shelf(1:nl)
      theta_i      = theta_i_shelf
      theta_w      = theta_w_shelf
    else
      ! initalize
      t_soil(0)      = T0
      t_soil(1:nl)   = t_skin_mean ! initialize soil temperature to be equal to skin temperature
      ! assume soil is saturated after ice/ocean shelf/lake are retreating
      ! if skin temperature below 0°C all frozen, otherwise all liquid
      if( t_skin_mean .lt. T0 ) then
        theta_i    = theta_sat
        theta_w    = 0._wp
      else
        theta_i    = 0._wp
        theta_w    = theta_sat
      endif
    endif

!    print *,f_shelf_old,f_veg
!    print *,frac_surf
!    print *,theta_w_shelf,theta_w
    w_w        = theta_w * dz(1:nl) * rho_w
    w_i        = theta_i * dz(1:nl) * rho_i
    theta      = w_w/(rho_w*dz(1:nl)) + w_i/(rho_i*dz(1:nl))

    ! initalize soil carbon
    litter_c   = 0._wp
    litter_c13 = 0._wp
    litter_c14 = 0._wp
    fast_c     = 0._wp
    fast_c13   = 0._wp
    fast_c14   = 0._wp
    slow_c     = 0._wp
    slow_c13   = 0._wp
    slow_c14   = 0._wp

    ! initialize peatlands
    litter_c_peat  = 0._wp
    litter_c13_peat= 0._wp
    litter_c14_peat= 0._wp
    acro_c         = 0._wp
    acro_c13       = 0._wp
    acro_c14       = 0._wp
    cato_c     = 0._wp
    cato_c13   = 0._wp
    cato_c14   = 0._wp
    if( peat_par%peat_area .or. peat_par%peat_carb ) then
      f_peat        = peat_par%f_peat_min*f_veg
      f_peat_pot    = peat_par%f_peat_min*f_veg
    else
      f_peat        = 0._wp
      f_peat_pot    = 0._wp
    endif
    w_table_min    = 0._wp
    w_table_peat   = 0._wp
    dCpeat_dt      = 0._wp

    do n=1,npft
      z0m(n)  = 0.1_wp * veg_h(n)
    enddo
    z0m(i_bare) = surf_par%z0m_bare

    return

  end subroutine init_cell_veg 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n i t _ c e l l _ i c e
  !   Purpose    :  Initialize prognostic variables for ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine init_cell_ice(frac_surf,t_skin,t_ice,w_snow,w_snow_max,h_snow,mask_snow)

    real(wp), dimension(:) :: t_skin, frac_surf
    real(wp) :: w_snow, w_snow_max, h_snow
    integer :: mask_snow
    real(wp), dimension(0:) :: t_ice

    integer :: n
    real(wp) :: t_skin_mean, sum_frac

     ! initialise prognostic variables

     ! compute grid cell mean skin temperature and use to initialize
     sum_frac = sum(frac_surf,mask=i_surf.ne.i_ice)
     if( sum_frac .gt. 0._wp ) then
      t_skin_mean = 0._wp
      do n=1,nsurf
       if( n .ne. i_ice ) t_skin_mean = t_skin_mean + t_skin(n)*frac_surf(n)/sum_frac
      enddo
     else ! all shelf, no skin temperature available
       t_skin_mean = T0
     endif
     t_skin(i_ice) = min(T0,t_skin_mean)

     w_snow         = 0._wp
     w_snow_max     = w_snow
     h_snow         = 0._wp
     mask_snow      = 0

     t_ice(0)      = T0
     t_ice(1:nl)   = min(T0,t_skin_mean) ! initialize ice temperature to be equal to skin temperature


     return

  end subroutine init_cell_ice


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n i t _ c e l l _ i c e _ g r d
  !   Purpose    :  Initialize prognostic variables for grounded ice
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine init_cell_ice_grd(litter_c_ice,litter_c13_ice,litter_c14_ice, &
                          fast_c_ice,fast_c13_ice,fast_c14_ice,slow_c_ice,slow_c13_ice,slow_c14_ice)

    real(wp), dimension(:) :: litter_c_ice,litter_c13_ice,litter_c14_ice
    real(wp), dimension(:) :: fast_c_ice,fast_c13_ice,fast_c14_ice,slow_c_ice,slow_c13_ice,slow_c14_ice


     ! initalize soil carbon
     litter_c_ice   = 0._wp
     litter_c13_ice = 0._wp
     litter_c14_ice = 0._wp
     fast_c_ice     = 0._wp
     fast_c13_ice   = 0._wp
     fast_c14_ice   = 0._wp
     slow_c_ice     = 0._wp
     slow_c13_ice   = 0._wp
     slow_c14_ice   = 0._wp


     return

  end subroutine init_cell_ice_grd


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n i t _ c e l l _ l a k e
  !   Purpose    :  Initialize prognostic variables for lake
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine init_cell_lake(frac_surf,t_skin,w_snow,w_snow_max,h_snow,mask_snow, &     
                            f_veg_old,t_soil,theta_i,theta_sat, &
                            t_sublake,w_w_sublake,w_i_sublake,theta_w_sublake,theta_i_sublake, &
                            t_lake,f_i_lake,f_lake_ice, &
                            litter_c_lake,litter_c13_lake,litter_c14_lake, &
                            fast_c_lake,fast_c13_lake,fast_c14_lake,slow_c_lake,slow_c13_lake,slow_c14_lake)

    real(wp), dimension(:) :: t_skin, frac_surf
    real(wp) :: w_snow, w_snow_max, h_snow
    integer :: mask_snow
    real(wp), intent(in) :: f_veg_old
    real(wp), dimension(0:), intent(in) :: t_soil
    real(wp), dimension(:), intent(in) :: theta_i, theta_sat
    real(wp), dimension(0:), intent(inout) :: t_sublake
    real(wp), dimension(0:), intent(inout) :: t_lake
    real(wp), dimension(:), intent(inout) :: w_w_sublake, w_i_sublake, theta_w_sublake, theta_i_sublake
    real(wp), dimension(:), intent(inout) :: f_i_lake 
    real(wp), intent(inout) :: f_lake_ice 
    real(wp), dimension(:) :: litter_c_lake,litter_c13_lake,litter_c14_lake
    real(wp), dimension(:) :: fast_c_lake,fast_c13_lake,fast_c14_lake,slow_c_lake,slow_c13_lake,slow_c14_lake

    integer :: n, k
    real(wp) :: t_skin_mean, sum_frac


     ! initialise prognostic variables

     ! compute grid cell mean skin temperature and use to initialize
     sum_frac = sum(frac_surf,mask=i_surf.ne.i_lake)
     if( sum_frac .gt. 0._wp ) then
      t_skin_mean = 0._wp
      do n=1,nsurf
       if( n .ne. i_lake ) t_skin_mean = t_skin_mean + t_skin(n)*frac_surf(n)/sum_frac
      enddo
     else ! all lake, no skin temperature available
       t_skin_mean = T0
     endif
     t_skin(i_lake) = t_skin_mean

     w_snow         = 0._wp
     w_snow_max     = w_snow
     h_snow         = 0._wp
     mask_snow      = 0

     t_lake(0)      = T0
     t_lake(1)      = t_skin_mean ! initialize top lake temperature to be equal to skin temperature
     t_lake(2:nl_l)   = max(t_skin_mean,T0+4._wp)  

     ! frozen water fraction
     do k=1,nl_l
       if (t_lake(k).ge.T0) then
         f_i_lake(k) = 0._wp 
       else
         f_i_lake(k) = 1._wp
       endif
     enddo
     f_lake_ice = f_i_lake(1)

    if (f_veg_old.gt.0._wp) then
     ! use temperature and water content from old vegetation soil column
     t_sublake(1:nl) = t_soil(1:nl)
     theta_i_sublake = theta_i
     theta_w_sublake = theta_sat - theta_i_sublake  ! fill the rest with liquid water
     w_w_sublake  = theta_w_sublake * dz(1:nl) * rho_w
     w_i_sublake  = theta_i_sublake * dz(1:nl) * rho_i
    else
     ! initialize
     t_sublake(1:nl) = T0 + 1._wp
     theta_i_sublake = 0.5_wp*theta_sat
     theta_w_sublake = 0.5_wp*theta_sat
     w_w_sublake  = theta_w_sublake * dz(1:nl) * rho_w
     w_i_sublake  = theta_i_sublake * dz(1:nl) * rho_i
    endif

    ! initalize soil carbon
    litter_c_lake   = 0._wp
    litter_c13_lake = 0._wp
    litter_c14_lake = 0._wp
    fast_c_lake     = 0._wp
    fast_c13_lake   = 0._wp
    fast_c14_lake   = 0._wp
    slow_c_lake     = 0._wp
    slow_c13_lake   = 0._wp
    slow_c14_lake   = 0._wp


     return

  end subroutine init_cell_lake


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n i t _ c e l l _ s h e l f
  !   Purpose    :  Initialize prognostic variables for shelf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine init_cell_shelf(f_veg_old,t_soil,theta_i,theta_sat, &
                            t_shelf,w_w_shelf,w_i_shelf,theta_w_shelf,theta_i_shelf, &
                            litter_c_shelf,litter_c13_shelf,litter_c14_shelf, &
                            fast_c_shelf,fast_c13_shelf,fast_c14_shelf,slow_c_shelf,slow_c13_shelf,slow_c14_shelf)

    implicit none    

    real(wp), intent(in) :: f_veg_old
    real(wp), dimension(0:), intent(in) :: t_soil
    real(wp), dimension(:), intent(in) :: theta_i, theta_sat
    real(wp), dimension(0:), intent(inout) :: t_shelf
    real(wp), dimension(:), intent(inout) :: w_w_shelf, w_i_shelf, theta_w_shelf, theta_i_shelf
    real(wp), dimension(:) :: litter_c_shelf,litter_c13_shelf,litter_c14_shelf
    real(wp), dimension(:) :: fast_c_shelf,fast_c13_shelf,fast_c14_shelf,slow_c_shelf,slow_c13_shelf,slow_c14_shelf


    if (f_veg_old.gt.0._wp) then
     ! use temperature and water content from old vegetation soil column
     t_shelf(0)    = t_soil(1)
     t_shelf(1:nl) = t_soil(1:nl)
     theta_i_shelf = theta_i
     theta_w_shelf = theta_sat - theta_i_shelf  ! fill the rest with liquid water
     w_w_shelf  = theta_w_shelf * dz(1:nl) * rho_w
     w_i_shelf  = theta_i_shelf * dz(1:nl) * rho_i
    else
     ! initialize
     t_shelf(0:nl) = T0 + 1._wp
     theta_i_shelf = 0.5_wp*theta_sat
     theta_w_shelf = 0.5_wp*theta_sat
     w_w_shelf  = theta_w_shelf * dz(1:nl) * rho_w
     w_i_shelf  = theta_i_shelf * dz(1:nl) * rho_i
    endif

    ! initalize soil carbon
    litter_c_shelf   = 0._wp
    litter_c13_shelf = 0._wp
    litter_c14_shelf = 0._wp
    fast_c_shelf     = 0._wp
    fast_c13_shelf   = 0._wp
    fast_c14_shelf   = 0._wp
    slow_c_shelf     = 0._wp
    slow_c13_shelf   = 0._wp
    slow_c14_shelf   = 0._wp

    return

  end subroutine init_cell_shelf

 
end module init_cell

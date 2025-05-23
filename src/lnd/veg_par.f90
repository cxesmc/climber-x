!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : v e g _ p a r _ m o d
!
!  Purpose : vegetation parameters
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
module veg_par_mod

  use precision, only : wp
  use timer, only : doy, mon, sec_year, nday_mon, nday_year, nmon_year, time_soy_lnd
  use constants, only : T0
  use lnd_grid, only : nl, z_int, npft, nsurf, flag_veg
  use lnd_params, only : dt, dt_day
  use lnd_params, only : pft_par, veg_par 

  implicit none

  private
  public :: phenology, dynveg_par, root_frac_update, litter_in_frac_update

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p h e n o l o g y
  !   Purpose    :  compute leaf area index, from LPJ
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine phenology(frac_surf,f_veg,t2m,t2m_min_mon,gdd5_temp,gdd5,gdd,phen_acc,phen,gamma_leaf, &
                      lai,lai_bal,i,j)

    implicit none

    real(wp), intent(in) :: f_veg, t2m_min_mon
    real(wp), dimension(:), intent(in) :: t2m, frac_surf
    real(wp), dimension(:), intent(in) :: lai_bal
    real(wp), intent(inout) :: gdd5_temp, gdd5
    real(wp), dimension(:), intent(inout) :: gdd, phen_acc, phen, gamma_leaf
    real(wp), dimension(:), intent(inout) :: lai

    integer :: n
    real(wp) :: lai_old, t2m_veg, f_deciduous, f_deciduous_T, f_deciduous_gdd5
    real(wp), parameter :: dT_dec = 5._wp ! degC
    real(wp), parameter :: dgdd5_dec = 200._wp 
    real(wp), parameter :: phen_max  = 250._wp !210._wp
    real(wp), parameter :: phen_min  = 10._wp

    integer :: i,j

    if( f_veg .gt. 0._wp ) then

      t2m_veg = 0._wp
      do n=1,nsurf
        if( flag_veg(n) .eq. 1 ) t2m_veg = t2m_veg + t2m(n)*frac_surf(n)/f_veg
      enddo
      if( time_soy_lnd ) gdd5_temp = 0._wp
      if( (t2m_veg - T0) .gt. 5._wp ) gdd5_temp = gdd5_temp + (t2m_veg - T0) * dt_day
      if( doy .eq. nday_year ) gdd5 = gdd5_temp

      do n=1,npft

        ! gdd above t_base_phen
        if( (t2m(n)-T0) .gt. pft_par%t_base_phen(n) ) then
          gdd(n) = gdd(n) + (t2m(n)-T0-pft_par%t_base_phen(n)) * dt_day
        endif

        ! factor for smooth transition between evergreen and deciduous
        f_deciduous_T = min(1._wp,max(0._wp,1._wp-1._wp/(2._wp*dT_dec)*(t2m_min_mon-T0-pft_par%t_cmon_phen(n)+dT_dec)))
        f_deciduous_gdd5 = min(1._wp,max(0._wp,1._wp-1._wp/(2._wp*dgdd5_dec)*(gdd5-pft_par%gdd5_phen(n)+dgdd5_dec)))
        f_deciduous = max(f_deciduous_T,f_deciduous_gdd5)

        ! deciduous
        if (f_deciduous.gt.0._wp) then

          ! save old lai
          lai_old = lai(n)

          ! budburst
          if( (t2m(n)-T0).gt.pft_par%t_base_phen(n) .and. phen_acc(n).lt.phen_max ) then
            if (phen(n).gt.gdd(n)/pft_par%ramp(n)) then
            else
              phen(n) = gdd(n) / pft_par%ramp(n)
            endif
            if( phen(n) .gt. 1._wp ) phen(n) = 1._wp
          endif

          ! senescence
          if( (t2m(n)-T0).lt.pft_par%t_base_phen(n) ) then
            phen(n) = phen(n) - phen(n)*veg_par%gamma_phen*dt
            if( phen(n) .lt. 0._wp ) phen(n) = 0._wp
            phen_acc(n) = 0._wp
            gdd(n) = 0._wp
          endif

          phen_acc(n) = phen_acc(n) + phen(n)

          ! evergreen
        else

          phen(n) = 1._wp

        endif

        lai(n) = (f_deciduous*phen(n)+(1._wp-f_deciduous)) * lai_bal(n) 

        gamma_leaf(n) = f_deciduous*1._wp/sec_year + (1._wp-f_deciduous)*pft_par%gamma_leaf(n)
        
      enddo

    else ! ice covered grid cell

      ! necessary?
      lai       = 0._wp
      phen      = 0._wp
      phen_acc  = 0._wp
      gdd       = 0._wp
      gdd5      = 0._wp

    endif

    return

  end subroutine phenology


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l i t t e r _ i n _ f r a c _ u p d a t e
  !   Purpose    :  adjust vertical litter input distribution so that all 
  !   input in the active layer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine litter_in_frac_update(alt,litter_in_frac)

    implicit none

    real(wp), intent(in) :: alt
    real(wp), dimension(:), intent(inout) :: litter_in_frac

    integer :: n, k


    ! first assign reference litter input distribution
    litter_in_frac = pft_par%litter_in_frac
    ! adjust for permafrost areas
    if (alt.gt.-1._wp) then
      ! apply active layer thickness constraint
      do k=1,nl
        if (z_int(k).gt.alt) litter_in_frac(k) = 0._wp
      enddo
      if (sum(litter_in_frac).eq.0._wp) litter_in_frac(1) = 1._wp
      ! normalize to 1
      litter_in_frac = litter_in_frac / sum(litter_in_frac)
    endif


    return

  end subroutine litter_in_frac_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r o o t _ f r a c _ u p d a t e
  !   Purpose    :  adjust vertical root distribution so that all roots in
  !   active layer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine root_frac_update(alt,root_frac)

    implicit none

    real(wp), intent(in) :: alt
    real(wp), dimension(:,:), intent(inout) :: root_frac

    integer :: n, k


    ! first assign reference root distributions
    root_frac = pft_par%root_frac
    ! adjust for permafrost areas
    if (alt.gt.-1._wp) then
      ! apply active layer thickness constraint
      do n=1,npft
        do k=1,nl
          if (z_int(k).gt.alt) root_frac(k,n) = 0._wp
        enddo
        if (sum(root_frac(:,n)).eq.0._wp) root_frac(1,n) = 1._wp
        ! normalize to 1
        root_frac(:,n) = root_frac(:,n) / sum(root_frac(:,n))
      enddo
    endif


    return

  end subroutine root_frac_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d y n v e g _ p a r
  !   Purpose    :  parameters for dynamic vegetation model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dynveg_par(disturbance,t2m_min_mon,gdd5,veg_c_above,theta_fire_cum,gamma_dist_cum)

    implicit none

    real(wp), dimension(:), intent(in) :: disturbance
    real(wp), intent(in) :: t2m_min_mon
    real(wp), intent(in) :: gdd5
    real(wp), intent(in) :: veg_c_above
    real(wp), intent(inout) :: theta_fire_cum
    real(wp), dimension(:), intent(inout) :: gamma_dist_cum

    integer :: n
    logical :: lim_bio
    real(wp) :: fac_lim_bio
    real(wp) :: gamma_fire, theta_fire


    if (mon.eq.1) gamma_dist_cum = 0._wp

    ! fire disturbance rate for vegetation dynamics, Reick 2013
    theta_fire = theta_fire_cum / real(nday_mon,wp) ! monthly average
    do n=1,npft
      ! increased disturbance rate if outside bioclimatic limit
      lim_bio = (t2m_min_mon-T0).ge.pft_par%t_cmon_min(n) .and. (t2m_min_mon-T0).le.pft_par%t_cmon_max(n) .and. gdd5.ge.pft_par%gdd5_min(n)
      if (.not.lim_bio) then
        fac_lim_bio = 100._wp
      else
        fac_lim_bio = 1._wp
      endif
      ! fire disturbance
      gamma_fire = 1._wp/pft_par%tau_fire(n) &
        * max(0._wp, (veg_par%theta_fire_crit-theta_fire)/veg_par%theta_fire_crit) !&  ! soil moisture factor
        !* max(0._wp, min(1._wp, (veg_c_above-veg_par%cveg_fire_low)/(veg_par%cveg_fire_high-veg_par%cveg_fire_low)))  ! fuel factor
      ! annual mean disturbance rate
      gamma_dist_cum(n) = gamma_dist_cum(n) &
        + (pft_par%gamma_dist_min(n)*fac_lim_bio &        ! minimum disturbance rate
        + gamma_fire &                        ! fire disturbance rate
        + disturbance(n)) &                   ! additional prescribed disturbance rate
        / real(nmon_year,wp)
    enddo

    ! reset cumulated variable
    theta_fire_cum = 0._wp

   return

  end subroutine dynveg_par


end module veg_par_mod

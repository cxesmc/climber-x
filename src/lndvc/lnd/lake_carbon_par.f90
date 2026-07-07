!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l a k e _ c a r b o n _ p a r _ m o d
!
!  Purpose : parameters for carbon in soil below lakes
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
module lndvc_lake_carbon_par_mod

  use precision, only : wp
  use timer, only : nstep_mon_lnd, mon
  use constants, only : T0, k_boltz
  use lnd_grid, only : z_int_c, z, z_c, dz_c, nl, nlc
  use lnd_params, only : soilc_par, ch4_par

  implicit none

  private
  public :: lake_carbon_par

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l a k e _ c a r b o n _ p a r
  !   Purpose    :  update lake carbon decomposition rate and diffusivity
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_carbon_par(t_lake_cum,theta_w_sublake_cum,theta_i_sublake_cum, &
                             k_litter_lake,k_fast_lake,k_slow_lake,diff_lakec,adv_lakec,ch4_frac_lake)

    implicit none

    real(wp), dimension(:), intent(inout) :: t_lake_cum, theta_w_sublake_cum, theta_i_sublake_cum 
    real(wp), dimension(:), intent(inout) :: k_litter_lake, k_fast_lake, k_slow_lake, diff_lakec, adv_lakec
    real(wp), dimension(:), intent(inout) :: ch4_frac_lake

    integer :: k
    real(wp), dimension(nl) :: t_lake_mean, theta_w_sublake_mean, theta_i_sublake_mean, ftemp, fmoist, fdepth, ftemp_ch4
    real(wp), dimension(nl) :: klitter_lake, kfast_lake, kslow_lake
    real(wp), dimension(nlc) :: diff, difflakec


    t_lake_mean       = t_lake_cum / nstep_mon_lnd
    theta_w_sublake_mean = theta_w_sublake_cum / nstep_mon_lnd
    theta_i_sublake_mean = theta_i_sublake_cum / nstep_mon_lnd

    ! temperature dependence of decomposition rate
    do k=1,nl
      if (soilc_par%iresp_temp.eq.1) then ! modified Arrhenius, as in LPJ
        if( t_lake_mean(k) .gt. 240._wp ) then
          ftemp(k) = exp(308.56_wp * (1._wp/56.02_wp - 1._wp/(46.02_wp+t_lake_mean(k)-T0)) )
        else
          ftemp(k) = 0._wp
        endif
      elseif (soilc_par%iresp_temp.eq.2) then ! Arrhenius
        if( t_lake_mean(k) .gt. 260._wp ) then
          ftemp(k) = exp(soilc_par%Ea*(1._wp/(k_boltz*283.15_wp)-1._wp/(k_boltz*t_lake_mean(k))))
        else
          ftemp(k) = 0._wp
        endif
      elseif (soilc_par%iresp_temp.eq.3) then ! Q10
        ftemp(k) = soilc_par%q10_c**((t_lake_mean(k)-283.15_wp)/10._wp)
      endif
    enddo

    ! temperature dependence of methane emission fraction of soil respiration
    do k=1,nl
      ftemp_ch4(k) = exp(ch4_par%Ea_ch4*(1._wp/(k_boltz*303.15_wp)-1._wp/(k_boltz*t_lake_mean(k))))
    enddo
    ch4_frac_lake = ch4_par%ch4_frac_lake * ftemp_ch4

    ! soil moisture dependence of decomposition rate
    do k=1,nl
      ! Porporato 2003
      if( theta_w_sublake_mean(k) .le. 0.3_wp ) then  ! 0.3 = theta at field capacity
        fmoist(k) = max(0._wp,theta_w_sublake_mean(k)/0.3_wp) ! linear increase below field capacity
      else
        fmoist(k) = 0.3_wp/theta_w_sublake_mean(k) ! hyperbolic decrease above field capacity
      endif
    enddo

    ! soil depth dependence of decomposition rate, Koven 2013
    do k=1,nl
      fdepth(k) = exp(-z(k)/soilc_par%z_tau)
    enddo

    ! decomposition rate, 1/s
    klitter_lake = soilc_par%k10_litter * ftemp*fmoist*fdepth 
    kfast_lake   = soilc_par%k10_fast   * ftemp*fmoist*fdepth 
    kslow_lake   = soilc_par%k10_slow   * ftemp*fmoist*fdepth 

    diff = soilc_par%diff_lake  ! m2/s
    ! carbon diffusivity at the soil levels (interfaces), m2/s
    if( soilc_par%diff_lake .gt. 0._wp ) then
      do k=1,nl-1
        difflakec(k) = diff(k) * diff(k+1) * (z_c(k+1) - z_c(k))  &
          / ( diff(k) * (z_c(k+1) - z_int_c(k)) + diff(k+1) * (z_int_c(k) - z_c(k)) )
      enddo
      difflakec(nl:nlc) = 0._wp
    else
      difflakec = 0._wp
    endif

    ! cumulate
    if (mon.eq.1) then
      k_litter_lake = 0._wp 
      k_fast_lake   = 0._wp 
      k_slow_lake   = 0._wp 
      diff_lakec    = 0._wp 
    endif

    k_litter_lake(1:nl) = k_litter_lake(1:nl) + klitter_lake
    k_fast_lake(1:nl)   = k_fast_lake(1:nl) + kfast_lake
    k_slow_lake(1:nl)   = k_slow_lake(1:nl) + kslow_lake
    ! no decomposition in burial layer
    k_litter_lake(nlc) = 0._wp
    k_fast_lake(nlc)   = 0._wp
    k_slow_lake(nlc)   = 0._wp

    diff_lakec    = diff_lakec + difflakec ! same value as for soil carbon in vegetated part
    adv_lakec     = 0._wp

    ! reset cumulated variables
    t_lake_cum  = 0._wp
    theta_w_sublake_cum = 0._wp
    theta_i_sublake_cum = 0._wp


    return

  end subroutine lake_carbon_par


end module lndvc_lake_carbon_par_mod

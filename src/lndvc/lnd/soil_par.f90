!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s o i l _ p a r _ m o d
!
!  Purpose : soil parameters
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
module lndvc_soil_par_mod
  
  use precision, only : wp
  use constants, only : lambda_i, lambda_w
  use constants, only : rho_w, rho_i, T0, cap_w, cap_i
  use lnd_grid, only : nl, nlc, z_int, z, dz, rdz
  use lnd_params, only : snow_par, soil_par, hydro_par, organic, soilc_par

  implicit none

  private
  public :: soil_par_update, soil_par_thermal, soil_par_hydro

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o i l _ p a r _ u p d a t e
  !   Purpose    :  update soil structure and related paramters
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_par_update(f_peat,litter_c,fast_c,slow_c,litter_c_peat,acro_c,cato_c, &
                            mineral_theta_sat,mineral_k_sat,mineral_psi_sat,mineral_Bi,mineral_lambda_s,mineral_lambda_dry, &
                            frac_soc,theta_sat,psi_sat,k_sat,k_exp,psi_exp,theta_field,theta_wilt,lambda_s,lambda_dry)

    implicit none

    real(wp), intent(in) :: f_peat, litter_c_peat, acro_c
    real(wp), dimension(:), intent(in) :: litter_c, fast_c, slow_c, cato_c
    real(wp), intent(in) :: mineral_theta_sat,mineral_k_sat,mineral_psi_sat,mineral_Bi,mineral_lambda_s,mineral_lambda_dry
    integer, dimension(:), intent(inout) :: k_exp, psi_exp
    real(wp), dimension(:), intent(inout) :: frac_soc
    real(wp), dimension(:), intent(inout) :: theta_sat, psi_sat, k_sat, theta_field, theta_wilt
    real(wp), dimension(:), intent(inout) :: lambda_s, lambda_dry

    integer :: k
    real(wp) :: Bi
    real(wp) :: soc


    do k=1,nl

     ! effect of soil carbon on soil properties, following partly Lawrence 2008, linear weighting

     if( k.eq.1 ) then
       soc = (1._wp-f_peat) * (litter_c(k) + fast_c(k) + slow_c(k)) &  ! kgC/m3
           + f_peat * (litter_c_peat*rdz(1) + acro_c*rdz(1))
     else
       soc = (1._wp-f_peat) * (litter_c(k) + fast_c(k) + slow_c(k)) &
           + f_peat * cato_c(k)
     endif

     ! fraction of soil organic carbon
     frac_soc(k) = min(1._wp, soc / soilc_par%rho_soc_max)

     ! update soil porosity
     if (.not.soil_par%constant_porosity) then
      ! porosity
      theta_sat(k) = (1._wp - frac_soc(k)) * mineral_theta_sat + frac_soc(k) * organic%theta_sat
     endif

     ! update soil hydrological parameters
     if (.not.soil_par%constant_soil_par_hydro) then
      ! hydraulic properties
      ! linear weighting, e.g. Lawrence 2008
      k_sat(k)   = (1._wp - frac_soc(k)) * mineral_k_sat + frac_soc(k) * organic%k_sat
      psi_sat(k) = (1._wp - frac_soc(k)) * mineral_psi_sat + frac_soc(k) * organic%psi_sat

      Bi = (1._wp - frac_soc(k)) * mineral_Bi + frac_soc(k) * organic%Bi
      k_exp(k) = 2*nint(Bi)+3
      psi_exp(k) = -nint(Bi)

      ! field capacity
      theta_field(k) = theta_sat(k) * (0.1_wp / (86400._wp * k_sat(k)))**(1._wp/k_exp(k))
      ! wilting point
      theta_wilt(k) = theta_sat(k) * (hydro_par%p_psi_min / psi_sat(k))**(1._wp/psi_exp(k))
     endif

     ! update soil thermal parameters
     if (.not.soil_par%constant_soil_par_therm) then
      ! thermal properties
      lambda_s(k)   = (1._wp - frac_soc(k)) * mineral_lambda_s + frac_soc(k) * organic%lambda_s
      lambda_dry(k) = (1._wp - frac_soc(k)) * mineral_lambda_dry + frac_soc(k) * organic%lambda_dry
     endif

    enddo

    return

  end subroutine soil_par_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o i l _ p a r _ t h e r m a l
  !   Purpose    :  update soil thermal properties
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_par_thermal(t_soil,theta_w,theta_i,theta,theta_sat, &
                             lambda_s,lambda_dry,h_snow, &
                             cap_soil,lambda_soil,lambda_int_soil)

    implicit none

    real(wp), dimension(0:), intent(in) :: t_soil
    real(wp), dimension(:), intent(in) :: theta_w, theta_i, theta, theta_sat
    real(wp), dimension(:), intent(in) :: lambda_s,lambda_dry
    real(wp), intent(in) :: h_snow
    real(wp), dimension(0:), intent(out) :: cap_soil, lambda_soil, lambda_int_soil

    integer :: k
    real(wp) :: lambda_sat
    real(wp) :: s_r, kersten
    real(wp), dimension(0:nl) :: z_loc


     z_loc(0) = -0.5_wp * h_snow
     z_loc(1:nl) = z(1:nl)

     ! snow heat capacity, J/m3/K
     cap_soil(0) = cap_i * snow_par%rho_snow

     ! snow thermal conductivity, W/m/K
     lambda_soil(0) = snow_par%lambda_snow

     do k=1,nl

      ! thermal properties
      lambda_sat = lambda_s(k)**(1._wp - theta_sat(k))               &
                 * lambda_w**(theta_w(k)/theta(k)*theta_sat(k)) &
                 * lambda_i**(theta_i(k)/theta(k)*theta_sat(k))
 
      s_r = theta(k)/theta_sat(k)
      if( t_soil(k) .ge. T0 ) then
       ! kersten = max(0._wp, log(s_r) + 1._wp)
!       kersten = max(0._wp, 1.428571_wp*(s_r-0.3_wp)) ! linear approximation of log version of kersten number
!       kersten = min(1._wp,max(0._wp, 1._wp/(1._wp-0.38_wp)*(s_r-0.33_wp)))
!       kersten = min(1._wp,max(0._wp, 1._wp/(1._wp-0.35_wp)*(s_r-0.3_wp)))
       kersten = min(1._wp,max(0._wp, 1._wp/(1._wp-0.35_wp)*(s_r-0.35_wp)))
      else
       kersten = s_r
      endif

      if( s_r .gt. 1.d-7 ) then 
       lambda_soil(k) = kersten * lambda_sat + (1._wp - kersten) * lambda_dry(k)
      else
       lambda_soil(k) = lambda_dry(k)
      endif

      ! heat capacity of soil, J/m3/K
      cap_soil(k) = (1._wp - theta_sat(k)) * soil_par%cap_s      &
                  + theta_w(k) * rho_w * cap_w  &
                  + theta_i(k) * rho_i * cap_i

     enddo

     ! soil thermal conductivity at the soil levels (interfaces), W/m/K
     do k=0,nl-1
      lambda_int_soil(k) = lambda_soil(k) * lambda_soil(k+1) * (z_loc(k+1) - z_loc(k))  &
                    / ( lambda_soil(k) * (z_loc(k+1) - z_int(k)) + lambda_soil(k+1) * (z_int(k) - z_loc(k)) )
     enddo
     lambda_int_soil(nl) = 0._wp


    return

  end subroutine soil_par_thermal


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o i l _ p a r _ h y d r o
  !   Purpose    :  update soil hydraulic properties
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine soil_par_hydro(theta_w,theta_sat,w_w,w_i,psi_sat,k_sat,k_exp,psi_exp, &
                           theta, &
                           psi,kappa_int)

    implicit none

    real(wp), dimension(:), intent(in) :: theta_w, theta_sat, w_w, w_i, psi_sat, k_sat
    integer, dimension(:), intent(in) :: psi_exp, k_exp
    real(wp), dimension(:), intent(inout) :: theta
    real(wp), dimension(:), intent(inout) :: psi, kappa_int

    integer :: k


    ! total volumetric water content (liquid + frozen)
    theta = w_w/(rho_w*dz(1:nl)) + w_i/(rho_i*dz(1:nl))

    ! psi, m
    if( hydro_par%i_cond_theta .eq. 1 ) then
     ! use liquid AND frozen water for metric potential formulation
     psi = psi_sat*(max(hydro_par%theta_min, theta/theta_sat))**psi_exp
    elseif( hydro_par%i_cond_theta .eq. 2 ) then
     ! use only liquid water, Swenson 2012
     psi = psi_sat*(max(hydro_par%theta_min, theta_w/theta_sat))**psi_exp
    elseif( hydro_par%i_cond_theta .eq. 3 ) then
     ! use only liquid water, Swenson 2012
     psi = psi_sat * ( max(hydro_par%theta_min, 0.5_wp*(theta+theta_w) / theta_sat) )**psi_exp
    endif
    psi = max(-1.e5_wp, psi)

    ! hydraulic conductivity, mm/s
    if( hydro_par%i_cond_theta .eq. 1 ) then 
     ! use liquid AND frozen water for conductivity formulation
     do k=1,nl-1
      kappa_int(k) = k_sat(k) &
                   * ( max(hydro_par%theta_min, (theta(k) + theta(k+1)) / (theta_sat(k) + theta_sat(k+1))) ) &
                   **k_exp(k)
     enddo
     k = nl
     kappa_int(k) = k_sat(k) &
                  * ( max(hydro_par%theta_min, theta(k) / theta_sat(k)) )**k_exp(k)
    elseif( hydro_par%i_cond_theta .eq. 2 ) then
     ! use only liquid water, recommended
     do k=1,nl-1
      kappa_int(k) = k_sat(k) &
                   * ( max(hydro_par%theta_min, (theta_w(k) + theta_w(k+1)) / (theta_sat(k) + theta_sat(k+1))) ) &
                   **k_exp(k)
     enddo
     k = nl
     kappa_int(k) = k_sat(k) &
                  * ( max(hydro_par%theta_min, theta_w(k) / theta_sat(k)) )**k_exp(k)
    elseif( hydro_par%i_cond_theta .eq. 3 ) then
     ! use only liquid water, recommended
     do k=1,nl-1
      kappa_int(k) = k_sat(k) &
                   * ( max(hydro_par%theta_min, 0.5_wp*(theta(k) + theta(k+1) + theta_w(k) + theta_w(k+1)) / (theta_sat(k) + theta_sat(k+1))) ) &
                   **k_exp(k)
     enddo
     k = nl
     kappa_int(k) = k_sat(k) &
                  * ( max(hydro_par%theta_min, 0.5_wp*(theta(k)+theta_w(k)) / theta_sat(k)) )**k_exp(k)
    endif


    return

  end subroutine soil_par_hydro

end module lndvc_soil_par_mod

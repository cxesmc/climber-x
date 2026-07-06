!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l a k e _ p a r _ m o d
!
!  Purpose : lake parameters
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
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
module lndvc_lake_par_mod

  use precision, only : wp
  use constants, only : lambda_a, lambda_i, cap_i, lambda_w, cap_w
  use constants, only : rho_w, rho_i, T0, karman, g
  use constants, only : lambda_i, lambda_w, cap_i, cap_w
  use constants, only : rho_w, rho_i
  use lnd_grid, only : z_int, z, z_l, nl, nl_l
  use lnd_params, only : snow_par
  use lnd_params, only : soil_par, K_eddy_lake_bg, K_eddy_lake_max
  use lake_rho_mod, only : lake_rho

  private
  public :: lake_par_thermal, sublake_par_thermal

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l a k e _ p a r _ t h e r m a l
  !   Purpose    :  update lake thermal properties
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_par_thermal(h_snow, h_lake, t_lake, f_i_lake, wind, lat, &
                             cap_lake, lambda_lake, lambda_int_lake)

    implicit none

    real(wp), intent(in) :: h_snow, h_lake
    real(wp), dimension(:), intent(in) :: t_lake
    real(wp), dimension(:), intent(in) :: f_i_lake
    real(wp), intent(in) :: wind 
    real(wp), intent(in) :: lat
    real(wp), dimension(0:), intent(out) :: cap_lake, lambda_lake, lambda_int_lake

    integer :: k
    real(wp) :: K_eddy, K_eddy_wind, lambdaw
    real(wp) :: w, k_star, exp_kstar_z, drho_dz, N2, Ri
    real(wp), dimension(0:nl_l) :: z_loc
    real(wp), dimension(0:nl_l) :: z_int_loc


     z_loc(0) = -0.5_wp*h_snow
     z_loc(1:nl_l-1) = z_l(1:nl_l-1)
     z_loc(nl_l) = (h_lake+0.5_wp*z_loc(nl_l-1))/1.5_wp 

     z_int_loc(0) = 0._wp ! snow - lake interface
     do k=1,nl_l-1
       z_int_loc(k) = 0.5_wp * ( z_loc(k) + z_loc(k+1) )
     enddo

     ! snow heat capacity, J/m3/K
     cap_lake(0) = cap_i*snow_par%rho_snow
     ! snow thermal conductivity, W/m/K
     lambda_lake(0) = snow_par%lambda_snow 

     ! friction velocity
     w = 0.0012_wp*max(0.1_wp,wind)
     ! latitudinally dependent parameter of the Ekman profile, Henderson-Sellers 1985
     k_star = 6.6_wp*max(0.1_wp,wind)**(-2)*sqrt(abs(lat))   

     do k=1,nl_l
       ! lake heat capacity, J/m3/K
       ! Note that the density of water is used for both ice and water fractions, as the thickness of the layer is fixed
       cap_lake(k) = (1._wp-f_i_lake(k))*rho_w*cap_w + f_i_lake(k)*rho_w*cap_i 
       ! lake thermal conductivity, including molecular and eddy diffusivities, W/m/K
       ! eddy diffusivity as sum of background value + wind driven eddy diffusivity
       if (k.ne.nl_l) then
         drho_dz = (lake_rho(t_lake(k+1))-lake_rho(t_lake(k))) / (z_loc(k+1)-z_loc(k))
       else
         drho_dz = (lake_rho(t_lake(k))-lake_rho(t_lake(k-1))) / (z_loc(k)-z_loc(k-1))
       endif
       drho_dz = max(0._wp,drho_dz)
       N2 = g/lake_rho(t_lake(k))*drho_dz
       exp_kstar_z = max(1.e-10_wp,exp(-k_star*z_loc(k)))
       Ri = (-1._wp+sqrt(1._wp+40._wp*N2*(karman*z_loc(k))**2/(w*exp_kstar_z)**2))/20._wp
       K_eddy_wind = (1._wp-f_i_lake(1)) * karman*w*z_loc(k)/(1._wp+37._wp*Ri**2)*exp_kstar_z   ! following CLM4.5 eq. 9.31
       K_eddy = K_eddy_lake_bg + K_eddy_wind     ! m2/s
       K_eddy = min(K_eddy, K_eddy_lake_max)
       lambdaw = lambda_w + K_eddy*rho_w*cap_w      ! W/m/K 
       lambda_lake(k) = lambdaw*lambda_i/(lambdaw*f_i_lake(k)+lambda_i*(1._wp-f_i_lake(k)))   ! CLM4.5 eq. 9.38
     enddo

     ! lake thermal conductivity at the interfaces, W/m/K
     do k=0,nl_l-1
       lambda_int_lake(k) = lambda_lake(k) * lambda_lake(k+1) * (z_loc(k+1) - z_loc(k))  &
         / ( lambda_lake(k) * (z_loc(k+1)-z_int_loc(k)) + lambda_lake(k+1) * (z_int_loc(k)-z_loc(k)) )
     enddo
     lambda_int_lake(nl_l) = 0._wp


    return

  end subroutine lake_par_thermal


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u b l a k e _ p a r _ t h e r m a l
  !   Purpose    :  update soil thermal properties below lake
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sublake_par_thermal(theta_w_sublake,theta_i_sublake,theta_sat,lambda_s, &
                              cap_sublake,lambda_int_sublake)

    implicit none

    real(wp), dimension(nl), intent(in) :: theta_w_sublake, theta_i_sublake, theta_sat, lambda_s
    real(wp), dimension(nl), intent(out) :: cap_sublake
    real(wp), dimension(0:nl), intent(out) :: lambda_int_sublake

    integer :: k
    real(wp), dimension(nl) :: lambda_sublake
    real(wp), dimension(0:nl) :: z_loc


     z_loc(0) = 0._wp
     z_loc(1:nl) = z(1:nl)

     do k=1,nl

      ! thermal properties
      lambda_sublake(k) = lambda_s(k)**(1._wp - theta_sat(k)) * lambda_w**theta_w_sublake(k) * lambda_i**theta_i_sublake(k)
 
      ! heat capacity of lake, J/m3/K
      cap_sublake(k) = (1._wp - theta_sat(k)) * soil_par%cap_s      &
                  + theta_w_sublake(k) * rho_w * cap_w  &
                  + theta_i_sublake(k) * rho_i * cap_i

     enddo

     ! lake thermal conductivity at the soil levels (interfaces), W/m/K
     lambda_int_sublake(0) = lambda_sublake(1)
     do k=1,nl-1
      lambda_int_sublake(k) = lambda_sublake(k) * lambda_sublake(k+1) * (z_loc(k+1) - z_loc(k))  &
                    / ( lambda_sublake(k) * (z_loc(k+1) - z_int(k)) + lambda_sublake(k+1) * (z_int(k) - z_loc(k)) )
     enddo
     lambda_int_sublake(nl) = 0._wp


    return

  end subroutine sublake_par_thermal

end module lndvc_lake_par_mod

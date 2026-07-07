!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : d u s t _ e m i s _ m o d
!
!  Purpose : dust emissions
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
module lndvc_dust_emis_mod

  use precision, only : wp
  use lnd_params, only : dust_par
  use lnd_grid, only : i_bare, i_grass, i_shrub

  implicit none

  private
  public :: dust_emission

contains

  subroutine dust_emission(frac_surf,z_veg_std,z_veg,z_veg_min,z_veg_max,f_snow,lai,sai,h_snow,tskin,tatm,theta_w,theta_i,wind, &
                      dust_emis_d, dust_emis_g, dust_emis_s, dust_emis)
!  Purpose: Computation of potential dust source area and dust emission

  implicit none

  real(wp), dimension(:), intent(in) :: frac_surf   !! surface type fractions [/]
  real(wp), intent(in) :: z_veg_std  !! standard deviation of surface topography [m]
  real(wp), intent(in) :: z_veg  !! mean surface elevation in grid cell [m]
  real(wp), intent(in) :: z_veg_min  !! min surface elevation in grid cell [m]
  real(wp), intent(in) :: z_veg_max  !! max surface elevation in grid cell [m]
  real(wp), intent(in) :: f_snow     !! snow cover fraction [/]
  real(wp), dimension(:), intent(in) :: lai, sai    !! leaf and stem area index [m2/m2]
  real(wp), intent(in) :: h_snow              !! snow thickness [m]
  real(wp), dimension(:), intent(in) :: tskin   !! skin temperature [K]
  real(wp), dimension(:), intent(in) :: tatm   !! atmospheric temperature [K]
  real(wp), intent(in) :: theta_w       !! liquid soil water content [m3/m3]
  real(wp), intent(in) :: theta_i       !! frozen soil water content [m3/m3]
  real(wp), intent(in) :: wind          !! surface wind speed [m/s]

  real(wp), intent(out) :: dust_emis_d, dust_emis_g, dust_emis_s, dust_emis

  integer :: i
  real(wp) :: theta, wind_eff
  real(wp) :: fac_sm, fac_wind, fac_stab, fac_lai, fac_snow, fac_topo
  real(wp), parameter :: h_snow_crit = 0.1_wp  ! m


  if (dust_par%l_dust_topo) then
    ! topographic erodibility after Ginoux 2001, see also Zender 2003
    if (z_veg_max.gt.z_veg_min) then
      fac_topo = (max(0._wp,z_veg_max-z_veg)/(z_veg_max-z_veg_min))**dust_par%topo_exp
    else
      fac_topo = 0._wp
    endif
  else
    fac_topo = 1._wp
  endif

  ! soil moisture factor for dust emission
  if (dust_par%i_theta.eq.1) then
    theta = theta_w+theta_i   ! assume both liquid and frozen water supress dust emissions
  else if (dust_par%i_theta.eq.2) then
    theta = theta_w   ! assume only liquid water supresses dust emissions
  endif
  fac_sm = 1._wp + tanh((theta - dust_par%sm_t) * dust_par%sm_n)
  ! wind factor for dust emission
  ! effective wind speed, accound for topographic roughness-dependent wind gusts (roughly derived from ERA5)
  wind_eff = wind * (1._wp+dust_par%wind_gust_fac*max(0._wp,z_veg_std-100._wp)/1000._wp)
  fac_wind = wind_eff**2*max(0._wp,wind_eff-dust_par%u0*fac_sm)
  ! snow factor for dust emission
  if (dust_par%i_fsnow.eq.1) then
    fac_snow = 1._wp-min(1._wp,h_snow/h_snow_crit)
  else if (dust_par%i_fsnow.eq.2) then
    fac_snow = 1._wp-f_snow
  endif

  ! Dust emission from DESERT area
  if (frac_surf(i_bare).gt.0._wp) then ! snowfree desert
    ! stratification factor
    if (dust_par%l_dust_stab) then
      fac_stab = max(0._wp,tskin(i_bare)-tatm(i_bare))
    else
      fac_stab = 1._wp
    endif
    dust_emis_d = dust_par%b0 * dust_par%qd * frac_surf(i_bare) * fac_wind * fac_stab * fac_snow * fac_topo
  else
    dust_emis_d = 0._wp
  endif

  ! Dust emission from GRASS area
  dust_emis_g = 0._wp
  do i=1,size(i_grass)
    if (frac_surf(i_grass(i)).gt.0._wp) then
      if (dust_par%l_dust_stab) then
        fac_stab = max(0._wp,tskin(i_grass(i))-tatm(i_grass(i)))
      else
        fac_stab = 1._wp
      endif
      fac_lai = 0.5_wp * (1._wp-tanh( (lai(i_grass(i))+sai(i_grass(i)) - dust_par%lai_t) * dust_par%lai_n) )
      dust_emis_g = dust_emis_g + dust_par%b0 * dust_par%qg * frac_surf(i_grass(i))*fac_lai * fac_wind * fac_stab * fac_snow * fac_topo
    endif
  enddo

  ! Dust emission from SHRUB area
  dust_emis_s = 0._wp
  do i=1,size(i_shrub)
    if (frac_surf(i_shrub(i)).gt.0._wp) then
      if (dust_par%l_dust_stab) then
        fac_stab = max(0._wp,tskin(i_shrub(i))-tatm(i_shrub(i)))
      else
        fac_stab = 1._wp
      endif
      fac_lai = 0.5_wp * (1._wp-tanh( (lai(i_shrub(i))+sai(i_shrub(i)) - dust_par%lai_t) * dust_par%lai_n) )
      dust_emis_s = dust_emis_s + dust_par%b0 * dust_par%qg * frac_surf(i_shrub(i))*fac_lai * fac_wind * fac_stab * fac_snow * fac_topo
    endif
  enddo

  ! total dust emission per grid cell area (kg/m2/s)
  dust_emis = dust_emis_d + dust_emis_g + dust_emis_s


  return

  end subroutine dust_emission

end module lndvc_dust_emis_mod

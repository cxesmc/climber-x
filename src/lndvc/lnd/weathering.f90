!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : w e a t h e r i n g _ m o d
!
!  Purpose : weathering
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
module lndvc_weathering_mod

  use precision, only : wp
  use constants, only : T0, c13_c12_std
  use control, only : d13c_weath
  use lnd_params, only : weath_gemco2_par, weath_uhh_par
  use lnd_params, only : i_weath_sc, kweath_scale, l_match_alk_weath

  implicit none

  private
  public :: weathering_gemco2, weathering_uhh

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  w e a t h e r i n g _ g e m c o 2
  !   Purpose    :  
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine weathering_gemco2(c13_c12_atm,c14_c_atm,weath_scale,f_veg,f_lit,runoff_ann, &
                        weath_carb,weath13_carb,weath14_carb, &
                        weath_sil,weath13_sil,weath14_sil,weath_loess)

    implicit none

    real(wp), intent(in) :: c13_c12_atm
    real(wp), intent(in) :: c14_c_atm
    real(wp), intent(in) :: weath_scale
    real(wp), intent(in) :: f_veg
    real(wp), dimension(:), intent(in) :: f_lit
    real(wp), intent(in) :: runoff_ann
    real(wp), intent(out) :: weath_carb, weath_sil, weath_loess
    real(wp), intent(out) :: weath13_carb, weath13_sil
    real(wp), intent(out) :: weath14_carb, weath14_sil

    integer :: n
    real(wp) :: factor, factor_temp, factor_npp
    real(wp) :: weath_scale_loc


    if (l_match_alk_weath) then
      weath_scale_loc = weath_scale
    else
      weath_scale_loc = 1._wp
    endif

    factor_temp = 1._wp ! todo
    factor_npp  = 1._wp ! todo
    factor = runoff_ann*factor_temp*factor_npp  ! kg/m2/yr

    ! bicarbonate (HCO3-) weathering rate
    weath_carb = 0._wp
    weath_sil  = 0._wp
    weath_loess = 0._wp
    ! sum over rock lithologies
    if (f_veg.gt.0._wp) then
      do n=1,weath_gemco2_par%nlit
        weath_carb = weath_carb + f_lit(n)*weath_gemco2_par%frac_carb(n)*weath_gemco2_par%kweath(n)*weath_scale_loc*kweath_scale * factor   ! mol C/m2/yr
        weath_sil  = weath_sil  + f_lit(n)*(1._wp-weath_gemco2_par%frac_carb(n))*weath_gemco2_par%kweath(n)*weath_scale_loc*kweath_scale * factor    ! mol C/m2/yr
      enddo
    endif

    weath13_carb = 0.5_wp*weath_carb*(d13c_weath/1000._wp+1._wp)*c13_c12_std + 0.5_wp*weath_carb*c13_c12_atm    ! half C from atm and half from rock
    weath13_sil  = weath_sil*c13_c12_atm    ! all C from atm

    weath14_carb = 0.5_wp*weath_carb*c14_c_atm    ! half C from atm and half from rock
    weath14_sil  = weath_sil*c14_c_atm    ! all C from atm


   end subroutine weathering_gemco2


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  w e a t h e r i n g _ u h h
  !   Purpose    :  weathering model from UHH (Börker and Hartmann) 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine weathering_uhh(c13_c12_atm,c14_c_atm,weath_scale,f_veg,f_lit,runoff_ann,t2m_ann_mean, &
                        weath_carb,weath13_carb,weath14_carb, &
                        weath_sil,weath13_sil,weath14_sil,weath_loess)

    implicit none

    real(wp), intent(in) :: c13_c12_atm
    real(wp), intent(in) :: c14_c_atm
    real(wp), intent(in) :: weath_scale
    real(wp), intent(in) :: f_veg
    real(wp), dimension(:), intent(in) :: f_lit
    real(wp), intent(in) :: runoff_ann
    real(wp), intent(in) :: t2m_ann_mean
    real(wp), intent(out) :: weath_carb, weath_sil, weath_loess
    real(wp), intent(out) :: weath13_carb, weath13_sil
    real(wp), intent(out) :: weath14_carb, weath14_sil

    integer :: n
    real(wp), parameter :: R = 8.3145_wp    ! J/mol/K, gas constant
    real(wp) :: weath_scale_loc


    if (l_match_alk_weath) then
      weath_scale_loc = weath_scale
    else
      weath_scale_loc = 1._wp
    endif

    ! bicarbonate (HCO3-) weathering rate

    weath_carb = 0._wp
    weath_sil  = 0._wp
    weath_loess = 0._wp
    if (f_veg.gt.0._wp) then
      ! sum over rock lithologies
      do n=1,weath_uhh_par%nlit
        if (f_lit(n).gt.0._wp) then
          if (n.eq.weath_uhh_par%i_sc) then
            ! carbonate sedimentary rocks
            if (i_weath_sc.eq.1) then
              ! Amiotte-Suchet and Probst (1995)
              weath_carb = weath_carb + f_lit(n)*3.1692_wp/1000._wp*runoff_ann *weath_scale_loc*kweath_scale    ! mol C/m2/yr
            else if (i_weath_sc.eq.2) then
              ! Romero-Mujalli et al. (2019)
              weath_carb = weath_carb + f_lit(n)*10._wp**(exp(-1.73_wp+0.28_wp*(t2m_ann_mean-T0)-0.0157_wp*(t2m_ann_mean-T0)**2))/1000._wp &
                *runoff_ann *weath_scale_loc*kweath_scale    ! mol C/m2/yr
            endif
          else if (n.eq.weath_uhh_par%i_loess) then
            ! loess, assume carbonate weathering behaviour
            weath_loess = f_lit(n)*10._wp**(9._wp/(1._wp+exp(-0.63_wp*(log10(max(1.e-10_wp,runoff_ann))+0.76_wp)))-2._wp)*1.e-6_wp *weath_scale_loc*kweath_scale ! mol C/m2/yr
            weath_carb = weath_carb + weath_loess 
          else
            ! other lithologies
            weath_carb = weath_carb + f_lit(n)*runoff_ann*weath_uhh_par%frac_carb(n) &
              * weath_uhh_par%b(n)*exp(1000._wp*weath_uhh_par%ca/R*(1._wp/284.2_wp-1._wp/t2m_ann_mean)) *weath_scale_loc*kweath_scale    ! mol C/m2/yr
            weath_sil  = weath_sil  + f_lit(n)*runoff_ann*(1._wp-weath_uhh_par%frac_carb(n)) &
              * weath_uhh_par%b(n)*exp(1000._wp*weath_uhh_par%sa(n)/R*(1._wp/284.2_wp-1._wp/t2m_ann_mean)) *weath_scale_loc*kweath_scale    ! mol C/m2/yr
          endif
        endif
      enddo
    endif

    ! TODO, 
!dC13 of carbonate minerals weathered: 1.6 permil
!dC13 of CO2 used to weather carbonates = dC13_atm - 20permil
!dC13 of CO2 used to weather silicates = dC13_atm - 20permil
!(it is assumed that the CO2 consumed by weathering processes stems from respiration in the soils)
    weath13_carb = 0.5_wp*weath_carb*(d13c_weath/1000._wp+1._wp)*c13_c12_std + 0.5_wp*weath_carb*c13_c12_atm    ! half C from atm and half from rock
    weath13_sil  = weath_sil*c13_c12_atm    ! all C from atm

    weath14_carb = 0.5_wp*weath_carb*c14_c_atm    ! half C from atm and half from rock, but rock assumed to be completely depleted in C14
    weath14_sil  = weath_sil*c14_c_atm    ! all C from atm


   end subroutine weathering_uhh


end module lndvc_weathering_mod


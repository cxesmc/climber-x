!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l a k e _ c o n v e c t i o n _ m o d
!
!  Purpose : convective adjustment in lake
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
module lake_convection_mod

  use precision, only : wp
  use constants, only : T0, rho_w, cap_w, cap_i
  use lnd_params, only : z_mix_lake_min
  use lnd_grid, only : nl_l
  use lake_rho_mod, only : lake_rho

  implicit none

  private
  public :: lake_convection

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l a k e _ c o n v e c t i o n
  !   Purpose    :  lake convective adjustment, for freshwater only
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_convection(t_lake, f_i_lake, dz, t_freeze, &
                             h_conv, h_mix)

    implicit none

    real(wp), dimension(:), intent(inout) :: t_lake
    real(wp), dimension(:), intent(inout) :: f_i_lake
    real(wp), dimension(:), intent(in) :: dz
    real(wp), intent(in) :: t_freeze
    real(wp), intent(out) :: h_conv 
    real(wp), intent(out) :: h_mix

    real(wp) :: tmx, tsm
    logical :: chk_la, chk_lb
    integer :: k, kk, k_mix
    integer :: kb, kt, la, lb
    real(wp) :: rl, ru, zsm
    real(wp) :: dz_mix, f_i_lake_avg, Q, t_frz, t_unfrz


    h_conv = 0._wp
    h_mix = 0._wp

    !-----------------------------------------------------------
    ! apply conventional convection scheme of Rahmstorf 1993 as implemented in MOM5
    ! for unfrozen layers

    ! search for unstable regions starting from the top unfrozen layer
    kt = 1
    do while (f_i_lake(kt)>0._wp .and. kt<nl_l)
      kt = kt+1
    enddo
    kb = kt+1
    do while (kt < nl_l)

      ! density 
      ru = lake_rho(t_lake(kt))
      rl = lake_rho(t_lake(kb))

      ! sum the first pair found in an unstable region
      if (ru > rl) then
        chk_la = .true.
        chk_lb = .true.
        zsm = dz(kt) + dz(kb)
        tsm = t_lake(kt)*dz(kt) + t_lake(kb)*dz(kb)
        tmx = tsm/zsm

        do while (chk_lb .or. chk_la)

          ! check for an unstable level (lb) below kb
          if (kb >= nl_l) chk_lb = .false.
          do while (chk_lb)
            chk_lb = .false.
            lb = kb + 1 
            ru = lake_rho(tmx)
            rl = lake_rho(t_lake(lb))

            if (ru > rl) then
              ! add new level to sums
              kb  = lb
              zsm = zsm + dz(kb)
              tsm = tsm + t_lake(kb)*dz(kb) 
              tmx = tsm/zsm

              chk_la = .true.
              if (kb < nl_l) chk_lb = .true.
            endif
          enddo

          ! check for an unstable level (la) above kt
          if (kt <= 1) chk_la = .false.
          do while (chk_la)
            chk_la = .false.
            la = kt - 1
            ru = lake_rho(t_lake(la))
            rl = lake_rho(tmx)

            if (ru > rl) then
              ! add new level to sums
              kt  = la
              zsm = zsm + dz(kt)
              tsm = tsm + t_lake(kt)*dz(kt)
              tmx = tsm/zsm

              chk_lb = .true.
              if (kt > 1) chk_la = .true.
            endif
          enddo

        enddo

        ! mix tracers from kt to kb
        do kk=kt,kb
          t_lake(kk) = tmx
        enddo

        ! convection diagnostics
        h_conv = max(h_conv,sum(dz(1:kb)))  ! maximum depth of convection

        kt = kb + 1

      else
        kt = kb
      endif

      ! continue the search for other unstable regions
      kb = kt + 1

    enddo   


    !-------------------------------------------------------
    ! additional mixing when ice in a layer that is below a layer which is not completely frozen
    ! following CLM4.5 section 9.5.9
    ! and mixing to a minimum depth

    ! find depth of mixing required to bring ice to the surface
    k_mix = 0
    do k=1,nl_l-1
      if ((f_i_lake(k)<1._wp .and. f_i_lake(k+1)>0._wp) .or. (sum(dz(1:k))<z_mix_lake_min)) then
      !if ((f_i_lake(k)<1._wp .and. f_i_lake(k+1)>0._wp) .or. (sum(dz(1:k))<((1._wp-f_i_lake(1))*z_mix_lake_min))) then
        k_mix = k+1
      endif
    enddo

    if (k_mix>1) then

      ! total thickness of layers to mix
      dz_mix = sum(dz(1:k_mix))
      h_mix = dz_mix

      ! average ice fraction over layers to be mixed
      f_i_lake_avg = sum(f_i_lake(1:k_mix)*dz(1:k_mix))/dz_mix

      ! compute total enthalpy
      Q = 0._wp
      do k=1,k_mix
        Q = Q + dz(k)*rho_w*(t_lake(k)-t_freeze)*((1._wp-f_i_lake(k))*cap_w + f_i_lake(k)*cap_i)
      enddo

      if (Q>0._wp) then
        t_frz   = t_freeze
        t_unfrz = Q/(rho_w*dz_mix*(1._wp-f_i_lake_avg)*cap_w) + t_freeze
      else
        t_frz   = Q/(rho_w*dz_mix*f_i_lake_avg*cap_i) + t_freeze
        t_unfrz = t_freeze
      endif

      do k=1,k_mix
        if (sum(dz(1:k)) < dz_mix*f_i_lake_avg) then
          ! layer is completely frozen and temperature at freezing point
          f_i_lake(k) = 1._wp
          t_lake(k) = t_frz
        else if (sum(dz(1:(k-1))) < dz_mix*f_i_lake_avg) then
          ! layer contains both ice and water 
          f_i_lake(k) = (dz_mix*f_i_lake_avg-sum(dz(1:(k-1))))/dz(k)
          t_lake(k) = (t_frz*f_i_lake(k)*cap_i + t_unfrz*(1._wp-f_i_lake(k))*cap_w) / (f_i_lake(k)*cap_i + (1._wp-f_i_lake(k))*cap_w) 
        else
          ! no ice in layer
          f_i_lake(k) = 0._wp
          t_lake(k) = t_unfrz
        endif
      enddo

    endif


   return
  
  end subroutine lake_convection

end module lake_convection_mod

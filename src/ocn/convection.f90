!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o n v e c t i o n _ m o d
!
!  Purpose : convective adjustment
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was partly ported from the original c-GOLDSTEIN model,
! see Edwards and Marsh (2005)
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
module convection_mod

  ! nconv counts occurences of mixing at each point not including the point at the top of each mixed region
  ! dconv is the maximum depth of convection
  ! kven is the integer depth of ventilation from the surface 
  ! dven is the depth of ventilation from the surface

  use precision, only : wp
  use ocn_grid
  use ocn_params, only : i_conv_shuffle, l_conv_shuffle_passive, l_mix_bgc_all, n_tracers_tot, n_tracers_ocn
  use eos_mod

  implicit none

  private
  public :: convection

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c o n v e c t i o n
  !   Purpose    :  convection code suitable for arbitrary functions rho(T,S)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine convection(l_trans_tracers,ts,rho,kbo,mask_coast,nconv,dconv,kven,dven)

    implicit none

    logical, dimension(:), intent(in) :: l_trans_tracers
    real(wp), dimension(:,:), intent(inout) :: ts
    real(wp), dimension(:), intent(inout) :: rho
    integer, intent(in) :: kbo
    integer, intent(in) :: mask_coast
    integer, intent(inout) :: kven, nconv
    real(wp), intent(inout) :: dven, dconv

    real(wp), dimension(3) :: tmx, tsm
    logical :: chk_la, chk_lb
    integer :: kk, l
    integer :: kb, kt, la, lb
    real(wp) :: rl, ru, zsm


    ! reset count of number of mixed layers 
    nconv = 0
    kven  = maxk
    dconv = 0._wp
    dven  = 0._wp

    if (i_conv_shuffle.eq.1 .or. (i_conv_shuffle.eq.2 .and. mask_coast.eq.1)) then
      ! Mueller convection scheme from Bern3D after Müller et al 
      ! 'ventilation time scales in an efficient 3-d ocean model'
      ! first apply coshuffle to mix directly down to density level
      ! don't comment out the call without ensuring kven is initialised
      call coshuffle(l_trans_tracers,kven,ts,rho,kbo)
    endif

    ! convection scheme of Rahmstorf 1993 as implemented in MOM5

    ! search for unstable regions starting from the top
    kt = maxk
    kb = maxk-1
    do while (kt > kbo)

      ! density referred to level interface depth between layers
      ru = eos(ts(kt,1),ts(kt,2),zw(kb))
      rl = eos(ts(kb,1),ts(kb,2),zw(kb))

      ! sum the first pair found in an unstable region
      if (ru > rl) then
        chk_la = .true.
        chk_lb = .true.
        zsm    = dz(kt) + dz(kb)
        tsm(1) = ts(kt,1)*dz(kt) + ts(kb,1)*dz(kb)
        tmx(1) = tsm(1)/zsm
        tsm(2) = ts(kt,2)*dz(kt) + ts(kb,2)*dz(kb)
        tmx(2) = tsm(2)/zsm

        do while (chk_lb .or. chk_la)

          ! check for an unstable level (lb) below kb
          if (kb <= kbo) chk_lb = .false.
          do while (chk_lb)
            chk_lb = .false.
            lb = kb - 1 
            ru = eos(tmx(1),tmx(2),zw(lb))
            rl = eos(ts(lb,1),ts(lb,2),zw(lb))

            if (ru > rl) then
              ! add new level to sums
              kb     = lb
              zsm    = zsm + dz(kb)
              tsm(1) = tsm(1) + ts(kb,1)*dz(kb) 
              tmx(1) = tsm(1)/zsm
              tsm(2) = tsm(2) + ts(kb,2)*dz(kb)
              tmx(2) = tsm(2)/zsm

              chk_la = .true.
              if (kb > kbo) chk_lb = .true.
            endif
          enddo

          ! check for an unstable level (la) above kt
          if (kt >= maxk) chk_la = .false.
          do while (chk_la)
            chk_la = .false.
            la = kt + 1
            ru = eos(ts(la,1),ts(la,2),zw(kt))
            rl = eos(tmx(1),tmx(2),zw(kt))

            if (ru > rl) then
              ! add new level to sums
              kt     = la
              zsm    = zsm + dz(kt)
              tsm(1) = tsm(1) + ts(kt,1)*dz(kt)
              tmx(1) = tsm(1)/zsm
              tsm(2) = tsm(2) + ts(kt,2)*dz(kt)
              tmx(2) = tsm(2)/zsm

              chk_lb = .true.
              if (kt < maxk) chk_la = .true.
            endif
          enddo

        enddo

        ! mix tracers from kt to kb
        do kk=kb,kt
          ts(kk,1) = tmx(1)
          ts(kk,2) = tmx(2)
          rho(kk) = eos(tmx(1),tmx(2),zro(kk))   ! update density of mixed layers
        enddo
        ! mix tracers other than temperature and salinity
        do l=3,n_tracers_tot
          if (l_mix_bgc_all) then
            tsm(3) = 0.0
            do kk=kt,kb
              tsm(3) = tsm(3) + ts(kk,l)*dz(kk) 
            enddo
            tmx(3) = tsm(3)/zsm 
            do kk=kt,kb
              ts(kk,l) = tmx(3)
            enddo
          else if (l_trans_tracers(l)) then
            tsm(3) = 0.0
            do kk=kt,kb
              tsm(3) = tsm(3) + ts(kk,l)*dz(kk) 
            enddo
            tmx(3) = tsm(3)/zsm 
            do kk=kt,kb
              ts(kk,l) = tmx(3)
            enddo
          endif
        enddo

        ! convection diagnostics
        nconv = nconv + kt - kb    ! number of mixed layers
        dconv = min(dconv,zw(kb-1))  ! maximum depth of convection
        if (kt .eq. maxk) then
          kven = min(kven,kb)   ! deepest level of ventilation from the surface
        endif

        kt = kb - 1

      else
        kt = kb
      endif

      ! continue the search for other unstable regions
      kb = kt - 1

    enddo   ! while (kt > kbo)

    dven = zw(kven-1) ! maximum depth of ventilation from the surface

   return
  
  end subroutine convection


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c o s h u f f l e
  !   Purpose    :  convection scheme from Bern3D after Müller et al
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine coshuffle(l_trans_tracers,kven,ts,rho,kbo)

    implicit none

    logical, dimension(:), intent(in) :: l_trans_tracers
    integer, intent(inout) :: kven  ! maximum integer depth of convection
    real(wp), dimension(:,:), intent(inout) :: ts
    real(wp), dimension(:), intent(inout) :: rho
    integer, intent(in) :: kbo

    integer :: k, k0, l, ipass, maxpass, n_tracers_mix
    real(wp) :: ts_top, rho_top, rho_k


    if (l_conv_shuffle_passive) then
      n_tracers_mix = n_tracers_tot
    else
      n_tracers_mix = 2 ! only active tracers (temperature and salinity)
    endif

    maxpass = maxk

    if (kbo.lt.maxk) then

      k0 = 0
      ipass = 0
      do while (k0.lt.maxk.and.ipass.lt.maxpass)
        ipass = ipass + 1

        k = maxk-1    ! start from second layer
        rho_top =  eos(ts(maxk,1),ts(maxk,2),zw(k))
        rho_k =  eos(ts(k,1),ts(k,2),zw(k))

        ! find level of neutral density
        do while((rho_top.gt.rho_k).and.(k.gt.kbo))
          k=k-1
          rho_top = eos(ts(maxk,1),ts(maxk,2),zw(k))   ! density of top water referenced to level k
          rho_k = eos(ts(k,1),ts(k,2),zw(k))
        enddo

        k0 = k+1
        if (k0.lt.maxk) then  ! unstable 
          do l = 1,n_tracers_mix
            if (l_trans_tracers(l)) then
              ! mix tracers
              ts_top = ts(maxk,l)
              do k = maxk,k0+1,-1
                ts(k,l) = ((dz(k)-dz(maxk))*ts(k,l)+dz(maxk)*ts(k-1,l))*rdz(k)
              enddo
              ts(k0,l) = ((dz(k0)-dz(maxk))*ts(k0,l)+dz(maxk) *ts_top)*rdz(k0)
            endif
          enddo
          do k = k0,maxk
            ! update density
            rho(k) = eos(ts(k,1),ts(k,2),zro(k))
          enddo
          kven = min(kven,k0)
        endif
      enddo

    endif

   return

  end subroutine coshuffle

  subroutine check_stability(unstable,rho,kbo,kto)

    implicit none

    logical,  intent(OUT) :: unstable
    real(wp), intent(IN)  :: rho(:)
    integer,  intent(IN)  :: kbo
    integer,  intent(IN)  :: kto 

    ! Local variables
    integer :: k, nz
    real(wp) :: mindens

    nz = size(rho)

    ! Start at the base
    mindens  = rho(kbo)
    unstable = .FALSE.
    do k = kbo+1, kto 
      if (rho(k) .gt. mindens) then
        unstable = .TRUE.
        exit 
      else
        mindens = rho(k)
      end if
    end do

    return

  end subroutine check_stability

end module convection_mod

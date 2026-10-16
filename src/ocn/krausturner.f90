!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : k r a u s t u r n e r _ m o d
!
!  Purpose : Kraus-Turner mixed layer scheme
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Kevin Oliver, Neil R. Edwards and Matteo Willeit
!
! This file is part of CLIMBER-X.
!
! This file was ported from the original c-GOLDSTEIN model,
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
module krausturner_mod

  ! Simplest Kraus-Turner mixed layer scheme, carried out after
  ! buoyancy-driven convection. The PE release from convection is
  ! calculated in convection.f90 (to be multiplied by an efficiency here).
  ! The KE input from wind is in the field kewind, pre-multiplied
  ! by an efficiency, which is not a function of mixed layer depth
  ! (mld). Internal KE sources are neglected. It is not necessary
  ! to determine mld, because the scheme works by partially mixing
  ! the last layer. However mld is output each time-step as a
  ! diagnostic

  use precision, only : wp
  use ocn_grid, only : maxk, zw, dz, dzg, rdzg, z2dzg
  use ocn_params, only : mlddec, mlddecd, n_tracers_tot
  use constants, only : g
  use eos_mod

  implicit none

  private
  public :: krausturner

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  k r a u s t u r n e r
  !   Purpose    :  Simplest Kraus-Turner mixed layer scheme
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine krausturner(pebuoy,ketau,k1, &
                         ts, &
                         mld,mldk)

    implicit none

    integer, intent(in) :: k1
    real(wp), intent(in) :: pebuoy      ! pe energy input
    real(wp), intent(in) :: ketau       ! ke energy input

    real(wp), intent(inout) :: ts(:,:)  ! tracers including T and S, ts(1:maxk,1:n_tracers_tot)

    integer, intent(out) :: mldk    ! index of layer containing mld
    real(wp), intent(out) :: mld    ! mixed layer depth (-ve value)

    integer :: k, l, n, partmix
    real(wp) :: em, empe, emke, eneed, emr
    real(wp) :: rhou, rhol, rhomix
    real(wp) :: smix, tmix, qmix(1:n_tracers_tot)
    real(wp) :: mlda, mldb, mldtadd


    ! initialise variables
    k = maxk
    empe = pebuoy           ! potential energy from buoyancy
    emke = ketau*mlddec(k)  ! kinetic energy from wind
    em = empe + emke        ! total energy available for mixing
    eneed = 0._wp
    tmix = ts(k,1)
    smix = ts(k,2)
    do l=1,n_tracers_tot
      qmix(l) = ts(k,l)
    enddo
    partmix = 1

    ! completely mix layers down to layer to be only partially mixed is found
    do while (eneed.lt.em .and. em.gt.0._wp)

      if (k.lt.maxk) then
        ! apply mixing from last step if complete
        qmix(1) = tmix
        qmix(2) = smix
        do l=3,n_tracers_tot
          qmix(l) = rdzg(maxk,k)*(qmix(l)*dzg(maxk,k+1) + ts(k,l)*dz(k))
        enddo
      endif

      if (k.eq.k1) then  ! bottom layer

        ! Mixed layer reaches bottom
        mldk = k
        mld = zw(k-1)
        partmix = 0
        em = -1.e-8_wp
        ! mix whole water column
        do l=1,n_tracers_tot
          do n=k,maxk
            ts(n,l) = qmix(l)
          enddo
        enddo

      else  ! not bottom layer

        ! Try to completely mix next layer
        k = k-1
        emr = (em - eneed)/em
        empe = empe*emr
        emke = emke*emr*mlddecd(k)
        em = empe + emke
        ! density of layer above
        rhou = eos(tmix,smix,zw(k))
        ! density of current layer
        rhol = eos(ts(k,1),ts(k,2),zw(k))
        ! compute density of layers if they would be mixed
        tmix = rdzg(maxk,k)*(tmix*dzg(maxk,k+1) + ts(k,1)*dz(k)) 
        smix = rdzg(maxk,k)*(smix*dzg(maxk,k+1) + ts(k,2)*dz(k))
        rhomix = eos(tmix,smix,zw(k))
        ! Energy that would be needed to completely mix layer
        eneed = 0.5_wp*g*(z2dzg(maxk,k+1)*rhou + z2dzg(k,k)*rhol - z2dzg(maxk,k)*rhomix)

      endif

    enddo

    if ((partmix.eq.1).and.(em.gt.0._wp)) then
      ! Partially mix next layer. This scheme is taken from Unified
      ! Model documentation paper No41. Ocean model mixed layer
      ! formulation (S. J. Foreman, 17 Sept 1990). Model Vers <2.0.
      ! (24) is the key equation, but it should read
      ! a=(h_(n-1)/h_n)*(M/E_mix).
      mldk = k ! k is now the index of the layer that is going to be only partially mixed 
      mlda = dzg(maxk,k+1)*rdzg(maxk,k)*em/eneed
      mldb = (dzg(maxk,k)*rdzg(maxk,k+1)-1._wp)*mlda
      do l=1,n_tracers_tot
        ts(maxk,l) = (1._wp-mldb)*qmix(l)+mldb*ts(k,l)
        do n=k+1,maxk-1
          ts(n,l) = ts(maxk,l)
        enddo
        ts(k,l) = mlda*qmix(l)+(1._wp-mlda)*ts(k,l)
      enddo
      ! Actual mixed layer depth. The key equation is (26), but 
      ! again the eqn is wrong (see above ref). It should read
      ! d = 2M / (g h_(n-1) (rho_n - rho_mix) )
      if (rhol.eq.rhou) then
        mldtadd = 0._wp ! avoid division by zero
      else
        mldtadd = 2._wp*em/(g*zw(k)*(rhol-rhou))
      endif
      if (-mldtadd.gt.(dzg(k,k)+1.e-3_wp)) then
        print *
        print*,'Warning: inconsistent result in mixed layer scheme partmixed layer: see krausturner.f'
        print *,-mldtadd,dzg(k,k),k,em,eneed
        print *,'rhol,rhou',rhol,rhou,rhol-rhou
        !print *,'mlda,mldb',mlda,mldb
        !print *,'dzg(maxk,k+1),rdzg(maxk,k)',dzg(maxk,k+1),rdzg(maxk,k)
        !print *,'dzg(maxk,k),rdzg(maxk,k+1)',dzg(maxk,k),rdzg(maxk,k+1)
        mldtadd=-dzg(k,k)
        ! Shouldn't happen, delete this line when confident it doesn't
      endif
      mld = zw(k)+mldtadd
    elseif ((partmix.eq.1).AND.(k.lt.maxk)) then
      mldk = k+1
      mld = zw(k) 
    endif 

    return

  end subroutine krausturner

end module krausturner_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : u 2 d _ m o d
!
!  Purpose : PBL and surface wind
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module u2d_mod

  use atm_params, only : wp
  use constants, only : g
  use atm_params, only : ra, i_kata_wind, h_kata, i_ugb_psi
  use atm_grid, only : im, imc, jm, jmc, nm, dxt, dxu, dy, aim
  use atm_grid, only : fcorg, fcorgu, fcorta, fcorua, signf
  !$use omp_lib

  implicit none
  
  private
  public :: u2d, usur

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u 2 d
  !   Purpose    :  computation of geostrophic and ageostrophic wind in 
  !                 planetary boundary layer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine u2d(slp, sin_cos_acbar, &
      ugb, vgb, ugbu, vgbv, psi_g, uab, vab)

    implicit none

    real(wp), intent(in ) :: slp(:,:)
    real(wp), intent(in ) :: sin_cos_acbar(:,:)

    real(wp), intent(out) :: ugb(:,:)
    real(wp), intent(out) :: vgb(:,:)
    real(wp), intent(out) :: ugbu(:,:)
    real(wp), intent(out) :: vgbv(:,:)
    real(wp), intent(out) :: psi_g(:,:)
    real(wp), intent(out) :: uab(:,:)
    real(wp), intent(out) :: vab(:,:)

    integer :: i, j, ipl, imi, jpl, jmi
    real(wp) :: dpdx, dpdy, dpdxa, dpdya, acbarb
    real(wp) :: slpz(jm), slpa(im,jm), psibar(jmc), ubarz(jm), slpc


    !$omp parallel do collapse(2) private(i,j,ipl,imi,jpl,jmi,dpdx,dpdy,dpdxa,dpdya,acbarb)
    do j=1,jm
      do i=1,im 

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1
        jpl = min(jm,j+1)
        jmi = max(1,j-1)

        !------------------------------------------------------
        ! PBL geostrophic wind components in T-points

        ! Horizontal SLP gradient in T-points      
        dpdx = 0.5_wp*(slp(ipl,j)-slp(imi,j))/dxt(j) 
        dpdy = 0.5_wp*(slp(i,jmi)-slp(i,jpl))/dy  

        ! fcorg is 1/fcort for i_fcorg=0 and the regular f/(f**2+fcormin**2) for i_fcorg=1, which
        ! takes ugb through zero at the equator instead of reversing it at full amplitude.
        ! See the fcorg block in atm_grid.f90.
        ugb(i,j) = -dpdy*fcorg(j)/ra
        vgb(i,j) =  dpdx*fcorg(j)/ra

        !------------------------------------------------------
        ! Ageostrophic wind components in PBL (on U-points)

        dpdxa = (slp(i,j)-slp(imi,j))/dxt(j)
        acbarb = 0.5_wp*(sin_cos_acbar(imi,j)+sin_cos_acbar(i,j))
        uab(i,j) = -1._wp/(fcorta(j)*ra)*(dpdxa*acbarb)

        if (j.gt.1) then
          dpdya = (slp(i,j-1)-slp(i,j))/dy
          acbarb = 0.5_wp*(sin_cos_acbar(i,j-1)+sin_cos_acbar(i,j))
          vab(i,j) = -1._wp/(fcorua(j)*ra)*(dpdya*acbarb)
        endif

        if (j.eq.jm) then
          vab(i,1)   = 0._wp
          vab(i,jmc) = 0._wp
        endif

        ! periodic boundary conditions
        if (i.eq.1) then
          uab(imc,j) = uab(1,j)
        endif

      enddo
    enddo 
    !$omp end parallel do

    !------------------------------------------------------
    ! Barotropic geostrophic wind on the faces, which is what the mass flux in
    ! u3d.f90 actually needs.  Building it here rather than averaging ugb/vgb
    ! inside u3d is what makes the streamfunction option possible: the curl of a
    ! corner streamfunction lives on the faces, and averaging it to T-points and
    ! back would apply a 1-2-1 filter and destroy the non-divergence.
    !------------------------------------------------------

    if (i_ugb_psi.eq.0) then

      ! Exactly the averaging u3d used to do internally, so this branch reproduces
      ! the previous code bit for bit.
      do j=1,jm
        do i=1,im
          imi = modulo(i - 2, im) + 1
          ugbu(i,j) = 0.5_wp*(ugb(imi,j)+ugb(i,j))
        enddo
        ugbu(imc,j) = ugbu(1,j)
      enddo
      vgbv(:,1)   = 0._wp
      vgbv(:,jmc) = 0._wp
      do j=2,jm
        do i=1,im
          vgbv(i,j) = 0.5_wp*(vgb(i,j-1)+vgb(i,j))
        enddo
      enddo
      psi_g(:,:) = 0._wp

    else if (i_ugb_psi.eq.1) then

      !----------------------------------------------------
      ! Streamfunction form.
      !
      ! psi is carried on the cell corners, and the face winds are its discrete
      ! curl in exactly the pairing that annihilates the column convergence
      ! stencil used in u3d.f90,
      !     conv(i,j) = Fx(i,j)-Fx(i+1,j) + Fy(i,j+1)-Fy(i,j) ,
      ! namely  Fx(i,j) ~ psi(i,j+1)-psi(i,j)  and  Fy(i,j) ~ psi(i+1,j)-psi(i,j).
      ! Substituting the two into conv cancels all eight terms identically, so the
      ! barotropic geostrophic velocity field carries no divergence at all - no
      ! beta term, and no contribution from the meridional variation of fcorgu or
      ! of any other latitude factor folded into psi.
      !
      ! The zonal mean is separated out and reinstated exactly.  This is essential:
      ! psi = fcorgu*p/ra with the FULL slp would add the term -p*d(fcorgu)/dy to
      ! the zonal wind, and with p ~ 1e5 Pa that is of order (p/ra)*beta/f**2, some
      ! 120 m/s at 45 deg.  Only the azonal slp, of order 1e3 Pa, may pass through
      ! fcorgu.  psibar is instead integrated straight from the zonal mean zonal
      ! wind of the i_ugb_psi=0 form, so that zonal mean is reproduced to round-off;
      ! being independent of i it produces no meridional wind and a zonal flux that
      ! telescopes away, so it adds no divergence either.
      !----------------------------------------------------

      ! zonal mean and azonal sea level pressure
      do j=1,jm
        slpz(j) = 0._wp
        do i=1,im
          slpz(j) = slpz(j) + slp(i,j)
        enddo
        slpz(j) = slpz(j)*aim
        do i=1,im
          slpa(i,j) = slp(i,j) - slpz(j)
        enddo
      enddo

      ! zonal mean zonal geostrophic wind, identical to the zonal mean of the
      ! ugb built above because the meridional difference is linear in slp
      do j=1,jm
        jpl = min(jm,j+1)
        jmi = max(1,j-1)
        ubarz(j) = -0.5_wp*(slpz(jmi)-slpz(jpl))/dy * fcorg(j)/ra
      enddo

      ! streamfunction of that zonal mean, from u = (psi(j+1)-psi(j))/dy
      psibar(1) = 0._wp
      do j=1,jm
        psibar(j+1) = psibar(j) + ubarz(j)*dy
      enddo

      ! full streamfunction on the corners. Both polar rows are held zonally
      ! constant, which is the no-flux-through-the-pole condition: the meridional
      ! flux there is psi(i+1,j)-psi(i,j) and dxu vanishes anyway.
      psi_g(:,1)   = psibar(1)
      psi_g(:,jmc) = psibar(jmc)
      do j=2,jm
        do i=1,im
          imi = modulo(i - 2, im) + 1
          slpc = 0.25_wp*(slpa(imi,j-1)+slpa(i,j-1)+slpa(imi,j)+slpa(i,j))
          psi_g(i,j) = psibar(j) + fcorgu(j)*slpc/ra
        enddo
      enddo

      ! face winds as the discrete curl
      do j=1,jm
        do i=1,im
          ugbu(i,j) = (psi_g(i,j+1)-psi_g(i,j))/dy
        enddo
        ugbu(imc,j) = ugbu(1,j)
      enddo
      vgbv(:,1)   = 0._wp
      vgbv(:,jmc) = 0._wp
      do j=2,jm
        do i=1,im
          ipl = modulo(i,im) + 1
          vgbv(i,j) = (psi_g(ipl,j)-psi_g(i,j))/dxu(j)
        enddo
      enddo

      ! T-point values, overwriting the f^-1*grad(p) form, so that the surface wind
      ! in usur and every diagnostic sees the same field the transport does
      do j=1,jm
        do i=1,im
          ipl = modulo(i,im) + 1
          ugb(i,j) = 0.5_wp*(ugbu(i,j)+ugbu(ipl,j))
          vgb(i,j) = 0.5_wp*(vgbv(i,j)+vgbv(i,j+1))
        enddo
      enddo

    else

      stop 'i_ugb_psi'

    endif

    return

  end subroutine u2d


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u s u r
  !   Purpose    :  computation of near-surface wind components
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine usur(ugb, vgb, epsa, cos_acbar, sin_acbar, t2a, tskina, cd0a, slope_x, slope_y, &
      usk, vsk, &
      us, vs)

    implicit none

    real(wp), intent(in   ) :: ugb(:,:)
    real(wp), intent(in   ) :: vgb(:,:)
    real(wp), intent(in   ) :: epsa(:,:,:)
    real(wp), intent(in   ) :: cos_acbar(:,:,:)
    real(wp), intent(in   ) :: sin_acbar(:,:,:)
    real(wp), intent(in   ) :: t2a(:,:)
    real(wp), intent(in   ) :: tskina(:,:)
    real(wp), intent(in   ) :: cd0a(:,:)
    real(wp), intent(in   ) :: slope_x(:,:)
    real(wp), intent(in   ) :: slope_y(:,:)

    real(wp), intent(inout) :: usk(:,:)
    real(wp), intent(inout) :: vsk(:,:)

    real(wp), intent(out  ) :: us(:,:,:)
    real(wp), intent(out  ) :: vs(:,:,:)

    integer :: i, j, n, ipl, imi
    real(wp) :: cd, theta0, theta1, dtheta

    real(wp), dimension(im,jm)  :: uk
    real(wp), dimension(im,jmc) :: vk


    !$omp parallel do collapse(2) private(i,j,n,ipl,imi,cd,dtheta,theta0,theta1)
    do j=1,jm
      do i=1,im 

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1

        if (j.gt.1 .and. j.lt.jm) then

          !------------------------------------------------------
          ! Near surface wind components in T-points
          do n=1,nm
            us(i,j,n) = epsa(i,j,n) * (ugb(i,j)*cos_acbar(i,j,n) - signf(j)*vgb(i,j)*sin_acbar(i,j,n))
            vs(i,j,n) = epsa(i,j,n) * (vgb(i,j)*cos_acbar(i,j,n) + signf(j)*ugb(i,j)*sin_acbar(i,j,n))
          enddo

        endif

        !------------------------------------------------------
        ! katabatic surface wind 

        if (i_kata_wind.ne.0) then

          ! compute katabatic winds from balance of bouyancy force and friction, ignoring Coriolis and background pressure gradient
          ! see e.g. Fedorovich and Shapiro 2009, Prandtl 1942 model

          ! zonal component on u-grid
          cd = 0.5_wp*(cd0a(i,j)+cd0a(imi,j))
          if (i_kata_wind.eq.1) then
            theta0 = 0.5_wp*(tskina(i,j)+tskina(imi,j))
            theta1 = 0.5_wp*(t2a(i,j)+t2a(imi,j))
          else if (i_kata_wind.eq.2) then
            ! use upstream values
            if (slope_x(i,j).gt.0._wp) then
              theta0 = tskina(i,j)
              theta1 = t2a(i,j)
            else
              theta0 = tskina(imi,j)
              theta1 = t2a(imi,j)
            endif
          endif
          dtheta = 2._wp*(theta1-theta0)
          if (dtheta.gt.0._wp) then
            uk(i,j) = sign(sqrt(g*h_kata/cd * dtheta/theta0 * abs(slope_x(i,j))), -slope_x(i,j))
          else
            uk(i,j) = 0._wp
          endif

          ! meridional component on v-grid
          if (j.gt.1) then

            cd = 0.5_wp*(cd0a(i,j-1)+cd0a(i,j))
            if (i_kata_wind.eq.1) then
              theta0 = 0.5_wp*(tskina(i,j-1)+tskina(i,j))
              theta1 = 0.5_wp*(t2a(i,j-1)+t2a(i,j))
            else if (i_kata_wind.eq.2) then
              ! use upstream values
              if (slope_y(i,j).gt.0._wp) then
                theta0 = tskina(i,j-1)
                theta1 = t2a(i,j-1)
              else
                theta0 = tskina(i,j)
                theta1 = t2a(i,j)
              endif
            endif
            dtheta = 2._wp*(theta1-theta0) 
            if (dtheta.gt.0._wp) then
              vk(i,j) = sign(sqrt(g*h_kata/cd * dtheta/theta0 * abs(slope_y(i,j))), -slope_y(i,j))
            else
              vk(i,j) = 0._wp
            endif

          else

            uk(i,j) = 0._wp
            vk(i,j) = 0._wp
            vk(i,jmc) = 0._wp

          endif

        else 

          uk(i,j) = 0._wp
          vk(i,j) = 0._wp
          vk(i,jmc) = 0._wp

        endif

      enddo
    enddo 
    !$omp end parallel do

    ! Wind near the poles
    us(:,1,:)  = 0.5_wp*us(:,2,:)       
    us(:,jm,:) = 0.5_wp*us(:,jm-1,:) 
    vs(:,1,:)  = 0.5_wp*vs(:,2,:)       
    vs(:,jm,:) = 0.5_wp*vs(:,jm-1,:)

    !$omp parallel do collapse(2) private(i,j,n,ipl,imi)
    do j=1,jm
      do i=1,im 

        ipl = modulo(i,im) + 1
        imi = modulo(i - 2, im) + 1

        ! interpolate to t-points and relax in time
        usk(i,j) = 0.1_wp * 0.5_wp*(uk(i,j)+uk(ipl,j)) + 0.9_wp*usk(i,j)
        vsk(i,j) = 0.1_wp * 0.5_wp*(vk(i,j)+vk(i,j+1)) + 0.9_wp*vsk(i,j)

        ! add katabatic wind to surface wind
        do n=1,nm
          us(i,j,n) = us(i,j,n) + usk(i,j)
          vs(i,j,n) = vs(i,j,n) + vsk(i,j)
        enddo

      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine usur

end module u2d_mod

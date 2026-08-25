!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : a d v _ f c t _ m o d
!
!  Purpose : flux corrected transport (Zalesak 1979) for the horizontal
!            advection of atmospheric tracers
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
!
! The scheme mirrors advection_fct in src/ocn/advection.f90, with two
! differences that the atmospheric grid forces:
!
! 1) The provisional low order solution is formed in ADVECTIVE, not flux,
!    form. The ocean velocity field is non-divergent cell by cell, so there
!    the donor cell update is automatically a convex combination of the
!    neighbouring tracer values, which is what Zalesak's bounds rest on.
!    The atmospheric mass fluxes are corrected for column mass balance only
!    (the psi/psi_topo Poisson correction in u3d), so LEVEL BY LEVEL they
!    are strongly divergent - the vertical mass flux w3 is precisely that
!    divergence accumulated over the column. Subtracting x*div(F_mass)
!    restores the convex combination:
!        xt = x + dt/m * sum_faces F_face*(x_donor - x)
!    which is bounded by the neighbouring values whenever the outgoing mass
!    per step does not exceed the cell mass. div(F_mass) does not depend on
!    the interpolation, so it is identical for the low and the high order
!    leg and the antidiffusive flux is unaffected.
!
! 2) What is returned is the LIMITER RATIO r in [0,1] per face, not the
!    limited antidiffusive flux. adifa multiplies the same face value by
!    three different mass fluxes (fax, fax_psi, fax_psi_topo) and relies on
!    the three products adding up exactly; handing back a face value
!        xface = xup + r*(0.5*(x_i+x_im1) - xup)
!    keeps that property by construction. It also avoids dividing by the
!    mass flux, which vanishes in every layer below the topography.
!
! There is no vertical term: adifa advects horizontally only, the vertical
! redistribution being implicit in the profile reconstruction of vesta.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module adv_fct_mod

  use atm_params, only : wp
  use atm_params, only : tstep
  use atm_grid, only : im, imc, jm, jmc, sqr
  !$ use omp_lib

  implicit none

  private
  public :: fct_ratios
  public :: fct_acc, fct_diag_reset

  ! --------------------------------------------------------------------
  ! Limiter activity diagnostics, accumulated per latitude when fct_ratios
  ! is called with idiag > 0. Slot 1 is dry static energy, slot 2 is water.
  !   1 : sum over faces and levels of |antidiffusive flux| * r
  !   2 : sum over faces and levels of |antidiffusive flux|
  !   3 : cells whose limiter was zeroed by the CFL fallback
  !   4 : cells tested
  ! The ratio (1)/(2) is the fraction of the available antidiffusion that
  ! actually gets through, that is, how much of the upstream numerical
  ! diffusion the scheme really removes. It is 1 for pure centred differencing
  ! and 0 when the scheme has degenerated back to plain upstream.
  ! --------------------------------------------------------------------
  real(wp) :: fct_acc(4,jm,2) = 0._wp

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f c t _ d i a g _ r e s e t
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fct_diag_reset

    implicit none

    fct_acc(:,:,:) = 0._wp

    return

  end subroutine fct_diag_reset

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f c t _ r a t i o s
  !   Purpose    :  antidiffusive flux limiter ratios for one tracer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fct_ratios(nlev, fmx, fmy, dm, x3, rx, ry, idiag)

    implicit none

    integer,  intent(in   ) :: nlev            ! number of levels to work on
    real(wp), intent(in   ) :: fmx(:,:,:)      ! (imc,jm,nlev) zonal mass flux on u-faces, kg/s, positive eastward
    real(wp), intent(in   ) :: fmy(:,:,:)      ! (im,jmc,nlev) meridional mass flux on v-faces, kg/s, positive southward (adifa convention)
    real(wp), intent(in   ) :: dm(:,:,:)       ! (im,jm,nlev) cell mass per unit area, kg/m2
    real(wp), intent(in   ) :: x3(:,:,:)       ! (im,jm,nlev) tracer at T-points

    real(wp), intent(  out) :: rx(:,:,:)       ! (imc,jm,nlev) limiter ratio on u-faces, [0,1]
    real(wp), intent(  out) :: ry(:,:,:)       ! (im,jmc,nlev) limiter ratio on v-faces, [0,1]

    integer,  intent(in   ) :: idiag           ! accumulate limiter diagnostics into fct_acc(:,:,idiag); 0 = do not

    integer :: i, j, k
    integer :: im1, im2, ip1, jm1, jm2, jp1
    real(wp) :: f, xup, a, m, netm, netf, fin, xmax, xmin, pp, pm, qp, qm

    ! level local work space
    real(wp) :: flx(imc,jm), fly(im,jmc)       ! low order (upstream) tracer fluxes
    real(wp) :: ax(imc,jm),  ay(im,jmc)        ! antidiffusive fluxes, prelimited
    logical  :: kx(imc,jm),  ky(im,jmc)        ! face kept by the prelimiter?
    real(wp) :: xt(im,jm)                      ! provisional low order solution
    real(wp) :: xa(im,jm), xb(im,jm)           ! local bounds
    real(wp) :: rp(im,jm), rn(im,jm)           ! Zalesak R+ / R-
    logical  :: lok(im,jm)                     ! cell has mass and satisfies the CFL condition

    real(wp), parameter :: big = huge(1._wp)

    ! per level partial sums, reduced after the loop; each k writes its own slice
    real(wp) :: dacc(4,jm,nlev)


    !$omp parallel do private(i,j,k,im1,im2,ip1,jm1,jm2,jp1,f,xup,a,m,netm,netf,fin,xmax,xmin,pp,pm,qp,qm) &
    !$omp private(flx,fly,ax,ay,kx,ky,xt,xa,xb,rp,rn,lok)
    do k=1,nlev

      !-----------------------------------
      ! low order (upstream) and antidiffusive fluxes
      !-----------------------------------

      ! zonal, face i separates cell i-1 (west) from cell i (east), as in adifa
      do j=1,jm
        do i=1,im
          im1 = modulo(i-2,im) + 1
          f = fmx(i,j,k)
          if (f.gt.0._wp) then
            xup = x3(im1,j,k)
          else
            xup = x3(i,j,k)
          endif
          flx(i,j) = f*xup
          ! centred minus upstream, which reduces to |F|/2 * (x_east - x_west)
          ax(i,j) = 0.5_wp*abs(f)*(x3(i,j,k)-x3(im1,j,k))
        enddo
      enddo
      ! periodic closure: face imc is the same physical face as face 1
      flx(imc,:) = flx(1,:)
      ax(imc,:)  = ax(1,:)

      ! meridional, face j separates cell j-1 (south) from cell j (north);
      ! fmy>0 means the flow goes from cell j to cell j-1
      fly(:,1)   = 0._wp
      fly(:,jmc) = 0._wp
      ay(:,1)    = 0._wp
      ay(:,jmc)  = 0._wp
      do j=2,jm
        do i=1,im
          f = fmy(i,j,k)
          if (f.gt.0._wp) then
            xup = x3(i,j,k)
          else
            xup = x3(i,j-1,k)
          endif
          fly(i,j) = f*xup
          ay(i,j) = 0.5_wp*abs(f)*(x3(i,j-1,k)-x3(i,j,k))
        enddo
      enddo

      ! antidiffusion available on this level, before any limiting; the
      ! denominator of the diagnostic below
      if (idiag.gt.0) then
        dacc(:,:,k) = 0._wp
        do j=1,jm
          do i=1,im
            dacc(2,j,k) = dacc(2,j,k) + abs(ax(i,j))
            if (j.gt.1) dacc(2,j,k) = dacc(2,j,k) + abs(ay(i,j))
          enddo
        enddo
      endif

      !-----------------------------------
      ! provisional low order solution, in advective form (see header)
      !-----------------------------------
      do j=1,jm
        do i=1,im
          m = dm(i,j,k)*sqr(i,j)   ! cell mass, kg
          if (m.gt.0._wp) then
            ! mass convergence and low order tracer convergence, same sign convention as adifa
            netm = fmx(i,j,k)-fmx(i+1,j,k) + fmy(i,j+1,k)-fmy(i,j,k)
            netf = flx(i,j)    -flx(i+1,j) + fly(i,j+1)  -fly(i,j)
            xt(i,j) = x3(i,j,k) + tstep/m*(netf - x3(i,j,k)*netm)
            ! In the advective form only the INFLOW faces carry a non zero
            ! (x_donor - x): on an outflow face the donor is the cell itself and
            ! the term drops out. The update is therefore
            !   xt = x*(1 - dt/m*fin) + dt/m * sum_inflow F*x_neighbour
            ! which is a convex combination of the cell and its neighbours, and so
            ! bounded by them, exactly while dt*fin <= m. The outflow does not enter.
            ! Testing the outflow instead would switch the limiter off wherever a
            ! level is strongly divergent - and the horizontal divergence of a single
            ! level IS the vertical mass flux, so that is most of the atmosphere.
            fin = max(0._wp,fmx(i,j,k)) + max(0._wp,-fmx(i+1,j,k)) &
                + max(0._wp,fmy(i,j+1,k)) + max(0._wp,-fmy(i,j,k))
            lok(i,j) = tstep*fin .le. m
          else
            ! layer entirely below the topography: nothing to advect, fall back to upstream
            xt(i,j) = x3(i,j,k)
            lok(i,j) = .false.
          endif
        enddo
      enddo

      !-----------------------------------
      ! prelimiting, eq. 14' of Zalesak 1979: drop antidiffusive fluxes that
      ! would push the provisional solution away from, instead of towards,
      ! the sharper profile
      !-----------------------------------
      do j=1,jm
        do i=1,im
          im1 = modulo(i-2,im) + 1
          im2 = modulo(i-3,im) + 1
          ip1 = modulo(i,im) + 1
          a = ax(i,j)
          ! positive a flows into cell i and out of cell im1
          kx(i,j) = .not.( a*(xt(i,j)-xt(im1,j)).lt.0._wp &
                     .and. ( a*(xt(ip1,j)-xt(i,j)).lt.0._wp .or. a*(xt(im1,j)-xt(im2,j)).lt.0._wp ) )
          if (.not.kx(i,j)) ax(i,j) = 0._wp
        enddo
      enddo
      kx(imc,:) = kx(1,:)
      ax(imc,:) = ax(1,:)

      ky(:,1)   = .false.
      ky(:,jmc) = .false.
      do j=2,jm
        do i=1,im
          jm2 = max(1,j-2)
          jp1 = min(jm,j+1)
          a = ay(i,j)
          ! positive a flows into cell j-1 and out of cell j
          ky(i,j) = .not.( a*(xt(i,j-1)-xt(i,j)).lt.0._wp &
                     .and. ( a*(xt(i,jm2)-xt(i,j-1)).lt.0._wp .or. a*(xt(i,j)-xt(i,jp1)).lt.0._wp ) )
          if (.not.ky(i,j)) ay(i,j) = 0._wp
        enddo
      enddo

      !-----------------------------------
      ! local bounds; cells without mass are pushed out of the max/min search
      !-----------------------------------
      do j=1,jm
        do i=1,im
          if (dm(i,j,k)*sqr(i,j).gt.0._wp) then
            xa(i,j) = max(x3(i,j,k),xt(i,j))
            xb(i,j) = min(x3(i,j,k),xt(i,j))
          else
            xa(i,j) = -big
            xb(i,j) =  big
          endif
        enddo
      enddo

      !-----------------------------------
      ! Zalesak R+ / R-
      !-----------------------------------
      do j=1,jm
        do i=1,im
          if (.not.lok(i,j)) then
            rp(i,j) = 0._wp
            rn(i,j) = 0._wp
            cycle
          endif
          im1 = modulo(i-2,im) + 1
          ip1 = modulo(i,im) + 1
          jm1 = max(1,j-1)
          jp1 = min(jm,j+1)
          xmax = max(xa(i,j),xa(im1,j),xa(ip1,j),xa(i,jm1),xa(i,jp1))
          xmin = min(xb(i,j),xb(im1,j),xb(ip1,j),xb(i,jm1),xb(i,jp1))
          ! antidiffusive inflow / outflow, same face-to-cell mapping as the convergence.
          ! Scaled by tstep because the fluxes here are rates (kg/s times tracer) while
          ! q below is an amount (kg times tracer); the ocean version carries the dt in
          ! the fluxes themselves instead.
          pp = tstep*( max(0._wp,ax(i,j)) - min(0._wp,ax(i+1,j)) + max(0._wp,ay(i,j+1)) - min(0._wp,ay(i,j)) )
          pm = tstep*( max(0._wp,ax(i+1,j)) - min(0._wp,ax(i,j)) + max(0._wp,ay(i,j)) - min(0._wp,ay(i,j+1)) )
          m = dm(i,j,k)*sqr(i,j)
          qp = (xmax-xt(i,j))*m
          qm = (xt(i,j)-xmin)*m
          ! no dimensional threshold on p: the ratio saturates at 1 anyway, and the
          ! tracers here span kg/kg to K
          if (pp.gt.0._wp) then
            rp(i,j) = min(1._wp,max(0._wp,qp)/pp)
          else
            rp(i,j) = 0._wp
          endif
          if (pm.gt.0._wp) then
            rn(i,j) = min(1._wp,max(0._wp,qm)/pm)
          else
            rn(i,j) = 0._wp
          endif
        enddo
      enddo

      !-----------------------------------
      ! face ratios: the smaller of the receiving cell's room to rise and the
      ! donating cell's room to fall
      !-----------------------------------
      do j=1,jm
        do i=1,im
          im1 = modulo(i-2,im) + 1
          if (.not.kx(i,j)) then
            rx(i,j,k) = 0._wp
          else if (ax(i,j).ge.0._wp) then
            rx(i,j,k) = min(rp(i,j),rn(im1,j))
          else
            rx(i,j,k) = min(rp(im1,j),rn(i,j))
          endif
        enddo
      enddo
      rx(imc,:,k) = rx(1,:,k)

      ry(:,1,k)   = 0._wp
      ry(:,jmc,k) = 0._wp
      do j=2,jm
        do i=1,im
          if (.not.ky(i,j)) then
            ry(i,j,k) = 0._wp
          else if (ay(i,j).ge.0._wp) then
            ry(i,j,k) = min(rp(i,j-1),rn(i,j))
          else
            ry(i,j,k) = min(rp(i,j),rn(i,j-1))
          endif
        enddo
      enddo

      ! antidiffusion that actually got through, and the CFL fallback count.
      ! ax/ay have been zeroed on prelimited faces, where r is zero as well, so
      ! the product is the same either way and the denominator above stays the
      ! full amount that was on offer.
      if (idiag.gt.0) then
        do j=1,jm
          do i=1,im
            dacc(1,j,k) = dacc(1,j,k) + abs(ax(i,j))*rx(i,j,k)
            if (j.gt.1) dacc(1,j,k) = dacc(1,j,k) + abs(ay(i,j))*ry(i,j,k)
            if (dm(i,j,k)*sqr(i,j).gt.0._wp) then
              dacc(4,j,k) = dacc(4,j,k) + 1._wp
              if (.not.lok(i,j)) dacc(3,j,k) = dacc(3,j,k) + 1._wp
            endif
          enddo
        enddo
      endif

    enddo
    !$omp end parallel do

    if (idiag.gt.0) then
      do k=1,nlev
        do j=1,jm
          fct_acc(:,j,idiag) = fct_acc(:,j,idiag) + dacc(:,j,k)
        enddo
      enddo
    endif

    return

  end subroutine fct_ratios

end module adv_fct_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i n v e r t _ m o d
!
!  Purpose : invert matrix for barotropic streamfunction
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
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
module invert_mod

  use precision, only : wp
  use ocn_grid, only : maxi, maxj, maxk, k1, c, cv, ds, dsv, dphi, rh, R_earth
  use ocn_params, only : fcor, fcorv, drag, drhcor_max, i_cor_form, nbw, qcor

  implicit none

  private
  public :: invert

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s t e n c i l _ c o e f
  !   Purpose    :  coefficients of psi(i+di,j+dj), di,dj = -1..1, of the discrete vorticity equation at psi point (i,j)
  !
  !   The equation is the discrete curl (dual-cell loop with the segment lengths of island.f90 / wind.f90) of
  !     drag*u_b + f x u_b = tau/(rho0 H) + JEBAR - grad p
  !   The drag part is the same for both i_cor_form. The Coriolis part is
  !     i_cor_form = 0 : advective form J(psi,f/H), central differences, drhcor_max limiter. NOT a discrete
  !                      divergence, so no line integrand matches it and the island constants of island.f90
  !                      depend on the integration path
  !     i_cor_form = 1 : flux form, J3 Jacobian d_phi(psi d_s q) - d_s(psi d_phi q) with q = f/H at tracer points;
  !                      a discrete divergence, so the matching path integral in island.f90 is path-independent.
  !                      Only psi averaging, and only the 4 cells of the psi point are used
  !   Both are 5-point stencils. Valid for 1 <= j <= maxj-1.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine stencil_coef(i,j,cf)

    implicit none

    integer, intent(in) :: i, j
    real(wp), dimension(-1:1,-1:1), intent(out) :: cf

    real(wp) :: tv, tv1, rh1, rh2, drh
    real(wp) :: dsq_e, dsq_w, dpq_n, dpq_s, fac_p, fac_s


    cf(:,:) = 0._wp

    ! drag part: curl of drag*u_b, symmetric 5-point stencil
    cf( 0,-1) = drag(1,i,j)*c(j)**2*rh(1,i,j)/(ds(j)*dsv(j)*R_earth**2)
    cf(-1, 0) = drag(2,i,j)*rh(2,i,j)/(cv(j)*dphi*R_earth)**2
    cf( 0, 0) = - (drag(2,i,j)*rh(2,i,j) + drag(2,i+1,j)*rh(2,i+1,j)) &
               /(cv(j)*dphi*R_earth)**2 &
               - (drag(1,i,j)*c(j)**2*rh(1,i,j)/ds(j) + drag(1,i,j+1)*c(j+1)**2*rh(1,i,j+1)/ds(j+1)) &
               /(dsv(j)*R_earth**2)
    cf( 1, 0) = drag(2,i+1,j)*rh(2,i+1,j)/(cv(j)*dphi*R_earth)**2
    cf( 0, 1) = drag(1,i,j+1)*c(j+1)**2*rh(1,i,j+1)/(ds(j+1)*dsv(j)*R_earth**2)

    if (i_cor_form.eq.0) then

      ! Coriolis terms, advective form, limit topography gradient
      rh1 = rh(1,i,j+1)
      rh2 = rh(1,i,j)
      drh = abs((rh1-rh2)/dphi)
      if (drh.gt.drhcor_max) then
        drh = drhcor_max
        rh2 = rh1 - sign(drh,(rh1-rh2)/dphi)*dphi
      endif
      tv  = (fcor(j+1)*rh1 - fcor(j)*rh2) /(2._wp*dphi*dsv(j)*R_earth**2)  ! 1/s/m3
      rh1 = rh(2,i+1,j)
      rh2 = rh(2,i,j)
      drh = abs((rh1-rh2)/dsv(j))
      if (drh.gt.drhcor_max) then
        drh = drhcor_max
        rh2 = rh1 - sign(drh,(rh1-rh2)/dsv(j))*dsv(j)
      endif
      tv1 = (fcorv(j)*rh1  - fcorv(j)*rh2)/(2._wp*dphi*dsv(j)*R_earth**2)

      cf( 0,-1) = cf( 0,-1) + tv1
      cf(-1, 0) = cf(-1, 0) - tv
      cf( 1, 0) = cf( 1, 0) + tv
      cf( 0, 1) = cf( 0, 1) - tv1

    else if (i_cor_form.eq.1) then

      ! Coriolis terms, J3 flux form: J(psi,q) = d_phi(psi d_s q) - d_s(psi d_phi q), q = f/H at tracer points.
      !   d_phi(psi d_s q) = [ (psi(i,j)+psi(i+1,j))/2 * dsq_e - (psi(i-1,j)+psi(i,j))/2 * dsq_w ] / dphi
      !   d_s  (psi d_phi q) = [ (psi(i,j)+psi(i,j+1))/2 * dpq_n - (psi(i,j-1)+psi(i,j))/2 * dpq_s ] / dsv(j)
      ! with the q differences taken on the dual-cell edges, i.e. over the four cells of this psi point only
      dsq_e = (qcor(i+1,j+1) - qcor(i+1,j))/dsv(j)
      dsq_w = (qcor(i  ,j+1) - qcor(i  ,j))/dsv(j)
      dpq_n = (qcor(i+1,j+1) - qcor(i  ,j+1))/dphi
      dpq_s = (qcor(i+1,j  ) - qcor(i  ,j  ))/dphi
      fac_p = 1._wp/(dphi*R_earth**2)
      fac_s = 1._wp/(dsv(j)*R_earth**2)

      cf( 1, 0) = cf( 1, 0) + 0.5_wp*fac_p*dsq_e
      cf(-1, 0) = cf(-1, 0) - 0.5_wp*fac_p*dsq_w
      cf( 0, 1) = cf( 0, 1) - 0.5_wp*fac_s*dpq_n
      cf( 0,-1) = cf( 0,-1) + 0.5_wp*fac_s*dpq_s
      cf( 0, 0) = cf( 0, 0) + 0.5_wp*fac_p*(dsq_e-dsq_w) - 0.5_wp*fac_s*(dpq_n-dpq_s)

    else

      stop 'stencil_coef: i_cor_form must be 0 or 1'

    endif

    return

  end subroutine stencil_coef


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i n v e r t
  !   Purpose    :  invert matrix for barotropic streamfunction
  !
  !   Band storage: gap(k, nbw+1+off) holds the coefficient of psi at row index k+off, with the half band width
  !   nbw = maxi+1 of the 5-point stencil (set in momentum_init).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine invert(gap,ratm)

    implicit none

    real(wp), dimension(:,:), intent(out) :: gap
    real(wp), dimension(:,:), intent(out) :: ratm


    integer i, j, k, l, n, m, im, di, dj, ic
    real(wp) :: rat
    real(wp), dimension(-1:1,-1:1) :: cf


    n = maxi
    m = maxj + 1

    ! Set equation at Psi points, assuming periodic b.c. in i.
    ! Cannot solve at both i=0 and i=maxi as periodicity => would
    ! have singular matrix. At dry points equation is trivial.

    gap(:,:) = 0._wp

    !$omp parallel do private(i,j,k,l,di,dj,ic,cf)
    do i=1,maxi
       do j=0,maxj
          k=i + j*n
          if (max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.maxk) then
             if (j.eq.0) stop 'j==0'
             if (j.eq.maxj) stop 'j==maxj'

             call stencil_coef(i,j,cf)

             ! scatter into the band: column index of psi(i+di,j+dj) with periodic wrap in i
             do dj=-1,1
               do di=-1,1
                 if (cf(di,dj).eq.0._wp) cycle
                 ic = i+di
                 if (ic.lt.1) ic = ic+maxi
                 if (ic.gt.maxi) ic = ic-maxi
                 l = ic + (j+dj)*n
                 gap(k,nbw+1+l-k) = gap(k,nbw+1+l-k) + cf(di,dj)
               enddo
             enddo

          else

             gap(k,nbw+1) = 1._wp 

          endif
       enddo
    enddo
    !$omp end parallel do

    ! now invert the thing (banded LU, half band width nbw)

    do i=1,n*m-1
       im = min(i+nbw,n*m)
       do j=i+1,im
          rat = gap(j,nbw+1-j+i)/gap(i,nbw+1)
          ratm(j,j-i) = rat
          if (rat.ne.0) then
             do k=nbw+1-j+i,2*nbw+1-j+i
                gap(j,k)=gap(j,k) - rat*gap(i,k+j-i)
             enddo
          endif
       enddo
    enddo

   return

  end subroutine invert

end module invert_mod

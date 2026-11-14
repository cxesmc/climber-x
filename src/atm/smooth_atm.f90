!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s m o o t h _ a t m _ m o d
!
!  Purpose : smoothing functions for atmosphere
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
module smooth_atm_mod

  use atm_params, only : wp
  use, intrinsic :: iso_c_binding

  implicit none

  ! include FFTW3 library for Fourier Transform (https://www.fftw.org/), polar_filter
  include 'fftw3.f03'

  private
  public :: smooth2, smooth2_m, smooth2eq, zona, zofil, shapiro2dx, limit_filt_wtr, polar_filter, polar_filter3

  ! FFTW plans of polar_filter, made on the first call
  integer :: n_pf_plan = 0
  type(C_PTR) :: plan_pf_r2c, plan_pf_c2r

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s h a p i r o 2 d x
  !   Purpose    :  conservative high-order 2-delta-x (Shapiro) filter of a 2D
  !                 field. The transfer function is  1 - s*[sin^2(k*dx/2)]^nord,
  !                 so at the grid scale (2dx) the amplitude is multiplied by
  !                 (1-s) for ANY order, while higher nord makes the filter
  !                 increasingly SELECTIVE: longer (resolved) waves are spared
  !                 ever more sharply (e.g. nord=1 still damps 4dx by ~0.75 per
  !                 pass, nord=8 leaves 4dx essentially untouched). 
  !                 Implemented as the order-nord negative Laplacian B applied
  !                 nord times:  X <- X - s*B^nord(X). B has eigenvalue
  !                 sin^2(k*dx/2) (=1 at 2dx, ~0 for smooth modes). Flux form
  !                 -> conserves the area-weighted integral sum(sqr*X) (energy
  !                 for tam, water mass for wcon). Zonal sweep (periodic) then
  !                 meridional sweep (no-flux poles).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine shapiro2dx(X, sqr, s, nord)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    real(wp), intent(in   ) :: sqr(:,:)
    real(wp), intent(in   ) :: s          ! filter strength in [0,1] (2dx damping = 1-s)
    integer,  intent(in   ) :: nord       ! filter order (selectivity); 1 = plain Laplacian

    integer :: i, j, n, imi, ipl, im, jm
    real(wp) :: P(size(X,1)), Q(size(X,1))         ! zonal work (row)
    real(wp) :: Pm(size(X,2)), Qm(size(X,2)), w(size(X,2))  ! merid work (column)
    real(wp) :: phim, phip


    im = size(X,1)
    jm = size(X,2)

    ! ---- zonal sweep (periodic; sqr is i-independent so the plain 2nd
    !      difference already conserves sum_i X) ----
    ! each latitude row j is independent (own private work buffers P,Q, writes the
    ! disjoint column X(:,j)) -> parallelise over j
    !$omp parallel do private(j,i,n,imi,ipl,P,Q)
    do j=1,jm
      do i=1,im
        P(i) = X(i,j)
      enddo
      ! apply the negative Laplacian B nord times: B(Z)_i = (2Z_i-Z_{i-1}-Z_{i+1})/4
      do n=1,nord
        do i=1,im
          imi = i-1
          if (imi.eq.0) imi = im
          ipl = i+1
          if (ipl.gt.im) ipl = 1
          Q(i) = 0.25_wp*(2._wp*P(i)-P(imi)-P(ipl))
        enddo
        do i=1,im
          P(i) = Q(i)
        enddo
      enddo
      ! X = X - s*B^nord(X)
      do i=1,im
        X(i,j) = X(i,j) - s*P(i)
      enddo
    enddo
    !$omp end parallel do

    ! ---- meridional sweep (conservative flux form, no flux through the poles) ----
    ! each longitude column i is independent (own private work buffers Pm,Qm,w, writes
    ! the disjoint row X(i,:)) -> parallelise over i
    !$omp parallel do private(i,j,n,Pm,Qm,w,phim,phip)
    do i=1,im
      do j=1,jm
        Pm(j) = X(i,j)
        w(j)  = sqr(i,j)
      enddo
      ! apply the negative Laplacian B nord times:
      !   B(Z)_j = -(phi_{j+1}-phi_j)/w_j , phi_j = 0.25*0.5*(w_{j-1}+w_j)*(Z_j-Z_{j-1})
      ! (eigenvalue sin^2 at 2dx, =1; conserves sum_j w_j*B(Z)_j = 0)
      do n=1,nord
        do j=1,jm
          if (j.gt.1) then
            phim = 0.25_wp*0.5_wp*(w(j-1)+w(j))*(Pm(j)-Pm(j-1))
          else
            phim = 0._wp
          endif
          if (j.lt.jm) then
            phip = 0.25_wp*0.5_wp*(w(j)+w(j+1))*(Pm(j+1)-Pm(j))
          else
            phip = 0._wp
          endif
          Qm(j) = -(phip-phim)/w(j)
        enddo
        do j=1,jm
          Pm(j) = Qm(j)
        enddo
      enddo
      ! X = X - s*B^nord(X)
      do j=1,jm
        X(i,j) = X(i,j) - s*Pm(j)
      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine shapiro2dx


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l i m i t _ f i l t _ w t r
  !   Purpose    :  keep the Shapiro filter of the moisture convergence from
  !                 draining a column. The filter is not positive definite: a
  !                 dry cell between strongly divergent neighbours (a 2dx
  !                 checkerboard, e.g. the Atacama) gets a filtered divergence
  !                 that removes more water in one step than the column holds,
  !                 and the qam floor then refills it, creating water that was
  !                 booked as negative precipitation.
  !                 Where the filtered convergence would take wcon below wmin
  !                 within tstep, its filter increment is reduced (never beyond
  !                 the unfiltered value). The withheld water is taken from the
  !                 positive filter increments within the filter's reach
  !                 (+-nord cells), so that sum(sqr*conv) is conserved.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine limit_filt_wtr(conv, conv_raw, sqr, wcon, wmin, evpa, tstep, nord)

    implicit none

    real(wp), intent(inout) :: conv(:,:)       ! filtered convergence, kg/m2/s
    real(wp), intent(in   ) :: conv_raw(:,:)   ! unfiltered convergence, kg/m2/s
    real(wp), intent(in   ) :: sqr(:,:)        ! cell area, m2
    real(wp), intent(in   ) :: wcon(:,:)       ! column water, kg/m2
    real(wp), intent(in   ) :: wmin(:,:)       ! lowest column water the step may leave, kg/m2
    real(wp), intent(in   ) :: evpa(:,:)       ! evaporation, kg/m2/s
    real(wp), intent(in   ) :: tstep           ! s
    integer,  intent(in   ) :: nord            ! filter order = reach of the filter in cells

    integer :: i, j, ii, jj, i2, im, jm
    real(wp) :: conv_min, conv_new, dw, pos, r
    real(wp) :: dinc(size(conv,1),size(conv,2))


    im = size(conv,1)
    jm = size(conv,2)

    ! filter increment
    dinc = conv - conv_raw

    ! serial: limited cells are rare, and the redistribution touches neighbours
    do j=1,jm
      do i=1,im
        if (dinc(i,j).ge.0._wp) cycle
        ! lowest convergence that leaves at least wmin after tstep
        conv_min = (wmin(i,j)-wcon(i,j))/tstep - evpa(i,j)
        if (conv(i,j).ge.conv_min) cycle
        conv_new = min(conv_raw(i,j), conv_min)
        if (conv_new.le.conv(i,j)) cycle
        dw = (conv_new-conv(i,j))*sqr(i,j)   ! kg/s withheld from this cell
        ! positive filter increments within reach
        pos = 0._wp
        do jj=max(1,j-nord),min(jm,j+nord)
          do ii=i-nord,i+nord
            i2 = modulo(ii-1,im)+1
            if (dinc(i2,jj).gt.0._wp) pos = pos + dinc(i2,jj)*sqr(i2,jj)
          enddo
        enddo
        if (pos.le.0._wp) cycle   ! nothing to take it from, leave this cell unlimited
        ! withhold at most what the neighbours received
        r = min(1._wp, dw/pos)
        do jj=max(1,j-nord),min(jm,j+nord)
          do ii=i-nord,i+nord
            i2 = modulo(ii-1,im)+1
            if (dinc(i2,jj).gt.0._wp) then
              conv(i2,jj) = conv(i2,jj) - r*dinc(i2,jj)
              dinc(i2,jj) = (1._wp-r)*dinc(i2,jj)
            endif
          enddo
        enddo
        conv(i,j) = conv(i,j) + r*pos/sqr(i,j)
        dinc(i,j) = conv(i,j) - conv_raw(i,j)
      enddo
    enddo

    return

  end subroutine limit_filt_wtr

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m o o t h 2
  !   Purpose    :  smoothing of 2D fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth2(X, niter)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: niter

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    ! smoothing is repeated niter times
    do iter=1,niter

      ! save field in temporary variable
      do j=1,jm
        do i=1,im
          Y(i,j) = X(i,j)
        enddo
      enddo        

      ! smooth with 4-cell neighbors
      do i=1,im
        imi = i-1
        if (imi.eq.0) imi = im
        ipl = i+1
        if (ipl.gt.im) ipl = 1

        j = 1
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
        do j=2,jm-1       
          X(i,j) = 0.2_wp*(Y(i,j)+Y(imi,j)+Y(ipl,j)+Y(i,j+1)+Y(i,j-1))
        enddo
        j = jm
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
      enddo

    enddo       

    deallocate(Y)

    return

  end subroutine smooth2


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m o o t h 2 _ m
  !   Purpose    :  smoothing of 2D fields
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth2_m(X,niter)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: niter

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    do iter=1,niter

      do j=1,jm
        do i=1,im
          Y(i,j) = X(i,j)
        enddo
      enddo        

      !$omp parallel do private(i,j,imi,ipl) 
      do i=1,im
        imi = i-1
        if (imi.eq.0) imi = im
        ipl = i+1
        if (ipl.gt.im) ipl = 1

        j = 1
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
        do j=2,jm-1       
          X(i,j) = 0.4_wp*Y(i,j) + 0.15_wp*(Y(imi,j)+Y(ipl,j)+Y(i,j+1)+Y(i,j-1))
        enddo
        j = jm
        X(i,j) = 1._wp/3._wp*(Y(i,j)+Y(imi,j)+Y(ipl,j))
      enddo
      !$omp end parallel do

    enddo       

    deallocate(Y)

    return

  end subroutine smooth2_m


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m o o t h 2 e q
  !   Purpose    :  smoothing of 2D fields around the equator
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth2eq(X, nj, niter)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: nj
    integer, intent(in) :: niter

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    do iter=1,niter

      do j=1,jm
        do i=1,im
          Y(i,j) = X(i,j)
        enddo
      enddo        

      do i=1,im
        imi = i-1
        if (imi.eq.0) imi = im
        ipl = i+1
        if (ipl.gt.im) ipl = 1

        do j=jm/2-nj,jm/2+1+nj       
          X(i,j) = 0.2_wp*(Y(i,j)+Y(imi,j)+Y(ipl,j)+Y(i,j+1)+Y(i,j-1))
        enddo
      enddo

    enddo       

    deallocate(Y)

    return

  end subroutine smooth2eq


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  z o n a
  !   Purpose    :  compute average of field over latitudinal belt
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function zona(X, cost, j1, j2)

    implicit none

    real(wp), intent(in) :: X(:)
    real(wp), intent(in) :: cost(:)
    integer, intent(in) :: j1, j2

    integer :: j
    real(wp) :: zona, weight


    zona = 0._wp
    weight = 0._wp

    do j=j1,j2
      zona = zona+cost(j)*X(j)  
      weight = weight+cost(j)
    enddo

    zona = zona/weight

    return

  end function zona


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  z o f i l
  !   Purpose    :  zonal filtering nite times for lat belt with index j 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine zofil(x,jf,nite)

    implicit none

    real(wp), intent(inout) :: X(:,:)
    integer, intent(in) :: jf, nite

    integer :: n, i, ipl, imi, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(X,1)
    jm = size(X,2)
    allocate(Y(im,jm))

    do n=1,nite

      do i=1,im
        Y(i,jf)=X(i,jf)
      enddo

      do i=1,im
        ipl=i+1
        if (ipl.gt.im) ipl=1
        imi=i-1
        if (imi.lt.1) imi=im 
        X(i,jf)=(Y(imi,jf)+Y(i,jf)+Y(ipl,jf))/3._wp
      enddo

    enddo

    deallocate(Y)

    return

  end subroutine zofil

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p o l a r _ f i l t e r
  !   Purpose    :  zonal Fourier truncation that keeps the resolution isotropic.
  !                 Row j keeps the zonal wavenumbers whose half wavelength,
  !                 pi*R*cos(lat)/m, is not shorter than the meridional spacing
  !                 dy = pi*R/jm, i.e. m <= jm*cos(lat), with at least m = 1.
  !                 Rows that resolve all wavenumbers (m_max >= im/2) are untouched.
  !                 Linear, so it commutes with zonal differences, and it keeps
  !                 the zonal mean of each row.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine polar_filter(x, cosl)

    implicit none

    real(wp), intent(inout) :: x(:,:)
    real(wp), intent(in) :: cosl(:)

    integer :: j

    call polar_filter_plan(size(x,1))

    !$omp parallel do private(j)
    do j=1,size(x,2)
      call polar_filter_row(x(:,j), size(x,2), cosl(j))
    enddo
    !$omp end parallel do

    return

  end subroutine polar_filter


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p o l a r _ f i l t e r 3
  !   Purpose    :  polar_filter on every level of a 3D field, (level,row) in parallel
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine polar_filter3(x, cosl)

    implicit none

    real(wp), intent(inout) :: x(:,:,:)
    real(wp), intent(in) :: cosl(:)

    integer :: j, k

    call polar_filter_plan(size(x,1))

    !$omp parallel do collapse(2) private(j,k)
    do k=1,size(x,3)
      do j=1,size(x,2)
        call polar_filter_row(x(:,j,k), size(x,2), cosl(j))
      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine polar_filter3


  ! FFTW plans for rows of length n, made once. The planner is not thread safe, so this is
  ! called outside the parallel loops; the rows are then transformed with the new-array execute
  ! interface (thread safe), FFTW_UNALIGNED because the row buffers are thread-private.
  subroutine polar_filter_plan(n)

    implicit none

    integer, intent(in) :: n

    real(C_DOUBLE) :: rowr(n)
    complex(C_DOUBLE_COMPLEX) :: rowc(n/2+1)

    !$omp critical (polar_filter_planner)
    if (n_pf_plan.ne.n) then
      if (n_pf_plan.gt.0) then
        call fftw_destroy_plan(plan_pf_r2c)
        call fftw_destroy_plan(plan_pf_c2r)
      endif
      plan_pf_r2c = fftw_plan_dft_r2c_1d(n, rowr, rowc, ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      plan_pf_c2r = fftw_plan_dft_c2r_1d(n, rowc, rowr, ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      n_pf_plan = n
    endif
    !$omp end critical (polar_filter_planner)

  end subroutine polar_filter_plan


  ! one row: keep the zonal wavenumbers 0..m_max, m_max = max(1, int(jm*cos(lat)))
  subroutine polar_filter_row(xr, jm, cosl)

    implicit none

    real(wp), intent(inout) :: xr(:)
    integer, intent(in) :: jm
    real(wp), intent(in) :: cosl

    integer :: n, mmax
    real(C_DOUBLE) :: rowr(size(xr))
    complex(C_DOUBLE_COMPLEX) :: rowc(size(xr)/2+1)

    n = size(xr)
    mmax = max(1, int(real(jm,wp)*cosl))
    if (2*mmax.ge.n) return

    rowr = xr
    call fftw_execute_dft_r2c(plan_pf_r2c, rowr, rowc)
    ! wavenumber m is rowc(m+1)
    rowc(mmax+2:) = (0._C_DOUBLE, 0._C_DOUBLE)
    call fftw_execute_dft_c2r(plan_pf_c2r, rowc, rowr)
    xr = rowr/real(n,wp)   ! FFTW is unnormalised

  end subroutine polar_filter_row

end module smooth_atm_mod

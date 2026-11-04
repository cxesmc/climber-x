!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : f r e e _ s u r f a c e _ m o d
!
!  Purpose : diagnose elevation of the free surface
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
module free_surface_mod

  use precision, only : wp
  use constants, only : g
  use ocn_params, only : rho0, fcor, fcorv, drag, drag_par
  use ocn_grid, only : maxi, maxj, maxk, k1, dz, dza, mask_ocn
  use ocn_grid, only : ku, rh, c, cv, rcv, ds, dsv, dphi, R_earth
  use island_mod, only : face_cb

  implicit none

  real(wp), parameter :: misval = -999._wp

  ! Cholesky factor of the face Laplacian, in envelope storage: row m of L spans columns env_first(m)..m and starts at
  ! lval(env_ptr(m)). The Laplacian only depends on the geometry, so it is factorized again only if the faces change.
  integer :: n_fac = -1
  integer :: nf_fac = -1
  integer :: nbasin = 0
  integer, allocatable :: f1_fac(:), f2_fac(:), basin(:), env_first(:), env_ptr(:)
  real(wp), allocatable :: w_fac(:), lval(:)

  private
  public :: free_surface

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f r e e _ s u r f a c e
  !   Purpose    :  diagnose elevation of the free surface
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine free_surface(rho, ub, psi, tau, f_ocn, ssh)

    implicit none

    real(wp), intent(in) :: rho(:,:,:)
    real(wp), dimension(:,0:,0:), intent(in) :: ub
    real(wp), dimension(0:,0:), intent(in) :: psi
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), intent(out) :: ssh(:,:)
    real(wp), allocatable :: ssh_tmp(:,:)
    integer :: i, j, nbr, ii, jj, iii, jjj
    real(wp) :: zs_sum


    call ssh_dynamic(rho, ub, psi, tau, f_ocn, ssh)

    allocate(ssh_tmp,source=ssh)

    ! extrapolate over land neighbors (needed for gradient computation)
    do j=1,maxj
      do i=1,maxi

        if (mask_ocn(i,j).eq.0) then
          nbr = 0
          zs_sum = 0._wp
          do ii=i-1,i+1
            do jj=j-1,j+1
              iii = ii
              if (iii.eq.0) iii = maxi
              if (iii.eq.maxi+1) iii = 1
              jjj = jj
              jjj = max(1,jjj)
              jjj = min(maxj,jjj)
              if (ssh_tmp(iii,jjj).ne.misval) then
                nbr = nbr+1
                zs_sum = zs_sum + ssh_tmp(iii,jjj)
              endif
            enddo
          enddo
          if (nbr.gt.0) then
            ssh(i,j) = zs_sum/nbr
          endif
        endif

      enddo
    enddo

    deallocate(ssh_tmp)

    return
  end subroutine free_surface


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s s h _ d y n a m i c
  !   Purpose    :  surface pressure of the rigid lid expressed as a height.
  !                 The depth-mean momentum balance on each velocity face, with the terms of the island integrand
  !                 (island.f90) and the velocity form of the Coriolis term, gives the difference of eta across the face:
  !                   g*d(eta) = -(drag*ub + cor + fcb - tau/(rho0*H))*L - <d(p_bc)>
  !                 p_bc = hydrostatic pressure of rho-rho0 from the surface, /rho0, <> = mean over the wet levels of the face.
  !                 eta is the area-weighted least squares fit of these differences, by direct solution of the normal
  !                 equations with a Cholesky factor kept between calls. Area mean over each connected basin = 0.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ssh_dynamic(rho, ub, psi, tau, f_ocn, ssh)

    implicit none

    real(wp), intent(in) :: rho(:,:,:)
    real(wp), dimension(:,0:,0:), intent(in) :: ub
    real(wp), dimension(0:,0:), intent(in) :: psi
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), dimension(:,:), intent(in) :: f_ocn
    real(wp), intent(out) :: ssh(:,:)

    logical, save :: first_call = .true.

    integer, allocatable :: id(:,:), f1(:), f2(:)
    real(wp), allocatable :: pbc(:,:,:), w(:), b(:)
    real(wp), allocatable :: x(:), area(:), eta_mean(:)
    integer :: i, j, k, ip1, n, nf, l
    real(wp) :: cor, fcb, len, force, dp, flx, misfit, target_norm


    ! wet cells
    allocate(id(maxi,maxj))
    id = 0
    n = 0
    do j=1,maxj
      do i=1,maxi
        if (mask_ocn(i,j).eq.1 .and. k1(i,j).le.maxk) then
          n = n+1
          id(i,j) = n
        endif
      enddo
    enddo

    ! hydrostatic pressure of rho-rho0 from the surface down at layer centres, /rho0 [m2/s2]
    allocate(pbc(maxi,maxj,maxk))
    pbc = 0._wp
    do j=1,maxj
      do i=1,maxi
        if (id(i,j).gt.0) then
          pbc(i,j,maxk) = g*(rho(i,j,maxk)-rho0)/rho0*0.5_wp*dz(maxk)
          do k=maxk-1,k1(i,j),-1
            pbc(i,j,k) = pbc(i,j,k+1) + g*(0.5_wp*(rho(i,j,k)+rho(i,j,k+1))-rho0)/rho0*dza(k)
          enddo
        endif
      enddo
    enddo

    ! difference of eta across the wet faces [m] and face weights (cross-face / along-face length)
    allocate(f1(2*maxi*maxj), f2(2*maxi*maxj), w(2*maxi*maxj), b(2*maxi*maxj))
    nf = 0
    do j=1,maxj
      do i=1,maxi
        if (id(i,j).eq.0) cycle
        ip1 = modulo(i,maxi) + 1

        ! u-face (1,i,j), cell (i,j) to cell (ip1,j)
        if (id(ip1,j).gt.0 .and. rh(1,i,j).gt.0._wp) then
          cor = -0.5_wp*(fcorv(j)*0.5_wp*(ub(2,i,j) + ub(2,i+1,j)) + fcorv(j-1)*0.5_wp*(ub(2,i,j-1) + ub(2,i+1,j-1)))
          fcb = 0._wp
          if (drag_par%drag_cb.gt.0._wp) fcb = face_cb(psi,1,i,j)
          len = R_earth*c(j)*dphi
          force = drag(1,i,j)*ub(1,i,j) + cor + fcb - tau(1,i,j)/rho0*rh(1,i,j)
          dp = 0._wp
          do k=ku(1,i,j),maxk
            dp = dp + dz(k)*(pbc(ip1,j,k)-pbc(i,j,k))
          enddo
          nf = nf+1
          f1(nf) = id(i,j)
          f2(nf) = id(ip1,j)
          b(nf) = (-force*len - dp*rh(1,i,j))/g
          w(nf) = R_earth*ds(j)/c(j)/len
        endif

        ! v-face (2,i,j), cell (i,j) to cell (i,j+1)
        if (j.lt.maxj) then
          if (id(i,j+1).gt.0 .and. rh(2,i,j).gt.0._wp) then
            cor = 0.5_wp*(fcor(j)*0.5_wp*(ub(1,i-1,j) + ub(1,i,j)) + fcor(j+1)*0.5_wp*(ub(1,i-1,j+1) + ub(1,i,j+1)))
            fcb = 0._wp
            if (drag_par%drag_cb.gt.0._wp) fcb = face_cb(psi,2,i,j)
            len = R_earth*rcv(j)*dsv(j)
            force = drag(2,i,j)*ub(2,i,j) + cor + fcb - tau(2,i,j)/rho0*rh(2,i,j)
            dp = 0._wp
            do k=ku(2,i,j),maxk
              dp = dp + dz(k)*(pbc(i,j+1,k)-pbc(i,j,k))
            enddo
            nf = nf+1
            f1(nf) = id(i,j)
            f2(nf) = id(i,j+1)
            b(nf) = (-force*len - dp*rh(2,i,j))/g
            w(nf) = R_earth*cv(j)*dphi/len
          endif
        endif
      enddo
    enddo

    ! normal equations: Laplacian(x) = rhs
    if (faces_changed(n, nf, f1, f2, w)) call factorize(n, nf, f1, f2, w)

    allocate(x(n))
    x = 0._wp
    do l=1,nf
      flx = w(l)*b(l)
      x(f1(l)) = x(f1(l)) - flx
      x(f2(l)) = x(f2(l)) + flx
    enddo
    call solve_factorized(n, x)

    ! area mean over each connected basin = 0, the balance only determines differences within a basin
    allocate(area(nbasin), eta_mean(nbasin))
    area = 0._wp
    eta_mean = 0._wp
    do j=1,maxj
      do i=1,maxi
        if (id(i,j).gt.0) then
          area(basin(id(i,j))) = area(basin(id(i,j))) + ds(j)*f_ocn(i,j)
          eta_mean(basin(id(i,j))) = eta_mean(basin(id(i,j))) + ds(j)*f_ocn(i,j)*x(id(i,j))
        endif
      enddo
    enddo
    where (area.gt.0._wp) eta_mean = eta_mean/area
    where (area.le.0._wp) eta_mean = 0._wp
    do l=1,n
      x(l) = x(l) - eta_mean(basin(l))
    enddo

    do j=1,maxj
      do i=1,maxi
        if (id(i,j).gt.0) then
          ssh(i,j) = x(id(i,j))
        else
          ssh(i,j) = misval
        endif
      enddo
    enddo

    if (first_call) then
      ! share of the face differences that is not a pure gradient
      misfit = 0._wp
      target_norm = 0._wp
      do l=1,nf
        misfit = misfit + w(l)*(x(f2(l))-x(f1(l))-b(l))**2
        target_norm = target_norm + w(l)*b(l)**2
      enddo
      print *,'ssh_dynamic: cells, faces, basins, gradient misfit ',n,nf,nbasin,sqrt(misfit/target_norm)
      first_call = .false.
    endif

    return
  end subroutine ssh_dynamic


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function   :  f a c e s _ c h a n g e d
  !   Purpose    :  check whether the Laplacian differs from the factorized one
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  logical function faces_changed(n, nf, f1, f2, w)

    implicit none

    integer, intent(in) :: n, nf
    integer, intent(in) :: f1(:), f2(:)
    real(wp), intent(in) :: w(:)

    integer :: l


    faces_changed = n.ne.n_fac .or. nf.ne.nf_fac
    if (faces_changed) return
    do l=1,nf
      if (f1(l).ne.f1_fac(l) .or. f2(l).ne.f2_fac(l) .or. w(l).ne.w_fac(l)) then
        faces_changed = .true.
        return
      endif
    enddo

    return
  end function faces_changed


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f a c t o r i z e
  !   Purpose    :  find the connected basins and compute the envelope Cholesky factor of the face Laplacian
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine factorize(n, nf, f1, f2, w)

    implicit none

    integer, intent(in) :: n, nf
    integer, intent(in) :: f1(:), f2(:)
    real(wp), intent(in) :: w(:)

    integer, allocatable :: parent(:), root_basin(:)
    logical, allocatable :: pinned(:)
    integer :: l, m, jc, kk, k0, r1, r2, lo, hi, row, col
    real(wp) :: tmp


    if (allocated(f1_fac)) deallocate(f1_fac, f2_fac, w_fac, basin, env_first, env_ptr, lval)
    allocate(f1_fac(nf), f2_fac(nf), w_fac(nf), basin(n), env_first(n), env_ptr(n+1))
    n_fac = n
    nf_fac = nf
    f1_fac = f1(1:nf)
    f2_fac = f2(1:nf)
    w_fac = w(1:nf)

    ! connected basins
    allocate(parent(n), root_basin(n))
    do m=1,n
      parent(m) = m
    enddo
    do l=1,nf
      r1 = f1(l)
      do while (parent(r1).ne.r1)
        r1 = parent(r1)
      enddo
      r2 = f2(l)
      do while (parent(r2).ne.r2)
        r2 = parent(r2)
      enddo
      if (r1.ne.r2) parent(r1) = r2
    enddo
    nbasin = 0
    root_basin = 0
    do m=1,n
      r1 = m
      do while (parent(r1).ne.r1)
        r1 = parent(r1)
      enddo
      if (root_basin(r1).eq.0) then
        nbasin = nbasin+1
        root_basin(r1) = nbasin
      endif
      basin(m) = root_basin(r1)
    enddo

    ! envelope of the lower triangle
    do m=1,n
      env_first(m) = m
    enddo
    do l=1,nf
      hi = max(f1(l),f2(l))
      lo = min(f1(l),f2(l))
      env_first(hi) = min(env_first(hi),lo)
    enddo
    env_ptr(1) = 1
    do m=1,n
      env_ptr(m+1) = env_ptr(m) + m-env_first(m)+1
    enddo

    ! assemble
    allocate(lval(env_ptr(n+1)-1))
    lval = 0._wp
    do l=1,nf
      hi = max(f1(l),f2(l))
      lo = min(f1(l),f2(l))
      lval(env_ptr(hi)+lo-env_first(hi)) = lval(env_ptr(hi)+lo-env_first(hi)) - w(l)
      lval(env_ptr(f1(l)+1)-1) = lval(env_ptr(f1(l)+1)-1) + w(l)
      lval(env_ptr(f2(l)+1)-1) = lval(env_ptr(f2(l)+1)-1) + w(l)
    enddo

    ! The Laplacian is singular, with a constant null vector on each basin. Extra weight on the diagonal of one cell per
    ! basin makes it positive definite: for a right hand side summing to zero over the basin the cell gets x = 0 and the
    ! other equations hold exactly. An isolated cell has no equation left, it gets x = 0 as well.
    allocate(pinned(nbasin))
    pinned = .false.
    do m=1,n
      if (.not.pinned(basin(m))) then
        pinned(basin(m)) = .true.
        if (lval(env_ptr(m+1)-1).gt.0._wp) then
          lval(env_ptr(m+1)-1) = 2._wp*lval(env_ptr(m+1)-1)
        else
          lval(env_ptr(m+1)-1) = 1._wp
        endif
      endif
    enddo

    ! Cholesky factorization, row by row
    do m=1,n
      row = env_ptr(m)-env_first(m)
      do jc=env_first(m),m
        col = env_ptr(jc)-env_first(jc)
        k0 = max(env_first(m),env_first(jc))
        tmp = lval(row+jc)
        do kk=k0,jc-1
          tmp = tmp - lval(row+kk)*lval(col+kk)
        enddo
        if (jc.lt.m) then
          lval(row+jc) = tmp/lval(col+jc)
        else
          lval(row+m) = sqrt(tmp)
        endif
      enddo
    enddo

    return
  end subroutine factorize


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s o l v e _ f a c t o r i z e d
  !   Purpose    :  forward and back substitution with the Cholesky factor, x is the right hand side on input
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine solve_factorized(n, x)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(inout) :: x(:)

    integer :: m, kk, row
    real(wp) :: tmp


    do m=1,n
      row = env_ptr(m)-env_first(m)
      tmp = x(m)
      do kk=env_first(m),m-1
        tmp = tmp - lval(row+kk)*x(kk)
      enddo
      x(m) = tmp/lval(row+m)
    enddo
    do m=n,1,-1
      row = env_ptr(m)-env_first(m)
      x(m) = x(m)/lval(row+m)
      tmp = x(m)
      do kk=env_first(m),m-1
        x(kk) = x(kk) - tmp*lval(row+kk)
      enddo
    enddo

    return
  end subroutine solve_factorized

end module free_surface_mod

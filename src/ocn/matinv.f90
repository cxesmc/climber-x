!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : m a t i n v _ m o d
!
!  Purpose : matrix inversion
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
module matinv_mod
  ! includes subroutines to solve a set of n linear equations by direct inversion
  ! solves amat*x = rhs, putting result into rhs
  ! split in to two parts to invert and multiply separately
  ! to save a lot of cpu if matrix is constant in time

  use precision, only : wp

  implicit none

  private
  public :: matinv, matmult

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m a t i n v
  !   Purpose    :  invert
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine matinv(nisl, &
                    amat,ipiv)

    implicit none

    integer, intent(in) :: nisl
    real(wp), intent(inout) :: amat(:,:)
    integer, intent(out) :: ipiv(:)

    integer :: i, j, k, p
    real(wp) :: tmp
    real(wp), parameter :: pivtol = 1.e-30_wp   ! singularity floor (exact-zero / denormal guard)

    ! LU factorization with partial (row) pivoting: P*amat = L*U.
    ! U is stored in the upper triangle (incl. diagonal), the unit-lower multipliers L in
    ! the strict lower triangle, and the row swaps in ipiv; matmult replays ipiv on the rhs
    ! and back-solves. Replaces the old fraction-free elimination, whose entries grew like
    ! the running product of pivots (overflow risk) and which could divide by a zero pivot
    ! without warning. Only the factorization columns (1..nisl) are touched; the rhs column
    ! held in amat(:,nisl+1) is left for matmult.

    do i=1,nisl
       ! find pivot row (largest magnitude in column i, rows i..nisl)
       p = i
       do j=i+1,nisl
          if (abs(amat(j,i)).gt.abs(amat(p,i))) p = j
       enddo
       ipiv(i) = p
       ! swap rows i and p over the factorization columns
       if (p.ne.i) then
          do k=1,nisl
             tmp = amat(i,k); amat(i,k) = amat(p,k); amat(p,k) = tmp
          enddo
       endif
       if (abs(amat(i,i)).lt.pivtol) then
          print *
          print *,'ERROR: singular island matrix in matinv'
          print *,'island',i,' pivot',amat(i,i)
          print *
          stop 'matinv: singular island matrix'
       endif
       ! eliminate below the pivot, storing the multipliers in the lower triangle
       do j=i+1,nisl
          amat(j,i) = amat(j,i)/amat(i,i)
          do k=i+1,nisl
             amat(j,k) = amat(j,k) - amat(j,i)*amat(i,k)
          enddo
       enddo
    enddo

   return

  end subroutine matinv


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m a t m u l t
  !   Purpose    :  multiply
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine matmult(nisl,amat,ipiv, &
                     rhs)

    implicit none

    integer, intent(in) :: nisl
    real(wp), intent(in) :: amat(:,:)
    integer, intent(in) :: ipiv(:)
    real(wp), intent(inout) :: rhs(:)

    integer :: i, j
    real(wp) :: tmp

    ! solve P*amat = L*U system for rhs using the factors stored by matinv

    ! apply the row permutation recorded during factorization
    do i=1,nisl
       if (ipiv(i).ne.i) then
          tmp = rhs(i); rhs(i) = rhs(ipiv(i)); rhs(ipiv(i)) = tmp
       endif
    enddo

    ! forward substitution (unit lower triangular L)
    do i=2,nisl
       do j=1,i-1
          rhs(i) = rhs(i) - amat(i,j)*rhs(j)
       enddo
    enddo

    ! back substitution (upper triangular U)
    do i=nisl,1,-1
       do j=i+1,nisl
          rhs(i) = rhs(i) - amat(i,j)*rhs(j)
       enddo
       rhs(i) = rhs(i)/amat(i,i)
    enddo

   return

  end subroutine matmult

end module matinv_mod

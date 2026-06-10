!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t r i d i a g
!
!  Purpose : solve tridiagonal system of equations
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
module tridiag

  use precision, only : wp

  implicit none

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t r i d i a g _ s o l v e
  !   Purpose    :  solve tridiagonal system of equations 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine tridiag_solve(a,b,c,r,x,n)

    implicit none

!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        r - right part
!        x - the answer
!        n - number of equations

        real(wp), intent(in), dimension(:) :: a, b, c, r
        real(wp), intent(out), dimension(:) :: x
        integer, intent(in) :: n

        integer :: i
        real(wp), dimension(n) :: cp, dp
        real(wp) :: m


! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = r(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (r(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

        return

   end subroutine tridiag_solve


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c y c l i c _ t r i d i a g _ s o l v e
  !   Purpose    :  solve a periodic (cyclic) tridiagonal system via the
  !                 Sherman-Morrison formula, reusing the Thomas solver.
  !                 a(1) is the (1,n) corner, c(n) is the (n,1) corner.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine cyclic_tridiag_solve(a, b, c, r, x, n)

    implicit none

    real(wp), intent(in   ) :: a(:), b(:), c(:), r(:)
    real(wp), intent(out  ) :: x(:)
    integer,  intent(in   ) :: n

    integer :: i
    real(wp) :: gamma, alpha, beta, fact
    real(wp) :: bb(n), u(n), z(n)


    ! Desired periodic corners: M(1,n) = a(1) (row-1 sub-diagonal wraps to col n),
    ! M(n,1) = c(n) (row-n super-diagonal wraps to col 1). The Sherman-Morrison
    ! rank-1 update u*v^T with u=(gamma,..,alpha), v=(1,..,beta/gamma) produces
    ! M(1,n)=beta and M(n,1)=alpha, so we must set alpha=c(n), beta=a(1) (NOT the
    ! naive alpha=a(1), beta=c(n), which swaps the corners and corrupts the i=1/i=n
    ! seam whenever a/=c, e.g. with upwind advection).
    alpha = c(n)    ! -> goes to M(n,1)
    beta  = a(1)    ! -> goes to M(1,n)
    gamma = -b(1)   ! arbitrary nonzero

    bb(1) = b(1) - gamma
    bb(n) = b(n) - alpha*beta/gamma
    do i=2,n-1
      bb(i) = b(i)
    enddo

    call tridiag_solve(a, bb, c, r, x, n)

    u(:) = 0._wp
    u(1) = gamma
    u(n) = alpha
    call tridiag_solve(a, bb, c, u, z, n)

    fact = (x(1) + beta*x(n)/gamma) / (1._wp + z(1) + beta*z(n)/gamma)
    do i=1,n
      x(i) = x(i) - fact*z(i)
    enddo

    return

  end subroutine cyclic_tridiag_solve


end module tridiag



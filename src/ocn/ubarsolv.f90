!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : u b a r s o l v _ m o d
!
!  Purpose : calculate barotropic velocity
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
module ubarsolv_mod

  use precision, only : wp
  use ocn_params, only : rho0
  use ocn_grid, only : maxi, maxj, c, rcv, rdphi, rds, rh, R_earth

  implicit none

  private
  public :: ubarsolv

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u b a r s o l v
  !   Purpose    :  calculate barotropic velocity on c grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ubarsolv(Lband,Uband,Udiag, &
                      gb, &
                      ub,psi)

    implicit none

    ! contiguous band storage of the LU factors (built once in momentum from ratm/gap):
    !   Lband(off,i) = ratm(i+off,off) ; Uband(off,i) = gap(i,n+2+off) ; Udiag(i) = gap(i,n+2)
    real(wp), dimension(:,:), intent(in) :: Lband
    real(wp), dimension(:,:), intent(in) :: Uband
    real(wp), dimension(:),   intent(in) :: Udiag

    real(wp), dimension(:), intent(inout) :: gb

    real(wp), dimension(:,0:,0:), intent(out) :: ub
    real(wp), dimension(0:,0:), intent(out) :: psi

    integer i, j, k, n, m, km, im, off


    n = maxi
    m = maxj + 1

    ! solve Psi equation

    ! forward substitution: gb(i+1:im) -= gb(i)*Lband(1:im-i,i) 
    do i=1,n*m-1
       im = min(i+n+1,n*m)
       do off=1,im-i
          gb(i+off) = gb(i+off) - Lband(off,i)*gb(i)
       enddo
    enddo
    ! back substitution: unit-stride dot product over the upper band
    gb(n*m) = gb(n*m)/Udiag(n*m)
    do i=n*m-1,1,-1
       km=min(n+1,n*m-i)
       do k=1,km
          gb(i)=gb(i) - Uband(k,i)*gb(i+k)
       enddo
       gb(i)=gb(i)/Udiag(i) ! kg/s
    enddo

    ! write to Psi for convenience
    do j=0,maxj
       do i=1,maxi
          k = i + j*n
          psi(i,j) = gb(k)  ! kg/s
       enddo
       psi(0,j) = psi(maxi,j)
    enddo

    ! derive barotropic velocity from streamfunction Psi
    do j=1,maxj
       do i=1,maxi
          ub(1,i,j) = -rh(1,i,j)*c(j)*(psi(i,j) - psi(i,j-1))*rds(j)/(rho0*R_earth) ! m/s
       enddo
    enddo
    do j=1,maxj-1
       do i=1,maxi
          ub(2,i,j) = rh(2,i,j)*(psi(i,j) - psi(i-1,j))*rcv(j)*rdphi/(rho0*R_earth) ! m/s
       enddo
    enddo

    ! set velocity to zero at N and S boundaries
    do i=1,maxi
       ub(2,i,maxj) = 0._wp
       ub(2,i,0) = 0._wp
    enddo

    ! periodic b.c. for ub(2) required only for island integral
    do j=1,maxj
       ub(2,maxi+1,j) = ub(2,1,j)
       ub(1,0,j) = ub(1,maxi,j)
       ub(1,maxi+1,j) = ub(1,1,j)
       ub(2,0,j) = ub(2,maxi,j)
    enddo
    ub(2,maxi+1,0) = ub(2,1,0)
    ub(2,0,0) = ub(2,maxi,0)


   return
         
  end subroutine ubarsolv

end module ubarsolv_mod

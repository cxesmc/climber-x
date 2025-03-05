!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b m b _ g r i d 
!
!  Purpose : grid definitions fro BMB model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module bmb_grid

  use precision, only : wp

  implicit none

  ! levels
  integer,  parameter :: nl = 4
  real(wp), dimension(0:nl) :: z = (/0._wp,0.2_wp,1._wp,5._wp,15._wp/)! z(0) is overwritten with 0.5*h_snow

  real(wp), dimension(0:nl) :: z_int, dz, rdz, rdz_pos, rdz_neg

contains

  subroutine bmb_grid_init

  implicit none

  integer :: k


    ! vertical grid 
    z(0) = 0._wp 
    dz(0) = 0._wp
    ! vertical layers thickness
    dz(1) = 0.5_wp * ( z(1) + z(2) )
    do k=2,nl-1
     dz(k) = 0.5_wp * ( z(k+1) - z(k-1) )
    enddo
    dz(nl) = z(nl) - z(nl-1)

    ! depth of vertical layer interfaces
    z_int(0) = 0._wp ! snow - soil interface
    do k=1,nl-1
     z_int(k) = 0.5_wp * ( z(k) + z(k+1) )
    enddo
    z_int(nl) = z(nl) + 0.5_wp * dz(nl)

    ! reciprocals to speed up fortran
    rdz(1:nl) = 1._wp/dz(1:nl)
    do k=1,nl-1
     rdz_pos(k) = 1._wp/(z(k+1)-z(k))
    enddo
    rdz_pos(nl) = 0._wp
    do k=2,nl
     rdz_neg(k) = 1._wp/(z(k)-z(k-1))
    enddo

  return

  end subroutine bmb_grid_init

end module bmb_grid




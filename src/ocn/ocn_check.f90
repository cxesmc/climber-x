!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ c h e c k
!
!  Purpose : check velocity range
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
module ocn_check

  use precision, only : wp
  use timer, only : year, doy
  use ocn_grid, only : k1, maxi, maxj, maxk, dx, dy
  use ocn_params, only : dt

  implicit none

  private
  public :: check_vel

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c h e c k _ v e l
  !   Purpose    :  check velocity range 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine check_vel(u, error)

    implicit none

    real(wp), dimension(1:,0:,0:,:), intent(in) :: u
    logical, intent(inout) :: error

    integer :: i, j, k
    real(wp), dimension(maxj) :: u_min_CFL
    real(wp), dimension(maxj) :: u_max_CFL
    real(wp) :: v_min_CFL
    real(wp) :: v_max_CFL
    real(wp), parameter :: vlim = 1._wp ! m/s


    do j=1,maxj
      !u_min_CFL(j) = max(-vlim,-dx(j)/dt)
      !u_max_CFL(j) = min( vlim, dx(j)/dt)
      u_min_CFL(j) = -dx(j)/dt
      u_max_CFL(j) =  dx(j)/dt
    enddo
    !v_min_CFL = max(-vlim,-dy/dt)
    !v_max_CFL = min( vlim, dy/dt)
    v_min_CFL = -dy/dt
    v_max_CFL =  dy/dt

    ! check CFL stability
    do j=1,maxj
      do i=1,maxi
        do k=k1(i,j),maxk
          if (u(1,i,j,k).gt.u_max_CFL(j)) then
            print *
            print *,'WARNING: violation of CFL criterium!'
            print *,'(i,j,k)',i,j,k
            print *,'doy',doy
            print *,'u',u(1,i,j,k)
            print *,'CFL: ',abs(u(1,i,j,k))*dt/dx(j)
            print *,'Try reducing the ocean time step dt_day_ocn in control.nml'
            if (year>1 .and. doy>1) error = .true.
          endif
          if (u(1,i,j,k).lt.u_min_CFL(j)) then
            print *
            print *,'WARNING: violation of CFL criterium!'
            print *,'(i,j,k)',i,j,k
            print *,'doy',doy
            print *,'u',u(1,i,j,k)
            print *,'CFL: ',abs(u(1,i,j,k))*dt/dx(j)
            print *,'Try reducing the ocean time step dt_day_ocn in control.nml'
            if (year>1 .and. doy>1) error = .true.
          endif
          if (u(2,i,j,k).gt.v_max_CFL) then
            print *
            print *,'WARNING: violation of CFL criterium!'
            print *,'(i,j,k)',i,j,k
            print *,'doy',doy
            print *,'v',u(2,i,j,k)
            print *,'CFL: ',abs(u(2,i,j,k))*dt/dy
            print *,'Try reducing the ocean time step dt_day_ocn in control.nml'
            if (year>1 .and. doy>1) error = .true.
          endif
          if (u(2,i,j,k).lt.v_min_CFL) then
            print *
            print *,'WARNING: violation of CFL criterium!'
            print *,'(i,j,k)',i,j,k
            print *,'doy',doy
            print *,'v',u(2,i,j,k)
            print *,'CFL: ',abs(u(2,i,j,k))*dt/dy
            print *,'Try reducing the ocean time step dt_day_ocn in control.nml'
            if (year>1 .and. doy>1) error = .true.
          endif
        enddo
      enddo
    enddo


    return

  end subroutine check_vel


end module ocn_check

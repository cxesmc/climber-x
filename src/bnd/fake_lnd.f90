!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f a k e _ l n d _ m o d
!
!  Purpose : fake land
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
module fake_lnd_mod

  use precision, only : wp
  use timer, only : doy, nday_year, monthly2daily
  use climber_grid, only : ni, nj
  use control, only : fake_lnd_const_file
  use ncio

  implicit none

  type fake_lnd_type
    real(wp), dimension(:,:), allocatable :: runoff, discharge
  end type fake_lnd_type

  real(wp), dimension(:,:,:), allocatable :: runoff, discharge

  integer, dimension(nday_year) :: m0, m1
  real(wp), dimension(nday_year) :: wtm0, wtm1
  
  private
  public :: fake_lnd_init, fake_lnd_update, fake_lnd_type


contains

  subroutine fake_lnd_init(lnd)

  implicit none

  type(fake_lnd_type) :: lnd


  ! allocate land type variables
  allocate(lnd%runoff(ni,nj))
  allocate(lnd%discharge(ni,nj))


  ! allocate local variables
  allocate(runoff(ni,nj,0:13))
  allocate(discharge(ni,nj,0:13))

  ! get weights for interpolation from monthly to daily
  call monthly2daily(m0,m1,wtm0,wtm1)

  call nc_read(fake_lnd_const_file,"runoff",runoff(:,:,1:12),start=[1,1,1],count=[ni,nj,12])
  call nc_read(fake_lnd_const_file,"discharge",discharge(:,:,1:12),start=[1,1,1],count=[ni,nj,12])

  ! periodic boundary conditions
  runoff(:,:,0) = runoff(:,:,12)
  discharge(:,:,0) = discharge(:,:,12)

  runoff(:,:,13) = runoff(:,:,1)
  discharge(:,:,13) = discharge(:,:,1)

  ! initial values
  lnd%discharge = discharge(:,:,1)
  lnd%runoff    = runoff(:,:,1)


  return

  end subroutine fake_lnd_init


  subroutine fake_lnd_update(lnd)

  implicit none

  type(fake_lnd_type), intent(inout) :: lnd


  lnd%discharge = wtm0(doy)*discharge(:,:,m0(doy)) + wtm1(doy)*discharge(:,:,m1(doy))
  lnd%runoff = wtm0(doy)*runoff(:,:,m0(doy)) + wtm1(doy)*runoff(:,:,m1(doy)) 


  return

  end subroutine fake_lnd_update


end module fake_lnd_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : n e p t u n e _ m o d
!
!  Purpose : Neptune parameterisation of eddy-topography interaction
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Neil R. Edwards and Matteo Willeit
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
! Eddy-topography interaction (Holloway 1992; Eby and Holloway 1994) drives the
! flow towards a non-zero statistical equilibrium rather than towards rest. The
! equilibrium depth-integrated volume transport streamfunction is
!
!    psi* = - f L^2 H
!
! with f the Coriolis parameter, H the water depth and L a length scale of the
! order of the eddy scale. The associated depth-averaged velocity u* follows from
! psi* exactly as the model's barotropic velocity follows from psi. This is the same
! formulation as in the original MOM implementation of Eby and Holloway, where
! pnep = -f*snep*snep*hnep and the velocity is the streamfunction gradient times 1/H.
!
! Eby and Holloway (1994) let L vary with latitude as a cosine fit between an equatorial
! and a polar value (12 km and 3 km respectively), mimicking the latitudinal variation of
! the deformation radius; see i_neptune_l==2. A uniform L is available as i_neptune_l==1.
!
! The parameterisation replaces the linear drag  lambda*u  by  lambda*(u - u*),
! i.e. it adds a body force  lambda*u*  to the depth-integrated momentum equation.
! That force enters the barotropic problem in exactly the same way as the surface
! wind stress does, so it is implemented here as an equivalent stress
!
!    tau* = rho0 * lambda * H * u*     [N/m2]
!
! which is simply added to the wind stress before the streamfunction is solved.
! This guarantees that the forcing and the island path integrals stay consistent,
! since both are computed from the same stress array.
!
! Because psi* is proportional to H it vanishes at the coast, which is the correct
! boundary condition for a transport streamfunction, and it varies most rapidly
! across the continental slope. The resulting flow is cyclonic around basins:
! equatorward along western continental slopes in the northern hemisphere, i.e. in
! the sense of the Labrador Current and the Oyashio.
!
! Note that this is a purely barotropic effect: the baroclinic shear computed in
! velc is deliberately left untouched.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module neptune_mod

  use precision, only : wp
  use constants, only : omega
  use ocn_grid, only : maxi, maxj, h, rh, c, cv, sv, rds, rcv, rdphi, R_earth
  use ocn_params, only : drag, rho0, neptune_par

  implicit none

  real(wp), allocatable :: psi_neptune(:,:)   !! Neptune volume transport streamfunction at psi points [m3/s]
  real(wp), allocatable :: u_neptune(:,:,:)   !! Neptune depth-averaged velocity on u/v grid [m/s]
  real(wp), allocatable :: tau_neptune(:,:,:) !! equivalent stress on u/v grid [N/m2]
  real(wp), allocatable :: lscale(:)          !! Neptune length scale [m]

  private
  public :: neptune_init, neptune_update
  public :: psi_neptune, u_neptune, tau_neptune

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  n e p t u n e _ i n i t
  !   Purpose    :  allocate variables and set the length scale
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine neptune_init

    implicit none

    integer :: j


    allocate(psi_neptune(0:maxi,0:maxj))
    allocate(u_neptune(2,maxi,maxj))
    allocate(tau_neptune(2,maxi,maxj))
    allocate(lscale(0:maxj))

    psi_neptune(:,:)   = 0._wp
    u_neptune(:,:,:)   = 0._wp
    tau_neptune(:,:,:) = 0._wp

    ! length scale at psi points (cell edges in latitude)
    do j=0,maxj
      if (neptune_par%i_l.eq.1) then
        ! uniform length scale
        lscale(j) = neptune_par%l
      else
        ! Eby and Holloway (1994): cosine fit between an equatorial and a polar value,
        ! mimicking the latitudinal variation of the deformation radius. Written there as
        ! L = L_pol + (L_eq-L_pol)*(0.5+0.5*cos(2*lat)), which is the same as the cos^2 form below.
        lscale(j) = neptune_par%l_pol + (neptune_par%l_eq - neptune_par%l_pol)*cv(j)*cv(j)
      endif
    enddo

    return

  end subroutine neptune_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  n e p t u n e _ u p d a t e
  !   Purpose    :  compute the Neptune equilibrium flow and the equivalent
  !                 stress needed to drive the model towards it.
  !                 Depends on geometry and drag only, so it has to be called
  !                 only when the ocean grid is updated.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine neptune_update

    implicit none

    integer :: i, j
    real(wp) :: fcorv_j, h_psi


    ! Neptune equilibrium streamfunction at psi points (grid cell corners).
    ! h is zero over land, so h_psi tapers to zero towards the coast and psi* with it.
    do j=0,maxj
      ! Coriolis parameter at psi points. Deliberately the true value, not the
      ! fcormin-limited one used for the geostrophic balance: fcormin exists to keep the
      ! momentum solver well behaved near the equator, whereas here f only sets the
      ! amplitude of the equilibrium flow and must be allowed to vanish at the equator.
      fcorv_j = 2._wp*omega*sv(j)
      do i=0,maxi
        h_psi = 0.25_wp*(h(3,i,j) + h(3,i+1,j) + h(3,i,j+1) + h(3,i+1,j+1))
        psi_neptune(i,j) = -fcorv_j*lscale(j)*lscale(j)*h_psi  ! m3/s
      enddo
    enddo

    ! depth-averaged Neptune velocity, derived from psi* exactly as ubarsolv derives
    ! the barotropic velocity from psi (psi* is a volume transport, hence no rho0 here).
    ! rh is zero wherever the velocity point is not wet, so u* vanishes there.
    u_neptune(:,:,:) = 0._wp
    do j=1,maxj
      do i=1,maxi
        u_neptune(1,i,j) = -rh(1,i,j)*c(j)*(psi_neptune(i,j) - psi_neptune(i,j-1))*rds(j)/R_earth  ! m/s
      enddo
    enddo
    do j=1,maxj-1
      do i=1,maxi
        u_neptune(2,i,j) = rh(2,i,j)*(psi_neptune(i,j) - psi_neptune(i-1,j))*rcv(j)*rdphi/R_earth  ! m/s
      enddo
    enddo

    ! limit the Neptune velocity. Needed because u* scales with 1/H and can become
    ! unreasonably large over very shallow topography.
    do j=1,maxj
      do i=1,maxi
        u_neptune(1,i,j) = sign(min(abs(u_neptune(1,i,j)),neptune_par%u_max),u_neptune(1,i,j))
        u_neptune(2,i,j) = sign(min(abs(u_neptune(2,i,j)),neptune_par%u_max),u_neptune(2,i,j))
      enddo
    enddo

    ! equivalent stress. Constructed so that tau*/rho0*rh = drag*u*, i.e. so that adding
    ! tau* to the wind stress turns the drag term  drag*u  into  drag*(u-u*), both in the
    ! streamfunction forcing (wind) and in the island path integrals (island).
    ! h is zero at dry velocity points, so tau* vanishes there.
    do j=1,maxj
      do i=1,maxi
        tau_neptune(1,i,j) = rho0*drag(1,i,j)*h(1,i,j)*u_neptune(1,i,j)  ! N/m2
        tau_neptune(2,i,j) = rho0*drag(2,i,j)*h(2,i,j)*u_neptune(2,i,j)
      enddo
    enddo

    return

  end subroutine neptune_update

end module neptune_mod

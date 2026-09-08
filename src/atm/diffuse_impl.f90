!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : d i f f u s e _ i m p l
!
!  Purpose : implicit horizontal diffusion (zonal + meridional, ADI) of the
!            prognostic atmospheric column variables (dry static energy /
!            temperature, column water, dust and CO2).
!
!            EXACT 2D EQUIVALENT. Because the diffusivity is vertically uniform,
!            the column-integrated zonal diffusive convergence is reproduced
!            exactly from 2D quantities by solving, per field,
!
!              (I - dt*fac*L) * dprog = dt*fac*conv_expl ,   prog := prog + dprog
!
!            where conv_expl = (fdx(i)-fdx(i+1))/sqr is adifa's EXACT zonal
!            diffusive flux convergence (the original layer-resolved physics; the
!            gad*z term cancels because adifa differences tp at common height
!            levels and zeroes sub-surface layers via the face mass dplx), and L
!            is a column Laplacian present ONLY on the LHS, so it does not change
!            WHAT is diffused - it only damps the stiff polar grid-scale modes.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module diffuse_impl_mod

  use atm_params, only : wp
  use atm_params, only : tstep, hatm, ra, cp, l_dust
  use atm_grid, only : im, jm, jmc
  use atm_grid, only : plx_trop, ply_trop, sqr, dxt, dxu, dy, cheat
  use tridiag, only : tridiag_solve, cyclic_tridiag_solve
  !$ use omp_lib

  implicit none

  private
  public :: diffuse_impl

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  d i f f u s e _ i m p l
  !   Purpose    :  implicit diffusion of the prognostic columns
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diffuse_impl(diffx, diffy, &  ! in
      ra2a, hdust, fdxdse, fdxwtr, fdxdst, fdxco2, fdydse, fdywtr, fdydst, fdyco2, &      ! in
      tam, wcon, dam, cam, &                                                              ! inout
      convwtr_dif)                                                                        ! out

    implicit none

    real(wp), intent(in   ) :: diffx(:,:), diffy(:,:)
    real(wp), intent(in   ) :: ra2a(:,:)
    real(wp), intent(in   ) :: hdust(:,:)
    real(wp), intent(in   ) :: fdxdse(:,:), fdydse(:,:)   ! adifa exact DSE diffusive fluxes (zonal/merid)
    real(wp), intent(in   ) :: fdxwtr(:,:), fdywtr(:,:)   ! adifa exact water diffusive fluxes
    real(wp), intent(in   ) :: fdxdst(:,:), fdydst(:,:)   ! adifa exact dust diffusive fluxes
    real(wp), intent(in   ) :: fdxco2(:,:), fdyco2(:,:)   ! adifa exact CO2 diffusive fluxes

    real(wp), intent(inout) :: tam(:,:)
    real(wp), intent(inout) :: wcon(:,:)
    real(wp), intent(inout) :: dam(:,:)
    real(wp), intent(inout) :: cam(:,:)

    real(wp), intent(out  ) :: convwtr_dif(:,:)

    integer :: i, j
    real(wp) :: fac(im,jm)
    real(wp) :: heff
    real(wp) :: wcon_old(im,jm)


    ! implicit diffusion (zonal sweep then meridional sweep) 
    ! (I - dt*fac*L)*dprog = dt*fac*conv_expl, prog += dprog,
    ! where conv_expl = (fdx(i)-fdx(i+1))/sqr is adifa's exact diffusive flux
    ! convergence. L (the column Laplacian) is
    ! on the LHS only and provides the unconditionally-stable implicit damping of the
    ! stiff polar grid-scale modes; it does not change what is diffused.

    !---------------------------------------------------------------
    ! dry static energy: dtam = dt * cp/cheat * conv.  The transport carries DSE (cp*T per
    ! unit mass) and the reservoir is the column enthalpy cheat=pzsa*amas*cp, so the ratio
    ! is simply 1/M, M being the local column mass.  Same conversion time_step applies to
    ! the advective convdse (deba/cheat), so the implicit and explicit diffusion paths and
    ! the column energy budget all use one and the same heat capacity.
    !---------------------------------------------------------------
    fac(:,:) = cp/cheat(:,:)
    call zonal_diffuse(tam, diffx, fac, .true., fdxdse)
    call merid_diffuse(tam, diffy, fac, .true., fdydse)

    !---------------------------------------------------------------
    ! column water content 
    !---------------------------------------------------------------
    wcon_old(:,:) = wcon(1:im,1:jm)
    fac(:,:) = 1._wp
    call zonal_diffuse(wcon, diffx, fac, .false., fdxwtr)
    call merid_diffuse(wcon, diffy, fac, .false., fdywtr)
    ! realized implicit diffusive moisture convergence (zonal+meridional, for precip)
    do j=1,jm
      do i=1,im
        convwtr_dif(i,j) = (wcon(i,j) - wcon_old(i,j))/tstep   ! kg/m2/s
      enddo
    enddo

    !---------------------------------------------------------------
    ! dust: ddam = dt * 1/(heff*ra) * conv_expl
    !---------------------------------------------------------------
    if (l_dust) then
      do j=1,jm
        do i=1,im
          heff = hdust(i,j)*hatm/(hdust(i,j)+hatm)   ! m
          fac(i,j) = 1._wp/(heff*ra)
        enddo
      enddo
      call zonal_diffuse(dam, diffx, fac, .true., fdxdst)
      call merid_diffuse(dam, diffy, fac, .true., fdydst)
    endif

    !---------------------------------------------------------------
    ! CO2: dcam = dt * 1/(hatm*ra2a) * conv
    !---------------------------------------------------------------
    do j=1,jm
      do i=1,im
        fac(i,j) = 1._wp/(hatm*ra2a(i,j))
      enddo
    enddo
    call zonal_diffuse(cam, diffx, fac, .true., fdxco2)
    call merid_diffuse(cam, diffy, fac, .true., fdyco2)

    return

  end subroutine diffuse_impl


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  z o n a l _ d i f f u s e
  !   Purpose    :  one implicit zonal diffusion step of a 2D column field x,
  !                 solving directly for the increment dprog:
  !                   (I - dt*fac*L) dprog = dt*fac*conv_expl ,  x := x + dprog
  !                 conv_expl = (fdx(i)-fdx(i+1))/sqr is adifa's zonal
  !                 diffusive flux convergence; L (the
  !                 column Laplacian, cyclic per latitude row) is on the LHS only
  !                 and provides the unconditionally-stable implicit damping of
  !                 the stiff polar grid-scale modes. Conservative (flux form).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine zonal_diffuse(x, diffx, fac, mass_weight, fdx)

    implicit none

    real(wp), intent(inout) :: x(:,:)
    real(wp), intent(in   ) :: diffx(:,:)
    real(wp), intent(in   ) :: fac(:,:)
    logical,  intent(in   ) :: mass_weight
    real(wp), intent(in   ) :: fdx(:,:)   ! (imc,jm) adifa exact zonal diffusive flux

    integer :: i, j, ip
    real(wp) :: g
    real(wp) :: dxf(im), wx(im)
    real(wp) :: a(im), b(im), c(im), r(im), sol(im)


    !$omp parallel do private(i,j,ip,g,dxf,wx,a,b,c,r,sol)
    do j=1,jm

      ! zonal face weights and diffusion coefficients (face i lies west of cell i)
      do i=1,im
        if (mass_weight) then
          wx(i) = plx_trop(i,j)
        else
          wx(i) = 1._wp
        endif
        dxf(i) = diffx(i,j)*dy*wx(i)/dxt(j)   ! kg/s per unit tracer gradient
      enddo

      do i=1,im
        ip = i+1
        if (ip.gt.im) ip = 1
        g = tstep*fac(i,j)/sqr(i,j)
        a(i) = -g*dxf(i)            ! couples to x_{i-1} (cyclic for i=1)
        c(i) = -g*dxf(ip)           ! couples to x_{i+1} (cyclic for i=im)
        b(i) = 1._wp + g*(dxf(i)+dxf(ip))
        r(i) = tstep*fac(i,j)*(fdx(i,j) - fdx(i+1,j))/sqr(i,j)   ! (I-D) dprog = dt*fac*conv_expl
      enddo

      call cyclic_tridiag_solve(a, b, c, r, sol, im)

      do i=1,im
        x(i,j) = x(i,j) + sol(i)    ! x + dprog
      enddo

    enddo
    !$omp end parallel do

    return

  end subroutine zonal_diffuse


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m e r i d _ d i f f u s e
  !   Purpose    :  one implicit meridional diffusion step of a 2D column field
  !                 x, solving directly for the increment dprog:
  !                   (I - dt*fac*L) dprog = dt*fac*conv_expl ,  x := x + dprog
  !                 conv_expl = (fdy(j+1)-fdy(j))/sqr is adifa's meridional
  !                 diffusive flux convergence; L (column Laplacian, tridiagonal
  !                 per longitude, no-flux poles) is on the LHS only, providing
  !                 the implicit stabilization. Conservative (flux form).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine merid_diffuse(x, diffy, fac, mass_weight, fdy)

    implicit none

    real(wp), intent(inout) :: x(:,:)
    real(wp), intent(in   ) :: diffy(:,:)
    real(wp), intent(in   ) :: fac(:,:)
    logical,  intent(in   ) :: mass_weight
    real(wp), intent(in   ) :: fdy(:,:)   ! (im,jmc) adifa exact meridional diffusive flux

    integer :: i, j, jp
    real(wp) :: g
    real(wp) :: dyf(jm+1), wy(jm)
    real(wp) :: a(jm), b(jm), c(jm), r(jm), sol(jm)


    !$omp parallel do private(i,j,jp,g,dyf,wy,a,b,c,r,sol)
    do i=1,im

      ! meridional face weights and diffusion coefficients (face j lies south of
      ! cell j); diffy is zero at j=1 and j=jmc -> automatic no-flux at the poles
      do j=1,jm
        if (mass_weight) then
          wy(j) = ply_trop(i,j)
        else
          wy(j) = 1._wp
        endif
        dyf(j) = diffy(i,j)*dxu(j)*wy(j)/dy
      enddo
      dyf(jm+1) = 0._wp   ! no flux through the pole face j=jmc

      do j=1,jm
        jp = j+1
        g = tstep*fac(i,j)/sqr(i,j)
        a(j) = -g*dyf(j)            ! couples to x_{j-1} (=0 at j=1)
        c(j) = -g*dyf(jp)           ! couples to x_{j+1} (=0 at j=jm)
        b(j) = 1._wp + g*(dyf(j)+dyf(jp))
        r(j) = tstep*fac(i,j)*(fdy(i,j+1) - fdy(i,j))/sqr(i,j)   ! (I-D) dprog = dt*fac*conv_expl
      enddo

      call tridiag_solve(a, b, c, r, sol, jm)

      do j=1,jm
        x(i,j) = x(i,j) + sol(j)    ! x + dprog
      enddo

    enddo
    !$omp end parallel do

    return

  end subroutine merid_diffuse

end module diffuse_impl_mod

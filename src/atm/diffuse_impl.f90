!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : d i f f u s e _ i m p l
!
!  Purpose : implicit horizontal diffusion (zonal + meridional, ADI) of the
!            prognostic atmospheric column variables (dry static energy /
!            temperature, column water, dust and CO2).
!
!            EXACT 2D EQUIVALENT. Because the diffusivity has a fixed vertical
!            profile (uniform, or exp(-z/h_diff) folded into the face masses),
!            the column-integrated zonal diffusive convergence is reproduced
!            exactly from 2D quantities by solving, per field,
!
!              (I - dt*fac*L) * dprog = dt*fac*conv_expl ,   prog := prog + dprog
!
!            where conv_expl = (fdx(i)-fdx(i+1))/sqr is adifa's EXACT zonal
!            diffusive flux convergence (the original layer-resolved physics; the
!            gad*z term cancels because adifa differences tp at common height
!            levels and zeroes sub-surface layers via the face mass dplx_dif), and L
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
  use atm_params, only : tstep, hatm, ra, cp, l_dust, f_diff_wtr, diff_stab
  use atm_grid, only : im, jm, jmc
  use atm_grid, only : plx_dif, ply_dif, sqr, dxt, dxu, dy, cheat
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
      ra2a, hdust, &                                                                      ! in
      fdxdse, fdxwtr, fdxdst, fdxco2, fdydse, fdywtr, fdydst, fdyco2, &                   ! inout
      tam, wcon, dam, cam, &                                                              ! inout
      convwtr_dif)                                                                        ! out

    implicit none

    real(wp), intent(in   ) :: diffx(:,:), diffy(:,:)
    real(wp), intent(in   ) :: ra2a(:,:)
    real(wp), intent(in   ) :: hdust(:,:)
    ! in: adifa's explicit diffusive fluxes of the current state; out: the fluxes actually
    ! realized by the implicit solve, whose convergence is exactly the applied increment
    real(wp), intent(inout) :: fdxdse(:,:), fdydse(:,:)   ! DSE diffusive fluxes (zonal/merid)
    real(wp), intent(inout) :: fdxwtr(:,:), fdywtr(:,:)   ! water diffusive fluxes
    real(wp), intent(inout) :: fdxdst(:,:), fdydst(:,:)   ! dust diffusive fluxes
    real(wp), intent(inout) :: fdxco2(:,:), fdyco2(:,:)   ! CO2 diffusive fluxes

    real(wp), intent(inout) :: tam(:,:)
    real(wp), intent(inout) :: wcon(:,:)
    real(wp), intent(inout) :: dam(:,:)
    real(wp), intent(inout) :: cam(:,:)

    real(wp), intent(out  ) :: convwtr_dif(:,:)

    integer :: i, j, n, nt
    real(wp), dimension(im,jm) :: fac_dse, fac_wtr, fac_dst, fac_co2
    real(wp) :: heff
    real(wp) :: wcon_old(im,jm)
    real(wp) :: diffx_wtr(size(diffx,1),size(diffx,2)), diffy_wtr(size(diffy,1),size(diffy,2))   ! water vapour diffusivities (EKE part x f_diff_wtr, as in adifa)


    ! implicit diffusion (zonal sweep then meridional sweep) 
    ! (I - dt*fac*L)*dprog = dt*fac*conv_expl, prog += dprog,
    ! where conv_expl = (fdx(i)-fdx(i+1))/sqr is adifa's exact diffusive flux
    ! convergence. L (the column Laplacian) is
    ! on the LHS only and provides the unconditionally-stable implicit damping of the
    ! stiff polar grid-scale modes; it does not change what is diffused.

    ! Conversion factors from the flux convergence to the tendency of each column variable.
    ! dry static energy: dtam = dt * cp/cheat * conv.  The transport carries DSE (cp*T per
    ! unit mass) and the reservoir is the column enthalpy cheat=pzsa*amas*cp, so the ratio
    ! is simply 1/M, M being the local column mass.  Same conversion time_step applies to
    ! the advective convdse (deba/cheat), so the diffusion, the advection and
    ! the column energy budget all use one and the same heat capacity.
    ! column water content: 1
    ! dust: ddam = dt * 1/(heff*ra) * conv_expl
    ! CO2: dcam = dt * 1/(hatm*ra2a) * conv
    !$omp parallel do private(i,j,heff)
    do j=1,jm
      do i=1,im
        fac_dse(i,j) = cp/cheat(i,j)
        fac_wtr(i,j) = 1._wp
        if (l_dust) then
          heff = hdust(i,j)*hatm/(hdust(i,j)+hatm)   ! m
          fac_dst(i,j) = 1._wp/(heff*ra)
        endif
        fac_co2(i,j) = 1._wp/(hatm*ra2a(i,j))
        wcon_old(i,j) = wcon(i,j)
      enddo
    enddo
    !$omp end parallel do

    diffx_wtr = f_diff_wtr*diffx
    diffy_wtr = f_diff_wtr*diffy

    ! The rows (zonal sweep) and then the columns (meridional sweep) of all tracers are
    ! independent, so each sweep is one parallel loop over tracers x rows (columns).
    ! Tracers: 1 dse, 2 water, 3 CO2, 4 dust (only if l_dust).
    nt = 3
    if (l_dust) nt = 4

    !$omp parallel do private(n,j)
    do n=1,nt*jm
      j = mod(n-1,jm)+1
      select case ((n-1)/jm+1)
      case (1)
        call zonal_row(j, tam, diffx, fac_dse, .true., fdxdse)
      case (2)
        call zonal_row(j, wcon, diffx_wtr, fac_wtr, .false., fdxwtr)
      case (3)
        call zonal_row(j, cam, diffx, fac_co2, .true., fdxco2)
      case (4)
        call zonal_row(j, dam, diffx, fac_dst, .true., fdxdst)
      end select
    enddo
    !$omp end parallel do

    !$omp parallel do private(n,i)
    do n=1,nt*im
      i = mod(n-1,im)+1
      select case ((n-1)/im+1)
      case (1)
        call merid_col(i, tam, diffy, fac_dse, .true., fdydse)
      case (2)
        call merid_col(i, wcon, diffy_wtr, fac_wtr, .false., fdywtr)
      case (3)
        call merid_col(i, cam, diffy, fac_co2, .true., fdyco2)
      case (4)
        call merid_col(i, dam, diffy, fac_dst, .true., fdydst)
      end select
    enddo
    !$omp end parallel do

    ! realized implicit diffusive moisture convergence (zonal+meridional, for precip)
    !$omp parallel do private(i,j)
    do j=1,jm
      do i=1,im
        convwtr_dif(i,j) = (wcon(i,j) - wcon_old(i,j))/tstep   ! kg/m2/s
      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine diffuse_impl


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  z o n a l _ r o w
  !   Purpose    :  one implicit zonal diffusion step of latitude row j of a 2D column field x,
  !                 solving directly for the increment dprog:
  !                   (I - dt*fac*L) dprog = dt*fac*conv_expl ,  x := x + dprog
  !                 conv_expl = (fdx(i)-fdx(i+1))/sqr is adifa's zonal
  !                 diffusive flux convergence; L (the
  !                 column Laplacian, cyclic per latitude row) is on the LHS only
  !                 and provides the unconditionally-stable implicit damping of
  !                 the stiff polar grid-scale modes. Conservative (flux form).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine zonal_row(j, x, diffx, fac, mass_weight, fdx)

    implicit none

    integer,  intent(in   ) :: j          ! latitude row
    real(wp), intent(inout) :: x(:,:)
    real(wp), intent(in   ) :: diffx(:,:)
    real(wp), intent(in   ) :: fac(:,:)
    logical,  intent(in   ) :: mass_weight
    real(wp), intent(inout) :: fdx(:,:)   ! (imc,jm) zonal diffusive flux: in explicit (adifa), out realized

    integer :: i, ip, imi
    real(wp) :: g
    real(wp) :: dxf(im), wx(im)
    real(wp) :: a(im), b(im), c(im), r(im), sol(im)


      ! zonal face weights and diffusion coefficients (face i lies west of cell i)
      do i=1,im
        if (mass_weight) then
          wx(i) = plx_dif(i,j)
        else
          wx(i) = 1._wp
        endif
        ! zonal diffusivity = EKE macro-diffusivity diffx + uniform stabilisation diffusivity diff_stab, as in adifa
        dxf(i) = (diffx(i,j)+diff_stab)*dy*wx(i)/dxt(j)   ! kg/s per unit tracer gradient
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

      ! realized flux through face i (west of cell i): the explicit flux plus the implicit
      ! correction of the LHS operator, F_i = fdx_i + dxf_i*(dprog_{i-1}-dprog_i), so that
      ! dprog_i = dt*fac*(F_i - F_{i+1})/sqr exactly
      do i=1,im
        imi = i-1
        if (imi.lt.1) imi = im
        fdx(i,j) = fdx(i,j) + dxf(i)*(sol(imi)-sol(i))
      enddo
      fdx(im+1,j) = fdx(1,j)

    return

  end subroutine zonal_row


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m e r i d _ c o l
  !   Purpose    :  one implicit meridional diffusion step of longitude column i of a 2D column field
  !                 x, solving directly for the increment dprog:
  !                   (I - dt*fac*L) dprog = dt*fac*conv_expl ,  x := x + dprog
  !                 conv_expl = (fdy(j+1)-fdy(j))/sqr is adifa's meridional
  !                 diffusive flux convergence; L (column Laplacian, tridiagonal
  !                 per longitude, no-flux poles) is on the LHS only, providing
  !                 the implicit stabilization. Conservative (flux form).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine merid_col(i, x, diffy, fac, mass_weight, fdy)

    implicit none

    integer,  intent(in   ) :: i          ! longitude column
    real(wp), intent(inout) :: x(:,:)
    real(wp), intent(in   ) :: diffy(:,:)
    real(wp), intent(in   ) :: fac(:,:)
    logical,  intent(in   ) :: mass_weight
    real(wp), intent(inout) :: fdy(:,:)   ! (im,jmc) meridional diffusive flux: in explicit (adifa), out realized

    integer :: j, jp
    real(wp) :: g
    real(wp) :: dyf(jm+1), wy(jm)
    real(wp) :: a(jm), b(jm), c(jm), r(jm), sol(jm)


      ! meridional face weights and diffusion coefficients (face j lies south of
      ! cell j); diffy is zero at j=1 and j=jmc -> automatic no-flux at the poles
      do j=1,jm
        if (mass_weight) then
          wy(j) = ply_dif(i,j)
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

      ! realized flux through face j (south of cell j, positive towards j-1),
      ! F_j = fdy_j + dyf_j*(dprog_j-dprog_{j-1}), so that dprog_j = dt*fac*(F_{j+1}-F_j)/sqr
      ! exactly; dyf(1) = dyf(jm+1) = 0 keeps the pole faces at zero
      do j=2,jm
        fdy(i,j) = fdy(i,j) + dyf(j)*(sol(j)-sol(j-1))
      enddo

    return

  end subroutine merid_col

end module diffuse_impl_mod

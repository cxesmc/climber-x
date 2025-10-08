!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : v e l c _ m o d
!
!  Purpose : 3D velocity field 
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
module velc_mod

  use precision, only : wp
  use ocn_grid, only : maxi, maxj, maxk, k1, f_pbl, c, rc, cv, rcv, rdsv, rds2, rdphi, dx, dxv, dy, dz, dza, rh, R_earth
  use ocn_params, only : fcor, fcorv, drag_bcl, rho0, rtv, rtv3, urelax, dt
  use constants, only : g

  implicit none

  private
  public :: velc

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  v e l c
  !   Purpose    :  determine the baroclinic solution of the velocity,
  !              :  add barotropic component and 
  !              :  compute vertical velocity from continuity equation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine velc(rho,dtau_dz2,dtav_dz2,ub, &
                  u)

    implicit none

    real(wp), dimension(:,:,:), intent(in) :: rho
    real(wp), dimension(:,:,:), intent(in) :: dtau_dz2, dtav_dz2
    real(wp), dimension(:,0:,0:), intent(in) :: ub

    real(wp), dimension(:,0:,0:,:), intent(inout) :: u

    real(wp), allocatable, dimension(:,:) :: du_dz    
    real(wp), allocatable, dimension(:,:) :: u1
    real(wp), allocatable, dimension(:,:) :: u_old

    integer :: i, j, k, l, ip1, im1
    real(wp) :: tv, tv1, tv2, tv4, tv5, tsum(2)
    real(wp) :: dy_dz_i, dy_dz_im1, dx_dz_j, dx_dz_jm1


    allocate(du_dz(2,maxk))
    allocate(u1(2,maxk))
    allocate(u_old(2,maxk))


    !$omp parallel do private(i,j,k,l,ip1,im1,tv1,tv2,tv4,tv5,tsum,du_dz,u1,u_old)
    do j=1,maxj

      do i=1,maxi
        ! save old velocity field to be used for relaxation
        u_old(1:2,:) = u(1:2,i,j,1:maxk)
        ! set velocity to zero outside domain
        do k=1,maxk
          if (k.lt.k1(i,j)) then
            u_old(1:2,k) = 0._wp
          endif
        enddo
        tsum(:) = 0._wp
        ip1 = i+1
        if (ip1.eq.maxi+1) ip1=1
        im1 = i-1
        if (im1.eq.0) im1=maxi
        do k=k1(i,j),maxk
          ! compute terms of zonal velocity component
          if (k1(ip1,j).gt.k) then
            tv1 = 0._wp
            tv2 = 0._wp
          else
            tv2 = g/R_earth*(rho(ip1,j,k) - rho(i,j,k))*rdphi*rc(j)  ! kg/s2/m3
            if (max(k1(i,j-1),k1(i,j+1),k1(ip1,j-1),k1(ip1,j+1)).le.k) then
              tv1 = g/R_earth*c(j)*(rho(ip1,j+1,k) - rho(ip1,j-1,k) &
                  + rho(i,j+1,k) - rho(i,j-1,k))*rds2(j)*0.25_wp
            elseif (max(k1(i,j-1),k1(ip1,j-1)).le.k) then
              tv1 = g/R_earth*c(j)*(rho(ip1,j,k) - rho(ip1,j-1,k) &
                  + rho(i,j,k) - rho(i,j-1,k))*rdsv(j-1)*0.5_wp
            elseif (max(k1(i,j+1),k1(ip1,j+1)).le.k) then
              tv1 = g/R_earth*c(j)*(rho(ip1,j+1,k) - rho(ip1,j,k) &
                  + rho(i,j+1,k) - rho(i,j,k))*rdsv(j)*0.5_wp
            else
              tv1 = 0._wp
            endif
          endif

          ! compute terms of meridional velocity component
          if (k1(i,j+1).gt.k) then
            tv4 = 0._wp
            tv5 = 0._wp
          else
            tv4 = g/R_earth*cv(j)*(rho(i,j+1,k) - rho(i,j,k))*rdsv(j)  ! kg/s2/m3
            if (max(k1(im1,j),k1(im1,j+1),k1(ip1,j),k1(ip1,j+1)).le.k) then
              tv5 = g/R_earth*(rho(ip1,j+1,k) - rho(im1,j+1,k) &
              + rho(ip1,j,k) - rho(im1,j,k))*rdphi*0.25*rcv(j)
              elseif (max(k1(im1,j),k1(im1,j+1)).le.k) then
              tv5 = g/R_earth*(rho(i,j+1,k) - rho(im1,j+1,k) &
              + rho(i,j,k) - rho(im1,j,k))*rdphi*0.5*rcv(j)
              elseif(max(k1(ip1,j),k1(ip1,j+1)).le.k)then
              tv5 = g/R_earth*(rho(ip1,j+1,k) - rho(i,j+1,k) &
              + rho(ip1,j,k) - rho(i,j,k))*rdphi*0.5*rcv(j)
            else
              tv5 = 0._wp
            endif
          endif

          ! add second vertical derivative of wind stress in the boundary layer only
          if (k1(ip1,j).le.k) then
            tv1 = tv1 + dtau_dz2(2,i,j)*f_pbl(k)  ! kg/s2/m3
            tv2 = tv2 + dtau_dz2(1,i,j)*f_pbl(k)
          endif
          if (k1(i,j+1).le.k) then
            tv4 = tv4 + dtav_dz2(2,i,j)*f_pbl(k)
            tv5 = tv5 + dtav_dz2(1,i,j)*f_pbl(k)
          endif

          ! vertical derivative of u and v
          du_dz(1,k) = (fcor(j)*tv1 + drag_bcl(1,i,j)*tv2)*rtv(i,j)/rho0  ! 1/s
          du_dz(2,k) = (drag_bcl(2,i,j)*tv4 - fcorv(j)*tv5)*rtv3(i,j)/rho0
          ! integrate
          do l=1,2
            if (k.eq.k1(i,j)) then
              u1(l,k) = 0._wp
            else
              u1(l,k) = u1(l,k-1) + dza(k-1)*(du_dz(l,k) + du_dz(l,k-1))*0.5_wp  ! m/s
              tsum(l) = tsum(l) + dz(k)*u1(l,k)  ! m2/s
            endif
          enddo 
        enddo

        ! add barotropic part and relax
        do k = k1(i,j),maxk
          if (k1(ip1,j).le.k) then
            u1(1,k) = u1(1,k) - tsum(1)*rh(1,i,j) + ub(1,i,j)  ! m/s
            u(1,i,j,k) = urelax*u_old(1,k) + (1._wp - urelax)*u1(1,k)
          else 
            u(1,i,j,k) = 0._wp
          endif
          if (k1(i,j+1).le.k) then
            u1(2,k) = u1(2,k) - tsum(2)*rh(2,i,j) + ub(2,i,j)
            u(2,i,j,k) = urelax*u_old(2,k) + (1._wp - urelax)*u1(2,k)
          else
            u(2,i,j,k) = 0._wp 
          endif
        enddo

      enddo

      ! set boundary conditions 
      do k=k1(1,j),maxk
        u(1,0,j,k) = u(1,maxi,j,k)
      enddo

    enddo
    !$omp end parallel do

    ! calculate vertical velocity w from continuity equation
    do j=1,maxj
      do i=1,maxi
        tv = 0._wp
        ip1 = i+1
        if (ip1.eq.maxi+1) ip1=1
        im1 = i-1
        if (im1.eq.0) im1=maxi
        do k=k1(i,j),maxk
          dy_dz_i = dy*dz(k)
          dy_dz_im1 = dy*dz(k)
          dx_dz_j = dxv(j)*dz(k)
          dx_dz_jm1 = dxv(j-1)*dz(k)
          tv1 = (u(1,i,j,k)*dy_dz_i - u(1,i-1,j,k)*dy_dz_im1) ! m3/s
          tv2 = (u(2,i,j,k)*dx_dz_j - u(2,i,j-1,k)*dx_dz_jm1) ! m3/s
          u(3,i,j,k) = tv - (tv1 + tv2) / (dx(j)*dy)
          tv = u(3,i,j,k)
        enddo
      enddo
    enddo

    deallocate(du_dz)
    deallocate(u1)
    deallocate(u_old)

   return

  end subroutine velc

end module velc_mod

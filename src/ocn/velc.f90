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
  use ocn_grid, only : maxi, maxj, maxk, k1, f_ek, c, rc, cv, rcv, rdsv, rds2, rdphi, dx, dxv, dy, dz, dza, rh, R_earth
  use ocn_params, only : fcor, fcorv, drag, drag_bcl, drag_par, rho0, rtv, rtv3, dt
  use constants, only : g, omega, pi

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
    real(wp), allocatable, dimension(:,:) :: u_ek

    integer :: i, j, k, l, ip1, im1, kb(2)
    real(wp) :: tv, tv1, tv2, tv4, tv5, tsum(2)
    real(wp) :: f0, ufr(2), hface
    real(wp) :: dy_dz_i, dy_dz_im1, dx_dz_j, dx_dz_jm1


    allocate(du_dz(2,maxk))
    allocate(u1(2,maxk))
    allocate(u_ek(2,maxk))

    ! equatorial regularisation of 1/f in the bottom return of the barotropic friction
    f0 = 2._wp*omega*sin(drag_par%lat_fric_ret*pi/180._wp)


    !$omp parallel do private(i,j,k,l,ip1,im1,tv1,tv2,tv4,tv5,tsum,du_dz,u1,u_ek,kb,ufr,hface)
    do j=1,maxj

      do i=1,maxi
        tsum(:) = 0._wp
        ip1 = modulo(i,maxi) + 1
        im1 = modulo(i - 2, maxi) + 1

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

          ! wind stress in the boundary layer: Ekman velocity from the exact layer mean of dtau/dz (dtau_dz2*f_ek),
          ! added to the level directly; sums to tau/(rho0*f) (frictional form with drag_bcl) over the boundary layer
          u_ek(:,k) = 0._wp
          if (k1(ip1,j).le.k) then
            u_ek(1,k) = (fcor(j)*dtau_dz2(2,i,j) + drag_bcl(1,i,j)*dtau_dz2(1,i,j))*f_ek(k)*rtv(i,j)/rho0  ! m/s
          endif
          if (k1(i,j+1).le.k) then
            u_ek(2,k) = (drag_bcl(2,i,j)*dtav_dz2(2,i,j) - fcorv(j)*dtav_dz2(1,i,j))*f_ek(k)*rtv3(i,j)/rho0
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

        ! add the Ekman velocity after the shear integration, so it does not accumulate upward
        do k=k1(i,j),maxk
          u1(:,k) = u1(:,k) + u_ek(:,k)
          tsum(:) = tsum(:) + dz(k)*u_ek(:,k)
        enddo

        ! i_fric_ret==1: the part of the barotropic velocity driven by the barotropic friction, r*vb/f (u-face) and
        ! r*ub/f (v-face), is moved from the whole column into the bottom level of the face (bottom Ekman layer), so
        ! that the friction returns the Ekman transport at the bottom, where topographic form stress acts in reality,
        ! instead of evenly over the column. The column transport and hence psi are unchanged.
        ! Momentum balances: f*ub + r*vb = Fy -> u_fric = -r*vb/f;  -f*vb + r*ub = Fx -> v_fric = r*ub/f,
        ! with the velocity of the other component averaged from the 4 neighbouring faces and 1/f -> f/(f^2+f0^2).
        ufr(:) = 0._wp
        kb(:) = maxk
        if (drag_par%i_fric_ret.eq.1) then
          if (rh(1,i,j).gt.0._wp) then
            ufr(1) = -drag(1,i,j)*fcor(j)/(fcor(j)**2+f0**2) &
                     * 0.25_wp*(ub(2,i,j)+ub(2,ip1,j)+ub(2,i,j-1)+ub(2,ip1,j-1))
            kb(1) = max(k1(i,j),k1(ip1,j))
          endif
          if (j.lt.maxj .and. rh(2,i,j).gt.0._wp) then
            ufr(2) = drag(2,i,j)*fcorv(j)/(fcorv(j)**2+f0**2) &
                     * 0.25_wp*(ub(1,i,j)+ub(1,im1,j)+ub(1,i,j+1)+ub(1,im1,j+1))
            kb(2) = max(k1(i,j),k1(i,j+1))
          endif
          do l=1,2
            if (ufr(l).ne.0._wp) then
              hface = sum(dz(kb(l):maxk))
              u1(l,kb(l)) = u1(l,kb(l)) + ufr(l)*hface/dz(kb(l))
              u1(l,kb(l):maxk) = u1(l,kb(l):maxk) - ufr(l)
            endif
          enddo
        endif

        ! add barotropic part
        do k = k1(i,j),maxk
          if (k1(ip1,j).le.k) then
            u1(1,k) = u1(1,k) - tsum(1)*rh(1,i,j) + ub(1,i,j)  ! m/s
            u(1,i,j,k) = u1(1,k)
          else 
            u(1,i,j,k) = 0._wp
          endif
          if (k1(i,j+1).le.k) then
            u1(2,k) = u1(2,k) - tsum(2)*rh(2,i,j) + ub(2,i,j)
            u(2,i,j,k) = u1(2,k)
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
    ! (columns are independent; the k-integration recurrence on tv stays within a column)
    !$omp parallel do private(i,j,k,tv,tv1,tv2,dy_dz_i,dy_dz_im1,dx_dz_j,dx_dz_jm1)
    do j=1,maxj
      do i=1,maxi
        tv = 0._wp
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
    !$omp end parallel do

    deallocate(du_dz)
    deallocate(u1)
    deallocate(u_ek)

   return

  end subroutine velc

end module velc_mod

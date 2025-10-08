!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b e r i n g _ m o d
!
!  Purpose : parameterisation of Bering Strait throughflow
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
module bering_mod

  use precision, only : wp
  use constants, only : g
  use ocn_params, only : c_bering, dt
  use ocn_grid, only : maxk
  use ocn_grid, only : dx, dy, dz, zw, zro

  !$  use omp_lib

  implicit none

  private
  public :: bering

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b e r i n g
  !   Purpose    :  compute Bering Strait throughflow and tracer transport 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bering(A_bering, f_ocn, ssh, saln0, &
                    sal, bering_tf, bering_fw)

    implicit none

    real(wp), intent(in) :: A_bering    ! Bering Strait cross-sectional area (m2)
    real(wp), dimension(:,:), intent(in) :: f_ocn   ! ocean fraction (1)
    real(wp), dimension(:,:), intent(in) :: ssh     ! sea surface height (m)
    real(wp), intent(in) :: saln0    ! reference salinity (psu)

    real(wp), dimension(:,:,:), intent(inout) :: sal

    real(wp), intent(out) :: bering_tf
    real(wp), intent(out) :: bering_fw

    real(wp), parameter :: h = 50._wp           ! indicative depth of Strait (m)
    real(wp), parameter :: f_cor = 1.4e-4_wp    ! Coriolis parameter

    integer :: i, j, k, kk, k1_h, nk
    integer, parameter :: ni = 3, nj_Arctic = 2, nj_Pacific = 2
    integer, dimension(ni) :: i_idx = (/1,2,3/)
    integer, dimension(nj_Arctic)  :: j_idx_Arctic = (/33,34/)
    integer, dimension(nj_Pacific) :: j_idx_Pacific = (/30,31/)
    real(wp) :: A_Arctic, A_Pacific, ssh_Arctic, ssh_Pacific, Qmax, fy, tmp, z_tot
    real(wp), dimension(:), allocatable :: sal_Arctic, sal_Pacific

    real(wp), parameter :: A_bering_present = 5383306._wp   ! present-day Bering Strait cross-sectional area (m2)


    k1_h = maxk
    tmp = 1000._wp
    do k=maxk,1,-1
      if (abs(zw(k-1)+h).lt.tmp) then
        k1_h = k
        tmp = abs(zw(k-1)+h)
      endif
    enddo
    nk = maxk-k1_h+1
    z_tot = -zw(k1_h-1) 

    ! compute average SSH and tracer concentrations over Arctic region
    allocate(sal_Arctic(nk))
    A_Arctic = 0._wp
    ssh_Arctic = 0._wp
    sal_Arctic(:) = 0._wp
    do i=1,ni
      do j=1,nj_Arctic
        if (f_ocn(i_idx(i),j_idx_Arctic(j)).gt.0._wp) then
          A_Arctic = A_Arctic+f_ocn(i_idx(i),j_idx_Arctic(j))
          ssh_Arctic = ssh_Arctic + ssh(i_idx(i),j_idx_Arctic(j))*f_ocn(i_idx(i),j_idx_Arctic(j))
          kk = 0
          do k=k1_h,maxk
            kk=kk+1
            sal_Arctic(kk) = sal_Arctic(kk) + sal(i_idx(i),j_idx_Arctic(j),k)*f_ocn(i_idx(i),j_idx_Arctic(j))
          enddo
        endif
      enddo
    enddo
    ! compute average SSH and tracer concentrations over Pacific region
    allocate(sal_Pacific(nk))
    A_Pacific = 0._wp
    ssh_Pacific = 0._wp
    sal_Pacific(:) = 0._wp
    do i=1,ni
      do j=1,nj_Pacific
        if (f_ocn(i_idx(i),j_idx_Pacific(j)).gt.0._wp) then
          A_Pacific = A_Pacific+f_ocn(i_idx(i),j_idx_Pacific(j))
          ssh_Pacific = ssh_Pacific + ssh(i_idx(i),j_idx_Pacific(j))*f_ocn(i_idx(i),j_idx_Pacific(j))
          kk = 0
          do k=k1_h,maxk
            kk=kk+1
            sal_Pacific(kk) = sal_Pacific(kk) + sal(i_idx(i),j_idx_Pacific(j),k)*f_ocn(i_idx(i),j_idx_Pacific(j))
          enddo
        endif
      enddo
    enddo

    if (A_Arctic.gt.0._wp .and. A_Pacific.gt.0._wp) then

      ssh_Arctic = ssh_Arctic/A_Arctic
      ssh_Pacific= ssh_Pacific/A_Pacific
      sal_Arctic = sal_Arctic/A_Arctic
      sal_Pacific= sal_Pacific/A_Pacific

      ! maximum water flow through the strait, after Goosse 1997
      Qmax = c_bering * A_bering/A_bering_present * g*h/f_cor * (ssh_Pacific-ssh_Arctic)  ! m/s2*m*s * m = m3/s
      !print *,'ssh_Arctic,ssh_Pacific,Qmax',ssh_Arctic,ssh_Pacific,Qmax*1e-6
      !print *,'sal_Arctic',sal_Arctic
      !print *,'sal_Pacific',sal_Pacific

      ! salinity transport 
      kk=0
      do k=k1_h,maxk
        kk=kk+1
        ! northward salt flux using upstream salinity
        if (Qmax.gt.0._wp) then
          fy = Qmax*dz(k)/z_tot * (sal_Pacific(kk)-saln0) * dt  ! m3/s * psu * s = m3 * psu
        else
          fy = Qmax*dz(k)/z_tot * (sal_Arctic(kk)-saln0) * dt  ! m3/s * psu * s = m3 * psu
        endif
        !if (kk.eq.3) then
        !  print *,'fy',fy
        !  print *,'sal Arctic before',sal(1:3,33,k)
        !  print *,'sal Pacific before',sal(1:3,31,k)
        !endif
        do i=1,ni
          do j=1,nj_Pacific
            if (f_ocn(i_idx(i),j_idx_Pacific(j)).gt.0._wp) then
              sal(i_idx(i),j_idx_Pacific(j),k) = sal(i_idx(i),j_idx_Pacific(j),k) &
                - fy*f_ocn(i_idx(i),j_idx_Pacific(j))/A_Pacific / (dx(j_idx_Pacific(j))*dy*dz(k)*f_ocn(i_idx(i),j_idx_Pacific(j)))
            endif
          enddo
        enddo
        do i=1,ni
          do j=1,nj_Arctic
            if (f_ocn(i_idx(i),j_idx_Arctic(j)).gt.0._wp) then
              sal(i_idx(i),j_idx_Arctic(j),k) = sal(i_idx(i),j_idx_Arctic(j),k) &
                + fy*f_ocn(i_idx(i),j_idx_Arctic(j))/A_Arctic / (dx(j_idx_Arctic(j))*dy*dz(k)*f_ocn(i_idx(i),j_idx_Arctic(j)))
            endif
          enddo
        enddo
        !if (kk.eq.3) then
        !  print *,'sal Arctic after',sal(1:3,33,k)
        !  print *,'sal Pacific after',sal(1:3,31,k)
        !endif
      enddo

    else

      Qmax = 0._wp

    endif

    bering_tf = Qmax*1.e-6_wp  ! Sv
    if (Qmax.gt.0._wp) then
      bering_fw = Qmax*(1._wp-sum(sal_Pacific*dz(k1_h:maxk))/z_tot/saln0)*1.e-6_wp  ! Sv
    else
      bering_fw = Qmax*(1._wp-sum(sal_Arctic*dz(k1_h:maxk))/z_tot/saln0)*1.e-6_wp  ! Sv
    endif


    return

  end subroutine bering

end module bering_mod 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : h y p s o _ t o p o _ m o d
!
!  Purpose : sub-grid hypsometry of the seafloor for the marine sediments.
!            For each coarse cell it builds an area-weighted histogram of the
!            high-resolution ocean depth (depth = -z_bed, so depth 0 is the
!            current sea surface since z_bed is referenced to sea level and
!            updated with sea-level change). Recomputed every geo update; the
!            coupler later re-bins this depth histogram onto the (time-varying)
!            ocean levels to obtain the per-level seafloor area fractions used
!            by the depth-resolved sediment model.
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
module hypso_topo_mod

  use precision, only : wp
  use constants, only : pi, r_earth
  use geo_grid, only : ni, nj, ni_topo, nj_topo
  use geo_grid, only : i0_topo, i1_topo, j0_topo, j1_topo, area, lon_topo, lat_topo

  implicit none

  ! fixed depth bins of the per-cell hypsometric histogram. Bin b (1..n_hypso)
  ! covers [ (b-1)*dz_hypso , b*dz_hypso ) with centre (b-0.5)*dz_hypso. Depths
  ! beyond the range are clamped into the deepest bin.
  real(wp), parameter :: depth_hypso_max = 6500._wp   !! deepest bin edge [m]
  real(wp), parameter :: dz_hypso        =   10._wp   !! bin width [m]
  integer,  parameter :: n_hypso = nint(depth_hypso_max/dz_hypso)

  ! coral photic zone: top z_coral metres, sampled in bins of dz_coral. Bin lev
  ! (1..n_coral_fine) covers depth ((lev-1)*dz_coral, lev*dz_coral]. The coral
  ! arrays on the bgc side are allocatable and sized from n_coral_fine (passed via
  ! the bgc grid), so dz_coral can be changed here without touching the bgc code.
  real(wp), parameter :: z_coral  = 50._wp   !! coral photic-zone depth [m]
  real(wp), parameter :: dz_coral =  5._wp   !! coral vertical sampling [m]
  integer,  parameter :: n_coral_fine = nint(z_coral/dz_coral)

  ! the coral topographic (slope) factor only drifts via slow bed changes (GIA),
  ! so it is recomputed only every n_topo_update geo updates to save time.
  integer, parameter :: n_topo_update = 10

  private
  public :: hypso_topo, hypso_topo_factor
  public :: dz_hypso, n_hypso, n_coral_fine, dz_coral, z_coral, n_topo_update, hypso_depth_centre

contains

  ! bin-centre depth [m] of histogram bin b
  elemental function hypso_depth_centre(b) result(z)
    integer, intent(in) :: b
    real(wp) :: z
    z = (real(b,wp)-0.5_wp)*dz_hypso
  end function hypso_depth_centre


  subroutine hypso_topo(z_bed, hypso_f_depth, hypso_f_fine)

    implicit none

    real(wp), intent(in)  :: z_bed(:,:)             !! high-res bed elevation [m] (0 = sea level)
    real(wp), intent(out) :: hypso_f_depth(:,:,:)   !! (ni,nj,n_hypso) area fraction per coarse depth bin (sediment classes)
    real(wp), intent(out) :: hypso_f_fine(:,:,:)    !! (ni,nj,n_coral_fine) area fraction per 1 m depth bin in the coral photic zone

    integer :: i, j, ii, jj, b, bf
    real(wp) :: depth, atot

    !$omp parallel do collapse(2) private(i,j,ii,jj,b,bf,depth,atot)
    do j=1,nj
      do i=1,ni
        hypso_f_depth(i,j,:) = 0._wp
        hypso_f_fine(i,j,:)  = 0._wp
        atot = 0._wp
        do jj=j0_topo(j),j1_topo(j)
          do ii=i0_topo(i),i1_topo(i)
            atot = atot + area(ii,jj)
            if (z_bed(ii,jj)<0._wp) then
              ! area-weighted histogram of ocean depth
              depth = -z_bed(ii,jj)
              ! coarse bins (sediment depth classes), full depth
              b = floor(depth/dz_hypso) + 1
              b = max(1,min(n_hypso,b))
              hypso_f_depth(i,j,b) = hypso_f_depth(i,j,b) + area(ii,jj)
              ! fine bins (coral photic zone), only the top z_coral metres
              if (depth<=z_coral) then
                bf = max(1,min(n_coral_fine,ceiling(depth/dz_coral)))   ! bin lev covers depth ((lev-1)*dz,lev*dz]
                hypso_f_fine(i,j,bf) = hypso_f_fine(i,j,bf) + area(ii,jj)
              endif
            endif
          enddo
        enddo
        ! normalise to fraction of total cell area (sum over bins = ocean fraction)
        if (atot>0._wp) then
          hypso_f_depth(i,j,:) = hypso_f_depth(i,j,:)/atot
          hypso_f_fine(i,j,:)  = hypso_f_fine(i,j,:)/atot
        endif
      enddo
    enddo
    !$omp end parallel do

  end subroutine hypso_topo


  subroutine hypso_topo_factor(z_bed, hypso_f_topo)
    !! Seabed-slope reef-suitability factor (Kleypas 1997), binned by depth below
    !! the current sea surface in the same 1 m bins as the coral area curve. The
    !! slope is a property of the bed, but its depth binning tracks sea level/GIA,
    !! so this is recomputed (throttled) from the current high-res bed topography.

    implicit none

    real(wp), intent(in)  :: z_bed(:,:)             !! high-res bed elevation [m] (0 = sea level)
    real(wp), intent(out) :: hypso_f_topo(:,:,:)    !! (ni,nj,n_coral_fine) mean topo factor per 1 m depth bin

    integer :: i, j, ii, jj, iii, jjj, bf
    real(wp) :: phi1, phi2, theta1, theta2, tmp, dist, alpha, depth
    real(wp), parameter :: deg2rad = pi/180._wp
    real(wp), dimension(:,:), allocatable :: tfac
    real(wp), dimension(n_coral_fine) :: asum

    allocate(tfac(ni_topo,nj_topo))

    ! per-pixel topographic (slope) factor from the 8-neighbour bed gradient
    !$omp parallel do collapse(2) private(i,j,ii,jj,iii,jjj,phi1,phi2,theta1,theta2,tmp,dist,alpha)
    do j=1,nj_topo
      do i=1,ni_topo
        alpha = 0._wp
        do jj=j-1,j+1
          do ii=i-1,i+1
            iii = ii
            if (iii.eq.0) iii = ni_topo
            if (iii.eq.ni_topo+1) iii = 1
            jjj = max(1,min(nj_topo,jj))
            if (i.eq.iii .and. j.eq.jjj) cycle
            phi1 = (90._wp - lat_topo(jjj))*deg2rad
            phi2 = (90._wp - lat_topo(j))*deg2rad
            theta1 = lon_topo(iii)*deg2rad
            theta2 = lon_topo(i)*deg2rad
            tmp = sin(phi1)*sin(phi2)*cos(theta1-theta2) + cos(phi1)*cos(phi2)
            dist = acos(max(-1._wp,min(1._wp,tmp)))*r_earth
            if (dist.gt.0._wp) alpha = alpha + atan((z_bed(iii,jjj)-z_bed(i,j))/dist)
          enddo
        enddo
        alpha = min(alpha,1.7_wp)
        alpha = max(alpha,0.013_wp)
        tfac(i,j) = log(alpha*100._wp)/5._wp
      enddo
    enddo
    !$omp end parallel do

    ! area-weighted mean factor per 1 m depth bin in each coarse cell
    !$omp parallel do collapse(2) private(i,j,ii,jj,bf,depth,asum)
    do j=1,nj
      do i=1,ni
        hypso_f_topo(i,j,:) = 0._wp
        asum = 0._wp
        do jj=j0_topo(j),j1_topo(j)
          do ii=i0_topo(i),i1_topo(i)
            if (z_bed(ii,jj)<0._wp) then
              depth = -z_bed(ii,jj)
              if (depth<=z_coral) then
                bf = max(1,min(n_coral_fine,ceiling(depth/dz_coral)))
                hypso_f_topo(i,j,bf) = hypso_f_topo(i,j,bf) + tfac(ii,jj)*area(ii,jj)
                asum(bf) = asum(bf) + area(ii,jj)
              endif
            endif
          enddo
        enddo
        where (asum>0._wp) hypso_f_topo(i,j,:) = hypso_f_topo(i,j,:)/asum
      enddo
    enddo
    !$omp end parallel do

    deallocate(tfac)

  end subroutine hypso_topo_factor

end module hypso_topo_mod

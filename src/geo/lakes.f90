!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l a k e _ m o d
!
!  Purpose : dynamic lake model
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
module lake_mod

  use precision, only : wp 
  use constants, only : rho_w
  use timer, only : sec_year, n_year_geo, n_accel
  use geo_params, only : i_lakes, lake_area_crit, lake_sea_z_crit
  use geo_def, only : n_lev, lake_type
  !$ use omp_lib

  implicit none

  real(wp), parameter :: lake_depth_min = 100._wp     ! m, minimum lake depth

  integer :: ni, nj
  integer, dimension(:,:), allocatable :: map

  private
  public :: lake_pot, update_lake_mask, lake_type, ocean_lake_mask_level

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l a k e _ p o t
  !   Purpose    :  determination of potential lakes 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lake_pot(z_topo, z_topo_fill, mask, area, &
      mask_lake_pot, lake)

    implicit none

    real(wp), intent(in) :: z_topo(:,:)
    real(wp), intent(in) :: z_topo_fill(:,:)
    integer, intent(in) :: mask(:,:)
    real(wp), intent(in) :: area(:,:)
    integer, intent(inout) :: mask_lake_pot(:,:)
    type(lake_type), allocatable, intent(inout) :: lake(:)

    integer :: i, j, n, k, n_old, n_new
    logical :: flag_ilakes
    integer, parameter :: n_lakes_max = 10000
    integer :: n_lakes, n_lakes_old
    integer :: ncells
    integer, parameter :: ncells_max = 400000
    integer, dimension(n_lakes_max) :: lake_i, lake_j 
    real(wp), dimension(n_lakes_max) :: lake_area, lake_vol
    real(wp), dimension(n_lakes_max) :: lake_area_pot, lake_vol_pot
    real(wp), dimension(n_lakes_max) :: lake_area_min, lake_vol_min
    real(wp), dimension(n_lakes_max,0:n_lev) :: lake_area_lev, lake_vol_lev, lake_z_lev
    real(wp), dimension(n_lakes_max) :: dh_p_e, dh_run, runoff_diag
    integer, dimension(:,:), allocatable :: mask_lake_pot_old
    logical, dimension(:,:), allocatable :: lake_sea
    integer, dimension(:,:), allocatable :: lake_mask
    real(wp), dimension(:,:), allocatable :: lake_depth
    integer, dimension(ncells_max) :: il, jl
    real(wp) :: depth, depth_max
    real(wp) :: zb, z_bot, z_k, dz, z, vol_k, area_k
    logical :: found_lake
    !$ real(wp) :: time1,time2


    ! get grid size
    ni = size(z_topo,1)
    nj = size(z_topo,2)

    ! allocate
    allocate(lake_sea(ni,nj))
    allocate(lake_mask(ni,nj))
    allocate(lake_depth(ni,nj))
    allocate(map(ni,nj))

    ! save old potential lake mask
    mask_lake_pot_old = mask_lake_pot

    !--------------------------------------------------------------
    ! identify and index lakes
    !--------------------------------------------------------------

    ! compute potential lake depth in each cell
    lake_depth = z_topo_fill-z_topo

    ! derive potential lake mask (exclude grounded ice points)
    where (mask.ne.0 .and. lake_depth.gt.0._wp)
      lake_mask = 1
    elsewhere
      lake_mask = 0
    endwhere

    ! identify 'lakes' below sea level 
    where (lake_mask.eq.1 .and. z_topo.le.0._wp) 
      lake_sea = .true.
    elsewhere
      lake_sea = .false.
    endwhere

    ! index all lakes based on lake mask
    where (lake_mask==1) 
      map = 1    ! lake
      mask_lake_pot = -1
    elsewhere
      map = 0     ! no lake
      mask_lake_pot = 0
    endwhere
    n_lakes = 0 ! initialise lake counter

    lake_area_min = 0._wp
    lake_vol_min = 0._wp
    lake_area_pot = 0._wp
    lake_vol_pot = 0._wp

    !!$ time1 = omp_get_wtime()
    do i=1,ni
      do j=1,nj
        if (mask_lake_pot(i,j).eq.-1) then   ! unclaimed lake,found new lake
          ! initialize
          ncells = 0
          !print *
          !print *,'point ',i,j
          call identify_lake(i,j,ncells,il,jl)   ! identify current lake extent
          if (ncells.gt.ncells_max) then
            print *,'ERROR: too many lake cells, increase ncells_max!'
            print *,'ncells, ncells_max',ncells,ncells_max
            stop
          endif
          if (i_lakes.eq.1) then
            ! all lakes active
            flag_ilakes = .true.
          else if (i_lakes.eq.2) then
            ! only 'sea lakes' active (below sea level)
            zb = 10000._wp
            do n=1,ncells
              zb = min(zb,z_topo(il(n),jl(n)))
            enddo
            flag_ilakes = zb.lt.lake_sea_z_crit
          endif
          if (flag_ilakes .and. ncells*area(i,j)>lake_area_crit) then    ! approximate lake area, area threshold for lake to be considered
            ! add new lake
            n_lakes = n_lakes+1 ! increase index
            depth = 0._wp
            lake_vol_min(n_lakes) = 0._wp
            lake_area_min(n_lakes) = 0._wp
            lake_vol_pot(n_lakes) = 0._wp
            lake_area_pot(n_lakes) = 0._wp
            do n=1,ncells
              mask_lake_pot(il(n),jl(n)) = n_lakes
              ! lake deepest point index
              if (lake_depth(il(n),jl(n))>depth) then
                depth = lake_depth(il(n),jl(n))
                lake_i(n_lakes) = il(n)
                lake_j(n_lakes) = jl(n)
              endif
              ! potential lake volume
              lake_vol_pot(n_lakes) = lake_vol_pot(n_lakes) + lake_depth(il(n),jl(n))*area(il(n),jl(n))
              ! potential lake area
              lake_area_pot(n_lakes) = lake_area_pot(n_lakes) + area(il(n),jl(n))
              if (lake_sea(il(n),jl(n))) then
                ! minimum lake volume
                !lake_vol_min(n_lakes) = lake_vol_min(n_lakes) + lake_depth(il(n),jl(n))*area(il(n),jl(n))
                lake_vol_min(n_lakes) = lake_vol_min(n_lakes) + (min(0._wp,z_topo_fill(il(n),jl(n)))-z_topo(il(n),jl(n)))*area(il(n),jl(n))
                ! minimum lake area
                lake_area_min(n_lakes) = lake_area_min(n_lakes) + area(il(n),jl(n))
              endif
            enddo
            ! lake_hypsometry
            depth_max = lake_depth(lake_i(n_lakes),lake_j(n_lakes))   ! depth of deepest point
            z_bot = z_topo(lake_i(n_lakes),lake_j(n_lakes)) ! elevation of deepest point
            dz = depth_max/real(n_lev,wp)
            !$omp parallel do private(k,n,z_k,area_k,vol_k,z)
            do k=1,n_lev
              z_k = z_bot + k*dz
              area_k = 0._wp
              vol_k = 0._wp
              do n=1,ncells
                z = z_topo(il(n),jl(n))
                if (z<z_k) then
                  vol_k = vol_k + (z_k-z)*area(il(n),jl(n))
                  area_k = area_k + area(il(n),jl(n))
                endif
              enddo
              lake_z_lev(n_lakes,k) = z_k
              lake_vol_lev(n_lakes,k) = vol_k
              lake_area_lev(n_lakes,k) = area_k
            enddo
            !$omp end parallel do
            lake_z_lev(n_lakes,0) = z_bot
            lake_vol_lev(n_lakes,0) = 0._wp
            lake_area_lev(n_lakes,0) = 0._wp
            ! minimum lake volume
            lake_vol_min(n_lakes) = lake_vol_min(n_lakes) + lake_vol_lev(n_lakes,1)
            ! minimum lake area
            lake_area_min(n_lakes) = lake_area_min(n_lakes) + lake_area_lev(n_lakes,1)
            !print *
            !print *,'lake number,lake_area,lake_vol',n_lakes,lake_area_pot(n_lakes),lake_vol_pot(n_lakes)
          else
            ! ignore small lake
            do n=1,ncells
              mask_lake_pot(il(n),jl(n)) = 0
            enddo
          endif
        endif
      enddo
    enddo
    !!$ time2 = omp_get_wtime()
    !!$ print *,'1',time2-time1

    print *,'n_lakes',n_lakes

    ! account for changing number of potential lakes when topography is changing
    !!$ time1 = omp_get_wtime()
    if (allocated(lake)) then
      n_lakes_old = size(lake)
      lake_area = 0._wp
      lake_vol = 0._wp
      dh_p_e = 0._wp
      dh_run = 0._wp
      runoff_diag = 0._wp
      ! old-new lake correspondence 
      do n_old=1,n_lakes_old
        found_lake = .false.
        do n_new=1,n_lakes
          !if (lake_i(n_new).eq.lake(n_old)%hires%i .and. lake_j(n_new).eq.lake(n_old)%hires%j) then
          ! fixme!?
          if (count(mask_lake_pot.eq.n_new .and. mask_lake_pot_old.eq.n_old).gt.0) then  ! check for overlap between old and new lake
            ! old potential lake still existing
            found_lake = .true.
            ! transfer area and volume
            lake_area(n_new) = lake(n_old)%area
            lake_vol(n_new) = lake(n_old)%vol
            ! transfer lake depth tendency diagnostics
            dh_p_e(n_new) = lake(n_old)%dh_p_e
            dh_run(n_new) = lake(n_old)%dh_run
            runoff_diag(n_new) = lake(n_old)%runoff_diag
            exit
          endif
        enddo
        if (.not.found_lake) then
          ! todo move old lake water to somewhere to conserve water
          ! lake(n_old)%vol
        endif
      enddo
      deallocate(lake)
      allocate(lake(n_lakes))
    else
      ! fixme, from restart!
      lake_area = 0._wp
      lake_vol = 0._wp
      dh_p_e = 0._wp
      dh_run = 0._wp
      runoff_diag = 0._wp
      allocate(lake(n_lakes))
    endif
    !!$ time2 = omp_get_wtime()
    !!$ print *,'2',time2-time1

    !$omp parallel do private(n)
    do n=1,n_lakes
      lake(n)%area_min = lake_area_min(n)
      lake(n)%vol_min  = lake_vol_min(n)
      lake(n)%area_pot = lake_area_pot(n)
      lake(n)%vol_pot  = lake_vol_pot(n)
      lake(n)%area = lake_area(n)
      lake(n)%area = max(lake(n)%area,lake(n)%area_min)
      lake(n)%vol  = lake_vol(n)
      lake(n)%vol  = max(lake(n)%vol,lake(n)%vol_min)
!      lake(n)%runoff_geo = max(0._wp,lake(n)%vol-lake(n)%vol_pot) * rho_w / (real(n_year_geo,wp)/real(n_accel,wp)*sec_year)   ! m3 * kg/m3 / s = kg/s
      lake(n)%runoff_geo = 0._wp
      lake(n)%vol  = min(lake(n)%vol,lake(n)%vol_pot)
      lake(n)%z_lev   = lake_z_lev(n,:)
      lake(n)%area_lev(0:)  = lake_area_lev(n,0:)
      lake(n)%vol_lev(0:)  = lake_vol_lev(n,0:)
      lake(n)%hires%i = lake_i(n)
      lake(n)%hires%j = lake_j(n)
      ! lake depth 
      lake(n)%depth_max = lake_depth(lake(n)%hires%i,lake(n)%hires%j)   ! depth of deepest point
      lake(n)%z_bot = z_topo(lake(n)%hires%i,lake(n)%hires%j) ! elevation of deepest point
      ! copy also diagnostics 
      lake(n)%dh_p_e  = dh_p_e(n)
      lake(n)%dh_run  = dh_run(n)
      lake(n)%runoff_diag  = runoff_diag(n)
    enddo
    !$omp end parallel do

    deallocate(lake_sea)
    deallocate(lake_mask)
    deallocate(lake_depth)
    deallocate(map)


   return

  end subroutine lake_pot


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  i d e n t i f y _ l a k e
  ! Purpose  :  identify lakes
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  recursive subroutine identify_lake(i,j,ncells,il,jl) 

    implicit none

    integer, intent(in) :: i, j
    integer, intent(inout) :: ncells
    integer, intent(inout) :: il(:), jl(:)

    integer :: im1, ip1, jm1, jp1


    if (map(i,j) == 1) then

      ! This is a contiguous cell - index it
      ncells = ncells+1
      il(ncells) = i
      jl(ncells) = j
      map(i,j) = 0

      ! Now recursively search neighbouring cells
      ! West
      im1 = i-1
      if (im1 > 0) then
        call identify_lake(im1,j,ncells,il,jl)
      endif
      ! East
      ip1 = i+1
      if (ip1 < (ni+1)) then
        call identify_lake(ip1,j,ncells,il,jl)
      endif
      ! South
      jm1 = j-1
      if (jm1 > 0) then
        call identify_lake(i,jm1,ncells,il,jl)
      endif
      ! North
      jp1 = j+1
      if (jp1 < (nj+1)) then
        call identify_lake(i,jp1,ncells,il,jl)
      endif
      ! North-west
      im1 = i-1
      jp1 = j+1
      if (im1 > 0 .and. jp1 < (nj+1)) then
        call identify_lake(im1,jp1,ncells,il,jl)
      endif
      ! South-west
      im1 = i-1
      jm1 = j-1
      if (im1 > 0 .and. jm1 > 0) then
        call identify_lake(im1,jm1,ncells,il,jl)
      endif
      ! North-east
      ip1 = i+1
      jp1 = j+1
      if (ip1 < (ni+1) .and. jp1 < (nj+1)) then
        call identify_lake(ip1,jp1,ncells,il,jl)
      endif
      ! South-east
      ip1 = i+1
      jm1 = j-1
      if (ip1 < (ni+1) .and. jm1 > 0) then
        call identify_lake(ip1,jm1,ncells,il,jl)
      endif

      ! wrap east-west
      ! Jump west
      if (i == 1) then
        call identify_lake(ni,j,ncells,il,jl)
      endif
      ! Jump east
      if (i == ni) then
        call identify_lake(1,j,ncells,il,jl)
      endif

    endif

  end subroutine identify_lake


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  u p d a t e _ l a k e _ m a s k
  !   Purpose    :  update lake mask based on actual lake water volume
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine update_lake_mask(z_topo, area, lake, mask_lake_pot, &
      mask_lake, z_lake)

    implicit none

    real(wp), intent(in) :: z_topo(:,:)
    real(wp), intent(in) :: area(:,:)
    integer, intent(in) :: mask_lake_pot(:,:)
    type(lake_type), intent(inout) :: lake(:)
    integer, intent(out) :: mask_lake(:,:)
    real(wp), intent(out) :: z_lake(:,:)

    integer :: n, k 
    integer :: n_lakes


    n_lakes = size(lake)

    mask_lake = 0
    z_lake = 0._wp

    !$omp parallel do private(n,k)
    do n=1,n_lakes
      ! volume considered for area is a weighted mean of old and new, to avoid lake on-off
      k = 0
      do while (lake(n)%vol_lev(k)<lake(n)%vol)
        k = k+1
        if (k.eq.n_lev) exit
      enddo
      if (k>0) then
        ! consider lake only if area larger than critical area
        !if (lake(n)%area_lev(k)>lake_area_crit) then
          where (mask_lake_pot==n .and. z_topo<lake(n)%z_lev(k))
            mask_lake = n    ! lake n
          endwhere          
          ! lake area
          lake(n)%area = lake(n)%area_lev(k) 
          ! lake surface elevation
          lake(n)%z = lake(n)%z_lev(k) 
          where (mask_lake.eq.n) z_lake = lake(n)%z
        !else
        !  lake(n)%area = 0._wp
        !endif
        ! compute average lake depth
        lake(n)%depth = max(lake_depth_min,lake(n)%vol/lake(n)%area)
      else
        ! no lake, but keep at least deepest potential lake point in the mask
        !mask_lake(lake(n)%hires%i,lake(n)%hires%j) = n
        lake(n)%area = lake(n)%area_min
        lake(n)%z = lake(n)%z_lev(0)
        lake(n)%depth = lake_depth_min
      endif

    enddo
    !omp end parallel do

!    lake%area = 0._wp
!    !$omp parallel do private(i,j,n,k)
!    do i=1,ni
!      do j=1,nj
!        n = mask_lake_pot(i,j)
!        k = 0
!        do while (lake(n)%vol_lev(k)<lake(n)%vol)
!          k = k+1
!          if (k.eq.n_lev) exit
!        enddo
!        if (k>0) then
!          if (z_topo(i,j)<lake(n)%z_lev(k)) then
!            mask_lake = n
!            lake(n)%area = lake(n)%area + area(i,j)
!          endif
!        endif
!      enddo
!    enddo
!    !omp end parallel do
!    print *,'arealake',lake(25)%area

    return

  end subroutine update_lake_mask


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  o c e a n _ l a k e _ m a s k _ l e v e l
  !   Purpose    :  update ocean/lake mask and sea/lake level
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocean_lake_mask_level(mask, z_lake, &
    mask_ocn_lake, z_sea_lake)

    implicit none

    integer, intent(in) :: mask(:,:)
    real(wp), intent(in) :: z_lake(:,:)
    integer, intent(out) :: mask_ocn_lake(:,:)
    real(wp), intent(out) :: z_sea_lake(:,:)

    integer :: i, j, im1, ip1, jm1, jp1, ni, nj, n_ocn, n_lake
    logical :: flag
    integer, allocatable :: mask_tmp(:,:)
    integer, parameter :: i_ocn = 0
    integer, parameter :: i_lake = 1

    where (mask.eq.2)
      ! ocean
      mask_ocn_lake = i_ocn
      z_sea_lake = 0._wp
    elsewhere (mask.eq.4)
      ! lake
      mask_ocn_lake = i_lake
      z_sea_lake = z_lake
    elsewhere
      ! to fill
      mask_ocn_lake = -1
      z_sea_lake = 0._wp
    endwhere

    ! get grid size
    ni = size(mask,1)
    nj = size(mask,2)

    ! allocate
    allocate(mask_tmp(ni,nj))
    flag = .true.
    mask_tmp = mask_ocn_lake
    do while (flag)
      flag = .false.
      do j = 1,nj
        do i = 1,ni
          if (mask_ocn_lake(i,j)==-1) then
            ip1 = i+1; if (ip1.eq.ni+1) ip1 = 1
            im1 = i-1; if (im1.eq.0) im1 = ni
            jp1 = min(j+1,nj)
            jm1 = max(j-1,1)
            n_ocn = 0
            if (mask_ocn_lake(im1,j)==i_ocn) then
              n_ocn = n_ocn+1
            endif
            if (mask_ocn_lake(ip1,j)==i_ocn) then
              n_ocn = n_ocn+1
            endif
            if (mask_ocn_lake(i,jm1)==i_ocn) then
              n_ocn = n_ocn+1
            endif
            if (mask_ocn_lake(i,jp1)==i_ocn) then
              n_ocn = n_ocn+1
            endif
            n_lake = 0
            if (mask_ocn_lake(im1,j)==i_lake) then
              n_lake = n_lake+1
            endif
            if (mask_ocn_lake(ip1,j)==i_lake) then
              n_lake = n_lake+1
            endif
            if (mask_ocn_lake(i,jm1)==i_lake) then
              n_lake = n_lake+1
            endif
            if (mask_ocn_lake(i,jp1)==i_lake) then
              n_lake = n_lake+1
            endif
            if (n_ocn>0 .and. n_ocn>=n_lake) then
              ! update mask 
              mask_tmp(i,j) = i_ocn
            else if (n_lake>0) then
              ! update mask 
              mask_tmp(i,j) = i_lake
            endif
            flag = .true.
          endif
        enddo
      enddo
      mask_ocn_lake = mask_tmp
    enddo
    deallocate(mask_tmp)

    return

  end subroutine ocean_lake_mask_level

end module lake_mod


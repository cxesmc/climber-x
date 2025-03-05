!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : h i r e s _ t o _ l o w r e s _ m o d
!
!  Purpose : aggregate to low coupler resolution 
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
module hires_to_lowres_mod

  use precision, only : wp
  use control, only : flag_lakes
  use constants, only : pi
  use climber_grid, only : lat
  use geo_params, only : f_crit, f_crit_eq
  use geo_params, only : l_ocn_below_shelf
  use geo_params, only : l_close_panama, l_close_bering
  use geo_params, only : i_z_min_max, z_ocn_max_quant
  use geo_grid, only : ni, nj, n_topo_sur, i_topo_sur, j_topo_sur, i0_topo, i1_topo, j0_topo, j1_topo

  implicit none

  private
  public :: hires_to_lowres

contains

  subroutine hires_to_lowres(n_lakes, hires_mask, hires_mask_lake, hires_z_topo, hires_z_topo_fil, hires_z_bed, hires_z_sur, &
    hires_map_runoff, hires_i_runoff, hires_j_runoff, &  ! in
    f_ocn, f_ocn2, f_lnd, f_ice, f_ice_grd, f_ice_flt, f_lake, f_lake_n, &   ! out
    z_sur, z_ocn, z_ocn_min, z_ocn_max, z_ocn_max_q, z_ice, z_lake, z_veg, z_veg_min, z_veg_max, z_bed, z_sur_std, z_sur_smooth_std, z_veg_std, z_sur_lnd_std, &    ! out
    f_drain_veg, f_drain_ice, i_runoff, j_runoff, i_runoff_veg, j_runoff_veg, i_runoff_ice, j_runoff_ice)   ! out

  implicit none

  integer, intent(in) :: n_lakes
  integer, intent(in) :: hires_mask(:,:)
  integer, intent(in) :: hires_mask_lake(:,:)
  real(wp), intent(in) :: hires_z_topo(:,:)
  real(wp), intent(in) :: hires_z_topo_fil(:,:)
  real(wp), intent(in) :: hires_z_bed(:,:)
  real(wp), intent(in) :: hires_z_sur(:,:)
  integer, intent(in) :: hires_map_runoff(:,:)
  integer, intent(in) :: hires_i_runoff(:,:)
  integer, intent(in) :: hires_j_runoff(:,:)

  real(wp), intent(out) :: f_ocn(:,:)
  real(wp), intent(out) :: f_ocn2(:,:)
  real(wp), intent(out) :: f_lnd(:,:)
  real(wp), intent(out) :: f_ice(:,:)
  real(wp), intent(out) :: f_ice_grd(:,:)
  real(wp), intent(out) :: f_ice_flt(:,:)
  real(wp), intent(out) :: f_lake(:,:)
  real(wp), intent(out) :: f_lake_n(:,:,:)
  real(wp), intent(out) :: z_sur(:,:)
  real(wp), intent(out) :: z_ocn(:,:)
  real(wp), intent(out) :: z_ocn_min(:,:)
  real(wp), intent(out) :: z_ocn_max(:,:)
  real(wp), intent(out) :: z_ocn_max_q(:,:)
  real(wp), intent(out) :: z_ice(:,:)
  real(wp), intent(out) :: z_lake(:,:)
  real(wp), intent(out) :: z_veg(:,:)
  real(wp), intent(out) :: z_veg_min(:,:)
  real(wp), intent(out) :: z_veg_max(:,:)
  real(wp), intent(out) :: z_bed(:,:)
  real(wp), intent(out) :: z_sur_std(:,:)
  real(wp), intent(out) :: z_sur_smooth_std(:,:)
  real(wp), intent(out) :: z_veg_std(:,:)
  real(wp), intent(out) :: z_sur_lnd_std(:,:)
  real(wp), intent(out) :: f_drain_veg(0:,:,:)
  real(wp), intent(out) :: f_drain_ice(0:,:,:)
  integer, intent(out) :: i_runoff(:,:)
  integer, intent(out) :: j_runoff(:,:)
  integer, intent(out) :: i_runoff_veg(:,:)
  integer, intent(out) :: j_runoff_veg(:,:)
  integer, intent(out) :: i_runoff_ice(:,:)
  integer, intent(out) :: j_runoff_ice(:,:)

  integer :: i, j, ii, jj, n, nocn, nice, nice_grd, nice_flt, nveg, nlake, nlaken, ncells, nrun_veg, nrun_ice, dir_loc, dir_loc_save, nrun_tot
  real(wp) :: fcrit
  real(wp) :: z_sur_lnd, z_bed_lnd

  logical, dimension(:), allocatable :: mask_ndir
  integer, dimension(:), allocatable :: i_runoff_cell
  integer, dimension(:), allocatable :: j_runoff_cell
  integer, dimension(:), allocatable :: dir_count
  integer, dimension(:), allocatable :: mask_cell
  real(wp), dimension(:), allocatable :: z_topo_cell
  real(wp), dimension(:), allocatable :: z_topo_fil_cell
  real(wp), dimension(:), allocatable :: z_sur_cell
  real(wp), dimension(:), allocatable :: z_bed_cell
  real(wp), dimension(:), allocatable :: z_ocn_tmp
  real(wp), dimension(:), allocatable :: z_ocn_tmp2
  integer, dimension(:), allocatable :: map_runoff_cell
  logical, dimension(:), allocatable :: mask_cell_roff


  allocate( mask_cell(n_topo_sur) )
  allocate( z_topo_cell(n_topo_sur) )
  allocate( z_topo_fil_cell(n_topo_sur) )
  allocate( z_sur_cell(n_topo_sur) )
  allocate( z_bed_cell(n_topo_sur) )
  allocate( map_runoff_cell(n_topo_sur) )
  allocate( mask_cell_roff(n_topo_sur) )
  allocate( i_runoff_cell(n_topo_sur))
  allocate( j_runoff_cell(n_topo_sur))
  allocate( mask_ndir(n_topo_sur))
  allocate( dir_count(n_topo_sur))

  allocate( z_ocn_tmp(n_topo_sur) )
  allocate( z_ocn_tmp2(n_topo_sur) )
  z_ocn_tmp2 = 0._wp

  !$omp parallel do private(i,j,ii,jj,n,fcrit,nocn,nice,nice_grd,nice_flt,nveg,nlake,nlaken,ncells,nrun_veg,nrun_ice) &
  !$omp private(mask_ndir,i_runoff_cell,j_runoff_cell,dir_count,dir_loc,dir_loc_save,nrun_tot) &
  !$omp private(mask_cell,z_topo_cell,z_topo_fil_cell,z_sur_cell,z_bed_cell,z_ocn_tmp,map_runoff_cell,mask_cell_roff,z_sur_lnd,z_bed_lnd)
  do j=1,nj
    do i=1,ni

      mask_cell = pack(hires_mask(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
      z_topo_cell = pack(hires_z_topo(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
      z_topo_fil_cell = pack(hires_z_topo_fil(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
      z_sur_cell = pack(hires_z_sur(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
      z_bed_cell = pack(hires_z_bed(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)
      map_runoff_cell = pack(hires_map_runoff(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)),.true.)

      ! ocean fraction
      if (l_ocn_below_shelf) then
        ! count number of ocean points and floating ice in grid cell
        nocn = count(mask_cell.eq.2 .or. mask_cell.eq.3) 
      else
        ! count number of ocean points in grid cell
        nocn = count(mask_cell.eq.2) 
      endif
      f_ocn(i,j) = real(nocn,wp)/real(n_topo_sur,wp)

      ! apply critical threshold for ocean/land
      if (j==1 .or. j.eq.nj) then
        fcrit = f_crit*2._wp  ! fixme?!
      else
        fcrit = f_crit + max(0._wp,(f_crit_eq-f_crit))*cos(lat(j)*pi/180._wp)**30
      endif
      if (f_ocn(i,j).lt.fcrit) then
        f_ocn(i,j) = 0._wp
      endif

      ! fix for Panama
      if (l_close_panama) then
        ! close Panama, ocean fraction=0
        f_ocn(16,23) = 0._wp
        f_ocn(17:18,22) = 0._wp
        f_ocn(19:20,21) = 0._wp
        f_ocn(21,20) = 0._wp
      endif
      ! option to close Bering 
      if (l_close_bering) then
        ! close Bering, ocean fraction=0
        f_ocn(1:4,32) = 0._wp
      endif

      ! land fraction
      f_lnd(i,j) = 1._wp - f_ocn(i,j)

      ! ocean fraction excluding floating ice
      nocn = count(mask_cell.eq.2) 
      f_ocn2(i,j) = real(nocn,wp)/real(n_topo_sur,wp)
      if (f_ocn2(i,j).gt.f_ocn(i,j)) f_ocn2(i,j)=f_ocn(i,j)     ! make sure it is not larger than total ocean fraction

      ! ice fraction
      ! count number of grounded ice sheet points in grid cell
      nice_grd = count(mask_cell.eq.0)
      f_ice_grd(i,j) = real(nice_grd,wp)/real(n_topo_sur,wp)
      ! count number of floating ice sheet points in grid cell
      nice_flt = count(mask_cell.eq.3)
      f_ice_flt(i,j) = real(nice_flt,wp)/real(n_topo_sur,wp)
      ! count number of ice sheet points (grounded and floating) in grid cell
      nice = nice_grd + nice_flt 
      if (f_ice_flt(i,j).gt.f_ocn(i,j)) then
        ! if floating ice fraction larger than ocean fraction add floating ice excess to grounded ice
        f_ice_grd(i,j) = f_ice_grd(i,j) + (f_ice_flt(i,j)-f_ocn(i,j))
        f_ice_flt(i,j) = f_ocn(i,j) 
      endif
      if (.not.l_ocn_below_shelf) then
        ! all ice is assumed grounded
        f_ice_grd(i,j) = f_ice_grd(i,j) + f_ice_flt(i,j)
        f_ice_flt(i,j) = 0._wp
      endif
      if (f_ocn2(i,j).gt.0._wp .and. f_ocn2(i,j).lt.fcrit) then
        ! ocean fraction below critical value, remove ocean and expand floating ice instead
        f_ice_flt(i,j) = f_ice_flt(i,j) + f_ocn2(i,j)
        f_ocn2(i,j) = 0._wp
      endif
      ! total ice fraction in grid cell
      f_ice(i,j) = f_ice_grd(i,j) + f_ice_flt(i,j)

      ! count number of ice-free and lake-free land points in grid cell
      nveg = count(mask_cell.eq.1) 
      if (f_lnd(i,j).eq.0._wp) nveg = 0

      ! count number of lake points in grid cell
      nlake = count(mask_cell.eq.4)

      if (flag_lakes) then

        ! lakes fraction
        do n=1,n_lakes
          ! count number of lake n points in grid cell
          nlaken = count(hires_mask_lake(i0_topo(i):i1_topo(i),j0_topo(j):j1_topo(j)).eq.n) 
          f_lake_n(n,i,j) = real(nlaken,wp)/real(n_topo_sur,wp)
        enddo
        ! total lake fraction
        f_lake(i,j) = sum(f_lake_n(1:n_lakes,i,j))
        ! lake fraction can not exceed land fraction
        if (f_lake(i,j).gt.f_lnd(i,j)) f_lake(i,j) = f_lnd(i,j)

        ! drainage fractions to lake n or ocean (0)
        ! count number of points draining into lake n
        if (nveg.gt.0._wp) then
          do n=0,n_lakes    ! index 0 is for direct drainage to ocean
            nrun_veg = count(map_runoff_cell.eq.n .and. mask_cell.eq.1) ! land points
            f_drain_veg(n,i,j) = real(nrun_veg,wp)/real(nveg)
          enddo
        else
          f_drain_veg(:,i,j) = 0._wp
        endif
        if (nice.gt.0._wp) then
          do n=0,n_lakes
            nrun_ice = count(map_runoff_cell.eq.n .and. (mask_cell.eq.0 .or. mask_cell.eq.3)) ! ice points
            f_drain_ice(n,i,j) = real(nrun_ice,wp)/real(nice)
          enddo
        else
          f_drain_ice(:,i,j) = 0._wp
        endif

      else ! flag_lakes==F

        f_lake(i,j) = real(nlake,wp)/real(n_topo_sur,wp)

        if (nveg.gt.0._wp) then
          f_drain_veg(0,i,j) = 1._wp ! all runoff directly routed to ocean
        else
          f_drain_veg(0,i,j) = 0._wp 
        endif
        if (nice.gt.0._wp) then
          f_drain_ice(0,i,j) = 1._wp 
        else
          f_drain_ice(0,i,j) = 0._wp 
        endif

      endif

      ! drainage fractions to most frequent ocean destination cell
      ! with filled lakes
      mask_cell_roff = mask_cell.ne.2   ! not ocean
      ncells = count(mask_cell_roff.eqv..true.)
      !print *,'ncells',i,j,ncells
      if (ncells.gt.0) then
        n = 0
        i_runoff_cell = 0
        j_runoff_cell = 0
        do jj=j0_topo(j),j1_topo(j)
          do ii=i0_topo(i),i1_topo(i)
            if (hires_mask(ii,jj).ne.2) then
              n = n+1
              i_runoff_cell(n) = hires_i_runoff(ii,jj)
              j_runoff_cell(n) = hires_j_runoff(ii,jj)
            endif
          enddo
        enddo
        ! mask out duplicate entries
        mask_ndir(1:ncells) = .true.
        mask_ndir(ncells+1:n_topo_sur) = .false.
        do n = ncells,2,-1
          mask_ndir(n) = .not.(any(i_runoff_cell(:n-1)==i_runoff_cell(n).and.j_runoff_cell(:n-1)==j_runoff_cell(n)))
        enddo
        ! count frequency of unique entries
        dir_count = 0
        nrun_tot = 0
        do n=1,ncells
          if (mask_ndir(n).eqv..true.) then
            dir_count(n) = count(i_runoff_cell==i_runoff_cell(n) .and. j_runoff_cell==j_runoff_cell(n))
            nrun_tot = nrun_tot + dir_count(n)
          endif
        enddo
        ! find most frequent entry
        dir_loc = maxloc(dir_count,1)
        ! set most frequent runoff destination point
        i_runoff(i,j) = i_runoff_cell(dir_loc)
        j_runoff(i,j) = j_runoff_cell(dir_loc)

        if (i_runoff(i,j).eq.0) then
          print *,'ERROR i_runoff==0'
          stop
        endif

      else

        i_runoff(i,j) = 0 
        j_runoff(i,j) = 0

      endif

      ! drainage fractions to most frequent ocean destination cell, veg fraction 
      ! only for points not draining to lakes
      mask_cell_roff = mask_cell.eq.1 .and. map_runoff_cell.eq.0
      ncells = count(mask_cell_roff.eqv..true.)
      !print *,'ncells',i,j,ncells
      if (ncells.gt.0) then
        n = 0
        i_runoff_cell = 0
        j_runoff_cell = 0
        do jj=j0_topo(j),j1_topo(j)
          do ii=i0_topo(i),i1_topo(i)
            if (hires_mask(ii,jj).eq.1 .and. hires_map_runoff(ii,jj).eq.0) then
              n = n+1
              i_runoff_cell(n) = hires_i_runoff(ii,jj)
              j_runoff_cell(n) = hires_j_runoff(ii,jj)
            endif
          enddo
        enddo
        ! mask out duplicate entries
        mask_ndir(1:ncells) = .true.
        mask_ndir(ncells+1:n_topo_sur) = .false.
        do n = ncells,2,-1
          mask_ndir(n) = .not.(any(i_runoff_cell(:n-1)==i_runoff_cell(n).and.j_runoff_cell(:n-1)==j_runoff_cell(n)))
        enddo
        ! count frequency of unique entries
        dir_count = 0
        do n=1,ncells
          if (mask_ndir(n).eqv..true.) then
            dir_count(n) = count(i_runoff_cell==i_runoff_cell(n) .and. j_runoff_cell==j_runoff_cell(n))
          endif
        enddo
        ! find most frequent entry
        dir_loc = maxloc(dir_count,1)
        ! set most frequent runoff destination points 
        i_runoff_veg(i,j) = i_runoff_cell(dir_loc)
        j_runoff_veg(i,j) = j_runoff_cell(dir_loc)

        !print *,''
        !print *,'i,j',i,j
        !print *,'ncells',ncells
        !print *,'f_ocn,f_lnd,f_ice,f_lake',f_ocn(i,j),1._wp-f_ocn(i,j),f_ice(i,j),f_lake(i,j)
        !print *,'i_runoff_cell,j_runoff_cell',i_runoff_cell(1:ncells),j_runoff_cell(1:ncells)
        !print *,'mask_ndir',mask_ndir(1:ncells)
        !print *,'dir_count',dir_count(1:ncells)
        !print *,'frac',dir_count(dir_loc)/real(ncells,wp)
        !print *,'dir_loc',dir_loc
        !print *,'i_runoff_veg',i_runoff_veg(i,j)
        !print *,'j_runoff_veg',j_runoff_veg(i,j)
        !print *,mask_cell

        if (i_runoff_veg(i,j).eq.0) then
          print *,'ERROR i_runoff_veg==0'
          stop
        endif

      else

        i_runoff_veg(i,j) = 0 
        j_runoff_veg(i,j) = 0

      endif

      ! drainage fractions to most frequent ocean destination cell, ice fraction 
      mask_cell_roff = (mask_cell.eq.0 .or. mask_cell.eq.3) .and. map_runoff_cell.eq.0
      ncells = count(mask_cell_roff.eqv..true.)
      if (ncells.gt.0) then
        n = 0
        i_runoff_cell = 0
        j_runoff_cell = 0
        do jj=j0_topo(j),j1_topo(j)
          do ii=i0_topo(i),i1_topo(i)
            if ((hires_mask(ii,jj).eq.0 .or. hires_mask(ii,jj).eq.3) .and. hires_map_runoff(ii,jj).eq.0) then
              n = n+1
              i_runoff_cell(n) = hires_i_runoff(ii,jj)
              j_runoff_cell(n) = hires_j_runoff(ii,jj)
            endif
          enddo
        enddo
        ! mask out duplicate entries
        mask_ndir(1:ncells) = .true.
        mask_ndir(ncells+1:n_topo_sur) = .false.
        do n = ncells,2,-1
          mask_ndir(n) = .not.(any(i_runoff_cell(:n-1)==i_runoff_cell(n).and.j_runoff_cell(:n-1)==j_runoff_cell(n)))
        enddo
        ! count frequency of unique entries
        dir_count = 0
        do n=1,ncells
          if (mask_ndir(n)) then
            dir_count(n) = count(i_runoff_cell==i_runoff_cell(n) .and. j_runoff_cell==j_runoff_cell(n))
          endif
        enddo
        ! find most frequent entries
        dir_loc = maxloc(dir_count,1)
        ! set most frequent runoff destination points and relative fraction
        i_runoff_ice(i,j) = i_runoff_cell(dir_loc)
        j_runoff_ice(i,j) = j_runoff_cell(dir_loc)

        !print *,''
        !print *,'i,j',i,j
        !print *,'f_ocn,f_veg,f_ice,f_lake',f_ocn(i,j),1._wp-f_ocn(i,j),f_ice(i,j),f_lake(i,j)
        !print *,'i_runoff_cell,j_runoff_cell',i_runoff_cell(1:ncells),j_runoff_cell(1:ncells)
        !print *,'mask_ndir',mask_ndir(1:ncells)
        !print *,'dir_count',dir_count(1:ncells)
        !print *,'frac',dir_count(dir_loc)/real(ncells,wp)
        !print *,'dir_loc',dir_loc
        !print *,'i_runoff_veg',i_runoff_ice(i,j)
        !print *,'j_runoff_veg',j_runoff_ice(i,j)

        if (i_runoff_ice(i,j).eq.0) then
          print *,'ERROR i_runoff_ice==0'
          stop
        endif
      else

        i_runoff_ice(i,j) = 0 
        j_runoff_ice(i,j) = 0

      endif

      ! grid-cell mean surface elevation
      z_sur(i,j) = sum(z_sur_cell)/real(n_topo_sur,wp)
      if (i_z_min_max.eq.1) then
        ! minimum and maximum elevation in grid cell
        z_veg_min(i,j) = minval(z_sur_cell)
        z_veg_max(i,j) = maxval(z_sur_cell)
      endif

      ! standard deviation of grid cell surface elevation 
      if (f_lnd(i,j).gt.0._wp) then
        z_sur_std(i,j) = sqrt(sum((max(0._wp,z_topo_cell)-z_sur(i,j))**2) / real(n_topo_sur-1,wp))
        z_sur_smooth_std(i,j) = sqrt(sum((max(0._wp,z_topo_cell)-max(0._wp,z_topo_fil_cell))**2) / real(n_topo_sur-1,wp))
      else 
        z_sur_std(i,j) = 0._wp
        z_sur_smooth_std(i,j) = 0._wp
      endif
      
      ! grid-cell mean ocean elevation 
      if (nocn.gt.0) then
        z_ocn(i,j) = sum(z_bed_cell, mask_cell.eq.2)/real(nocn,wp)
        z_ocn_min(i,j) = minval(z_bed_cell, mask_cell.eq.2)
        z_ocn_max(i,j) = maxval(z_bed_cell, mask_cell.eq.2)
        z_ocn_tmp = pack(z_bed_cell, mask_cell.eq.2,z_ocn_tmp2)
        call bubble_sort(z_ocn_tmp)
        z_ocn_max_q(i,j) = z_ocn_tmp(max(1,nint(real(nocn,wp)*z_ocn_max_quant/100._wp)))
        !print *,'nocn,z,z_min,z_max,z_q10',nocn,z_ocn(i,j),z_ocn_min(i,j),maxval(z_bed_cell, mask_cell.eq.2),z_ocn_max_q(i,j)
      else
        z_ocn(i,j) = 0._wp
        z_ocn_min(i,j) = 0._wp 
        z_ocn_max(i,j) = 0._wp
        z_ocn_max_q(i,j) = 0._wp
      endif

      ! grid-cell mean ice-free and land-free land elevation 
      if (nveg.gt.0) then
        z_veg(i,j) = sum(z_sur_cell, mask_cell.eq.1)/real(nveg,wp)
        if (i_z_min_max.eq.2) then
          z_veg_min(i,j) = minval(z_sur_cell, mask_cell.eq.1)
          z_veg_max(i,j) = maxval(z_sur_cell, mask_cell.eq.1)
        endif
      else
        z_veg(i,j) = 0._wp
        if (i_z_min_max.eq.2) then
          z_veg_min(i,j) = 0._wp 
          z_veg_max(i,j) = 0._wp
        endif
      endif

      ! grid-cell mean ice-free standard deviation, including lakes 
      if ((nveg+nlake-1).gt.0) then
        z_veg_std(i,j) = sqrt(sum((max(0._wp,z_topo_cell)-z_veg(i,j))**2, mask_cell.eq.1 .or. mask_cell.eq.4) / real(nveg+nlake-1,wp))
      else
        z_veg_std(i,j) = 0._wp
      endif

      ! grid-cell mean surface ice sheet elevation
      if (nice.gt.0) then 
        z_ice(i,j) = sum(z_sur_cell, mask_cell.eq.0 .or. mask_cell.eq.3)/real(nice,wp)
      else
        z_ice(i,j) = 0._wp
      endif

      ! grid-cell mean surface lakes elevation
      if (nlake.gt.0) then 
        z_lake(i,j) = sum(z_sur_cell, mask_cell.eq.4)/real(nlake,wp)
      else
        z_lake(i,j) = 0._wp
      endif

      ! grid-cell mean bedrock elevation
      z_bed(i,j) = sum(z_bed_cell)/real(n_topo_sur,wp)

      ! grid-cell standard deviation of bedrock elevation over land
      if (nveg+nlake+nice_grd.gt.1) then
        z_sur_lnd = sum(z_topo_cell, mask_cell.eq.0 .or. mask_cell.eq.1 .or. mask_cell.ne.4)/real(nveg+nlake+nice_grd,wp)
        z_sur_lnd_std(i,j) = sqrt(sum((max(0._wp,z_topo_cell)-max(0._wp,z_sur_lnd))**2, mask_cell.eq.0 .or. mask_cell.eq.1 .or. mask_cell.ne.4) &
          / real(nveg+nlake+nice_grd-1,wp))
        z_bed_lnd = sum(z_bed_cell, mask_cell.eq.0 .or. mask_cell.eq.1 .or. mask_cell.ne.4)/real(nveg+nlake+nice_grd,wp)
      else
        z_sur_lnd_std(i,j) = 0._wp
      endif

    enddo
  enddo
  !$omp end parallel do

  deallocate( mask_cell )
  deallocate( z_topo_cell )
  deallocate( z_topo_fil_cell )
  deallocate( z_bed_cell )
  deallocate( map_runoff_cell )
  deallocate( mask_cell_roff )
  deallocate( i_runoff_cell)
  deallocate( j_runoff_cell)
  deallocate( mask_ndir)
  deallocate( dir_count)


  end subroutine hires_to_lowres

  SUBROUTINE Bubble_Sort(a)
  REAL(wp), INTENT(in out), DIMENSION(:) :: a
  REAL(wp) :: temp
  INTEGER :: i, j
  LOGICAL :: swapped

  DO j = SIZE(a)-1, 1, -1
    swapped = .FALSE.
    DO i = 1, j
      IF (a(i) > a(i+1)) THEN
        temp = a(i)
        a(i) = a(i+1)
        a(i+1) = temp
        swapped = .TRUE.
      END IF
    END DO
    IF (.NOT. swapped) EXIT
  END DO
END SUBROUTINE Bubble_Sort

end module hires_to_lowres_mod

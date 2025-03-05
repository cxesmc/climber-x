!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : b m b _ m o d e l
!
!  Purpose : main BMB model
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
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
module bmb_model

  use nml
  use ncio
  use precision, only : wp, dp
  use constants, only : rho_i, rho_sw, Lf, cap_w
  use timer, only : nmon_year, nstep_year_bmb, time_soy_bmb, time_eoy_bmb
  use control, only : out_dir, restart_in_dir, bmb_restart
  use bmb_grid, only : bmb_grid_init
  use bmb_params, only : bmb_params_init, dt, i_bmb, bmb_const, k_1, k_2, l_fix_depth, fix_depth, depth_disc, l_depth_scale, zl_ref, i_bmb_lake
  use bmb_def, only : bmb_class
  use coord, only : grid_allocate, grid_class
  use coord, only : map_scrip_init, map_scrip_class, map_scrip_field

  implicit none


  private
  public :: bmb_init, bmb_update, bmb_end, bmb_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b m b _ u p d a t e
  !   Purpose    :  update basal mass balance
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_update(bmb)

  !$  use omp_lib

  implicit none

  type(bmb_class) :: bmb

  integer :: i, j, im1, ip1, jm1, jp1, n, k, k1, k2
  integer :: nk_ocn, k_fix_ocn, k_disc_ocn
  integer :: nk_lake, k_fix_lake, k_disc_lake
  logical :: flag
  real(wp) :: t_sum, s_sum
  real(wp) :: t_npol, s_npol, t_spol, s_spol
  real(wp) :: w1, w2


  !---------------------------
  ! ocean

  nk_ocn = size(bmb%z_ocn_in)

  ! extrapolate ocean temperature and salinity to cells where they are undefined
  do k = 1,nk_ocn
    flag = .true.
    bmb%mask_ocn_in_tmp = bmb%mask_ocn_in(:,:,k)
    do while (flag)
      flag = .false.
      do j = 1,bmb%grid_in%G%ny
        do i = 1,bmb%grid_in%G%nx
          if (bmb%mask_ocn_in(i,j,k)==0) then
            n = 0
            t_sum = 0._wp
            s_sum = 0._wp
            ip1 = i+1; if (ip1.eq.bmb%grid_in%G%nx+1) ip1 = 1
            im1 = i-1; if (im1.eq.0) im1 = bmb%grid_in%G%nx
            jp1 = min(j+1,bmb%grid_in%G%ny)
            jm1 = max(j-1,1)
            if (bmb%mask_ocn_in(im1,j,k)==1) then
              n = n+1
              t_sum = t_sum + bmb%t_ocn_in(im1,j,k)
              s_sum = s_sum + bmb%s_ocn_in(im1,j,k)
            endif
            if (bmb%mask_ocn_in(ip1,j,k)==1) then
              n = n+1
              t_sum = t_sum + bmb%t_ocn_in(ip1,j,k)
              s_sum = s_sum + bmb%s_ocn_in(ip1,j,k)
            endif
            if (bmb%mask_ocn_in(i,jm1,k)==1) then
              n = n+1
              t_sum = t_sum + bmb%t_ocn_in(i,jm1,k)
              s_sum = s_sum + bmb%s_ocn_in(i,jm1,k)
            endif
            if (bmb%mask_ocn_in(i,jp1,k)==1) then
              n = n+1
              t_sum = t_sum + bmb%t_ocn_in(i,jp1,k)
              s_sum = s_sum + bmb%s_ocn_in(i,jp1,k)
            endif
            if (n>0) then
              flag = .true.
              bmb%t_ocn_in(i,j,k) = t_sum/real(n,wp)
              bmb%s_ocn_in(i,j,k) = s_sum/real(n,wp)
              ! update ocean mask 
              bmb%mask_ocn_in_tmp(i,j) = 1
            endif
          endif
        enddo
      enddo
      bmb%mask_ocn_in(:,:,k) = bmb%mask_ocn_in_tmp
    enddo
  enddo

  ! polar smoothing
  do k = 1,nk_ocn
    t_npol = 0._wp
    s_npol = 0._wp
    t_spol = 0._wp
    s_spol = 0._wp
    do i = 1,bmb%grid_in%G%nx
      t_npol = t_npol + bmb%t_ocn_in(i,bmb%grid_in%G%ny,k) / bmb%grid_in%G%nx
      s_npol = s_npol + bmb%s_ocn_in(i,bmb%grid_in%G%ny,k) / bmb%grid_in%G%nx
      t_spol = t_spol + bmb%t_ocn_in(i,1,k) / bmb%grid_in%G%nx
      s_spol = s_spol + bmb%s_ocn_in(i,1,k) / bmb%grid_in%G%nx
    enddo
    bmb%t_ocn_in(:,1,k) = t_spol
    bmb%s_ocn_in(:,1,k) = s_spol
    bmb%t_ocn_in(:,bmb%grid_in%G%ny,k) = t_npol
    bmb%s_ocn_in(:,bmb%grid_in%G%ny,k) = s_npol
  enddo

  if (.not.allocated(bmb%t_ocn)) allocate(bmb%t_ocn(bmb%grid%G%nx,bmb%grid%G%ny,nk_ocn))
  if (.not.allocated(bmb%s_ocn)) allocate(bmb%s_ocn(bmb%grid%G%nx,bmb%grid%G%ny,nk_ocn))

  ! map shelf temperature and salinity to ice sheet grid
  !$omp parallel do private(k)
  do k=1,nk_ocn
    call map_scrip_field(bmb%maps_cmn_to_ice,"t_ocn",bmb%t_ocn_in(:,:,k), bmb%t_ocn(:,:,k),method="mean",missing_value=-9999._dp)
    call map_scrip_field(bmb%maps_cmn_to_ice,"s_ocn",bmb%s_ocn_in(:,:,k), bmb%s_ocn(:,:,k),method="mean",missing_value=-9999._dp)
  enddo
  !$omp end parallel do 


  !---------------------------
  ! lake

  if (i_bmb_lake.ne.0) then

    nk_lake = size(bmb%z_lake_in)

    ! extrapolate lake temperature and salinity to cells where they are undefined
    flag = .true.
    bmb%mask_lake_in_tmp = bmb%mask_lake_in
    do while (flag)
      flag = .false.
      do j = 1,bmb%grid_in%G%ny
        do i = 1,bmb%grid_in%G%nx
          if (bmb%mask_lake_in(i,j)==0) then
            ip1 = i+1; if (ip1.eq.bmb%grid_in%G%nx+1) ip1 = 1
            im1 = i-1; if (im1.eq.0) im1 = bmb%grid_in%G%nx
            jp1 = min(j+1,bmb%grid_in%G%ny)
            jm1 = max(j-1,1)
            do k = 1,nk_lake
              n = 0
              t_sum = 0._wp
              s_sum = 0._wp
              if (bmb%mask_lake_in(im1,j)==1) then
                n = n+1
                t_sum = t_sum + bmb%t_lake_in(im1,j,k)
                s_sum = s_sum + bmb%s_lake_in(im1,j,k)
              endif
              if (bmb%mask_lake_in(ip1,j)==1) then
                n = n+1
                t_sum = t_sum + bmb%t_lake_in(ip1,j,k)
                s_sum = s_sum + bmb%s_lake_in(ip1,j,k)
              endif
              if (bmb%mask_lake_in(i,jm1)==1) then
                n = n+1
                t_sum = t_sum + bmb%t_lake_in(i,jm1,k)
                s_sum = s_sum + bmb%s_lake_in(i,jm1,k)
              endif
              if (bmb%mask_lake_in(i,jp1)==1) then
                n = n+1
                t_sum = t_sum + bmb%t_lake_in(i,jp1,k)
                s_sum = s_sum + bmb%s_lake_in(i,jp1,k)
              endif
              if (n>0) then
                bmb%t_lake_in(i,j,k) = t_sum/real(n,wp)
                bmb%s_lake_in(i,j,k) = s_sum/real(n,wp)
                ! update ocean mask 
                bmb%mask_lake_in_tmp(i,j) = 1
              endif
            enddo
            flag = .true.
          endif
        enddo
      enddo
      bmb%mask_lake_in = bmb%mask_lake_in_tmp
    enddo

    if (.not.allocated(bmb%t_lake)) allocate(bmb%t_lake(bmb%grid%G%nx,bmb%grid%G%ny,nk_lake))
    if (.not.allocated(bmb%s_lake)) allocate(bmb%s_lake(bmb%grid%G%nx,bmb%grid%G%ny,nk_lake))

    ! map lake temperature and salinity to ice sheet grid
    !$omp parallel do private(k)
    do k=1,nk_lake
      call map_scrip_field(bmb%maps_cmn_to_ice,"t_lake",bmb%t_lake_in(:,:,k), bmb%t_lake(:,:,k),method="mean",missing_value=-9999._dp)
      call map_scrip_field(bmb%maps_cmn_to_ice,"s_lake",bmb%s_lake_in(:,:,k), bmb%s_lake(:,:,k),method="mean",missing_value=-9999._dp)
    enddo
    !$omp end parallel do 

  endif

  ! index of fix_depth 
  k_fix_ocn = minloc(abs(bmb%z_ocn_in-fix_depth),1)
  k_disc_ocn = minloc(abs(bmb%z_ocn_in-depth_disc),1)
  k_fix_lake = minloc(abs(bmb%z_lake_in-fix_depth),1)
  k_disc_lake = minloc(abs(bmb%z_lake_in-depth_disc),1)

  !$omp parallel do collapse(2) private(i,j,k,k1,k2,w1,w2)
  do i = 1,bmb%grid%G%nx
    do j = 1,bmb%grid%G%ny

      if (time_soy_bmb) then
        bmb%bmb_ann(i,j) = 0._wp
        bmb%t_disc(i,j) = 0._wp
        bmb%s_disc(i,j) = 0._wp
      endif

      if (bmb%mask_ocn_lake(i,j).eq.0 .or. i_bmb_lake.eq.0) then
        ! use ocean values

        ! annual mean values for small-scale basal melt in ice sheet model
        bmb%t_disc(i,j) = bmb%t_disc(i,j) + bmb%t_ocn(i,j,k_disc_ocn)/real(nstep_year_bmb,wp)
        bmb%s_disc(i,j) = bmb%s_disc(i,j) + bmb%s_ocn(i,j,k_disc_ocn)/real(nstep_year_bmb,wp) 

        ! select which level to use
        if (l_fix_depth) then
          ! use temperature and salinity at fixed depth
          bmb%t_bmb(i,j) = bmb%t_ocn(i,j,k_fix_ocn)
          bmb%s_bmb(i,j) = bmb%s_ocn(i,j,k_fix_ocn) 
        else
          ! interpolate temperature and salinity to depth of ice shelf base
          k = 1
          do while (bmb%z_ocn_in(k).lt.-bmb%zb(i,j) .and. k.le.nk_ocn)
            k = k+1
          enddo
          k1 = max(1,k-1)
          k2 = min(nk_ocn,k)
          if (k1.eq.k2) then
            w1 = 1._wp
          else
            w1  = 1._wp - (-bmb%zb(i,j)-bmb%z_ocn_in(k1)) / (bmb%z_ocn_in(k2)-bmb%z_ocn_in(k1))
          endif
          w2  = 1._wp-w1
          !print *
          !print *,'zb',-bmb%zb(i,j)
          !print *,'z_ocn_in',bmb%z_ocn_in
          !print *,'k1,k2',k1,k2
          !print *,'w1,w2',w1,w2
          !print *,'t_ocn',bmb%t_ocn(i,j,:)
          !print *,'s_ocn',bmb%s_ocn(i,j,:)
          bmb%t_bmb(i,j) = w1*bmb%t_ocn(i,j,k1) + w2*bmb%t_ocn(i,j,k2)
          bmb%s_bmb(i,j) = w1*bmb%s_ocn(i,j,k1) + w2*bmb%s_ocn(i,j,k2)
          !print *,'t_bmb',bmb%t_bmb(i,j)
          !print *,'s_bmb',bmb%s_bmb(i,j)
        endif

      else if (bmb%mask_ocn_lake(i,j).eq.1) then
        ! use lake values

        ! annual mean values for small-scale basal melt in ice sheet model
        bmb%t_disc(i,j) = bmb%t_disc(i,j) + bmb%t_lake(i,j,k_disc_lake)/real(nstep_year_bmb,wp)
        bmb%s_disc(i,j) = bmb%s_disc(i,j) + bmb%s_lake(i,j,k_disc_lake)/real(nstep_year_bmb,wp) 

        ! select which level to use
        if (l_fix_depth) then
          ! use temperature and salinity at fixed depth
          bmb%t_bmb(i,j) = bmb%t_lake(i,j,k_fix_lake)
          bmb%s_bmb(i,j) = bmb%s_lake(i,j,k_fix_lake) 
        else
          ! interpolate temperature and salinity to depth of ice shelf base
          k = 1
          do while (bmb%z_lake_in(k).lt.-bmb%zb(i,j) .and. k.le.nk_lake)
            k = k+1
          enddo
          k1 = max(1,k-1)
          k2 = min(nk_lake,k)
          if (k1.eq.k2) then
            w1 = 1._wp
          else
            w1  = 1._wp - (-bmb%zb(i,j)-bmb%z_lake_in(k1)) / (bmb%z_lake_in(k2)-bmb%z_lake_in(k1))
          endif
          w2  = 1._wp-w1
          !print *
          !print *,'zb',-bmb%zb(i,j)
          !print *,'z_lake_in',bmb%z_lake_in
          !print *,'k1,k2',k1,k2
          !print *,'w1,w2',w1,w2
          !print *,'t_lake',bmb%t_lake(i,j,:)
          !print *,'s_lake',bmb%s_lake(i,j,:)
          bmb%t_bmb(i,j) = w1*bmb%t_lake(i,j,k1) + w2*bmb%t_lake(i,j,k2)
          bmb%s_bmb(i,j) = w1*bmb%s_lake(i,j,k1) + w2*bmb%s_lake(i,j,k2)
          !print *,'t_bmb',bmb%t_bmb(i,j)
          !print *,'s_bmb',bmb%s_bmb(i,j)
        endif

      endif

      ! freezing temperature after Beckmann & Goose 2003 (eq. 2), zb<0
      bmb%t_freeze(i,j) = 0.0939_wp - 0.057_wp*bmb%s_bmb(i,j) + 7.64e-4_wp*bmb%zb(i,j) ! degC

      if (i_bmb==0) then
        ! constant and uniform basal mass balance
        bmb%bmb(i,j) = bmb_const  ! kg/m2/s
      else if (i_bmb==1) then
        ! Beckmann & Goose 2003
        bmb%bmb(i,j) = - k_1 * rho_sw*cap_w/(rho_i*Lf) * (bmb%t_bmb(i,j)-bmb%t_freeze(i,j))  ! kg/m2/s
      else if (i_bmb==2) then
        ! Pollard & DeConto 2012, eq 17
        bmb%bmb(i,j) = - k_2 * rho_sw*cap_w/(rho_i*Lf) * abs(bmb%t_bmb(i,j)-bmb%t_freeze(i,j))*(bmb%t_bmb(i,j)-bmb%t_freeze(i,j))  ! kg/m2/s
      endif

      if (bmb%mask_ocn_lake(i,j).eq.1 .and. i_bmb_lake.eq.2) then
        ! zero basal mass balance in lakes
        bmb%bmb(i,j) = 0._wp
      endif

      ! additional dependence on ocean depth, reduce melt in shallow water and enhance in deep water
      if (l_depth_scale) then
        if (bmb%bmb(i,j).lt.0._wp) then ! only for melt case, exclude refreezing
          bmb%bmb(i,j) = bmb%bmb(i,j) * (1._wp+max(0._wp,max(0._wp,-bmb%zl_fil(i,j))-zl_ref)/zl_ref)
        endif
      endif

      ! exclude melt in Hudson Bay
!      where (bmb%grid%lat.gt.50._wp .and. bmb%grid%lat.lt.70._wp .and. bmb%grid%lon.gt.-100._wp .and. bmb%grid%lon.lt.-60._wp)
!        bmb%bmb = 0._wp
!      endwhere

      bmb%bmb_ann(i,j) = bmb%bmb_ann(i,j) + bmb%bmb(i,j)*dt   ! bmb integrated over the year, kg/m2

    enddo
  enddo
  !$omp end parallel do


  return

  end subroutine bmb_update
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b m b _ i n i t
  !   Purpose    :  initialize basal mass balance 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_init(bmb,grid,cmn_grid)

    implicit none

    type(bmb_class) :: bmb
    type(grid_class), intent(in) :: grid 
    type(grid_class), intent(in) :: cmn_grid 


    call bmb_params_init
    call bmb_grid_init

    ! store grid information
    bmb%grid = grid

    ! Generate mapping from cmn to bmb/ice
    call map_scrip_init(bmb%maps_cmn_to_ice,cmn_grid,bmb%grid,method="bil",fldr="maps",load=.TRUE.,clean=.FALSE.)

    bmb%grid_in = cmn_grid

    ! Allocate arrays 
    call grid_allocate(bmb%grid_in, bmb%mask_ocn_in_tmp )
    call grid_allocate(bmb%grid_in, bmb%mask_lake_in_tmp )
    call grid_allocate(bmb%grid, bmb%mask_ocn_lake)
    call grid_allocate(bmb%grid, bmb%mask_ice_shelf )
    call grid_allocate(bmb%grid, bmb%zb )
    call grid_allocate(bmb%grid, bmb%t_bmb )
    call grid_allocate(bmb%grid, bmb%s_bmb )
    call grid_allocate(bmb%grid, bmb%t_disc )
    call grid_allocate(bmb%grid, bmb%s_disc )
    call grid_allocate(bmb%grid, bmb%t_freeze )
    call grid_allocate(bmb%grid, bmb%bmb )
    call grid_allocate(bmb%grid, bmb%bmb_ann )

    if (bmb_restart) then

      ! read restart file 
      call bmb_read_restart(trim(restart_in_dir)//"/bmb_"//trim(bmb%grid%name)//"_restart.nc",bmb)
      print *,'read restart file ',trim(restart_in_dir)//"/bmb_"//trim(bmb%grid%name)//"_restart.nc"

    else

      ! initialize for new run
      bmb%bmb     = 0._wp
      bmb%bmb_ann = 0._wp

    endif

    print*
    print*,'======================================================='
    print*,' Initialisation of bmb ', trim(bmb%grid%name),' complete'
    print*,'======================================================='
    print*

  return

  end subroutine bmb_init
      
     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ e n d 
  ! Purpose  :  end bmb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_end(bmb)

    implicit none

    type(bmb_class) :: bmb


    ! Deallocate all state variables to free memory
    call bmb_dealloc(bmb)


    return

  end subroutine bmb_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ d e a l l o c 
  ! Purpose  :  deallocate bmb variables
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_dealloc(bmb)

    implicit none 

    type(bmb_class) :: bmb


    ! Deallocate all state variables to free memory
    deallocate(bmb%mask_ocn_in)
    deallocate(bmb%mask_ocn_in_tmp)
    deallocate(bmb%t_ocn_in )
    deallocate(bmb%s_ocn_in )
    deallocate(bmb%mask_lake_in)
    deallocate(bmb%mask_lake_in_tmp)
    deallocate(bmb%t_lake_in )
    deallocate(bmb%s_lake_in )
    deallocate(bmb%t_ocn )
    deallocate(bmb%s_ocn )
    deallocate(bmb%t_lake )
    deallocate(bmb%s_lake )
    deallocate(bmb%mask_ocn_lake)
    deallocate(bmb%mask_ice_shelf )
    deallocate(bmb%zb )
    deallocate(bmb%t_bmb )
    deallocate(bmb%s_bmb )
    deallocate(bmb%t_disc )
    deallocate(bmb%s_disc )
    deallocate(bmb%t_freeze )
    deallocate(bmb%bmb )
    deallocate(bmb%bmb_ann )

    return

  end subroutine bmb_dealloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_write_restart(fnm,bmb)

    use dim_name, only: dim_x, dim_y

    implicit none

    type(bmb_class) :: bmb

    character (len=*) :: fnm
    integer :: ncid
    integer :: nx, ny


    nx = bmb%grid%G%nx
    ny = bmb%grid%G%ny

    ! Create the netcdf file and the dimension variables
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_x,x=bmb%grid%G%x0,dx=bmb%grid%G%dx,nx=nx,axis="x",units="m",ncid=ncid)
    call nc_write_dim(fnm,dim_y,x=bmb%grid%G%y0,dx=bmb%grid%G%dy,nx=ny,axis="y",units="m",ncid=ncid)

    call nc_write(fnm,"bmb_ann", bmb%bmb_ann, dims=[dim_x,dim_y],start=[1,1],count=[nx,ny], &
      long_name="floating ice basal mass balance",grid_mapping="polar_stereographic",units="kg/m2/a",ncid=ncid)    

    call nc_close(ncid)


  end subroutine bmb_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ r e a d _ r e s t a r t
  ! Purpose  :  read bmb restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_read_restart(fnm,bmb)

    implicit none

    type(bmb_class), intent(inout) :: bmb

    character (len=*) :: fnm


    call nc_read(fnm,"bmb_ann", bmb%bmb_ann)


    return

  end subroutine bmb_read_restart

end module bmb_model

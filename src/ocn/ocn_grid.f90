!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : o c n _ g r i d
!
!  Purpose : ocean grid definition
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
module ocn_grid

  use precision, only : wp, dp
  use ncio
  use timer, only : year
  use constants, only : pi, r_earth
  use climber_grid, only: ni, nj, lon, lat, basin_mask, i_atlantic
  use control, only: in_dir, out_dir
  use ocn_params, only : shelf_depth, nlayers, i_smooth, smooth_fac, zw_in
  use ocn_params, only : z_mix_brines, dbl
  use ocn_params, only : i_isl, l_isl_ant, l_isl_aus, l_isl_grl, l_isl_ame, n_isl, lat_isl, lon_isl

  implicit none

  integer :: maxi     !! number of gridcells in longitudinal direction []
  integer :: maxj     !! number of gridcells in latitudinal direction []
  integer :: maxk     !! number of gridcells in vertical direction []
  integer, parameter :: maxisles = 10     !! maximum allowed number of islands []
  integer :: n_isles   !! actual number of islands []
  integer :: mpi 
  integer :: mpxi 
  integer :: mpxj 

  ! ocean masks, vary in time
  integer, allocatable :: mask_ocn(:,:)  !! 2D ocean mask []
  integer, allocatable :: mask_c(:,:,:)  !! 3D ocean mask on c-grid (tracer grid) []
  integer, allocatable :: mask_u(:,:,:)  !! 3D ocean mask on u-grid []
  integer, allocatable :: mask_v(:,:,:)  !! 3D ocean mask on v-grid []
  integer, allocatable :: mask_w(:,:,:)  !! 3D ocean mask on w-grid []

  real(wp), dimension(:,:), allocatable :: topo   !! topography [m]
  real(wp), dimension(:,:), allocatable :: bathy   !! bathymetry [m]
  integer, dimension(:,:), allocatable :: k1_pot   !! potential index of the first wet ocean grid cell, counting starts from bottom []
  integer, dimension(:,:), allocatable :: k1       !! index of the first wet ocean grid cell, counting starts from bottom []
  integer :: k1_shelf
  integer :: k1_1000
  integer :: k1_3000
  integer :: k_mix_brines
  integer, dimension(:,:,:), allocatable :: ku       !! index of the first wet ocean grid cell on velocity grid (u,v) []
  logical, dimension(:,:), allocatable :: getj     !!  array to avoid J term in flat regions

  integer, dimension(:,:), allocatable :: map, edge, islands
  integer, dimension(:,:), allocatable :: map_isles
  integer, dimension(:,:), allocatable :: map_isles_ext
  integer, dimension(:,:,:), allocatable :: map_edge
  integer, dimension(:,:,:), allocatable :: map_path
  real(wp), allocatable :: psiles(:)
  integer, dimension(:,:), allocatable :: lpisl, ipisl, jpisl !! islands path integral arrays
  integer, dimension(:), allocatable :: npi    !! islands path integral array

  integer, dimension(20) :: i_idx_isl       !! i-index of islands
  integer, dimension(20) :: j_idx_isl       !! j-index of islands

  real(wp) :: depth  !! maximum ocean depth, updated during runtime! [m]
  real(wp) :: phi0   !! easternmost longitude [radians] 
  real(wp) :: dphi   !! longitudinal resolution [radians] 
  real(wp) :: rdphi  !! reverse of dphi [1/radians] 
  real(wp) :: dtheta !! latitudinal resolution [radians]
  real(wp) :: dy     !! latitudinal grid-cell width [m]
  real(wp) :: rdy    !! reverse of dy [1/m]
  real(wp), dimension(:), allocatable :: s   !! sine of latitude at cell centers []
  real(wp), dimension(:), allocatable :: sv  !! sine of latitude at cell edges []
  real(wp), dimension(:), allocatable :: ds  !! Delta(s) []
  real(wp), dimension(:), allocatable :: dsv !! Delta(sv) []
  real(wp), dimension(:), allocatable :: rds !! reverse of ds
  real(wp), dimension(:), allocatable :: rdsv!! reverse of dsv
  real(wp), dimension(:), allocatable :: rds2!! 2.0/(dsv(j)+dsv(j-1))
  real(wp), dimension(:), allocatable :: c   !! cosine of latitude at cell centers []
  real(wp), dimension(:), allocatable :: cv  !! cosine of latitude at cell edges []
  real(wp), dimension(:), allocatable :: cv2 !! cv(j)*cv(j)*rdsv(j)
  real(wp), dimension(:), allocatable :: rc  !! reverse of c
  real(wp), dimension(:), allocatable :: rc2 !! rc*rc*rdphi
  real(wp), dimension(:), allocatable :: rcv !! reverse of cv
  real(wp), dimension(:), allocatable :: dx  !! longitudinal grid-cell width at cell centers [m] 
  real(wp), dimension(:), allocatable :: rdx !! reverse of dx
  real(wp), dimension(:), allocatable :: dxv !! longitudinal grid-cell width at cell edges [m]
  real(wp), dimension(:), allocatable :: rdxv!! reverse of dxv
  real(wp), dimension(:), allocatable :: zro !! depth of layer center, zero at the surface, decreasing [m]
  real(wp), dimension(:), allocatable :: zw  !! depth of layer edges, zero at the surface, decreasing [m]
  real(wp), dimension(:), allocatable :: dz  !! thickness of layers [m] 
  real(wp), dimension(:), allocatable :: dza !! thickness between layer centers [m]
  real(wp), dimension(:), allocatable :: rdz  !! reverse of dz [1/m] 
  real(wp), dimension(:), allocatable :: rdza !! reverse of dza [1/m]
  real(wp), dimension(:), allocatable :: f_pbl  !! fraction of layers in boundary layer [] 

  real(wp) :: dzz
  real(wp), dimension(:,:), allocatable :: dzg !! 2D "depth" grid of difference between each pair of levels [m]
  real(wp), dimension(:,:), allocatable :: rdzg!! reverse of dzg [1/m]
  real(wp), dimension(:,:), allocatable :: z2dzg!! 2D "depth" grid of difference between squares of levels [m2]
  real(wp), dimension(:,:,:), allocatable :: h !! seabed depth [m]
  real(wp), dimension(:,:,:), allocatable :: rh !! reverse of seabed depth [1/m]

  real(dp) :: ocn_area_tot   !! total surface area of ocean [m2]
  real(dp) :: ocn_vol_tot    !! total ocean volume [m3]
  real(wp), dimension(:,:), allocatable :: ocn_area  !! horizonzal area of ocean cells (ocean fraction is accounted for) [m2]
  real(wp), dimension(:,:,:), allocatable :: ocn_vol  !! volume of ocean cells (ocean fraction is accounted for) [m2]

  type grid_class
    integer :: ni, nj, nk
    real(wp), dimension(:), allocatable :: lat  !! latitude [deg]
    real(wp), dimension(:), allocatable :: lon  !! longitude [deg]
    integer, dimension(:,:), allocatable :: k1  !! ocean first layer (bottom) index []
    integer, dimension(:,:), allocatable :: k1_pot   !! potential index of the first wet ocean grid cell, counting starts from bottom []
    integer, dimension(:,:), allocatable :: mask_ocn  !! ocean mask []
    real(dp) :: ocn_area_tot   !! total surface area of ocean [m2]
    real(wp), dimension(:,:), allocatable :: ocn_area  !! horizonzal area of ocean cells (ocean fraction is accounted for) [m2]
    real(wp), dimension(:,:), allocatable :: ocn_area_old  !! horizonzal area of ocean cells (ocean fraction is accounted for) [m2]
    real(dp) :: ocn_vol_tot    !! total ocean volume [m3]
    real(dp) :: ocn_vol_tot_real    !! real total ocean volume derived from high resolution bathymetry and sea level [m3]
    real(wp), dimension(:,:,:), allocatable :: ocn_vol  !! volume of ocean cells (ocean fraction is accounted for) [m2]
    real(wp), dimension(:,:,:), allocatable :: ocn_vol_old  !! volume of ocean cells (ocean fraction is accounted for) [m2]
    integer, allocatable :: mask_c(:,:,:)  !! 3D ocean mask on c-grid (tracer grid) []
    real(wp), dimension(:), allocatable :: zro !! depth of layer center, zero at the surface, decreasing [m]
    real(wp), dimension(:), allocatable :: zw  !! depth of layer edges, zero at the surface, decreasing [m]
    real(wp), dimension(:), allocatable :: dz  !! thickness of layers [m] 
    real(wp), dimension(:), allocatable :: dza !! thickness between layer centers [m]
    logical :: l_large_vol_change    
  end type

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ g r i d _ i n i t
  ! Purpose  :  initialize ocean grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_grid_init(f_ocn,topo_in,z_ocn_max,mask_coast,ocn_vol_tot_real,grid)

    implicit none
 
    real(wp), intent(in) :: f_ocn(:,:)  !! ocean fracŧion
    real(wp), intent(in) :: topo_in(:,:)  !! 'real' ocean bathymetry [m]
    real(wp), intent(in) :: z_ocn_max(:,:)  !! 'real' ocean bathymetry [m]
    integer, intent(in) :: mask_coast(:,:)  !! 'real' ocean bathymetry [m]
    real(dp), intent(in) :: ocn_vol_tot_real
    type(grid_class), intent(inout) :: grid

    integer :: i, j, k, n
    real(wp) :: th0, th1, s0, s1, theta, thv, dth, deg_to_rad
    real(wp) :: tv1
    real(wp) :: zfac
    real(wp) :: dist_isl


    maxi = ni   
    maxj = nj  
    mpi = 2*(maxi+maxj)
    mpxi = maxi
    mpxj = maxj+1
    ! number of ocean layers
    maxk = nlayers

    ! allocate variables
    allocate(grid%k1(maxi,maxj))
    allocate(grid%k1_pot(maxi,maxj))
    allocate(grid%mask_ocn(maxi,maxj))
    allocate(grid%ocn_area(maxi,maxj))
    allocate(grid%ocn_vol(maxi,maxj,maxk))
    allocate(grid%ocn_area_old(maxi,maxj))
    allocate(grid%ocn_vol_old(maxi,maxj,maxk))
    allocate(grid%mask_c(maxi,maxj,maxk))
    allocate(grid%zro(maxk))
    allocate(grid%zw(0:maxk))
    allocate(grid%dz(maxk))
    allocate(grid%dza(maxk))
    allocate(grid%lat(maxj))
    allocate(grid%lon(maxi))

    allocate(s(0:maxj))
    allocate(sv(0:maxj))
    allocate(ds(maxj))
    allocate(dsv(1:maxj-1))
    allocate(rds(maxj))
    allocate(rds2(2:maxj-1))
    allocate(rdsv(1:maxj-1))
    allocate(c(0:maxj))
    allocate(cv(0:maxj))
    allocate(cv2(1:maxj-1))
    allocate(rc(0:maxj))
    allocate(rc2(0:maxj)) 
    allocate(rcv(1:maxj-1))
    allocate(dx(maxj))
    allocate(rdx(maxj))
    allocate(dxv(0:maxj))
    allocate(rdxv(0:maxj))
    allocate(zro(maxk))
    allocate(zw(0:maxk))
    allocate(dz(maxk))
    allocate(dza(maxk))
    allocate(dzg(maxk,maxk))
    allocate(rdz(maxk))
    allocate(rdza(maxk))
    allocate(rdzg(maxk,maxk))
    allocate(z2dzg(maxk,maxk))
    allocate(f_pbl(maxk))
    allocate(h(3,0:maxi+1,0:maxj+1))
    allocate(rh(3,0:maxi+1,0:maxj+1))
    allocate(ocn_area(maxi,maxj))
    allocate(ocn_vol(maxi,maxj,maxk))
    allocate(topo(maxi,maxj))
    allocate(bathy(maxi,maxj))
    allocate(k1_pot(maxi,maxj))
    allocate(k1(0:maxi+1,0:maxj+1))
    allocate(ku(2,maxi,maxj))
    allocate(getj(maxi,maxj))
    allocate(map(maxi,maxj))
    allocate(edge(maxi,maxj))
    allocate(islands(maxi,maxj))
    allocate(map_edge(maxi,maxj,maxisles))
    allocate(map_path(maxi,maxj,maxisles))
    allocate(map_isles(maxi,maxj))
    allocate(map_isles_ext(0:maxi+1,0:maxj+1))
    allocate(psiles(mpxi*mpxj))
    allocate(lpisl(mpi,maxisles))
    allocate(ipisl(mpi,maxisles))
    allocate(jpisl(mpi,maxisles))
    allocate(npi(maxisles))
    allocate(mask_ocn(maxi,maxj))
    allocate(mask_c(maxi,maxj,maxk))
    allocate(mask_u(maxi,maxj,maxk))
    allocate(mask_v(maxi,maxj,maxk))
    allocate(mask_w(maxi,maxj,maxk))
  
    !-------------------------------------------------------------
    ! set up horizontal grid
    !-------------------------------------------------------------

    th0 = - pi/2._wp
    th1 = pi/2._wp
    s0 = sin(th0)
    s1 = sin(th1)
    deg_to_rad = pi/180._wp

    phi0 = -180.0*deg_to_rad 
    dphi = 2._wp*pi/maxi
    rdphi = 1._wp/dphi

    sv(0) = s0
    cv(0) = cos(th0)
    ! set up const dlat grid
    dth = (th1 - th0)/maxj
    dtheta = dth
    do j=1,maxj
      thv = th0 + j*dth
      theta = thv - 0.5_wp*dth
      sv(j) = sin(thv)
      s(j) = sin(theta)
      cv(j) = cos(thv)
    enddo

    do j=1,maxj
      ds(j) = sv(j) - sv(j-1)
      rds(j) = 1._wp/ds(j)
      c(j) = sqrt(1._wp - s(j)*s(j))
      rc(j) = 1._wp/c(j)
      rc2(j) = rc(j)*rc(j)*rdphi
      if (j.lt.maxj) then
        dsv(j) = s(j+1) - s(j)
        rdsv(j) = 1._wp/dsv(j)
        cv2(j) = cv(j)*cv(j)*rdsv(j)
        rcv(j) = 1._wp/cv(j)
        if(j.gt.1) rds2(j) = 2._wp/(dsv(j)+dsv(j-1))
      endif
    enddo

    ! meridional grid resolution
    dy = dtheta*R_earth ! m
    rdy = 1._wp/dy
    ! zonal grid resolution at cell center (tracer points)
    do j=1,maxj
      dx(j) = dphi*R_earth*c(j)
      rdx(j) = 1._wp/dx(j)
    enddo
    ! zonal grid resolution on v-grid (meridional cell edges)
    do j=0,maxj
      dxv(j) = dphi*R_earth*cv(j)
      rdxv(j) = 1._wp/dxv(j)
    enddo

    !-------------------------------------------------------------
    ! set up vertical grid
    !-------------------------------------------------------------

    do k=1,maxk
      dz(maxk-k+1) = (zw_in(k+1)-zw_in(k))/zw_in(maxk+1)
    enddo
    depth = zw_in(maxk+1)

    zw(maxk) = 0._wp
    do k=maxk,1,-1
      zw(k-1) = zw(k) - dz(k)
      zro(k) = 0.5_wp*(zw(k) + zw(k-1))
    enddo
    dza(maxk) = 0._wp ! never referenced
    do k=maxk-1,1,-1
      dza(k) = zro(k+1)-zro(k)
    enddo

    zro = zro*depth  ! m
    zw  = zw*depth  ! m
    dz  = dz*depth  ! m
    dza = dza*depth  ! m

    ! define fraction of layers in boundary layer
    do k=maxk,1,-1
      f_pbl(k) = max(0.,min(1.,(zw(k)+dbl)/(zw(k)-zw(k-1))))
    enddo

    print*,'layer #, layer top, layer mid, layer thick, dz layer mid, f_pbl'
    write(6,'(i7,5f12.4)')(k,zw(k),zro(k),dz(k),dza(k),f_pbl(k),k=maxk,1,-1)

    ! define shelf depth, needed? fixme
    k1_shelf = maxk
    tv1 = 1000._wp
    do k=maxk,1,-1
      if (abs(zw(k-1)+shelf_depth).lt.tv1) then
        k1_shelf = k
        tv1 = abs(zw(k-1)+shelf_depth)
      endif
    enddo
    print *,'k1_shelf',k1_shelf
    print *,'depth shelf',zw(k1_shelf-1),zro(k1_shelf)

    ! define 1000 m level depth
    k1_1000 = maxk
    tv1 = 1000._wp
    do k=maxk,1,-1
      if (abs(zw(k-1)+1000._wp).lt.tv1) then
        k1_1000 = k
        tv1 = abs(zw(k-1)+1000._wp)
      endif
    enddo
    print *,'k1_1000',k1_1000
    print *,'depth 1000',zw(k1_1000-1),zro(k1_1000)

    ! define 3000 m level depth
    k1_3000 = maxk
    tv1 = 3000._wp
    do k=maxk,1,-1
      if (abs(zw(k-1)+3000._wp).lt.tv1) then
        k1_3000 = k
        tv1 = abs(zw(k-1)+3000._wp)
      endif
    enddo
    print *,'k1_3000',k1_3000
    print *,'depth 3000',zw(k1_3000-1),zro(k1_3000)

    ! define depth of mixing brines 
    k_mix_brines = maxk
    tv1 = 5000._wp
    do k=maxk,1,-1
      if (abs(zw(k-1)+z_mix_brines).lt.tv1) then
        k_mix_brines = k
        tv1 = abs(zw(k-1)+z_mix_brines)
      endif
    enddo

    !-------------------------------------------------------------
    ! set bathymetry
    !-------------------------------------------------------------

    ! first derive 'potential' bathymetry (k1_pot)
    topo = topo_in

    ! smooth real bathymetry
    if (i_smooth.eq.1) then
      call smooth_topo(topo, int(smooth_fac))
    else if (i_smooth.eq.2) then
      call smooth_topo2(topo, int(smooth_fac))
    else if (i_smooth.eq.3) then
      ! smooth using FFT filter
      call smooth_topo_fft(topo, smooth_fac)
    else if (i_smooth.eq.4) then
      call smooth_topo(topo, int(smooth_fac))
      do j=1,maxj
        do i=1,maxi
          if (j.lt.7 .and.mask_coast(i,j).eq.1) then
            topo(i,j) = z_ocn_max(i,j)
          endif
        enddo
      enddo
    else if (i_smooth.eq.5) then
      call smooth_topo(topo, int(smooth_fac))
      do j=1,maxj
        do i=1,maxi
          if (mask_coast(i,j).eq.1) then
            topo(i,j) = z_ocn_max(i,j)
          endif
        enddo
      enddo
    endif
      

    ! truncate to depth levels and get index of bottom layer in k1_pot
    call truncate_topo(topo,k1_pot)

    ! minimum depth <-> shelf depth
    where (k1_pot.gt.k1_shelf) k1_pot = k1_shelf
    !where (f_ocn.lt.0.5_wp) k1_pot = k1_shelf

    ! derive actual bathymetry applying ocean mask
    k1(1:maxi,1:maxj) = k1_pot(1:maxi,1:maxj)
    ! apply land/sea mask
    where (f_ocn.eq.0._wp)
      k1(1:maxi,1:maxj) = 99
    endwhere

    ! periodic boundary conditions
    k1(0,1:maxj) = k1(maxi,1:maxj)
    k1(maxi+1,1:maxj) = k1(1,1:maxj)
    ! North Pole and South Pole 'islands'
    k1(:,0) = 99
    k1(:,maxj+1) = 99

    ! actual total ocean volume
    ocn_vol_tot = 0._wp
    do j=1,maxj
      do i=1,maxi
        do k=1,maxk
          if (k.ge.k1(i,j)) then
            ocn_vol_tot = ocn_vol_tot + dx(j)*dy*dz(k)*f_ocn(i,j)
          endif
        enddo
      enddo
    enddo

    ! scale depth to get total ocean volume right (matching high resolution bathymetry)
    zfac = ocn_vol_tot_real/ocn_vol_tot

    ! update layer depth and thickness
    zro = zro * zfac  ! m
    zw  = zw  * zfac  ! m
    dz  = dz  * zfac  ! m
    dza = dza * zfac  ! m

    ! truncated bathymetry for output
    do j=1,maxj
      do i=1,maxi
        bathy(i,j) = zw(min(maxk,k1(i,j)-1))
      enddo
    enddo

    ! ocean area and volume
    do j=1,maxj
      do i=1,maxi
        if (k1(i,j).le.maxk) then
          ocn_area(i,j) = dx(j)*dy*f_ocn(i,j)
        else
          ocn_area(i,j) = 0._wp
        endif
        do k=1,maxk
          if (k.ge.k1(i,j)) then
            ocn_vol(i,j,k) = dx(j)*dy*dz(k)*f_ocn(i,j)
          else
            ocn_vol(i,j,k) = 0._wp
          endif
        enddo
      enddo
    enddo

    ! copy to grid type
    grid%ni = maxi
    grid%nj = maxj
    grid%nk = maxk
    grid%lat = lat
    grid%lon = lon

    grid%zro    = zro     
    grid%zw     = zw     
    grid%dz     = dz     
    grid%dza    = dza    
    grid%k1     = k1(1:maxi,1:maxj)
    grid%k1_pot = k1_pot 

    grid%ocn_area = ocn_area
    grid%ocn_vol  = ocn_vol

    ocn_area_tot = 0._wp
    do j=1,maxj
      do i=1,maxi
        if (k1(i,j).le.maxk) then
          ocn_area(i,j) = dx(j)*dy*f_ocn(i,j)
          ocn_area_tot = ocn_area_tot + ocn_area(i,j)
        else
          ocn_area(i,j) = 0._wp
        endif
      enddo
    enddo
    grid%ocn_area_tot = ocn_area_tot

    ! compute total ocean volume for new surface area
    ocn_vol_tot = 0._wp
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (k.ge.k1(i,j)) then
            ocn_vol_tot = ocn_vol_tot + dx(j)*dy*dz(k)*f_ocn(i,j)
          endif
        enddo
      enddo
    enddo
    grid%ocn_vol_tot = ocn_vol_tot

    ! find index of island grid cells
    if (i_isl.eq.2) then
      do n=1,n_isl
        i_idx_isl(n) = 0
        dist_isl = 999.
        do i=1,maxi
          if (modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_isl(n)/180.0)),2.0*pi).lt.dist_isl) then
            i_idx_isl(n) = i
            dist_isl = modulo(abs(phi0+(i-0.5)*dphi-(pi*lon_isl(n)/180.0)),2.0*pi)
          endif
        enddo
        j_idx_isl(n) = 0
        dist_isl = 999.
        do j=1,maxj
          if (abs(s(j)-sin(pi*lat_isl(n)/180.0)).lt.dist_isl) then
            j_idx_isl(n) = j
            dist_isl = abs(s(j)-sin(pi*lat_isl(n)/180.0))
          endif
        enddo
      enddo
    endif

    return

  end subroutine ocn_grid_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  o c n _ g r i d _ u p d a t e
  ! Purpose  :  update ocean grid
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine ocn_grid_update(f_ocn,grid)

    use dim_name, only: dim_lon, dim_lat, dim_time
    use climber_grid, only : lon, lat

    implicit none

    real(wp), dimension(:,:), intent(in) :: f_ocn
    type(grid_class), intent(inout) :: grid

    integer :: i, j, k, kk, ip1, jp1, ncells, n
    integer :: spole, npole
    integer :: flag, sum1, sum2
    integer :: flagi, sum1i, sum2i
    logical :: is_mainland
    logical :: is_island
    real(wp) :: zfac
    real(dp) :: ocn_vol_tot_1000
    real(dp) :: ocn_vol_tot_old

    integer :: ncid
    character (len=256) :: fnm
    real(8) :: empty_time(0)
    logical, save :: l_write_isl = .false.

    grid%ocn_area_old = grid%ocn_area
    grid%ocn_vol_old  = grid%ocn_vol

    !--------------------------------------------------------
    ! update ocean mask
    !--------------------------------------------------------
    mask_ocn = grid%mask_ocn

    !--------------------------------------------------------
    ! update k1
    !--------------------------------------------------------
    do j=1,maxj
      do i=1,maxi
        if (mask_ocn(i,j).eq.1 .and. k1(i,j).gt.maxk) then
          ! new ocean (shelf) cell forming
          k1(i,j) = k1_pot(i,j)
        else if (mask_ocn(i,j).eq.0 .and. k1(i,j).le.maxk) then
          ! ocean cell disappearing
          k1(i,j) = 99
        endif
      enddo
    enddo
    ! periodic boundary conditions
    k1(0,1:maxj) = k1(maxi,1:maxj)
    k1(maxi+1,1:maxj) = k1(1,1:maxj)
    ! North Pole and South Pole 'islands'
    k1(:,0) = 99
    k1(:,maxj+1) = 99

    grid%k1 = k1(1:maxi,1:maxj)

    !--------------------------------------------------------
    ! update ocean surface area and volume
    !--------------------------------------------------------

    ocn_area_tot = 0._wp
    do j=1,maxj
      do i=1,maxi
        if (k1(i,j).le.maxk) then
          ocn_area(i,j) = dx(j)*dy*f_ocn(i,j)
          ocn_area_tot = ocn_area_tot + ocn_area(i,j)
        else
          ocn_area(i,j) = 0._wp
        endif
      enddo
    enddo
    grid%ocn_area_tot = ocn_area_tot
    grid%ocn_area = ocn_area

    ! compute total ocean volume for new surface area
    ocn_vol_tot = 0._wp
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (k.ge.k1(i,j)) then
            ocn_vol_tot = ocn_vol_tot + dx(j)*dy*dz(k)*f_ocn(i,j)
          endif
        enddo
      enddo
    enddo

    !print *,'ocn_vol_tot_real',grid%ocn_vol_tot_real,'ocn_vol_tot before',ocn_vol_tot

    !--------------------------------------------------------
    ! scale depth
    !--------------------------------------------------------

    !-------------------------------------------------------------
    ! set up vertical grid
    !-------------------------------------------------------------

    ! compute total ocean volume below 1000 m
    ocn_vol_tot_1000 = 0._wp
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (k.ge.k1(i,j) .and. k.le.k1_1000) then
            ocn_vol_tot_1000 = ocn_vol_tot_1000 + dx(j)*dy*dz(k)*f_ocn(i,j)
          endif
        enddo
      enddo
    enddo

    ! compute depth factor needed to get required volume change
    zfac = 1._wp + (grid%ocn_vol_tot_real-ocn_vol_tot)/ocn_vol_tot_1000 

    dz(1:k1_1000) = dz(1:k1_1000) * zfac

    zw(maxk) = 0._wp
    do k=maxk,1,-1
      zw(k-1) = zw(k) - dz(k)
      zro(k) = 0.5_wp*(zw(k) + zw(k-1))
    enddo
    dza(maxk) = 0._wp ! never referenced
    do k=maxk-1,1,-1
      dza(k) = zro(k+1)-zro(k)
    enddo

    grid%zro = zro
    grid%zw  = zw 
    grid%dz  = dz 
    grid%dza = dza

    ! update depth dependent variables
    do k=1,maxk-1
     rdz(k) = 1._wp/dz(k)
     rdza(k) = 1._wp/dza(k)
    enddo
    rdz(maxk) = 1._wp/dz(maxk)

    dzz = dz(maxk)*dza(maxk-1)
    do k=maxk,1,-1
     do kk=maxk,1,-1
      dzg(k,kk) = zw(k)-zw(kk-1)
      z2dzg(k,kk) = -zw(k)*zw(k) + zw(kk-1)*zw(kk-1)
      if (k.ne.(kk-1)) then
       rdzg(k,kk) = 1._wp/dzg(k,kk)
      else
       ! This number should be inifinite. Large number used instead
       rdzg(k,kk) = 1.e20_wp
      endif
     enddo
    enddo

    ! truncated bathymetry for output
    do j=1,maxj
      do i=1,maxi
        bathy(i,j) = zw(min(maxk,k1(i,j)-1))
      enddo
    enddo

    !--------------------------------------------------------
    ! update ocean volume
    !--------------------------------------------------------
    ocn_vol_tot = 0._wp
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (k.ge.k1(i,j)) then
            ocn_vol(i,j,k) = dx(j)*dy*dz(k)*f_ocn(i,j)
            ocn_vol_tot = ocn_vol_tot + ocn_vol(i,j,k)
          else
            ocn_vol(i,j,k) = 0._wp
          endif
        enddo
      enddo
    enddo
    grid%ocn_vol_tot = ocn_vol_tot
    grid%ocn_vol = ocn_vol

    ocn_vol_tot_old = sum(grid%ocn_vol_old)

    if (abs(ocn_vol_tot-ocn_vol_tot_old)/ocn_vol_tot.gt.1e-2_wp) then
      grid%l_large_vol_change = .true.
      print *,'WARNING: large change in ocean volume!'
      print *,'use of horizontal diffusion enforced for all tracers for the current year'
      print *,'old ocean volume',ocn_vol_tot_old
      print *,'new ocean volume',ocn_vol_tot
    else
      grid%l_large_vol_change = .false.
    endif

    !print *,'ocn_vol_tot_real',grid%ocn_vol_tot_real,'ocn_vol_tot after',ocn_vol_tot

    !--------------------------------------------------------
    ! update ocean k1 masks
    !--------------------------------------------------------

    ! update ocean mask on c-grid
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          if (k.ge.k1(i,j)) then
            mask_c(i,j,k) = 1
          else
            mask_c(i,j,k) = 0
          endif
        enddo
      enddo
    enddo
    grid%mask_c = mask_c
    ! update ocean mask on velocity grids
    do k=1,maxk
      do j=1,maxj
        do i=1,maxi
          ! update ocean mask on u-grid
          ip1 = i+1
          if (ip1.eq.maxi+1) ip1=1
          if (k.ge.max(k1(i,j),k1(ip1,j))) then
            mask_u(i,j,k) = 1
          else
            mask_u(i,j,k) = 0
          endif
          ! update ocean mask on v-grid
          jp1 = min(maxj,j+1)
          if (j.lt.maxj .and. k.ge.max(k1(i,j),k1(i,jp1))) then
            mask_v(i,j,k) = 1
          else
            mask_v(i,j,k) = 0
          endif
          ! update ocean mask on w-grid
          if (k.ge.k1(i,j) .and. k.lt.maxk) then
            mask_w(i,j,k) = 1
          else
            mask_w(i,j,k) = 0
          endif
        enddo
      enddo
    enddo

    ! update seabed depth h
    h = 0._wp
    rh = 0._wp
    do j=maxj+1,0,-1
     do i=0,maxi+1
      if(k1(i,j).le.maxk) then
       do k=k1(i,j),maxk
        h(3,i,j) = h(3,i,j) + dz(k)
       enddo
       rh(3,i,j) = 1._wp/h(3,i,j)
      endif
     enddo
    enddo
    do j=0,maxj+1
     do i=0,maxi
      h(1,i,j) = min(h(3,i,j),h(3,i+1,j))
      if (max(k1(i,j),k1(i+1,j)).le.maxk) rh(1,i,j) = 1._wp/h(1,i,j)
     enddo
    enddo
    do j=0,maxj
     do i=0,maxi+1
      h(2,i,j) = min(h(3,i,j),h(3,i,j+1))
      if (max(k1(i,j),k1(i,j+1)).le.maxk) rh(2,i,j) = 1._wp/h(2,i,j)
     enddo
    enddo

    ! derive mask on velocity grid
    do j=1,maxj
     do i=1,maxi
      ku(1,i,j) = max(k1(i,j),k1(i+1,j))
      ku(2,i,j) = max(k1(i,j),k1(i,j+1))
     enddo
    enddo

    ! array to avoid J term in flat regions
    ! For non-trivial coasts essential to avoid adding J term at non-Psi points.
    ! Hence first condition ensures (i,j) is a wet Psi point, 2nd that bottom is not flat.
    do j=1,maxj
       do i=1,maxi
          if ((max(k1(i,j),k1(i+1,j),k1(i,j+1),k1(i+1,j+1)).le.maxk) .and. &
             (k1(i,j).ne.k1(i,j+1).or.k1(i,j).ne.k1(i+1,j) .or.k1(i,j).ne.k1(i+1,j+1))) then
             getj(i,j) = .true.
          else
             getj(i,j) = .false.
          endif
       enddo
    enddo

    !--------------------------------------------------------
    ! update islands 
    !--------------------------------------------------------

    where (k1(1:maxi,1:maxj).gt.maxk) 
      map = 1    ! land
      map_isles = -1
    elsewhere
      map = 0     ! ocean
      map_isles = 0
    endwhere

    ! try to connect lonely land cells to larger islands close by
    call lonely_cell(f_ocn,map,map_isles)

    ! identify isles
    if (i_isl.eq.0) then
      n_isles = 0
    else
      n_isles = 1
    endif
    do i=1,maxi
      do j=1,maxj
        if (map_isles(i,j).eq.-1) then   ! unclaimed land,found new island
          ! initialize
          islands = 0
          edge = 0
          call fill_island(i,j)   ! fill current island
          ! count number of cells in island
          ncells = count(islands==1)
          is_island = .false.
          is_mainland = .false.

          if (i_isl.eq.1) then
            ! manual island determination
            if (ncells>10) then
              if (islands(53,27).eq.1) then
                ! Eurasia, set to mainland, island 1
                is_island = .true.
                is_mainland = .true.
              endif
              if (l_isl_ant .and. islands(36,1).eq.1) then
                ! Antarctica
                is_island = .true.
              endif
              if (l_isl_aus .and. islands(63,13).eq.1) then
                ! Australia
                is_island = .true.
              endif
              if (l_isl_grl .and. islands(28,33).eq.1) then
                ! Greenland
                is_island = .true.
              endif
              if (l_isl_ame .and. islands(16,28).eq.1) then
                ! Greenland
                is_island = .true.
              endif
            endif
            if (is_island) then
              ! create new island
              if (is_mainland) then
                ! Mainland, island 1
                map_edge(:,:,1) = edge
                where (islands.eq.1)
                  map_isles = 1
                endwhere
              else
                n_isles = n_isles+1 ! increase index
                map_edge(:,:,n_isles) = edge
                where (islands.eq.1)
                  map_isles = n_isles
                endwhere
              endif
            else
              ! add island to mainland (island 1)
              where (islands.eq.1)
                map_isles = 1
              endwhere
            endif
          endif

          if (i_isl.eq.2) then
            ! manual island determination
            do n=1,n_isl
              if (islands(i_idx_isl(n),j_idx_isl(n)).eq.1) then
                is_island = .true.
                if (n.eq.1) is_mainland = .true.
              endif
            enddo
            if (is_island) then
              ! create new island
              if (is_mainland) then
                ! Mainland, island 1
                map_edge(:,:,1) = edge
                where (islands.eq.1)
                  map_isles = 1
                endwhere
              else
                n_isles = n_isles+1 ! increase index
                map_edge(:,:,n_isles) = edge
                where (islands.eq.1)
                  map_isles = n_isles
                endwhere
              endif
            else
              ! add island to mainland (island 1)
              where (islands.eq.1)
                map_isles = 1
              endwhere
            endif
          endif

          if (i_isl.eq.0) then
            ! automatic island determination
            if (ncells>50) then 
              ! create new island
              n_isles = n_isles+1 ! increase index
              map_edge(:,:,n_isles) = edge
              where (islands.eq.1)
                map_isles = n_isles
              endwhere
            else
              ! add island to mainland (island 1)
              where (islands.eq.1)
                map_isles = 1
              endwhere
            endif
          endif

        endif
      enddo
    enddo

    ! number of isles, remove the mainland
    n_isles = n_isles-1

    ! east-west axis is mirrored and polar values are land (this bit is tricky and may
    ! need hand-editting if there are complex polar islands)
    map_isles_ext(1:maxi,1:maxj) = map_isles
    spole = maxval(map_isles_ext(:,1:2))
    npole = maxval(map_isles_ext(:,maxj-1:maxj))
    map_isles_ext(:,0) = spole
    map_isles_ext(:,maxj+1) = npole
    map_isles_ext(0,:) = map_isles_ext(maxi,:)
    map_isles_ext(maxi+1,:) = map_isles_ext(1,:)

    ! islands on psi grid
    do i=1,maxi
      do j=0,maxj
        psiles(i + j*maxi) = max(map_isles_ext(i,j),map_isles_ext(i+1,j),map_isles_ext(i,j+1),map_isles_ext(i+1,j+1))
      enddo
    enddo

    !--------------------------------------------------------
    ! process island edges
    !--------------------------------------------------------

    do i=1,n_isles+1
      flag = 0
      do while (flag==0)
        sum1 = sum(map_edge(:,:,i))
        flagi = 0
        ! get rid of dangling edge points with only one neighbor in its four cell neighborhood
        do while (flagi==0)
          sum1i = sum(map_edge(:,:,i))
          call lonely_cell_edge(map_edge(:,:,i))
          sum2i = sum(map_edge(:,:,i))
          if (sum2i==sum1i) flagi = 1
        enddo
        flagi = 0
        ! clean edge, fix cells with too many neighbors
        do while (flagi==0)
          sum1i = sum(map_edge(:,:,i))
          call clean_edge(map_edge(:,:,i))
          sum2i = sum(map_edge(:,:,i))
          if (sum2i==sum1i) flagi = 1
        enddo
        sum2 = sum(map_edge(:,:,i))
        if (sum2==sum1) flag = 1
      enddo
    enddo

    !--------------------------------------------------------
    ! create paths from edges
    !--------------------------------------------------------

    do i = 2,n_isles+1
      call create_path(map_edge(:,:,i), npi(i-1), ipisl(:,i-1), jpisl(:,i-1), lpisl(:,i-1), map_path(:,:,i))
    enddo

    do i=1,n_isles
      if(npi(i).gt.mpi) then
        print *,'path integral around island too long'
        l_write_isl = .true.
      endif
      do j=1,npi(i)
        if (abs(lpisl(j,i)).ne.1.and.abs(lpisl(j,i)).ne.2) then
          print *,'bad path'
          print *,j,i
          print *,lpisl(j,i)
          print *,k1(ipisl(j,i),jpisl(j,i)),maxk
        l_write_isl = .true.
        endif
        if(ipisl(j,i).gt.maxi.or.ipisl(j,i).lt.0) then
          print *,'bad path' 
          print *,j,i
          print *,ipisl(j,i)
          print *,k1(ipisl(j,i),jpisl(j,i)),maxk
        l_write_isl = .true.
        end if
        if(jpisl(j,i).gt.maxj.or.jpisl(j,i).lt.0) then
          print *,'bad path' 
          print *,j,i
          print *,jpisl(j,i)
          print *,k1(ipisl(j,i),jpisl(j,i)),maxk
        l_write_isl = .true.
        end if
        if(k1(ipisl(j,i),jpisl(j,i)).gt.maxk) then
          print *,'dry path',j,i,ipisl(j,i),jpisl(j,i),k1(ipisl(j,i),jpisl(j,i)),maxk
        l_write_isl = .true.
        endif
      enddo
    enddo

    if (l_write_isl) then
    fnm = trim(out_dir)//"/check_islands.nc"
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,dim_time, x=empty_time, units="years BP", unlimited=.TRUE.,ncid=ncid)
    call nc_write_dim(fnm,dim_lon,x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,dim_lat,x=lat,axis="y",ncid=ncid)
    call nc_write_dim(fnm,"nisles",x=[(i, i=1,maxisles)],axis="z",ncid=ncid)
    call nc_close(ncid)

    call nc_open(fnm,ncid)
    call nc_write(fnm,dim_time,    real(year+1,wp), dim1=dim_time,start=[year+1],count=[1],ncid=ncid)    
    call nc_write(fnm,"k1_pot",     k1_pot,dims=[dim_lon,dim_lat,dim_time],start=[1,1,year+1],count=[maxi,maxj,1],long_name="density",units="?",ncid=ncid)
    call nc_write(fnm,"k1",     k1(1:maxi,1:maxj),dims=[dim_lon,dim_lat,dim_time],start=[1,1,year+1],count=[maxi,maxj,1],long_name="density",units="?",ncid=ncid)
    call nc_write(fnm,"mask_ocn",     real(mask_ocn,wp),dims=[dim_lon,dim_lat,dim_time],start=[1,1,year+1],count=[maxi,maxj,1],long_name="density",units="?",ncid=ncid)
    call nc_write(fnm,"map_isles",     map_isles,dims=[dim_lon,dim_lat,dim_time],start=[1,1,year+1],count=[maxi,maxj,1],long_name="density",units="?",ncid=ncid)
    call nc_write(fnm,"map_edge",     map_edge,dims=[dim_lon,dim_lat,"nisles",dim_time],start=[1,1,1,year+1],count=[maxi,maxj,maxisles,1],long_name="density",units="?",ncid=ncid)
    call nc_write(fnm,"map_path",     map_path,dims=[dim_lon,dim_lat,"nisles",dim_time],start=[1,1,1,year+1],count=[maxi,maxj,maxisles,1],long_name="density",units="?",ncid=ncid)
    call nc_close(ncid)

    stop
  endif


  end subroutine ocn_grid_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s m o o t h _ t o p o
  ! Purpose  :  smooth bathymetry
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth_topo(topo, nsmooth)

  implicit none

    real(wp), intent(inout) :: topo(:,:)
    integer, intent(in) :: nsmooth

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp), allocatable :: Y(:,:)


    im = size(topo,1)
    jm = size(topo,2)
    allocate(Y(im,jm))

    do iter=1,nsmooth

      do j=1,jm
        do i=1,im
          Y(i,j)=topo(i,j)
        enddo
      enddo        

      do j=2,jm-1       
        do i=1,im
          imi=i-1
          if (imi.eq.0)imi=im
          ipl=i+1
          if (ipl.gt.im)ipl=1
          topo(i,j)=0.2*(Y(i,j)+Y(imi,j)+Y(ipl,j)+Y(i,j+1)+Y(i,j-1))
        enddo
      enddo

    enddo       

    deallocate(Y)

  end subroutine smooth_topo

  subroutine smooth_topo2(topo, nsmooth)

  implicit none

    real(wp), intent(inout) :: topo(:,:)
    integer, intent(in) :: nsmooth

    integer :: iter, i, j, imi, ipl, im, jm
    real(wp) :: wx, wy, wsum
    real(wp), allocatable :: Y(:,:)


    im = size(topo,1)
    jm = size(topo,2)
    allocate(Y(im,jm))

    do iter=1,nsmooth

      do j=1,jm
        do i=1,im
          Y(i,j)=topo(i,j)
        enddo
      enddo        

      j=1
      do i=1,im
        imi=i-1
        if (imi.eq.0)imi=im
        ipl=i+1
        if (ipl.gt.im)ipl=1
        topo(i,j)=0.25*(Y(i,j)+Y(imi,j)+Y(ipl,j)+Y(i,j+1))
      enddo

      do j=2,jm-1       
        do i=1,im
          imi=i-1
          if (imi.eq.0)imi=im
          ipl=i+1
          if (ipl.gt.im)ipl=1
          wx = 1._wp/dx(j)
          wy = 1._wp/dy
          wsum = wx+wy
          topo(i,j)=0.2_wp*Y(i,j) + 0.8_wp*(wx*0.5_wp*(Y(imi,j)+Y(ipl,j))+wy*0.5_wp*(Y(i,j+1)+Y(i,j-1)))/wsum
        enddo
      enddo

      j=jm
      do i=1,im
        imi=i-1
        if (imi.eq.0)imi=im
        ipl=i+1
        if (ipl.gt.im)ipl=1
        topo(i,j)=0.25*(Y(i,j)+Y(imi,j)+Y(ipl,j)+Y(i,j-1))
      enddo

    enddo       

    deallocate(Y)

  end subroutine smooth_topo2


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s m o o t h _ t o p o _ f f t
  ! Purpose  :  smooth bathymetry using filter in fourier transform space
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smooth_topo_fft(topo, smooth_fac)

  use, intrinsic :: iso_c_binding 

  implicit none

  include 'fftw3.f03'


    real(wp), intent(inout) :: topo(:,:)
    real(wp), intent(in) :: smooth_fac

    integer :: i, j, nx, ny, nx2, ny2
    integer*8 :: plan
    real(wp) :: X, Y, R, R0, lambda
    real(wp), allocatable :: topo_ext(:,:)
    real(wp), allocatable, dimension(:,:) :: M
    complex(wp), allocatable :: topo_fft(:,:)

    character(len=256) :: fnm
    integer :: ncid

    nx = size(topo,1)
    ny = size(topo,2)

    nx2 = nx/2
    ny2 = ny/2

    allocate(topo_ext(nx*2,ny*2))
    topo_ext = 0._wp
    ! Put data matrix into centre of frame
    topo_ext((nx2+1):(nx2+nx),(ny2+1):(ny2+ny)) = topo(1:nx,1:ny)
    ! Repeat data matrix either side (i.e. where it wraps east-west) 
    topo_ext(1:nx2,(ny2+1):(ny2+ny)) = topo(nx2+1:nx,1:ny)
    topo_ext((nx2+nx+1):(2*nx),(ny2+1):(ny2+ny)) = topo(1:nx2,1:ny)

    !
    ! Now, need to do something clever to remove low amplitude/high frequency stuff
    ! Here I'm creating a circular mask
    allocate(M(nx*2,ny*2))
    ! R0 is the radius of the circle within which FFT components are destroyed
    R0 = nx * smooth_fac
    lambda = (5._wp/360._wp)*nx
    do i=1,2*nx
      do j=1,2*ny
        X = i - nx - 0.5_wp
        Y = j - ny - 0.5_wp
        R = sqrt(X**2 + Y**2)
        M(i,j) = (tanh((R - R0)/lambda)+1._wp)/2._wp
      enddo
    enddo

    !    fnm = "test_fft.nc"
    !    call nc_create(fnm)
    !    call nc_open(fnm,ncid)
    !    call nc_write_dim(fnm,"nx",x=1._wp,dx=1._wp,nx=nx,ncid=ncid)
    !    call nc_write_dim(fnm,"nx1",x=1._wp,dx=1._wp,nx=nx+1,ncid=ncid)
    !    call nc_write_dim(fnm,"nx2",x=1._wp,dx=1._wp,nx=nx*2,ncid=ncid)
    !    call nc_write_dim(fnm,"ny",x=1._wp,dx=1._wp,nx=ny,ncid=ncid)
    !    call nc_write_dim(fnm,"ny2",x=1._wp,dx=1._wp,nx=ny*2,ncid=ncid)
    !    call nc_write(fnm,"topo",     topo,dims=["nx","ny"],start=[1,1],count=[nx,ny],ncid=ncid)
    !    call nc_write(fnm,"topo_ext",     topo_ext,dims=["nx2","ny2"],start=[1,1],count=[nx*2,ny*2],ncid=ncid)
    !    call nc_write(fnm,"M",     M,dims=["nx2","ny2"],start=[1,1],count=[nx*2,ny*2],ncid=ncid)
    !    call nc_close(ncid)


    allocate(topo_fft(nx+1, 2*ny))

    !  Make a plan for the FFT, and forward transform the data.
    !
    call dfftw_plan_dft_r2c_2d(plan,nx*2,ny*2,topo_ext,topo_fft,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, topo_ext, topo_fft)
    call dfftw_destroy_plan(plan)

    !    call nc_open(fnm,ncid)
    !    call nc_write(fnm,"topo_fft_re",     real(topo_fft),dims=["nx1","ny2"],start=[1,1],count=[nx+1,ny*2],ncid=ncid)
    !    call nc_write(fnm,"topo_fft_im",     imag(topo_fft),dims=["nx1","ny2"],start=[1,1],count=[nx+1,ny*2],ncid=ncid)
    !    call nc_write(fnm,"topo_fft",     real(topo_fft)**2+imag(topo_fft)**2,dims=["nx1","ny2"],start=[1,1],count=[nx+1,ny*2],ncid=ncid)
    !    call nc_close(ncid)

    ! Apply this to the real and imaginary components
    do i=1,nx
      do j=1,2*ny
        topo_fft(i,j) = topo_fft(i,j) * M(i,j)
      enddo
    enddo

    !    call nc_open(fnm,ncid)
    !    call nc_write(fnm,"topo_fft_fil",     real(topo_fft)**2+imag(topo_fft)**2,dims=["nx1","ny2"],start=[1,1],count=[nx+1,ny*2],ncid=ncid)
    !    call nc_close(ncid)

    !  Make a plan for the backward FFT, and back transform.
    !
    call dfftw_plan_dft_c2r_2d(plan, nx*2, ny*2, topo_fft, topo_ext, FFTW_ESTIMATE )
    call dfftw_execute_dft_c2r(plan,topo_fft,topo_ext )
    call dfftw_destroy_plan(plan)

    ! De-frame the data matrix
    topo = topo_ext((nx2+1):(nx2+nx),(ny2+1):(ny2+ny))

    ! normalize
    topo = topo/real(nx*2*ny*2,wp)


    !    call nc_open(fnm,ncid)
    !    call nc_write(fnm,"topo_inv",     topo,dims=["nx","ny"],start=[1,1],count=[nx,ny],ncid=ncid)
    !    call nc_write(fnm,"topo_ext_inv",     topo_ext,dims=["nx2","ny2"],start=[1,1],count=[nx*2,ny*2],ncid=ncid)
    !    call nc_close(ncid)

    deallocate(topo_fft)
    deallocate(M)

  end subroutine smooth_topo_fft


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  t r u n c a t e _ t o p o
  ! Purpose  :  truncate bathymetry to depth levels
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine truncate_topo(topo,k1)

    implicit none

    real(wp), intent(in) :: topo(:,:)
    integer, intent(out) :: k1(:,:)
    
    integer :: k


    k1 = 99
    do k=maxk,1,-1
      where (topo < zro(k)) k1 = k
    enddo
    where (topo > zro(maxk) .and. topo < 0) k1 = maxk


  end subroutine truncate_topo


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  f i l l _ i s l a n d
  ! Purpose  :  identify islands
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  recursive subroutine fill_island(i,j) 

    implicit none

    integer, intent(in) :: i, j

    integer :: im1, ip1, jm1, jp1


    if (map(i,j) == 1 .and. islands(i,j) == 0) then

      ! This is a contiguous cell - index it
      islands(i,j) = 1

      ! Check its perimeters for cells with a different value
      if (i > 1) then
        if (map((i-1),j).ne.1) edge((i-1),j) = 1 
      endif
      if (i > 1 .and. j > 1) then
        if (map((i-1),(j-1)) .ne. 1) edge((i-1),(j-1)) = 1
      endif
      if (i > 1 .and. j < maxj) then
        if (map((i-1),(j+1)) .ne. 1) edge((i-1),(j+1)) = 1
      endif
      if (i < maxi) then
        if (map(i+1,j) .ne. 1) edge(i+1,j) = 1
      endif
      if (i < maxi .and. j > 1) then
        if (map((i+1),(j-1)) .ne. 1) edge((i+1),(j-1)) = 1
      endif
      if (i < maxi .and. j < maxj) then
        if (map((i+1),(j+1)) .ne. 1) edge((i+1),(j+1)) = 1
      endif
      if (j > 1) then
        if (map(i,j-1) .ne. 1) edge(i,j-1) = 1
      endif
      if (j < maxj) then
        if (map(i,j+1) .ne. 1) edge(i,j+1) = 1
      endif

      ! And again for east-west boundary 
      if (i == 1) then
        if (map(maxi,j) .ne. 1) edge(maxi,j) = 1
        if (j > 1) then
          if (map(maxi,(j-1)) .ne. 1) edge(maxi,(j-1)) = 1
        endif
        if (j < maxj) then
          if (map(maxi,(j+1)) .ne. 1) edge(maxi,(j+1)) = 1
        endif
      endif
      if (i == maxi) then
        if (map(1,j) .ne. 1) edge(1,j) = 1
        if (j > 1) then
          if (map(1,(j-1)) .ne. 1) edge(1,(j-1)) = 1
        endif
        if (j < maxj) then
          if (map(1,(j+1)) .ne. 1) edge(1,(j+1)) = 1
        endif
      endif

      ! Now recursively search neighbouring cells
      ! West
      im1 = i-1
      if (im1 > 0) then
        call fill_island(im1,j)
      endif
      ! East
      ip1 = i+1
      if (ip1 < (maxi+1)) then
        call fill_island(ip1,j)
      endif
      ! South
      jm1 = j-1
      if (jm1 > 0) then
        call fill_island(i,jm1)
      endif
      ! North
      jp1 = j+1
      if (jp1 < (maxj+1)) then
        call fill_island(i,jp1)
      endif
      ! North-west
      im1 = i-1
      jp1 = j+1
      if (im1 > 0 .and. jp1 < (maxj+1)) then
        call fill_island(im1,jp1)
      endif
      ! South-west
      im1 = i-1
      jm1 = j-1
      if (im1 > 0 .and. jm1 > 0) then
        call fill_island(im1,jm1)
      endif
      ! North-east
      ip1 = i+1
      jp1 = j+1
      if (ip1 < (maxi+1) .and. jp1 < (maxj+1)) then
        call fill_island(ip1,jp1)
      endif
      ! South-east
      ip1 = i+1
      jm1 = j-1
      if (ip1 < (maxi+1) .and. jm1 > 0) then
        call fill_island(ip1,jm1)
      endif

      ! wrap east-west
      ! Jump west
      if (i == 1) then
        call fill_island(maxi,j)
        jm1 = j-1
        if (jm1 > 0) then
          call fill_island(maxi,jm1)
        endif
        jp1 = j+1
        if (jp1 < (maxj+1)) then
          call fill_island(maxi,jp1)
        endif
      endif
      ! Jump east
      if (i == maxi) then
        call fill_island(1,j)
        jm1 = j-1
        if (jm1 > 0) then
          call fill_island(1,jm1)
        endif
        jp1 = j+1
        if (jp1 < (maxj+1)) then
          call fill_island(1,jp1)
        endif
      endif

    endif

  end subroutine fill_island


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l o n e l y _ c e l l
  ! Purpose  :  find lonely land cells and try to connect them to larger islands nearby
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lonely_cell(f_ocn,map,map_isles)

    implicit none

    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(inout) :: map(:,:)
    integer, intent(inout) :: map_isles(:,:)

    integer :: i, j, ii, jj, iii, jjj, i1, j1, ii1, jj1, iii1, jjj1, iiii1, jjjj1
    integer :: ncell_ocn, n, nfill, ncell_lnd(8), idx(8), jdx(8)
    real(wp) :: focn(8)


    ! Loop through cells counting 9 cell neighbours
    do i=1,maxi
      do j=1,maxj
        if (map(i,j).eq.1) then ! land point
          ! count number of ocean neighbors
          ncell_ocn = 0
          do ii=i-1,i+1
            do jj=j-1,j+1
              iii = ii
              if (iii.eq.0) iii = maxi
              if (iii.eq.maxi+1) iii = 1
              jjj = jj
              jjj = max(1,jjj)
              jjj = min(maxj,jjj)
              if (map(iii,jjj).eq.0) then !  ocean
                ncell_ocn = ncell_ocn+1
              endif
            enddo
          enddo
          if (ncell_ocn.eq.8) then  ! isolated land point, all neighbors are ocean
            ! count number of land point neighbors of neighbors
            n=0
            ncell_lnd(:) = 0
            focn(:) = 0._wp
            do i1=i-1,i+1
              do j1=j-1,j+1
                ii1 = i1
                if (ii1.eq.0) ii1 = maxi
                if (ii1.eq.maxi+1) ii1 = 1
                jj1 = j1
                jj1 = max(1,jj1)
                jj1 = min(maxj,jj1)
                if (.not.(ii1.eq.i .and. jj1.eq.j)) then
                  n=n+1
                  idx(n) = ii1
                  jdx(n) = jj1
                  focn(n) = f_ocn(ii1,jj1)
                  do iii1=ii1-1,ii1+1
                    do jjj1=jj1-1,jj1+1
                      iiii1 = iii1
                      if (iiii1.eq.0) iiii1 = maxi
                      if (iiii1.eq.maxi+1) iiii1 = 1
                      jjjj1 = jjj1
                      jjjj1 = max(1,jjjj1)
                      jjjj1 = min(maxj,jjjj1)
                      if (map(iiii1,jjjj1).eq.1) then !  land
                        ncell_lnd(n) = ncell_lnd(n)+1
                        ncell_lnd(n) = min(2,ncell_lnd(n))
                      endif
                    enddo
                  enddo
                endif
              enddo
            enddo
            ! find index of cell with most land neighbors and consider also land fraction
            nfill = maxloc(ncell_lnd+(1._wp-focn),1)
            ! turn ocean cell to land to connect isolated land point
            map(idx(nfill),jdx(nfill)) = 1
            map_isles(idx(nfill),jdx(nfill)) = -1
          endif
        endif
      enddo
    enddo


  end subroutine lonely_cell


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  l o n e l y _ c e l l _ e d g e
  ! Purpose  :  find and remove lonely cells in island edges
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lonely_cell_edge(map_edge)

    implicit none

    integer, intent(inout) :: map_edge(:,:)

    integer, allocatable :: map_edge_new(:,:)

    integer :: i, j, nbors


    allocate(map_edge_new, MOLD=map_edge)

    ! Set up output array
    map_edge_new = map_edge

    ! Loop through cells counting four cell neighbours
    do i=1,maxi
      do j=1,maxj
        nbors = 0
        ! Standard four cell neighbourhood
        if (j < maxj) then
          if (map_edge(i,j+1) == 1) nbors = nbors + 1
        endif
        if (i < maxi) then
          if (map_edge(i+1,j) == 1) nbors = nbors + 1
        endif
        if (j > 1) then
          if (map_edge(i,j-1) == 1) nbors = nbors + 1
        endif
        if (i > 1) then
          if (map_edge(i-1,j) == 1) nbors = nbors + 1
        endif
        ! East-west boundary cells
        if (i == 1) then
          if (map_edge(maxi,j) == 1) nbors = nbors + 1
        endif
        if (i == maxi) then 
          if (map_edge(1,j) == 1) nbors = nbors + 1
        endif
        ! If no or one neighbour, cell dies
        if (nbors<=1) map_edge_new(i,j) = 0
      enddo
    enddo

    map_edge = map_edge_new

    deallocate(map_edge_new)

  end subroutine lonely_cell_edge


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c l e a n _ e d g e
  ! Purpose  :  clean island edges with too many neighbors
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine clean_edge(map_edge)

    implicit none

    integer, intent(inout) :: map_edge(:,:)

    integer, allocatable :: map_nbors(:,:)
    integer, allocatable :: map_sqr(:,:)

    integer :: i, j, ip1, nbors


    allocate(map_nbors, MOLD=map_edge)
    allocate(map_sqr, MOLD=map_edge)

    ! Loop through cells counting four cell neighbours
    do i=1,maxi
      do j=1,maxj
        if (map_edge(i,j).eq.1) then
          nbors = 0
          ! Standard four cell neighbourhood
          if (j < maxj) then
            if (map_edge(i,j+1) == 1) nbors = nbors + 1
          endif
          if (i < maxi) then
            if (map_edge(i+1,j) == 1) nbors = nbors + 1
          endif
          if (j > 1) then
            if (map_edge(i,j-1) == 1) nbors = nbors + 1
          endif
          if (i > 1) then
            if (map_edge(i-1,j) == 1) nbors = nbors + 1
          endif
          ! East-west boundary cells
          if (i == 1) then
            if (map_edge(maxi,j) == 1) nbors = nbors + 1
          endif
          if (i == maxi) then 
            if (map_edge(1,j) == 1) nbors = nbors + 1
          endif
          ! store neighbor count 
          map_nbors(i,j) = nbors
        else
          map_nbors(i,j) = 0
        endif
      enddo
    enddo

    do i=1,maxi
      ip1 = i+1
      if (ip1==maxi+1) ip1=1
      do j=1,maxj-1
        if (map_edge(i,j)==1 .and. map_edge(ip1,j)==1 .and. map_edge(i,j+1)==1 .and. map_edge(ip1,j+1)==1) then
          ! square of edge points found
          ! remove from edge those points in the square that have less than 3 neighbors
          if (map_nbors(i,  j  )<3) map_edge(i,  j  ) = 0
          if (map_nbors(ip1,j  )<3) map_edge(ip1,j  ) = 0
          if (map_nbors(i,  j+1)<3) map_edge(i,  j+1) = 0
          if (map_nbors(ip1,j+1)<3) map_edge(ip1,j+1) = 0
        endif
      enddo
    enddo

    deallocate(map_sqr, map_nbors)

  end subroutine clean_edge


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  c r e a t e _ p a t h
  ! Purpose  :  check island paths and save
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine create_path(map_edge, &
                        npi, ipisl, jpisl, lpisl, map_path)

    implicit none

    integer, intent(in) :: map_edge(:,:)

    integer, intent(out) :: npi
    integer, dimension(:), intent(out) :: ipisl, jpisl, lpisl
    integer, dimension(:,:), intent(out) :: map_path

    integer :: i, j, n, firstx, firsty, lastx, lasty
    integer :: flag, endflag, lastcell
    integer, dimension(4) :: nbors


    ipisl = 0
    jpisl = 0
    lpisl = 0

    ! Work out most NW edge cell
    flag = 0
    do i=1,maxi
      do j=maxj,1,-1
        if (flag == 0 .and. map_edge(i,j).eq.1) then
          firstx = i 
          firsty = j
          flag = 1
        endif
      enddo
    enddo

    ! Initialise path arrays
    lastx = firstx
    lasty = firsty
    n = 1
    endflag = 0

    ! Initialise path variables
    map_path = 0
    lastcell = 1

    ! Loop through cells finding four cell neighbours
    do while (endflag == 0)
      i = lastx 
      j = lasty
      map_path(lastx, lasty) = n
      nbors = 0
      ! Standard four cell neighbourhood
      if (i < maxi) then
        if (map_edge(i+1,j) == 1) nbors(2) = 1
      endif
      if (j < maxj) then
        if (map_edge(i,j+1) == 1) nbors(1) = 1
      endif
      if (i > 1) then
        if (map_edge(i-1,j) == 1) nbors(4) = 1
      endif
      if (j > 1) then 
        if (map_edge(i,j-1) == 1) nbors(3) = 1
      endif
      ! East-west boundary cells
      if (i == 1) then
        if (map_edge(maxi,j) == 1) nbors(4) = 1
      endif
      if (i == maxi) then
        if (map_edge(1,j) == 1) nbors(2) = 1
      endif
      ! Blank neighbours if this is the first cell
      if (n == 1) then
        nbors(1) = 0
        nbors(4) = 0
      endif
      ! Blank last neighbour
      nbors(lastcell) = 0
      ! Check that only one neighbour remains
      if (sum(nbors) > 1) then
        nbors(4) = 0
      endif
      if (sum(nbors) > 1) then
        nbors(1) = 0
      endif
      if (sum(nbors) > 1) then
        nbors(2) = 0
      endif
      if (sum(nbors) > 1) stop
      ! Act on neighbourhood information
      ! Need to note the interface rather than the cell!
      if (nbors(1) == 1) then
        ! North
        ipisl(n) = lastx
        jpisl(n) = lasty
        lpisl(n) = +2
        lastcell = 3
        lastx = lastx
        lasty = lasty + 1
      else if (nbors(2) == 1) then
        ! East
        ipisl(n) = lastx
        jpisl(n) = lasty
        lpisl(n) = +1
        lastcell = 4
        lastx = lastx + 1
        lasty = lasty 
      else if (nbors(3) == 1) then
        ! South
        ipisl(n) = lastx
        jpisl(n) = lasty - 1
        lpisl(n) = -2
        lastcell = 1
        lastx = lastx
        lasty = lasty - 1
      else if (nbors(4) == 1) then
        ! West
        ipisl(n) = lastx - 1
        jpisl(n) = lasty
        lpisl(n) = -1
        lastcell = 2
        lastx = lastx - 1
        lasty = lasty
      endif
      ! Check if this means the cell has cross the east-west boundary
      if (lastx == maxi + 1) then
        lastx = 1
      else if (lastx == 0) then
        lastx = maxi
      endif
      if (ipisl(n) == maxi + 1) then
        ipisl(n) = 1
      else if (ipisl(n) == 0) then
        ipisl(n) = maxi
      endif
      ! Check if we're back at the start cell
      if (lastx == firstx .and. lasty == firsty) then
        endflag = 1
      endif
      n = n + 1

    enddo

    npi = n-1

  end subroutine create_path


end module ocn_grid

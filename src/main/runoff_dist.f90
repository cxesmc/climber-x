!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : r u n o f f _ d i s t _ m o d
!
!  Purpose : ocean cells over which the freshwater (and weathering)
!            fluxes reaching a coastal destination cell are spread
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2026 Potsdam Institute for Climate Impact Research,
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
module runoff_dist_mod

  ! For every coastal destination cell (ir,jr) a list of receiving ocean cells
  ! n=1..nc with a fraction r(n) of the ocean area of each cell is built.
  ! A flux F (per unit time, e.g. kg/s) arriving at (ir,jr) is applied as a
  ! uniform flux density over the receiving ocean area A_sum = sum(r*area*f_ocn):
  !
  !   fld(i(n),j(n)) += F * r(n) / A_sum      [F/m2], n=1..nc
  !
  ! which conserves F exactly.
  !
  ! i_mode = 1 : the first n_cells entries of the coastal neighbour list from
  !              coast_cells (geo), r = 1 (the original CLIMBER-X scheme)
  ! i_mode = 2 : a fixed receiving ocean area A_dist [m2], independent of the
  !              grid resolution and of the ocean fraction of coastal cells.
  !              Ocean cells are added in order of increasing distance from
  !              the destination cell, measured along paths through ocean
  !              cells (Dijkstra on the 8-neighbour grid, diagonal steps
  !              only when one of the two orthogonal cells is ocean, so the
  !              flux cannot jump across land into another basin). The last
  !              cell is used only partially (r<1) so that A_sum = A_dist
  !              exactly. If the connected ocean is smaller than A_dist, all
  !              of it is used.

  use precision, only : wp
  use constants, only : pi, R_earth
  use climber_grid, only : ni, nj, lon, lat, area

  implicit none

  type dist_type
    integer :: n = 0                        !! number of receiving cells
    integer, allocatable :: i(:), j(:)      !! indices of receiving cells
    real(wp), allocatable :: r(:)           !! fraction of the ocean area of each receiving cell used []
    real(wp) :: A_sum = 0._wp               !! total receiving ocean area [m2]
  end type dist_type

  private
  public :: dist_type, runoff_dist_init, dist_add

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  r u n o f f _ d i s t _ i n i t
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine runoff_dist_init(i_mode, n_cells, A_dist, f_ocn, coast_nbr, i_coast_nbr, j_coast_nbr, dist, label, verbose)

    implicit none

    integer, intent(in) :: i_mode
    integer, intent(in) :: n_cells          !! number of coastal cells, i_mode=1
    real(wp), intent(in) :: A_dist          !! receiving ocean area [m2], i_mode=2
    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(in) :: coast_nbr(:,:)
    integer, intent(in) :: i_coast_nbr(:,:,:)
    integer, intent(in) :: j_coast_nbr(:,:,:)
    type(dist_type), intent(inout) :: dist(:,:)
    character(len=*), intent(in) :: label
    logical, intent(in) :: verbose

    integer :: i, j, n, nc, ndest
    integer, allocatable :: ic(:), jc(:)
    real(wp), allocatable :: rc(:)
    real(wp) :: A_sum, A_min, A_max, A_mean, nc_mean

    allocate(ic(ni*nj), jc(ni*nj), rc(ni*nj))

    ndest = 0
    A_min = huge(1._wp)
    A_max = 0._wp
    A_mean = 0._wp
    nc_mean = 0._wp

    do j=1,nj
      do i=1,ni

        nc = 0
        A_sum = 0._wp

        if (coast_nbr(i,j).gt.0) then

          if (i_mode.eq.1) then
            nc = min(coast_nbr(i,j),n_cells)
            do n=1,nc
              ic(n) = i_coast_nbr(i,j,n)
              jc(n) = j_coast_nbr(i,j,n)
              rc(n) = 1._wp
              A_sum = A_sum + area(ic(n),jc(n))*f_ocn(ic(n),jc(n))  ! m2
            enddo
          else if (i_mode.eq.2) then
            call dist_area(i, j, A_dist, f_ocn, nc, ic, jc, rc, A_sum)
          else
            print *,'ERROR in runoff_dist_init: i_runoff_dist must be 1 or 2, got ',i_mode
            stop
          endif

          if (A_sum.le.0._wp) then
            print *,'WARNING runoff_dist_init (',trim(label),'): no receiving ocean for cell',i,j
            nc = 0
          else
            ndest = ndest+1
            A_min = min(A_min,A_sum)
            A_max = max(A_max,A_sum)
            A_mean = A_mean + A_sum
            nc_mean = nc_mean + real(nc,wp)
          endif

        endif

        if (allocated(dist(i,j)%i)) deallocate(dist(i,j)%i)
        if (allocated(dist(i,j)%j)) deallocate(dist(i,j)%j)
        if (allocated(dist(i,j)%r)) deallocate(dist(i,j)%r)
        allocate(dist(i,j)%i(nc), dist(i,j)%j(nc), dist(i,j)%r(nc))
        dist(i,j)%n = nc
        dist(i,j)%i = ic(1:nc)
        dist(i,j)%j = jc(1:nc)
        dist(i,j)%r = rc(1:nc)
        dist(i,j)%A_sum = A_sum

      enddo
    enddo

    if (verbose .and. ndest.gt.0) then
      print '(a,a,a,i2,a,i5,a,3es10.2,a,f6.1)', ' runoff_dist ',trim(label),': mode ',i_mode,', ',ndest, &
        ' destination cells, receiving area min/mean/max [m2]',A_min,A_mean/ndest,A_max, &
        ', mean number of cells',nc_mean/ndest
    endif

    deallocate(ic, jc, rc)

    return

  end subroutine runoff_dist_init


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  d i s t _ a r e a
  !   receiving cells covering ocean area A_dist around (i0,j0), ordered by
  !   path distance through the ocean (Dijkstra)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dist_area(i0, j0, A_dist, f_ocn, nc, ic, jc, rc, A_sum)

    implicit none

    integer, intent(in) :: i0, j0
    real(wp), intent(in) :: A_dist
    real(wp), intent(in) :: f_ocn(:,:)
    integer, intent(out) :: nc
    integer, intent(inout) :: ic(:), jc(:)
    real(wp), intent(inout) :: rc(:)
    real(wp), intent(out) :: A_sum

    real(wp), allocatable :: d(:,:)
    logical, allocatable :: done(:,:)
    integer, allocatable :: fi(:), fj(:)
    real(wp), allocatable :: fd(:)
    integer :: nf, k, kmin, i, j, ii, jj, di, dj
    real(wp) :: dmin, a, take, dnew
    logical :: ok

    allocate(d(ni,nj), done(ni,nj))
    allocate(fi(9*ni*nj), fj(9*ni*nj), fd(9*ni*nj))

    d = huge(1._wp)
    done = .false.
    nc = 0
    A_sum = 0._wp

    ! frontier (lazy priority queue, linear search: the frontier stays small)
    nf = 1
    fi(1) = i0
    fj(1) = j0
    fd(1) = 0._wp
    d(i0,j0) = 0._wp

    do while (nf.gt.0)

      ! pop the frontier cell closest to the destination
      kmin = 1
      dmin = fd(1)
      do k=2,nf
        if (fd(k).lt.dmin) then
          dmin = fd(k)
          kmin = k
        endif
      enddo
      i = fi(kmin)
      j = fj(kmin)
      fi(kmin) = fi(nf)
      fj(kmin) = fj(nf)
      fd(kmin) = fd(nf)
      nf = nf-1
      if (done(i,j)) cycle
      done(i,j) = .true.

      ! add the ocean area of this cell
      a = area(i,j)*f_ocn(i,j)
      if (a.gt.0._wp) then
        take = min(a, A_dist-A_sum)
        nc = nc+1
        ic(nc) = i
        jc(nc) = j
        rc(nc) = take/a
        A_sum = A_sum + take
        if (A_sum.ge.A_dist) exit
      endif

      ! relax the ocean neighbours
      do dj=-1,1
        jj = j+dj
        if (jj.lt.1 .or. jj.gt.nj) cycle
        do di=-1,1
          if (di.eq.0 .and. dj.eq.0) cycle
          ii = modulo(i+di-1,ni)+1
          if (f_ocn(ii,jj).le.0._wp .or. done(ii,jj)) cycle
          if (di.ne.0 .and. dj.ne.0) then
            ! diagonal step only through an ocean corner
            ok = f_ocn(ii,j).gt.0._wp .or. f_ocn(i,jj).gt.0._wp
            if (.not.ok) cycle
          endif
          dnew = d(i,j) + gc_dist(lon(i),lat(j),lon(ii),lat(jj))
          if (dnew.lt.d(ii,jj)) then
            d(ii,jj) = dnew
            nf = nf+1
            fi(nf) = ii
            fj(nf) = jj
            fd(nf) = dnew
          endif
        enddo
      enddo

    enddo

    deallocate(d, done, fi, fj, fd)

    return

  end subroutine dist_area


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  g c _ d i s t
  !   great-circle distance between two points given in degrees [m]
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function gc_dist(lon1, lat1, lon2, lat2) result(dist)

    implicit none

    real(wp), intent(in) :: lon1, lat1, lon2, lat2
    real(wp) :: dist

    real(wp) :: p1, p2, dp, dl, h

    p1 = lat1*pi/180._wp
    p2 = lat2*pi/180._wp
    dp = p2-p1
    dl = (lon2-lon1)*pi/180._wp
    h = sin(0.5_wp*dp)**2 + cos(p1)*cos(p2)*sin(0.5_wp*dl)**2
    dist = 2._wp*R_earth*asin(min(1._wp,sqrt(h)))

    return

  end function gc_dist


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  d i s t _ a d d
  !   spread flux F arriving at a destination cell over its receiving cells
  !   fld [F/m2] += F * r / A_sum
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine dist_add(fld, d, flux)

    implicit none

    real(wp), intent(inout) :: fld(:,:)
    type(dist_type), intent(in) :: d
    real(wp), intent(in) :: flux

    integer :: n

    if (flux.eq.0._wp) return

    do n=1,d%n
      fld(d%i(n),d%j(n)) = fld(d%i(n),d%j(n)) + flux*d%r(n)/d%A_sum
    enddo

    return

  end subroutine dist_add

end module runoff_dist_mod

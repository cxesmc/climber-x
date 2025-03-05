!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c o a s t _ c e l l s _ m o d
!
!  Purpose : find coastal cells
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
module coast_cells_mod

  use precision, only : wp
  use geo_params, only : n_coast_cells

  implicit none

  private
  public :: coast_cells

contains

  subroutine coast_cells(f_lnd, i_runoff, j_runoff, &
      mask_coast, mask_coast2, i_coast_nbr, j_coast_nbr, coast_nbr)

    implicit none

    real(wp), intent(in) :: f_lnd(:,:)
    integer, intent(in) :: i_runoff(:,:)
    integer, intent(in) :: j_runoff(:,:)
    integer, intent(out) :: mask_coast(:,:)
    integer, intent(out) :: mask_coast2(:,:)
    integer, intent(out) :: i_coast_nbr(:,:,:)
    integer, intent(out) :: j_coast_nbr(:,:,:)
    integer, intent(out) :: coast_nbr(:,:)

    integer :: i, j, ii, jj, iii, jjj, ni, nj, n, no, nbr, n_nbr_lnd, n_nbr_ocn, n_nbr_lnd2
    integer :: ip1, im1, jp1, jm1, ip2, im2, jp2, jm2, ip3, im3, ip4, im4
    integer, dimension(9) :: idx1, jdx1
    integer, dimension(28) :: idx2, jdx2
    integer, dimension(:,:), allocatable :: nbr_ocn

    ni = size(f_lnd,1)
    nj = size(f_lnd,2)

    allocate(nbr_ocn(ni,nj))

    ! find coastal cells
    do i=1,ni
      do j=1,nj
        mask_coast(i,j) = 0
        mask_coast2(i,j) = 0
        if (f_lnd(i,j).lt.1._wp) then
          n_nbr_lnd = 0
          n_nbr_lnd2 = 0
          n_nbr_ocn = 0
          ! look for land in 9-cell neighborhood
          do ii=i-1,i+1
            do jj=j-1,j+1
              iii = ii
              if (iii.eq.0) iii = ni
              if (iii.eq.ni+1) iii = 1
              jjj = jj
              jjj = max(1,jjj)
              jjj = min(nj,jjj)
              if (f_lnd(iii,jjj).gt.0._wp) then
                n_nbr_lnd = n_nbr_lnd+1
              endif
              if (f_lnd(iii,jjj).eq.1._wp) then
                n_nbr_lnd2 = n_nbr_lnd2+1
              endif
              if (f_lnd(iii,jjj).lt.1._wp) then
                n_nbr_ocn = n_nbr_ocn+1
              endif
            enddo
          enddo
          nbr_ocn(i,j) = n_nbr_ocn
          if (n_nbr_lnd>0 .and. n_nbr_ocn<9) then
            ! coastal cell 
            mask_coast(i,j) = 1
          endif
          if (n_nbr_lnd2>0) then
            ! coastal cell 
            mask_coast2(i,j) = 1
          endif
        endif
      enddo
    enddo

    ! make sure runoff destination cells are part of coastal mask
    do i=1,ni
      do j=1,nj
        if (i_runoff(i,j).ne.0 .and. j_runoff(i,j).ne.0) then
          mask_coast(i_runoff(i,j),j_runoff(i,j)) = 1
        endif
      enddo
    enddo

    ! index neighbors of each coastal cell
    do i=1,ni
      do j=1,nj
        if (mask_coast(i,j).eq.1) then
          ! find n_coast_cells coastal cells in 25-cell neighborhood
          ip1 = i+1
          if (ip1.eq.ni+1) ip1 = 1
          im1 = i-1
          if (im1.eq.0) im1 = ni
          jp1 = j+1
          if (jp1.eq.nj+1) jp1 = nj
          jm1 = j-1
          if (jm1.eq.0) jm1 = 1
          ip2 = i+2
          if (ip2.gt.ni) ip2 = ip2-ni
          im2 = i-2
          if (im2.lt.1) im2 = ni+im2
          ip3 = i+3
          if (ip3.gt.ni) ip3 = ip3-ni
          im3 = i-3
          if (im3.lt.1) im3 = ni+im3
          ip4 = i+4
          if (ip4.gt.ni) ip4 = ip4-ni
          im4 = i-4
          if (im4.lt.1) im4 = ni+im4
          jp2 = j+2
          if (jp2.gt.nj) jp2 = nj
          jm2 = j-2
          if (jm2.lt.1) jm2 = 1

          idx1 = (/i, i,   i,   ip1, im1, im1, im1, ip1, ip1/)
          jdx1 = (/j, jp1, jm1, j,   j,   jm1, jp1, jm1, jp1/)
          idx2 = (/i,   i,   ip2, im2, im2, im2, ip2, ip2, im1, ip1, im1, ip1, im2, im2, ip2, ip2, &
            ip3, im3, ip4, im4, ip3, im3, ip4, im4, ip3, im3, ip4, im4/)
          jdx2 = (/jp2, jm2, j,   j,   jm1, jp1, jm1, jp1, jm2, jm2, jp2, jp2, jm2, jp2, jm2, jp2, &
            j,   j,   j,   j,   jm1, jm1, jm1, jm1, jp1, jp1, jp1, jp1/)

          ! first add the grid cell itself
          nbr = 1
          i_coast_nbr(i,j,nbr) = i
          j_coast_nbr(i,j,nbr) = j
          ! then look for neighbors in 3x3 neighbourhood
          ! add first neighbors that have low number of ocean neighbors
          do no=1,8
            n = 1
            do while (nbr.lt.n_coast_cells .and. n.lt.size(idx1))
              n = n+1
              if (nbr_ocn(idx1(n),jdx1(n)).eq.no .and. mask_coast(idx1(n),jdx1(n)).eq.1) then
                nbr = nbr+1
                i_coast_nbr(i,j,nbr) = idx1(n)
                j_coast_nbr(i,j,nbr) = jdx1(n)
              endif
            enddo
          enddo
          ! then look into more distant neighbors, if needed
          n = 0
          do while (nbr.lt.n_coast_cells .and. n.lt.size(idx2))
            n = n+1
            if (mask_coast(idx2(n),jdx2(n)).eq.1) then
              nbr = nbr+1
              i_coast_nbr(i,j,nbr) = idx2(n)
              j_coast_nbr(i,j,nbr) = jdx2(n)
            endif
          enddo

          coast_nbr(i,j) = nbr
        else
          i_coast_nbr(i,j,:) = 0 
          j_coast_nbr(i,j,:) = 0
          coast_nbr(i,j) = 0
        endif
      enddo
    enddo


    return

  end subroutine coast_cells

end module coast_cells_mod


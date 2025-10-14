!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f r e e _ s u r f a c e _ m o d
!
!  Purpose : diagnose elevation of the free surface
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
module free_surface_mod

  use precision, only : wp
  use constants, only : g
  use ocn_params, only : rho0
  use ocn_grid, only : maxi, maxj, maxk, k1, dz, zw, zro, mask_ocn, mask_c

  implicit none

  real(wp), parameter :: depth_ref = 1500._wp

  private
  public :: free_surface

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f r e e _ s u r f a c e
  !   Purpose    :  diagnose elevation of the free surface 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine free_surface(rho, ssh, error)

    implicit none

    real(wp), intent(in) :: rho(:,:,:)
    real(wp), intent(out) :: ssh(:,:)
    logical, intent(inout) :: error

    real(wp), allocatable :: rho_tmp(:,:,:)
    real(wp), allocatable :: rho_tmp1(:,:,:)
    real(wp), allocatable :: ssh_tmp(:,:)

    integer :: i, j, n, k, k_ref, nbr, ii, jj, iii, jjj
    real(wp) :: rhotmp, rho_sum, zs_sum
    logical :: flag
    real(wp) :: p_ref, z_ref, integral
    
    real(wp), parameter :: misval = -999._wp


    ! local copy of rho which will be filled
    allocate(rho_tmp,source=rho)
    where (mask_c.eq.0) rho_tmp = misval  ! missing value
    allocate(rho_tmp1,source=rho_tmp)

    ! find level closest to reference depth
    k_ref = minloc(abs(-zw-depth_ref),1)
    z_ref = -zw(k_ref-1)

    ! reference pressure at depth z_ref
    p_ref = rho0*g*z_ref

    flag = .true.

    n = 0

    loop : do while (flag)

      n=n+1
      flag = .false.

      ! fill rho values in regions shallower than z_ref
      do j=1,maxj
        do i=1,maxi

          if (mask_ocn(i,j).eq.1 .and. k1(i,j).gt.k_ref) then ! fixme, check

            do k=k_ref,k1(i,j)-1
              if (rho_tmp(i,j,k).eq.misval) then  
                if (n.gt.1e5) then
                  print *
                  print *,'ERROR: stuck in ocean free_surface loop'
                  print *,'i,j,k',i,j,k
                  print *,'k1,k_ref',k1(i,j),k_ref
                  print *,'mask_ocn',mask_ocn(i,j)
                  print *,'rho',rho_tmp(i,j,k)
                  print *
                  print *,'mask'
                  print '(5i4)', mask_ocn(max(1,i-2):min(maxi,i+2),min(maxj,j+2))
                  print '(5i4)', mask_ocn(max(1,i-2):min(maxi,i+2),min(maxj,j+1))
                  print '(5i4)', mask_ocn(max(1,i-2):min(maxi,i+2),j)
                  print '(5i4)', mask_ocn(max(1,i-2):min(maxi,i+2),max(1,j-1))
                  print '(5i4)', mask_ocn(max(1,i-2):min(maxi,i+2),max(1,j-2))
                  print *
                  print *,'k1'
                  print '(5i4)', k1(max(1,i-2):min(maxi,i+2),min(maxj,j+2))
                  print '(5i4)', k1(max(1,i-2):min(maxi,i+2),min(maxj,j+1))
                  print '(5i4)', k1(max(1,i-2):min(maxi,i+2),j)
                  print '(5i4)', k1(max(1,i-2):min(maxi,i+2),max(1,j-1))
                  print '(5i4)', k1(max(1,i-2):min(maxi,i+2),max(1,j-2))
                  print *
                  print *,'rho'
                  print '(5F7.1)', rho_tmp(max(1,i-2):min(maxi,i+2),min(maxj,j+2),k)
                  print '(5F7.1)', rho_tmp(max(1,i-2):min(maxi,i+2),min(maxj,j+1),k)
                  print '(5F7.1)', rho_tmp(max(1,i-2):min(maxi,i+2),j,k)
                  print '(5F7.1)', rho_tmp(max(1,i-2):min(maxi,i+2),max(1,j-1),k)
                  print '(5F7.1)', rho_tmp(max(1,i-2):min(maxi,i+2),max(1,j-2),k)
                  !stop
                  error = .true.
                  exit loop
                endif
                flag = .true.
                ! mean density from neighbor wet cells at same level
                nbr = 0
                rho_sum = 0._wp
                do ii=i-1,i+1
                  do jj=j-1,j+1
                    iii = ii
                    if (iii.eq.0) iii = maxi
                    if (iii.eq.maxi+1) iii = 1
                    jjj = jj
                    jjj = max(1,jjj)
                    jjj = min(maxj,jjj)
                    rhotmp = rho_tmp(iii,jjj,k)
                    if (rhotmp.ne.misval) then
                      nbr = nbr+1
                      rho_sum = rho_sum + rhotmp
                    endif
                  enddo
                enddo
                if (nbr.gt.0) then
                  rho_tmp1(i,j,k) = rho_sum/nbr
                endif
                if (n.eq.100) then
                  ! probably isolated ocean points, fill with global mean density 
                  nbr = 0
                  rho_sum = 0._wp
                  do ii=1,maxi
                    do jj=1,maxj
                      rhotmp = rho_tmp(ii,jj,k)
                      if (rhotmp.ne.misval) then
                        nbr = nbr+1
                        rho_sum = rho_sum + rhotmp
                      endif
                    enddo
                  enddo
                  rho_tmp1(i,j,k) = rho_sum/nbr
                endif
              endif
            enddo

          endif

        enddo
      enddo

      rho_tmp = rho_tmp1

    enddo loop

    ! derive free surface elevation
    do j=1,maxj
      do i=1,maxi

        if (mask_ocn(i,j).eq.1) then

          integral = 0._wp
          do k=k_ref,maxk
            integral = integral + g*(rho_tmp(i,j,k)-rho0)*dz(k)
          enddo

          ssh(i,j) = (p_ref-integral)/(rho0*g) - z_ref

        else

          ssh(i,j) = misval

        endif
      enddo
    enddo

    allocate(ssh_tmp,source=ssh)

    ! extrapolate over land neighbors (needed for gradient computation)
    do j=1,maxj
      do i=1,maxi

        if (mask_ocn(i,j).eq.0) then
          nbr = 0
          zs_sum = 0._wp
          do ii=i-1,i+1
            do jj=j-1,j+1
              iii = ii
              if (iii.eq.0) iii = maxi
              if (iii.eq.maxi+1) iii = 1
              jjj = jj
              jjj = max(1,jjj)
              jjj = min(maxj,jjj)
              if (ssh_tmp(iii,jjj).ne.misval) then
                nbr = nbr+1
                zs_sum = zs_sum + ssh_tmp(iii,jjj)
              endif
            enddo
          enddo
          if (nbr.gt.0) then
            ssh(i,j) = zs_sum/nbr 
          endif
        endif

      enddo
    enddo

    deallocate(rho_tmp)
    deallocate(rho_tmp1)
    deallocate(ssh_tmp)

    return
  end subroutine free_surface

end module free_surface_mod

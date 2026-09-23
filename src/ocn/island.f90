!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i s l a n d _ m o d
!
!  Purpose : calculate path integral around islands
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
module island_mod

  use precision, only : wp
  use ocn_grid, only : maxi, maxj, maxk, c, dphi, rcv, dsv, dz, rh, R_earth, ku
  use ocn_grid, only : npi, lpisl, ipisl, jpisl
  use ocn_params, only : fcor, fcorv, drag, rho0, i_cor_form, qcor
  !$use om_lib

  implicit none

  private
  public :: island

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  i s l a n d
  !   Purpose    :  calculate path integral around island
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine island(ubloc,tau,bp, &
                   isl,indj, &
                   erisl1, psiloc)

    implicit none

    real(wp), dimension(:,0:,0:), intent(in) :: ubloc
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), intent(in) :: bp(:,:,:)
    integer, intent(in) :: isl, indj
  
    real(wp), intent(out) :: erisl1
    ! streamfunction of the same solution as ubloc, required by the J3 flux form (i_cor_form=1), whose Coriolis
    ! integrand is psi itself rather than the velocity
    real(wp), dimension(0:,0:), intent(in), optional :: psiloc

    integer :: i, k, lpi, ipi, jpi, im1
    real(wp) :: cor, tv1, tv2


    erisl1 = 0._wp
    !!$omp parallel do private (i,k,lpi,ipi,jpi,cor,tv1,tv2) reduction(+:erisl1)
    do i=1,npi(isl)
       lpi = lpisl(i,isl)
       ipi = ipisl(i,isl)
       jpi = jpisl(i,isl)
       if (abs(lpi).eq.1) then
          !cor = - fcor(jpi)*0.25_wp*(ubloc(2,ipi,jpi) + ubloc(2,ipi+1,jpi) &  ! m/s2
          !    + ubloc(2,ipi,jpi-1) + ubloc(2,ipi+1,jpi-1))
          cor = - 0.5_wp * (fcorv(jpi)*0.5_wp*(ubloc(2,ipi,jpi) + ubloc(2,ipi+1,jpi))  &  ! m/s2
              + fcorv(jpi-1)*0.5_wp*(ubloc(2,ipi,jpi-1) + ubloc(2,ipi+1,jpi-1)))
       else
          !cor = fcorv(jpi)*0.25_wp*(ubloc(1,ipi-1,jpi) + ubloc(1,ipi,jpi) &
          !    + ubloc(1,ipi-1,jpi+1) + ubloc(1,ipi,jpi+1))
          cor = 0.5_wp * (fcor(jpi)*0.5_wp*(ubloc(1,ipi-1,jpi) + ubloc(1,ipi,jpi)) &
              + fcor(jpi+1)*0.5_wp*(ubloc(1,ipi-1,jpi+1) + ubloc(1,ipi,jpi+1)))
       endif

       if (i_cor_form.eq.1) then
          ! J3 flux form: the Coriolis flux through a face is 0.5*(psi at the two ends of the face)*(q across the face),
          ! with the metric factors of the segment length cancelling, so it is added here already multiplied by its length
          ! (and divided below by the length that the common factor applies)
          cor = 0._wp
       endif

       if (jpi.lt.maxj) then ! to avoid acessing rcv(maxj) and dsv(maxj), out of bound!
          erisl1 = erisl1 + sign(1,lpi)*(drag(abs(lpi),ipi,jpi) &  ! m/s2
              * ubloc(abs(lpi),ipi,jpi) + cor &
              - indj*tau(abs(lpi),ipi,jpi)/rho0*rh(abs(lpi),ipi,jpi)) &
              * (c(jpi)*dphi*(2._wp - abs(lpi)) &
              + rcv(jpi)*dsv(jpi)*(abs(lpi) - 1._wp))  
       else
          erisl1 = erisl1 + sign(1,lpi)*(drag(abs(lpi),ipi,jpi) &  ! m/s2
              * ubloc(abs(lpi),ipi,jpi) + cor &
              - indj*tau(abs(lpi),ipi,jpi)/rho0*rh(abs(lpi),ipi,jpi)) &
              * (c(jpi)*dphi*(2._wp - abs(lpi)))
       endif

       if (i_cor_form.eq.1) then
          ! flux form J3: sign * 0.5*(psi at the two ends of the face) * (q difference across the face) / (rho0 R).
          ! Its area sum over the enclosed psi points is exactly the J3 Coriolis term of stencil_coef, so the value of
          ! the path integral does not depend on the path taken around the island
          im1 = ipi-1
          if (im1.lt.0) im1 = im1+maxi
          if (abs(lpi).eq.1) then
             ! u-face (1,ipi,jpi): spans psi(ipi,jpi-1)..psi(ipi,jpi), separates cells (ipi,jpi) and (ipi+1,jpi)
             erisl1 = erisl1 + sign(1,lpi)*0.5_wp*(psiloc(ipi,jpi-1) + psiloc(ipi,jpi)) &
                    * (qcor(ipi+1,jpi) - qcor(ipi,jpi))/(rho0*R_earth)
          else
             ! v-face (2,ipi,jpi): spans psi(ipi-1,jpi)..psi(ipi,jpi), separates cells (ipi,jpi) and (ipi,jpi+1)
             erisl1 = erisl1 + sign(1,lpi)*0.5_wp*(psiloc(im1,jpi) + psiloc(ipi,jpi)) &
                    * (qcor(ipi,jpi+1) - qcor(ipi,jpi))/(rho0*R_earth)
          endif
       endif

       ! calc tricky bits and add to source term for path integral round 
       ! islands, all sums have at least one element
     
       if (indj.eq.1) then
          if (abs(lpi).eq.1) then
             tv1 = 0._wp
             do k=ku(1,ipi,jpi),maxk
                tv1 = tv1 + bp(ipi+1,jpi,k)*dz(k)  ! kg/s2
             enddo
             tv2 = 0._wp
             do k=ku(1,ipi,jpi),maxk
                tv2 = tv2 + bp(ipi,jpi,k)*dz(k)
             enddo
             erisl1 = erisl1 + (tv1-tv2) &  ! m/s2
                    * sign(1,lpi)*rh(1,ipi,jpi)/(rho0*R_earth)
          else
             tv1 = 0._wp
             do k=ku(2,ipi,jpi),maxk
                tv1 = tv1 + bp(ipi,jpi+1,k)*dz(k)
             enddo
             tv2 = 0._wp
             do k=ku(2,ipi,jpi),maxk
                tv2 = tv2 + bp(ipi,jpi,k)*dz(k)
             enddo
             erisl1 = erisl1 + (tv1 - tv2) &
                    * sign(1,lpi)*rh(2,ipi,jpi)/(rho0*R_earth)
          endif
       endif
    enddo
    !!$omp end parallel do

   return

  end subroutine island

end module island_mod

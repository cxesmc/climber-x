!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : a d i f a _ m o d
!
!  Purpose : compute advective and diffusive fluxes of energy, 
!            water and dust
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Andrey Ganopolski and Matteo Willeit
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
module adifa_mod

  use atm_params, only : wp
  use atm_params, only : cp
  use atm_grid, only : im, imc, jm, jmc, km
  use atm_grid, only : dplx, dply, dy, dxt, dxu, sqr
  !$ use omp_lib

  implicit none
  
  private
  public :: adifa
  
contains
    
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a d i f a
  !   Purpose    :  compute advective and diffusive fluxes of energy, 
  !              :  water and dust
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine adifa(fax, fay, tp, q3, qam, d3, cam, diffxdse, diffydse, diffxwtr, diffywtr, diffxdst, diffydst, &
    convdse, convwtr, convdst, convco2, faxdse, faxwtr, faxdst, faxco2, faydse, faywtr, faydst, fayco2, &
    fdxdse, fdxwtr, fdxdst, fdxco2, fdydse, fdywtr, fdydst, fdyco2) 

    implicit none

    real(wp), intent(in   ) :: fax(:,:,:)
    real(wp), intent(in   ) :: fay(:,:,:)
    real(wp), intent(in   ) :: tp(:,:,:)
    real(wp), intent(in   ) :: q3(:,:,:)
    real(wp), intent(in   ) :: qam(:,:,:)
    real(wp), intent(in   ) :: d3(:,:,:)
    real(wp), intent(in   ) :: cam(:,:)
    real(wp), intent(in   ) :: diffxdse(:,:)
    real(wp), intent(in   ) :: diffydse(:,:)
    real(wp), intent(in   ) :: diffxwtr(:,:)
    real(wp), intent(in   ) :: diffywtr(:,:)
    real(wp), intent(in   ) :: diffxdst(:,:)
    real(wp), intent(in   ) :: diffydst(:,:)
    
    real(wp), intent(inout) :: convdse(:,:)
    real(wp), intent(inout) :: convwtr(:,:,:)
    real(wp), intent(inout) :: convdst(:,:)
    real(wp), intent(inout) :: convco2(:,:)

    real(wp), intent(out  ) :: faxdse(:,:)
    real(wp), intent(out  ) :: faxwtr(:,:,:)
    real(wp), intent(out  ) :: faxdst(:,:)
    real(wp), intent(out  ) :: faxco2(:,:)
    real(wp), intent(out  ) :: faydse(:,:)
    real(wp), intent(out  ) :: faywtr(:,:,:)
    real(wp), intent(out  ) :: faydst(:,:)
    real(wp), intent(out  ) :: fayco2(:,:)
    real(wp), intent(out  ) :: fdxdse(:,:)
    real(wp), intent(out  ) :: fdxwtr(:,:,:)
    real(wp), intent(out  ) :: fdxdst(:,:)
    real(wp), intent(out  ) :: fdxco2(:,:)
    real(wp), intent(out  ) :: fdydse(:,:)
    real(wp), intent(out  ) :: fdywtr(:,:,:)
    real(wp), intent(out  ) :: fdydst(:,:)
    real(wp), intent(out  ) :: fdyco2(:,:)

    integer :: i, j, k, imi, jmi
    real(wp) :: tpup, qup, dup, cup
    real(wp) :: dpl_x, dpl_y
    real(wp) :: tp_ijk, tp_i1jk, tp_ij1k
    real(wp) :: q3_ijk, q3_i1jk, q3_ij1k
    real(wp) :: d3_ijk, d3_i1jk, d3_ij1k
    real(wp) :: c3_ij, c3_i1j, c3_ij1
    real(wp) :: fax_ijk, fay_ijk
    real(wp) :: fro16_fwtr, fro18_fwtr


    !$omp parallel do private(i, j, k, imi, jmi, tpup, qup, dup, cup, dpl_x, dpl_y) &
    !$omp private (tp_ijk, tp_i1jk, tp_ij1k, q3_ijk, q3_i1jk, q3_ij1k, d3_ijk, d3_i1jk, d3_ij1k, c3_ij, c3_i1j, c3_ij1, fax_ijk, fay_ijk, fro16_fwtr, fro18_fwtr)
    do j=1,jm

      jmi=max(1,j-1)

      do i=1,im

        imi=i-1
        if (imi.eq.0) imi=im

        ! initialize vertically integrated fluxes
        faxdse(i,j)   = 0._wp      
        faxwtr(i,j,:) = 0._wp
        faxdst(i,j)   = 0._wp
        faxco2(i,j)   = 0._wp

        faydse(i,j)   = 0._wp      
        faywtr(i,j,:) = 0._wp
        faydst(i,j)   = 0._wp
        fayco2(i,j)   = 0._wp

        fdxdse(i,j)   = 0._wp      
        fdxwtr(i,j,:) = 0._wp       
        fdxdst(i,j)   = 0._wp       
        fdxco2(i,j)   = 0._wp       

        fdydse(i,j)   = 0._wp      
        fdywtr(i,j,:) = 0._wp
        fdydst(i,j)   = 0._wp 
        fdyco2(i,j)   = 0._wp 

        c3_ij  = cam(i,j)
        c3_i1j = cam(imi,j)
        c3_ij1 = cam(i,jmi)

        ! integrate fluxes vertically
        do k=1,km

          ! for efficiency, to avoid accessing same element of 3d arrays several times
          tp_ijk  = tp(i,j,k)
          tp_i1jk = tp(imi,j,k)
          tp_ij1k = tp(i,jmi,k)
          q3_ijk  = q3(i,j,k)
          q3_i1jk = q3(imi,j,k)
          q3_ij1k = q3(i,jmi,k)
          d3_ijk  = d3(i,j,k)
          d3_i1jk = d3(imi,j,k)
          d3_ij1k = d3(i,jmi,k)

          fax_ijk = fax(i,j,k)
          fay_ijk = fay(i,j,k)

          !-----------------------------------
          ! advective fluxes
          !-----------------------------------

          !-----------------------------------
          ! zonal components

          ! Upstream values
          if (fax_ijk.gt.0._wp) then
            tpup = tp_i1jk
            qup  = q3_i1jk
            dup  = d3_i1jk
            cup  = c3_i1j
          else    
            tpup = tp_ijk
            qup  = q3_ijk
            dup  = d3_ijk
            cup  = c3_ij
          endif
          faxdse(i,j)   = faxdse(i,j)   + fax_ijk*tpup ! kg/s * K
          faxwtr(i,j,1) = faxwtr(i,j,1) + fax_ijk*qup  ! kg/s * kg/kg
          faxdst(i,j)   = faxdst(i,j)   + fax_ijk*dup  
          faxco2(i,j)   = faxco2(i,j)   + fax_ijk*cup  ! kg/s * kgCO2/kg = kgCO2/s

          !-----------------------------------
          ! meridional components

          ! Upstream values
          if (fay_ijk.gt.0._wp) then
            tpup = tp_ijk
            qup  = q3_ijk
            dup  = d3_ijk
            cup  = c3_ij
          else    
            tpup = tp_ij1k
            qup  = q3_ij1k
            dup  = d3_ij1k
            cup  = c3_ij1
          endif 
          faydse(i,j)   = faydse(i,j)   + fay_ijk*tpup ! kg/s * K
          faywtr(i,j,1) = faywtr(i,j,1) + fay_ijk*qup  ! kg/s * kg/kg
          faydst(i,j)   = faydst(i,j)   + fay_ijk*dup
          fayco2(i,j)   = fayco2(i,j)   + fay_ijk*cup

          !-----------------------------------
          ! diffusive fluxes
          !-----------------------------------

          if (k.le.km-2) then   ! limit to troposphere

            !-----------------------------------
            ! zonal diffusive fluxes
            dpl_x = dplx(i,j,k)
            fdxdse(i,j)   = fdxdse(i,j)   + diffxdse(i,j)*dy*dpl_x*(tp_i1jk-tp_ijk)/dxt(j) ! m2/s * K * kg/m2 = kg/s * K
            fdxwtr(i,j,1) = fdxwtr(i,j,1) + diffxwtr(i,j)*dy*dpl_x*(q3_i1jk-q3_ijk)/dxt(j) 
            fdxdst(i,j)   = fdxdst(i,j)   + diffxdst(i,j)*dy*dpl_x*(d3_i1jk-d3_ijk)/dxt(j)
            fdxco2(i,j)   = fdxco2(i,j)   + diffxdst(i,j)*dy*dpl_x*(c3_i1j-c3_ij)/dxt(j)

            !-----------------------------------
            ! meridional diffusive fluxes
            dpl_y = dply(i,j,k)
            fdydse(i,j)   = fdydse(i,j)   + diffydse(i,j)*dxu(j)*dpl_y*(tp_ijk-tp_ij1k)/dy
            fdywtr(i,j,1) = fdywtr(i,j,1) + diffywtr(i,j)*dxu(j)*dpl_y*(q3_ijk-q3_ij1k)/dy
            fdydst(i,j)   = fdydst(i,j)   + diffydst(i,j)*dxu(j)*dpl_y*(d3_ijk-d3_ij1k)/dy
            fdyco2(i,j)   = fdyco2(i,j)   + diffydst(i,j)*dxu(j)*dpl_y*(c3_ij-c3_ij1)/dy

          endif

        enddo

        !-------------------------------------
        ! oxygen isotope fractions in vertically integrated water fluxes
        !-------------------------------------

        ! isotope fractions in near-surface specific humidity
        fro16_fwtr = 0.1_wp*qam(i,j,2)/qam(i,j,1)  ! qam(i,j,2) is in g/m2/day, fro16 in percent
        fro18_fwtr = 0.1_wp*qam(i,j,3)/qam(i,j,1)  ! ditto
        ! advective fluxes
        faxwtr(i,j,2) = faxwtr(i,j,1)*10._wp*fro16_fwtr  ! factor 1000 between water/isotope tracers
        faxwtr(i,j,3) = faxwtr(i,j,1)*10._wp*fro18_fwtr
        faywtr(i,j,2) = faywtr(i,j,1)*10._wp*fro16_fwtr
        faywtr(i,j,3) = faywtr(i,j,1)*10._wp*fro18_fwtr
        ! diffusive fluxes
        fdxwtr(i,j,2) = fdxwtr(i,j,1)*10._wp*fro16_fwtr  ! ditto
        fdxwtr(i,j,3) = fdxwtr(i,j,1)*10._wp*fro18_fwtr
        fdywtr(i,j,2) = fdywtr(i,j,1)*10._wp*fro16_fwtr
        fdywtr(i,j,3) = fdywtr(i,j,1)*10._wp*fro18_fwtr

      enddo

      ! no-flux condition at the poles
      if (j.eq.jm) then
        faydse(:,jmc) = 0._wp
        faywtr(:,jmc,1) = 0._wp
        faywtr(:,jmc,2) = 0._wp
        faywtr(:,jmc,3) = 0._wp
        faydst(:,jmc) = 0._wp
        fayco2(:,jmc) = 0._wp
        fdydse(:,jmc) = 0._wp
        fdywtr(:,jmc,1) = 0._wp
        fdywtr(:,jmc,2) = 0._wp
        fdywtr(:,jmc,3) = 0._wp
        fdydst(:,jmc) = 0._wp
        fdyco2(:,jmc) = 0._wp
      endif

      ! Cycling
      faxdse(imc,j) = faxdse(1,j)
      faxwtr(imc,j,1) = faxwtr(1,j,1)
      faxwtr(imc,j,2) = faxwtr(1,j,2)
      faxwtr(imc,j,3) = faxwtr(1,j,3)
      faxdst(imc,j) = faxdst(1,j)
      faxco2(imc,j) = faxco2(1,j)
      fdxdse(imc,j) = fdxdse(1,j)
      fdxwtr(imc,j,1) = fdxwtr(1,j,1)
      fdxwtr(imc,j,2) = fdxwtr(1,j,2)
      fdxwtr(imc,j,3) = fdxwtr(1,j,3)
      fdxdst(imc,j) = fdxdst(1,j)
      fdxco2(imc,j) = fdxco2(1,j)

    enddo
    !$omp end parallel do


    !-----------------------------------
    ! fluxes convergency
    !-----------------------------------

    !$omp parallel do collapse(2) private(i,j)
    do j=1,jm
      do i=1,im 

        !-----------------------------------
        ! dry static energy
        convdse(i,j)= &
                       (faxdse(i,j  )-faxdse(i+1,j) &  
                       +faydse(i,j+1)-faydse(i  ,j) &
                       +fdxdse(i,j  )-fdxdse(i+1,j) &
                       +fdydse(i,j+1)-fdydse(i  ,j)) &
                       /sqr(i,j) * cp  ! K * kg/s / m2 * J/kg/K = J/m2/s = W/m2

        !-----------------------------------
        ! water
        convwtr(i,j,:)= 0.9_wp*convwtr(i,j,:) + 0.1_wp * &  ! relax in time
                         (faxwtr(i,j  ,:)-faxwtr(i+1,j,:) &
                         +faywtr(i,j+1,:)-faywtr(i  ,j,:) &
                         +fdxwtr(i,j  ,:)-fdxwtr(i+1,j,:) &
                         +fdywtr(i,j+1,:)-fdywtr(i  ,j,:)) &
                         /sqr(i,j)       ! kg/kg * kg/s / m2 = kg/m2/s 

        !-----------------------------------
        ! dust
        convdst(i,j)= &
                       (faxdst(i,j  )-faxdst(i+1,j) &
                       +faydst(i,j+1)-faydst(i  ,j) &
                       +fdxdst(i,j  )-fdxdst(i+1,j) &
                       +fdydst(i,j+1)-fdydst(i  ,j)) &
                       /sqr(i,j)

        !-----------------------------------
        ! carbon
        convco2(i,j)= &
                       (faxco2(i,j  )-faxco2(i+1,j) &
                       +fayco2(i,j+1)-fayco2(i  ,j) &
                       +fdxco2(i,j  )-fdxco2(i+1,j) &
                       +fdyco2(i,j+1)-fdyco2(i  ,j)) &
                       /sqr(i,j)        ! kgCO2/s/m2

      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine adifa

end module adifa_mod

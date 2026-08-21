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
  use atm_params, only : tstep
  use atm_params, only : l_diff_impl
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
  subroutine adifa(fax, fay, fax_psi, fay_psi, fax_psi_topo, fay_psi_topo, &
    tp, q3, d3, cam, diffxdse, diffydse, diffxwtr, diffywtr, diffxdst, diffydst, &
    convdse, convwtr_adv, convwtr_dif, convdst, convco2, faxdse, faxwtr, faxdst, faxco2, faydse, faywtr, faydst, fayco2, &
    fdxdse, fdxwtr, fdxdst, fdxco2, fdydse, fdywtr, fdydst, fdyco2, &
    convdse_psi, convdse_psi_topo)

    implicit none

    real(wp), intent(in   ) :: fax(:,:,:)
    real(wp), intent(in   ) :: fay(:,:,:)
    ! the column mass correction, split into its grad(psi) and grad(psi_topo) parts
    real(wp), intent(in   ) :: fax_psi(:,:,:)
    real(wp), intent(in   ) :: fay_psi(:,:,:)
    real(wp), intent(in   ) :: fax_psi_topo(:,:,:)
    real(wp), intent(in   ) :: fay_psi_topo(:,:,:)
    real(wp), intent(in   ) :: tp(:,:,:)
    real(wp), intent(in   ) :: q3(:,:,:)
    real(wp), intent(in   ) :: d3(:,:,:)
    real(wp), intent(in   ) :: cam(:,:)
    real(wp), intent(in   ) :: diffxdse(:,:)
    real(wp), intent(in   ) :: diffydse(:,:)
    real(wp), intent(in   ) :: diffxwtr(:,:)
    real(wp), intent(in   ) :: diffywtr(:,:)
    real(wp), intent(in   ) :: diffxdst(:,:)
    real(wp), intent(in   ) :: diffydst(:,:)
    
    real(wp), intent(inout) :: convdse(:,:)
    real(wp), intent(inout) :: convwtr_adv(:,:)   ! advective moisture convergence
    real(wp), intent(  out) :: convwtr_dif(:,:)   ! explicit diffusive moisture convergence (0 when l_diff_impl, then set by diffuse_impl)
    real(wp), intent(inout) :: convdst(:,:)
    real(wp), intent(inout) :: convco2(:,:)

    real(wp), intent(out  ) :: faxdse(:,:)
    real(wp), intent(out  ) :: faxwtr(:,:)
    real(wp), intent(out  ) :: faxdst(:,:)
    real(wp), intent(out  ) :: faxco2(:,:)
    real(wp), intent(out  ) :: faydse(:,:)
    real(wp), intent(out  ) :: faywtr(:,:)
    real(wp), intent(out  ) :: faydst(:,:)
    real(wp), intent(out  ) :: fayco2(:,:)
    real(wp), intent(out  ) :: fdxdse(:,:)
    real(wp), intent(out  ) :: fdxwtr(:,:)
    real(wp), intent(out  ) :: fdxdst(:,:)
    real(wp), intent(out  ) :: fdxco2(:,:)
    real(wp), intent(out  ) :: fdydse(:,:)
    real(wp), intent(out  ) :: fdywtr(:,:)
    real(wp), intent(out  ) :: fdydst(:,:)
    real(wp), intent(out  ) :: fdyco2(:,:)
    ! advective DSE convergence carried by each part of the mass correction, W/m2.
    ! Formed with the same upstream values as convdse, so the three add up exactly:
    ! the advective part of convdse is the sum of these two and the uncorrected flux.
    real(wp), intent(out  ) :: convdse_psi(:,:)
    real(wp), intent(out  ) :: convdse_psi_topo(:,:)

    integer :: i, j, k, imi, jmi
    real(wp) :: tpup, qup, dup, cup
    real(wp) :: dpl_x, dpl_y
    real(wp) :: tp_ijk, tp_i1jk, tp_ij1k
    real(wp) :: q3_ijk, q3_i1jk, q3_ij1k
    real(wp) :: d3_ijk, d3_i1jk, d3_ij1k
    real(wp) :: c3_ij, c3_i1j, c3_ij1
    real(wp) :: fax_ijk, fay_ijk
    real(wp) :: fdivdse, fdivwtr, fdivdst, fdivco2
    real(wp) :: cmass, ctp, mc_psi, mc_psi_topo, m_k
    integer :: ipl
    real(wp), dimension(imc,jm) :: faxdse_psi, faxdse_psi_topo
    ! tropospheric mass weighted mean of tp, the reference the parts are measured against
    real(wp), dimension(im,jm) :: tp_col
    real(wp), dimension(im,jmc) :: faydse_psi, faydse_psi_topo


    !$omp parallel do private(i, j, k, imi, jmi, tpup, qup, dup, cup, dpl_x, dpl_y) &
    !$omp private (tp_ijk, tp_i1jk, tp_ij1k, q3_ijk, q3_i1jk, q3_ij1k, d3_ijk, d3_i1jk, d3_ij1k, c3_ij, c3_i1j, c3_ij1, fax_ijk, fay_ijk) &
    !$omp private (ipl, cmass, ctp, m_k) 
    do j=1,jm

      jmi=max(1,j-1)

      do i=1,im

        imi=i-1
        if (imi.eq.0) imi=im

        ! initialize vertically integrated fluxes
        faxdse(i,j) = 0._wp      
        faxdse_psi(i,j) = 0._wp
        faxdse_psi_topo(i,j) = 0._wp
        faxwtr(i,j) = 0._wp
        faxdst(i,j) = 0._wp
        faxco2(i,j) = 0._wp

        faydse(i,j) = 0._wp
        faydse_psi(i,j) = 0._wp
        faydse_psi_topo(i,j) = 0._wp
        faywtr(i,j) = 0._wp
        faydst(i,j) = 0._wp
        fayco2(i,j) = 0._wp

        fdxdse(i,j) = 0._wp
        fdxwtr(i,j) = 0._wp       
        fdxdst(i,j) = 0._wp       
        fdxco2(i,j) = 0._wp       

        fdydse(i,j) = 0._wp      
        fdywtr(i,j) = 0._wp
        fdydst(i,j) = 0._wp 
        fdyco2(i,j) = 0._wp 

        ipl = modulo(i,im) + 1
        cmass = 0._wp
        ctp   = 0._wp
        do k=1,km-2                                    ! troposphere, as for the diffusion
          m_k   = 0.5_wp*(dplx(i,j,k)+dplx(ipl,j,k))   ! cell layer mass
          cmass = cmass + m_k
          ctp   = ctp   + m_k*tp(i,j,k)
        enddo
        if (cmass.gt.0._wp) then
          tp_col(i,j) = ctp/cmass
        else
          tp_col(i,j) = tp(i,j,1)
        endif

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

          ! upstream zonal advection: upstream is the west cell (imi)
          ! for fax>0, the current cell (i) for fax<=0.
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
          faxdse(i,j) = faxdse(i,j) + fax_ijk*tpup ! kg/s * K
          faxdse_psi(i,j)      = faxdse_psi(i,j)      + fax_psi(i,j,k)*tpup
          faxdse_psi_topo(i,j) = faxdse_psi_topo(i,j) + fax_psi_topo(i,j,k)*tpup
          faxwtr(i,j) = faxwtr(i,j) + fax_ijk*qup  ! kg/s * kg/kg
          faxdst(i,j) = faxdst(i,j) + fax_ijk*dup
          faxco2(i,j) = faxco2(i,j) + fax_ijk*cup  ! kg/s * kgCO2/kg = kgCO2/s

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
          faydse(i,j) = faydse(i,j) + fay_ijk*tpup ! kg/s * K
          faydse_psi(i,j)      = faydse_psi(i,j)      + fay_psi(i,j,k)*tpup
          faydse_psi_topo(i,j) = faydse_psi_topo(i,j) + fay_psi_topo(i,j,k)*tpup
          ! the same mass flux and upstream DSE the advection uses, kept per level so that the
          ! zonal mean part can be subtracted after the i loop
          faywtr(i,j) = faywtr(i,j) + fay_ijk*qup  ! kg/s * kg/kg
          faydst(i,j) = faydst(i,j) + fay_ijk*dup
          fayco2(i,j) = fayco2(i,j) + fay_ijk*cup

          !-----------------------------------
          ! diffusive fluxes
          !-----------------------------------

          if (k.le.km-2) then   ! limit to troposphere

            !-----------------------------------
            ! zonal diffusive fluxes
            dpl_x = dplx(i,j,k)
            fdxdse(i,j) = fdxdse(i,j) + diffxdse(i,j)*dy*dpl_x*(tp_i1jk-tp_ijk)/dxt(j) ! m2/s * K * kg/m2 = kg/s * K
            fdxwtr(i,j) = fdxwtr(i,j) + diffxwtr(i,j)*dy*dpl_x*(q3_i1jk-q3_ijk)/dxt(j) 
            fdxdst(i,j) = fdxdst(i,j) + diffxdst(i,j)*dy*dpl_x*(d3_i1jk-d3_ijk)/dxt(j)
            fdxco2(i,j) = fdxco2(i,j) + diffxdst(i,j)*dy*dpl_x*(c3_i1j-c3_ij)/dxt(j)

            !-----------------------------------
            ! meridional diffusive fluxes
            dpl_y = dply(i,j,k)
            fdydse(i,j) = fdydse(i,j) + diffydse(i,j)*dxu(j)*dpl_y*(tp_ijk-tp_ij1k)/dy
            fdywtr(i,j) = fdywtr(i,j) + diffywtr(i,j)*dxu(j)*dpl_y*(q3_ijk-q3_ij1k)/dy
            fdydst(i,j) = fdydst(i,j) + diffydst(i,j)*dxu(j)*dpl_y*(d3_ijk-d3_ij1k)/dy
            fdyco2(i,j) = fdyco2(i,j) + diffydst(i,j)*dxu(j)*dpl_y*(c3_ij-c3_ij1)/dy

          endif

        enddo

      enddo

      ! no-flux condition at the poles
      if (j.eq.jm) then
        faydse(:,jmc) = 0._wp              
        faydse_psi(:,jmc) = 0._wp
        faydse_psi_topo(:,jmc) = 0._wp
        faywtr(:,jmc) = 0._wp          
        faydst(:,jmc) = 0._wp          
        fayco2(:,jmc) = 0._wp          
        fdydse(:,jmc) = 0._wp              
        fdywtr(:,jmc) = 0._wp        
        fdydst(:,jmc) = 0._wp 
        fdyco2(:,jmc) = 0._wp 
      endif

      ! Cycling
      faxdse(imc,j) = faxdse(1,j)              
      faxdse_psi(imc,j) = faxdse_psi(1,j)
      faxdse_psi_topo(imc,j) = faxdse_psi_topo(1,j)
      faxwtr(imc,j) = faxwtr(1,j)          
      faxdst(imc,j) = faxdst(1,j)          
      faxco2(imc,j) = faxco2(1,j)          
      fdxdse(imc,j) = fdxdse(1,j)              
      fdxwtr(imc,j) = fdxwtr(1,j)        
      fdxdst(imc,j) = fdxdst(1,j)        
      fdxco2(imc,j) = fdxco2(1,j)        

    enddo
    !$omp end parallel do


    !-----------------------------------
    ! fluxes convergency
    !-----------------------------------

    !$omp parallel do collapse(2) private(i,j,k,fdivdse,fdivwtr,fdivdst,fdivco2,mc_psi,mc_psi_topo)
    do j=1,jm
      do i=1,im

        ! Diffusive flux divergence, only added if l_diff_impl==false
        if (l_diff_impl) then
          fdivdse = 0._wp
          fdivwtr = 0._wp
          fdivdst = 0._wp
          fdivco2 = 0._wp
        else
          fdivdse = fdxdse(i,j)-fdxdse(i+1,j) + fdydse(i,j+1)-fdydse(i,j)
          fdivwtr = fdxwtr(i,j)-fdxwtr(i+1,j) + fdywtr(i,j+1)-fdywtr(i,j)
          fdivdst = fdxdst(i,j)-fdxdst(i+1,j) + fdydst(i,j+1)-fdydst(i,j)
          fdivco2 = fdxco2(i,j)-fdxco2(i+1,j) + fdyco2(i,j+1)-fdyco2(i,j)
        endif

        !-----------------------------------
        ! dry static energy (advection [+ explicit diffusion])
        convdse(i,j)= &
                       (faxdse(i,j)  -faxdse(i+1,j) &
                       +faydse(i,j+1)-faydse(i,j) &
                       +fdivdse) &
                       /sqr(i,j) * cp  ! K * kg/s / m2 * J/kg/K = J/m2/s = W/m2

        ! The part of it carried by each half of the column mass correction.
        !
        ! Neither half is mass conserving on its own - grad(psi) removes the dynamic
        ! share of the spurious column convergence and grad(psi_topo) the topographic
        ! share, and only their sum with the uncorrected flux closes the column. Each
        ! is therefore reported against a reference: the SAME column mass correction
        ! spread through the troposphere in proportion to layer mass, which is the
        ! least perturbing way of removing it. What is written out is the excess over
        ! that reference,
        !     sum_k mc_k*(tp_k - tp_trop) * cp / area ,
        ! so it answers how much the CHOICE OF LEVEL for the compensation adds to the
        ! local dry static energy budget. It is in W/m2 and directly comparable with
        ! convdse. A column that exports mass, mc<0, from levels above the tropospheric
        ! mean reads negative: the compensation is cooling that column, and it cools it
        ! more the higher the band sits.
        !
        ! Note that taken raw, without any reference, each part would instead carry a
        ! term (column mass imbalance)*tp of order fac*cp*theta/area - several hundred
        ! W/m2, an order of magnitude above the placement signal, cancelling only when
        ! the three are added back up.
        !
        ! The reference is the TROPOSPHERIC mean, not the whole column, because both
        ! bands lie in the troposphere and because the two stratospheric layers carry a
        ! fifth of the column pressure at theta of 380-520 K. That drags a full column
        ! mean to ~334 K at 67N, above the tropopause value itself, so every possible
        ! tropospheric placement would read as a large positive anomaly and a band
        ! sitting at the tropopause would read as harmless.
        !
        ! Subtracting a reference does not break the decomposition: the three column
        ! imbalances sum to zero, so the subtracted pieces do too and the parts still
        ! add up to the advective convdse. Only the split is reference dependent - the
        ! sum, and the difference between two runs, are not.
        mc_psi      = 0._wp
        mc_psi_topo = 0._wp
        do k=1,km
          mc_psi      = mc_psi      + fax_psi(i,j,k)     -fax_psi(i+1,j,k) &
                                    + fay_psi(i,j+1,k)   -fay_psi(i,j,k)
          mc_psi_topo = mc_psi_topo + fax_psi_topo(i,j,k)-fax_psi_topo(i+1,j,k) &
                                    + fay_psi_topo(i,j+1,k)-fay_psi_topo(i,j,k)
        enddo
        convdse_psi(i,j) = &
                       (faxdse_psi(i,j)  -faxdse_psi(i+1,j) &
                       +faydse_psi(i,j+1)-faydse_psi(i,j) &
                       -mc_psi*tp_col(i,j)) &
                       /sqr(i,j) * cp
        convdse_psi_topo(i,j) = &
                       (faxdse_psi_topo(i,j)  -faxdse_psi_topo(i+1,j) &
                       +faydse_psi_topo(i,j+1)-faydse_psi_topo(i,j) &
                       -mc_psi_topo*tp_col(i,j)) &
                       /sqr(i,j) * cp

        !-----------------------------------
        ! water (advection [+ explicit diffusion])
        convwtr_adv(i,j)= (faxwtr(i,j)  -faxwtr(i+1,j) &
                          +faywtr(i,j+1)-faywtr(i,j)) &
                          /sqr(i,j)       ! kg/kg * kg/s / m2 = kg/m2/s
        convwtr_dif(i,j)= fdivwtr/sqr(i,j)

        !-----------------------------------
        ! dust (advection [+ explicit diffusion])
        convdst(i,j)= &
                       (faxdst(i,j)  -faxdst(i+1,j) &
                       +faydst(i,j+1)-faydst(i,j) &
                       +fdivdst) &
                       /sqr(i,j)

        !-----------------------------------
        ! carbon (advection [+ explicit diffusion])
        convco2(i,j)= &
                       (faxco2(i,j)  -faxco2(i+1,j) &
                       +fayco2(i,j+1)-fayco2(i,j) &
                       +fdivco2) &
                       /sqr(i,j)        ! kgCO2/s/m2

      enddo
    enddo
    !$omp end parallel do

    return

  end subroutine adifa

end module adifa_mod

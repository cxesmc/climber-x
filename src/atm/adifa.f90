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
  use atm_params, only : l_diff_impl, l_dust, i_adv, c_adv_polar, lat_adv_polar
  use atm_params, only : l_write_timer
  use atm_grid, only : im, imc, jm, jmc, km
  use atm_grid, only : dplt, dplx, dply, dplx_dif, dply_dif, dy, dxt, dxu, sqr, fit
  use constants, only : pi
  use adv_fct_mod, only : fct_level, fct_diag_add
  !$ use omp_lib

  implicit none

  private
  public :: adifa

  ! flux corrected transport (i_adv=2) work space, allocated on first use.
  ! rx/ry are the antidiffusive limiter ratios on the u- and v-faces; the CO2
  ! ones carry a single level because cam is uniform over the column.
  real(wp), allocatable, dimension(:,:,:) :: rx_dse, rx_wtr, rx_dst, rx_co2
  real(wp), allocatable, dimension(:,:,:) :: ry_dse, ry_wtr, ry_dst, ry_co2
  real(wp), allocatable, dimension(:,:,:) :: fmx_col, fmy_col, dm_col, cam_col

contains
    
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  a d i f a
  !   Purpose    :  compute advective and diffusive fluxes of energy, 
  !              :  water and dust
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine adifa(fax, fay, fax_psi, fay_psi, fax_psi_topo, fay_psi_topo, &
    tp, q3, d3, cam, diffx, diffy, &
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
    real(wp), intent(in   ) :: diffx(:,:)
    real(wp), intent(in   ) :: diffy(:,:)
    
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
    integer :: idg_dse, idg_wtr
    integer :: n, nfct
    real(wp), dimension(4,jm,km) :: dacc_dse, dacc_wtr   ! per level limiter diagnostics
    real(wp), dimension(4,jm,1)  :: dacc_nul             ! not used, the diagnostics are off for these
    real(wp), dimension(imc,jm) :: faxdse_psi, faxdse_psi_topo
    ! tropospheric mass weighted mean of tp, the reference the parts are measured against
    real(wp), dimension(im,jm) :: tp_col
    real(wp), dimension(im,jmc) :: faydse_psi, faydse_psi_topo


    !-----------------------------------
    ! flux corrected transport: antidiffusive limiter ratios
    !-----------------------------------
    ! Computed here, before the flux loop, because the limiter needs the whole
    ! horizontal field of a level at once. They are then applied face by face
    ! below, where they shift the face value from upstream towards centred.
    if (i_adv.eq.2) then

      if (.not.allocated(rx_dse)) then
        allocate(rx_dse(imc,jm,km), ry_dse(im,jmc,km))
        allocate(rx_wtr(imc,jm,km), ry_wtr(im,jmc,km))
        allocate(rx_dst(imc,jm,km), ry_dst(im,jmc,km))
        allocate(rx_co2(imc,jm,1),  ry_co2(im,jmc,1))
        allocate(fmx_col(imc,jm,1), fmy_col(im,jmc,1))
        allocate(dm_col(im,jm,1),   cam_col(im,jm,1))
        rx_dst = 0._wp
        ry_dst = 0._wp
      endif

      ! the limiter activity diagnostics ride along with the rest of the timer
      ! output; without it fct_ratios does no book keeping at all
      idg_dse = 0
      idg_wtr = 0
      if (l_write_timer) then
        idg_dse = 1
        idg_wtr = 2
      endif

      ! CO2 is uniform over the column, so its limiter is a single level problem
      ! carried by the column integrated mass flux and the column mass
      !$omp parallel do private(i,j,k)
      do j=1,jmc
        do i=1,imc
          if (j.le.jm) then
            fmx_col(i,j,1) = 0._wp
            do k=1,km
              fmx_col(i,j,1) = fmx_col(i,j,1) + fax(i,j,k)
            enddo
          endif
          if (i.le.im) then
            fmy_col(i,j,1) = 0._wp
            do k=1,km
              fmy_col(i,j,1) = fmy_col(i,j,1) + fay(i,j,k)
            enddo
          endif
          if (i.le.im .and. j.le.jm) then
            dm_col(i,j,1) = 0._wp
            do k=1,km
              dm_col(i,j,1) = dm_col(i,j,1) + dplt(i,j,k)
            enddo
            cam_col(i,j,1) = cam(i,j)
          endif
        enddo
      enddo
      !$omp end parallel do

      ! The limiter of each level is independent of the others, and so is that of each tracer:
      ! all (tracer, level) pairs are done in one parallel loop, as there are only km levels
      ! per tracer. Items 1:km are dse, km+1:2km water, 2km+1 CO2 and 2km+2:3km+1 dust.
      nfct = 2*km+1
      if (l_dust) nfct = nfct+km
      !$omp parallel do schedule(dynamic,1) private(n)
      do n=1,nfct
        if (n.le.km) then
          call fct_level(n, fax, fay, dplt, tp, rx_dse, ry_dse, dacc_dse(:,:,n), idg_dse)
        else if (n.le.2*km) then
          call fct_level(n-km, fax, fay, dplt, q3, rx_wtr, ry_wtr, dacc_wtr(:,:,n-km), idg_wtr)
        else if (n.eq.2*km+1) then
          call fct_level(1, fmx_col, fmy_col, dm_col, cam_col, rx_co2, ry_co2, dacc_nul(:,:,1), 0)
        else
          call fct_level(n-2*km-1, fax, fay, dplt, d3, rx_dst, ry_dst, dacc_nul(:,:,1), 0)
        endif
      enddo
      !$omp end parallel do
      if (idg_dse.gt.0) call fct_diag_add(km, dacc_dse, idg_dse)
      if (idg_wtr.gt.0) call fct_diag_add(km, dacc_wtr, idg_wtr)

    endif


    !$omp parallel do collapse(2) private(i, j, k, imi, jmi, tpup, qup, dup, cup, dpl_x, dpl_y) &
    !$omp private (tp_ijk, tp_i1jk, tp_ij1k, q3_ijk, q3_i1jk, q3_ij1k, d3_ijk, d3_i1jk, d3_ij1k, c3_ij, c3_i1j, c3_ij1, fax_ijk, fay_ijk) &
    !$omp private (ipl, cmass, ctp, m_k) 
    do j=1,jm
      do i=1,im

        jmi=max(1,j-1)
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
          if (i_adv.eq.2) then
            ! move the face value from upstream towards centred by the limiter
            ! ratio; the same face value multiplies fax, fax_psi and fax_psi_topo
            ! below, so the three still add up to the total advective flux
            tpup = tpup + rx_dse(i,j,k)*(0.5_wp*(tp_ijk+tp_i1jk)-tpup)
            qup  = qup  + rx_wtr(i,j,k)*(0.5_wp*(q3_ijk+q3_i1jk)-qup)
            dup  = dup  + rx_dst(i,j,k)*(0.5_wp*(d3_ijk+d3_i1jk)-dup)
            cup  = cup  + rx_co2(i,j,1)*(0.5_wp*(c3_ij+c3_i1j)-cup)
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
          if (i_adv.eq.2) then
            tpup = tpup + ry_dse(i,j,k)*(0.5_wp*(tp_ijk+tp_ij1k)-tpup)
            qup  = qup  + ry_wtr(i,j,k)*(0.5_wp*(q3_ijk+q3_ij1k)-qup)
            dup  = dup  + ry_dst(i,j,k)*(0.5_wp*(d3_ijk+d3_ij1k)-dup)
            cup  = cup  + ry_co2(i,j,1)*(0.5_wp*(c3_ij+c3_ij1)-cup)
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
            ! zonal diffusive fluxes; the face mass carries the vertical profile of the diffusivity (h_diff)
            dpl_x = dplx_dif(i,j,k)
            fdxdse(i,j) = fdxdse(i,j) + diffx(i,j)*dy*dpl_x*(tp_i1jk-tp_ijk)/dxt(j) ! m2/s * K * kg/m2 = kg/s * K
            fdxwtr(i,j) = fdxwtr(i,j) + diffx(i,j)*dy*dpl_x*(q3_i1jk-q3_ijk)/dxt(j) 
            fdxdst(i,j) = fdxdst(i,j) + diffx(i,j)*dy*dpl_x*(d3_i1jk-d3_ijk)/dxt(j)
            fdxco2(i,j) = fdxco2(i,j) + diffx(i,j)*dy*dpl_x*(c3_i1j-c3_ij)/dxt(j)

            !-----------------------------------
            ! meridional diffusive fluxes
            dpl_y = dply_dif(i,j,k)
            fdydse(i,j) = fdydse(i,j) + diffy(i,j)*dxu(j)*dpl_y*(tp_ijk-tp_ij1k)/dy
            fdywtr(i,j) = fdywtr(i,j) + diffy(i,j)*dxu(j)*dpl_y*(q3_ijk-q3_ij1k)/dy
            fdydst(i,j) = fdydst(i,j) + diffy(i,j)*dxu(j)*dpl_y*(d3_ijk-d3_ij1k)/dy
            fdyco2(i,j) = fdyco2(i,j) + diffy(i,j)*dxu(j)*dpl_y*(c3_ij-c3_ij1)/dy

          endif

        enddo

      enddo
    enddo
    !$omp end parallel do

    ! zonal sub-cycling in rows that exceed the zonal Courant limit: overwrites the zonal advective
    ! fluxes of those rows only, after (and outside) the main loop
    call polar_zonal_subcycle(fax, fax_psi, fax_psi_topo, tp, q3, d3, cam, &
      faxdse, faxdse_psi, faxdse_psi_topo, faxwtr, faxdst, faxco2)

    ! no-flux condition at the poles
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

    ! Cycling
    faxdse(imc,:) = faxdse(1,:)              
    faxdse_psi(imc,:) = faxdse_psi(1,:)
    faxdse_psi_topo(imc,:) = faxdse_psi_topo(1,:)
    faxwtr(imc,:) = faxwtr(1,:)          
    faxdst(imc,:) = faxdst(1,:)          
    faxco2(imc,:) = faxco2(1,:)          
    fdxdse(imc,:) = fdxdse(1,:)              
    fdxwtr(imc,:) = fdxwtr(1,:)        
    fdxdst(imc,:) = fdxdst(1,:)        
    fdxco2(imc,:) = fdxco2(1,:)        


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


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p o l a r _ z o n a l _ s u b c y c l e
  !   Purpose    :  zonal tracer advection in rows that violate the zonal Courant condition
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Applied only poleward of lat_adv_polar. Near the poles dxt is a few tens of km, and a few m/s of azonal geostrophic wind give a zonal
  ! Courant number above 1 on the fast step. There the FCT limiter falls back to upstream, which
  ! is itself unstable for C>1 (output/cflchk: C_x 2.3 at 87.5N, moisture goes negative, crash).
  ! In each row whose zonal inflow-mass Courant number C exceeds c_adv_polar, every level and
  ! tracer is advected along the row with the zonal mass fluxes only, in nsub = ceil(C/c_adv_polar)
  ! upstream sub-steps of the advective form used by the FCT low order leg,
  !   x <- x + dt/nsub/m * sum_inflow F*(x_donor - x) ,
  ! each a convex combination. The time mean of the upstream face value over the sub-steps
  ! replaces the zonal face value of the row, and the zonal fluxes are rebuilt from it. The mass
  ! fluxes are untouched, so fax, fax_psi and fax_psi_topo still add up and the column budgets
  ! stay consistent with the meridional fluxes; tracers are conserved because the fluxes are on
  ! shared faces. CO2 uses the column mass flux, as its FCT limiter does.
  subroutine polar_zonal_subcycle(fax, fax_psi, fax_psi_topo, tp, q3, d3, cam, &
      faxdse, faxdse_psi, faxdse_psi_topo, faxwtr, faxdst, faxco2)

    implicit none

    real(wp), intent(in   ) :: fax(:,:,:), fax_psi(:,:,:), fax_psi_topo(:,:,:)
    real(wp), intent(in   ) :: tp(:,:,:), q3(:,:,:), d3(:,:,:), cam(:,:)
    real(wp), intent(inout) :: faxdse(:,:), faxdse_psi(:,:), faxdse_psi_topo(:,:)
    real(wp), intent(inout) :: faxwtr(:,:), faxdst(:,:), faxco2(:,:)

    integer :: i, j, k, n, r, nsub, nrow
    integer :: jrow(jm), nsub_row(jm)
    real(wp) :: c, cmax, m
    real(wp) :: fcol(imc), mcol(im)
    real(wp), allocatable :: xb_dse(:,:,:), xb_wtr(:,:,:), xb_dst(:,:,:), xb_co2(:,:)
    integer, save :: nsub_warned(jm) = 1


    ! select the rows to treat and their number of sub-steps
    nrow = 0
    do j=1,jm

      ! polar rows only
      if (abs(fit(j))*180._wp/pi.le.lat_adv_polar) cycle

      ! zonal inflow-mass Courant number of the row, over all levels
      cmax = 0._wp
      do i=1,im
        do k=1,km
          m = dplt(i,j,k)*sqr(i,j)
          if (m.gt.0._wp) then
            c = tstep*(max(0._wp,fax(i,j,k)) + max(0._wp,-fax(i+1,j,k)))/m
            cmax = max(cmax,c)
          endif
        enddo
      enddo
      if (cmax.le.c_adv_polar) cycle

      nsub = ceiling(cmax/c_adv_polar)
      if (l_write_timer .and. nsub.gt.nsub_warned(j)) then
        print '(a,i3,a,f7.2,a,i3)', ' atm polar zonal sub-cycling: j =',j,' , zonal Courant =',cmax,' , nsub =',nsub
        nsub_warned(j) = nsub
      endif

      nrow = nrow + 1
      jrow(nrow) = j
      nsub_row(nrow) = nsub

    enddo
    if (nrow.eq.0) return

    allocate(xb_dse(im,km,nrow), xb_wtr(im,km,nrow), xb_dst(im,km,nrow), xb_co2(im,nrow))

    ! each (row, tracer, level) is independent. Items 1:km are dse, km+1:2km water, 2km+1:3km dust
    ! and 3km+1 CO2.
    !$omp parallel do collapse(2) schedule(dynamic,1) private(r, n, i, j, k, m, fcol, mcol)
    do r=1,nrow
      do n=1,3*km+1
        j = jrow(r)
        if (n.le.km) then
          k = n
          call subcycle_row(nsub_row(r), fax(:,j,k), dplt(:,j,k)*sqr(:,j), tp(:,j,k), xb_dse(:,k,r))
        else if (n.le.2*km) then
          k = n-km
          call subcycle_row(nsub_row(r), fax(:,j,k), dplt(:,j,k)*sqr(:,j), q3(:,j,k), xb_wtr(:,k,r))
        else if (n.le.3*km) then
          k = n-2*km
          call subcycle_row(nsub_row(r), fax(:,j,k), dplt(:,j,k)*sqr(:,j), d3(:,j,k), xb_dst(:,k,r))
        else
          ! CO2 uses the column mass flux and mass
          do i=1,im
            mcol(i) = 0._wp
            do k=1,km
              m = dplt(i,j,k)*sqr(i,j)
              mcol(i) = mcol(i) + m
            enddo
          enddo
          do i=1,imc
            fcol(i) = 0._wp
            do k=1,km
              fcol(i) = fcol(i) + fax(i,j,k)
            enddo
          enddo
          call subcycle_row(nsub_row(r), fcol, mcol, cam(:,j), xb_co2(:,r))
        endif
      enddo
    enddo
    !$omp end parallel do

    ! rebuild the zonal advective fluxes of the rows from the time mean face values
    !$omp parallel do private(r, i, j, k)
    do r=1,nrow
      j = jrow(r)
      do i=1,im
        faxdse(i,j) = 0._wp
        faxdse_psi(i,j) = 0._wp
        faxdse_psi_topo(i,j) = 0._wp
        faxwtr(i,j) = 0._wp
        faxdst(i,j) = 0._wp
        faxco2(i,j) = 0._wp
        do k=1,km
          faxdse(i,j) = faxdse(i,j) + fax(i,j,k)*xb_dse(i,k,r)
          faxdse_psi(i,j) = faxdse_psi(i,j) + fax_psi(i,j,k)*xb_dse(i,k,r)
          faxdse_psi_topo(i,j) = faxdse_psi_topo(i,j) + fax_psi_topo(i,j,k)*xb_dse(i,k,r)
          faxwtr(i,j) = faxwtr(i,j) + fax(i,j,k)*xb_wtr(i,k,r)
          faxdst(i,j) = faxdst(i,j) + fax(i,j,k)*xb_dst(i,k,r)
          faxco2(i,j) = faxco2(i,j) + fax(i,j,k)*xb_co2(i,r)
        enddo
      enddo
    enddo
    !$omp end parallel do

    deallocate(xb_dse, xb_wtr, xb_dst, xb_co2)

    return

  end subroutine polar_zonal_subcycle


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s u b c y c l e _ r o w
  !   Purpose    :  nsub upstream sub-steps of one tracer along one periodic row,
  !                 returns the time mean upstream value on the zonal faces
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine subcycle_row(nsub, f, m, x0, xbar)

    implicit none

    integer,  intent(in ) :: nsub
    real(wp), intent(in ) :: f(:)      ! (imc) zonal mass flux, face i between cells i-1 and i, positive eastward, kg/s
    real(wp), intent(in ) :: m(:)      ! (im) cell mass, kg
    real(wp), intent(in ) :: x0(:)     ! (im) tracer at the start of the step
    real(wp), intent(out) :: xbar(:)   ! (im) time mean upstream value on face i

    integer :: i, s, iw, ie
    real(wp) :: x(im), xn(im), dts


    dts = tstep/real(nsub,wp)
    x(:) = x0(1:im)
    xbar(:) = 0._wp
    do s=1,nsub
      do i=1,im
        iw = modulo(i-2,im) + 1
        if (f(i).gt.0._wp) then
          xbar(i) = xbar(i) + x(iw)
        else
          xbar(i) = xbar(i) + x(i)
        endif
      enddo
      do i=1,im
        iw = modulo(i-2,im) + 1
        ie = modulo(i,im) + 1
        if (m(i).gt.0._wp) then
          xn(i) = x(i) + dts/m(i)*( max(0._wp,f(i))*(x(iw)-x(i)) + max(0._wp,-f(i+1))*(x(ie)-x(i)) )
        else
          xn(i) = x(i)
        endif
      enddo
      x(:) = xn(:)
    enddo
    xbar(:) = xbar(:)/real(nsub,wp)

    return

  end subroutine subcycle_row

end module adifa_mod

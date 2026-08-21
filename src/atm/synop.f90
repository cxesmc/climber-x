!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s y n o p _ m o d
!
!  Purpose : synoptic processes
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
module synop_mod

  use atm_params, only : wp
  use constants, only : g, omega
  use atm_grid, only : im, jm, nm, imc, jmc, km, dxt, dy, zl, k850, k700, k500, i_ocn, sint, cost 
  use atm_params, only : tstep, ra
  use atm_params, only : c_syn_1, c_syn_2, c_syn_3, c_syn_4, c_syn_5, c_syn_6, c_syn_7, c_syn_8, windmin, synsurmin, c_wind_ele
  use atm_params, only : tau_fac
  use atm_params, only : c_diff_dse, i_diff_wtr, c_diff_wtr, l_diff_impl, c_diffx_pol
  use smooth_atm_mod, only : zofil
  use tridiag, only : tridiag_solve, cyclic_tridiag_solve
  use timer, only : dt_atm
  !$ use omp_lib

  implicit none

  private
  public :: synop

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s y n o p 
  !   Purpose    :  1) compute macro-turbulent horizontal diffusion coefficients
  !                 2) compute synoptic wind and module of surface wind
  !                 3) compute synoptic vertical velocity on cloudiness level
  !                 4) compute zonal surface wind stress over the ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine synop(frst, zs, ut3f, vt3f, u700, v700, us, vs, tp, zsa, cda, cd, epsa, sam, cdif, &
      synprod, syndiss, synadv, syndif, synsur, winda, wind, taux, tauy, diffxdse, diffydse, diffxwtr, diffywtr, diffxdst, diffydst, wsyn)

    implicit none

    real(wp), intent(in   ) :: frst(:,:,:)
    real(wp), intent(in   ) :: zs(:,:,:)
    real(wp), intent(in   ) :: ut3f(:,:,:)
    real(wp), intent(in   ) :: vt3f(:,:,:)
    real(wp), intent(in   ) :: u700(:,:)
    real(wp), intent(in   ) :: v700(:,:)
    real(wp), intent(in   ) :: us(:,:,:)
    real(wp), intent(in   ) :: vs(:,:,:)
    real(wp), intent(in   ) :: tp(:,:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: cda(:,:)
    real(wp), intent(in   ) :: cd(:,:,:)
    real(wp), intent(in   ) :: epsa(:,:,:)

    real(wp), intent(inout) :: sam(:,:)
    real(wp), intent(inout) :: cdif(:,:)

    real(wp), intent(out  ) :: synprod(:,:)
    real(wp), intent(out  ) :: syndiss(:,:)
    real(wp), intent(out  ) :: synadv(:,:)
    real(wp), intent(out  ) :: syndif(:,:)
    real(wp), intent(out  ) :: synsur(:,:,:)
    real(wp), intent(out  ) :: winda(:,:)
    real(wp), intent(out  ) :: wind(:,:,:)
    real(wp), intent(out  ) :: taux(:,:,:)
    real(wp), intent(out  ) :: tauy(:,:,:)
    real(wp), intent(out  ) :: diffxdse(:,:)
    real(wp), intent(out  ) :: diffydse(:,:)
    real(wp), intent(out  ) :: diffxwtr(:,:)
    real(wp), intent(out  ) :: diffywtr(:,:)
    real(wp), intent(out  ) :: diffxdst(:,:)
    real(wp), intent(out  ) :: diffydst(:,:)
    real(wp), intent(out  ) :: wsyn(:,:)

    integer :: i, j, k, n, imi, ipl, jmi, jpl
    real(wp) :: dsdt, sqsam
    real(wp) :: sxadv, syadv, sxdif, sydif
    real(wp) :: uef, vef
    real(wp) :: dudz, dvdz
    real(wp) :: ugrad, Nfreq
    real(wp) :: diffxmx(jm)
    real(wp) :: sam_sqrt(im,jm)
    ! implicit (ADI) solve of the synoptic-energy (sam) equation
    real(wp) :: kdiss(im,jm), samx(im,jm)
    real(wp) :: gx, gy
    real(wp) :: a(im), b(im), c(im), r(im), sol(im)
    real(wp) :: aj(jm), bj(jm), cj(jm), rj(jm), solj(jm)


    do j=1,jm
      ! maximum zonal diffusion coefficient for CFL stability 
      diffxmx(j) = 0.5_wp*dxt(j)**2/tstep 
    enddo 

    !-----------------------------------------------
    ! synoptic kinetic energy
    !-----------------------------------------------

    !$omp parallel do collapse(2) private(i, j, k, uef, vef, imi, ipl, jmi, jpl, ugrad, dudz, dvdz, Nfreq, sxadv, syadv, sxdif, sydif)
    do i=1,im
      do j=2,jm-1

        imi=i-1
        if (imi.lt.1) imi=im
        ipl=i+1
        if (ipl.gt.im) ipl=1
        jmi=j-1
        if (jmi.lt.1) jmi=1
        jpl=j+1
        if (jpl.gt.jm) jpl=jm        

        !-----------------------------------------------
        ! synoptic energy production

        ! proportional to max Eady model baroclinic growth rate (e.g. Hoskins and Valdes 1990)

        Nfreq = sqrt(g/tp(i,j,k700)*(tp(i,j,k500)-tp(i,j,k850))/(zl(k500)-zl(k850)))   ! Brunt-Vaisala frequency
        dudz = (ut3f(i,j,k500)-ut3f(i,j,k850))/(zl(k500)-zl(k850))
        dvdz = (vt3f(i,j,k500)-vt3f(i,j,k850))/(zl(k500)-zl(k850)) * cost(j)
        ugrad = sqrt(dudz**2+dvdz**2)
        synprod(i,j) = c_syn_1 + c_syn_2 * 2._wp*omega*abs(sint(j)) / Nfreq * ugrad * (1._wp-c_syn_8*zsa(i,j)/3000._wp)
        synprod(i,j) = max(0._wp,synprod(i,j))

        !-----------------------------------------------
        ! synoptic energy dissipation  

        syndiss(i,j) = (c_syn_3+c_syn_4*cda(i,j))*sam(i,j)*sqrt(sam(i,j))

        !-----------------------------------------------
        ! synoptic energy advection         

        ! effective horizontal velocities
        ! use steering wind at ~700 hPa for advection
        uef = u700(i,j)   
        vef = v700(i,j)
        if (uef.gt.0._wp) then
          sxadv = uef*(sam(imi,j)-sam(i,j))/dxt(j)
        else 
          sxadv = uef*(sam(i,j)-sam(ipl,j))/dxt(j)
        endif
        if (vef.gt.0._wp) then
          syadv = vef*(sam(i,jpl)-sam(i,j))/dy
        else         
          syadv = vef*(sam(i,j)-sam(i,jmi))/dy 
        endif 
        synadv(i,j) = sxadv+syadv

        !-----------------------------------------------
        ! synoptic energy diffusion

        sxdif = min(cdif(i,j),diffxmx(j))*(sam(imi,j)+sam(ipl,j)-2._wp*sam(i,j))/dxt(j)**2
        sydif = cdif(i,j)*(sam(i,jmi)+sam(i,jpl)-2._wp*sam(i,j))/dy**2    
        syndif(i,j) = sxdif+sydif

      enddo
    enddo
    !$omp end parallel do

    !-----------------------------------------------
    ! synoptic energy evolution - FULLY IMPLICIT (advection + diffusion + semi-implicit dissipation).
    ! Solve (I + dt_atm*(kdiss + Adv_upwind - cdif*Lap)) sam_new = sam + dt_atm*synprod
    ! by ADI: cyclic zonal sweep, then meridional sweep (Dirichlet poles). 

    do j=2,jm-1
      do i=1,im
        kdiss(i,j) = (c_syn_3+c_syn_4*cda(i,j))*sqrt(sam(i,j))
      enddo
    enddo

    ! zonal implicit sweep (advection + diffusion; periodic -> cyclic tridiagonal)
    !$omp parallel do private(i,j,uef,gx,a,b,c,r,sol)
    do j=2,jm-1
      do i=1,im
        uef = u700(i,j)
        gx  = dt_atm*cdif(i,j)/dxt(j)**2
        ! diffusion + upwind advection (sxadv: imi donor for uef>0, ipl for uef<0)
        a(i) = -gx - dt_atm*max(uef,0._wp)/dxt(j)
        c(i) = -gx + dt_atm*min(uef,0._wp)/dxt(j)
        ! full linearized dissipation applied once here
        b(i) = 1._wp + dt_atm*kdiss(i,j) - a(i) - c(i)
        r(i) = sam(i,j) + dt_atm*synprod(i,j)
      enddo
      call cyclic_tridiag_solve(a, b, c, r, sol, im)
      do i=1,im
        samx(i,j) = sol(i)
      enddo
    enddo
    !$omp end parallel do

    ! meridional implicit sweep (advection + diffusion; Dirichlet poles), per longitude
    !$omp parallel do private(i,j,vef,gy,aj,bj,cj,rj,solj)
    do i=1,im
      do j=2,jm-1
        vef = v700(i,j)
        gy  = dt_atm*cdif(i,j)/dy**2
        ! diffusion + upwind advection (syadv: jpl donor for vef>0, jmi for vef<0)
        aj(j) = -gy + dt_atm*min(vef,0._wp)/dy
        cj(j) = -gy - dt_atm*max(vef,0._wp)/dy
        bj(j) = 1._wp - aj(j) - cj(j)
        rj(j) = samx(i,j)
      enddo
      ! Dirichlet poles: move the j=1 / j=jm couplings (diffusion+advection) to the RHS
      rj(2)    = rj(2)    - aj(2)*sam(i,1)
      rj(jm-1) = rj(jm-1) - cj(jm-1)*sam(i,jm)
      aj(2)    = 0._wp
      cj(jm-1) = 0._wp
      call tridiag_solve(aj(2:jm-1), bj(2:jm-1), cj(2:jm-1), rj(2:jm-1), solj(2:jm-1), jm-2)
      do j=2,jm-1
        sam(i,j) = max(solj(j),1._wp)
      enddo
    enddo
    !$omp end parallel do

    ! polar fill and zonal smoothing

    do i=1,im
      sam(i,1)  = sam(i,2)
      sam(i,jm) = sam(i,jm-1)
    enddo 
    call zofil(sam,1,12)
    call zofil(sam,2,8)
    call zofil(sam,3,4)
    call zofil(sam,4,2)
    call zofil(sam,5,1)
    call zofil(sam,jm-4,1)
    call zofil(sam,jm-3,2)
    call zofil(sam,jm-2,4)
    call zofil(sam,jm-1,8)
    call zofil(sam,jm,12)

    
    !-----------------------------------------------
    ! Products: macro diffusion, wind, wind stress
    !-----------------------------------------------

    !$omp parallel do collapse(2) private(i,j,n,sqsam)
    do i=1,im
      do j=1,jm

        sqsam = sqrt(sam(i,j)) 
        sam_sqrt(i,j) = sqsam

        ! effective macrodiffusive coefficient          

        cdif(i,j) = c_syn_5*sqsam

        do n=1,nm

          ! synoptic surface wind
          synsur(i,j,n) = c_syn_6*sqsam*epsa(i,j,n)
          synsur(i,j,n) = max(synsur(i,j,n),synsurmin)

          ! total surface wind
          wind(i,j,n) = sqrt(us(i,j,n)**2+vs(i,j,n)**2+synsur(i,j,n)**2) + c_wind_ele*zs(i,j,n)
          if (n.eq.i_ocn) then
            wind(i,j,n) = max(wind(i,j,n),windmin)
          endif

          ! wind stress 
          taux(i,j,n) = tau_fac*cd(i,j,n)*ra*us(i,j,n)*wind(i,j,n)
          tauy(i,j,n) = tau_fac*cd(i,j,n)*ra*vs(i,j,n)*wind(i,j,n)

        enddo

        winda(i,j) = sum(wind(i,j,:)*frst(i,j,:))

        ! synoptic vertical velocity at ~700 hPa
        wsyn(i,j) = c_syn_7*sqsam

      enddo
    enddo
    !$omp end parallel do

    !-----------------------------------------------
    ! macro diffusivities for energy, water and dust

    do i=1,im
      imi=i-1
      if (imi.lt.1) imi=im
      do j=1,jm
        ! diffusivity for dry static energy, proportional to sqrt(EKE)
        diffxdse(i,j) = c_diff_dse * 0.5_wp*(sam_sqrt(imi,j)+sam_sqrt(i,j))
        ! diffusivity for water vapor, proportional to EKE (e.g. Caballero & Hanley, 2012)
        if (i_diff_wtr.eq.1) then
          diffxwtr(i,j) = c_diff_wtr * 0.5_wp*(sam(imi,j)+sam(i,j)) 
        else if (i_diff_wtr.eq.2) then
          diffxwtr(i,j) = c_diff_wtr * 0.5_wp*(sam_sqrt(imi,j)+sam_sqrt(i,j))
        endif
        ! diffusivity for dust
        diffxdst(i,j) = c_diff_dse * 0.5_wp*(sam_sqrt(imi,j)+sam_sqrt(i,j))
        ! limit the zonal diffusivities. Explicit scheme: CFL stability requires the
        ! zonal diffusion number diffx*tstep/dxt^2 < 0.5, i.e. diffx <= diffxmx.
        ! Implicit scheme: unconditionally stable, but near the pole (dxt->0) the
        ! zonal conductance ~1/dxt becomes very large => limit
        if (.not.l_diff_impl) then
          diffxdse(i,j) = min(diffxdse(i,j),diffxmx(j))
          diffxwtr(i,j) = min(diffxwtr(i,j),diffxmx(j))
          diffxdst(i,j) = min(diffxdst(i,j),diffxmx(j))
        else
          diffxdse(i,j) = min(diffxdse(i,j),2._wp*c_diffx_pol*diffxmx(j))
          diffxwtr(i,j) = min(diffxwtr(i,j),2._wp*c_diffx_pol*diffxmx(j))
          diffxdst(i,j) = min(diffxdst(i,j),2._wp*c_diffx_pol*diffxmx(j))
        endif
      enddo
    enddo 
    do j=1,jm
      diffxdse(imc,j) = diffxdse(1,j)
      diffxwtr(imc,j) = diffxwtr(1,j)
      diffxdst(imc,j) = diffxdst(1,j)
    enddo

    do j=2,jm
      jmi=j-1
      if (jmi.lt.1) jmi=1
      do i=1,im
        ! diffusivity for dry static energy, proportional to sqrt(EKE)
        diffydse(i,j) = c_diff_dse * 0.5_wp*(sam_sqrt(i,jmi)+sam_sqrt(i,j))
        ! diffusivity for water vapor, proportional to EKE (e.g. Caballero & Hanley, 2012)
        if (i_diff_wtr.eq.1) then
          diffywtr(i,j) = c_diff_wtr * 0.5_wp*(sam(i,jmi)+sam(i,j)) 
        else if (i_diff_wtr.eq.2) then
          diffywtr(i,j) = c_diff_wtr * 0.5_wp*(sam_sqrt(i,jmi)+sam_sqrt(i,j))
        endif
        ! diffusivity for dust 
        diffydst(i,j) = c_diff_dse * 0.5_wp*(sam_sqrt(i,jmi)+sam_sqrt(i,j))
      enddo
    enddo
    diffydse(:,1)   = 0._wp
    diffydse(:,jmc) = 0._wp
    diffywtr(:,1)   = 0._wp
    diffywtr(:,jmc) = 0._wp
    diffydst(:,1)   = 0._wp
    diffydst(:,jmc) = 0._wp


    return

  end subroutine synop

end module synop_mod


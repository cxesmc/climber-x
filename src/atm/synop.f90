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
  use constants, only : g, omega, pi, z_sfl
  use atm_grid, only : im, jm, nm, imc, jmc, km, dxt, dy, z3, k850, k700, k500, i_ocn, sint, cost, fit, cosu, dxtg 
  use atm_params, only : tstep, ra
  use atm_params, only : c_syn_1, c_syn_2, c_syn_3, c_syn_4, c_syn_5, c_syn_6, c_syn_7, windmin, synsurmin
  use atm_params, only : i_syn_diss, syn_e0, c_eke_slope, hgams, dz_syn_min
  use atm_params, only : tau_fac, c_tau_syn, l_tau_stab, f_tau_Ri_unstab, f_tau_Ri_stab
  use atm_params, only : c_diff
  use smooth_atm_mod, only : polar_filter
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
  subroutine synop(frst, zs, ut3f, vt3f, uadv, vadv, us, vs, tp, zsa, cda, cd, epsa, tskin, t2, sam, cdif, &
      synprod, synsur, winda, wind, taux, tauy, diffx, diffy, wsyn)

    implicit none

    real(wp), intent(in   ) :: frst(:,:,:)
    real(wp), intent(in   ) :: zs(:,:,:)
    real(wp), intent(in   ) :: ut3f(:,:,:)
    real(wp), intent(in   ) :: vt3f(:,:,:)
    real(wp), intent(in   ) :: uadv(:,:)
    real(wp), intent(in   ) :: vadv(:,:)
    real(wp), intent(in   ) :: us(:,:,:)
    real(wp), intent(in   ) :: vs(:,:,:)
    real(wp), intent(in   ) :: tp(:,:,:)
    real(wp), intent(in   ) :: zsa(:,:)
    real(wp), intent(in   ) :: cda(:,:)
    real(wp), intent(in   ) :: cd(:,:,:)
    real(wp), intent(in   ) :: epsa(:,:,:)
    real(wp), intent(in   ) :: tskin(:,:,:)
    real(wp), intent(in   ) :: t2(:,:,:)

    real(wp), intent(inout) :: sam(:,:)
    real(wp), intent(inout) :: cdif(:,:)

    real(wp), intent(out  ) :: synprod(:,:)
    real(wp), intent(out  ) :: synsur(:,:,:)
    real(wp), intent(out  ) :: winda(:,:)
    real(wp), intent(out  ) :: wind(:,:,:)
    real(wp), intent(out  ) :: taux(:,:,:)
    real(wp), intent(out  ) :: tauy(:,:,:)
    real(wp), intent(out  ) :: diffx(:,:)
    real(wp), intent(out  ) :: diffy(:,:)
    real(wp), intent(out  ) :: wsyn(:,:)

    integer :: i, j, k, n, imi, jmi, ipl
    integer :: iop, i_n, j_n, i_s, j_s
    real(wp) :: kn, ks
    real(wp) :: sqsam, wtau, richardson, f_stab

    real(wp), parameter :: sin45 = 0.7071067811865476_wp
    real(wp) :: uef, vef
    real(wp) :: dudz, dvdz
    real(wp) :: ugrad, Nfreq, dz58
    real(wp) :: f_slope, N2, dzsdx, dzsdy, s_b, s_th
    real(wp) :: vtf(im,jm,km)  ! polar-filtered meridional thermal wind [m/s]
    real(wp) :: zb, zt         ! m, bottom and top heights of the shear layer
    real(wp) :: f_depth        ! depth of the shear layer relative to z500-z850
    real(wp) :: zsf(im,jm)     ! terrain height for the slope term [m]
    real(wp) :: kdif(im,jm)    ! macro-diffusivity at cell centres [m2/s]
    ! implicit (ADI) solve of the synoptic-energy (sam) equation
    real(wp) :: kdiss(im,jm), samx(im,jm)
    real(wp) :: gx
    real(wp) :: a(im), b(im), c(im), r(im), sol(im)
    real(wp) :: aj(jm), bj(jm), cj(jm), rj(jm), solj(jm)


    !-----------------------------------------------
    ! synoptic kinetic energy
    !-----------------------------------------------

    ! The meridional thermal wind carries the zonal temperature gradient and the slope term the
    ! zonal terrain gradient, whose grid-scale part is amplified by 1/cos(lat) toward the poles.
    ! The polar filter truncates both to the zonal wavenumbers resolved at the meridional spacing
    ! (smooth_atm), on the layers that can enter the shear.
    vtf = 0._wp
    do k=k850,km
      vtf(:,:,k) = vt3f(:,:,k)
      call polar_filter(vtf(:,:,k), cost)
    enddo
    zsf  = zsa
    call polar_filter(zsf, cost)

    !$omp parallel do collapse(2) private(i, j, k, zb, zt, f_depth, ugrad, dudz, dvdz, Nfreq, dz58, f_slope, N2, ipl, imi, dzsdx, dzsdy, s_b, s_th, iop, i_n, j_n, i_s, j_s)
    do i=1,im
      do j=1,jm

        !-----------------------------------------------
        ! synoptic energy production

        ! proportional to max Eady model baroclinic growth rate (e.g. Hoskins and Valdes 1990)

        ! Shear and stability between a bottom at the 850 hPa layer centre, raised to the top of the
        ! surface layer (zsa+hgams) where that is higher, and a top at the 500 hPa layer centre,
        ! raised to keep at least dz_syn_min above the bottom (high plateaus). Neither the boundary
        ! layer nor the extrapolated sub-surface column enters. Values at zb and zt are interpolated
        ! linearly in geometric height between the layer centres (vesta), consistent with
        ! tp = t + g*z/cp. The thermal wind (u3d) is the full one without polar damping, so the
        ! difference over the range only involves the temperatures in it.
        zb = max(z3(i,j,k850), zsa(i,j)+hgams)
        zt = max(z3(i,j,k500), zb+dz_syn_min)
        dz58 = zt-zb
        N2 = g/(0.5_wp*(zinterp(z3(i,j,:),tp(i,j,:),zb)+zinterp(z3(i,j,:),tp(i,j,:),zt))) &
          * (zinterp(z3(i,j,:),tp(i,j,:),zt)-zinterp(z3(i,j,:),tp(i,j,:),zb))/dz58
        Nfreq = sqrt(max(N2,1.e-6_wp))   ! Brunt-Vaisala frequency, floor 1e-3 1/s
        dudz = (zinterp(z3(i,j,:),ut3f(i,j,:),zt)-zinterp(z3(i,j,:),ut3f(i,j,:),zb))/dz58
        dvdz = (zinterp(z3(i,j,:),vtf(i,j,:),zt)-zinterp(z3(i,j,:),vtf(i,j,:),zb))/dz58
        ugrad = sqrt(dudz**2+dvdz**2)
        ! Depth of the baroclinic layer relative to the free 850-500 hPa range (z850 extrapolated
        ! below ground over high terrain). The production is the growth rate times a reference eddy
        ! energy, and a shallow baroclinic layer holds smaller and weaker eddies (deformation radius
        ! N*H over f, velocity difference shear*H), so the production scales with the layer depth.
        f_depth = min(1._wp, dz58/(z3(i,j,k500)-z3(i,j,k850)))
        ! Sloping lower boundary (Blumsack & Gierasch 1972): terrain rising toward the cold side,
        ! i.e. in the same sense as the isentropes, weakens the baroclinic instability, and in the
        ! Eady model the flow is stable once the terrain slope reaches the isentrope slope (delta>=1).
        ! The cold-side direction follows from the thermal wind, -grad(T) ~ sign(f)*(-dv/dz, du/dz)
        ! (y northward, as in u3d), and the isentrope slope is |f|*|du/dz|/N^2.
        f_slope = 1._wp
        if (c_eke_slope.gt.0._wp) then
          N2 = Nfreq**2
          if (ugrad.gt.0._wp .and. N2.gt.0._wp) then
            ipl = i+1; if (ipl.gt.im) ipl = 1
            imi = i-1; if (imi.lt.1)  imi = im
            dzsdx = (zsf(ipl,j)-zsf(imi,j))/(2._wp*dxtg(j))
            ! meridional neighbours, beyond the pole the cell of the same row half way round
            iop = modulo(i-1+im/2,im) + 1
            i_n = i ; j_n = j-1
            i_s = i ; j_s = j+1
            if (j.eq.1)  then ; i_n = iop ; j_n = 1  ; endif
            if (j.eq.jm) then ; i_s = iop ; j_s = jm ; endif
            dzsdy = (zsf(i_n,j_n)-zsf(i_s,j_s))/(2._wp*dy)   ! northward (j increases southward)
            s_b  = sign(1._wp,sint(j))*(dzsdy*dudz-dzsdx*dvdz)/ugrad
            s_th = 2._wp*omega*abs(sint(j))*ugrad/N2
            f_slope = max(0._wp, 1._wp-c_eke_slope*s_b/s_th)
            f_slope = min(1._wp, f_slope)
          endif
        endif
        synprod(i,j) = c_syn_1 + c_syn_2 * 2._wp*omega*abs(sint(j)) / Nfreq * ugrad * f_slope * f_depth
        synprod(i,j) = max(0._wp,synprod(i,j))

        !-----------------------------------------------
        ! synoptic energy dissipation rate (1/s), linearized about the old sam for the implicit solve

        if (i_syn_diss.eq.1) then
          ! quadratic (cascade/drag) dissipation, (c_syn_3+c_syn_4*cda)*sam^1.5, local balance gives sam ~ synprod^(2/3)
          kdiss(i,j) = (c_syn_3+c_syn_4*cda(i,j))*sqrt(sam(i,j))
        else if (i_syn_diss.eq.2) then
          ! linear damping on a fixed ~1 day time scale, local balance gives sam ~ synprod, which
          ! matches the ERA-Interim seasonal and meridional scaling of EKE with the Eady growth rate
          kdiss(i,j) = (c_syn_3+c_syn_4*cda(i,j))*sqrt(syn_e0)
        endif

      enddo
    enddo
    !$omp end parallel do

    !-----------------------------------------------
    ! synoptic energy evolution - FULLY IMPLICIT (advection + diffusion + semi-implicit dissipation).
    ! Solve (I + dt_atm*(kdiss + Adv_upwind - cdif*Lap)) sam_new = sam + dt_atm*synprod
    ! by ADI: cyclic zonal sweep, then meridional sweep (no flux through the poles). 

    ! zonal implicit sweep (advection + diffusion; periodic -> cyclic tridiagonal)
    !$omp parallel do private(i,j,uef,gx,a,b,c,r,sol)
    do j=1,jm
      do i=1,im
        uef = uadv(i,j)
        gx  = dt_atm*cdif(i,j)/dxt(j)**2
        ! diffusion + upwind advection (imi donor for uef>0, ipl for uef<0)
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

    ! meridional implicit sweep over all rows, per longitude. Diffusion in flux form on the sphere,
    ! (1/cos) d/dy (cos K dS/dy), with K and cos on the cell faces; the polar faces have zero
    ! length, so nothing crosses the pole and the zonal-mean sam is conserved by the diffusion.
    ! Upwind advection in advective form, with no donor beyond the pole.
    !$omp parallel do private(i,j,vef,kn,ks,aj,bj,cj,rj,solj)
    do i=1,im
      do j=1,jm
        vef = vadv(i,j)
        kn = 0._wp
        ks = 0._wp
        if (j.gt.1)  kn = 0.5_wp*(cdif(i,j-1)+cdif(i,j))*cosu(j)
        if (j.lt.jm) ks = 0.5_wp*(cdif(i,j)+cdif(i,j+1))*cosu(j+1)
        aj(j) = -dt_atm*kn/(cost(j)*dy**2)
        cj(j) = -dt_atm*ks/(cost(j)*dy**2)
        if (j.gt.1)  aj(j) = aj(j) + dt_atm*min(vef,0._wp)/dy
        if (j.lt.jm) cj(j) = cj(j) - dt_atm*max(vef,0._wp)/dy
        bj(j) = 1._wp - aj(j) - cj(j)
        rj(j) = samx(i,j)
      enddo
      call tridiag_solve(aj, bj, cj, rj, solj, jm)
      do j=1,jm
        sam(i,j) = max(solj(j),1._wp)
      enddo
    enddo
    !$omp end parallel do


    
    !-----------------------------------------------
    ! Products: macro diffusion, wind, wind stress
    !-----------------------------------------------

    !$omp parallel do collapse(2) private(i,j,n,sqsam,wtau,richardson,f_stab)
    do i=1,im
      do j=1,jm

        sqsam = sqrt(sam(i,j)) 

        ! macro-diffusivity of the synoptic energy itself, proportional to EKE,
        ! K = c*V'^2 (fixed eddy time scale c), as for energy, water and dust
        cdif(i,j) = c_syn_5*sam(i,j)

        do n=1,nm

          ! synoptic surface wind
          synsur(i,j,n) = c_syn_6*sqsam*epsa(i,j,n)
          synsur(i,j,n) = max(synsur(i,j,n),synsurmin)

          ! total surface wind
          wind(i,j,n) = sqrt(us(i,j,n)**2+vs(i,j,n)**2+synsur(i,j,n)**2) 
          if (n.eq.i_ocn) then
            wind(i,j,n) = max(wind(i,j,n),windmin)
          endif

          ! wind stress
          ! The mean stress is rho*cd*<|V| V>, not rho*cd*V_mean*<|V|>: the strong winds blow preferentially along the
          ! mean direction, so the speed and direction fluctuations covary. For isotropic Gaussian fluctuations of
          ! per-component variance s^2 about the mean wind U,
          !   <|V| u> -> 1.5*U*<|V'|>   (U << s) ,     <|V| u> -> U*(U^2 + 1.5 s^2)/U   (U >> s) ,
          ! while wind (the mean speed used for the scalar fluxes) carries synsur ~ <|V'|> (U << s) or synsur^2 ~ s^2 (U >> s).
          ! Both limits are reproduced by an effective speed for the stress
          !   wtau = sqrt(us^2 + vs^2 + c_tau_syn*synsur^2) ,   c_tau_syn = 2.25 (U << s) .. 3 (U >> s) .
          wtau = sqrt(us(i,j,n)**2+vs(i,j,n)**2+c_tau_syn*synsur(i,j,n)**2)
          if (n.eq.i_ocn) then
            wtau = max(wtau,windmin)
          endif
          ! Stability dependence of the momentum drag. 
          ! Same BATS bulk Richardson factor as for the ocean heat and moisture exchange in surface_par_sic.
          f_stab = 1._wp
          if (l_tau_stab .and. frst(i,j,n).gt.0._wp .and. t2(i,j,n).gt.0._wp) then
            richardson = g*z_sfl*(1._wp-tskin(i,j,n)/t2(i,j,n))/max(wind(i,j,n),1._wp)**2
            if (richardson.lt.0._wp) then
              f_stab = 1._wp - f_tau_Ri_unstab*richardson
            else
              f_stab = 1._wp/(1._wp + f_tau_Ri_stab*richardson)
            endif
          endif
          taux(i,j,n) = tau_fac*f_stab*cd(i,j,n)*ra*us(i,j,n)*wtau
          tauy(i,j,n) = tau_fac*f_stab*cd(i,j,n)*ra*vs(i,j,n)*wtau

        enddo

        winda(i,j) = sum(wind(i,j,:)*frst(i,j,:))

        ! Synoptic vertical velocity at ~700 hPa. Continuity gives w' ~ V'*H/L for an eddy of horizontal scale L and depth H, 
        ! and L is the deformation radius N*H/|f|, so the conversion from the eddy speed V'=sqrt(sam) to the eddy vertical velocity scales with |f| 
        wsyn(i,j) = c_syn_7*sqsam*abs(sint(j))/sin45

      enddo
    enddo
    !$omp end parallel do

    !-----------------------------------------------
    ! macro diffusivities for energy, water and dust

    ! proportional to EKE, K = c*V'^2 (mixing length ~ V'*c, i.e. eddy time scale c; Caballero & Hanley 2012).
    kdif = c_diff*sam

    do i=1,im
      imi=i-1
      if (imi.lt.1) imi=im
      do j=1,jm
        diffx(i,j) = 0.5_wp*(kdif(imi,j)+kdif(i,j))
      enddo
    enddo 
    do j=1,jm
      diffx(imc,j) = diffx(1,j)
    enddo

    do j=2,jm
      jmi=j-1
      if (jmi.lt.1) jmi=1
      do i=1,im
        diffy(i,j) = 0.5_wp*(kdif(i,jmi)+kdif(i,j))
      enddo
    enddo
    diffy(:,1)   = 0._wp
    diffy(:,jmc) = 0._wp


    return

  end subroutine synop


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function :  z i n t e r p
  !   Purpose  :  linear interpolation of a column field f given at the layer centre
  !               heights zc (increasing upward) to the height z; constant beyond the ends
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  pure function zinterp(zc, f, z) result(fz)

    implicit none

    real(wp), intent(in) :: zc(:), f(:), z
    real(wp) :: fz
    integer :: k, n

    n = size(zc)
    if (z.le.zc(1)) then
      fz = f(1)
      return
    endif
    do k=1,n-1
      if (z.le.zc(k+1)) then
        fz = f(k) + (f(k+1)-f(k))*(z-zc(k))/(zc(k+1)-zc(k))
        return
      endif
    enddo
    fz = f(n)

  end function zinterp

end module synop_mod


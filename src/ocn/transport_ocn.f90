!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t r a n s p o r t _ o c n _ m o d
!
!  Purpose : time step of ocean tracers
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
module transport_ocn_mod

  use precision, only : wp
  use dim_name, only: dim_lon, dim_lat, dim_depth, dim_time, dim_depth1, dim_lon1, dim_lat1
  use ncio

  use timer, only : year, doy, time_soy_ocn
  use control, only: out_dir
  use climber_grid, only : lon, lat
  use constants, only : g, pi, fcoriolis
  use ocn_grid, only : mask_ocn, mask_c, mask_u, mask_v, mask_w, k1
  use ocn_grid, only : maxi, maxj, maxk, mask_v, mask_w
  use ocn_grid, only : zro, zw, dz, rdza, z2dzg, depth, dx, dxv, rdx, dy, rdy, s, sv
  use ocn_params, only : dt, rho0
  use ocn_params, only : n_tracers_tot, n_tracers_ocn, n_tracers_trans, idx_tracers_trans
  use ocn_params, only : i_diff_dia
  use ocn_params, only : diff_iso, diff_gm
  use ocn_params, only : l_mle, c_mle, mle_fac, tau_mle, h_mle_max, cfl_mle
  use ocn_params, only : diff_dia_ref, diff_dia, diff_dia_zref, diff_dia_min, diff_dia_max
  use ocn_params, only : l_diff_dia_strat, brunt_vaisala_ref, alpha_strat
  use ocn_params, only : slope_max, slope_crit, drho_hbl
  use ocn_params, only : diffx_max, diffy_max
  use ocn_params, only : mlddec, mlddecd, ke_wind_dec, pe_buoy_coeff
  use ocn_params, only : l_downslope

  use advection_mod, only : advection_fct
  use diffusion_mod, only : diffusion_coef, diffusion, slope_taper
  use eos_mod
  use convection_mod, only : convection
  use krausturner_mod
  use downslope_mod, only : downslope

  implicit none

  real(wp), dimension(:,:,:), allocatable :: drho_dx, drho_dy, drho_dz, Ri
  real(wp), dimension(:,:,:), allocatable :: v_gm   !! Gent-McWilliams eddy-induced meridional velocity on the v-grid, diagnostic only [m/s]
  real(wp), dimension(:,:,:), allocatable :: u_mle, v_mle, w_mle   !! mixed-layer eddy (Fox-Kemper) velocities on the u-, v- and w-grid [m/s]
  real(wp), dimension(:,:), allocatable :: hbl     !! surface boundary layer depth for the isoneutral diffusion and GM (positive) [m]
  real(wp), dimension(:,:), allocatable :: sx_bl, sy_bl   !! tapered isoneutral slopes at the base of the boundary layer []
  real(wp), dimension(:,:), allocatable :: z_topo_std_dl  !! standard deviation of bathymetry, for i_diff_dia==4 [m]

  private
  public :: transport, transport_init, transport_coef_update, drho_dx, drho_dy, drho_dz, Ri, v_gm, v_mle

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t r a n s p o r t
  !   Purpose    :  transport of tracers
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine transport(l_tracers_trans, &
                      u,ke_tau,flx_sur,flx_bot,f_ocn,mask_coast,z_ocn_max, &
                      ts,rho,nconv,dconv,kven,dven,dconv_int,hconv_int,conv_pe,fds,zds, &
                      mld,fdx,fdy,fdz,fax,fay,faz,dts_dt_adv,dts_dt_diff, error)

    !$ use omp_lib

    implicit none

    logical, dimension(:), intent(in) :: l_tracers_trans
    real(wp), dimension(1:,0:,0:,:), intent(in) :: u
    real(wp), dimension(:,:), intent(in) :: ke_tau
    real(wp), dimension(:,:,:), intent(in) :: flx_sur
    real(wp), dimension(:,:,:), intent(in) :: flx_bot
    real(wp), dimension(:,:), intent(in) :: f_ocn
    integer, dimension(:,:), intent(in) :: mask_coast
    real(wp), dimension(:,:), intent(in) :: z_ocn_max

    real(wp), dimension(:,:,:), intent(inout) :: rho
    real(wp), dimension(:,:,:,:), intent(inout) :: ts
    integer, dimension(:,:), intent(inout) :: kven, nconv
    real(wp), dimension(:,:), intent(inout) :: dven, dconv
    real(wp), dimension(:,:), intent(inout) :: dconv_int, hconv_int
    real(wp), dimension(:,:), intent(inout) :: conv_pe
    real(wp), dimension(:,:), intent(inout) :: fds
    real(wp), dimension(:,:), intent(inout) :: zds

    real(wp), dimension(:,:), intent(inout) :: mld   ! Kraus-Turner mixed layer depth (negative) [m]
    real(wp), dimension(:,:,:,:), intent(out) :: fdx
    real(wp), dimension(:,:,:,:), intent(out) :: fdy
    real(wp), dimension(:,:,:,:), intent(out) :: fdz
    real(wp), dimension(:,:,:,:), intent(out) :: fax
    real(wp), dimension(:,:,:,:), intent(out) :: fay
    real(wp), dimension(:,:,:,:), intent(out) :: faz
    real(wp), dimension(:,:,:,:), intent(out) :: dts_dt_adv
    real(wp), dimension(:,:,:,:), intent(out) :: dts_dt_diff

    logical, intent(inout) :: error

    integer :: i, j, k, l, ll, ip1, k1_max
    logical :: par_tracers
    integer, allocatable, dimension(:,:) :: mldk
    real(wp), allocatable, dimension(:,:,:,:) :: uadv   ! advection velocity: resolved + mixed-layer eddy (Fox-Kemper)
    real(wp), dimension(maxk) :: rho_old
    real(wp), allocatable, dimension(:,:) :: pe_layer1
    real(wp) :: pe_conv, pe_buoy, e_mix
    real(wp) :: mldtstmp(2), mldrhotmp
    real(wp) :: brunt_vaisala
    real(wp) :: diff_dia_tmp
    real(wp) :: rho1, rho2, rho_c, drhodz, dudz2, tv1
    integer :: kb, knp, nnp, i0, i1, jm, j0
    real(wp) :: slope_x, slope_y, taper, mc_x, yfac
    real(wp), parameter :: eps = 1.e-15_wp

    !$ logical, parameter :: print_omp = .false.
    !$ real(wp) :: time1,time2


    allocate(mldk(maxi,maxj))
    allocate(pe_layer1(maxi,maxj))

    !$ time1 = omp_get_wtime()
    ! before calculating fluxes, calculate energy consumed
    ! or released in mixing surface forcing over top layer. Needed for
    ! mld calculation, especially if mixed layer is <1 cell thick.
    ! NOTE: all energies in this scheme are calculated in units of (Energy/area). Still true?
    !$omp parallel do collapse(2) private(i,j,mldtstmp,mldrhotmp)
    do j=1,maxj
       do i=1,maxi
         if (mask_ocn(i,j).eq.1) then
           ! apply surface fluxes and compute new virtual temperature and salinity
           mldtstmp(1) = ts(i,j,maxk,1)-flx_sur(i,j,1)*dt/dz(maxk)  ! C
           mldtstmp(2) = ts(i,j,maxk,2)-flx_sur(i,j,2)*dt/dz(maxk)  ! psu
           mldtstmp(2) = max(0._wp,mldtstmp(2))
           ! compute new virtual density
           mldrhotmp = eos(mldtstmp(1),mldtstmp(2),zro(maxk))
           ! potential energy change induced by mixing the top layer
           pe_layer1(i,j) = 0.5_wp*g*(mldrhotmp-rho(i,j,maxk))*z2dzg(maxk,maxk) ! J/m2 or kg/s2
         endif
       enddo
    enddo
    !$omp end parallel do
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'transport: mldini ',time2-time1


    ! tracer-independent part of isoneutral diffusion
    !$ time1 = omp_get_wtime()

      ! check density with reference values
      !rho1 = eos(5._wp,0._wp,0._wp)
      !print *,rho1
      !rho1 = eos(5._wp,0._wp,-10000._wp)
      !print *,rho1
      !rho1 = eos(25._wp,0._wp,0._wp)
      !print *,rho1
      !rho1 = eos(25._wp,0._wp,-10000._wp)
      !print *,rho1
      !rho1 = eos(5._wp,35._wp,0._wp)
      !print *,rho1
      !rho1 = eos(5._wp,35._wp,-10000._wp)
      !print *,rho1
      !rho1 = eos(25._wp,35._wp,0._wp)
      !print *,rho1
      !rho1 = eos(25._wp,35._wp,-10000._wp)
      !print *,rho1
      !print *
      !rho1 = eos(3._wp,35.5_wp,-3000._wp)
      !print *,rho1,1041.83267

       ! compute x,y,z-density gradient needed for isoneutral diffusion
       !$omp parallel do collapse(2) &
       !$omp private ( i,j,k ,brunt_vaisala,diff_dia_tmp,ip1,rho1,rho2,rho_c)
       do k=1,maxk
         do j=1,maxj
           do i=1,maxi
             if (mask_w(i,j,K).eq.1) then
               ! get vertical density gradient on w-grid
               rho1 = eos(ts(i,j,k,1),ts(i,j,k,2),zw(K))
               rho2 = eos(ts(i,j,k+1,1),ts(i,j,k+1,2),zw(K))
               drho_dz(i,j,K) = (rho2-rho1)*rdza(K)
               !dudz2 = 0.5_wp * ( &
               !      + (u(1,I-1,j,k+1) - u(1,I-1,j,k))**2 &
               !      + (u(1,I,j,k+1)   - u(1,I,j,k))**2   &
               !      + (u(2,i,J-1,k+1) - u(2,i,J-1,k))**2 &
               !      + (u(2,i,J,k+1)   - u(2,i,J,k))**2 ) * rdza(K)**2
               !drhodz = (rho2-rho1)*rdza(k)
               !brunt_vaisala = sqrt(-g/rho0*min(0._wp,drhodz)) ! Brunt-Vaisala frequency
               !Ri(i,j,K) = brunt_vaisala**2 / dudz2
               if (l_diff_dia_strat) then
                 ! stratification dependent diapycnal diffusivity following Marzeion 2007
                 if (drho_dz(i,j,K).lt.-1.e-12_wp) then
                   brunt_vaisala = sqrt(-g/rho0*drho_dz(i,j,K)) ! Brunt-Vaisala frequency
                   diff_dia_tmp = diff_dia_ref*(brunt_vaisala/brunt_vaisala_ref)**(-alpha_strat)
                   diff_dia_tmp = max(diff_dia_min,diff_dia_tmp)
                   diff_dia(i,j,K) = min(diff_dia_max,diff_dia_tmp)
                 else
                   diff_dia(i,j,K) = diff_dia_max
                 endif
               endif
             endif
             if (mask_c(i,j,k).eq.1) then
               ! central-cell density at zro(k), shared by the zonal and meridional gradients below
               rho_c = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
             endif
             if (mask_u(i,j,k).eq.1) then
               ! zonal density gradient on u-grid
               ip1 = modulo(i,maxi) + 1
               rho2 = eos(ts(ip1,j,k,1),ts(ip1,j,k,2),zro(k))
               drho_dx(i,j,k) = (rho2-rho_c)*rdx(j)
             endif
             if (mask_v(i,j,k).eq.1) then
               ! meridional density gradient on v-grid
               rho2 = eos(ts(i,j+1,k,1),ts(i,j+1,k,2),zro(k))
               drho_dy(i,j,k) = (rho2-rho_c)*rdy
             endif
           enddo
         enddo
       enddo
      !$omp end parallel do

      ! surface boundary layer of the isoneutral diffusion and GM: depth and tapered slopes at its base,
      ! averaged over the 4 corners of the w-point as in the vertical isoneutral flux in diffusion.f90
      !$omp parallel do collapse(2) private(i,j,k,kb,knp,nnp,i0,i1,jm,j0,drhodz,slope_x,slope_y,taper,mc_x,yfac)
      do j=1,maxj
        do i=1,maxi
          hbl(i,j) = 1._wp
          sx_bl(i,j) = 0._wp
          sy_bl(i,j) = 0._wp
          if (mask_ocn(i,j).eq.1) then
            hbl(i,j) = hbl_density(ts(i,j,:,1),ts(i,j,:,2),k1(i,j),drho_hbl)
            ! first w-level below the boundary layer, none if it reaches the bottom
            kb = 0
            do k=maxk-1,k1(i,j),-1
              if (zw(k).lt.-hbl(i,j)) then
                kb = k
                exit
              endif
            enddo
            if (kb.gt.0) then
              drhodz = min(-eps,drho_dz(i,j,kb))
              do knp=0,1
                do nnp=0,1
                  i0 = modulo(i-2+2*nnp, maxi) + 1
                  i1 = modulo(i+nnp-2, maxi) + 1
                  mc_x = real(mask_c(i0,j,kb+knp),wp)
                  slope_x = mc_x * drho_dx(i1,j,kb+knp)/drhodz
                  j0 = j-1+2*nnp
                  jm = max(j+nnp-1, 1)
                  yfac = merge(real(mask_c(i,min(max(j0,1),maxj),kb+knp),wp), 0._wp, j0.ge.1 .and. j0.le.maxj)
                  slope_y = yfac * drho_dy(i,jm,kb+knp)/drhodz
                  taper = slope_taper(slope_x**2+slope_y**2+eps,slope_crit(j,kb))
                  sx_bl(i,j) = sx_bl(i,j) + 0.25_wp*taper*slope_x
                  sy_bl(i,j) = sy_bl(i,j) + 0.25_wp*taper*slope_y
                enddo
              enddo
            endif
          endif
        enddo
      enddo
      !$omp end parallel do

      ! eddy-induced meridional velocity, needed only for the overturning and heat transport diagnostics
      call eddy_induced_velocity

      ! mixed-layer eddy restratification (Fox-Kemper et al. 2008, 2011), from the Kraus-Turner mixed layer of the
      ! previous step and the density gradients above
      if (l_mle) call mle_velocity(mld)

      ! tracer-independent coefficients of the diffusion, shared by all tracers
      call diffusion_coef(diff_iso,diff_dia,drho_dx,drho_dy,drho_dz,slope_crit,hbl,sx_bl,sy_bl)

    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'transport: diffini ',time2-time1

    ! advection velocity: resolved + mixed-layer eddy velocities (divergence-free, zero net transport per column face)
    allocate(uadv(3,0:maxi,0:maxj,maxk))
    uadv = u
    if (l_mle) then
      uadv(1,1:maxi,1:maxj,:) = uadv(1,1:maxi,1:maxj,:) + u_mle
      uadv(2,1:maxi,1:maxj,:) = uadv(2,1:maxi,1:maxj,:) + v_mle
      uadv(3,1:maxi,1:maxj,:) = uadv(3,1:maxi,1:maxj,:) + w_mle
      uadv(1,0,:,:) = uadv(1,maxi,:,:)
    endif

    ! tracer loop; parallel over tracers if there are enough of them to keep the threads busy (e.g. with BGC),
    ! otherwise serial with advection and diffusion parallel over the grid
    !$ time1 = omp_get_wtime()
    ! (two explicit branches rather than an if clause on the directive: with an inactive outer region the grid loops
    ! become nested regions, which start a second thread pool that competes with the other components)
    par_tracers = .false.
    !$ par_tracers = 2*n_tracers_trans.ge.omp_get_max_threads()
    if (par_tracers) then
      !$omp parallel do private(ll)
      do ll=1,n_tracers_trans
        call advdiff_tracer(idx_tracers_trans(ll))
      enddo
      !$omp end parallel do
    else
      do ll=1,n_tracers_trans
        call advdiff_tracer(idx_tracers_trans(ll))
      enddo
    endif
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'transport: advdiff ',time2-time1

    ! convection and mixed layer scheme
    !$ time1 = omp_get_wtime()
    !$omp parallel do collapse(2) private(i,j,k,rho_old,pe_conv,pe_buoy,e_mix,tv1,k1_max)
    do j=1,maxj
      do i=1,maxi
        if (mask_ocn(i,j).eq.1) then
          ! check for NaNs
          if (ts(i,j,maxk,1).ne.ts(i,j,maxk,1) .or. ts(i,j,maxk,2).ne.ts(i,j,maxk,2)) then
            print *,'NaN produced in ocn'
            print *,'i,j',i,j
            print *,'t',ts(i,j,maxk,1)
            print *,'s',ts(i,j,maxk,2)
            error=.true.
          endif
          ! update density
          do k=k1(i,j),maxk
            ts(i,j,k,2) = max(0._wp,ts(i,j,k,2))  ! make sure that salinity is positive
            rho(i,j,k) = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
            ! remember old density needed to calculate PE change. 
            rho_old(k) = rho(i,j,k)
          enddo
          ! convection 
          k1_max = k1(i,j)
          tv1 = 5000._wp
          do k=maxk,1,-1
            if (abs(z_ocn_max(i,j)-zw(k-1)).lt.tv1) then
              k1_max = k
              tv1 = abs(z_ocn_max(i,j)-zw(k-1))
            endif
          enddo
          k1_max = max(k1_max,k1(i,j))
          call convection(l_tracers_trans,ts(i,j,:,:),rho(i,j,:),k1(i,j),k1_max,mask_coast(i,j),nconv(i,j),dconv(i,j),kven(i,j),dven(i,j),dconv_int(i,j),hconv_int(i,j),i,j) 
          ! calculate potential energy released by convection (only layers that are ventilated by the surface!)
          pe_conv = 0._wp
          do k=kven(i,j),maxk 
            pe_conv = pe_conv+0.5_wp*g*(rho(i,j,k)-rho_old(k))*z2dzg(k,k)  ! J/m2 or kg/s2
          enddo
          conv_pe(i,j) = pe_conv   ! J/m2
          ! add energy consumed or released in mixing surface forcing over top layer
          pe_buoy = pe_conv+pe_layer1(i,j)
          ! multiply by efficiency of recycling of potential energy
          if (pe_buoy.gt.0._wp) then
            pe_buoy = pe_buoy*pe_buoy_coeff
          endif
          ! Add wind energy to get total energy available for mixing 
          e_mix = pe_buoy+ke_tau(i,j)*mlddec(maxk)  ! J/m2
          if (e_mix.gt.0._wp) then
            ! Kraus Turner mixed layer scheme. Static instability driven convection has already been done in convection, 
            ! and will differ from standard KT if l_conv_shuffle.eq.T. 
            ! krausturner uses PE released in convection and KE from the wind to deepen the mixed layer further.
            call krausturner(pe_buoy,ke_tau(i,j),k1(i,j), & 
              ts(i,j,1:maxk,1:n_tracers_tot), &
              mld(i,j),mldk(i,j))
          else
            ! Not enough energy even to homogenise first layer. The first layer *is* still homogeneous for all tracers, 
            ! but an mld shallower than the first cell is output as a diagnostic
            mldk(i,j) = maxk
            if (pe_layer1(i,j).lt.0._wp) then
              !          mldtadd = 2._wp*em/(g*zw(K)*(rhol-rhou))
              !          mldt = zw(K)+mldtadd
              mld(i,j) = zw(maxk-1)*(1._wp-e_mix/pe_layer1(i,j))
            else
              mld(i,j) = zw(maxk-1)
            endif
            ! in rare cases the resulting mld can be positive due to termobaricity effects (see email discussion with Edwards and Oliver)
            mld(i,j) = min(0._wp,mld(i,j))
          endif
          !if(mld(i,j).gt.0._wp) then
          !  print *,''
          !  print *,'WARNING: mld>0', i,j
          !  print *,'mld',mld(i,j),mldk(i,j)
          !  print *,'e_mix',e_mix
          !  print *,'pe_layer1',pe_layer1(i,j)
          !  print *,'pe_conv',pe_conv
          !  print *,'pe_buoy',pe_buoy
          !  print *,'ke_tau',ke_tau(i,j)*mlddec(maxk)
          !  !stop 'mld > 0'
          !endif
          ! update density
          do k=k1(i,j),maxk
            rho(i,j,k) = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
          enddo
        endif
      enddo
    enddo
    !$omp end parallel do
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'transport: mld ',time2-time1

    ! downsloping flow of dense shelf water (after convection, so that the shelf column is already homogenised)
    if (l_downslope) then
      call downslope(l_tracers_trans,f_ocn,ts,fds,zds)
      ! update density
      do j=1,maxj
        do i=1,maxi
          if (mask_ocn(i,j).eq.1) then
            do k=k1(i,j),maxk
              rho(i,j,k) = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
            enddo
          endif
        enddo
      enddo
    else
      fds = 0._wp
      zds = 0._wp
    endif

    deallocate(mldk)
    deallocate(uadv)
    deallocate(pe_layer1)

   return

  contains

    subroutine advdiff_tracer(l)

    ! advection and diffusion of tracer l

    implicit none

    integer, intent(in) :: l
    integer :: i, j, k


       ! advection, flux corrected transport scheme following Zalesak 1979
       call advection_fct(f_ocn,uadv(:,1:maxi,1:maxj,1:maxk),ts(:,:,:,l),flx_sur(:,:,l),flx_bot(:,:,l), &
                         fax(:,:,:,l), fay(:,:,:,l), faz(:,:,:,l), dts_dt_adv(:,:,:,l))

       ! diffusion
       call diffusion(f_ocn,ts(:,:,:,l),diff_iso, &
         fdx(:,:,:,l), fdy(:,:,:,l), fdz(:,:,:,l), dts_dt_diff(:,:,:,l))

       ! apply advection and diffusion to tracer field
       !$omp parallel do collapse(2) private(i,j,k)
       do k=1,maxk
         do j=1,maxj
           do i=1,maxi
             if (mask_c(i,j,k).eq.1) ts(i,j,k,l) = ts(i,j,k,l) + dt*(dts_dt_adv(i,j,k,l)+dts_dt_diff(i,j,k,l))
           enddo
         enddo
       enddo
       !$omp end parallel do

       !if (l.gt.1 .and. count(ts(:,:,:,l).lt.0._wp).gt.0) then
       !  print *,'WARNING: negative tracer concentration for tracer #', l-2,' in # points: ',count(ts(:,:,:,l).lt.0._wp)
       !endif

       !do j=1,maxj
       !  do i=1,maxi
       !    if (l.gt.1 .and. any(ts(i,j,k1(i,j):maxk,l).lt.0._wp)) then
       !      print *
       !      print *,'after diffusion+advection'
       !      print *,'WARNING: negative tracer concentration for tracer #', l-2,' in cell (i,j) ',(i,j)
       !      print *,doy
       !      print *,ts(i,j,k1(i,j):maxk,l)
       !      print *,'set tracer concentration to 0'
       !      !ts(i,j,k1(i,j):maxk,l) = max(0._wp,ts(i,j,k1(i,j):maxk,l))
       !    endif
       !  enddo
       !enddo

   return

    end subroutine advdiff_tracer

  end subroutine transport


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  e d d y _ i n d u c e d _ v e l o c i t y
  !   Purpose    :  Gent-McWilliams eddy-induced (bolus) meridional velocity, diagnostic only.
  !                 GM enters the model as a skew flux in diffusion.f90; its advective
  !                 equivalent is v* = d(psi)/dz with the streamfunction psi = diff_gm*taper*slope_y
  !                 on the w-levels (same slope, taper and sign convention as the skew flux),
  !                 zero at the surface, the bottom and closed faces, so that v* integrates to
  !                 zero over each column. Used to put the overturning streamfunctions and the
  !                 overturning/gyre heat transport split on the CMIP definition (resolved plus
  !                 parameterized advection).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine eddy_induced_velocity

    implicit none

    integer :: i, j, k, n
    real(wp) :: drhody, slope_y, taper
    real(wp) :: psi(0:maxk)
    real(wp), parameter :: eps = 1.e-15_wp


    v_gm(:,maxj,:) = 0._wp   ! no v-points on the last row
    !$omp parallel do collapse(2) private(i,j,k,n,drhody,slope_y,taper,psi)
    do j=1,maxj-1
      do i=1,maxi
        psi(:) = 0._wp
        do k=1,maxk-1
          ! interface between levels k and k+1 of the v-column
          if (mask_v(i,j,k).eq.1 .and. mask_v(i,j,k+1).eq.1) then
            drhody = 0.5_wp*(drho_dy(i,j,k)+drho_dy(i,j,k+1))
            ! average over the two adjacent tracer columns
            do n=0,1
              if (mask_w(i,j+n,k).eq.1) then
                if (zw(k).ge.-hbl(i,j+n)) then
                  ! boundary layer: slope linear in depth from its value at the base, no taper (as in diffusion.f90)
                  slope_y = -zw(k)/hbl(i,j+n)*sy_bl(i,j+n)
                  taper = 1._wp
                else
                  slope_y = drhody/min(-eps,drho_dz(i,j+n,k))
                  taper = slope_taper(slope_y**2,slope_crit(j+n,k))
                endif
                psi(k) = psi(k) + 0.5_wp*diff_gm*taper*slope_y   ! m2/s
              endif
            enddo
          endif
        enddo
        do k=1,maxk
          if (mask_v(i,j,k).eq.1) then
            v_gm(i,j,k) = (psi(k)-psi(k-1))/dz(k)   ! m/s
          else
            v_gm(i,j,k) = 0._wp
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine eddy_induced_velocity


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m l e _ v e l o c i t y
  !   Purpose    :  mixed-layer eddy (submesoscale) restratification, Fox-Kemper et al. (2008, 2011).
  !                 On each u- and v-face the eddy-induced streamfunction
  !                   psi = c_mle*mle_fac * H^2 * db/dx_ML / sqrt(f^2+1/tau_mle^2) * mu(z)
  !                 (db/dy_ML on v-faces) is defined on the w-levels inside the mixed layer H, with
  !                   mu(z) = [1-(2z/H+1)^2]*[1+5/21*(2z/H+1)^2],  zero at the surface and at z = -H,
  !                 and the face velocity is d(psi)/dz (z upward), so that it integrates to zero over the
  !                 column of each face (no net transport); w follows from continuity as in velc.f90.
  !                 Light water moves towards the dense side at the top of the mixed layer and dense water
  !                 under it: always restratifying. H is the Kraus-Turner mixed layer depth averaged over the
  !                 two cells of the face and capped at h_mle_max, db/dx_ML = -g/rho0 * thickness-mean of drho_dx over the face levels
  !                 whose centre lies within H (at least the top level). mle_fac is the effective
  !                 Delta_s/L_f of Fox-Kemper et al. (2011) at 5 deg resolution, scaled with the local grid
  !                 spacing Delta_s = min(dx,dy) of the face relative to dy (FK11), so that it shrinks in the
  !                 narrow polar cells. Each face profile is then limited to the Courant number cfl_mle
  !                 (horizontal velocity and the vertical velocity it induces in the adjacent cells).
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine mle_velocity(mld)

    implicit none

    real(wp), dimension(:,:), intent(in) :: mld   ! Kraus-Turner mixed layer depth (negative) [m]

    integer :: i, j, k, ip1, im1
    real(wp) :: h, gsum, hsum, coef, fe, sk, tv, tv1, tv2
    real(wp) :: psi(0:maxk)
    real(wp) :: rtau2

    rtau2 = 1._wp/(tau_mle*86400._wp)**2

    !$omp parallel do collapse(2) private(i,j,k,ip1,h,gsum,hsum,coef,fe,sk,psi)
    do j=1,maxj
      do i=1,maxi

        ! zonal face between i and i+1
        u_mle(i,j,:) = 0._wp
        ip1 = modulo(i,maxi) + 1
        if (mask_u(i,j,maxk).eq.1) then
          h = min(h_mle_max, max(0._wp, -0.5_wp*(mld(i,j)+mld(ip1,j))))
          if (h.gt.0._wp) then
            gsum = 0._wp
            hsum = 0._wp
            do k=maxk,1,-1
              if (mask_u(i,j,k).eq.0) exit
              if (k.lt.maxk .and. zro(k).lt.-h) exit
              gsum = gsum + drho_dx(i,j,k)*dz(k)
              hsum = hsum + dz(k)
            enddo
            fe = sqrt((fcoriolis*s(j))**2 + rtau2)
            coef = c_mle*mle_fac*min(dx(j),dy)/dy * h**2 * (-g/rho0*gsum/hsum) / fe   ! m2/s
            psi(:) = 0._wp
            do k=1,maxk-1   ! w-level zw(k) = interface between levels k and k+1
              if (mask_u(i,j,k).eq.1 .and. mask_u(i,j,k+1).eq.1 .and. zw(k).gt.-h) then
                sk = 2._wp*zw(k)/h + 1._wp
                psi(k) = coef*(1._wp-sk**2)*(1._wp+5._wp/21._wp*sk**2)
              endif
            enddo
            call mle_cfl_limit(psi, dx(j), dx(j))
            do k=1,maxk
              if (mask_u(i,j,k).eq.1) u_mle(i,j,k) = (psi(k)-psi(k-1))/dz(k)   ! m/s
            enddo
          endif
        endif

        ! meridional face between j and j+1
        v_mle(i,j,:) = 0._wp
        if (j.lt.maxj) then
          if (mask_v(i,j,maxk).eq.1) then
            h = min(h_mle_max, max(0._wp, -0.5_wp*(mld(i,j)+mld(i,j+1))))
            if (h.gt.0._wp) then
              gsum = 0._wp
              hsum = 0._wp
              do k=maxk,1,-1
                if (mask_v(i,j,k).eq.0) exit
                if (k.lt.maxk .and. zro(k).lt.-h) exit
                gsum = gsum + drho_dy(i,j,k)*dz(k)
                hsum = hsum + dz(k)
              enddo
              fe = sqrt((fcoriolis*sv(j))**2 + rtau2)
              coef = c_mle*mle_fac*min(dxv(j),dy)/dy * h**2 * (-g/rho0*gsum/hsum) / fe   ! m2/s
              psi(:) = 0._wp
              do k=1,maxk-1
                if (mask_v(i,j,k).eq.1 .and. mask_v(i,j,k+1).eq.1 .and. zw(k).gt.-h) then
                  sk = 2._wp*zw(k)/h + 1._wp
                  psi(k) = coef*(1._wp-sk**2)*(1._wp+5._wp/21._wp*sk**2)
                endif
              enddo
              call mle_cfl_limit(psi, dy, min(dx(j),dx(j+1))*dy/dxv(j))
              do k=1,maxk
                if (mask_v(i,j,k).eq.1) v_mle(i,j,k) = (psi(k)-psi(k-1))/dz(k)   ! m/s
              enddo
            endif
          endif
        endif

      enddo
    enddo
    !$omp end parallel do

    ! vertical velocity from continuity, integrated up from the bottom as in velc.f90
    !$omp parallel do collapse(2) private(i,j,k,im1,tv,tv1,tv2)
    do j=1,maxj
      do i=1,maxi
        w_mle(i,j,:) = 0._wp
        im1 = modulo(i-2,maxi) + 1
        tv = 0._wp
        do k=k1(i,j),maxk
          tv1 = (u_mle(i,j,k) - u_mle(im1,j,k))*dy*dz(k)
          tv2 = v_mle(i,j,k)*dxv(j)*dz(k)
          if (j.gt.1) tv2 = tv2 - v_mle(i,j-1,k)*dxv(j-1)*dz(k)
          w_mle(i,j,k) = tv - (tv1 + tv2)/(dx(j)*dy)
          tv = w_mle(i,j,k)
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine mle_velocity


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  m l e _ c f l _ l i m i t
  !   Purpose    :  scale the MLE streamfunction of one face so that the face velocity d(psi)/dz satisfies
  !                 |u|*dt/dl <= cfl_mle and the vertical velocity it induces in an adjacent cell,
  !                 |psi|/dlc, satisfies |w|*dt/dz <= cfl_mle on both levels around each w-level.
  !                 dl is the grid spacing normal to the face, dlc = cell area / face length.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine mle_cfl_limit(psi, dl, dlc)
    implicit none
    real(wp), intent(inout) :: psi(0:)
    real(wp), intent(in) :: dl, dlc
    integer :: k
    real(wp) :: r
    r = 0._wp
    do k=1,maxk
      r = max(r, abs(psi(k)-psi(k-1))/dz(k)*dt/dl)
    enddo
    do k=1,maxk-1
      r = max(r, abs(psi(k))/dlc*dt/min(dz(k),dz(k+1)))
    enddo
    if (r.gt.cfl_mle) psi = psi*cfl_mle/r
  end subroutine mle_cfl_limit


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function :  h b l _ d e n s i t y
  !   Purpose  :  boundary layer depth from a potential density criterion: depth where potential density
  !               is drho above its value at 10 m, linear interpolation between layer centres, minimum 10 m,
  !               full depth if the column is mixed to the bottom (as the mldst03 diagnostic in ocn_out.f90)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function hbl_density(t_col, s_col, kbot, drho) result(h)

    implicit none

    real(wp), dimension(:), intent(in) :: t_col, s_col
    integer, intent(in) :: kbot
    real(wp), intent(in) :: drho
    real(wp) :: h

    integer :: k, k10
    real(wp) :: rho10, rho_k, rho_ml, z_ml


    h = -zw(kbot-1)   ! mixed to the bottom
    if (kbot.lt.maxk) then
      k10 = maxk-1
      do while (zro(k10).gt.-10._wp .and. k10.gt.kbot)
        k10 = k10-1
      enddo
      rho_k = eos(t_col(k10),s_col(k10),0._wp)
      if (zro(k10).gt.-10._wp) then
        rho10 = rho_k
      else
        rho_ml = eos(t_col(k10+1),s_col(k10+1),0._wp)
        rho10 = rho_ml + (rho_k-rho_ml)*(-10._wp-zro(k10+1))/(zro(k10)-zro(k10+1))
      endif
      z_ml = -10._wp
      rho_ml = rho10
      do k=k10,kbot,-1
        rho_k = eos(t_col(k),s_col(k),0._wp)
        if (rho_k-rho10.ge.drho) then
          h = -(z_ml + (zro(k)-z_ml)*(rho10+drho-rho_ml)/(rho_k-rho_ml))
          exit
        endif
        z_ml = zro(k)
        rho_ml = rho_k
      enddo
      h = max(h,10._wp)
    endif

  end function hbl_density


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t r a n s p o r t _ i n i t
  !   Purpose    :  initialize transport
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine transport_init

    implicit none

    character (len=256) :: fnm

    allocate(slope_crit(maxj,maxk))
    allocate(mlddec(maxk))
    allocate(mlddecd(maxk))
    allocate(diff_dia(maxi,maxj,maxk))
    allocate(diffx_max(maxj))
    allocate(diffy_max(maxj))

    allocate(drho_dx(maxi,maxj,maxk))
    allocate(drho_dy(maxi,maxj,maxk))
    allocate(drho_dz(maxi,maxj,maxk))
    allocate(Ri(maxi,maxj,maxk))
    allocate(v_gm(maxi,maxj,maxk))
    allocate(u_mle(maxi,maxj,maxk), source=0._wp)
    allocate(v_mle(maxi,maxj,maxk), source=0._wp)
    allocate(w_mle(maxi,maxj,maxk), source=0._wp)
    allocate(hbl(maxi,maxj), source=1._wp)
    allocate(sx_bl(maxi,maxj), source=0._wp)
    allocate(sy_bl(maxi,maxj), source=0._wp)

    drho_dz = 0._wp
    drho_dx = 0._wp
    drho_dy = 0._wp
    Ri = 0._wp
    v_gm = 0._wp

    if (i_diff_dia.eq.4) then
      ! read standard deviation of bathymetry
      allocate(z_topo_std_dl(maxi,maxj))
      fnm = "input/RTopo-2.0.1_2min_std_30.nc"
      call nc_read(fnm,"topo_r",z_topo_std_dl)
    endif

    call transport_coef_update

    return

  end subroutine transport_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  t r a n s p o r t _ c o e f _ u p d a t e
  ! Purpose    :  coefficients that depend on the vertical grid (zro, zw, dz, k1);
  !               recomputed after every ocean grid update, so that they follow the
  !               current grid both in a continuous run and after a restart
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine transport_coef_update

    implicit none

    integer :: i, j, k
    real(wp) :: r, K_b, h0

    ! mld scheme - calculate wind decay efficiency
    do k=maxk,1,-1
       mlddec(k) = exp(zro(k)/ke_wind_dec)
       if (k.lt.maxk) then
         if (mlddec(k+1).gt.0._wp) then
           mlddecd(k) = mlddec(k)/mlddec(k+1) 
         else
           mlddecd(k) = 0._wp
         endif
       else
          mlddecd(maxk) = mlddec(maxk)
       endif
    enddo


    !******************************************************
    ! define isopycnal and diapycnal diffusivities

    if (i_diff_dia.eq.0) then
      diff_dia = diff_dia_min
    else if (i_diff_dia.eq.1) then
      ! Use Bryan & Lewis (1979) type profile
      do i=1,maxi
        do j=1,maxj
          do k=1,maxk
            diff_dia(i,j,K) = 8.e-5 + 2.e-4/pi*atan(-2e-3*(zw(K)+1000.))
          enddo
        enddo
      enddo
    else if (i_diff_dia.eq.2) then
      ! Use Bryan & Lewis (1979) type profile as used in CLIMBER-2
      do i=1,maxi
        do j=1,maxj
          do k=1,maxk
            diff_dia(i,j,K) = 8e-5 + 1.05e-4_wp/pi*atan(-4.5e-3_wp*(zw(K)+2500._wp))
          enddo
        enddo
      enddo
    else if (i_diff_dia.eq.3) then
      ! Use Bryan & Lewis (1979) type profile
      do i=1,maxi
        do j=1,maxj
          do k=1,maxk
            diff_dia(i,j,K) = diff_dia_min + (atan(-(zw(K)+diff_dia_zref)/1000._wp)-atan(-diff_dia_zref/1000._wp)) &
              / ((atan(-(-5000._wp+diff_dia_zref)/1000._wp)-atan(-diff_dia_zref/1000._wp))) &
              * (diff_dia_max-diff_dia_min)
          enddo
        enddo
      enddo
    else if (i_diff_dia.eq.4) then
      ! vertical profile dependent on roughness at the ocean floor, Decloedt & Luther 2010
      do i=1,maxi
        do j=1,maxj
          r = z_topo_std_dl(i,j)
          if (r.lt.830._wp) then
            K_b = 1.87e-5_wp*exp(3.e-8*r**3-5.8e-5*r**2+0.0325*r)
          else
            K_b = 1.8e-3_wp
          endif
          if (r.lt.540._wp) then
            h0 = -2.9e-6*r**3+0.0046*r**2-2.5896*r+670._wp
          else
            h0 = 170._wp
          endif
          do k=1,maxk
            if (k1(i,j).le.k) then
              diff_dia(i,j,K) = diff_dia_min + K_b*(1._wp+(zw(K)-zw(k1(i,j)))/h0)**(-2)
            else
              diff_dia(i,j,K) = diff_dia_min + K_b
            endif
          enddo
        enddo
      enddo
    endif

        
    ! check that diff_dia is not negative
    if (minval(diff_dia).lt.0._wp) stop 'ERROR: negative diapycnal diffusivity'

    ! stability criterion for isoneutral diffusion, eq. A.2 in Gerdes 1991 (eq. C2 in Griffies 1998?)
    do j=1,maxj
      do k=1,maxk
        ! A33 is implicit: only the limits of the explicit off-diagonal terms remain
        slope_crit(j,k) = min( dx(j)*dz(k)/(4._wp*diff_iso*dt), dy*dz(k)/(4._wp*diff_iso*dt) )
        slope_crit(j,k) = min(slope_max,slope_crit(j,k))
      enddo
    enddo

    ! maximum horizontal diffusivity (theoretical limit is 0.125*dx^2/dt). 
    do j=1,maxj
      diffx_max(j) = 0.1_wp*dx(j)**2/dt
      diffy_max(j) = 0.1_wp*dy**2/dt
    enddo

    
    return

  end subroutine transport_coef_update

end module transport_ocn_mod

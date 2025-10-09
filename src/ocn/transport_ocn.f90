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

  use timer, only : year, doy, time_soy_ocn, sec_day
  use control, only: out_dir
  use climber_grid, only : lon, lat
  use constants, only : g, pi
  use ocn_grid, only : mask_ocn, mask_c, mask_u, mask_v, mask_w, k1, k1_pot
  use ocn_grid, only : maxi, maxj, maxk
  use ocn_grid, only : zro, zw, dz, rdza, z2dzg, depth, dx, dxv, rdx, dy, rdy
  use ocn_params, only : dt, rho0
  use ocn_params, only : n_tracers_tot, n_tracers_ocn, n_tracers_trans, idx_tracers_trans
  use ocn_params, only : i_advection
  use ocn_params, only : i_diff, i_diff_dia, l_diff33_impl
  use ocn_params, only : diff_iso
  use ocn_params, only : diff_dia_ref, diff_dia, diff_dia_bgc, diff_dia_zref, diff_dia_min, diff_dia_bgc_min, diff_dia_max
  use ocn_params, only : l_diff_dia_strat, brunt_vaisala_ref, alpha_strat
  use ocn_params, only : slope_max, slope_crit
  use ocn_params, only : diffx_max, diffy_max
  use ocn_params, only : l_mld, mlddec, mlddecd, ke_wind_dec, pe_buoy_coeff

  use advection_mod, only : advection_upstream, advection_fct
  use diffusion_mod, only : diffusion, diffusion_33
  use eos_mod
  use convection_mod, only : convection
  use krausturner_mod

  implicit none

  real(wp), dimension(:,:,:), allocatable :: drho_dx, drho_dy, drho_dz, Ri, slope2_w

  private
  public :: transport, transport_init, drho_dx, drho_dy, drho_dz, Ri

contains


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t r a n s p o r t
  !   Purpose    :  transport of tracers
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine transport(l_tracers_trans,l_tracer_dic,l_tracers_isodiff,l_large_vol_change, &
                      u,ke_tau,flx_sur,flx_bot,f_ocn,mask_coast,z_ocn_max, &
                      ts,rho,nconv,dconv,kven,dven,conv_pe, &
                      mld,fdx,fdy,fdz,fax,fay,faz,dts_dt_adv,dts_dt_diff, error)

    !$ use omp_lib

    implicit none

    logical, dimension(:), intent(in) :: l_tracers_trans
    logical, dimension(:), intent(in) :: l_tracer_dic
    logical, dimension(:), intent(in) :: l_tracers_isodiff
    logical, intent(in) :: l_large_vol_change
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
    real(wp), dimension(:,:), intent(inout) :: conv_pe

    real(wp), dimension(:,:), intent(out) :: mld
    real(wp), dimension(:,:,:,:), intent(out) :: fdx
    real(wp), dimension(:,:,:,:), intent(out) :: fdy
    real(wp), dimension(:,:,:,:), intent(out) :: fdz
    real(wp), dimension(:,:,:,:), intent(out) :: fax
    real(wp), dimension(:,:,:,:), intent(out) :: fay
    real(wp), dimension(:,:,:,:), intent(out) :: faz
    real(wp), dimension(:,:,:,:), intent(out) :: dts_dt_adv
    real(wp), dimension(:,:,:,:), intent(out) :: dts_dt_diff

    logical, intent(inout) :: error

    integer :: i, j, k, l, ll, ip1, idiff, k1_max
    integer, allocatable, dimension(:,:) :: mldk
    real(wp), dimension(maxk) :: rho_old
    real(wp), allocatable, dimension(:,:) :: pe_layer1
    real(wp) :: pe_conv, pe_buoy, e_mix
    real(wp) :: mldtstmp(2), mldrhotmp
    real(wp) :: brunt_vaisala
    real(wp) :: diff_dia_tmp
    real(wp) :: rho1, rho2, drhodz, dudz2, tv1

    !$ logical, parameter :: print_omp = .false.
    !$ real(wp) :: time1,time2


    allocate(mldk(maxi,maxj))
    allocate(pe_layer1(maxi,maxj))

    if (l_mld) then
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
    endif


    ! tracer-independent part of isoneutral diffusion
    if (i_diff.eq.1) then
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
       !$omp private ( i,j,k ,brunt_vaisala,diff_dia_tmp,ip1,rho1,rho2)
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
                 if (drho_dz(i,j,k).lt.-1.e-12_wp) then
                   !rho1 = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
                   !rho2 = eos(ts(i,j,k+1,1),ts(i,j,k+1,2),zro(k+1))
                   drhodz = (rho2-rho1)*rdza(k)
                   brunt_vaisala = sqrt(-g/rho0*drhodz) ! Brunt-Vaisala frequency
                   diff_dia_tmp = diff_dia_ref*(brunt_vaisala/brunt_vaisala_ref)**(-alpha_strat)                   
                   diff_dia_tmp = max(diff_dia_min,diff_dia_tmp)
                   diff_dia(i,j,K) = diff_dia_tmp !min(diff_dia_max,diff_dia_tmp)
                   diff_dia_bgc(i,j,K) = max(diff_dia_bgc_min,diff_dia(i,j,K)) 
                 else
                   diff_dia(i,j,K) = diff_dia_max
                   diff_dia_bgc(i,j,K) = diff_dia_max
                 endif
               endif
             endif
             if (mask_u(i,j,k).eq.1) then
               ! zonal density gradient on u-grid
               ip1 = modulo(i,maxi) + 1
               rho1 = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
               rho2 = eos(ts(ip1,j,k,1),ts(ip1,j,k,2),zro(k))
               drho_dx(i,j,k) = (rho2-rho1)*rdx(j)
             endif
             if (mask_v(i,j,k).eq.1) then
               ! meridional density gradient on v-grid
               rho1 = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
               rho2 = eos(ts(i,j+1,k,1),ts(i,j+1,k,2),zro(k))
               drho_dy(i,j,k) = (rho2-rho1)*rdy
             endif
           enddo
         enddo
       enddo
      !$omp end parallel do
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'transport: diffini ',time2-time1
    endif

    ! tracer loop
    !$ time1 = omp_get_wtime()
    !$omp parallel do private ( ll, l , idiff)
    do ll=1,n_tracers_trans
      l = idx_tracers_trans(ll)
       !!$ print *,l,ll,omp_get_thread_num(),'/',omp_get_num_threads()

       ! advection
       if (i_advection.eq.1) then
          ! Fiadeiro and Veronis 1977 weighted upstream/centred differences scheme (standard in Goldstein)
          call advection_upstream(f_ocn,u(:,1:maxi,1:maxj,1:maxk),ts(:,:,:,l),flx_sur(:,:,l),flx_bot(:,:,l), &
                                 fax(:,:,:,l), fay(:,:,:,l), faz(:,:,:,l), dts_dt_adv(:,:,:,l))
       elseif (i_advection.eq.2) then
          ! Flux corrected transport advection scheme following Zalesak 1979
          call advection_fct(f_ocn,u(:,1:maxi,1:maxj,1:maxk),ts(:,:,:,l),flx_sur(:,:,l),flx_bot(:,:,l), &
                            fax(:,:,:,l), fay(:,:,:,l), faz(:,:,:,l), dts_dt_adv(:,:,:,l))
       endif

       if (i_diff.eq.1 .and. l_tracers_isodiff(l) .and. .not.l_large_vol_change) then
         idiff = 1      ! use isopycnal diffusion
       else
         idiff = 0      ! use horizontal diffusion
       endif

       ! diffusion
!       if (ll.eq.1 .or. ll.eq.2) then
!         ! active tracers temperature and salinity
!         call diffusion(idiff,f_ocn,ts(:,:,:,l),diff_iso,diff_dia,drho_dx,drho_dy,drho_dz,slope_crit, &
!           fdx(:,:,:,l), fdy(:,:,:,l), fdz(:,:,:,l), slope2_w, dts_dt_diff(:,:,:,l),l)
!       else
!         ! passive tracers 
!         call diffusion(idiff,f_ocn,ts(:,:,:,l),diff_iso,diff_dia_bgc,drho_dx,drho_dy,drho_dz,slope_crit, &
!           fdx(:,:,:,l), fdy(:,:,:,l), fdz(:,:,:,l), slope2_w, dts_dt_diff(:,:,:,l),l)
!       endif

       if (l_tracer_dic(l)) then
         ! DIC tracers 
         call diffusion(idiff,f_ocn,ts(:,:,:,l),diff_iso,diff_dia_bgc,drho_dx,drho_dy,drho_dz,slope_crit, &
           fdx(:,:,:,l), fdy(:,:,:,l), fdz(:,:,:,l), slope2_w, dts_dt_diff(:,:,:,l))
       else
         call diffusion(idiff,f_ocn,ts(:,:,:,l),diff_iso,diff_dia,drho_dx,drho_dy,drho_dz,slope_crit, &
           fdx(:,:,:,l), fdy(:,:,:,l), fdz(:,:,:,l), slope2_w, dts_dt_diff(:,:,:,l))
       endif

       ! apply advection and diffusion to tracer field
       where (mask_c.eq.1)
         ts(:,:,:,l) = ts(:,:,:,l) + dt*(dts_dt_adv(:,:,:,l)+dts_dt_diff(:,:,:,l))
       endwhere

       if (l_diff33_impl) then
         ! treat A33 term of diffusion matrix implicitely
         call diffusion_33(idiff,ts(:,:,:,l),diff_iso,diff_dia,slope2_w,slope_max)
       endif

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

    enddo ! tracer loop
    !$omp end parallel do
    !$ time2 = omp_get_wtime()
    !$ if(print_omp) print *,'transport: advdiff ',time2-time1

    ! convection and mixed layer scheme
    if (l_mld) then
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
            call convection(l_tracers_trans,ts(i,j,:,:),rho(i,j,:),k1(i,j),k1_max,mask_coast(i,j),nconv(i,j),dconv(i,j),kven(i,j),dven(i,j),i,j) 
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
              call krausturner(l_tracers_trans,pe_buoy,ke_tau(i,j),k1(i,j), & 
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

    else

      ! Not applying mixed layer scheme. Just call convection scheme
      !$omp parallel do private(i,j)
      do j=1,maxj
        do i=1,maxi
          if (mask_ocn(i,j).eq.1) then
            ! update density
            do k=k1(i,j),maxk
              ts(i,j,k,2) = max(0._wp,ts(i,j,k,2))  ! make sure that salinity is positive
              rho(i,j,k) = eos(ts(i,j,k,1),ts(i,j,k,2),zro(k))
            enddo
            call convection(l_tracers_trans,ts(:,i,j,:),rho(i,j,:),k1(i,j),k1(i,j),mask_coast(i,j),nconv(i,j),dconv(i,j),kven(i,j),dven(i,j),i,j) 
          endif
        enddo
      enddo
      !$omp end parallel do

    endif

    deallocate(mldk)
    deallocate(pe_layer1)

   return

  end subroutine transport


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t r a n s p o r t _ i n i t
  !   Purpose    :  initialize transport
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine transport_init

    implicit none

    integer :: i, j, k
    real(wp) :: r, K_b, h0
    character (len=256) :: fnm

    real(wp), dimension(maxi,maxj) :: z_topo_std

    allocate(slope_crit(maxj,maxk))
    allocate(mlddec(maxk))
    allocate(mlddecd(maxk))
    allocate(diff_dia(maxi,maxj,maxk))
    allocate(diff_dia_bgc(maxi,maxj,maxk))
    allocate(diffx_max(maxj))
    allocate(diffy_max(maxj))

    allocate(drho_dx(maxi,maxj,maxk))
    allocate(drho_dy(maxi,maxj,maxk))
    allocate(drho_dz(maxi,maxj,maxk))
    allocate(Ri(maxi,maxj,maxk))
    allocate(slope2_w(maxi,maxj,maxk))

    drho_dz = 0._wp
    drho_dx = 0._wp
    drho_dy = 0._wp
    Ri = 0._wp
    slope2_w= 0._wp
  
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
      diff_dia_bgc = diff_dia_bgc_min
    else if (i_diff_dia.eq.1) then
      ! Use Bryan & Lewis (1979) type profile
      do i=1,maxi
        do j=1,maxj
          do k=1,maxk
            diff_dia(i,j,K) = 8.e-5 + 2.e-4/pi*atan(-2e-3*(zw(K)+1000.))
            diff_dia_bgc(i,j,K) = max(diff_dia_bgc_min,diff_dia(i,j,K))
          enddo
        enddo
      enddo
    else if (i_diff_dia.eq.2) then
      ! Use Bryan & Lewis (1979) type profile as used in CLIMBER-2
      do i=1,maxi
        do j=1,maxj
          do k=1,maxk
            diff_dia(i,j,K) = 8e-5 + 1.05e-4_wp/pi*atan(-4.5e-3_wp*(zw(K)+2500._wp))
            diff_dia_bgc(i,j,K) = max(diff_dia_bgc_min,diff_dia(i,j,K))
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
            diff_dia_bgc(i,j,K) = max(diff_dia_bgc_min,diff_dia(i,j,K))
          enddo
        enddo
      enddo
    else if (i_diff_dia.eq.4) then
      ! vertical profile dependent on roughness at the ocean floor, Decloedt & Luther 2010

      ! read standard deviation of bathymetry
      fnm = "input/RTopo-2.0.1_2min_std_30.nc"
      call nc_read(fnm,"topo_r",z_topo_std)

      do i=1,maxi
        do j=1,maxj
          r = z_topo_std(i,j)
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
            if (k1_pot(i,j).le.k) then
              diff_dia(i,j,K) = diff_dia_min + K_b*(1._wp+(zw(K)-zw(k1_pot(i,j)))/h0)**(-2)
            else
              diff_dia(i,j,K) = diff_dia_min + K_b
            endif
            diff_dia_bgc(i,j,K) = max(diff_dia_bgc_min,diff_dia(i,j,K))
          enddo
        enddo
      enddo
    endif

        
    ! check that diff_dia is not negative
    if (minval(diff_dia).lt.0._wp) stop 'ERROR: negative diapycnal diffusivity'

    ! stability criterion for isoneutral diffusion, eq. A.2 in Gerdes 1991 (eq. C2 in Griffies 1998?)
    do j=1,maxj
      do k=1,maxk
        if (l_diff33_impl) then
          slope_crit(j,k) = min( dx(j)*dz(k)/(4._wp*diff_iso*dt), dy*dz(k)/(4._wp*diff_iso*dt) )
        else
          slope_crit(j,k) = min( dz(k)/sqrt(4._wp*diff_iso*dt), dx(j)*dz(k)/(4._wp*diff_iso*dt), dy*dz(k)/(4._wp*diff_iso*dt) )
        endif
        slope_crit(j,k) = min(slope_max,slope_crit(j,k))
      enddo
    enddo

    ! maximum horizontal diffusivity (theoretical limit is 0.125*dx^2/dt). 
    do j=1,maxj
      diffx_max(j) = 0.1_wp*dx(j)**2/sec_day
      diffy_max(j) = 0.1_wp*dy**2/sec_day
    enddo

    
    return

  end subroutine transport_init

end module transport_ocn_mod

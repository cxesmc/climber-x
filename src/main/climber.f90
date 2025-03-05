!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Program : c l i m b e r
!
!  Purpose : main climber program
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit, Andrey Ganopolski and
!                         Alexander Robinson
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
program climber

  use precision, only : wp, dp, sp
  use timer, only : timer_init, timer_update
  use timer, only : step, nstep, year_ini, year, year_now, mon, time_soy, time_eoy, time_eoy_smb, time_soy_smb, time_soy_bmb, time_eoy_bmb
  use timer, only : time_call_atm, time_call_ocn, time_call_lnd, time_call_sic, time_call_smb, time_call_bmb, time_call_ice, time_call_geo, time_call_clim
  use timer, only : time_call_bgc, doy
  use timer, only : time_spinup_cc_1, time_spinup_cc_2
  use timer, only : time_out_ice 
  use timer, only : time_write_restart
  use timer, only : n_year_ice
  use control, only: in_dir, out_dir, control_load, args, l_debug_main_loop, l_write_timer
  use control, only: flag_atm, flag_co2, flag_ch4, flag_ocn, flag_lnd, flag_sic, flag_smb, flag_bmb, flag_bgc, flag_geo, ifake_geo, flag_lakes
  use control, only: flag_ice, ice_model_name, ice_domain_name, n_ice_domain, ice_restart
  use control, only: l_aquaplanet
  use control, only : l_spinup_cc, l_daily_input_save_ocn, l_daily_input_save_bgc
  use coord, only : grid_class

  use climber_grid, only: climber_grid_init, ice_grid_init

  use coupler, only: cmn_class, cmn_init, &
    &                atm_to_cmn, cmn_to_atm, &
    &                cmn_to_bgc, bgc_to_cmn, &
    &                sic_to_cmn, cmn_to_sic, &
    &                lnd_to_cmn, cmn_to_lnd, &
    &                ocn_to_cmn, cmn_to_ocn, &
    &                smb_to_cmn, cmn_to_smb, &
    &                ice_to_smb, ice_to_cmn, smb_to_ice, &
    &                bmb_to_cmn, cmn_to_bmb, &
    &                ice_to_bmb, bmb_to_ice, &
    &                cmn_to_geo, geo_to_cmn, &
    &                ice_to_geo, bnd_to_geo, geo_to_ice, &
    &                geo_to_smb, geo_to_bmb, &
    &                co2_to_cmn, cmn_to_co2, &
    &                ch4_to_cmn, cmn_to_ch4, &
    &                bnd_to_cmn, &
    &                lakes_update, runoff_merge, runoff_to_ocn, &
    &                aquaplanet, aqua_init, aqua_end

  use cmn_out

  use atm_model, only: atm_init, atm_update, atm_end, atm_write_restart
  use atm_def, only: atm_class
  use atm_out, only: atm_diag, atm_diag_init

  use ocn_model, only: ocn_init, ocn_update, ocn_end, ocn_write_restart
  use ocn_def, only: ocn_class
  use ocn_out, only: ocn_diag, ocn_diag_init
  use ocn_params, only: n_tracers_tot, n_tracers_bgc, n_tracers_ocn, n_tracers_trans, idx_tracers_trans

  use bgc_model, only: bgc_ini, bgc_update, bgc_end, bgc_write_restart
  use bgc_def, only : bgc_class, bgc_ntra
  use bgc_out, only: bgc_diag, bgc_diag_init
  use bgc_params, only : l_spinup_bgc, l_spinup_sed

  use sic_model, only: sic_init, sic_update, sic_end, sic_write_restart
  use sic_def, only: sic_class
  use sic_out, only: sic_diag, sic_diag_init

  use lnd_model, only: lnd_init, lnd_update_wrapper, lnd_end, lnd_write_restart
  use lnd_def, only: lnd_class
  use lnd_out, only: lnd_diag, lnd_diag_init

  use ice_model, only : ice_init_domains, ice_init, ice_update, ice_write_restart
  use ice_def, only : ice_class

  use smb_model, only: smb_init, smb_update, smb_end, smb_write_restart
  use smb_def, only : smb_in_class, smb_class 
  use smb_out, only: smb_diag, smb_diag_init

  use bmb_model, only: bmb_init, bmb_update, bmb_end, bmb_write_restart
  use bmb_def, only : bmb_class
  use bmb_out, only: bmb_diag, bmb_diag_init

  use co2_model, only: co2_init, co2_update, co2_end, co2_write_restart
  use co2_def, only : co2_class
  use co2_out, only : co2_diag, co2_diag_init

  use ch4_model, only: ch4_init, ch4_update, ch4_end, ch4_write_restart
  use ch4_def, only: ch4_class
  use ch4_out, only : ch4_diag, ch4_diag_init

  use bnd_mod, only: bnd_init, bnd_update
  use bnd_mod, only: bnd_class

  use geo_mod, only: geo_init, geo_update, geo_end, geo_write_restart
  use geo_def, only: geo_class
  use geo_out, only: geo_diag, geo_diag_init

  use gitversion, only: git_commit_hash

  !$  use omp_lib

  implicit none

  type(cmn_class) :: cmn
  type(atm_class) :: atm
  type(ocn_class) :: ocn
  type(bgc_class) :: bgc
  type(sic_class) :: sic
  type(lnd_class) :: lnd
  type(ice_class), allocatable :: ice(:)
  type(smb_in_class) :: smb_in
  type(smb_class), allocatable :: smb(:) 
  type(bmb_class), allocatable :: bmb(:) 
  type(co2_class) :: co2
  type(ch4_class) :: ch4
  type(geo_class) :: geo
  type(bnd_class) :: bnd
  type(grid_class), allocatable :: ice_grid(:)

  integer :: n, l
  logical :: l_error
  !$ real(dp) :: time_ini, time_end
  !$ real(dp) :: time_atm
  !$ real(dp) :: time_ocn
  !$ real(dp) :: time_sic
  !$ real(dp) :: time_lnd
  !$ real(dp) :: time_bgc
  !$ real(dp) :: time_bnd
  !$ real(dp) :: time_geo
  !$ real(dp) :: time_ice
  !$ real(dp) :: time_smb
  !$ real(dp) :: time_bmb

  write(*,*) 'climber-x git_commit_hash', git_commit_hash
   
    ! ==============================
    ! Initialize model components
    ! ==============================

    ! input directory
    in_dir = "input/"

    ! get output folder name
    call args(out_dir)

    ! read control namelist
    call control_load(trim(out_dir)//"/control.nml")

    !$OMP PARALLEL
    !$OMP MASTER
    !$ print *,'Using ',omp_get_num_threads(),' OpenMP threads.'
    !$OMP END MASTER
    !$OMP END PARALLEL

    ! initialize timer
    call timer_init(trim(out_dir)//"/control.nml")

    ! initialize restart directories
    call restart_init(trim(out_dir))

    ! initialize climber grid
    call climber_grid_init

    ! initialize coupler
    call cmn_init(cmn)

    ! initalize boundary conditions
    call bnd_init(cmn%grid, bnd)

    ! Initialise geography 
    call geo_init(geo, &
      bnd%geo%grid, bnd%geo%z_bed, bnd%geo%z_bed_ref, bnd%ice%grid, bnd%ice%h_ice)  ! used only if not from restart
    call geo_diag_init(geo)

    ! define ice grid object for all domains, if required
    if (flag_ice .or. flag_smb .or. flag_bmb) then
      allocate(ice_grid(n_ice_domain))
      do n=1,n_ice_domain
        call ice_grid_init(ice_domain_name(n), ice_grid(n))
      enddo
    endif

    ! Initialize ice model
    if (flag_ice) then
      allocate(ice(n_ice_domain))
      call ice_init_domains(n_ice_domain, ice_model_name)
      ! loop over ice domains
      do n=1,n_ice_domain
        call ice_init(ice(n),n,ice_model_name,real(year_ini,wp),ice_restart, ice_grid(n), cmn%grid, geo%hires%grid, &
          geo%hires%z_bed, geo%hires%z_bed_rel, geo%hires%h_ice, geo%hires%q_geo, geo%hires%h_sed)  
        where (ice(n)%grid_ice_to_cmn%ncells.gt.0)
          cmn%mask_ice = 1
        endwhere
      enddo
      call ice_to_geo(n_ice_domain,ice,geo)
    endif

    ! update geography
    call geo_update(geo)
    call geo_diag(geo)

    if (flag_ice) then
      call geo_to_ice(n_ice_domain,geo,ice)
    endif

    ! Initialise atmosphere model
    if (flag_atm) then
       if (l_aquaplanet) then
         call aqua_init(cmn)
         geo%f_lnd = 0._wp
         geo%z_sur = 0._wp
       endif
       call atm_init(atm,geo%f_lnd,geo%z_sur,geo%z_sur_std,bnd%co2) 
       call atm_diag_init
    endif
    
    ! Initialise ocean model
    if (flag_ocn) then
       ! initialize tracer number based on biogeochemistry
       if (flag_bgc) then
         n_tracers_bgc = BGC_NTRA
       else
         n_tracers_bgc = 0
       endif
       call ocn_init(ocn,geo%f_ocn,geo%z_bed,geo%z_ocn_max_q,geo%mask_coast2,geo%ocn_vol_tot,geo%A_bering,l_daily_input_save_ocn)
       call ocn_diag_init
    endif

    ! initialize ocean biogeochemistry model
    if (flag_bgc) then
      call bgc_ini(bgc, ocn%grid%ni, ocn%grid%nj, ocn%grid%lon, ocn%grid%lat, &
        ocn%grid%dz(ocn%grid%nk:1:-1), -ocn%grid%zro(ocn%grid%nk:1:-1), -ocn%grid%zw(ocn%grid%nk:0:-1), -ocn%grid%k1+ocn%grid%nk+1, &
        ocn%grid%ocn_area,l_daily_input_save_bgc)
      call bgc_to_cmn(bgc, cmn,ocn) ! initialize tracers in ocean module
      call bgc_diag_init(ocn%grid%ni, ocn%grid%nj, ocn%grid%nk, ocn%grid%zro)
      ! flag for bgc tracers transport
      ocn%l_tracers_trans(n_tracers_ocn+1:) = bgc%l_trans_tracers
      ocn%l_tracer_dic(n_tracers_ocn+1:) = bgc%l_tracer_dic
      ocn%l_tracers_isodiff(n_tracers_ocn+1:) = bgc%l_isodiff_tracers
    endif

    ! index of transported tracers
    n_tracers_trans = 0
    do l=1,n_tracers_tot
      if (ocn%l_tracers_trans(l)) then
        n_tracers_trans = n_tracers_trans+1
        idx_tracers_trans(n_tracers_trans) = l
      endif
    enddo
    print *,'n_tracers_trans',n_tracers_trans

    ! Initialise sea ice model
    if (flag_sic) then
       if (allocated(ocn%grid%dz)) then
         call sic_init(sic,geo%f_ocn2,ocn%grid%dz(ocn%grid%nk))
       else
         call sic_init(sic,geo%f_ocn2,10._wp)
       endif
       call sic_diag_init
    endif

    ! Initialise land model
    if (flag_lnd) then
       call lnd_init(lnd,geo%f_lnd,geo%f_ocn,geo%f_ice,geo%f_ice_grd,geo%f_lake)
       call lnd_diag_init
    endif

    ! Initialize smb
    if (flag_smb) then
      allocate(smb(n_ice_domain))
      do n=1,n_ice_domain
        call smb_init(smb_in,smb(n),n,ice_grid(n),cmn%grid,geo%hires%z_bed_1min,geo%hires%lon_1min,geo%hires%lat_1min)  
        call smb_diag_init(smb(n))
        where (smb(n)%grid_smb_to_cmn%ncells.gt.0)
          cmn%mask_smb = 1
        endwhere
      enddo
    endif

    ! Initialize bmb
    if (flag_bmb) then
      allocate(bmb(n_ice_domain))
      do n=1,n_ice_domain
        call bmb_init(bmb(n),ice_grid(n),cmn%grid) 
        call bmb_diag_init(bmb(n))
      enddo
    endif

    ! Initialise atmospheric CO2
    if (flag_co2) then
       call co2_init(co2)
       call co2_diag_init
    endif
    
    ! Initialise atmospheric CH4
    if (flag_ch4) then
       call ch4_init(ch4)
       call ch4_diag_init
    endif
    
    call cmn_diag_init

    ! transfer initial state to coupler
    call geo_to_cmn(geo,cmn)
    call bnd_to_cmn(bnd,cmn,1)
    if (flag_atm) call atm_to_cmn(atm,cmn)
    if (flag_co2) call co2_to_cmn(co2,cmn)
    if (flag_ch4) call ch4_to_cmn(ch4,cmn)
    if (flag_ocn) call ocn_to_cmn(ocn,cmn)
    if (flag_sic) call sic_to_cmn(sic,cmn) 
    if (flag_bgc) call bgc_to_cmn(bgc,cmn,ocn)
    if (flag_smb) then
      do n=1,n_ice_domain
        call smb_to_cmn(smb(n),cmn)
      enddo
    endif
    if (flag_lnd) call lnd_to_cmn(lnd,cmn)
    if (flag_bmb) then
      do n=1,n_ice_domain
        call bmb_to_cmn(bmb(n),cmn)
      enddo
    endif
    if (flag_ice) then
      do n=1,n_ice_domain
        call ice_to_cmn(ice(n),cmn)
      enddo
    endif

    print*
    print*,'*******************************************************'
    print*,' Initialisation complete, simulation starting'
    print*,'*******************************************************'


    do step=1,nstep

      call timer_update

      if (time_soy) print *
      if (l_debug_main_loop) print*, 'timer (year, mon, doy): ',year,mon,doy

      !$ if (time_soy) then
      !$   time_atm = 0._dp  
      !$   time_ocn = 0._dp
      !$   time_sic = 0._dp
      !$   time_lnd = 0._dp
      !$   time_bgc = 0._dp
      !$   time_bnd = 0._dp
      !$   time_geo = 0._dp
      !$   time_ice = 0._dp
      !$   time_smb = 0._dp
      !$   time_bmb = 0._dp
      !$ endif

      ! for carbon cycle spinup
      if (l_spinup_cc) then
        if (year.eq.1 .and. doy.eq.1) then
          ! set weathering scale to 1, will be re-computed later during spinup
          lnd%l0d%weath_scale = 1._wp
        endif
        if (time_spinup_cc_1) then
          ! continue with ocean and bgc only with prescribed daily saved input
          flag_atm = .false.
          flag_lnd = .false.
          flag_sic = .false.
        endif
        if (time_spinup_cc_2) then
          ! end bgc spinup phase 
          l_spinup_bgc = .false.
          ! end sediments spinup
          l_spinup_sed = .false.
        endif
      endif

      ! update boundary
      !$ time_ini = omp_get_wtime()
      call bnd_update(bnd)
      call bnd_to_cmn(bnd,cmn)
      !$ time_end = omp_get_wtime()
      !$ time_bnd = time_bnd + (time_end-time_ini)

      if (l_debug_main_loop) print*, 'bnd exchange done.'


      !---------------------------------------------------------
      ! update geo model
      if ((flag_geo .or. ifake_geo.eq.1) .and. time_call_geo) then 
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update geo model'
        call bnd_to_geo(bnd,geo)
        if (flag_ice) then
          call ice_to_geo(n_ice_domain,ice,geo)
        endif
        call cmn_to_geo(cmn,geo)
        call geo_update(geo)
        call geo_diag(geo)
        if (flag_ice) then
          call geo_to_ice(n_ice_domain,geo,ice)
        endif
        call geo_to_cmn(geo,cmn)
        !$ time_end = omp_get_wtime()
        !$ time_geo = time_geo + (time_end-time_ini)
      endif

      !---------------------------------------------------------
      ! update atmosphere model
      if (flag_atm .and. time_call_atm) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update atm model'
        if (l_aquaplanet) call aquaplanet(cmn)
        call cmn_to_atm(cmn,atm)
        call atm_update(atm)
        call atm_diag(atm)
        call atm_to_cmn(atm,cmn)
        !$ time_end = omp_get_wtime()
        !$ time_atm = time_atm + (time_end-time_ini)
        if (atm%error) exit
      endif

      !---------------------------------------------------------
      ! update surface mass balance model
      if (flag_smb .and. time_call_smb) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update smb model'
        call cmn_to_smb(cmn,smb_in)
        do n=1,n_ice_domain
          if (time_soy_smb) then
            if (flag_ice) then
              call ice_to_smb(ice(n),smb(n))
            else
              call geo_to_smb(geo, smb(n))
            endif
          endif
          call smb_update(smb_in,smb(n))
          call smb_diag(smb(n))
          call smb_to_cmn(smb(n),cmn)
          if (flag_ice .and. time_eoy_smb) then
            call smb_to_ice(smb(n),ice(n))
          endif
        enddo
        !$ time_end = omp_get_wtime()
        !$ time_smb = time_smb + (time_end-time_ini)
      endif

      !---------------------------------------------------------
      ! update ice melt to ocean model
      if (flag_bmb .and. time_call_bmb) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update bmb model'
        do n=1,n_ice_domain
          call cmn_to_bmb(cmn,bmb(n))
          if (time_soy_bmb) then
            if (flag_ice) then
              call ice_to_bmb(ice(n),bmb(n))
            else
              call geo_to_bmb(geo, bmb(n))
            endif
          endif
          call bmb_update(bmb(n))
          call bmb_diag(bmb(n))
          if (time_eoy_bmb) then
            call bmb_to_cmn(bmb(n),cmn)
            if (flag_ice) then
              call bmb_to_ice(bmb(n),ice(n))
            endif
          endif
        enddo
        !$ time_end = omp_get_wtime()
        !$ time_bmb = time_bmb + (time_end-time_ini)
      endif

      !---------------------------------------------------------
      ! update ice model
      if (flag_ice .and. time_call_ice) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update ice model'
        l_error = .false.
        do n=1,n_ice_domain
          call ice_update(ice(n),n,ice_model_name,n_year_ice,real(year_now,wp),time_out_ice)
          call ice_to_cmn(ice(n),cmn)
          if (ice(n)%error) l_error = l_error .or. ice(n)%error
        enddo
        if (l_error) exit
        !$ time_end = omp_get_wtime()
        !$ time_ice = time_ice + (time_end-time_ini)
      endif

      !---------------------------------------------------------
      ! update land model
      if (flag_lnd .and. time_call_lnd) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update lnd model'
        call cmn_to_lnd(cmn,lnd)
        call lnd_update_wrapper(lnd)
        call runoff_merge(cmn)
        if (flag_lakes) then
          call lakes_update(cmn)
        endif
        call lnd_diag(lnd%l2d,lnd%l0d)
        call lnd_to_cmn(lnd,cmn)
        !$ time_end = omp_get_wtime()
        !$ time_lnd = time_lnd + (time_end-time_ini)
      endif

      !---------------------------------------------------------
      ! update sea ice model
      if (flag_sic .and. time_call_sic) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update sic model'
        call cmn_to_sic(cmn,sic)
        call sic_update(sic)
        call sic_diag(sic)
        call sic_to_cmn(sic,cmn)
        !$ time_end = omp_get_wtime()
        !$ time_sic = time_sic + (time_end-time_ini)
        if (sic%error) exit
      endif

      !---------------------------------------------------------
      ! update ocean model
      if (flag_ocn .and. time_call_ocn) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update ocn model'
        call runoff_to_ocn(cmn)
        !print *,'before'
        !print *,sum(ocn%ts(11,1:ocn%grid%ni,1:ocn%grid%nj,1:ocn%grid%nk)*ocn%grid%ocn_vol)
        call cmn_to_ocn(cmn,ocn)
        call ocn_update(ocn)
        call ocn_diag(ocn)
        call ocn_to_cmn(ocn,cmn)
        !print *,'after'
        !print *,sum(ocn%ts(11,1:ocn%grid%ni,1:ocn%grid%nj,1:ocn%grid%nk)*ocn%grid%ocn_vol)
        !print *
        !$ time_end = omp_get_wtime()
        !$ time_ocn = time_ocn + (time_end-time_ini)
        if (ocn%error) exit
      endif

      !---------------------------------------------------------
      ! update ocean biogeochemistry
      if (flag_bgc .and. time_call_bgc) then
        !$ time_ini = omp_get_wtime()
        if (l_debug_main_loop) print*, 'update bgc model'
        call cmn_to_bgc(cmn,ocn, bgc)
        call bgc_update(bgc)
        call bgc_diag(bgc)
        call bgc_to_cmn(bgc, cmn,ocn)
        !$ time_end = omp_get_wtime()
        !$ time_bgc = time_bgc + (time_end-time_ini)
      endif

      !---------------------------------------------------------
      ! update atmospheric CH4
      if (flag_ch4 .and. time_eoy) then
        call cmn_to_ch4(cmn,ch4)
        call ch4_update(ch4,real(year_now,kind=wp))
        call ch4_diag(ch4)
        call ch4_to_cmn(ch4,cmn)
      endif

      !---------------------------------------------------------
      ! update atmospheric CO2 
      if (flag_co2 .and. time_eoy) then
        call cmn_to_co2(cmn,co2)
        call co2_update(co2,real(year_now,kind=wp))
        call co2_diag(co2)
        call co2_to_cmn(co2,cmn)
      endif

      if (time_call_clim) then
        call cmn_diag(cmn)
      endif

      !$ if (time_eoy .and. l_write_timer) then
      !$   print *,'time atm',time_atm   
      !$   print *,'time ocn',time_ocn 
      !$   print *,'time sic',time_sic 
      !$   print *,'time lnd',time_lnd 
      !$   print *,'time bgc',time_bgc 
      !$   print *,'time bnd',time_bnd
      !$   print *,'time geo',time_geo 
      !$   print *,'time ice',time_ice 
      !$   print *,'time smb',time_smb 
      !$   print *,'time bmb',time_bmb 
      !$ endif

      !---------------------------------------------------------
      ! write restart files if required 
      if (time_write_restart) then
        call write_restart(trim(out_dir)//"/restart_out/", year_now)
      endif

    enddo


    ! check if loop was aborted because errors were encountered
    l_error = .false.
    if (flag_atm) l_error = l_error .or. atm%error
    if (flag_ocn) l_error = l_error .or. ocn%error
    if (flag_sic) l_error = l_error .or. sic%error
    if (flag_ice) then
      do n=1,n_ice_domain
        l_error = l_error .or. ice(n)%error
      enddo
    endif

    if (.not.l_error) then

      print*
      print*,'*******************************************************'
      print*,' Simulation complete, shutdown starting'
      print*,'*******************************************************'
      print*

      ! End geo
      call geo_end(geo)

      ! End atmosphere model
      if (flag_atm) call atm_end(atm)

      ! End land model
      if (flag_lnd) call lnd_end(lnd%l2d)

      ! End ocean biogeochemistry model
      if (flag_bgc) call bgc_end(bgc)

      ! End ocean model
      if (flag_ocn) call ocn_end(ocn)

      ! End sea ice model
      if (flag_sic) call sic_end(sic)

      ! End smb
      if (flag_smb) then
        do n=1,n_ice_domain
          if (flag_smb) call smb_end(smb_in, smb(n))
        enddo
      endif

      ! End bmb
      !    if (flag_bmb) call bmb_end(bmb)

      if (l_aquaplanet) call aqua_end(cmn)

      print*
      print*,'*******************************************************'
      print*,' Shutdown complete; home time'
      print*,'*******************************************************'
      print*

    else

      print *,'run aborted because of instabilities:'
      print *,'atm error : ',atm%error
      print *,'ocn error : ',ocn%error
      print *,'sic error : ',sic%error
      if (flag_ice) then
        do n=1,n_ice_domain
          print *,'ice n, error : ',n,ice(n)%error
        enddo
      endif

      call write_restart(trim(out_dir)//"/restart_crash/", 9999)

      stop 'ERROR: run aborted because of instabilities, restart files written to restart_crash directory'

    endif


contains

subroutine restart_init(folder)

  implicit none

  character(len=*), intent(in) :: folder

  integer :: cstat, estat
  character(len=10) :: year_now_str
  character(len=256) :: rest_dir, cmsg

  print *, '== restart_init ==============================='
  rest_dir = trim(folder)//"/restart_crash//year_9999"
  call execute_command_line('mkdir -p ' // adjustl(trim(rest_dir)), &
  exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)

  if (cstat == 0) then
    print *, 'mkdir restart_crash/year_9999 with status', estat
  else
    print *, 'command execution failed with error ', trim(cmsg)
  end if
  do step=1,nstep

    call timer_update

    if (time_write_restart) then
      ! cast year_now to string
      write (year_now_str, '(I10)')  year_now

      rest_dir = adjustl(trim(folder))//"/restart_out/"//"/year_"//adjustl(trim(year_now_str))

      ! make directory 
      call execute_command_line('mkdir -p ' // adjustl(trim(rest_dir)), &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)

      if (cstat == 0) then
        print *, 'mkdir ',adjustl(trim(rest_dir)),' with status', estat
      else
        print *, 'command execution failed with error ', trim(cmsg)
      end if

      rest_dir = adjustl(trim(folder))//"/restart_out/"//"/year_"//trim(adjustl(year_now_str))//"/vilma"

      ! make directory 
      call execute_command_line('mkdir -p ' // adjustl(trim(rest_dir)), &
      exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)

      if (cstat == 0) then
        print *, 'mkdir ',adjustl(trim(rest_dir)),' with status', estat
      else
        print *, 'command execution failed with error ', trim(cmsg)
      end if
    end if

  end do

  ! reset timer
  doy = 0
  mon = 0
  year = 0
  step = 1
  call timer_update

  return

end subroutine restart_init


subroutine write_restart(restart_out_dir, year_now)

  implicit none

  character(len=*), intent(in) :: restart_out_dir
  integer, intent(in) :: year_now

  character(len=10) :: year_now_str
  character(len=256) :: rest_dir


  ! cast year_now to string
  write (year_now_str, '(I10)')  year_now

  rest_dir = adjustl(trim(restart_out_dir))//"/year_"//adjustl(trim(year_now_str))

  ! write restart files
  call geo_write_restart(trim(rest_dir),trim(rest_dir)//"/geo_restart.nc",geo)
  if (flag_atm) call atm_write_restart(trim(rest_dir)//"/atm_restart.nc",atm)
  if (flag_ocn) call ocn_write_restart(trim(rest_dir)//"/ocn_restart.nc",ocn)
  if (flag_bgc) call bgc_write_restart(trim(rest_dir)//"/bgc_restart.nc",bgc)
  if (flag_sic) call sic_write_restart(trim(rest_dir)//"/sic_restart.nc",sic)
  if (flag_lnd) call lnd_write_restart(trim(rest_dir)//"/lnd_restart.nc",lnd%l2d,lnd%l0d)
  if (flag_smb) then
    do n=1,n_ice_domain
      call smb_write_restart(trim(rest_dir)//"/smb_"//trim(smb(n)%grid%name)//"_restart.nc",smb(n))
    enddo
  endif
  if (flag_bmb) then
    do n=1,n_ice_domain
      call bmb_write_restart(trim(rest_dir)//"/bmb_"//trim(bmb(n)%grid%name)//"_restart.nc",bmb(n))
    enddo
  endif
  if (flag_ice) then
    do n=1,n_ice_domain
      call ice_write_restart(ice_model_name,ice(n),n,real(year_now,wp),trim(rest_dir))
    enddo
  endif
  if (flag_co2) call co2_write_restart(trim(rest_dir)//"/co2_restart.nc",co2)
  if (flag_ch4) call ch4_write_restart(trim(rest_dir)//"/ch4_restart.nc",ch4)

  print *,'restart files written at year: ', year_now, ' in ', trim(rest_dir)

  return

end subroutine write_restart

end program climber


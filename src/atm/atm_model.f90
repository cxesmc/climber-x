!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : a t m _ m o d e l
!
!  Purpose : main of atmospheric model SESAM
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
module atm_model

    use atm_params, only : wp
    use ncio
    use constants, only : T0, fqsat, pi, frac_vu
    use timer, only: nday_year, doy, time_soy_atm, time_eoy_atm
    use timer, only : time_feedback_save, time_feedback_analysis
    use timer, only : monthly2daily
    use climber_grid, only: lon, lat 
    use control, only: restart_in_dir, atm_restart, flag_dust
    use control, only : l_feedbacks

    use atm_params, only : atm_params_init
    use atm_params, only : nstep_fast
    use atm_params, only : p0, ra, nsmooth_cda
    use atm_params, only : ars_ot, ars_im, l_dust, l_dust_rad
    use atm_params, only : i_rbstr
    use atm_params, only : l_write_timer
    use atm_params, only : tam_init
    use atm_params, only : ecs_scale, ecs_scale_dT
    use atm_grid, only : atm_grid_init, atm_grid_update
    use atm_grid, only : im, imc, jm, jmc, km, kmc, k700, nm, cost, pl, zl
    use atm_grid, only : i_ocn, i_sic, i_lnd, i_ice, i_lake
    use atm_def, only : atm_class

    use lw_radiation_mod, only : lw_radiation
    use sw_radiation_mod, only : sw_radiation
    use wvel_mod, only : wvel
    use clouds_mod, only : clouds
    use vesta_mod, only : hscales, vesta, tropoheight
    use crisa_mod, only : crisa
    use slp_mod, only : azslp, zslp
    use u2d_mod, only : u2d, usur
    use u3d_mod, only : u3d
    use synop_mod, only : synop
    use adifa_mod, only : adifa
    use dust_mod, only : dust
    use time_step_mod, only : time_step    
    use feedbacks_mod, only : feedback_type, feedback_init, feedback_save, feedback_analysis, feedback_write
    use rad_kernels_mod, only : rad_kernels_type, rad_kernels_init, rad_kernels, rad_kernels_write
    use smooth_atm_mod, only : smooth2
    !$  use omp_lib

    implicit none
    
    type(feedback_type) :: fb 
    type(rad_kernels_type) :: rk

    real(wp) :: eke_mon(im,jm,0:13)
    integer, dimension(nday_year) :: m0, m1
    real(wp), dimension(nday_year) :: wtm0, wtm1

    private
    public :: atm_init, atm_update, atm_end
    public :: atm_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ u p d a t e
  ! Purpose  :  update atmosphere
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_update(atm)

    implicit none

    type(atm_class) :: atm

    integer :: niter
    real(wp) :: ecs_scale_now
    !$ real(wp) :: time1,time2

    real(wp), dimension(im,jm) :: prc_save

    
    atm%eke = wtm0(doy)*eke_mon(:,:,m0(doy)) + wtm1(doy)*eke_mon(:,:,m1(doy)) 

    !$ if(l_write_timer) print *

    !-------------------------------------------------
    ! update topography
    !-------------------------------------------------
    if (time_soy_atm) then
      call atm_grid_update(atm%zs, atm%frst, atm%zsa, &      ! in
        atm%zsa_smooth, atm%slope, atm%slope_x, atm%slope_y, atm%ra2, atm%ra2a, atm%pzsa0, atm%pzsa, atm%psa, atm%ps)   ! out
    endif

    !-------------------------------------------------
    ! cross isobar angle
    !-------------------------------------------------
    if (doy.eq.1 .or. mod(doy,10).eq.0) then
      !$ time1 = omp_get_wtime()
      call crisa(atm%frst, atm%z0m, atm%zoro, &
        atm%cd, atm%cda, atm%cd0, atm%cd0a, atm%acbar, atm%sin_cos_acbar, atm%cos_acbar, atm%sin_acbar, atm%epsa)
      call smooth2(atm%cda,nsmooth_cda)
      !$ time2 = omp_get_wtime()
      !$ if(l_write_timer) print *,'crisa',time2-time1
    endif

    !-------------------------------------------------
    ! longwave radiation
    !-------------------------------------------------
    !$ time1 = omp_get_wtime()
    ecs_scale_now = ecs_scale + ecs_scale_dT*atm%dt2m_glob_ann_cum
    call lw_radiation(ecs_scale_now, atm%frst, atm%zsa, atm%zs, atm%htrop, atm%hcld, atm%ra2, &   ! in
      atm%gams, atm%gamb, atm%gamt, atm%tam, atm%ram, atm%hrm, atm%ttrop, atm%cld, atm%clot, &
      atm%co2, atm%ch4, atm%n2o, atm%cfc11, atm%cfc12, atm%co2e, atm%o3, &  ! in
      atm%flwr_up_sur, &  ! in
      atm%lwr_sur, atm%flwr_dw_sur, atm%flwr_dw_sur_cs, atm%flwr_dw_sur_cld, &    ! out
      atm%lwr_top, atm%lwr_top_cs, atm%lwr_top_cld, atm%lwr_tro, atm%lwr_cld)    ! out
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'lwr',time2-time1

    if (flag_dust .and. l_dust .and. l_dust_rad) then
      atm%aerosol_ot(:,:) = atm%dust_ot(:,:)
    else
      atm%aerosol_ot(:,:) = ars_ot 
    endif
    atm%aerosol_im(:,:) = ars_im

    !-------------------------------------------------
    ! shortwave radiation
    !-------------------------------------------------
    !$ time1 = omp_get_wtime()
    call sw_radiation(atm%frst, atm%swr_dw_top, atm%coszm, atm%cld, atm%q2, atm%ra2, &    ! in
      atm%alb_vu_s, atm%alb_vu_c, atm%alb_ir_s, atm%alb_ir_c, atm%clot, atm%hcld, atm%hqeff, atm%aerosol_ot, atm%aerosol_im, atm%so4, &    ! in
      atm%alb_cld, atm%swr_top, atm%swr_top_cs, atm%swr_top_cld, atm%swr_sur, atm%fswr_sur, atm%fswr_sur_cs, atm%fswr_sur_cld, & ! out
      atm%dswd_dalb_vu_cs, atm%dswd_dalb_ir_cs, atm%dswd_dalb_vu_cld, atm%dswd_dalb_ir_cld, atm%dswd_dz_ir_cs, atm%dswd_dz_ir_cld, &  ! out
      atm%swr_dw_sur_vis_cs, atm%swr_dw_sur_nir_cs, atm%swr_dw_sur_vis_cld, atm%swr_dw_sur_nir_cld)
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'swr',time2-time1

    !-------------------------------------------------
    ! total radiation
    !-------------------------------------------------
    atm%rb_top = atm%swr_top + atm%lwr_top
    atm%rb_sur = atm%swr_sur + atm%lwr_sur
    atm%rb_atm = atm%rb_top  - atm%rb_sur 
    if (i_rbstr.eq.1) then
      ! longwave radiation balance only
      atm%rb_str = atm%lwr_top - atm%lwr_tro
    else if (i_rbstr.eq.2) then
      ! include also solar radiation absorbtion by stratospheric O3
      atm%rb_str = atm%lwr_top - atm%lwr_tro + frac_vu*0.02*atm%swr_dw_top    ! 0.02=1-ITF_O3
    endif

    !-------------------------------------------------
    ! update prognostic tropopause height
    !-------------------------------------------------
    !$ time1 = omp_get_wtime()
    call tropoheight(atm%had_fi, atm%had_width, atm%rb_str, atm%hcld, &      ! in
      atm%htrop, atm%ptrop)      ! out
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'htrop',time2-time1

    !-------------------------------------------------
    ! feedback analysis
    !-------------------------------------------------
    if (l_feedbacks .and. time_feedback_analysis) then 
      ! feedback analysis
      call feedback_analysis(fb, atm%frst, atm%zs, atm%zsa, atm%htrop, atm%hcld, atm%tskin, atm%t2, atm%ra2, & 
        atm%gams, atm%gamb, atm%gamt, atm%tam, atm%ram, atm%hrm, atm%hqeff, atm%q2, atm%ttrop, atm%cld, &
        atm%co2, atm%ch4, atm%n2o, atm%cfc11, atm%cfc12, atm%o3, atm%flwr_up_sur, &  
        atm%swr_dw_top, atm%coszm, atm%alb_vu_s, atm%alb_vu_c, atm%alb_ir_s, atm%alb_ir_c, atm%clot, atm%aerosol_ot, atm%aerosol_im, atm%so4, &
        atm%had_fi, atm%had_width)
      ! radiative kernels 
      call rad_kernels(rk, atm%frst, atm%zs, atm%zsa, atm%htrop, atm%hcld, atm%ra2, & 
        atm%gams, atm%gamb, atm%gamt, atm%tam, atm%ram, atm%hrm, atm%hqeff, atm%q2, atm%ttrop, atm%cld, &
        atm%co2, atm%ch4, atm%n2o, atm%cfc11, atm%cfc12, atm%o3, atm%flwr_up_sur, &  
        atm%swr_dw_top, atm%coszm, atm%alb_vu_s, atm%alb_vu_c, atm%alb_ir_s, atm%alb_ir_c, atm%clot, atm%aerosol_ot, atm%aerosol_im, atm%so4)
      if (time_eoy_atm) then
        call feedback_write(fb)
        call rad_kernels_write(rk)
      endif
    endif

    !-------------------------------------------------
    ! dynamics 
    !-------------------------------------------------

    !-------------------------------------------------
    ! azonal sea level pressure
    !$ time1 = omp_get_wtime()
    call azslp(atm%frst, atm%tsksl, atm%htrop, atm%zsa, atm%uz500, &  ! in
      atm%aslp, &   ! inout
      atm%aslp_temp, atm%aslp_topo, atm%dz500, atm%atsl) ! out
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'azslp',(time2-time1)

    !-------------------------------------------------
    ! zonal sea level pressure and total sea level pressure
    !$ time1 = omp_get_wtime()
    call zslp(atm%zsa, atm%sin_cos_acbar, atm%tsl, atm%aslp, & ! in
      atm%slp, atm%had_fi, atm%had_width)  ! out
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'zslp',(time2-time1)

    !-------------------------------------------------
    ! 2d geostrophic and ageostrophic wind components in the PBL
    !$ time1 = omp_get_wtime()
    call u2d(atm%slp, atm%sin_cos_acbar, &  ! in
      atm%ugb, atm%vgb, atm%ugbf, atm%vgbf, atm%uab, atm%vab) ! out
    !$ time2 = omp_get_wtime()
    !$ if(l_write_timer) print *,'u2d',(time2-time1)


    ! fast time step

    ! initialize precipitation to be cumulated over fast time steps
    prc_save = atm%prc
    atm%prc  = 0._wp 
    atm%prcw = 0._wp 
    atm%prcs = 0._wp 

    do niter=1,nstep_fast

      !-------------------------------------------------
      ! 3d velocity field
      !-------------------------------------------------
      !$ time1 = omp_get_wtime()
      call u3d(niter, atm%pzsa, atm%ptrop, atm%ugb, atm%vgb, atm%uab, atm%vab, atm%t3, & ! in
        atm%ua, atm%va, atm%uter, atm%vter, atm%uterf, atm%vterf, atm%u3, atm%v3, atm%w3, atm%uz500, &  ! out   
        atm%fax, atm%faxo, atm%fay, atm%fayo, atm%fac)  ! out
      !$ time2 = omp_get_wtime()
      !$ if(l_write_timer .and. niter.eq.1) print *,'u3d1',(time2-time1)*nstep_fast
      !$ if(l_write_timer .and. niter.eq.2) print *,'u3d2',(time2-time1)*nstep_fast

      !-------------------------------------------------
      ! 2d near surface (10 m) wind components
      !-------------------------------------------------
      if (niter.eq.1) then
        !$ time1 = omp_get_wtime()
        call usur(atm%ugbf, atm%vgbf, atm%epsa, atm%cos_acbar, atm%sin_acbar, atm%t2a, atm%tskina, atm%cd0a, atm%slope_x, atm%slope_y, &
          atm%usk, atm%vsk, & 
          atm%us, atm%vs)
        !$ time2 = omp_get_wtime()
        !$ if(l_write_timer) print *,'usur',(time2-time1)
      endif

      !-------------------------------------------------
      ! effective vertical velocities
      !-------------------------------------------------
      if (niter.eq.1) then
        !$ time1 = omp_get_wtime()
        call wvel(atm%w3, atm%wsyn, atm%winda, atm%sigoro, &    ! in
          atm%wcld, atm%woro, atm%weff)   ! out
        !$ time2 = omp_get_wtime()
        !$ if(l_write_timer) print *,'wvel',time2-time1
      endif

      !-------------------------------------------------
      ! clouds
      !-------------------------------------------------
      if (niter.eq.1) then
        !$ time1 = omp_get_wtime()
        call clouds(atm%frst, atm%weff, atm%wcld, atm%zsa, atm%t2a, atm%ram, atm%qam, atm%rskina, atm%wcon, atm%htrop, atm%so4, atm%sam2, &    ! in
          atm%fweff, atm%cld_rh, atm%cld_low, atm%cld, atm%hcld, atm%clot) ! out
        !$ time2 = omp_get_wtime()
        !$ if(l_write_timer) print *,'cld',time2-time1
      endif


      !-------------------------------------------------
      ! vertical structure of atmosphere
      !-------------------------------------------------

      ! lapse rate and height scales of moisture and dust
      !$ time1 = omp_get_wtime()
      call hscales(atm%frst, atm%f_ice_lake, atm%ra2a, atm%rb_sur, atm%tam, atm%tskin, atm%qam, atm%wcon, atm%wcld, &  ! in
        atm%had_fi, atm%had_width, &   ! in
        atm%gams, atm%gamb, atm%gamt, atm%hrm, &    ! inout
        atm%hqeff, atm%hdust)    ! out
      !$ time2 = omp_get_wtime()
      !$ if(l_write_timer .and. niter.eq.1) print *,'hscales',(time2-time1)*nstep_fast

      ! vertical profiles of temperature, humidity and dust
      !$ time1 = omp_get_wtime()
      call vesta(atm%zsa, atm%tam, atm%gams, atm%gamb, atm%gamt, atm%htrop, atm%ram, atm%hrm, atm%dam, atm%hdust, &  ! in
        atm%wcon, atm%t3, atm%q3, atm%tp, atm%d3, atm%ttrop)    ! out
      !$ time2 = omp_get_wtime()
      !$ if(l_write_timer .and. niter.eq.1) print *,'vesta',(time2-time1)*nstep_fast

      !-------------------------------------------------
      ! dust
      !-------------------------------------------------
      if (flag_dust .and. l_dust .and. niter.eq.1) then
        ! compute dust load, deposition and optical thickness
        !$ time1 = omp_get_wtime()
        call dust(atm%dam, prc_save, atm%ra2a, atm%hdust, &
          atm%dust_load, atm%dust_dep_dry, atm%dust_dep_wet, atm%dust_dep, atm%dust_ot)
        !$ time2 = omp_get_wtime()
        !$ if(l_write_timer .and. niter.eq.1) print *,'dust',time2-time1
      endif

      !-------------------------------------------------
      ! synoptic processes
      !-------------------------------------------------
      if (niter.eq.1) then
        !$ time1 = omp_get_wtime()
        call synop(atm%frst, atm%zs, atm%uter, atm%vter, atm%uterf, atm%vterf, atm%u3(:,:,k700), atm%v3(:,:,k700), atm%us, atm%vs, atm%tp, &    ! in
          atm%zsa, atm%sigoro, atm%cda, atm%cd, atm%epsa, atm%cos_acbar, &    ! in
          atm%sam, atm%sam2, atm%cdif, &    ! inout
          atm%synprod, atm%syndiss, atm%synadv, atm%syndif, atm%synsur, atm%winda, atm%wind, atm%taux, atm%tauy, &  ! out 
          atm%diffxdse, atm%diffydse, atm%diffxwtr, atm%diffywtr, atm%diffxdst, atm%diffydst, atm%wsyn)    ! out
        !$ time2 = omp_get_wtime()
        !$ if(l_write_timer) print *,'synop',time2-time1
      endif

      !-------------------------------------------------
      ! advection-diffusion
      !-------------------------------------------------
      !$ time1 = omp_get_wtime()
      call adifa(atm%fax, atm%fay, atm%tp, atm%q3, atm%d3, atm%cam, &
        atm%diffxdse, atm%diffydse, atm%diffxwtr, atm%diffywtr, atm%diffxdst, atm%diffydst,  &   ! in
        atm%convdse, atm%convwtr, atm%convdst, atm%convco2, &   ! out
        atm%faxdse, atm%faxwtr, atm%faxdst, atm%faxco2, &  ! out
        atm%faydse, atm%faywtr, atm%faydst, atm%fayco2, &  ! out
        atm%fdxdse, atm%fdxwtr, atm%fdxdst, atm%fdxco2, &  ! out
        atm%fdydse, atm%fdywtr, atm%fdydst, atm%fdyco2)   ! out
      !$ time2 = omp_get_wtime()
      !$ if(l_write_timer .and. niter.eq.1) print *,'adifa',(time2-time1)*nstep_fast

      !-------------------------------------------------
      ! time step, prognostic equations for temperature, humidity and dust
      !-------------------------------------------------
      !$ time1 = omp_get_wtime()
      call time_step(atm%frst, atm%zs, atm%zsa, atm%ps, atm%psa, atm%ra2a, atm%slope, atm%evpa, atm%convwtr, atm%wcon, atm%hqeff, atm%sam, atm%eke, atm%sam2, &   ! in
        atm%tskin, atm%convdse, atm%rb_atm, atm%rb_sur, atm%sha, atm%gams, atm%gamb, atm%gamt, &     ! in
        atm%convdst, atm%dust_emis, atm%dust_dep, atm%hdust, &     ! in
        atm%convco2, atm%co2flx, &     ! in
        atm%tam, atm%qam, atm%dam, atm%cam, atm%prc, atm%prcw, atm%prcs, atm%prc_conv, atm%prc_wcon, atm%prc_over, &   ! inout
        atm%q2, atm%q2a, atm%ram, atm%r2, atm%r2a, atm%rskina, atm%tsl, atm%tsksl, atm%t2, atm%t2a, atm%tskina, atm%error)  ! out
      !$ time2 = omp_get_wtime()
      !$ if(l_write_timer .and. niter.eq.1) print *,'time_Step',(time2-time1)*nstep_fast

      if (atm%error) exit

    enddo

    ! mean precipitation over fast time steps
    atm%prc  = atm%prc/nstep_fast
    atm%prcw = atm%prcw/nstep_fast
    atm%prcs = atm%prcs/nstep_fast

    !-------------------------------------------------
    ! save for feedback analysis
    !-------------------------------------------------
    if (l_feedbacks .and. time_feedback_save) then
      ! save fields needed for feedback analysis in derived type fb
      call feedback_save(atm%co2, atm%tam, atm%cld, atm%hcld, atm%clot, atm%gams, atm%gamb, atm%gamt, &
        atm%htrop, atm%ttrop, atm%ram, atm%hrm, atm%hqeff, atm%q2, atm%aerosol_ot, atm%aerosol_im, atm%so4, &
        atm%frst, atm%tskin, atm%t2, atm%alb_vu_s, atm%alb_vu_c, atm%alb_ir_s, atm%alb_ir_c, atm%flwr_up_sur, &
        fb)
    endif

   return

  end subroutine atm_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ i n i t
  ! Purpose  :  initialize atm
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_init(atm,frlnd,zsa,sigoro,co2)

    implicit none

    type(atm_class) :: atm
    real(wp), intent(in) :: frlnd(:,:)
    real(wp), intent(in) :: zsa(:,:)
    real(wp), intent(in) :: sigoro(:,:)
    real(wp), intent(in) :: co2

    integer :: i, j, ja, n


    ! initialize/read atmosphere model parameters
    call atm_params_init

    ! initialize atmospheric grid
    call atm_grid_init

    ! allocate 
    call atm_alloc(atm)
    
    atm%pl = pl(:)
    atm%zl = zl(:)

    ! CO2 mass mixing ratio
    atm%cam(:,:) = co2*1.e-6_wp * 44.0095_wp/28.97_wp  ! kgCO2/kg

    if (atm_restart) then
      ! read restart file

      call atm_read_restart(trim(restart_in_dir)//"/atm_restart.nc",atm)

    else
      ! initial atmosperic state without restart      
      
      atm%had_fi   = 0._wp
      atm%had_width = pi/3._wp

       do i=1,im
        do j=1,jm
         ja = jm-j+1

         atm%zsa(i,j) = max(0._wp,zsa(i,ja))
         atm%sigoro(i,j) = sigoro(i,ja)
         atm%frlnd(i,j) = frlnd(i,ja)
         atm%f_ice_lake(i,j) = 0._wp
         atm%frst(i,j,:) = 0._wp   ! fixme
         atm%ra2a(i,j) = ra

         !atm%tam(i,j) = 30._wp*cost(j)+T0
         atm%tam(i,j) = tam_init*cost(j)+T0
         atm%t2a(i,j) = atm%tam(i,j)
         atm%tskina(i,j) = atm%tam(i,j)
         atm%tsl(i,j) = atm%tam(i,j)
         atm%tsksl(i,j) = atm%tam(i,j)
         atm%gams(i,j) = 5.e-3_wp
         atm%gamb(i,j) = 5.e-3_wp
         atm%gamt(i,j) = 8.e-3_wp
         atm%hrm(i,j) = 2000._wp
         atm%qam(i,j) = 0.8_wp*fqsat(atm%tam(i,j),p0)
         atm%q2a(i,j) = atm%qam(i,j)
         atm%ram(i,j) = 0.8_wp
         atm%hqeff(i,j) = 2000._wp
         atm%wcon(i,j) = 1._wp
         atm%prc(i,j) = 0._wp
         atm%prcs(i,j,:) = 0._wp
         atm%prcw(i,j,:) = 0._wp
         atm%evpa(i,j) = 0._wp
         atm%cld(i,j) = 0.5_wp
         atm%cld_rh(i,j) = 0.5_wp
         atm%cld_low(i,j) = 0._wp
         atm%hcld(i,j) = 4000._wp
         atm%clot(i,j) = 1._wp
         atm%htrop(i,j) = 12.e3_wp*(1._wp+0.5_wp*cost(j))
         atm%rb_sur(i,j) = 0._wp
         atm%rb_atm(i,j) = 0._wp
         atm%sam(i,j) = 0._wp
         atm%sam2(i,j) = 0._wp
         atm%cdif(i,j) = 0._wp
         atm%ps(i,j,:) = p0
         atm%psa(i,j) = p0
         atm%slp(i,j) = p0
         atm%aslp(i,j) = 0._wp      
         atm%aslp_topo(i,j) = 0._wp      
         atm%dz500(i,j) = 0._wp      

         atm%winda(i,j) = 1._wp
         atm%wind(i,j,:) = 1._wp
         atm%us(i,j,:) = 0._wp
         atm%vs(i,j,:) = 0._wp
         atm%usk(i,j) = 0._wp
         atm%vsk(i,j) = 0._wp
         atm%ugb(i,j) = 0._wp
         atm%vgb(i,j) = 0._wp
         atm%uab(i,j) = 0._wp
         atm%vab(i,j) = 0._wp
         atm%wsyn(i,j) = 0._wp
         atm%wcld(i,j) = 0._wp
         atm%weff(i,j) = 0._wp
         atm%woro(i,j) = 0._wp
         atm%fweff(i,j) = 0._wp
         atm%uz500(j) = 0._wp

         atm%convwtr(i,j) = 0._wp
         atm%convdse(i,j) = 0._wp

         atm%dam(i,j) = 0._wp
         atm%hdust(i,j) = 2000._wp
         atm%dust_ot(i,j) = 0._wp
         atm%dust_load(i,j) = 0._wp 
         atm%dust_dep(i,j) = 0._wp   
         atm%dust_dep_dry(i,j) = 0._wp 
         atm%dust_dep_wet(i,j) = 0._wp 

         atm%so4(i,j) = 0._wp

         do n=1,nm
           atm%tskin(i,j,n) = atm%tam(i,j)
           atm%t2(i,j,n) = atm%tam(i,j)
           atm%q2(i,j,n) = atm%qam(i,j)
           atm%r2(i,j,n) = atm%ram(i,j)
         enddo

        enddo
       enddo        

         atm%u3 = 0._wp  
         atm%v3 = 0._wp
         atm%w3 = 0._wp

         atm%uter = 0._wp  
         atm%vter = 0._wp
         atm%uterf = 0._wp  
         atm%vterf = 0._wp

         call vesta(atm%zsa, atm%tam, atm%gams, atm%gamb, atm%gamt, atm%htrop, atm%ram, atm%hrm, atm%dam, atm%hdust, &  ! in
           atm%wcon, atm%t3, atm%q3, atm%tp, atm%d3, atm%ttrop)    ! out

    endif

    if (l_feedbacks) then
      call feedback_init(fb)
      call rad_kernels_init(rk)
    endif


    call monthly2daily(m0,m1,wtm0,wtm1)
    call nc_read("input/era_interim_monthly_eke_1981_2010_5x5.nc","eke_2_6_int",eke_mon(:,jm:1:-1,1:12) )
    eke_mon(:,:,0)   = eke_mon(:,:,12)
    eke_mon(:,:,13)  = eke_mon(:,:,1)

    atm%t2m_glob_ann = 0._wp
    atm%dt2m_glob_ann_cum = 0._wp

    atm%error = .false.

    print*
    print*,'======================================================='
    print*,' Initialisation of ATMOSPHERE complete'
    print*,'======================================================='
    print*

  return

  end subroutine atm_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ a l l o c
  ! Purpose  :  allocate all state variables 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_alloc(atm)

    implicit none 

    type(atm_class) :: atm

     allocate(atm%pl(kmc))
     allocate(atm%zl(kmc))

     allocate(atm%solar(nday_year,24,jm))
     allocate(atm%solarm(nday_year,jm))
     allocate(atm%cosz(nday_year,24,jm))
     allocate(atm%coszm(nday_year,jm))

     allocate(atm%cam(im,jm))  
     allocate(atm%co2flx(im,jm))
     allocate(atm%C13flx(im,jm))
     allocate(atm%C14flx(im,jm))

     allocate(atm%idivide_pac_atl(jm))
     allocate(atm%idivide_atl_indpac(jm))

     allocate(atm%zs(im,jm,nm))
     allocate(atm%zsa(im,jm))
     allocate(atm%zsa_smooth(im,jm))
     allocate(atm%slope(im,jm))
     allocate(atm%slope_x(im,jm))
     allocate(atm%slope_y(im,jm))
     allocate(atm%pzsa0(im,jm))
     allocate(atm%pzsa(im,jm))
     allocate(atm%frlnd(im,jm))
     allocate(atm%frocn(im,jm))
     allocate(atm%f_ice_lake(im,jm))
     allocate(atm%sigoro(im,jm))

     allocate(atm%tam(im,jm))  
     allocate(atm%qam(im,jm))  
     allocate(atm%ram(im,jm))  
     allocate(atm%gams(im,jm)) 
     allocate(atm%gamb(im,jm)) 
     allocate(atm%gamt(im,jm)) 
     allocate(atm%dam(im,jm))  
     allocate(atm%hrm(im,jm))  
     allocate(atm%hqeff(im,jm)) 
     allocate(atm%wcon(im,jm)) 
     allocate(atm%cld_rh(im,jm)) 
     allocate(atm%cld_low(im,jm)) 
     allocate(atm%cld(im,jm)) 
     allocate(atm%cld_dat(im,jm)) 
     allocate(atm%cld_day_dat(im,jm)) 
     allocate(atm%prc(im,jm)) 
     allocate(atm%prcw(im,jm,nm))
     allocate(atm%prcs(im,jm,nm))
     allocate(atm%prc_conv(im,jm))
     allocate(atm%prc_wcon(im,jm))
     allocate(atm%prc_over(im,jm))
     allocate(atm%hcld(im,jm))
     allocate(atm%clot(im,jm)) 
     allocate(atm%alb_cld(im,jm)) 
     allocate(atm%htrop(im,jm)) 
     allocate(atm%ptrop(jm)) 
     allocate(atm%ttrop(im,jm))  
     allocate(atm%winda(im,jm))   
     allocate(atm%wind(im,jm,nm))   
     allocate(atm%aerosol_ot(im,jm)) 
     allocate(atm%aerosol_im(im,jm)) 
     allocate(atm%so4(im,jm))
     allocate(atm%o3(im,jm,kmc)) 

     allocate(atm%hdust(im,jm)) 
     allocate(atm%dust_load(im,jm)) 
     allocate(atm%dust_emis(im,jm)) 
     allocate(atm%dust_dep(im,jm)) 
     allocate(atm%dust_dep_dry(im,jm)) 
     allocate(atm%dust_dep_wet(im,jm)) 
     allocate(atm%dust_ot(im,jm)) 

     allocate(atm%frst(im,jm,nm))
     allocate(atm%tskin(im,jm,nm)) 
     allocate(atm%t2(im,jm,nm)) 
     allocate(atm%ra2(im,jm,nm))
     allocate(atm%q2(im,jm,nm))  
     allocate(atm%r2(im,jm,nm))   
     allocate(atm%alb_vu_s(im,jm,nm))
     allocate(atm%alb_vu_c(im,jm,nm))
     allocate(atm%alb_ir_s(im,jm,nm))
     allocate(atm%alb_ir_c(im,jm,nm))
     allocate(atm%cd(im,jm,nm)) 
     allocate(atm%cd0(im,jm,nm)) 
     allocate(atm%z0m(im,jm,nm)) 
     allocate(atm%zoro(im,jm)) 

     allocate(atm%ra2a(im,jm))
     allocate(atm%cda(im,jm))
     allocate(atm%cd0a(im,jm))
     allocate(atm%sha(im,jm))
     allocate(atm%lha(im,jm))
     allocate(atm%evpa(im,jm))
     allocate(atm%tskina(im,jm))
     allocate(atm%t2a(im,jm))
     allocate(atm%q2a(im,jm))
     allocate(atm%r2a(im,jm))
     allocate(atm%rskina(im,jm))
 
     allocate(atm%t3(im,jm,km))
     allocate(atm%q3(im,jm,km))
     allocate(atm%tp(im,jm,km)) 
     allocate(atm%d3(im,jm,km))
 
     allocate(atm%acbar(im,jm))
     allocate(atm%sin_cos_acbar(im,jm))
     allocate(atm%cos_acbar(im,jm,nm))
     allocate(atm%sin_acbar(im,jm,nm))
     allocate(atm%epsa(im,jm,nm))
     allocate(atm%slp(im,jm))
     allocate(atm%slp_dat(im,jm))
     allocate(atm%tsl_dat(im,jm))
     allocate(atm%ps(im,jm,nm))
     allocate(atm%psa(im,jm))
     allocate(atm%atsl(im,jm))
     allocate(atm%aslp(im,jm))
     allocate(atm%aslp_temp(im,jm))
     allocate(atm%aslp_topo(im,jm))
     allocate(atm%dz500(im,jm))
     allocate(atm%us(im,jm,nm))
     allocate(atm%vs(im,jm,nm))
     allocate(atm%usk(im,jm))
     allocate(atm%vsk(im,jm))
     allocate(atm%ugb(im,jm))
     allocate(atm%vgb(im,jm))
     allocate(atm%ugbf(im,jm))
     allocate(atm%vgbf(im,jm))
     allocate(atm%uab(imc,jm))
     allocate(atm%vab(im,jmc))
     allocate(atm%taux(im,jm,nm))
     allocate(atm%tauy(im,jm,nm))
     allocate(atm%uz500(jm))
     allocate(atm%wcld(im,jm))
     allocate(atm%woro(im,jm))
     allocate(atm%wsyn(im,jm))
     allocate(atm%weff(im,jm))   
     allocate(atm%fweff(im,jm))   

     allocate(atm%ua(im,jm,km))
     allocate(atm%va(im,jm,km))
     allocate(atm%u3(im,jm,km))
     allocate(atm%v3(im,jm,km))
     allocate(atm%uter(im,jm,km))
     allocate(atm%vter(im,jm,km))
     allocate(atm%uterf(im,jm,km))
     allocate(atm%vterf(im,jm,km))
     allocate(atm%fax(imc,jm,km))
     allocate(atm%faxo(imc,jm,km))
     allocate(atm%fay(im,jmc,km))
     allocate(atm%fayo(im,jmc,km))
     allocate(atm%fac(im,jm))
     allocate(atm%w3(im,jm,kmc))

     allocate(atm%convdse(im,jm))
     allocate(atm%convwtr(im,jm))
     allocate(atm%convdst(im,jm))
     allocate(atm%convco2(im,jm))
     allocate(atm%faxdse(imc,jm))
     allocate(atm%faxwtr(imc,jm))
     allocate(atm%faxdst(imc,jm))
     allocate(atm%faxco2(imc,jm))
     allocate(atm%faydse(im,jmc))
     allocate(atm%faywtr(im,jmc))
     allocate(atm%faydst(im,jmc))
     allocate(atm%fayco2(im,jmc))
     allocate(atm%fdxdse(imc,jm))
     allocate(atm%fdxwtr(imc,jm))
     allocate(atm%fdxdst(imc,jm))
     allocate(atm%fdxco2(imc,jm))
     allocate(atm%fdydse(im,jmc))
     allocate(atm%fdywtr(im,jmc))
     allocate(atm%fdydst(im,jmc))
     allocate(atm%fdyco2(im,jmc))

     allocate(atm%fswr_sur(im,jm,nm))
     allocate(atm%fswr_sur_cs(im,jm,nm))
     allocate(atm%fswr_sur_cld(im,jm,nm))
     allocate(atm%flwr_dw_sur(im,jm,nm))
     allocate(atm%flwr_dw_sur_cs(im,jm,nm))
     allocate(atm%flwr_dw_sur_cld(im,jm,nm))
     allocate(atm%flwr_up_sur(im,jm,nm))

     allocate(atm%dswd_dalb_vu_cs(im,jm))
     allocate(atm%dswd_dalb_ir_cs(im,jm))
     allocate(atm%dswd_dalb_vu_cld(im,jm))
     allocate(atm%dswd_dalb_ir_cld(im,jm))
     allocate(atm%dswd_dz_ir_cs(im,jm))
     allocate(atm%dswd_dz_ir_cld(im,jm))
     allocate(atm%swr_dw_sur_vis_cs(im,jm))
     allocate(atm%swr_dw_sur_nir_cs(im,jm))
     allocate(atm%swr_dw_sur_vis_cld(im,jm))
     allocate(atm%swr_dw_sur_nir_cld(im,jm))

     allocate(atm%rb_top(im,jm))
     allocate(atm%rb_sur(im,jm))
     allocate(atm%rb_atm(im,jm))
     allocate(atm%rb_str(im,jm))     
     allocate(atm%swr_dw_top(im,jm))
     allocate(atm%swr_top(im,jm))
     allocate(atm%swr_top_cs(im,jm))
     allocate(atm%swr_top_cld(im,jm))
     allocate(atm%swr_sur(im,jm))
     allocate(atm%lwr_top(im,jm))
     allocate(atm%lwr_top_cs(im,jm))
     allocate(atm%lwr_top_cld(im,jm))
     allocate(atm%lwr_sur(im,jm))
     allocate(atm%lwr_tro(im,jm))
     allocate(atm%lwr_cld(im,jm))

     allocate(atm%tsl(im,jm))
     allocate(atm%tsksl(im,jm))

     allocate(atm%eke(im,jm))
     allocate(atm%sam(im,jm))
     allocate(atm%sam2(im,jm))
     allocate(atm%synprod(im,jm))
     allocate(atm%syndiss(im,jm))
     allocate(atm%synadv(im,jm))
     allocate(atm%syndif(im,jm))
     allocate(atm%synsur(im,jm,nm))
     allocate(atm%cdif(im,jm) )    
     allocate(atm%diffxdse(imc,jm))
     allocate(atm%diffydse(im,jmc))
     allocate(atm%diffxwtr(imc,jm))
     allocate(atm%diffywtr(im,jmc))
     allocate(atm%diffxdst(imc,jm))
     allocate(atm%diffydst(im,jmc))

   return

  end subroutine atm_alloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ e n d
  ! Purpose  :  deallocate all state variables 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_end(atm)

    implicit none 

    type(atm_class) :: atm


    call atm_dealloc(atm)


   return

  end subroutine atm_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ d e a l l o c
  ! Purpose  :  deallocate all state variables 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_dealloc(atm)

    implicit none 

    type(atm_class) :: atm

     deallocate(atm%solar)
     deallocate(atm%solarm)
     deallocate(atm%cosz)
     deallocate(atm%coszm)

     deallocate(atm%cam)  
     deallocate(atm%co2flx)
     deallocate(atm%C13flx)
     deallocate(atm%C14flx)

     deallocate(atm%idivide_pac_atl)
     deallocate(atm%idivide_atl_indpac)

     deallocate(atm%zs)
     deallocate(atm%zsa)
     deallocate(atm%zsa_smooth)
     deallocate(atm%slope)
     deallocate(atm%slope_x)
     deallocate(atm%slope_y)
     deallocate(atm%pzsa0)
     deallocate(atm%pzsa)
     deallocate(atm%frlnd)
     deallocate(atm%frocn)
     deallocate(atm%f_ice_lake)
     deallocate(atm%sigoro)

     deallocate(atm%tam)  
     deallocate(atm%qam)  
     deallocate(atm%ram)  
     deallocate(atm%gams) 
     deallocate(atm%gamb) 
     deallocate(atm%gamt) 
     deallocate(atm%dam)  
     deallocate(atm%hrm)  
     deallocate(atm%hqeff) 
     deallocate(atm%wcon) 
     deallocate(atm%cld_rh) 
     deallocate(atm%cld_low) 
     deallocate(atm%cld) 
     deallocate(atm%cld_dat) 
     deallocate(atm%cld_day_dat) 
     deallocate(atm%prc) 
     deallocate(atm%prcw)
     deallocate(atm%prcs)
     deallocate(atm%prc_conv)
     deallocate(atm%prc_wcon)
     deallocate(atm%prc_over)
     deallocate(atm%hcld)
     deallocate(atm%clot) 
     deallocate(atm%alb_cld) 
     deallocate(atm%htrop) 
     deallocate(atm%ptrop) 
     deallocate(atm%ttrop)  
     deallocate(atm%winda)   
     deallocate(atm%wind)   
     deallocate(atm%aerosol_ot) 
     deallocate(atm%aerosol_im) 
     deallocate(atm%so4) 
     deallocate(atm%o3) 

     deallocate(atm%hdust) 
     deallocate(atm%dust_load) 
     deallocate(atm%dust_emis) 
     deallocate(atm%dust_dep) 
     deallocate(atm%dust_dep_dry) 
     deallocate(atm%dust_dep_wet) 
     deallocate(atm%dust_ot) 

     deallocate(atm%frst)
     deallocate(atm%t2) 
     deallocate(atm%tskin) 
     deallocate(atm%ra2)
     deallocate(atm%q2)  
     deallocate(atm%r2)   
     deallocate(atm%alb_vu_s)
     deallocate(atm%alb_vu_c)
     deallocate(atm%alb_ir_s)
     deallocate(atm%alb_ir_c)
     deallocate(atm%cd) 
     deallocate(atm%cd0)
     deallocate(atm%z0m) 
     deallocate(atm%zoro) 

     deallocate(atm%ra2a)
     deallocate(atm%cda)
     deallocate(atm%cd0a)
     deallocate(atm%sha)
     deallocate(atm%lha)
     deallocate(atm%evpa)
     deallocate(atm%tskina)
     deallocate(atm%t2a)
     deallocate(atm%q2a)
     deallocate(atm%r2a)
     deallocate(atm%rskina)
 
     deallocate(atm%t3)
     deallocate(atm%q3)
     deallocate(atm%tp) 
     deallocate(atm%d3)
 
     deallocate(atm%acbar)
     deallocate(atm%sin_cos_acbar)
     deallocate(atm%cos_acbar)
     deallocate(atm%sin_acbar)
     deallocate(atm%epsa)
     deallocate(atm%slp)
     deallocate(atm%slp_dat)
     deallocate(atm%tsl_dat)
     deallocate(atm%ps)
     deallocate(atm%psa)
     deallocate(atm%atsl)
     deallocate(atm%aslp)
     deallocate(atm%aslp_temp)
     deallocate(atm%aslp_topo)
     deallocate(atm%dz500)
     deallocate(atm%us)
     deallocate(atm%vs)
     deallocate(atm%usk)
     deallocate(atm%vsk)
     deallocate(atm%ugb)
     deallocate(atm%vgb)
     deallocate(atm%ugbf)
     deallocate(atm%vgbf)
     deallocate(atm%uab)
     deallocate(atm%vab)
     deallocate(atm%taux)
     deallocate(atm%tauy)
     deallocate(atm%uz500)
     deallocate(atm%wcld)
     deallocate(atm%woro)
     deallocate(atm%wsyn)
     deallocate(atm%weff)   
     deallocate(atm%fweff)   

     deallocate(atm%ua)
     deallocate(atm%va)
     deallocate(atm%u3)
     deallocate(atm%v3)
     deallocate(atm%w3)
     deallocate(atm%uter)
     deallocate(atm%vter)
     deallocate(atm%uterf)
     deallocate(atm%vterf)
     deallocate(atm%fax)
     deallocate(atm%faxo)
     deallocate(atm%fay)
     deallocate(atm%fayo)
     deallocate(atm%fac)
     deallocate(atm%diffxdse)
     deallocate(atm%diffydse)
     deallocate(atm%diffxwtr)
     deallocate(atm%diffywtr)
     deallocate(atm%diffxdst)
     deallocate(atm%diffydst)

     deallocate(atm%convdse)
     deallocate(atm%convwtr)
     deallocate(atm%convdst)
     deallocate(atm%convco2)
     deallocate(atm%faxdse)
     deallocate(atm%faxwtr)
     deallocate(atm%faxdst)
     deallocate(atm%faxco2)
     deallocate(atm%faydse)
     deallocate(atm%faywtr)
     deallocate(atm%faydst)
     deallocate(atm%fayco2)
     deallocate(atm%fdxdse)
     deallocate(atm%fdxwtr)
     deallocate(atm%fdxdst)
     deallocate(atm%fdxco2)
     deallocate(atm%fdydse)
     deallocate(atm%fdywtr)
     deallocate(atm%fdydst)
     deallocate(atm%fdyco2)

     deallocate(atm%fswr_sur)
     deallocate(atm%fswr_sur_cs)
     deallocate(atm%fswr_sur_cld)
     deallocate(atm%flwr_dw_sur)
     deallocate(atm%flwr_dw_sur_cs)
     deallocate(atm%flwr_dw_sur_cld)
     deallocate(atm%flwr_up_sur)

     deallocate(atm%dswd_dalb_vu_cs)
     deallocate(atm%dswd_dalb_ir_cs)
     deallocate(atm%dswd_dalb_vu_cld)
     deallocate(atm%dswd_dalb_ir_cld)
     deallocate(atm%dswd_dz_ir_cs)
     deallocate(atm%dswd_dz_ir_cld)
     deallocate(atm%swr_dw_sur_vis_cs)
     deallocate(atm%swr_dw_sur_nir_cs)
     deallocate(atm%swr_dw_sur_vis_cld)
     deallocate(atm%swr_dw_sur_nir_cld)

     deallocate(atm%rb_top)
     deallocate(atm%rb_sur)
     deallocate(atm%rb_atm)
     deallocate(atm%rb_str)     
     deallocate(atm%swr_dw_top)
     deallocate(atm%swr_top)
     deallocate(atm%swr_top_cs)
     deallocate(atm%swr_top_cld)
     deallocate(atm%swr_sur)
     deallocate(atm%lwr_top)
     deallocate(atm%lwr_top_cs)
     deallocate(atm%lwr_top_cld)
     deallocate(atm%lwr_sur)
     deallocate(atm%lwr_tro)
     deallocate(atm%lwr_cld)

     deallocate(atm%tsl)
     deallocate(atm%tsksl)

     deallocate(atm%eke)
     deallocate(atm%sam)
     deallocate(atm%sam2)
     deallocate(atm%synprod)
     deallocate(atm%syndiss)
     deallocate(atm%synadv)
     deallocate(atm%syndif)
     deallocate(atm%synsur)
     deallocate(atm%cdif)    

   return

  end subroutine atm_dealloc 


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_write_restart(fnm,atm)

    implicit none

    character (len=*) :: fnm
    type(atm_class) :: atm


    call nc_create(fnm)
    call nc_write_dim(fnm,"lon",x=lon,axis="x")
    call nc_write_dim(fnm,"lat",x=lat,axis="y")
    call nc_write_dim(fnm,"nm",x=1,dx=1,nx=nm)
    call nc_write_dim(fnm,"plev",x=pl(1:km),axis="z")
    call nc_write_dim(fnm,"plevc",x=pl(1:kmc),axis="z")
    call nc_write_dim(fnm,"x",x=[1])

    call nc_write(fnm,"had_fi   ",     atm%had_fi   ,     dim1="x",  long_name="",units="")
    call nc_write(fnm,"had_width",     atm%had_width,     dim1="x",  long_name="",units="")

    call nc_write(fnm,"tam      ",     atm%tam      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"t2a      ",     atm%t2a      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"tskina   ",     atm%tskina   ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"t2       ",     atm%t2       ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"tskin    ",     atm%tskin    ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"tsl      ",     atm%tsl      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"tsksl     ",     atm%tsksl     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"gams     ",     atm%gams     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"gamb     ",     atm%gamb     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"gamt     ",     atm%gamt     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"hrm      ",     atm%hrm      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"qam      ",     atm%qam      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"q2a      ",     atm%q2a      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"q2       ",     atm%q2       ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"ram      ",     atm%ram      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"rskina   ",     atm%rskina   ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"r2       ",     atm%r2       ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"r2a      ",     atm%r2a      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"hqeff    ",     atm%hqeff    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"wcon     ",     atm%wcon     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"prc      ",     atm%prc      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"cld      ",     atm%cld      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"cld_rh   ",     atm%cld_rh   ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"cld_low  ",     atm%cld_low  ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"hcld     ",     atm%hcld     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"clot     ",     atm%clot     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"dam      ",     atm%dam      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"hdust    ",     atm%hdust    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"dust_ot  ",     atm%dust_ot  ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"so4      ",     atm%so4      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"htrop    ",     atm%htrop    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"ttrop    ",     atm%ttrop    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"rb_sur   ",     atm%rb_sur   ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"sam      ",     atm%sam      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"sam2     ",     atm%sam2     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"cdif     ",     atm%cdif     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"slp      ",     atm%slp      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"ps       ",     atm%ps       ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"aslp     ",     atm%aslp     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"aslp_topo",     atm%aslp_topo,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"dz500    ",     atm%dz500    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"winda    ",     atm%winda    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"wind     ",     atm%wind     ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"us       ",     atm%us       ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"vs       ",     atm%vs       ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"usk      ",     atm%usk      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"vsk      ",     atm%vsk      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"uz500    ",     atm%uz500    ,     dims=["lat"],long_name="",units="")
    call nc_write(fnm,"wsyn     ",     atm%wsyn     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"woro     ",     atm%woro     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"wcld     ",     atm%wcld     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"weff     ",     atm%weff     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"fweff    ",     atm%fweff    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"convdse  ",     atm%convdse  ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"convwtr  ",     atm%convwtr  ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"ra2a     ",     atm%ra2a     ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"frlnd    ",     atm%frlnd    ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"f_ice_lake ",   atm%f_ice_lake,    dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"frst     ",     atm%frst     ,     dims=["lon","lat","nm "],long_name="",units="")
    call nc_write(fnm,"zsa      ",     atm%zsa      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"psa      ",     atm%psa      ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"sigoro   ",     atm%sigoro   ,     dims=["lon","lat"],long_name="",units="")
    call nc_write(fnm,"u3      ",     atm%u3      ,     dims=["lon ","lat ","plev"],long_name="",units="")
    call nc_write(fnm,"v3      ",     atm%v3      ,     dims=["lon ","lat ","plev"],long_name="",units="")
    call nc_write(fnm,"w3      ",     atm%w3      ,     dims=["lon  ","lat  ","plevc"],long_name="",units="")
    call nc_write(fnm,"t3      ",     atm%t3      ,     dims=["lon ","lat ","plev"],long_name="",units="")
    call nc_write(fnm,"tp      ",     atm%tp      ,     dims=["lon ","lat ","plev"],long_name="",units="")
    call nc_write(fnm,"q3      ",     atm%q3      ,     dims=["lon ","lat ","plev"],long_name="",units="")
    call nc_write(fnm,"d3      ",     atm%d3      ,     dims=["lon ","lat ","plev"],long_name="",units="")


   return

  end subroutine atm_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  a t m _ r e a d _ r e s t a r t
  ! Purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine atm_read_restart(fnm,atm)

    implicit none

    character (len=*) :: fnm
    type(atm_class) :: atm


    call nc_read(fnm,"had_fi   ",     atm%had_fi   )
    call nc_read(fnm,"had_width",     atm%had_width)

    call nc_read(fnm,"tam      ",     atm%tam     ) 
    call nc_read(fnm,"t2a      ",     atm%t2a     ) 
    call nc_read(fnm,"t2       ",     atm%t2      ) 
    call nc_read(fnm,"tskina   ",     atm%tskina  ) 
    call nc_read(fnm,"tskin    ",     atm%tskin   ) 
    call nc_read(fnm,"tsl      ",     atm%tsl     ) 
    call nc_read(fnm,"tsksl     ",     atm%tsksl    ) 
    call nc_read(fnm,"gams     ",     atm%gams    ) 
    call nc_read(fnm,"gamb     ",     atm%gamb    ) 
    call nc_read(fnm,"gamt     ",     atm%gamt    ) 
    call nc_read(fnm,"hrm      ",     atm%hrm     ) 
    call nc_read(fnm,"qam      ",     atm%qam     ) 
    call nc_read(fnm,"q2a      ",     atm%q2a     ) 
    call nc_read(fnm,"q2       ",     atm%q2      ) 
    call nc_read(fnm,"ram      ",     atm%ram     ) 
    call nc_read(fnm,"rskina   ",     atm%rskina  ) 
    call nc_read(fnm,"r2       ",     atm%r2      ) 
    call nc_read(fnm,"r2a      ",     atm%r2a     ) 
    call nc_read(fnm,"hqeff    ",     atm%hqeff   ) 
    call nc_read(fnm,"wcon     ",     atm%wcon    ) 
    call nc_read(fnm,"prc      ",     atm%prc     ) 
    call nc_read(fnm,"cld      ",     atm%cld     ) 
    call nc_read(fnm,"cld_rh   ",     atm%cld_rh  ) 
    call nc_read(fnm,"cld_low  ",     atm%cld_low ) 
    call nc_read(fnm,"hcld     ",     atm%hcld    ) 
    call nc_read(fnm,"clot     ",     atm%clot    ) 
    call nc_read(fnm,"dam      ",     atm%dam     ) 
    call nc_read(fnm,"hdust    ",     atm%hdust   ) 
    call nc_read(fnm,"dust_ot  ",     atm%dust_ot ) 
    call nc_read(fnm,"so4      ",     atm%so4     ) 
    call nc_read(fnm,"htrop    ",     atm%htrop   ) 
    call nc_read(fnm,"ttrop    ",     atm%ttrop   ) 
    call nc_read(fnm,"rb_sur   ",     atm%rb_sur  ) 
    call nc_read(fnm,"sam      ",     atm%sam     ) 
    call nc_read(fnm,"sam2     ",     atm%sam2    ) 
    call nc_read(fnm,"cdif     ",     atm%cdif    ) 
    call nc_read(fnm,"slp      ",     atm%slp     ) 
    call nc_read(fnm,"ps       ",     atm%ps      ) 
    call nc_read(fnm,"aslp     ",     atm%aslp    ) 
    call nc_read(fnm,"aslp_topo",     atm%aslp_topo) 
    call nc_read(fnm,"dz500    ",     atm%dz500   ) 
    call nc_read(fnm,"winda    ",     atm%winda   ) 
    call nc_read(fnm,"wind     ",     atm%wind    ) 
    call nc_read(fnm,"us       ",     atm%us      ) 
    call nc_read(fnm,"vs       ",     atm%vs      ) 
    call nc_read(fnm,"usk      ",     atm%usk     ) 
    call nc_read(fnm,"vsk      ",     atm%vsk     ) 
    call nc_read(fnm,"uz500    ",     atm%uz500   ) 
    call nc_read(fnm,"wsyn     ",     atm%wsyn    ) 
    call nc_read(fnm,"woro     ",     atm%woro    ) 
    call nc_read(fnm,"wcld     ",     atm%wcld    ) 
    call nc_read(fnm,"weff     ",     atm%weff    ) 
    call nc_read(fnm,"fweff    ",     atm%fweff   ) 
    call nc_read(fnm,"convdse  ",     atm%convdse ) 
    call nc_read(fnm,"convwtr  ",     atm%convwtr ) 
    call nc_read(fnm,"ra2a     ",     atm%ra2a    ) 
    call nc_read(fnm,"frlnd    ",     atm%frlnd   ) 
    call nc_read(fnm,"f_ice_lake   ",     atm%f_ice_lake  ) 
    call nc_read(fnm,"frst     ",     atm%frst    ) 
    call nc_read(fnm,"zsa      ",     atm%zsa     ) 
    call nc_read(fnm,"psa      ",     atm%psa     ) 
    call nc_read(fnm,"sigoro   ",     atm%sigoro  ) 
    call nc_read(fnm,"u3       ",     atm%u3      ) 
    call nc_read(fnm,"v3       ",     atm%v3      ) 
    call nc_read(fnm,"w3       ",     atm%w3      ) 
    call nc_read(fnm,"t3       ",     atm%t3      ) 
    call nc_read(fnm,"tp       ",     atm%tp      ) 
    call nc_read(fnm,"q3       ",     atm%q3      ) 
    call nc_read(fnm,"d3       ",     atm%d3      ) 

    print *,'read restart file ',fnm

   return

  end subroutine atm_read_restart


end module atm_model



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : f e e d b a c k s _ m o d
!
!  Purpose : feedback analysis
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
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
module feedbacks_mod

  use atm_params, only : wp
  use constants, only : sigma
  use control, only : out_dir
  use timer, only : doy, nday_year, time_eoy_atm
  use climber_grid, only : lon, lat
  use ncio
  use atm_grid, only : im, jm, nm, sqr, esqr
  use lw_radiation_mod, only : lw_radiation
  use sw_radiation_mod, only : sw_radiation
  use vesta_mod, only : tropoheight, t_prof

  implicit none

  integer, parameter :: i_control = 0
  integer, parameter :: i_pl  = 1
  integer, parameter :: i_wv  = 2
  integer, parameter :: i_cld = 3
  integer, parameter :: i_lr  = 4
  integer, parameter :: i_alb = 5
  integer, parameter :: i_cld_frac = 6
  integer, parameter :: i_cld_clot = 7
  integer, parameter :: i_cld_hcld = 8
  integer, parameter :: i_lr_gam  = 9
  integer, parameter :: i_lr_tam  = 10
  integer, parameter :: i_temp = 11
  integer, parameter :: i_all = 12
  integer, parameter :: nfb = 12

  type feedback_type
    real(wp) :: co2
    real(wp), dimension(:,:,:), allocatable :: tam
    real(wp), dimension(:,:,:), allocatable :: cld
    real(wp), dimension(:,:,:), allocatable :: hcld
    real(wp), dimension(:,:,:), allocatable :: clot
    real(wp), dimension(:,:,:), allocatable :: gams
    real(wp), dimension(:,:,:), allocatable :: gamb
    real(wp), dimension(:,:,:), allocatable :: gamt
    real(wp), dimension(:,:,:), allocatable :: htrop
    real(wp), dimension(:,:,:), allocatable :: ttrop
    real(wp), dimension(:,:,:), allocatable :: ram
    real(wp), dimension(:,:,:), allocatable :: hrm
    real(wp), dimension(:,:,:), allocatable :: hqeff
    real(wp), dimension(:,:,:), allocatable :: aerosol_ot 
    real(wp), dimension(:,:,:), allocatable :: aerosol_im
    real(wp), dimension(:,:,:), allocatable :: so4
    real(wp), dimension(:,:,:,:), allocatable :: frst
    real(wp), dimension(:,:,:,:), allocatable :: tskin
    real(wp), dimension(:,:,:,:), allocatable :: t2
    real(wp), dimension(:,:,:,:), allocatable :: q2
    real(wp), dimension(:,:,:,:), allocatable :: alb_vu_s
    real(wp), dimension(:,:,:,:), allocatable :: alb_vu_c
    real(wp), dimension(:,:,:,:), allocatable :: alb_ir_s
    real(wp), dimension(:,:,:,:), allocatable :: alb_ir_c
    real(wp), dimension(:,:,:,:), allocatable :: flwr_up_sur

    real(wp), dimension(2) :: tg
    real(wp) :: delta_t
    real(wp) :: rf_top_ave
    real(wp) :: rf_trop_ave
    real(wp), dimension(:,:), allocatable :: rf_top
    real(wp), dimension(:,:), allocatable :: rf_trop
    real(wp), dimension(:,:), allocatable :: dhtrop_rf
    real(wp), dimension(:,:,:), allocatable :: flwr_top
    real(wp), dimension(:,:,:), allocatable :: fswr_top
    real(wp), dimension(:,:,:), allocatable :: d_flwr_top
    real(wp), dimension(:,:,:), allocatable :: d_fswr_top
    real(wp), dimension(:,:,:), allocatable :: d_f_top
    real(wp), dimension(0:nfb) :: flwr_top_ave
    real(wp), dimension(0:nfb) :: fswr_top_ave
    real(wp), dimension(nfb) :: d_flwr_top_ave
    real(wp), dimension(nfb) :: d_fswr_top_ave
    real(wp), dimension(nfb) :: d_f_top_ave
  end type 

  private
  public :: feedback_type, feedback_init, feedback_save, feedback_analysis, feedback_write

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f e e d b a c k _ i n i t
  !   Purpose    :  initialisation for feedback analysis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine feedback_init(fb)

    implicit none

    type(feedback_type), intent(inout) :: fb


    allocate(fb%tam(im,jm,nday_year))
    allocate(fb%cld(im,jm,nday_year))
    allocate(fb%hcld(im,jm,nday_year))
    allocate(fb%clot(im,jm,nday_year))
    allocate(fb%gams(im,jm,nday_year))
    allocate(fb%gamb(im,jm,nday_year))
    allocate(fb%gamt(im,jm,nday_year))
    allocate(fb%htrop(im,jm,nday_year))
    allocate(fb%ttrop(im,jm,nday_year))
    allocate(fb%ram(im,jm,nday_year))
    allocate(fb%hrm(im,jm,nday_year))
    allocate(fb%hqeff(im,jm,nday_year))
    allocate(fb%aerosol_ot(im,jm,nday_year))
    allocate(fb%aerosol_im(im,jm,nday_year))
    allocate(fb%so4(im,jm,nday_year))
    allocate(fb%frst(im,jm,nm,nday_year))
    allocate(fb%tskin(im,jm,nm,nday_year))
    allocate(fb%t2(im,jm,nm,nday_year))
    allocate(fb%q2(im,jm,nm,nday_year))
    allocate(fb%alb_vu_s(im,jm,nm,nday_year))
    allocate(fb%alb_vu_c(im,jm,nm,nday_year))
    allocate(fb%alb_ir_s(im,jm,nm,nday_year))
    allocate(fb%alb_ir_c(im,jm,nm,nday_year))
    allocate(fb%flwr_up_sur(im,jm,nm,nday_year))

    allocate(fb%rf_top(im,jm))
    allocate(fb%rf_trop(im,jm))
    allocate(fb%dhtrop_rf(im,jm))
    allocate(fb%flwr_top(im,jm,0:nfb))
    allocate(fb%fswr_top(im,jm,0:nfb))
    allocate(fb%d_flwr_top(im,jm,nfb))
    allocate(fb%d_fswr_top(im,jm,nfb))
    allocate(fb%d_f_top(im,jm,nfb))

    fb%rf_top   = 0._wp
    fb%rf_trop  = 0._wp
    fb%dhtrop_rf  = 0._wp
    fb%flwr_top = 0._wp
    fb%fswr_top = 0._wp

    fb%tg = 0._wp

    return

  end subroutine feedback_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f e e d b a c k _ s a v e 
  !   Purpose    :  save variables needed for feedback analysis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine feedback_save(co2, tam, cld, hcld, clot, gams, gamb, gamt, htrop, ttrop, ram, hrm, hqeff, q2, aerosol_ot, aerosol_im, so4, &
      frst, tskin, t2, alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, flwr_up_sur, &
      fb)

    implicit none

    real(wp), intent(in) :: co2
    real(wp), dimension(:,:),   intent(in) :: tam
    real(wp), dimension(:,:),   intent(in) :: cld
    real(wp), dimension(:,:),   intent(in) :: hcld
    real(wp), dimension(:,:),   intent(in) :: clot
    real(wp), dimension(:,:),   intent(in) :: gams
    real(wp), dimension(:,:),   intent(in) :: gamb
    real(wp), dimension(:,:),   intent(in) :: gamt
    real(wp), dimension(:,:),   intent(in) :: htrop
    real(wp), dimension(:,:),   intent(in) :: ttrop
    real(wp), dimension(:,:),   intent(in) :: ram
    real(wp), dimension(:,:),   intent(in) :: hrm
    real(wp), dimension(:,:),   intent(in) :: hqeff
    real(wp), dimension(:,:,:), intent(in) :: q2
    real(wp), dimension(:,:),   intent(in) :: aerosol_ot 
    real(wp), dimension(:,:),   intent(in) :: aerosol_im
    real(wp), dimension(:,:),   intent(in) :: so4 
    real(wp), dimension(:,:,:), intent(in) :: frst
    real(wp), dimension(:,:,:), intent(in) :: tskin
    real(wp), dimension(:,:,:), intent(in) :: t2
    real(wp), dimension(:,:,:), intent(in) :: alb_vu_s
    real(wp), dimension(:,:,:), intent(in) :: alb_vu_c
    real(wp), dimension(:,:,:), intent(in) :: alb_ir_s
    real(wp), dimension(:,:,:), intent(in) :: alb_ir_c
    real(wp), dimension(:,:,:), intent(in) :: flwr_up_sur

    type(feedback_type), intent(inout) :: fb

    integer :: i, j, n, ncid
    character (len=256) :: fnm


    do i=1,im
      do j=1,jm
        fb%tam(i,j,doy)           = tam(i,j)
        fb%cld(i,j,doy)           = cld(i,j)
        fb%hcld(i,j,doy)          = hcld(i,j)
        fb%clot(i,j,doy)          = clot(i,j)
        fb%gams(i,j,doy)          = gams(i,j)
        fb%gamb(i,j,doy)          = gamb(i,j)
        fb%gamt(i,j,doy)          = gamt(i,j)
        fb%htrop(i,j,doy)         = htrop(i,j)
        fb%ttrop(i,j,doy)         = ttrop(i,j)
        fb%ram(i,j,doy)           = ram(i,j)
        fb%hrm(i,j,doy)           = hrm(i,j)
        fb%hqeff(i,j,doy)         = hqeff(i,j)
        fb%aerosol_ot(i,j,doy)    = aerosol_ot(i,j)
        fb%aerosol_im(i,j,doy)    = aerosol_im(i,j)
        fb%so4(i,j,doy)           = so4(i,j)
        do n=1,nm
          fb%frst(i,j,n,doy)        = frst(i,j,n)
          fb%tskin(i,j,n,doy)       = tskin(i,j,n)
          fb%t2(i,j,n,doy)          = t2(i,j,n)
          fb%q2(i,j,n,doy)          = q2(i,j,n)
          fb%alb_vu_s(i,j,n,doy)    = alb_vu_s(i,j,n)
          fb%alb_vu_c(i,j,n,doy)    = alb_vu_c(i,j,n)
          fb%alb_ir_s(i,j,n,doy)    = alb_ir_s(i,j,n)
          fb%alb_ir_c(i,j,n,doy)    = alb_ir_c(i,j,n)
          fb%flwr_up_sur(i,j,n,doy) = flwr_up_sur(i,j,n)
        enddo
      enddo
    enddo

    fb%tg(1) = fb%tg(1) + sum(sum(t2*frst,3)*sqr)/(esqr*nday_year)

    if (time_eoy_atm) then
      fb%co2         = co2          

      fnm = trim(out_dir)//"/feedback_save.nc"
      call nc_create(fnm)
      call nc_open(fnm,ncid)
      call nc_write_dim(fnm,"doy",x=1._wp,dx=1._wp,nx=nday_year,ncid=ncid)
      call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
      call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)
      call nc_write_dim(fnm,"st",x=1._wp,dx=1._wp,nx=nm,ncid=ncid)

      call nc_write(fnm,"tam        ", sngl(fb%tam         ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"cld        ", sngl(fb%cld         ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"hcld       ", sngl(fb%hcld        ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"clot       ", sngl(fb%clot        ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"gams       ", sngl(fb%gams        ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"gamb       ", sngl(fb%gamb        ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"gamt       ", sngl(fb%gamt        ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"htrop      ", sngl(fb%htrop       ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"ttrop      ", sngl(fb%ttrop       ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"ram        ", sngl(fb%ram         ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"hrm        ", sngl(fb%hrm         ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"hqeff      ", sngl(fb%hqeff       ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"aerosol_ot ", sngl(fb%aerosol_ot  ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"aerosol_im ", sngl(fb%aerosol_im  ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"so4        ", sngl(fb%so4         ), dims=["lon","lat","doy"],start=[1,1,1],count=[im,jm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"frst       ", sngl(fb%frst        ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"tskin      ", sngl(fb%tskin       ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"t2         ", sngl(fb%t2          ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"alb_vu_s   ", sngl(fb%alb_vu_s    ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"alb_vu_c   ", sngl(fb%alb_vu_c    ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"alb_ir_s   ", sngl(fb%alb_ir_s    ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"alb_ir_c   ", sngl(fb%alb_ir_c    ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)
      call nc_write(fnm,"flwr_up_sur", sngl(fb%flwr_up_sur ), dims=["lon","lat","st ","doy"],start=[1,1,1,1],count=[im,jm,nm,nday_year],long_name="",units="",ncid=ncid)

      call nc_close(ncid)

    endif

    return

  end subroutine feedback_save


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f e e d b a c k _ a n a l y s i s 
  !   Purpose    :  perform feedback analysis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine feedback_analysis(fb, frst, zs, zsa, htrop, hcld, tskin, t2, ra2, & 
      gams, gamb, gamt, tam, ram, hrm, hqeff, q2, ttrop, cld, co2, ch4, n2o, cfc11, cfc12, o3, flwr_up_sur, &  
      swr_dw_top, coszm, alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, aerosol_ot, aerosol_im, so4, had_fi, had_width)

    implicit none

    type(feedback_type), intent(inout) :: fb

    real(wp), dimension(:,:,:), intent(in) :: zs
    real(wp), dimension(:,:),   intent(in) :: zsa
    real(wp), dimension(:,:,:), intent(in) :: ra2
    real(wp), dimension(:,:),   intent(in) :: swr_dw_top
    real(wp), dimension(:,:),   intent(in) :: coszm
    real(wp), intent(in) :: co2
    real(wp), intent(in) :: ch4
    real(wp), intent(in) :: n2o
    real(wp), intent(in) :: cfc11
    real(wp), intent(in) :: cfc12
    real(wp), dimension(:,:,:), intent(in) :: o3
    real(wp), dimension(:,:),   intent(in) :: tam
    real(wp), dimension(:,:),   intent(in) :: cld
    real(wp), dimension(:,:),   intent(in) :: hcld
    real(wp), dimension(:,:),   intent(in) :: clot
    real(wp), dimension(:,:),   intent(in) :: gams
    real(wp), dimension(:,:),   intent(in) :: gamb
    real(wp), dimension(:,:),   intent(in) :: gamt
    real(wp), dimension(:,:),   intent(in) :: htrop
    real(wp), dimension(:,:),   intent(in) :: ttrop
    real(wp), dimension(:,:),   intent(in) :: ram
    real(wp), dimension(:,:),   intent(in) :: hrm
    real(wp), dimension(:,:),   intent(in) :: hqeff
    real(wp), dimension(:,:,:), intent(in) :: q2
    real(wp), dimension(:,:),   intent(in) :: aerosol_ot 
    real(wp), dimension(:,:),   intent(in) :: aerosol_im
    real(wp), dimension(:,:),   intent(in) :: so4 
    real(wp), dimension(:,:,:), intent(in) :: frst
    real(wp), dimension(:,:,:), intent(in) :: tskin
    real(wp), dimension(:,:,:), intent(in) :: t2
    real(wp), dimension(:,:,:), intent(in) :: alb_vu_s
    real(wp), dimension(:,:,:), intent(in) :: alb_vu_c
    real(wp), dimension(:,:,:), intent(in) :: alb_ir_s
    real(wp), dimension(:,:,:), intent(in) :: alb_ir_c
    real(wp), dimension(:,:,:), intent(in) :: flwr_up_sur
    real(wp), intent(in) :: had_fi
    real(wp), intent(in) :: had_width

    real(wp), dimension(:,:), allocatable :: lwr_sur
    real(wp), dimension(:,:,:), allocatable :: flwr_dw_sur
    real(wp), dimension(:,:,:), allocatable :: flwr_dw_sur_cs
    real(wp), dimension(:,:,:), allocatable :: flwr_dw_sur_cld
    real(wp), dimension(:,:), allocatable :: lwr_top
    real(wp), dimension(:,:), allocatable :: lwr_top_cs
    real(wp), dimension(:,:), allocatable :: lwr_top_cld
    real(wp), dimension(:,:), allocatable :: lwr_tro
    real(wp), dimension(:,:), allocatable :: lwr_cld
    real(wp), dimension(:,:),   allocatable :: alb_cld
    real(wp), dimension(:,:), allocatable :: swr_top
    real(wp), dimension(:,:), allocatable :: swr_top_cs
    real(wp), dimension(:,:), allocatable :: swr_top_cld
    real(wp), dimension(:,:), allocatable :: swr_sur
    real(wp), dimension(:,:,:), allocatable :: fswr_sur
    real(wp), dimension(:,:,:), allocatable :: fswr_sur_cs
    real(wp), dimension(:,:,:), allocatable :: fswr_sur_cld
    real(wp), dimension(:,:), allocatable :: flwr_trop_control
    real(wp), dimension(:,:), allocatable :: flwr_top_control
    real(wp), dimension(:,:), allocatable :: rb_str 
    real(wp), dimension(:,:), allocatable :: fb_htrop
    real(wp), dimension(:), allocatable :: fb_ptrop
    real(wp), dimension(:,:), allocatable :: fb_ttrop
    real(wp), dimension(:,:,:), allocatable :: fb_flwr_up_sur

    integer :: i, j, n
    real(wp) :: co2e


    allocate(lwr_sur(im,jm))
    allocate(flwr_dw_sur(im,jm,nm))
    allocate(flwr_dw_sur_cs(im,jm,nm))
    allocate(flwr_dw_sur_cld(im,jm,nm))
    allocate(lwr_top(im,jm))
    allocate(lwr_top_cs (im,jm))
    allocate(lwr_top_cld (im,jm))
    allocate(lwr_tro(im,jm))
    allocate(lwr_cld(im,jm))
    allocate(alb_cld    (im,jm))
    allocate(swr_top   (im,jm))
    allocate(swr_top_cs (im,jm))
    allocate(swr_top_cld (im,jm))
    allocate(swr_sur   (im,jm))
    allocate(fswr_sur   (im,jm,nm))
    allocate(fswr_sur_cs (im,jm,nm))
    allocate(fswr_sur_cld (im,jm,nm))
    allocate(flwr_trop_control(im,jm))
    allocate(flwr_top_control(im,jm))
    allocate(rb_str(im,jm))
    allocate(fb_htrop(im,jm))
    allocate(fb_ptrop(jm))
    allocate(fb_ttrop(im,jm))
    allocate(fb_flwr_up_sur(im,jm,nm))

    !-------------------------------
    ! longwave radiation
    !-------------------------------

    ! control 
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    fb%flwr_top(:,:,i_control) = fb%flwr_top(:,:,i_control) + lwr_top/nday_year
    flwr_top_control = lwr_top     ! net flux at TOA, positive down
    flwr_trop_control = lwr_tro    ! net flux at tropopause, positive down

    ! radiative forcing
    ! longwave radiation with 2xCO2 to diagnose stratosphere fluxes 
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, &    ! in
      fb%co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    ! radiative forcing at the top of atmosphere 
    fb%rf_top(:,:) = fb%rf_top(:,:) + (lwr_top-flwr_top_control)/nday_year
    ! stratosphere radiative balance for 2xCO2
    rb_str = lwr_top-lwr_tro
    ! compute implied new tropopause height
    fb_htrop = htrop
    do n=1,10
      call tropoheight(had_fi, had_width, rb_str, hcld, fb_htrop, fb_ptrop)
    enddo
    fb%dhtrop_rf(:,:) = fb%dhtrop_rf(:,:) + (fb_htrop-htrop)/nday_year
    ! compute ttrop for new tropopause height
    do i=1,im
      do j=1,jm
        fb_ttrop(i,j) = t_prof(zsa(i,j), fb_htrop(i,j), tam(i,j), gams(i,j), gamb(i,j), gamt(i,j), fb_htrop(i,j), 1)
      enddo
    enddo
    ! longwave radiation with 2xCO2 and adjusted stratospheric temperature (htrop and ttrop)
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, fb_ttrop, cld, clot, & ! in
      fb%co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &   ! out
      gams, gamb, gamt, tam, ttrop, htrop) ! optional input arguments for feedback (moisture)
    fb%rf_trop(:,:) = fb%rf_trop(:,:) + (lwr_tro-flwr_trop_control)/nday_year

    ! Planck feedback
    fb_flwr_up_sur = sigma*((flwr_up_sur/sigma)**0.25 + (fb%tskin(:,:,:,doy)-tskin))**4
    ! uniform tropospheric warming (same warming as near surface air)
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam+sum((fb%t2(:,:,:,doy)-t2)*frst,3), ram, hrm, ttrop, cld, clot, &  ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, fb_flwr_up_sur, &  ! in, account also for surface emission changes
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &    ! out
      gams, gamb, gamt, tam, ttrop, htrop) ! optional input arguments for feedback (moisture)
    fb%flwr_top(:,:,i_pl) = fb%flwr_top(:,:,i_pl) + lwr_top/nday_year

    ! water vapor feedback
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, fb%ram(:,:,doy), fb%hrm(:,:,doy), ttrop, cld, clot, &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &    ! out
      fb%gams(:,:,doy), fb%gamb(:,:,doy), fb%gamt(:,:,doy), fb%tam(:,:,doy), ttrop, htrop) ! optional input arguments for feedback (moisture)
    fb%flwr_top(:,:,i_wv) = fb%flwr_top(:,:,i_wv) + lwr_top/nday_year

    ! cloud feedback
    call lw_radiation(1._wp, frst, zsa, zs, htrop, fb%hcld(:,:,doy), ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, fb%cld(:,:,doy), fb%clot(:,:,doy), &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    fb%flwr_top(:,:,i_cld) = fb%flwr_top(:,:,i_cld) + lwr_top/nday_year

    ! cloud fraction feedback
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, fb%cld(:,:,doy), clot, &  ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    fb%flwr_top(:,:,i_cld_frac) = fb%flwr_top(:,:,i_cld_frac) + lwr_top/nday_year

    ! cloud top height feedback
    call lw_radiation(1._wp, frst, zsa, zs, htrop, fb%hcld(:,:,doy), ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, &  ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    fb%flwr_top(:,:,i_cld_hcld) = fb%flwr_top(:,:,i_cld_hcld) + lwr_top/nday_year

    ! cloud optical thickness feedback
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, fb%clot(:,:,doy), &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    fb%flwr_top(:,:,i_cld_clot) = fb%flwr_top(:,:,i_cld_clot) + lwr_top/nday_year

    ! lapse rate feedback
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      fb%gams(:,:,doy), fb%gamb(:,:,doy), fb%gamt(:,:,doy), fb%tam(:,:,doy)-sum((fb%t2(:,:,:,doy)-t2)*frst,3), ram, hrm, ttrop, cld, clot, & ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &    ! out
      gams, gamb, gamt, tam, ttrop, htrop) ! optional input arguments for feedback analysis (moisture profile)
    fb%flwr_top(:,:,i_lr) = fb%flwr_top(:,:,i_lr) + lwr_top/nday_year

    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      fb%gams(:,:,doy), fb%gamb(:,:,doy), fb%gamt(:,:,doy), tam, ram, hrm, ttrop, cld, clot, &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &    ! out
      gams, gamb, gamt, tam, ttrop, htrop) ! optional input arguments for feedback analysis (moisture profile)
    fb%flwr_top(:,:,i_lr_gam) = fb%flwr_top(:,:,i_lr_gam) + lwr_top/nday_year

    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, fb%tam(:,:,doy)-sum((fb%t2(:,:,:,doy)-t2)*frst,3), ram, hrm, ttrop, cld, clot, &  ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &    ! out
      gams, gamb, gamt, tam, ttrop, htrop) ! optional input arguments for feedback analysis (moisture profile)
    fb%flwr_top(:,:,i_lr_tam) = fb%flwr_top(:,:,i_lr_tam) + lwr_top/nday_year

    ! total temperature feedback
    fb_flwr_up_sur = sigma*((flwr_up_sur/sigma)**0.25 + (fb%tskin(:,:,:,doy)-tskin))**4 ! assuming unit emissivity
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      fb%gams(:,:,doy), fb%gamb(:,:,doy), fb%gamt(:,:,doy), fb%tam(:,:,doy), ram, hrm, ttrop, cld, clot, &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, fb_flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &    ! out
      gams, gamb, gamt, tam, ttrop, htrop) ! optional input arguments for feedback analysis (moisture profile)
    fb%flwr_top(:,:,i_temp) = fb%flwr_top(:,:,i_temp) + lwr_top/nday_year

    ! albedo feedback
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    fb%flwr_top(:,:,i_alb) = fb%flwr_top(:,:,i_alb) + lwr_top/nday_year

    ! all longwave feedbacks
    fb_flwr_up_sur = sigma*((flwr_up_sur/sigma)**0.25 + (fb%tskin(:,:,:,doy)-tskin))**4 ! assuming unit emissivity
    call lw_radiation(1._wp, frst, zsa, zs, htrop, fb%hcld(:,:,doy), ra2, &   ! in
      fb%gams(:,:,doy), fb%gamb(:,:,doy), fb%gamt(:,:,doy), fb%tam(:,:,doy), fb%ram(:,:,doy), fb%hrm(:,:,doy), ttrop, fb%cld(:,:,doy), fb%clot(:,:,doy), &    ! in
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, fb_flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &    ! out
      fb%gams(:,:,doy), fb%gamb(:,:,doy), fb%gamt(:,:,doy), fb%tam(:,:,doy), ttrop, htrop) ! optional input arguments for feedback (moisture)
    fb%flwr_top(:,:,i_all) = fb%flwr_top(:,:,i_all) + lwr_top/nday_year


    !-------------------------------
    ! shortwave radiation
    !-------------------------------

    ! control
    call sw_radiation(frst, swr_dw_top, coszm, cld, q2, ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, hcld, hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_control) = fb%fswr_top(:,:,i_control) + swr_top/nday_year

    ! Planck feedback 
    fb%fswr_top(:,:,i_pl) = fb%fswr_top(:,:,i_control) 

    ! water vapor feedback 
    call sw_radiation(frst, swr_dw_top, coszm, cld, fb%q2(:,:,:,doy), ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, hcld, fb%hqeff(:,:,doy), aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_wv) = fb%fswr_top(:,:,i_wv) + swr_top/nday_year

    ! cloud feedback
    call sw_radiation(frst, swr_dw_top, coszm, fb%cld(:,:,doy), q2, ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, fb%clot(:,:,doy), hcld, hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_cld) = fb%fswr_top(:,:,i_cld) + swr_top/nday_year

    ! cloud fraction feedback
    call sw_radiation(frst, swr_dw_top, coszm, fb%cld(:,:,doy), q2, ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, hcld, hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_cld_frac) = fb%fswr_top(:,:,i_cld_frac) + swr_top/nday_year

    ! cloud optical thickness feedback
    call sw_radiation(frst, swr_dw_top, coszm, cld, q2, ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, fb%clot(:,:,doy), hcld, hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_cld_clot) = fb%fswr_top(:,:,i_cld_clot) + swr_top/nday_year

    ! cloud top height feedback
    call sw_radiation(frst, swr_dw_top, coszm, cld, q2, ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, fb%hcld(:,:,doy), hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_cld_hcld) = fb%fswr_top(:,:,i_cld_hcld) + swr_top/nday_year

    ! lapse rate feedback
    fb%fswr_top(:,:,i_lr) = fb%fswr_top(:,:,i_control) 
    fb%fswr_top(:,:,i_lr_gam) = fb%fswr_top(:,:,i_control) 
    fb%fswr_top(:,:,i_lr_tam) = fb%fswr_top(:,:,i_control) 

    ! total temperature feedback
    fb%fswr_top(:,:,i_temp) = fb%fswr_top(:,:,i_control) 

    ! albedo feedback
    call sw_radiation(fb%frst(:,:,:,doy), swr_dw_top, coszm, cld, q2, ra2, &    ! in
      fb%alb_vu_s(:,:,:,doy), fb%alb_vu_c(:,:,:,doy), fb%alb_ir_s(:,:,:,doy), fb%alb_ir_c(:,:,:,doy), clot, hcld, hqeff, aerosol_ot, aerosol_im, so4, & ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_alb) = fb%fswr_top(:,:,i_alb) + swr_top/nday_year

    ! all SW feedbacks
    call sw_radiation(fb%frst(:,:,:,doy), swr_dw_top, coszm, fb%cld(:,:,doy), fb%q2(:,:,:,doy), ra2, &    ! in
      fb%alb_vu_s(:,:,:,doy), fb%alb_vu_c(:,:,:,doy), fb%alb_ir_s(:,:,:,doy), fb%alb_ir_c(:,:,:,doy), &  ! in
      fb%clot(:,:,doy), fb%hcld(:,:,doy), fb%hqeff(:,:,doy), aerosol_ot, aerosol_im, so4, & ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    fb%fswr_top(:,:,i_all) = fb%fswr_top(:,:,i_all) + swr_top/nday_year



    ! global mean temperature
    fb%tg(2) = fb%tg(2) + sum(sum(t2*frst,3)*sqr)/(esqr*nday_year)


    if (time_eoy_atm) then

      fb%delta_t = fb%tg(1)-fb%tg(2)

      fb%rf_top_ave = sum(fb%rf_top(:,:)*sqr)/esqr
      fb%rf_trop_ave = sum(fb%rf_trop(:,:)*sqr)/esqr

      do n=0,nfb
        fb%flwr_top_ave(n) = sum(fb%flwr_top(:,:,n)*sqr)/esqr
        fb%fswr_top_ave(n) = sum(fb%fswr_top(:,:,n)*sqr)/esqr
      enddo

      do n=1,nfb
        fb%d_flwr_top_ave(n) = (fb%flwr_top_ave(n)-fb%flwr_top_ave(i_control)) / fb%delta_t
        fb%d_fswr_top_ave(n) = (fb%fswr_top_ave(n)-fb%fswr_top_ave(i_control)) / fb%delta_t
        fb%d_f_top_ave(n) = fb%d_flwr_top_ave(n) + fb%d_fswr_top_ave(n)
      enddo

      do n=1,nfb
        fb%d_flwr_top(:,:,n) = (fb%flwr_top(:,:,n)-fb%flwr_top(:,:,i_control)) / fb%delta_t
        fb%d_fswr_top(:,:,n) = (fb%fswr_top(:,:,n)-fb%fswr_top(:,:,i_control)) / fb%delta_t
        fb%d_f_top(:,:,n)    = fb%d_flwr_top(:,:,n)+fb%d_fswr_top(:,:,n) 
      enddo

    endif

    deallocate(lwr_sur)
    deallocate(flwr_dw_sur)
    deallocate(flwr_dw_sur_cs)
    deallocate(flwr_dw_sur_cld)
    deallocate(lwr_top)
    deallocate(lwr_top_cs )
    deallocate(lwr_top_cld )
    deallocate(lwr_tro)
    deallocate(lwr_cld)
    deallocate(swr_top   )
    deallocate(swr_top_cs )
    deallocate(swr_top_cld )
    deallocate(swr_sur   )
    deallocate(fswr_sur   )
    deallocate(fswr_sur_cs   )
    deallocate(fswr_sur_cld   )
    deallocate(flwr_trop_control)
    deallocate(flwr_top_control)
    deallocate(rb_str)
    deallocate(fb_htrop)
    deallocate(fb_ttrop)

    return

  end subroutine feedback_analysis


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  f e e d b a c k _ w r i t e
  !   Purpose    :  write results of feedback analysis
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine feedback_write(fb)

    implicit none

    type(feedback_type), intent(in) :: fb

    integer :: ncid
    character (len=256) :: fnm


    ! write results of feedback analysis
    fnm = trim(out_dir)//"/feedbacks.nc"
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,"c",x=1,ncid=ncid)
    call nc_write_dim(fnm,"fb",x=1,dx=1,nx=nfb,ncid=ncid)
    call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)

    call nc_write(fnm,'i_pl      ', i_pl      ,dims=["c"],ncid=ncid)    
    call nc_write(fnm,'i_wv      ', i_wv      ,dims=["c"],ncid=ncid) 
    call nc_write(fnm,'i_cld     ', i_cld     ,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_lr      ', i_lr      ,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_alb     ', i_alb     ,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_cld_frac', i_cld_frac,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_cld_clot', i_cld_clot,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_cld_hcld', i_cld_hcld,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_lr_gam  ', i_lr_gam  ,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_lr_tam  ', i_lr_tam  ,dims=["c"],ncid=ncid)
    call nc_write(fnm,'i_temp    ', i_temp    ,dims=["c"],ncid=ncid) 
    call nc_write(fnm,'i_all     ', i_all     ,dims=["c"],ncid=ncid) 

    call nc_write(fnm,'ecs',fb%delta_t,dims=["c"],long_name="equilibrium climate sensitivity",units="K",ncid=ncid)
    call nc_write(fnm,'rf_trop',fb%rf_trop_ave,dims=["c"],long_name="radiative forcing at the tropopause for 2xCO2",units="W/m2",ncid=ncid)
    call nc_write(fnm,'rf_top',fb%rf_top_ave,dims=["c"],long_name="radiative forcing at the top of atmosphere for 2xCO2",units="W/m2",ncid=ncid)
    call nc_write(fnm,'fb_glob',fb%d_f_top_ave(1:nfb),dims=["fb"],long_name="global feedback factors [PLANCK,WV,CLD,LR,ALB]",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,'fb_glob_lw',fb%d_flwr_top_ave(1:nfb),dims=["fb"], &
      long_name="global longwave feedback factors [PLANCK,WV,CLD,LR,ALB]",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,'fb_glob_sw',fb%d_fswr_top_ave(1:nfb),dims=["fb"], &
      long_name="global shortwave feedback factors [PLANCK,WV,CLD,LR,ALB]",units="W/m2/K",ncid=ncid)

    call nc_write(fnm,"dhtrop_rf_2d", sngl(fb%dhtrop_rf(:,jm:1:-1)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="change in tropopause height due to CO2 radiative forcing",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rf_trop_2d", sngl(fb%rf_trop(:,jm:1:-1)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="radiative forcing at the tropopause",units="W/m2",ncid=ncid)
    call nc_write(fnm,"rf_top_2d", sngl(fb%rf_top(:,jm:1:-1)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="radiative forcing at the top of atmosphere",units="W/m2",ncid=ncid)

    call nc_write(fnm,"fb2d_pl", sngl(fb%d_f_top(:,jm:1:-1,i_pl)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="Planck feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_wv", sngl(fb%d_f_top(:,jm:1:-1,i_wv)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="water vapor feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_wv_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_wv)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave water vapor feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_wv_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_wv)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave water vapor feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld", sngl(fb%d_f_top(:,jm:1:-1,i_cld)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="cloud feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_cld)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave cloud feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_frac_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_cld_frac)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave cloud feedback (cloud fraction only)",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_clot_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_cld_clot)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave cloud feedback (optical thickness only)",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_hcld_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_cld_hcld)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave cloud feedback (cloud height only)",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_cld)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave cloud feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_frac_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_cld_frac)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave cloud feedback (cloud fraction only)",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_clot_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_cld_clot)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave cloud feedback (optical thickness only)",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_cld_hcld_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_cld_hcld)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave cloud feedback (cloud height only)",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr", sngl(fb%d_f_top(:,jm:1:-1,i_lr)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="lapse rate feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_lr)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave lapse rate feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_lr)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave lapse rate feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_gam", sngl(fb%d_f_top(:,jm:1:-1,i_lr_gam)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="lapse rate feedback from changes in gamma only",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_gam_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_lr_gam)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave lapse rate feedback from changes in gamma only",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_gam_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_lr_gam)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave lapse rate feedback from changes in gamma only",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_tam", sngl(fb%d_f_top(:,jm:1:-1,i_lr_tam)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="lapse rate feedback from changes in tam only",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_tam_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_lr_tam)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave lapse rate feedback from changes in tam only",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_lr_tam_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_lr_tam)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave lapse rate feedback from changes in tam only",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_temp", sngl(fb%d_f_top(:,jm:1:-1,i_temp)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="total temperature feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_alb", sngl(fb%d_f_top(:,jm:1:-1,i_alb)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="albedo feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_alb_lw", sngl(fb%d_flwr_top(:,jm:1:-1,i_alb)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="long wave albedo feedback",units="W/m2/K",ncid=ncid)
    call nc_write(fnm,"fb2d_alb_sw", sngl(fb%d_fswr_top(:,jm:1:-1,i_alb)), dims=["lon","lat"],start=[1,1],count=[im,jm], &
      long_name="short wave albedo feedback",units="W/m2/K",ncid=ncid)

    call nc_close(ncid)

    return

  end subroutine feedback_write

end module feedbacks_mod


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : r a d _ k e r n e l s _ m o d
!
!  Purpose : computation of radiative kernels
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
module rad_kernels_mod

  use atm_params, only : wp
  use constants, only : sigma
  use control, only : out_dir
  use timer, only : doy, nday_year, time_eoy_atm
  use climber_grid, only : lon, lat
  use atm_grid, only : im, jm, nm, sqr, esqr
  use lw_radiation_mod, only : lw_radiation
  use sw_radiation_mod, only : sw_radiation
  use ncio

  implicit none

  integer, parameter :: i_sw_control = 0
  integer, parameter :: i_sw_wv      = 1
  integer, parameter :: i_sw_alb     = 2
  integer, parameter :: n_sw = 2

  integer, parameter :: i_lw_control = 0
  integer, parameter :: i_lw_wv      = 1
  integer, parameter :: i_lw_Ta      = 2
  integer, parameter :: i_lw_Ts      = 3
  integer, parameter :: i_lw_co2     = 4
  integer, parameter :: n_lw = 4

  type rad_kernels_type
    real(wp), dimension(:,:,:), allocatable :: flwr_top
    real(wp), dimension(:,:,:), allocatable :: fswr_top
    real(wp), dimension(:,:,:), allocatable :: flwr_top_cs
    real(wp), dimension(:,:,:), allocatable :: fswr_top_cs
    real(wp), dimension(:,:,:), allocatable :: flwr_sur 
    real(wp), dimension(:,:,:), allocatable :: fswr_sur 
    real(wp), dimension(:,:,:), allocatable :: flwr_sur_cs
    real(wp), dimension(:,:,:), allocatable :: fswr_sur_cs
  end type 

  private
  public :: rad_kernels_type, rad_kernels_init, rad_kernels, rad_kernels_write

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r a d _ k e r n e l s _ i n i t
  !   Purpose    :  initialisation for radiative kernel computation 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rad_kernels_init(rk)

    implicit none

    type(rad_kernels_type), intent(inout) :: rk


    allocate(rk%flwr_top(im,jm,0:n_lw))
    allocate(rk%fswr_top(im,jm,0:n_sw))
    allocate(rk%flwr_top_cs(im,jm,0:n_lw))
    allocate(rk%fswr_top_cs(im,jm,0:n_sw))
    allocate(rk%flwr_sur (im,jm,0:n_lw))
    allocate(rk%fswr_sur (im,jm,0:n_sw))
    allocate(rk%flwr_sur_cs(im,jm,0:n_lw))
    allocate(rk%fswr_sur_cs(im,jm,0:n_sw))

    rk%flwr_top    = 0._wp
    rk%fswr_top    = 0._wp
    rk%flwr_top_cs = 0._wp
    rk%fswr_top_cs = 0._wp
    rk%flwr_sur    = 0._wp
    rk%fswr_sur    = 0._wp
    rk%flwr_sur_cs = 0._wp
    rk%fswr_sur_cs = 0._wp


    return

  end subroutine rad_kernels_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r a d _ k e r n e l s 
  !   Purpose    :  derive radiative kernels
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rad_kernels(rk, frst, zs, zsa, htrop, hcld, ra2, & 
      gams, gamb, gamt, tam, ram, hrm, hqeff, q2, ttrop, cld, co2, ch4, n2o, cfc11, cfc12, o3, flwr_up_sur, &  
      swr_dw_top, coszm, alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, aerosol_ot, aerosol_im, so4)

    implicit none

    type(rad_kernels_type), intent(inout) :: rk

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
    real(wp), dimension(:,:,:), intent(in) :: alb_vu_s
    real(wp), dimension(:,:,:), intent(in) :: alb_vu_c
    real(wp), dimension(:,:,:), intent(in) :: alb_ir_s
    real(wp), dimension(:,:,:), intent(in) :: alb_ir_c
    real(wp), dimension(:,:,:), intent(in) :: flwr_up_sur

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
    real(wp), dimension(:,:,:), allocatable :: rk_flwr_up_sur

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
    allocate(rk_flwr_up_sur(im,jm,nm))

    !------------------------
    ! longwave radiation
    !------------------------

    ! control 
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld)    ! out
    rk%flwr_top(:,:,i_lw_control)    = rk%flwr_top(:,:,i_lw_control)    + lwr_top  / nday_year
    rk%flwr_sur(:,:,i_lw_control)    = rk%flwr_sur(:,:,i_lw_control)    + lwr_sur   / nday_year
    rk%flwr_top_cs(:,:,i_lw_control) = rk%flwr_top_cs(:,:,i_lw_control) + lwr_top_cs / nday_year
    rk%flwr_sur_cs(:,:,i_lw_control) = rk%flwr_sur_cs(:,:,i_lw_control) + sum((flwr_dw_sur_cs-flwr_up_sur)*frst,3) / nday_year

    ! water vapor, increase corresponding to 1K temperature increase
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, & ! out
      gams, gamb, gamt, tam+1._wp, ttrop, htrop) ! optional input arguments for feedback (moisture)
    rk%flwr_top(:,:,i_lw_wv)    = rk%flwr_top(:,:,i_lw_wv)    + lwr_top  / nday_year
    rk%flwr_sur(:,:,i_lw_wv)    = rk%flwr_sur(:,:,i_lw_wv)    + lwr_sur   / nday_year
    rk%flwr_top_cs(:,:,i_lw_wv) = rk%flwr_top_cs(:,:,i_lw_wv) + lwr_top_cs / nday_year
    rk%flwr_sur_cs(:,:,i_lw_wv) = rk%flwr_sur_cs(:,:,i_lw_wv) + sum((flwr_dw_sur_cs-flwr_up_sur)*frst,3) / nday_year

    ! air temperature, 1 K temperature increase
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam+1._wp, ram, hrm, ttrop, cld, clot, co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, & ! out
      gams, gamb, gamt, tam, ttrop, htrop) ! optional input arguments for feedback (moisture)
    rk%flwr_top(:,:,i_lw_Ta)    = rk%flwr_top(:,:,i_lw_Ta)    + lwr_top  / nday_year
    rk%flwr_sur(:,:,i_lw_Ta)    = rk%flwr_sur(:,:,i_lw_Ta)    + lwr_sur   / nday_year
    rk%flwr_top_cs(:,:,i_lw_Ta) = rk%flwr_top_cs(:,:,i_lw_Ta) + lwr_top_cs / nday_year
    rk%flwr_sur_cs(:,:,i_lw_Ta) = rk%flwr_sur_cs(:,:,i_lw_Ta) + sum((flwr_dw_sur_cs-flwr_up_sur)*frst,3) / nday_year

    ! surface temperature, 1 K temperature increase (assuming unit surface emissivity)
    rk_flwr_up_sur = sigma*((flwr_up_sur/sigma)**0.25 + 1._wp)**4
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, co2, ch4, n2o, cfc11, cfc12, co2e, o3, rk_flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld) ! out
    rk%flwr_top(:,:,i_lw_Ts)    = rk%flwr_top(:,:,i_lw_Ts)    + lwr_top  / nday_year
    rk%flwr_sur(:,:,i_lw_Ts)    = rk%flwr_sur(:,:,i_lw_Ts)    + lwr_sur   / nday_year
    rk%flwr_top_cs(:,:,i_lw_Ts) = rk%flwr_top_cs(:,:,i_lw_Ts) + lwr_top_cs / nday_year
    rk%flwr_sur_cs(:,:,i_lw_Ts) = rk%flwr_sur_cs(:,:,i_lw_Ts) + sum((flwr_dw_sur_cs-rk_flwr_up_sur)*frst,3) / nday_year

    ! CO2, doubling 
    call lw_radiation(1._wp, frst, zsa, zs, htrop, hcld, ra2, &   ! in
      gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, 2._wp*co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &  ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld) ! out
    rk%flwr_top(:,:,i_lw_co2)    = rk%flwr_top(:,:,i_lw_co2)    + lwr_top  / nday_year
    rk%flwr_sur(:,:,i_lw_co2)    = rk%flwr_sur(:,:,i_lw_co2)    + lwr_sur   / nday_year
    rk%flwr_top_cs(:,:,i_lw_co2) = rk%flwr_top_cs(:,:,i_lw_co2) + lwr_top_cs / nday_year
    rk%flwr_sur_cs(:,:,i_lw_co2) = rk%flwr_sur_cs(:,:,i_lw_co2) + sum((flwr_dw_sur_cs-flwr_up_sur)*frst,3) / nday_year


    !------------------------
    ! shortwave radiation
    !------------------------

    ! control
    call sw_radiation(frst, swr_dw_top, coszm, cld, q2, ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, hcld, hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    rk%fswr_top(:,:,i_sw_control)    = rk%fswr_top(:,:,i_sw_control)    + swr_top   / nday_year
    rk%fswr_sur(:,:,i_sw_control)    = rk%fswr_sur(:,:,i_sw_control)    + sum(fswr_sur*frst,3)   / nday_year
    rk%fswr_top_cs(:,:,i_sw_control) = rk%fswr_top_cs(:,:,i_sw_control) + swr_top_cs / nday_year
    rk%fswr_sur_cs(:,:,i_sw_control) = rk%fswr_sur_cs(:,:,i_sw_control) + sum(fswr_sur_cs*frst,3) / nday_year

    ! water vapor, increase water content corresponding to 1K warming, ~7% from Clausius-Klapeyron
    call sw_radiation(frst, swr_dw_top, coszm, cld, q2*1.07, ra2, &    ! in
      alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, clot, hcld, hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    rk%fswr_top(:,:,i_sw_wv)    = rk%fswr_top(:,:,i_sw_wv)    + swr_top   / nday_year
    rk%fswr_sur(:,:,i_sw_wv)    = rk%fswr_sur(:,:,i_sw_wv)    + sum(fswr_sur*frst,3)   / nday_year
    rk%fswr_top_cs(:,:,i_sw_wv) = rk%fswr_top_cs(:,:,i_sw_wv) + swr_top_cs / nday_year
    rk%fswr_sur_cs(:,:,i_sw_wv) = rk%fswr_sur_cs(:,:,i_sw_wv) + sum(fswr_sur_cs*frst,3) / nday_year

    ! albedo, 1% increase
    call sw_radiation(frst, swr_dw_top, coszm, cld, q2, ra2, &    ! in
      alb_vu_s+0.01, alb_vu_c+0.01, alb_ir_s+0.01, alb_ir_c+0.01, clot, hcld, hqeff, aerosol_ot, aerosol_im, so4, &    ! in
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld) ! out
    rk%fswr_top(:,:,i_sw_alb)    = rk%fswr_top(:,:,i_sw_alb)    + swr_top   / nday_year
    rk%fswr_sur(:,:,i_sw_alb)    = rk%fswr_sur(:,:,i_sw_alb)    + sum(fswr_sur*frst,3)   / nday_year
    rk%fswr_top_cs(:,:,i_sw_alb) = rk%fswr_top_cs(:,:,i_sw_alb) + swr_top_cs / nday_year
    rk%fswr_sur_cs(:,:,i_sw_alb) = rk%fswr_sur_cs(:,:,i_sw_alb) + sum(fswr_sur_cs*frst,3) / nday_year


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
    deallocate(rk_flwr_up_sur)

    return

  end subroutine rad_kernels


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  r a d _ k e r n e l s _ w r i t e
  !   Purpose    :  write radiative kernels
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine rad_kernels_write(rk)

    implicit none

    type(rad_kernels_type), intent(in) :: rk

    integer :: ncid
    character (len=256) :: fnm


    ! write radiative kernels to file

    fnm = trim(out_dir)//"/radiative_kernels.nc"
    call nc_create(fnm)
    call nc_open(fnm,ncid)
    call nc_write_dim(fnm,"lon",x=lon,axis="x",ncid=ncid)
    call nc_write_dim(fnm,"lat",x=lat,axis="y",ncid=ncid)

    ! write 2D

    call nc_write(fnm,"K_2d_toa_sw_wv",    sngl(rk%fswr_top(:,jm:1:-1,i_sw_wv)   -rk%fswr_top(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_toa_cs_sw_wv", sngl(rk%fswr_top_cs(:,jm:1:-1,i_sw_wv)-rk%fswr_top_cs(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA clear-sky shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_sw_wv",    sngl(rk%fswr_sur(:,jm:1:-1,i_sw_wv)   -rk%fswr_sur(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_cs_sw_wv", sngl(rk%fswr_sur_cs(:,jm:1:-1,i_sw_wv)-rk%fswr_sur_cs(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface clear-sky shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_2d_toa_sw_alb",    sngl(rk%fswr_top(:,jm:1:-1,i_sw_alb)   -rk%fswr_top(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA shortwave albedo radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_toa_cs_sw_alb", sngl(rk%fswr_top_cs(:,jm:1:-1,i_sw_alb)-rk%fswr_top_cs(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA clear-sky shortwave albedo radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_sw_alb",    sngl(rk%fswr_sur(:,jm:1:-1,i_sw_alb)   -rk%fswr_sur(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface shortwave albedo radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_cs_sw_alb", sngl(rk%fswr_sur_cs(:,jm:1:-1,i_sw_alb)-rk%fswr_sur_cs(:,jm:1:-1,i_sw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface clear-sky shortwave albedo radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_2d_toa_lw_wv",    sngl(rk%flwr_top(:,jm:1:-1,i_lw_wv)   -rk%flwr_top(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA longwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_toa_cs_lw_wv", sngl(rk%flwr_top_cs(:,jm:1:-1,i_lw_wv)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA longwave clear-sky water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_lw_wv",    sngl(rk%flwr_sur(:,jm:1:-1,i_lw_wv)   -rk%flwr_sur(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface longwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_cs_lw_wv", sngl(rk%flwr_sur_cs(:,jm:1:-1,i_lw_wv)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface clear sky longwave water vapor radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_2d_toa_lw_Ta",    sngl(rk%flwr_top(:,jm:1:-1,i_lw_Ta)   -rk%flwr_top(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA longwave air temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_toa_cs_lw_Ta", sngl(rk%flwr_top_cs(:,jm:1:-1,i_lw_Ta)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA clear-sky longwave air temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_lw_Ta",    sngl(rk%flwr_sur(:,jm:1:-1,i_lw_Ta)   -rk%flwr_sur(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface longwave air temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_cs_lw_Ta", sngl(rk%flwr_sur_cs(:,jm:1:-1,i_lw_Ta)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface clear-sky longwave air temperature radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_2d_toa_lw_Ts",    sngl(rk%flwr_top(:,jm:1:-1,i_lw_Ts)   -rk%flwr_top(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_toa_cs_lw_Ts", sngl(rk%flwr_top_cs(:,jm:1:-1,i_lw_Ts)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA clear-sky longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_lw_Ts",    sngl(rk%flwr_sur(:,jm:1:-1,i_lw_Ts)   -rk%flwr_sur(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_cs_lw_Ts", sngl(rk%flwr_sur_cs(:,jm:1:-1,i_lw_Ts)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface clear-sky longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_2d_toa_lw_co2",    sngl(rk%flwr_top(:,jm:1:-1,i_lw_co2)   -rk%flwr_top(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA longwave CO2 radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_toa_cs_lw_co2", sngl(rk%flwr_top_cs(:,jm:1:-1,i_lw_co2)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="TOA clear-sky longwave CO2 radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_lw_co2",    sngl(rk%flwr_sur(:,jm:1:-1,i_lw_co2)   -rk%flwr_sur(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface longwave CO2 radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_2d_sur_cs_lw_co2", sngl(rk%flwr_sur_cs(:,jm:1:-1,i_lw_co2)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control)), &
      dims=["lon","lat"],start=[1,1],count=[im,jm], long_name="surface clear-sky longwave CO2 radiative kernel",units="W/m2",ncid=ncid)


    ! write zonal mean

    call nc_write(fnm,"K_zon_toa_sw_wv",    sngl(sum(rk%fswr_top(:,jm:1:-1,i_sw_wv)   -rk%fswr_top(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_toa_cs_sw_wv", sngl(sum(rk%fswr_top_cs(:,jm:1:-1,i_sw_wv)-rk%fswr_top_cs(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA clear-sky shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_sw_wv",    sngl(sum(rk%fswr_sur(:,jm:1:-1,i_sw_wv)   -rk%fswr_sur(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_cs_sw_wv", sngl(sum(rk%fswr_sur_cs(:,jm:1:-1,i_sw_wv)-rk%fswr_sur_cs(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface clear-sky shortwave water vapor radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_zon_toa_sw_alb",    sngl(sum(rk%fswr_top(:,jm:1:-1,i_sw_alb)   -rk%fswr_top(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA shortwave albedo radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_toa_cs_sw_alb", sngl(sum(rk%fswr_top_cs(:,jm:1:-1,i_sw_alb)-rk%fswr_top_cs(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA clear-sky shortwave albedo radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_sw_alb",    sngl(sum(rk%fswr_sur(:,jm:1:-1,i_sw_alb)   -rk%fswr_sur(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface shortwave albedo radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_cs_sw_alb", sngl(sum(rk%fswr_sur_cs(:,jm:1:-1,i_sw_alb)-rk%fswr_sur_cs(:,jm:1:-1,i_sw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface clear-sky shortwave albedo radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_zon_toa_lw_wv",    sngl(sum(rk%flwr_top(:,jm:1:-1,i_lw_wv)   -rk%flwr_top(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA longwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_toa_cs_lw_wv", sngl(sum(rk%flwr_top_cs(:,jm:1:-1,i_lw_wv)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA longwave clear-sky water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_lw_wv",    sngl(sum(rk%flwr_sur(:,jm:1:-1,i_lw_wv)   -rk%flwr_sur(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface longwave water vapor radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_cs_lw_wv", sngl(sum(rk%flwr_sur_cs(:,jm:1:-1,i_lw_wv)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface clear sky longwave water vapor radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_zon_toa_lw_Ta",    sngl(sum(rk%flwr_top(:,jm:1:-1,i_lw_Ta)   -rk%flwr_top(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA longwave air temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_toa_cs_lw_Ta", sngl(sum(rk%flwr_top_cs(:,jm:1:-1,i_lw_Ta)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA clear-sky longwave air temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_lw_Ta",    sngl(sum(rk%flwr_sur(:,jm:1:-1,i_lw_Ta)   -rk%flwr_sur(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface longwave air temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_cs_lw_Ta", sngl(sum(rk%flwr_sur_cs(:,jm:1:-1,i_lw_Ta)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface clear-sky longwave air temperature radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_zon_toa_lw_Ts",    sngl(sum(rk%flwr_top(:,jm:1:-1,i_lw_Ts)   -rk%flwr_top(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_toa_cs_lw_Ts", sngl(sum(rk%flwr_top_cs(:,jm:1:-1,i_lw_Ts)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA clear-sky longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_lw_Ts",    sngl(sum(rk%flwr_sur(:,jm:1:-1,i_lw_Ts)   -rk%flwr_sur(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_cs_lw_Ts", sngl(sum(rk%flwr_sur_cs(:,jm:1:-1,i_lw_Ts)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface clear-sky longwave surface temperature radiative kernel",units="W/m2",ncid=ncid)

    call nc_write(fnm,"K_zon_toa_lw_co2",    sngl(sum(rk%flwr_top(:,jm:1:-1,i_lw_co2)   -rk%flwr_top(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA longwave CO2 radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_toa_cs_lw_co2", sngl(sum(rk%flwr_top_cs(:,jm:1:-1,i_lw_co2)-rk%flwr_top_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="TOA clear-sky longwave CO2 radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_lw_co2",    sngl(sum(rk%flwr_sur(:,jm:1:-1,i_lw_co2)   -rk%flwr_sur(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface longwave CO2 radiative kernel",units="W/m2",ncid=ncid)
    call nc_write(fnm,"K_zon_sur_cs_lw_co2", sngl(sum(rk%flwr_sur_cs(:,jm:1:-1,i_lw_co2)-rk%flwr_sur_cs(:,jm:1:-1,i_lw_control),1)/im), &
      dims=["lat"],start=[1],count=[jm], long_name="surface clear-sky longwave CO2 radiative kernel",units="W/m2",ncid=ncid)

    call nc_close(ncid)

    return

  end subroutine rad_kernels_write

end module rad_kernels_mod


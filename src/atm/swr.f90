!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : s w _ r a d i a t i o n _ m o d
!
!  Purpose : shortwave radiation
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
module sw_radiation_mod

  use atm_params, only : wp
  use constants, only : T0, frac_vu
  use timer, only : doy
  use atm_params, only : eps, f_ice_pow, hatm, r_scat, a1_w, a2_w, b1_w, b2_w, c_itf_c, c_itf_cc
  use atm_params, only : l_so4_de, beta_so4, sigma_so4
  use atm_params, only : i_swr_ray, tau_ray, i_swr_gas, a_gas_vu, b_gas_vu, a_gas_ir, b_gas_ir
  use atm_params, only : hpbl, c_hcld_dz_min, n_wp, p_cosz_sw
  use atm_grid, only : im ,jm, nm, i_ice
  !$ use omp_lib

  implicit none

  real(wp), parameter :: p_1 = -1.97_wp 
  real(wp), parameter :: p_2 = 0.82_wp 
  real(wp), parameter :: p_3 = 0.35_wp
  real(wp), parameter :: p_4 = 0.67_wp         
  real(wp), parameter :: alf_1 = 7.73e-2_wp 
  real(wp), parameter :: alf_2 = 2.39e-2_wp 
  real(wp), parameter :: alf_3 = 1.51e2_wp 
  real(wp), parameter :: gam_ar_1 = 2.75_wp 
  real(wp), parameter :: gam_ar_2 = 0.636_wp 
  real(wp), parameter :: gl_c = 0.14_wp 
  real(wp), parameter :: c_itf_o = 0.98_wp 

  private
  public :: sw_radiation, sw_cosz_eff

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s w _ c o s z _ e f f
  !   Purpose    :  effective cos z of the single daily SW column: the insolation-weighted power mean
  !                 (sum_h cosz**(1+p) / sum_h cosz)**(1/p), p = p_cosz_sw, over the hourly cos z of the day.
  !                 p = -1 is the daytime mean (exact for absorption linear in the path, 1/cosz), p = 1 the
  !                 insolation-weighted mean (exact for terms linear in cosz); gas/water absorption (~ path**0.3)
  !                 and the cloud albedo (~ cosz**-p_3) lie in between
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sw_cosz_eff(cosz_h, coszm, cosz_eff)

    implicit none

    real(wp), intent(in ) :: cosz_h(:,:)   ! (hour, lat) hourly cos z of the day
    real(wp), intent(in ) :: coszm(:)      ! (lat) daytime-mean cos z, used where no hourly sample is sunlit
    real(wp), intent(out) :: cosz_eff(:)   ! (lat)

    integer :: j
    real(wp) :: s1

    do j=1,size(cosz_eff)
      s1 = sum(cosz_h(:,j), mask=cosz_h(:,j).gt.0._wp)
      if (s1.gt.0._wp) then
        if (abs(p_cosz_sw).lt.1.e-6_wp) then
          cosz_eff(j) = exp(sum(cosz_h(:,j)*log(max(cosz_h(:,j),eps)), mask=cosz_h(:,j).gt.0._wp)/s1)
        else
          cosz_eff(j) = (sum(cosz_h(:,j)**(1._wp+p_cosz_sw), mask=cosz_h(:,j).gt.0._wp)/s1)**(1._wp/p_cosz_sw)
        endif
      else
        cosz_eff(j) = coszm(j)
      endif
    enddo

    return

  end subroutine sw_cosz_eff


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s w _ r a d i a t i o n 
  !   Purpose    :  driver for short-wave radiation at the top and bottom of the atmosphere 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sw_radiation(frst, swr_dw_top, coszm, cld, wcon, zsa, alb_vu_s, alb_vu_c, alb_ir_s, alb_ir_c, &
      clot, hcld, hqeff, aerosol_ot, aerosol_im, so4, &
      alb_cld, swr_top, swr_top_cs, swr_top_cld, swr_sur, fswr_sur, fswr_sur_cs, fswr_sur_cld, &
      dswd_dalb_vu_cs, dswd_dalb_ir_cs, dswd_dalb_vu_cld, dswd_dalb_ir_cld, dswd_dz_ir_cs, dswd_dz_ir_cld, &
      swr_dw_sur_vis_cs, swr_dw_sur_nir_cs, swr_dw_sur_vis_cld, swr_dw_sur_nir_cld, zs)

    implicit none

    real(wp), intent(in ) :: frst(:,:,:)
    real(wp), intent(in ) :: swr_dw_top(:,:)
    real(wp), intent(in ) :: coszm(:,:)
    real(wp), intent(in ) :: cld(:,:)
    real(wp), intent(in ) :: wcon(:,:)   ! column water vapour (kg/m2), at the grid-cell mean elevation zsa
    real(wp), intent(in ) :: zsa(:,:)    ! grid-cell mean surface elevation (m)
    real(wp), intent(in ) :: alb_vu_s(:,:,:)
    real(wp), intent(in ) :: alb_vu_c(:,:,:)
    real(wp), intent(in ) :: alb_ir_s(:,:,:)
    real(wp), intent(in ) :: alb_ir_c(:,:,:)
    real(wp), intent(in ) :: clot(:,:)
    real(wp), intent(in ) :: hcld(:,:)
    real(wp), intent(in ) :: hqeff(:,:)
    real(wp), intent(in ) :: aerosol_ot(:,:)
    real(wp), intent(in ) :: aerosol_im(:,:)
    real(wp), intent(in ) :: so4(:,:)

    real(wp), intent(out) :: alb_cld(:,:)
    real(wp), intent(out) :: swr_top(:,:)
    real(wp), intent(out) :: swr_top_cs(:,:)
    real(wp), intent(out) :: swr_top_cld(:,:)
    real(wp), intent(out) :: swr_sur(:,:)
    real(wp), intent(out) :: fswr_sur(:,:,:)
    real(wp), intent(out) :: fswr_sur_cs(:,:,:)
    real(wp), intent(out) :: fswr_sur_cld(:,:,:)

    real(wp), intent(out), optional :: dswd_dalb_vu_cs(:,:)
    real(wp), intent(out), optional :: dswd_dalb_ir_cs(:,:)
    real(wp), intent(out), optional :: dswd_dalb_vu_cld(:,:)
    real(wp), intent(out), optional :: dswd_dalb_ir_cld(:,:)
    real(wp), intent(out), optional :: dswd_dz_ir_cs(:,:)
    real(wp), intent(out), optional :: dswd_dz_ir_cld(:,:)
    real(wp), intent(out), optional :: swr_dw_sur_vis_cs(:,:)  
    real(wp), intent(out), optional :: swr_dw_sur_nir_cs(:,:)  
    real(wp), intent(out), optional :: swr_dw_sur_vis_cld(:,:) 
    real(wp), intent(out), optional :: swr_dw_sur_nir_cld(:,:) 
    real(wp), intent(in ) :: zs(:,:,:)   ! surface elevation of the surface types (m); sets the cloud geometry and
                                         ! scales the well-mixed gas path (i_swr_gas==1)

    integer :: i, j, n
    real(wp) :: fps, zsn, h_cld_top, h_cld_gt, w_col
    logical :: l_dswd_dalb
    real(wp) :: solar_top_up, solar_top_up_s, solar_top_up_c, solar_sur, solar_sur_s, solar_sur_c

    real(wp), allocatable, dimension(:) :: fst
    real(wp), allocatable, dimension(:) :: fswr_top
    real(wp), allocatable, dimension(:) :: fswr_top_cs
    real(wp), allocatable, dimension(:) :: fswr_top_cld


    if (present(dswd_dalb_vu_cs)) then
      l_dswd_dalb = .true.
    else
      l_dswd_dalb = .false.
    endif

    allocate(fst(nm))
    allocate(fswr_top(nm))
    allocate(fswr_top_cs(nm))
    allocate(fswr_top_cld(nm))

    !$omp parallel do collapse(2) private(i,j,n,solar_top_up,solar_sur,solar_top_up_s,solar_top_up_c,solar_sur_s,solar_sur_c) &
    !$omp private(fst,fswr_top,fswr_top_cs,fswr_top_cld,fps,zsn,h_cld_top,h_cld_gt,w_col)
    do j=1,jm
      do i=1,im

        ! option to increase ice-albedo feedback by increasing 'effective' ice
        ! fraction in grid cell, only for shortwave radiation
        fst(i_ice) = frst(i,j,i_ice)**f_ice_pow
        do n=1,nm
          if (n.ne.i_ice) then
            if (fst(i_ice).lt.1._wp) then
              fst(n) = frst(i,j,n)*(1._wp-fst(i_ice))/(1._wp-frst(i,j,i_ice))
            else
              fst(n) = 0._wp
            endif
          endif
        enddo

        do n=1,nm
          if (fst(n).gt.0._wp .and. swr_dw_top(i,j).gt.eps) then 

            zsn = zs(i,j,n)
            fps = exp(-zsn/hatm)
            ! cloud geometry: the cloud reflects at its top, hcld (above sea level, as in the longwave
            ! column), and fills the layer below it down to the longwave cloud base, hpbl above the
            ! surface (at least c_hcld_dz_min thick)
            h_cld_top = max(hcld(i,j)-zsn, 0._wp)
            h_cld_gt  = max(hcld(i,j)-zsn-hpbl, c_hcld_dz_min)
            ! water vapour above the surface type: the column wcon (consistent with hqeff = wcon/(ra2a*qam))
            ! moved from the cell-mean elevation to the elevation of the surface type with the scale height hqeff
            w_col = wcon(i,j)*exp(-(zsn-zsa(i,j))/hqeff(i,j))
            call sw_radiation_col(swr_dw_top(i,j), coszm(doy,j), cld(i,j), w_col, &
              alb_vu_s(i,j,n), alb_ir_s(i,j,n), alb_vu_c(i,j,n), alb_ir_c(i,j,n), &
              aerosol_ot(i,j), aerosol_im(i,j), so4(i,j), clot(i,j), h_cld_top, hqeff(i,j), h_cld_gt, &
              alb_cld(i,j), solar_top_up, solar_top_up_s, solar_top_up_c, solar_sur, solar_sur_s, solar_sur_c, &
              fps=fps)

            ! net shortwave radiation at TOA
            fswr_top(n)     = swr_dw_top(i,j)-solar_top_up
            fswr_top_cs(n)  = swr_dw_top(i,j)-solar_top_up_s
            fswr_top_cld(n) = swr_dw_top(i,j)-solar_top_up_c
            ! net shortwave radiation at the surface
            fswr_sur(i,j,n)     = solar_sur
            fswr_sur_cs(i,j,n)  = solar_sur_s
            fswr_sur_cld(i,j,n) = solar_sur_c

          else

            fswr_top(n)         = 0._wp
            fswr_top_cs(n)      = 0._wp
            fswr_top_cld(n)     = 0._wp
            fswr_sur(i,j,n)     = 0._wp
            fswr_sur_cs(i,j,n)  = 0._wp
            fswr_sur_cld(i,j,n) = 0._wp

          endif 

        enddo

        ! grid cell averages
        swr_top(i,j)     = sum(fswr_top(:)*fst(:))
        swr_top_cs(i,j)  = sum(fswr_top_cs(:)*fst(:))
        swr_top_cld(i,j) = sum(fswr_top_cld(:)*fst(:))
        swr_sur(i,j)     = sum(fswr_sur(i,j,:)*fst(:))

        ! compute partial derivatives of downward surface solar radiation with
        ! respect to surface albedo and surface elevation 
        if (l_dswd_dalb) then
          if (swr_dw_top(i,j).gt.eps) then

            zsn = sum(zs(i,j,:)*fst(:))
            fps = exp(-zsn/hatm)
            h_cld_top = max(hcld(i,j)-zsn, 0._wp)
            h_cld_gt  = max(hcld(i,j)-zsn-hpbl, c_hcld_dz_min)
            w_col = wcon(i,j)*exp(-(zsn-zsa(i,j))/hqeff(i,j))
            call sw_radiation_col(swr_dw_top(i,j), coszm(doy,j), cld(i,j), w_col, &
              sum(alb_vu_s(i,j,:)*fst(:)), sum(alb_ir_s(i,j,:)*fst(:)), sum(alb_vu_c(i,j,:)*fst(:)), sum(alb_ir_c(i,j,:)*fst(:)), &
              aerosol_ot(i,j), aerosol_im(i,j), so4(i,j), clot(i,j), h_cld_top, hqeff(i,j), h_cld_gt, &
              alb_cld(i,j), solar_top_up, solar_top_up_s, solar_top_up_c, solar_sur, solar_sur_s, solar_sur_c, &  ! out
              dswd_dalb_vu_cs(i,j), dswd_dalb_ir_cs(i,j), dswd_dalb_vu_cld(i,j), dswd_dalb_ir_cld(i,j), dswd_dz_ir_cs(i,j), dswd_dz_ir_cld(i,j), & ! out
              swr_dw_sur_vis_cs(i,j), swr_dw_sur_nir_cs(i,j), swr_dw_sur_vis_cld(i,j), swr_dw_sur_nir_cld(i,j), &
              fps=fps)

          else

            dswd_dalb_vu_cs(i,j)  = 0._wp
            dswd_dalb_ir_cs(i,j)  = 0._wp
            dswd_dalb_vu_cld(i,j) = 0._Wp
            dswd_dalb_ir_cld(i,j) = 0._wp
            dswd_dz_ir_cs(i,j)    = 0._wp
            dswd_dz_ir_cld(i,j)   = 0._wp
            swr_dw_sur_vis_cs(i,j)  = 0._wp
            swr_dw_sur_nir_cs(i,j)  = 0._wp
            swr_dw_sur_vis_cld(i,j) = 0._wp
            swr_dw_sur_nir_cld(i,j) = 0._wp

          endif
        endif

      enddo
    enddo 
    !$omp end parallel do

    deallocate(fst)
    deallocate(fswr_top)
    deallocate(fswr_top_cs)
    deallocate(fswr_top_cld)

    return

  end subroutine sw_radiation


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s w _ r a d i a t i o n _ c o l
  !   Purpose    :  calculation of SW radiation at the top and bottom
  !                 of the atmosphere in a single column using two-stream
  !                 approximation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sw_radiation_col(solar_top, coszen, cld, w_col, &
      alb_sur_vu_s, alb_sur_ir_s, alb_sur_vu_c, alb_sur_ir_c, &
      aerosol_ot, aerosol_im, so4, cld_ot, h_c, h_q, cld_gt, &
      alb_cld, solar_top_up, solar_top_up_s, solar_top_up_c, solar_sur, solar_sur_s, solar_sur_c, &
      dswd_dalb_vu_cs, dswd_dalb_ir_cs, dswd_dalb_vu_cld, dswd_dalb_ir_cld, dswd_dz_ir_cs, dswd_dz_ir_cld, &
      swr_dw_sur_vis_cs, swr_dw_sur_nir_cs, swr_dw_sur_vis_cld, swr_dw_sur_nir_cld, fps)

    implicit none

    real(wp), intent(in ) :: solar_top
    real(wp), intent(in ) :: coszen
    real(wp), intent(in ) :: cld
    real(wp), intent(in ) :: w_col    ! column water vapour above the surface (kg/m2)
    real(wp), intent(in ) :: alb_sur_vu_s
    real(wp), intent(in ) :: alb_sur_ir_s
    real(wp), intent(in ) :: alb_sur_vu_c
    real(wp), intent(in ) :: alb_sur_ir_c
    real(wp), intent(in ) :: aerosol_ot
    real(wp), intent(in ) :: aerosol_im
    real(wp), intent(in ) :: so4 
    real(wp), intent(in ) :: cld_ot
    real(wp), intent(in ) :: h_c      ! height of the cloud top above the surface (m); the cloud occupies [h_c-cld_gt, h_c]
    real(wp), intent(in ) :: h_q      ! water vapour scale height (m)
    real(wp), intent(in ) :: cld_gt   ! geometrical thickness of the cloud (m)

    real(wp), intent(out) :: alb_cld
    real(wp), intent(out) :: solar_top_up
    real(wp), intent(out) :: solar_top_up_s
    real(wp), intent(out) :: solar_top_up_c
    real(wp), intent(out) :: solar_sur
    real(wp), intent(out) :: solar_sur_s
    real(wp), intent(out) :: solar_sur_c

    real(wp), intent(out), optional :: dswd_dalb_vu_cs
    real(wp), intent(out), optional :: dswd_dalb_ir_cs
    real(wp), intent(out), optional :: dswd_dalb_vu_cld
    real(wp), intent(out), optional :: dswd_dalb_ir_cld
    real(wp), intent(out), optional :: dswd_dz_ir_cs
    real(wp), intent(out), optional :: dswd_dz_ir_cld
    real(wp), intent(out), optional :: swr_dw_sur_vis_cs  
    real(wp), intent(out), optional :: swr_dw_sur_nir_cs  
    real(wp), intent(out), optional :: swr_dw_sur_vis_cld 
    real(wp), intent(out), optional :: swr_dw_sur_nir_cld 
    real(wp), intent(in ), optional :: fps   ! surface pressure / p0, scales the well-mixed gas path (i_swr_gas==1)

    ! local variables
    logical :: l_dswd_dalb
    real(wp) :: exp_hc_hq, exp_hb_hq
    real(wp) :: fp
    real(wp) :: cos_zen, cos_zen_o
    real(wp) :: f_1, f_2, f_3
    real(wp) :: b_ar, icos, b_c
    real(wp) :: alb_a_s, alb_a_c
    real(wp) :: alb_sa_vu_s, alb_sa_ir_s, alb_sa_vu_c, alb_sa_ir_c
    real(wp) :: alb_sct_ir, alb_sct_vu, alb_cld_ir, alb_cld_vu
    real(wp) :: alb_sct_ir_0, alb_sct_vu_0
    real(wp) :: alb_sa_ir_0, alb_sa_vu_0
    real(wp) :: alb_atm_ir_s, alb_atm_vu_s, alb_atm_ir_c, alb_atm_vu_c
    real(wp) :: b_arb, b_arb_d1, b_arb_d2
    real(wp) :: rqh
    real(wp) :: itf_w_vu_s, itf_w_ir_s, itf_a_vu_s, itf_a_ir_s
    real(wp) :: itf_o_vu_s, itf_o_ir_s
    real(wp) :: itf_w_vu_c, itf_w_ir_c, itf_a_vu_c, itf_a_ir_c
    real(wp) :: itf_o_vu_c, itf_o_ir_c, itf_c_vu, itf_c_ir
    real(wp) :: itf_c_vu_d1, itf_c_ir_d1, itf_c_vu_d2, itf_c_ir_d2
    real(wp) :: itf_o_vu_d1, itf_o_vu_d2, itf_o_ir_d1, itf_o_ir_d2
    real(wp) :: itf_w_vu_s_d1, itf_w_ir_s_d1, itf_w_vu_s_d2, itf_w_ir_s_d2
    real(wp) :: itf_w_vu_c_d1, itf_w_ir_c_d1, itf_w_vu_c_d2, itf_w_ir_c_d2
    real(wp) :: itf_a_vu_d1, itf_a_ir_d1, itf_a_vu_d2, itf_a_ir_d2
    real(wp) :: itf_atm_vu_s, itf_atm_ir_s, itf_atm_vu_c, itf_atm_ir_c
    real(wp) :: m_w_s_d1, m_w_s_d2, m_w_c_d1, m_w_c_d2, exp_1, exp_2
    real(wp) :: m_w_c, m_w_s
    real(wp) :: drqh_dz
    real(wp) :: dm_w_s_d1_dz, dm_w_s_d2_dz, ditf_w_ir_s_d1_dz, ditf_w_ir_s_d2_dz
    real(wp) :: dm_w_c_d1_dz, dm_w_c_d2_dz, ditf_w_ir_c_d1_dz, ditf_w_ir_c_d2_dz

    !***********************************************************
    !*     itf_    integral transmission function
    !*        _atm_ entire atmosphere
    !*        _w_   water
    !*        _a_   aerosol
    !*        -o_   ozone
    !*        _c_   clouds
    !*          _vu_  visible+UV
    !*          _ir_  infrared 
    !*            _s   clear sky
    !*            -c   cloudy conditions
    !*
    !*     alb_    albedo
    !*        _sur_ surface
    !*        _sct_ scattering
    !*        _cld_ clouds
    !*
    !*    frac_vu  fraction of VU radiation in total radiation
    !*    cld_ot  optical thickness of clouds
    !*    cld_gt   geometrical thickness of clouds
    !*    aerosol_ot  optical thickness of erosol
    !*    aerosol_im  imaginary part of aerosol refractive index
    !*   
    !************************************************************************

    if (present(dswd_dalb_vu_cs)) then
      l_dswd_dalb = .true.
    else
      l_dswd_dalb = .false.
    endif

    ! fractions of the column water vapour above the cloud top (exp_hc_hq) and above the cloud base (exp_hb_hq)
    exp_hc_hq = exp(-h_c/h_q)
    exp_hb_hq = exp(-max(h_c-cld_gt,0._wp)/h_q)

    fp = 1._wp
    if (present(fps)) fp = fps

    cos_zen = max(coszen,0.1_wp)
    cos_zen_o = 1._wp/1.66_wp
    icos = 1._wp/cos_zen+1._wp/cos_zen_o

    if (l_so4_de) then
      ! direct sulfate aerosol effect

      ! sulfate aerosol albedo (eq. 6 in Bauer et al., 2008)
      ! clear sky
      alb_a_s = beta_so4*sigma_so4*so4/cos_zen
      ! cloudy sky
      alb_a_c = beta_so4*sigma_so4*so4/cos_zen_o

      ! combined surface-aerosol albedo following Bauer et al., 2008 (eq. 5)
      alb_sa_vu_s = alb_a_s + (1._wp-alb_a_s)**2*alb_sur_vu_s/(1._wp-alb_a_s*alb_sur_vu_s)
      alb_sa_ir_s = alb_a_s + (1._wp-alb_a_s)**2*alb_sur_ir_s/(1._wp-alb_a_s*alb_sur_ir_s)
      alb_sa_vu_c = alb_a_c + (1._wp-alb_a_c)**2*alb_sur_vu_c/(1._wp-alb_a_c*alb_sur_vu_c)
      alb_sa_ir_c = alb_a_c + (1._wp-alb_a_c)**2*alb_sur_ir_c/(1._wp-alb_a_c*alb_sur_ir_c)

    else

      alb_sa_vu_s = alb_sur_vu_s
      alb_sa_ir_s = alb_sur_ir_s
      alb_sa_vu_c = alb_sur_vu_c
      alb_sa_ir_c = alb_sur_ir_c

    endif

    ! albedo of the surface (with sulfate aerosol) for the diffuse, surface-reflected light that the clear
    ! atmosphere scatters back down: the diffuse surface albedo
    alb_sa_vu_0 = alb_sa_vu_c
    alb_sa_ir_0 = alb_sa_ir_c

    !-----------------------------------------------
    ! 1) SW flux at the top of atmosphere   
    !-----------------------------------------------

    !-----------------------------------------------
    ! 1.1) clear sky conditions

    ! scatteres atmospheric albedo  VU and IR

    b_ar = 0.55_wp*aerosol_ot 
    b_arb = b_ar*icos
    f_1 = cos_zen**p_1
    f_2 = b_ar**p_2
    f_3 = alf_1-alf_2*log(1._wp+alf_3*aerosol_im)

    alb_sct_vu = 1._wp-(1._wp-r_scat_mu(cos_zen))*exp(-f_1*f_2*f_3)
    alb_sct_ir = 1._wp-exp(-f_1*f_2*f_3)

    ! albedo of the clear atmosphere for the diffuse light reflected by the surface (from below, at the diffuse angle)
    f_1 = cos_zen_o**p_1
    alb_sct_vu_0 = 1._wp-(1._wp-r_scat_mu(cos_zen_o))*exp(-f_1*f_2*f_3)
    alb_sct_ir_0 = 1._wp-exp(-f_1*f_2*f_3)

    ! integral transmission function (ITF) for water

    rqh = 0.1_wp*w_col    ! column water content in g/cm2
    ! effective water path: x (ps/p0)**n_wp; at fixed column water RRTMG absorbs more SW in columns with a lower
    ! surface pressure (n_wp < 0, fitted on elevated Antarctic/Greenland/Tibet columns)
    rqh = rqh*fp**n_wp
    m_w_s = rqh*icos

    itf_w_ir_s = a1_w*exp(-b1_w*m_w_s)+a2_w*exp(-b2_w*m_w_s)
    itf_w_vu_s = 1._wp

    ! integral transmission function (ITF) for aerosol

    itf_a_vu_s = exp(-gam_ar_1*b_arb*aerosol_im**gam_ar_2)
    itf_a_ir_s = itf_a_vu_s 

    ! integral transmission function (ITF) for O3 (and other non-water gases if i_swr_gas==1)

    itf_o_vu_s = itf_gas_vu(icos)
    itf_o_ir_s = itf_gas_ir(fp*(icos))

    ! planetary albedo: the (diffuse) light reflected by the surface escapes with 1-alb_sct_*_0 or is scattered
    ! back down with alb_sct_*_0 onto the diffuse surface albedo (as in the surface flux below)

    alb_atm_vu_s = &
      (alb_sct_vu+(1._wp-alb_sct_vu)*(1._wp-alb_sct_vu_0)*alb_sa_vu_s/ &
      (1._wp-alb_sct_vu_0*alb_sa_vu_0))* &
      itf_w_vu_s*itf_a_vu_s*itf_o_vu_s

    alb_atm_ir_s = &
      (alb_sct_ir+(1._wp-alb_sct_ir)*(1._wp-alb_sct_ir_0)*alb_sa_ir_s/ &
      (1._wp-alb_sct_ir_0*alb_sa_ir_0))* &
      itf_w_ir_s*itf_a_ir_s*itf_o_ir_s

    ! SW flux at the top of atmosphere, clear sky

    solar_top_up_s = solar_top*(frac_vu*alb_atm_vu_s+(1.-frac_vu)*alb_atm_ir_s)

    !-----------------------------------------------
    ! 1.2) cloudy conditions

    ! integral transmission function (ITF) for clouds

    itf_c_ir = c_itf_c  
    itf_c_vu = c_itf_c 

    ! integral transmission function (ITF) for water

    ! slant path above the cloud top, plus one vertical path through the cloud layer
    m_w_c = rqh*(exp_hc_hq*icos+(exp_hb_hq-exp_hc_hq))

    itf_w_ir_c = a1_w*exp(-b1_w*m_w_c)+a2_w*exp(-b2_w*m_w_c)
    itf_w_vu_c = 1._wp

    ! integral transmission function (ITF) for aerosol

    b_arb = b_ar*(exp_hc_hq*icos+(exp_hb_hq-exp_hc_hq))

    itf_a_vu_c = exp(-gam_ar_1*b_arb*aerosol_im**gam_ar_2)
    itf_a_ir_c = itf_a_vu_c 

    ! integral transmission function (ITF) for O3 (and other non-water gases if i_swr_gas==1);
    ! the well-mixed IR absorbers only above the cloud top

    itf_o_vu_c = itf_gas_vu(icos)
    itf_o_ir_c = itf_gas_ir(fp*(icos*exp(-h_c/hatm)))

    ! cloud albedo   

    b_c = gl_c/cos_zen**p_3 
    alb_cld_vu = 1._wp-(1._wp-alb_sct_vu)*exp(-b_c*cld_ot**p_4)
    alb_cld_ir = 1._wp-(1._wp-alb_sct_ir)*exp(-b_c*cld_ot**p_4)
    alb_cld = alb_cld_vu

    ! planetary albedo

    alb_atm_ir_c = &
      (alb_cld_ir+((1._wp-alb_cld_ir)**2*alb_sa_ir_c)/ &
      (1._wp-alb_cld_ir*alb_sa_ir_s))* &
      itf_w_ir_c*itf_c_ir*itf_a_ir_c*itf_o_ir_c

    alb_atm_vu_c = &
      (alb_cld_vu+((1._wp-alb_cld_vu)**2*alb_sa_vu_c)/ &
      (1._wp-alb_cld_vu*alb_sa_vu_s))* &
      itf_w_vu_c*itf_c_vu*itf_a_vu_c*itf_o_vu_c

    ! SW flux at the top of atmosphere cloudy sky 

    solar_top_up_c = solar_top*(frac_vu*alb_atm_vu_c+(1._wp-frac_vu)*alb_atm_ir_c)


    ! SW flux at the top of atmosphere 

    solar_top_up = (1._wp-cld)*solar_top_up_s+cld*solar_top_up_c


    !-----------------------------------------------
    ! 2)    SW flux at the surface   
    !-----------------------------------------------

    !-----------------------------------------------
    ! 2.1) clear sky conditions

    ! integral transmission function (ITF) for O3

    itf_o_vu_d1 = itf_gas_vu(1._wp/cos_zen)
    itf_o_ir_d1 = itf_gas_ir(fp*(1._wp/cos_zen))
    itf_o_vu_d2 = itf_gas_vu(1._wp/cos_zen+(1._wp-0.7788008_wp)*2._wp/cos_zen_o)
    itf_o_ir_d2 = itf_gas_ir(fp*(1._wp/cos_zen+(1._wp-0.7788008_wp)*2._wp/cos_zen_o))

    ! integral transmission function (ITF) for water

    m_w_s_d1 = rqh/cos_zen
    m_w_s_d2 = m_w_s_d1+rqh*(1._wp-0.7788008_wp)*2._wp/cos_zen_o    ! 0.7788008 = exp(-0.25)

    itf_w_vu_s_d1 = 1._wp
    itf_w_vu_s_d2 = 1._wp
    itf_w_ir_s_d1 = a1_w*exp(-b1_w*m_w_s_d1)+a2_w*exp(-b2_w*m_w_s_d1)
    itf_w_ir_s_d2 = a1_w*exp(-b1_w*m_w_s_d2)+a2_w*exp(-b2_w*m_w_s_d2)

    ! integral transmission function (ITF) for aerosol

    b_arb_d1 = b_ar/cos_zen
    b_arb_d2 = b_arb_d1+b_ar*(1._wp-0.7788008_wp)*2._wp/cos_zen_o     ! 0.7788008 = exp(-0.25)

    itf_a_vu_d1 = exp(-gam_ar_1*b_arb_d1*aerosol_im**gam_ar_2)
    itf_a_ir_d1 = itf_a_vu_d1 

    itf_a_vu_d2 = exp(-gam_ar_1*b_arb_d2*aerosol_im**gam_ar_2)
    itf_a_ir_d2 = itf_a_vu_d2 

    ! atmospheric ITF

    itf_atm_vu_s = (1._wp-alb_sct_vu)*(1._wp-alb_sa_vu_s)* &
      itf_w_vu_s_d1*itf_a_vu_d1*itf_o_vu_d1+ &
      (1._wp-alb_sct_vu)*alb_sa_vu_s*alb_sct_vu_0*(1._wp-alb_sa_vu_0)/ &
      (1._wp-alb_sct_vu_0*alb_sa_vu_0)* &
      itf_w_vu_s_d2*itf_a_vu_d2*itf_o_vu_d2

    itf_atm_ir_s = (1._wp-alb_sct_ir)*(1._wp-alb_sa_ir_s)* &
      itf_w_ir_s_d1*itf_a_ir_d1*itf_o_ir_d1+ &
      (1._wp-alb_sct_ir)*alb_sa_ir_s*alb_sct_ir_0*(1._wp-alb_sa_ir_0)/ &
      (1._wp-alb_sct_ir_0*alb_sa_ir_0)* &
      itf_w_ir_s_d2*itf_a_ir_d2*itf_o_ir_d2

    solar_sur_s = solar_top*(frac_vu*itf_atm_vu_s+(1._wp-frac_vu)*itf_atm_ir_s)

    if (l_dswd_dalb) then

      ! downward flux: transmitted beam plus the surface-reflected light scattered back down
      swr_dw_sur_vis_cs = solar_top*(1._wp-alb_sct_vu)*(itf_w_vu_s_d1*itf_a_vu_d1*itf_o_vu_d1 &
        + alb_sa_vu_s*alb_sct_vu_0/(1._wp-alb_sct_vu_0*alb_sa_vu_0)*itf_w_vu_s_d2*itf_a_vu_d2*itf_o_vu_d2)
      swr_dw_sur_nir_cs = solar_top*(1._wp-alb_sct_ir)*(itf_w_ir_s_d1*itf_a_ir_d1*itf_o_ir_d1 &
        + alb_sa_ir_s*alb_sct_ir_0/(1._wp-alb_sct_ir_0*alb_sa_ir_0)*itf_w_ir_s_d2*itf_a_ir_d2*itf_o_ir_d2)

      ! derivative of the downward flux with respect to the surface albedo (direct and diffuse changed alike)
      dswd_dalb_vu_cs = solar_top*(1._wp-alb_sct_vu)*alb_sct_vu_0 &
        * (1._wp/(1._wp-alb_sct_vu_0*alb_sa_vu_0) + alb_sa_vu_s*alb_sct_vu_0/(1._wp-alb_sct_vu_0*alb_sa_vu_0)**2) &
        * itf_w_vu_s_d2*itf_a_vu_d2*itf_o_vu_d2

      dswd_dalb_ir_cs = solar_top*(1._wp-alb_sct_ir)*alb_sct_ir_0 &
        * (1._wp/(1._wp-alb_sct_ir_0*alb_sa_ir_0) + alb_sa_ir_s*alb_sct_ir_0/(1._wp-alb_sct_ir_0*alb_sa_ir_0)**2) &
        * itf_w_ir_s_d2*itf_a_ir_d2*itf_o_ir_d2

      ! d/dz (per cm) of the water path: w_col ~ exp(-z/h_q), times fp**n_wp ~ exp(-n_wp*z/hatm)
      drqh_dz = -rqh*(1._wp/h_q + n_wp/hatm)/100._wp
      dm_w_s_d1_dz = 1._wp/cos_zen * drqh_dz 
      dm_w_s_d2_dz = dm_w_s_d1_dz + (1._wp-0.7788008_wp)*2._wp/cos_zen_o * drqh_dz 
      ditf_w_ir_s_d1_dz = a1_w*(-b1_w)*exp(-b1_w*m_w_s_d1)*dm_w_s_d1_dz +a2_w*(-b2_w)*exp(-b2_w*m_w_s_d1)*dm_w_s_d1_dz
      ditf_w_ir_s_d2_dz = a1_w*(-b1_w)*exp(-b1_w*m_w_s_d2)*dm_w_s_d2_dz +a2_w*(-b2_w)*exp(-b2_w*m_w_s_d2)*dm_w_s_d2_dz
      dswd_dz_ir_cs = solar_top * ((1._wp-alb_sct_ir)*itf_a_ir_d1*itf_o_ir_d1 * ditf_w_ir_s_d1_dz + &
        (1._wp-alb_sct_ir)*alb_sa_ir_s*alb_sct_ir_0/(1._wp-alb_sct_ir_0*alb_sa_ir_0)*itf_a_ir_d2*itf_o_ir_d2 * ditf_w_ir_s_d2_dz)
      dswd_dz_ir_cs = dswd_dz_ir_cs*100._wp  ! convert from W/m2/cm to W/m2/m

    endif

    !-----------------------------------------------
    ! 2.2) cloudy conditions

    ! integral transmission function (ITF) for O3

    ! well-mixed IR absorbers: direct beam above the cloud, diffuse below
    itf_o_vu_d1 = itf_gas_vu(1._wp/cos_zen)
    itf_o_ir_d1 = itf_gas_ir(fp*(exp(-h_c/hatm)/cos_zen+(1._wp-exp(-h_c/hatm))/cos_zen_o))
    itf_o_vu_d2 = itf_gas_vu(1._wp/cos_zen+(1._wp-0.7788008_wp)*2._wp/cos_zen_o)
    itf_o_ir_d2 = itf_gas_ir(fp*(exp(-h_c/hatm)/cos_zen+3._wp*(1._wp-exp(-h_c/hatm))/cos_zen_o))

    ! integral transmission function (ITF) for clouds

    itf_c_vu_d1 = c_itf_cc 
    itf_c_vu_d2 = c_itf_cc 
    itf_c_ir_d1 = c_itf_cc 
    itf_c_ir_d2 = c_itf_cc 

    ! integral transmission function (ITF) for water

    ! water in the cloud layer (vertical path)
    exp_1 = exp_hb_hq-exp_hc_hq
    ! diffuse path through the water below the cloud base, (1-exp(-h_b/h_q))*1.66; was 1-exp(-h_c/h_q)*1.66
    ! (operator-precedence typo inherited from CLIMBER-2), which is negative for low clouds
    exp_2 = (1._wp-exp_hb_hq)/cos_zen_o

    m_w_c_d1 = rqh*(exp_hc_hq/cos_zen+exp_1+exp_2)
    m_w_c_d2 = m_w_c_d1+rqh*(2._wp*exp_2+exp_1)

    itf_w_vu_c_d1 = 1._wp
    itf_w_vu_c_d2 = 1._wp
    itf_w_ir_c_d1 = a1_w*exp(-b1_w*m_w_c_d1)+a2_w*exp(-b2_w*m_w_c_d1)
    itf_w_ir_c_d2 = a1_w*exp(-b1_w*m_w_c_d2)+a2_w*exp(-b2_w*m_w_c_d2)

    ! integral transmission function (ITF) for aerosol

    b_arb_d1 = b_ar*(exp_hc_hq/cos_zen+exp_1+exp_2)
    b_arb_d2 = b_arb_d1+b_ar*(exp_1+2._wp*exp_2)

    itf_a_vu_d1 = exp(-gam_ar_1*b_arb_d1*aerosol_im**gam_ar_2)
    itf_a_ir_d1 = itf_a_vu_d1 

    itf_a_vu_d2 = exp(-gam_ar_1*b_arb_d2*aerosol_im**gam_ar_2)
    itf_a_ir_d2 = itf_a_vu_d2 

    ! atmospheric ITF

    itf_atm_vu_c = (1._wp-alb_cld_vu)*(1._wp-alb_sa_vu_c)* &
      itf_c_vu_d1*itf_w_vu_c_d1*itf_a_vu_d1*itf_o_vu_d1+ &
      (1._wp-alb_cld_vu)*alb_sa_vu_c*alb_cld_vu* &
      (1._wp-alb_sa_vu_c)/ &
      (1._wp-alb_cld_vu*alb_sa_vu_c)* &
      itf_c_vu_d2*itf_w_vu_c_d2*itf_a_vu_d2*itf_o_vu_d2

    itf_atm_ir_c = (1._wp-alb_cld_ir)*(1._wp-alb_sa_ir_c)* &
      itf_c_ir_d1*itf_w_ir_c_d1*itf_a_ir_d1*itf_o_ir_d1+ &
      (1._wp-alb_cld_ir)*alb_sa_ir_c*alb_cld_ir* &
      (1._wp-alb_sa_ir_c)/ &
      (1._wp-alb_cld_ir*alb_sa_ir_c)* &
      itf_c_ir_d2*itf_w_ir_c_d2*itf_a_ir_d2*itf_o_ir_d2

    solar_sur_c = solar_top*(frac_vu*itf_atm_vu_c+(1._wp-frac_vu)*itf_atm_ir_c)

    if (l_dswd_dalb) then

      swr_dw_sur_vis_cld = solar_top*itf_atm_vu_c / (1._wp-alb_sa_vu_c)
      swr_dw_sur_nir_cld = solar_top*itf_atm_ir_c / (1._wp-alb_sa_ir_c)

      dswd_dalb_vu_cld = solar_top * ((1._wp-alb_cld_vu)*alb_cld_vu*(1._wp-alb_cld_vu*alb_sa_vu_c) &
        - (1._wp-alb_cld_vu)*alb_sa_vu_c*alb_cld_vu*(-1._wp)*alb_cld_vu) / (1._wp-alb_cld_vu*alb_sa_vu_c)**2 &
        * itf_c_vu_d2*itf_w_vu_c_d2*itf_a_vu_d2*itf_o_vu_d2

      dswd_dalb_ir_cld = solar_top * ((1._wp-alb_cld_ir)*alb_cld_ir*(1._wp-alb_cld_ir*alb_sa_ir_c) &
        - (1._wp-alb_cld_ir)*alb_sa_ir_c*alb_cld_ir*(-1._wp)*alb_cld_ir) / (1._wp-alb_cld_ir*alb_sa_ir_c)**2 &
        * itf_c_ir_d2*itf_w_ir_c_d2*itf_a_ir_d2*itf_o_ir_d2

      ! d/dz (per cm) of the water path: w_col ~ exp(-z/h_q), times fp**n_wp ~ exp(-n_wp*z/hatm)
      drqh_dz = -rqh*(1._wp/h_q + n_wp/hatm)/100._wp
      dm_w_c_d1_dz = (exp_hc_hq/cos_zen+exp_1+exp_2) * drqh_dz 
      dm_w_c_d2_dz = dm_w_c_d1_dz + (2._wp*exp_2+exp_1) * drqh_dz 
      ditf_w_ir_c_d1_dz = a1_w*(-b1_w)*exp(-b1_w*m_w_c_d1)*dm_w_c_d1_dz +a2_w*(-b2_w)*exp(-b2_w*m_w_c_d1)*dm_w_c_d1_dz
      ditf_w_ir_c_d2_dz = a1_w*(-b1_w)*exp(-b1_w*m_w_c_d2)*dm_w_c_d2_dz +a2_w*(-b2_w)*exp(-b2_w*m_w_c_d2)*dm_w_c_d2_dz
      dswd_dz_ir_cld = solar_top * ((1._wp-alb_cld_ir)*itf_c_ir_d1*itf_a_ir_d1*itf_o_ir_d1 * ditf_w_ir_c_d1_dz + &
        (1._wp-alb_cld_ir)*alb_sa_ir_s*alb_cld_ir/(1._wp-alb_cld_ir*alb_sa_ir_s)*itf_c_ir_d2*itf_a_ir_d2*itf_o_ir_d2 * ditf_w_ir_c_d2_dz)
      dswd_dz_ir_cld = dswd_dz_ir_cld*100._wp  ! convert from W/m2/cm to W/m2/m

    endif

    ! final calculation

    solar_sur = (1._wp-cld)*solar_sur_s+cld*solar_sur_c


    return

  end subroutine sw_radiation_col


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function :  r _ s c a t _ m u
  !   Purpose  :  Rayleigh scattering albedo (VU) for cosine of zenith angle mu;
  !               i_swr_ray==1: thin-layer zenith dependence, r_scat is the value at mu=0.5
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function r_scat_mu(mu) result(r)

    implicit none

    real(wp), intent(in) :: mu
    real(wp) :: r

    if (i_swr_ray.eq.1) then
      r = r_scat*(1._wp-exp(-tau_ray/mu))/(1._wp-exp(-tau_ray/0.5_wp))
    else
      r = r_scat
    endif

  end function r_scat_mu


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Functions :  i t f _ g a s _ v u ,  i t f _ g a s _ i r
  !   Purpose   :  integral transmission function of the non-water gases for a slant path m
  !                (in units of the vertical column); i_swr_gas==0: constant O3 transmission
  !                c_itf_o in VU and none in IR
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function itf_gas_vu(m) result(itf)

    implicit none

    real(wp), intent(in) :: m
    real(wp) :: itf

    if (i_swr_gas.eq.1) then
      itf = max(0._wp, 1._wp-a_gas_vu*m**b_gas_vu)
    else
      itf = c_itf_o
    endif

  end function itf_gas_vu

  function itf_gas_ir(m) result(itf)

    implicit none

    real(wp), intent(in) :: m
    real(wp) :: itf

    if (i_swr_gas.eq.1) then
      itf = max(0._wp, 1._wp-a_gas_ir*m**b_gas_ir)
    else
      itf = 1._wp
    endif

  end function itf_gas_ir

end module sw_radiation_mod



      

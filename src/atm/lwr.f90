!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l w _ r a d i a t i o n _ m o d
!
!  Purpose : longwave radiation 
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
module lw_radiation_mod

  use atm_params, only : wp
  use constants, only : fqsat, sigma
  use control, only : co2_ref, ch4_ref, n2o_ref
  use atm_params, only : p0, ra, hatm, hpbl, c_hcld_dz_min
  use atm_params, only : i_lw_cld, c_lw_clot
  use atm_params, only : ak_co2, beta_co2
  use atm_params, only : ak_wv, a_vap, beta_vap, a2_vap, beta2_vap, a3_vap, rh_strat
  use atm_params, only : l_o3
  use atm_params, only : i_lwr_band, nu_win1, nu_win2, nu_co21, nu_co22, k_cont, c_cont_T, k_wvw, ak_wvw, f_wv_co2
  use atm_params, only : ak_co2_3b, beta_co2_3b, a0_co2_3b, a1_co2_3b
  use atm_params, only : ak_wv_3b, a_vap_3b, beta_vap_3b, a2_vap_3b, beta2_vap_3b, a3_vap_3b
  use atm_params, only : ak_o3_3b, a_o3_3b, beta_o3_3b
  use atm_params, only : k_co2w, k_co2r, b_co2w, c_co2_hi, z_top_lw
  use atm_grid, only : im, jm, kmc, nm, i_ocn
  use atm_grid, only : zl, llwr, llwr1, llwr2, llwr3, llwr4, nlwr1, nlwr2, nlwr3, nlwr4
  use vesta_mod, only : t_prof, rh_prof
  !$ use omp_lib

  implicit none

  real(wp), parameter :: emis=1._wp      !! atmospheric emissivity
  real(wp) :: h0    !! height scale [cm]
  real(wp), parameter :: beta0=1.66_wp
  ! CO2 parameters
  real(wp), parameter :: a0_co2=0.247_wp
  real(wp), parameter :: a1_co2=0.755_wp
  real(wp) :: q_co2
  ! ozone
  real(wp), parameter :: ak_o3=0.6_wp   ! from tuning of total LW contribution by O3
  real(wp), parameter :: a_o3=8.246_wp
  real(wp), parameter :: beta_o3=0.539_wp


  logical :: fb_q  ! flag for feedback analysis

  ! three-band LW (i_lwr_band==1): bands 1 = window (nu_win1-nu_win2), 2 = CO2 (nu_co21-nu_co22),
  ! 3 = rest (water vapour); fraction of sigma*T^4 in each band tabulated as a function of temperature
  integer, parameter :: nband = 3
  real(wp), parameter :: tpf_min = 100._wp, dtpf = 0.5_wp
  integer, parameter :: ntpf = 601
  real(wp) :: pfrac_tab(ntpf,nband)
  logical :: l_pfrac_init = .false.

  private
  public :: lw_radiation

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l w _ r a d i a t i o n 
  !   Purpose    :  driver for LW radiation 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lw_radiation(ecs_scale, frst, zsa, zs, htrop, hcld, ra2, gams, gamb, gamt, tam, ram, hrm, ttrop, cld, clot, &
      co2, ch4, n2o, cfc11, cfc12, co2e, o3, flwr_up_sur, &    ! in
      lwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, lwr_top, lwr_top_cs, lwr_top_cld, lwr_tro, lwr_cld, &     ! out
      gams_q, gamb_q, gamt_q, tam_q, ttrop_q, htrop_q, & ! optional input arguments for feedback analysis
      lwr_tro_cs) ! optional output, clear-sky flux at the tropopause

    implicit none

    real(wp), intent(in ) :: ecs_scale
    real(wp), intent(in ) :: frst(:,:,:)
    real(wp), intent(in ) :: zsa(:,:)
    real(wp), intent(in ) :: zs(:,:,:)
    real(wp), intent(in ) :: htrop(:,:)
    real(wp), intent(in ) :: hcld(:,:)
    real(wp), intent(in ) :: ra2(:,:,:)
    real(wp), intent(in ) :: gams(:,:)
    real(wp), intent(in ) :: gamb(:,:)
    real(wp), intent(in ) :: gamt(:,:)
    real(wp), intent(in ) :: tam(:,:)
    real(wp), intent(in ) :: ram(:,:)
    real(wp), intent(in ) :: hrm(:,:)
    real(wp), intent(in ) :: ttrop(:,:)
    real(wp), intent(in ) :: cld(:,:)
    real(wp), intent(in ) :: clot(:,:)
    real(wp), intent(in ) :: co2
    real(wp), intent(in ) :: ch4
    real(wp), intent(in ) :: n2o
    real(wp), intent(in ) :: cfc11
    real(wp), intent(in ) :: cfc12
    real(wp), intent(out) :: co2e
    real(wp), intent(in ) :: o3(:,:,:)
    real(wp), intent(in ) :: flwr_up_sur(:,:,:)

    real(wp), intent(out) :: lwr_sur(:,:)
    real(wp), intent(out) :: flwr_dw_sur(:,:,:)
    real(wp), intent(out) :: flwr_dw_sur_cs(:,:,:)
    real(wp), intent(out) :: flwr_dw_sur_cld(:,:,:)
    real(wp), intent(out) :: lwr_top(:,:)
    real(wp), intent(out) :: lwr_top_cs(:,:)
    real(wp), intent(out) :: lwr_top_cld(:,:)
    real(wp), intent(out) :: lwr_tro(:,:)
    real(wp), intent(out) :: lwr_cld(:,:)

    ! optional arguments needed for feedback analysis
    real(wp), intent(in ), optional :: gams_q(:,:)
    real(wp), intent(in ), optional :: gamb_q(:,:)
    real(wp), intent(in ), optional :: gamt_q(:,:)
    real(wp), intent(in ), optional :: tam_q(:,:)
    real(wp), intent(in ), optional :: ttrop_q(:,:)
    real(wp), intent(in ), optional :: htrop_q(:,:)

    ! optional output needed to diagnose the clear-sky radiative forcing at the tropopause
    real(wp), intent(out), optional :: lwr_tro_cs(:,:)

    ! local variables
    integer :: i, j, n
    real(wp), allocatable, dimension(:) :: zlwr
    real(wp), allocatable, dimension(:) :: tlwr
    real(wp), allocatable, dimension(:) :: qlwr
    real(wp), allocatable, dimension(:) :: O3lwr
    real(wp), allocatable, dimension(:,:) :: DCS
    real(wp), allocatable, dimension(:,:) :: DCL
    real(wp), allocatable, dimension(:) :: BSB
    real(wp), allocatable, dimension(:) :: fcl_up
    real(wp), allocatable, dimension(:) :: fcl_dw
    real(wp), allocatable, dimension(:) :: fcs_up
    real(wp), allocatable, dimension(:) :: fcs_dw
    real(wp), allocatable, dimension(:) :: fbl_up
    real(wp), allocatable, dimension(:) :: fbl_dw
    real(wp), allocatable, dimension(:) :: fbs_up
    real(wp), allocatable, dimension(:) :: fbs_dw
    real(wp), allocatable, dimension(:,:) :: BSB3
    real(wp), allocatable, dimension(:,:,:) :: DCS3
    real(wp), allocatable, dimension(:,:,:) :: DCL3
    integer :: ib
    real(wp) :: tsur_eff
    real(wp), dimension(nband) :: pf_sur
    integer :: n_col
    logical :: l_new_col

    real(wp), dimension(nm) :: fst
    real(wp), dimension(nm) :: flwr_sur
    real(wp), dimension(nm) :: flwr_top
    real(wp), dimension(nm) :: flwr_top_cs
    real(wp), dimension(nm) :: flwr_top_cld
    real(wp), dimension(nm) :: flwr_tro
    real(wp), dimension(nm) :: flwr_tro_cs
    real(wp), dimension(nm) :: flwr_cld

    real(wp) :: co2_bar
    real(wp) :: ch4_bar
    real(wp) :: n2o_bar
    real(wp) :: rf_co2
    real(wp) :: rf_ch4
    real(wp) :: rf_n2o
    real(wp) :: rf_cfc11
    real(wp) :: rf_cfc12

    real(wp), parameter :: a1 = -2.4e-7_wp     ! W/m2/ppm
    real(wp), parameter :: b1 = 7.2e-4_wp      ! W/m2/ppm
    real(wp), parameter :: c1 = -2.1e-4_wp     ! W/m2/ppb
    real(wp), parameter :: a2 = -8.0e-6_wp     ! W/m2/ppm
    real(wp), parameter :: b2 = 4.2e-6_wp      ! W/m2/ppb
    real(wp), parameter :: c2 = -4.9e-6_wp     ! W/m2/ppb
    real(wp), parameter :: a3 = -1.3e-6_wp     ! W/m2/ppb  
    real(wp), parameter :: b3 = -8.2e-6_wp     ! W/m2/ppb


    ! check wether feedback analysis has to be performed
    if (present(gams_q)) then
      fb_q = .true.
    else
      fb_q = .false.
    endif

    ! allocate
    allocate(zlwr(llwr))
    allocate(tlwr(llwr))
    allocate(qlwr(llwr))
    allocate(O3lwr(llwr))
    allocate(DCS(llwr,llwr))
    allocate(DCL(llwr,llwr))
    allocate(BSB(llwr))
    allocate(fcl_up(llwr))
    allocate(fcl_dw(llwr))
    allocate(fcs_up(llwr))
    allocate(fcs_dw(llwr))
    allocate(fbl_up(llwr))
    allocate(fbl_dw(llwr))
    allocate(fbs_up(llwr))
    allocate(fbs_dw(llwr))
    allocate(BSB3(llwr,nband))
    allocate(DCS3(llwr,llwr,nband))
    allocate(DCL3(llwr,llwr,nband))

    h0 = hatm*100._wp   ! cm

    if (i_lwr_band.eq.1 .and. .not.l_pfrac_init) call pfrac_init

    !-----------------------------------------------------------
    ! compute effective CO2 concentration for longwave radiation
    !-----------------------------------------------------------

    ! compute radiative forcings following Table 1 in Etminan et al., 2016 
    co2_bar = 0.5_wp*(co2+co2_ref)
    ch4_bar = 0.5_wp*(ch4+ch4_ref)
    n2o_bar = 0.5_wp*(n2o+n2o_ref)
    rf_co2 = (a1*(co2-co2_ref)**2+b1*abs(co2-co2_ref)+c1*n2o_bar+5.36_wp) * log(co2/co2_ref)
    rf_n2o = (a2*co2_bar+b2*n2o_bar+c2*ch4_bar+0.117_wp) * (sqrt(n2o)-sqrt(n2o_ref))
    rf_ch4 = (a3*min(5000._wp,ch4_bar)+b3*n2o_bar+0.043_wp) * (sqrt(ch4)-sqrt(ch4_ref))
    rf_cfc11 = 0.25_wp*1.e-3_wp*cfc11   ! Myhre et al., 1998, Table 3
    rf_cfc12 = 0.33_wp*1.e-3_wp*cfc12   ! Myhre et al., 1998, Table 3

    ! convert greenhouse gases concentrations into an effective CO2 concentration for longwave radiation
    co2e = co2_ref*exp((rf_co2+rf_ch4+rf_n2o+rf_cfc11+rf_cfc12)/(a1*(co2-co2_ref)**2+b1*abs(co2-co2_ref)+c1*n2o_bar+5.36_wp))   ! ppmv

    ! additional scaling of equivalent CO2 to mimic different climate sensitivity
    co2e = exp(ecs_scale*log(co2e)+(1._wp-ecs_scale)*log(co2_ref))

    ! convert from ppm to mass mixing ratio
    q_co2 = co2e*1.e-6_wp * 44.0095_wp/28.97_wp ! kg/kg

    !$omp parallel do collapse(2) private(i, j, n, zlwr, tlwr, qlwr, O3lwr, DCS, DCL, BSB, fcl_up, fcl_dw, fcs_up, fcs_dw) &
    !$omp private (fst, flwr_sur, flwr_top, flwr_top_cs, flwr_top_cld, flwr_tro, flwr_tro_cs, flwr_cld) &
    !$omp private (fbl_up, fbl_dw, fbs_up, fbs_dw, BSB3, DCS3, DCL3, ib, tsur_eff, pf_sur, n_col, l_new_col)
    do i=1,im
      do j=1,jm 

        n_col = 0   ! last surface type for which the column and the transmission functions were computed
        do n=1,nm
          fst(n) = frst(i,j,n)
          if (fst(n).gt.0._wp .or. n.eq.i_ocn) then     ! always compute it over the ocean (for SEMI)

            ! the column and the transmission functions depend on the surface type only through zs and ra2,
            ! reuse them if these are the same as for the previous surface type (e.g. sea ice and ocean)
            l_new_col = .true.
            if (n_col.gt.0) l_new_col = zs(i,j,n).ne.zs(i,j,n_col) .or. ra2(i,j,n).ne.ra2(i,j,n_col)
            if (l_new_col) n_col = n

            !-----------------------------------------------------------
            ! Atmospheric characteristics at vertical levels
            if (.not.l_new_col) then
              ! same as for n_col
            else if (.not.fb_q) then
              call lwr_column(zsa(i,j), zs(i,j,n), htrop(i,j), hcld(i,j), &    ! in
                gams(i,j), gamb(i,j), gamt(i,j), tam(i,j), ram(i,j), hrm(i,j), ttrop(i,j), o3(i,j,:), &   ! in
                zlwr, tlwr, qlwr, O3lwr)    ! out
            else
              call lwr_column(zsa(i,j), zs(i,j,n), htrop(i,j), hcld(i,j), &    ! in
                gams(i,j), gamb(i,j), gamt(i,j), tam(i,j), ram(i,j), hrm(i,j), ttrop(i,j), o3(i,j,:), &   ! in
                zlwr, tlwr, qlwr, O3lwr, &    ! out
                gams_q(i,j), gamb_q(i,j), gamt_q(i,j), tam_q(i,j), ttrop_q(i,j), htrop_q(i,j)) ! optional input arguments for feedback analysis
            endif

            if (i_lwr_band.eq.0) then

              !-----------------------------------------------------------
              ! Computation of clear-sky and cloudy LWR transmission functions
              if (l_new_col) call lwr_transfer(ra2(i,j,n), clot(i,j), zlwr, tlwr, qlwr, O3lwr, &
                BSB, DCS, DCL)

              !-----------------------------------------------------------
              ! LWR fluxes for clear-sky atmosphere  
              call lwr_clear_sky(flwr_up_sur(i,j,n), BSB, DCS, &
                fcs_up, fcs_dw)

              !-----------------------------------------------------------
              ! LWR fluxes for cloudy conditon
              call lwr_clouds(flwr_up_sur(i,j,n), BSB, DCL, &
                fcl_up, fcl_dw)

            else

              !-----------------------------------------------------------
              ! three Planck-weighted bands, fluxes summed over the bands;
              ! the surface emission is split with the Planck fraction at the surface radiative temperature
              tsur_eff = (flwr_up_sur(i,j,n)/sigma)**0.25_wp
              call pfrac3(tsur_eff, pf_sur)
              if (l_new_col) call lwr_transfer_3b(ra2(i,j,n), clot(i,j), zlwr, tlwr, qlwr, O3lwr, &
                BSB3, DCS3, DCL3)
              fcs_up = 0._wp
              fcs_dw = 0._wp
              fcl_up = 0._wp
              fcl_dw = 0._wp
              do ib=1,nband
                call lwr_clear_sky(flwr_up_sur(i,j,n)*pf_sur(ib), BSB3(:,ib), DCS3(:,:,ib), &
                  fbs_up, fbs_dw)
                call lwr_clouds(flwr_up_sur(i,j,n)*pf_sur(ib), BSB3(:,ib), DCL3(:,:,ib), &
                  fbl_up, fbl_dw)
                fcs_up = fcs_up + fbs_up
                fcs_dw = fcs_dw + fbs_dw
                fcl_up = fcl_up + fbl_up
                fcl_dw = fcl_dw + fbl_dw
              enddo

            endif

            !-----------------------------------------------------------
            ! total LWR fluxes 
            call lwr_total(cld(i,j), fcs_up, fcs_dw, fcl_up, fcl_dw, &
              flwr_sur(n), flwr_dw_sur(i,j,n), flwr_dw_sur_cs(i,j,n), flwr_dw_sur_cld(i,j,n), &
              flwr_top(n), flwr_top_cs(n), flwr_top_cld(n), flwr_tro(n), flwr_tro_cs(n), flwr_cld(n))

          else

              flwr_sur(n)            = 0._wp
              flwr_dw_sur(i,j,n)     = 0._wp
              flwr_dw_sur_cs(i,j,n)  = 0._wp
              flwr_dw_sur_cld(i,j,n) = 0._wp
              flwr_top(n)            = 0._wp
              flwr_top_cs(n)         = 0._wp
              flwr_top_cld(n)        = 0._wp
              flwr_tro(n)            = 0._wp
              flwr_tro_cs(n)         = 0._wp
              flwr_cld(n)            = 0._wp

          endif
        enddo              

        ! grid cell averages
        lwr_sur(i,j)     = sum(flwr_sur(:)*fst(:))
        lwr_top(i,j)     = sum(flwr_top(:)*fst(:))
        lwr_top_cs(i,j)  = sum(flwr_top_cs(:)*fst(:))
        lwr_top_cld(i,j) = sum(flwr_top_cld(:)*fst(:))
        lwr_tro(i,j)     = sum(flwr_tro(:)*fst(:))
        if (present(lwr_tro_cs)) lwr_tro_cs(i,j) = sum(flwr_tro_cs(:)*fst(:))
        lwr_cld(i,j)     = sum(flwr_cld(:)*fst(:))

      enddo
    enddo 
    !$omp end parallel do

    deallocate(zlwr)
    deallocate(tlwr)
    deallocate(qlwr)
    deallocate(O3lwr)
    deallocate(DCS)
    deallocate(DCL)
    deallocate(BSB)
    deallocate(fcl_up)
    deallocate(fcl_dw)
    deallocate(fcs_up)
    deallocate(fcs_dw)
    deallocate(fbl_up)
    deallocate(fbl_dw)
    deallocate(fbs_up)
    deallocate(fbs_dw)
    deallocate(BSB3)
    deallocate(DCS3)
    deallocate(DCL3)

    return

  end subroutine lw_radiation


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l w r _ t o t a l
  !   Purpose    :  combine clear sky and cloudy longwave radiation
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lwr_total(cld, fcs_up, fcs_dw, fcl_up, fcl_dw, &
      flwr_sur, flwr_dw_sur, flwr_dw_sur_cs, flwr_dw_sur_cld, &
      flwr_top, flwr_top_cs, flwr_top_cld, flwr_tro, flwr_tro_cs, flwr_cld)

    implicit none

    real(wp), intent(in ) :: cld
    real(wp), intent(in ) :: fcs_up(:)  
    real(wp), intent(in ) :: fcs_dw(:)
    real(wp), intent(in ) :: fcl_up(:)
    real(wp), intent(in ) :: fcl_dw(:)

    real(wp), intent(out) :: flwr_sur
    real(wp), intent(out) :: flwr_dw_sur
    real(wp), intent(out) :: flwr_dw_sur_cs
    real(wp), intent(out) :: flwr_dw_sur_cld
    real(wp), intent(out) :: flwr_top
    real(wp), intent(out) :: flwr_top_cs
    real(wp), intent(out) :: flwr_top_cld
    real(wp), intent(out) :: flwr_tro
    real(wp), intent(out) :: flwr_tro_cs
    real(wp), intent(out) :: flwr_cld

    integer :: k
    real(wp) :: flwr_dw_tro, flwr_up_tro
    real(wp) :: flwr_dw_cld, flwr_up_cld

    ! fluxes at the surface
    flwr_sur = (1._wp-cld)*(fcs_dw(1)-fcs_up(1)) + cld*(fcl_dw(1)-fcl_up(1))
    flwr_dw_sur = (1._wp-cld)*fcs_dw(1) + cld*fcl_dw(1) 
    flwr_dw_sur_cs  = fcs_dw(1) 
    flwr_dw_sur_cld = fcl_dw(1) 

    ! net fluxes at TOA, positive down
    flwr_top = -((1._wp-cld)*fcs_up(llwr) + cld*fcl_up(llwr))
    flwr_top_cs  = -fcs_up(llwr)
    flwr_top_cld = -fcl_up(llwr)

    ! fluxes at the tropopause
    k = llwr3       
    flwr_up_tro = (1._wp-cld)*fcs_up(k) + cld*fcl_up(k)
    flwr_dw_tro = (1._wp-cld)*fcs_dw(k) + cld*fcl_dw(k)
    flwr_tro = flwr_dw_tro-flwr_up_tro
    ! clear-sky net flux at the tropopause, needed for the clear-sky radiative forcing
    flwr_tro_cs = fcs_dw(k)-fcs_up(k)

    ! fluxes at cloud base
    k = llwr1
    flwr_up_cld = (1._wp-cld)*fcs_up(k) + cld*fcl_up(k)
    flwr_dw_cld = (1._wp-cld)*fcs_dw(k) + cld*fcl_dw(k)
    flwr_cld = flwr_dw_cld-flwr_up_cld

    return

  end subroutine lwr_total
    

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l w r _ c o l u m n
  !   Purpose    :  derive temperature, humidity and ozone at lwr levels
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lwr_column(zsa, zs, htrop, hcld, gams, gamb, gamt, tam, ram, hrm, ttrop, O3, &
      zlwr, tlwr, qlwr, O3lwr, &    
      gams_q, gamb_q, gamt_q, tam_q, ttrop_q, htrop_q) ! optional input arguments for feedback analysis
    
    implicit none

    real(wp), intent(in ) :: zsa
    real(wp), intent(in ) :: zs
    real(wp), intent(in ) :: htrop
    real(wp), intent(in ) :: hcld
    real(wp), intent(in ) :: gams
    real(wp), intent(in ) :: gamb
    real(wp), intent(in ) :: gamt
    real(wp), intent(in ) :: tam
    real(wp), intent(in ) :: ram
    real(wp), intent(in ) :: hrm
    real(wp), intent(in ) :: ttrop
    real(wp), intent(in ) :: O3(:)

    real(wp), intent(out), dimension(:) :: zlwr
    real(wp), intent(out), dimension(:) :: tlwr
    real(wp), intent(out), dimension(:) :: qlwr
    real(wp), intent(out), dimension(:) :: O3lwr

    ! optional arguments needed for feedback analysis
    real(wp), intent(in), optional :: gams_q
    real(wp), intent(in), optional :: gamb_q
    real(wp), intent(in), optional :: gamt_q
    real(wp), intent(in), optional :: tam_q
    real(wp), intent(in), optional :: ttrop_q
    real(wp), intent(in), optional :: htrop_q

    integer :: k, kk
    real(wp) :: z_cld_bot, z_cld_top
    real(wp) :: dz_l1, dz_l2, dz_l3, dz_l4
    real(wp) :: tamz, tlwr_q, rqlwr, qtrop
    real(wp) :: w


    ! Cloud parameters

    z_cld_bot = zs+hpbl
    z_cld_top = max(hcld,z_cld_bot+c_hcld_dz_min)
    z_cld_top = min(z_cld_top,htrop-1000._wp)

    ! Layers thickness

    dz_l1 = (z_cld_bot-zs)/(nlwr1-1._wp)
    dz_l2 = (z_cld_top-z_cld_bot)/nlwr2
    dz_l3 = (htrop-z_cld_top)/nlwr3      
    dz_l4 = (z_top_lw-htrop)/nlwr4

    ! levels 

    ! k=1 surface
    ! k=llwr1 cloud bottom
    ! k=llwr2 cloud top
    ! k=llwr3 tropopause
    ! k=llwr4 top of the atmosphere

    zlwr(1)=zs  

    do k=2,llwr1
      zlwr(k) = zlwr(k-1)+dz_l1
    enddo 

    do k=llwr1+1,llwr2
      zlwr(k) = zlwr(k-1)+dz_l2
    enddo 

    do k=llwr2+1,llwr3
      zlwr(k) = zlwr(k-1)+dz_l3
    enddo 

    do k=llwr3+1,llwr4
      zlwr(k) = zlwr(k-1)+dz_l4
    enddo  

    ! Temperature and humidity at levels

    ! Surface

    ! tam adjusted for surface elevation difference from the mean
    tamz    = t_prof(zsa, zlwr(1), tam, gams, gamb, gamt, htrop, 0)
    tlwr(1) = tamz 
    qlwr(1) = fqsat(tlwr(1),p0*exp(-zlwr(1)/hatm))*ram

    ! Troposphere

    do k=2,llwr3
      tlwr(k) = t_prof(zs, zlwr(k), tamz, gams, gamb, gamt, htrop, 1)
      rqlwr   = rh_prof(zs, zlwr(k), ram, hrm, htrop)
      qlwr(k) = fqsat(tlwr(k),p0*exp(-zlwr(k)/hatm))*rqlwr 
    enddo

    ! Stratosphere      

    qtrop = fqsat(ttrop,p0*exp(-htrop/hatm))*rh_strat
    do k=llwr3+1,llwr4
      tlwr(k) = ttrop
      qlwr(k) = qtrop 
    enddo

    ! for feedback analysis, temperature used to compute moisture profile can differ from actual temperature profile
    if (fb_q) then
      tamz = t_prof(zsa, zlwr(1), tam_q, gams_q, gamb_q, gamt_q, htrop_q, 0)
      qtrop = FQSAT(ttrop_q,p0*exp(-htrop_q/hatm))*rh_strat
      do k=2,llwr4
        ! fixme possible problems at tropopause because of non-adaptive layers?
        if (zlwr(k).le.(htrop_q+10.)) then
          ! troposphere
          tlwr_q = t_prof(zs, zlwr(k), tamz, gams_q, gamb_q, gamt_q, htrop_q, 1)
          rqlwr  = rh_prof(zs, zlwr(k), ram, hrm, htrop_q)
          qlwr(k)= FQSAT(tlwr_q,p0*exp(-zlwr(k)/hatm))*rqlwr 
        else
          ! stratosphere
          qlwr(k) = qtrop
        endif
      enddo
    endif

    if (l_o3) then
      ! interpolate O3 to lwr levels

      do k=1,llwr4
        if (zlwr(k)<=zl(1)) then
          w = (zlwr(k)-zl(1)) / (zl(2)-zl(1))
          o3lwr(k) = (1._wp-w)*o3(1) + w*o3(2)
        else if (zlwr(k)>=zl(kmc)) then
          ! above the top atmospheric level (30 km) O3 is kept at its value there
          o3lwr(k) = o3(kmc)
        else
          do kk=2,kmc
            if (zlwr(k)>=zl(kk-1) .and. zlwr(k)<=zl(kk)) then
              w = (zlwr(k)-zl(kk-1)) / (zl(kk)-zl(kk-1))
              o3lwr(k) = (1._wp-w)*o3(kk-1) + w*o3(kk)
              exit
            endif
          enddo
        endif
      enddo

    else
      O3lwr(:) = 0._wp
    endif

    return

  end subroutine lwr_column


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l w r _ t r a n s f e r
  !   Purpose    :  computation of tranfer functions for lwr
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lwr_transfer(ra2, clot, zlwr, tlwr, qlwr, O3lwr, &
      BSB, DCS, DCL)

    implicit none
    
    real(wp), intent(in) :: ra2
    real(wp), intent(in) :: clot
    real(wp), intent(in) :: zlwr(:)
    real(wp), intent(in) :: tlwr(:)
    real(wp), intent(in) :: qlwr(:)
    real(wp), intent(in) :: O3lwr(:)

    real(wp), intent(out) :: BSB(:)
    real(wp), intent(out) :: DCS(:,:)
    real(wp), intent(out) :: DCL(:,:)

    integer :: k, l

    ! size of these arrays has to be larger than llwr!
    real(wp), dimension(20) :: expc
    real(wp), dimension(20) :: am_cld
    real(wp), dimension(20) :: am_o3
    real(wp), dimension(20) :: am_wv
    real(wp), dimension(20) :: am_co2
    real(wp), dimension(20) :: cum_cld, cum_o3, cum_wv, cum_co2   ! absorber amounts from level k up to level l

    real(wp) :: rhos, zsur
    real(wp) :: zlsk 
    real(wp) :: Z1,Z2,Zm,dZ,hq
    real(wp) :: ql1,ql2
    real(wp) :: kappa_co2, kappa_wv
    real(wp) :: am_cld_kl, am_o3_kl, am_wv_kl, am_co2_kl
    real(wp) :: d_co2, d_vap, d_o3, d_cld

    ! Conversion to CGS units

    rhos = ra2*0.001_wp    ! g/cm3

    zsur = zlwr(1)*100._wp

    kappa_co2 = (ak_co2+1._wp)/h0  ! 1/cm

    ! 1) Calculations of intermediate values at the levels

    do k=1,llwr

      zlsk = zlwr(k)*100._wp    ! cm
      expc(k) = exp(-kappa_co2*zlsk)
      BSB(k) = emis*sigma*tlwr(k)**4

    enddo 

    ! 2) Calculations of absorptions in the layers

    do k=1,llwr-1

      Z1 = zlwr(k)*100._wp    ! cm
      Z2 = zlwr(k+1)*100._wp  ! cm
      Zm = 0.5_wp*(Z1+Z2)     ! cm
      dZ = Z2-Z1              ! cm

      ! 2.1)  Influence of water vapour

      ql1 = qlwr(k)
      ql2 = qlwr(k+1)
      if (ql1.gt.ql2 .and. ql2.gt.0._wp) then
        ! derive local moisture height scale assuming exponential profile between z1 and z2
        hq = dZ/log(ql1/ql2)    ! cm
        hq = min(hq,h0)
      else
        hq = h0
      endif       
      kappa_wv = (ak_wv+1._wp)/h0 + 1._wp/hq
      am_wv(k) = rhos*ql1*exp(z1/hq+zsur/h0)/kappa_wv*(exp(-kappa_wv*z1)-exp(-kappa_wv*z2))   ! g/cm3 * g/g * cm = g/cm2 

      ! 2.2)  Influence of CO2

      am_co2(k) = q_co2/kappa_co2*(expc(k)-expc(k+1))  ! cm

      ! 2.3) Influence of ozone

      if (l_o3) then
        am_o3(k) = ra*1.e-3_wp*exp(-Zm*(ak_o3+1._wp)/h0)*0.5_wp*(O3lwr(k)+O3lwr(k+1))*dZ  ! g/cm2
      else
        am_o3(k) = 0._wp
      endif       

      ! 2.4) Influence of clouds

      if (i_lw_cld.eq.1) then
        ! set to zero, clouds are considered later
        am_cld(k) = 0._wp
      else if (i_lw_cld.eq.2) then
        ! account for clouds directly in the transmission function
        if (k>=llwr1 .and. k<llwr2) then
          ! in cloud layer
          am_cld(k) = clot/nlwr2 
        else
          ! outside clouds
          am_cld(k) = 0._wp
        endif
      endif

    enddo
   
    ! 3) Calculation of the integral transmision functions

    do k=1,llwr-1

      ! absorber amounts between levels k and l, accumulated first so that the loop over l
      ! below has no dependency and can be vectorised (it is dominated by the power laws)
      am_cld_kl = 0._wp
      am_o3_kl  = 0._wp
      am_wv_kl  = 0._wp
      am_co2_kl = 0._wp
      do l=k+1,llwr
        am_cld_kl = am_cld_kl + am_cld(l-1)
        am_o3_kl  = am_o3_kl  + am_o3(l-1)
        am_wv_kl  = am_wv_kl  + am_wv(l-1)
        am_co2_kl = am_co2_kl + am_co2(l-1)
        cum_cld(l) = am_cld_kl
        cum_o3(l)  = am_o3_kl
        cum_wv(l)  = am_wv_kl
        cum_co2(l) = am_co2_kl
      enddo

      do l=k+1,llwr

        am_cld_kl = cum_cld(l)
        am_o3_kl  = cum_o3(l)
        am_wv_kl  = cum_wv(l)
        am_co2_kl = cum_co2(l)

        ! equation (6.7) in PIK report 81
        d_o3 = 1.0_wp-a_o3*(am_o3_kl**beta_o3)      

        ! equation (6.5) in PIK report 81
        d_vap = 1.0_wp/(1.0_wp+a_vap*((beta0*am_wv_kl)**beta_vap)+a2_vap*((beta0*am_wv_kl)**beta2_vap)+a3_vap*((beta0*am_wv_kl)**3))

        ! equation (6.6) in PIK report 81, valid for a broad range of CO2 concentrations (up to 20 times present day)
        ! modified with additional factor (1.-0.1*(am_co2_kl/1000.)^2) to get increasing CO2 radiative forcing for increasing CO2 
        ! (in agreement with Hansen 2005, Table 1 and Colman & McAvaney 2009, Fig. 1)
        d_co2 = (1._wp-min(0.2_wp,0.1_wp*(am_co2_kl/1000._wp)**2))*(1._wp+a0_co2*a1_co2*((beta0*am_co2_kl)**beta_co2))/(1._wp+a0_co2*((beta0*am_co2_kl)**beta_co2))

        ! clouds, using optical thickness
        d_cld = exp(-c_lw_clot*am_cld_kl)

        ! clear sky transmission function
        DCS(l,k) = d_vap*d_co2*d_o3

        ! cloudy sky transmission function
        DCL(l,k) = d_vap*d_co2*d_o3*d_cld

      enddo 
    enddo

    ! symmetric matrix
    do l=1, llwr-1
      do k=l+1, llwr
        DCS(l,k) = DCS(k,l)
        DCL(l,k) = DCL(k,l)
      enddo
    enddo

    ! diagonal terms
    do k=1, llwr
      DCS(k,k) = 1.0_wp
      DCL(k,k) = 1.0_wp
    enddo 

    return

  end subroutine lwr_transfer


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p f r a c _ i n i t
  !   Purpose    :  tabulate the fraction of sigma*T^4 emitted in each band
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine pfrac_init

    implicit none

    integer :: it
    real(wp) :: T

    do it=1,ntpf
      T = tpf_min + (it-1)*dtpf
      pfrac_tab(it,1) = planck_frac(T, nu_win1, nu_win2)
      pfrac_tab(it,2) = planck_frac(T, nu_co21, nu_co22)
      pfrac_tab(it,3) = 1._wp - pfrac_tab(it,1) - pfrac_tab(it,2)
    enddo
    l_pfrac_init = .true.

    return

  end subroutine pfrac_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Function :  p l a n c k _ f r a c
  !   Purpose  :  fraction of sigma*T^4 between wavenumbers nu1 and nu2 (cm-1)
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function planck_frac(T, nu1, nu2) result(pf)

    implicit none

    real(wp), intent(in) :: T, nu1, nu2
    real(wp) :: pf

    integer, parameter :: n = 200
    real(wp), parameter :: c1 = 1.191042e-8_wp   ! 2hc^2, W/m2/sr/(cm-1)^4
    real(wp), parameter :: c2 = 1.4388_wp        ! hc/k, cm K
    integer :: i
    real(wp) :: nu, dnu, w

    dnu = (nu2-nu1)/n
    pf = 0._wp
    do i=0,n
      nu = nu1+i*dnu
      w = 1._wp
      if (i.eq.0 .or. i.eq.n) w = 0.5_wp
      pf = pf + w*nu**3/(exp(c2*nu/T)-1._wp)
    enddo
    pf = 3.14159265358979_wp*c1*pf*dnu/(sigma*T**4)

  end function planck_frac


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  p f r a c 3
  !   Purpose    :  band fractions of sigma*T^4 from the table, all bands
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine pfrac3(T, pf)

    implicit none

    real(wp), intent(in) :: T
    real(wp), intent(out) :: pf(nband)

    integer :: it
    real(wp) :: x

    x = (min(max(T,tpf_min),tpf_min+(ntpf-1)*dtpf)-tpf_min)/dtpf
    it = min(int(x)+1,ntpf-1)
    x = x-(it-1)
    pf(:) = (1._wp-x)*pfrac_tab(it,:) + x*pfrac_tab(it+1,:)

  end subroutine pfrac3


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l w r _ t r a n s f e r _ 3 b
  !   Purpose    :  Planck-weighted emission and transmission functions for the
  !                 three LW bands (i_lwr_band==1), all bands in one pass:
  !                 1 window:  water vapour self-continuum, weak water vapour lines, O3
  !                 2 CO2:     CO2 + water vapour (absorber amount scaled by f_wv_co2)
  !                 3 rest:    water vapour
  !                 Same absorber amounts and transmission forms as lwr_transfer, with the
  !                 band parameters *_3b. The powers of the absorber amounts are evaluated
  !                 as exp(b*log(u)) with one log per absorber shared by all bands, and the
  !                 cloud transmission is band independent
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lwr_transfer_3b(ra2, clot, zlwr, tlwr, qlwr, O3lwr, &
      BSB, DCS, DCL)

    implicit none
    
    real(wp), intent(in) :: ra2
    real(wp), intent(in) :: clot
    real(wp), intent(in) :: zlwr(:)
    real(wp), intent(in) :: tlwr(:)
    real(wp), intent(in) :: qlwr(:)
    real(wp), intent(in) :: O3lwr(:)

    real(wp), intent(out) :: BSB(:,:)
    real(wp), intent(out) :: DCS(:,:,:)
    real(wp), intent(out) :: DCL(:,:,:)

    integer :: k, l, ib, ip, np

    real(wp), dimension(20) :: expc
    real(wp), dimension(20) :: am_cld, am_o3, am_wv, am_co2, am_cnt, am_wvw
    real(wp), dimension(190) :: cum_cld, cum_o3, cum_wv, cum_co2, cum_cnt, cum_wvw
    real(wp), dimension(190) :: d_win, d_co2, d_wvr, d_cld
    real(wp), dimension(nband) :: pf

    real(wp) :: rhos, zsur
    real(wp) :: Z1,Z2,Zm,dZ,hq
    real(wp) :: ql1,ql2,qm,tm
    real(wp) :: kappa_co2, kappa_wv, kappa_wvw
    real(wp) :: u_cld, u_o3, u_wv, u_co2, u_cnt, u_wvw
    real(wp) :: x_wv, lx_wv, e1_wv, e2_wv, e3_wv, lx_co2, pw_co2, pc_co2
    real(wp) :: f1_wv_co2, f2_wv_co2, f3_wv_co2


    rhos = ra2*0.001_wp    ! g/cm3
    zsur = zlwr(1)*100._wp
    kappa_co2 = (ak_co2_3b+1._wp)/h0  ! 1/cm

    ! water vapour of the CO2 band: (f_wv_co2*u)**b = f_wv_co2**b * u**b
    f1_wv_co2 = f_wv_co2**beta_vap_3b
    f2_wv_co2 = f_wv_co2**beta2_vap_3b
    f3_wv_co2 = f_wv_co2**3

    ! 1) band emission at the levels
    do k=1,llwr
      expc(k) = exp(-kappa_co2*zlwr(k)*100._wp)
      call pfrac3(tlwr(k), pf)
      BSB(k,:) = emis*sigma*tlwr(k)**4*pf(:)
    enddo 

    ! 2) absorber amounts in the layers
    do k=1,llwr-1

      Z1 = zlwr(k)*100._wp    ! cm
      Z2 = zlwr(k+1)*100._wp  ! cm
      Zm = 0.5_wp*(Z1+Z2)     ! cm
      dZ = Z2-Z1              ! cm

      ql1 = qlwr(k)
      ql2 = qlwr(k+1)
      if (ql1.gt.ql2 .and. ql2.gt.0._wp) then
        hq = dZ/log(ql1/ql2)    ! cm
        hq = min(hq,h0)
      else
        hq = h0
      endif       

      ! window: weak water vapour lines, own pressure scaling
      kappa_wvw = (ak_wvw+1._wp)/h0 + 1._wp/hq
      am_wvw(k) = rhos*ql1*exp(z1/hq+zsur/h0)/kappa_wvw*(exp(-kappa_wvw*z1)-exp(-kappa_wvw*z2))   ! g/cm2
      ! window: self-continuum, proportional to q*e, temperature dependent
      tm = 0.5_wp*(tlwr(k)+tlwr(k+1))
      qm = sqrt(max(ql1,1.e-12_wp)*max(ql2,1.e-12_wp))
      am_cnt(k) = rhos*exp(-(Zm-zsur)/h0)*qm*(qm*exp(-Zm/h0)/0.622_wp)*exp(c_cont_T*(296._wp/tm-1._wp))*dZ   ! g/cm2
      ! window: ozone
      if (l_o3) then
        am_o3(k) = ra*1.e-3_wp*exp(-Zm*(ak_o3_3b+1._wp)/h0)*0.5_wp*(O3lwr(k)+O3lwr(k+1))*dZ  ! g/cm2
      else
        am_o3(k) = 0._wp
      endif
      ! water vapour band (the CO2 band sees f_wv_co2 times this)
      kappa_wv = (ak_wv_3b+1._wp)/h0 + 1._wp/hq
      am_wv(k) = rhos*ql1*exp(z1/hq+zsur/h0)/kappa_wv*(exp(-kappa_wv*z1)-exp(-kappa_wv*z2))   ! g/cm2
      ! CO2: main band (ib=2) and, through k_co2w/k_co2r, the band wings in the window and in the water vapour band
      am_co2(k) = q_co2/kappa_co2*(expc(k)-expc(k+1))  ! cm

      ! clouds
      if (i_lw_cld.eq.2 .and. k>=llwr1 .and. k<llwr2) then
        am_cld(k) = clot/nlwr2 
      else
        am_cld(k) = 0._wp
      endif

    enddo
   
    ! 3) absorber amounts between all pairs of levels k<l
    np = 0
    do k=1,llwr-1
      u_cld = 0._wp
      u_o3  = 0._wp
      u_wv  = 0._wp
      u_co2 = 0._wp
      u_cnt = 0._wp
      u_wvw = 0._wp
      do l=k+1,llwr
        u_cld = u_cld + am_cld(l-1)
        u_o3  = u_o3  + am_o3(l-1)
        u_wv  = u_wv  + am_wv(l-1)
        u_co2 = u_co2 + am_co2(l-1)
        u_cnt = u_cnt + am_cnt(l-1)
        u_wvw = u_wvw + am_wvw(l-1)
        np = np+1
        cum_cld(np) = u_cld
        cum_o3(np)  = u_o3
        cum_wv(np)  = u_wv
        cum_co2(np) = u_co2
        cum_cnt(np) = u_cnt
        cum_wvw(np) = u_wvw
      enddo
    enddo

    ! 4) integral transmission functions, one branch-free (vectorisable) loop over all pairs;
    ! zero absorber amounts are floored at tiny so that their powers vanish
    do ip=1,np

      ! powers of the water vapour amount, shared by the CO2 and water vapour bands
      x_wv  = max(beta0*cum_wv(ip),tiny(1._wp))
      lx_wv = log(x_wv)
      e1_wv = a_vap_3b*exp(beta_vap_3b*lx_wv)
      e2_wv = a2_vap_3b*exp(beta2_vap_3b*lx_wv)
      e3_wv = a3_vap_3b*x_wv**3

      ! powers of the CO2 amount, shared by all bands
      lx_co2 = log(max(beta0*cum_co2(ip),tiny(1._wp)))
      pw_co2 = exp(b_co2w*lx_co2)
      pc_co2 = exp(beta_co2_3b*lx_co2)

      ! 1 window
      d_win(ip) = exp(-beta0*(k_cont*cum_cnt(ip)+k_wvw*cum_wvw(ip)) - k_co2w*pw_co2) &
        * (1.0_wp-a_o3_3b*exp(beta_o3_3b*log(max(cum_o3(ip),tiny(1._wp)))))

      ! 2 CO2
      d_co2(ip) = 1.0_wp/(1.0_wp+f1_wv_co2*e1_wv+f2_wv_co2*e2_wv+f3_wv_co2*e3_wv) &
        * (1._wp-min(0.2_wp,c_co2_hi*(cum_co2(ip)/1000._wp)**2)) &
        * (1._wp+a0_co2_3b*a1_co2_3b*pc_co2)/(1._wp+a0_co2_3b*pc_co2)

      ! 3 water vapour
      d_wvr(ip) = exp(-k_co2r*pw_co2)/(1.0_wp+e1_wv+e2_wv+e3_wv)

      ! clouds, band independent
      d_cld(ip) = exp(-c_lw_clot*cum_cld(ip))

    enddo

    np = 0
    do k=1,llwr-1
      do l=k+1,llwr
        np = np+1
        DCS(l,k,1) = d_win(np)
        DCS(l,k,2) = d_co2(np)
        DCS(l,k,3) = d_wvr(np)
        DCL(l,k,1) = d_win(np)*d_cld(np)
        DCL(l,k,2) = d_co2(np)*d_cld(np)
        DCL(l,k,3) = d_wvr(np)*d_cld(np)
      enddo
    enddo

    do ib=1,nband
      ! symmetric matrix
      do l=1, llwr-1
        do k=l+1, llwr
          DCS(l,k,ib) = DCS(k,l,ib)
          DCL(l,k,ib) = DCL(k,l,ib)
        enddo
      enddo
      ! diagonal terms
      do k=1, llwr
        DCS(k,k,ib) = 1.0_wp
        DCL(k,k,ib) = 1.0_wp
      enddo 
    enddo

    return

  end subroutine lwr_transfer_3b


      
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l w r _ c l e a r _ s k y
  !   Purpose    :  computation of LWR fluxes for clear sky conditions
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lwr_clear_sky(flwr_up_sur, BSB, DCS, &
      fcs_up, fcs_dw)

    implicit none

    real(wp), intent(in ) :: flwr_up_sur
    real(wp), intent(in ) :: BSB(:)
    real(wp), intent(in ) :: DCS(:,:)

    real(wp), intent(out) :: fcs_up(:)
    real(wp), intent(out) :: fcs_dw(:)

    integer :: k, l

     
    ! Upward LWR fluxes
    fcs_up(1) = flwr_up_sur
    do k=2, llwr
      fcs_up(k) = BSB(k)+(fcs_up(1)-BSB(1))*DCS(k,1)
      do l=1, k-1
        fcs_up(k) = fcs_up(k)-(BSB(l+1)-BSB(l))*0.5_wp*(DCS(k,l+1)+DCS(k,l))
      enddo
    enddo

    ! Downward LWR fluxes  
    fcs_dw(llwr) = 0._wp
    do k=llwr-1, 1, -1
      fcs_dw(k) = BSB(k)-BSB(llwr)*DCS(k,llwr)
      do l=k, llwr-1
        fcs_dw(k) = fcs_dw(k)+(BSB(l+1)-BSB(l))*0.5_wp*(DCS(k,l)+DCS(k,l+1))
      enddo
    enddo

    return

  end subroutine lwr_clear_sky

     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  l w r _ c l o u d s
  !   Purpose    :  computation of LWR fluxes for cloudy conditions
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine lwr_clouds(flwr_up_sur, BSB, DCL, &
      fcl_up, fcl_dw)

    implicit none

    real(wp), intent(in ) :: flwr_up_sur
    real(wp), intent(in ) :: BSB(:)
    real(wp), intent(in ) :: DCL(:,:)

    real(wp), intent(out) :: fcl_up(:)
    real(wp), intent(out) :: fcl_dw(:)

    integer :: k, l


    if (i_lw_cld.eq.1) then
      ! consider clouds as black bodies

      ! Upward LWR fluxes
      fcl_up(1) = flwr_up_sur
      do k=2, llwr1
        fcl_up(k) = BSB(k)+(fcl_up(1)-BSB(1))*DCL(k,1)
        do l=1, k-1
          fcl_up(k) = fcl_up(k)-(BSB(l+1)-BSB(l))*0.5_wp*(DCL(k,l+1)+DCL(k,l))
        end do
      end do
      do k=llwr1+1,llwr2
        fcl_up(k) = BSB(k)
      enddo  
      do k=llwr2+1,llwr
        fcl_up(k) = BSB(k)
        do l=llwr2, k-1
          fcl_up(k) = fcl_up(k)-(BSB(l+1)-BSB(l))*0.5_wp*(DCL(k,l+1)+DCL(k,l))
        end do
      end do

      ! Downward LWR fluxes  
      fcl_dw(llwr) = 0._wp
      do k=llwr-1, llwr2, -1
        fcl_dw(k) = BSB(k)-BSB(llwr)*DCL(k,llwr)
        do l=k, llwr-1
          fcl_dw(k) = fcl_dw(k)+(BSB(l+1)-BSB(l))*0.5_wp*(DCL(k,l)+DCL(k,l+1))
        end do
      end do
      do k=llwr2-1,llwr1,-1
        fcl_dw(k) = BSB(k)
      enddo  
      do k=llwr1-1,1,-1
        fcl_dw(k) = BSB(k)
        do l=k,llwr1-1
          fcl_dw(k) = fcl_dw(k)+(BSB(l+1)-BSB(l))*0.5_wp*(DCL(k,l)+DCL(k,l+1))
        end do
      end do

    else if (i_lw_cld.eq.2) then
      ! consider optical thickness of clouds

      ! Upward LWR fluxes
      fcl_up(1) = flwr_up_sur
      do k=2, llwr
        fcl_up(k) = BSB(k)+(fcl_up(1)-BSB(1))*DCL(k,1)
        do l=1, k-1
          fcl_up(k) = fcl_up(k)-(BSB(l+1)-BSB(l))*0.5_wp*(DCL(k,l+1)+DCL(k,l))
        enddo
      enddo

      ! Downward LWR fluxes  
      fcl_dw(llwr) = 0._wp
      do k=llwr-1, 1, -1
        fcl_dw(k) = BSB(k)-BSB(llwr)*DCL(k,llwr)
        do l=k, llwr-1
          fcl_dw(k) = fcl_dw(k)+(BSB(l+1)-BSB(l))*0.5_wp*(DCL(k,l)+DCL(k,l+1))
        enddo
      enddo

    endif

    return

  end subroutine lwr_clouds

end module lw_radiation_mod


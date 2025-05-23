!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : l n d _ d e f
!
!  Purpose : definition of land model class 
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
module lnd_def

    use precision, only : wp
    use lnd_grid, only : nl_l 

    implicit none

    type lnd_0d_class
      real(wp) :: co2, c13_c12_atm, c14_c_atm
      real(wp) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd
      real(wp) :: Cflx_avg
      real(wp) :: Cflx_burial, C13flx_burial, C14flx_burial
      real(wp) :: ch4_emis
      real(wp) :: landc, landc13, landc14
      real(wp) :: burc, burc13, burc14
      real(wp) :: weath_scale
      real(wp) :: weath_carb_avg, weath_sil_avg
    end type 

    type lnd_2d_class 

     integer :: mask_lnd
     real(wp) :: Cflx_atm_lnd, C13flx_atm_lnd, C14flx_atm_lnd
     real(wp) :: f_land, f_land0, f_ice, f_ice_old, f_ice_grd, f_ice_grd_old, f_ice_nbr, f_shelf, f_shelf_old, f_lake, f_lake_old, f_veg, f_veg_old
     real(wp) :: z_veg_std, z_veg, z_veg_min, z_veg_max
     real(wp) :: f_crop, f_pasture
     real(wp) :: t2m_min_mon
     real(wp) :: t2m_ann_mean
     real(wp) :: infiltration, w_table, w_table_peat
     real(wp) :: f_wet, f_wet_cum, f_wet_max, f_wetland, w_table_cum, w_table_min, cti_lim
     real(wp) :: t_skin_veg, flx_g_veg, dflxg_dT_veg, flx_melt_veg, t_2m
     real(wp) :: f_lake_ice
     real(wp) :: h_lake
     real(wp) :: h_lake_conv
     real(wp) :: h_lake_mix
     real(wp) :: lake_water_tendency
     real(wp) :: gdd5, gdd5_temp, npp_real, npp13_real, npp14_real
     real(wp) :: veg_c_above, veg_c13_above, veg_c14_above
     real(wp) :: theta_fire_cum
     real(wp) :: alt
     real(wp) :: k_litter_peat, k_acro, k_litter_peat_anox, k_acro_anox, f_oxic_peat
     real(wp) :: litter_c_peat, acro_c, litter_c13_peat, acro_c13, litter_c14_peat, acro_c14
     real(wp) :: f_peat, f_peat_pot, acro_h, cato_h, peat_c_ini_year, dCpeat_dt
     real(wp) :: ch4_emis_wetland, ch4_emis_shelf, ch4_emis_peat, ch4_emis_lake
     real(wp) :: c13h4_emis_wetland, c13h4_emis_shelf, c13h4_emis_peat, c13h4_emis_lake
     real(wp) :: dust_emis_d, dust_emis_g, dust_emis_s, dust_emis
     real(wp) :: dust_dep
     real(wp) :: runoff_ann
     real(wp) :: f_carb
     real(wp) :: weath_carb, weath_sil, weath_loess
     real(wp) :: weath13_carb, weath13_sil
     real(wp) :: weath14_carb, weath14_sil
     real(wp) :: poc_export, poc13_export, poc14_export
     real(wp) :: doc_export, doc13_export, doc14_export
     real(wp) :: energy_cons_soil, energy_cons_ice, energy_cons_shelf, energy_cons_lake
     real(wp) :: carbon_cons_veg, carbon13_cons_veg, carbon14_cons_veg

     real(wp), allocatable, dimension(:) :: frac_surf
     real(wp), allocatable, dimension(:) :: disturbance
     real(wp), allocatable, dimension(:) :: coszm, daylength
     real(wp), allocatable, dimension(:) :: tatm, t2m, qatm, q2m, lwdown, swnet, swnet_min
     real(wp), allocatable, dimension(:) :: rain, snow
     real(wp), allocatable, dimension(:) :: wind
     real(wp), allocatable, dimension(:) :: pressure
     real(wp), allocatable, dimension(:) :: rough_m, rough_h, Ch, z0m, Ri
     real(wp), allocatable, dimension(:) :: r_a, r_s, beta_s
     real(wp), allocatable, dimension(:) :: r_a_can, r_s_can, beta_s_can
     real(wp), allocatable, dimension(:) :: alb_snow_vis_dir, alb_snow_vis_dif, alb_snow_nir_dir, alb_snow_nir_dif
     real(wp), allocatable, dimension(:) :: albedo, alb_vis_dir, alb_vis_dif, alb_nir_dir, alb_nir_dif
     integer,  allocatable, dimension(:) :: mask_snow
     real(wp), allocatable, dimension(:) :: f_snow, h_snow, w_snow, w_snow_max, w_snow_old, snowmelt, icemelt, snow_grain, dust_con
     real(wp), allocatable, dimension(:) :: runoff, runoff_sur, calving, drainage, water_cons
     real(wp), allocatable, dimension(:) :: rain_ground, evap_can, snow_ground, subl_can
     real(wp), allocatable, dimension(:) :: w_can, w_can_old, s_can, s_can_old, f_wat_can, f_snow_can
     real(wp), allocatable, dimension(:) :: transpiration, evap_surface, et
     real(wp), allocatable, dimension(:) :: f_wet_mon, w_table_mon
     real(wp), allocatable, dimension(:) :: f_wet_long
     real(wp), allocatable, dimension(:) :: flx_sh, flx_lh, flx_g, dflxg_dT, flx_melt, flx_lwu, lwnet
     real(wp), allocatable, dimension(:) :: t_skin, t_skin_old, t_skin_amp
     real(wp), allocatable, dimension(:) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
     real(wp), allocatable, dimension(:) :: f_sh, f_e, f_t, f_le, f_lt, f_lw, lh_ecan, qsat_e, dqsatdT_e, qsat_t, dqsatdT_t
     real(wp), allocatable, dimension(:) :: ci, g_can, gpp, npp, npp13, npp14, aresp
     real(wp), allocatable, dimension(:) :: discrimination
     real(wp), allocatable, dimension(:) :: lai, sai, phen, phen_acc, gdd, gamma_leaf, lambda, lai_bal
     real(wp), allocatable, dimension(:) :: npp_cum, npp13_cum, npp14_cum
     real(wp), allocatable, dimension(:) :: npp_ann, npp13_ann, npp14_ann
     real(wp), allocatable, dimension(:) :: veg_c, veg_h, pft_frac, seed_frac
     real(wp), allocatable, dimension(:) :: veg_c_below, veg_c13_below, veg_c14_below
     real(wp), allocatable, dimension(:) :: leaf_c, stem_c, root_c, veg_c13, veg_c14
     real(wp), allocatable, dimension(:) :: gamma_luc, gamma_ice, gamma_dist, gamma_dist_cum
     real(wp), allocatable, dimension(:) :: lambda_soil, lambda_int_soil, cap_soil
     real(wp), allocatable, dimension(:) :: lambda_ice, lambda_int_ice, cap_ice
     real(wp), allocatable, dimension(:) :: lambda_lake, lambda_int_lake, cap_lake 
     real(wp), allocatable, dimension(:) :: lambda_sublake, lambda_int_sublake, cap_sublake
     real(wp), allocatable, dimension(:) :: cap_shelf
     real(wp), allocatable, dimension(:) :: lambda_int_shelf
     real(wp), allocatable, dimension(:) :: lambda_s, lambda_dry
     real(wp), allocatable, dimension(:) :: kappa_int
     real(wp), allocatable, dimension(:) :: theta_sat, k_sat, psi_sat, theta_field, theta_wilt
     real(wp), allocatable, dimension(:) :: t_soil, t_soil_old, t_soil_max
     real(wp), allocatable, dimension(:) :: t_ice, t_ice_old
     real(wp), allocatable, dimension(:) :: t_shelf, t_shelf_old, t_shelf_max
     real(wp), allocatable, dimension(:) :: t_lake, t_lake_old
     real(wp), allocatable, dimension(:) :: t_sublake
     real(wp), allocatable, dimension(:) :: theta_w, theta_i, theta, w_w, w_i, w_w_old, w_i_old, w_w_phase, w_i_phase
     real(wp), allocatable, dimension(:) :: theta_w_shelf, theta_i_shelf, w_w_shelf, w_i_shelf
     real(wp), allocatable, dimension(:) :: theta_w_sublake, theta_i_sublake, w_w_sublake, w_i_sublake
     real(wp), allocatable, dimension(:) :: w_w_lake, w_i_lake, f_i_lake
     real(wp), allocatable, dimension(:) :: t_soil_cum, theta_w_cum, theta_i_cum
     real(wp), allocatable, dimension(:) :: t_shelf_cum, theta_w_shelf_cum, theta_i_shelf_cum
     real(wp), allocatable, dimension(:) :: t_sublake_cum, theta_w_sublake_cum, theta_i_sublake_cum
     real(wp), allocatable, dimension(:) :: psi
     integer,  allocatable, dimension(:) :: k_exp, psi_exp
     real(wp), allocatable, dimension(:) :: ftemp, fmoist, fdepth
     real(wp), allocatable, dimension(:) :: k_litter, k_fast, k_slow, k_litter_wet, k_fast_wet, k_slow_wet, diff_soilc, adv_soilc
     real(wp), allocatable, dimension(:) :: k_litter_shelf, k_fast_shelf, k_slow_shelf, diff_shelfc, adv_shelfc
     real(wp), allocatable, dimension(:) :: k_litter_lake, k_fast_lake, k_slow_lake, diff_lakec, adv_lakec
     real(wp), allocatable, dimension(:) :: k_litter_ice, k_fast_ice, k_slow_ice, diff_icec, adv_icec
     real(wp), allocatable, dimension(:) :: k_cato, ch4_frac_wet, ch4_frac_peat, ch4_frac_shelf, ch4_frac_lake
     real(wp), allocatable, dimension(:) :: frac_soc
     real(wp), allocatable, dimension(:) :: litter_c, fast_c, slow_c, litter_c13, fast_c13, slow_c13, litter_c14, fast_c14, slow_c14
     real(wp), allocatable, dimension(:) :: litter_c_shelf, fast_c_shelf, slow_c_shelf, litter_c13_shelf, fast_c13_shelf, slow_c13_shelf, litter_c14_shelf, fast_c14_shelf, slow_c14_shelf
     real(wp), allocatable, dimension(:) :: litter_c_lake, fast_c_lake, slow_c_lake, litter_c13_lake, fast_c13_lake, slow_c13_lake, litter_c14_lake, fast_c14_lake, slow_c14_lake
     real(wp), allocatable, dimension(:) :: litter_c_ice, fast_c_ice, slow_c_ice, litter_c13_ice, fast_c13_ice, slow_c13_ice, litter_c14_ice, fast_c14_ice, slow_c14_ice
     real(wp), allocatable, dimension(:) :: cato_c, cato_c13, cato_c14
     real(wp), allocatable, dimension(:) :: soil_c_tot, soil_resp, soil_c13_tot, soil_resp13, soil_c14_tot, soil_resp14
     real(wp), allocatable, dimension(:,:) :: litterfall, litterfall13, litterfall14
     real(wp), allocatable, dimension(:) :: litter_in_frac
     real(wp), allocatable, dimension(:,:) :: wilt
     real(wp), allocatable, dimension(:,:) :: root_frac
     real(wp), allocatable, dimension(:,:) :: soil_resp_l

     real(wp), allocatable, dimension(:) :: lithology_gemco2
     real(wp), allocatable, dimension(:) :: lithology_uhh
     real(wp), allocatable, dimension(:) :: lithology_shelf_uhh

     real(wp), allocatable, dimension(:) :: energy_cons_surf1, energy_cons_surf2
     real(wp), allocatable, dimension(:) :: carbon_cons_soil, carbon13_cons_soil, carbon14_cons_soil

    end type 

    type lnd_class
      integer :: ncells
      integer, allocatable :: id_map(:,:)
      integer, allocatable :: ij_1d(:,:)
      real(wp), dimension(nl_l) :: z_lake
      type(lnd_2d_class), allocatable :: l2d(:,:)
      type(lnd_0d_class) :: l0d
    end type

    private
    public :: lnd_class, lnd_0d_class, lnd_2d_class

end module

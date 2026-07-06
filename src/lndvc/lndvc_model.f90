module lndvc_model
    ! Driver / orchestrator for the virtual-cell (lndvc) framework.
    !
    ! Owns the per-domain control (init/update/end), the loop over leaf virtual
    ! cells, dispatch to single-column physics by surface class, and the
    ! reduction back to coarse-cell state for coupling. Physics is ported from
    ! lnd/smb as single-column, single-class modules and called from here.
    ! See docs/design/virtual-cells.md.
    !
    ! Phase 1: allocation + decomposition + aggregation plumbing, wired into
    ! climber.x behind flag_lndvc (default off). Physics dispatch is a no-op
    ! until the port (Phase 2); numerical identity with the reference land model
    ! is therefore gated on Phase 2.

    use ncio
    use precision, only : sp, dp, wp

    use lndvc_def
    use lndvc_grid
    use lndvc_decomp,    only : lndvc_decompose, lndvc_remap_state
    use lndvc_downscale, only : lndvc_downscale_forcing
    use lndvc_aggregate, only : lndvc_aggregate_weights, lndvc_aggregate_cell, lndvc_conservation_check

    ! ported single-column physics
    use semi_m,          only : semi
    use smb_par_m,       only : p0, h_atm, gamma, prc_par, surf_par
    ! ported single-column lake physics
    use lndvc_surface_par_lnd, only : resist_aer_lake, snow_albedo_lake, surface_albedo_lake, resist_sur_lake
    use lndvc_lake_par_mod,    only : lake_par_thermal
    use lndvc_ebal_lake_mod,   only : ebal_lake, update_tskin_lake
    use lndvc_lake_temp_mod,   only : lake_temp
    use lndvc_hydrology_mod,   only : surface_hydrology_lake
    use lndvc_init_cell_mod,   only : lndvc_init_cell_veg
    use lnd_params,            only : lnd_surf_par => surf_par
    use lnd_params,            only : soil_par, hydro_par
    use wiso_params,           only : l_wiso, nwiso, i_o18, Rstd

    implicit none

    private
    public :: lndvc_init
    public :: lndvc_update
    public :: lndvc_init_land
    public :: lndvc_end

contains

    ! === Whole-domain control =================================================

    subroutine lndvc_init(lnd, nx, ny)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        integer,           intent(in)    :: nx, ny

        call lndvc_grid_init()
        call lndvc_alloc(lnd, nx, ny)

        return
    end subroutine lndvc_init

    subroutine lndvc_update(lnd)
        ! One coupling step: update every leaf virtual cell, then aggregate to
        ! coarse-cell state for coupling. The decomposition is assumed current
        ! (refreshed by the coupler via lndvc_decompose before this call).

        implicit none

        type(lndvc_class), intent(inout) :: lnd

        integer :: n, i, j, k
        real(wp), allocatable :: wt(:)

        allocate(wt(lnd%n_vc))

        ! TODO: OMP over the compressed active-leaf list (lnd%leaf_1d).
        do n = 1, lnd%ncells
            i = lnd%ij_1d(1,n)
            j = lnd%ij_1d(2,n)
            do k = 1, lnd%n_vc
                if (lnd%vc(i,j,k)%desc%class == 0) cycle
                call lndvc_update_vc(lnd%vc(i,j,k))
            end do
            call lndvc_aggregate_weights(lnd%vc(i,j,:), wt)
            call lndvc_aggregate_cell(lnd%cell(i,j), lnd%vc(i,j,:), wt)
        end do

        deallocate(wt)

        call lndvc_conservation_check(lnd)

        return

    end subroutine lndvc_update

    subroutine lndvc_end(lnd)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        call lndvc_dealloc(lnd)
        return
    end subroutine lndvc_end

    ! === Allocation ===========================================================

    subroutine lndvc_alloc(lnd, nx, ny)
        implicit none
        type(lndvc_class), intent(inout) :: lnd
        integer,           intent(in)    :: nx, ny

        lnd%n_vc   = 4          ! identity baseline: land, lake, ice (+shelf)
        lnd%ncells = 0

        if (allocated(lnd%vc))     deallocate(lnd%vc)
        if (allocated(lnd%cell))   deallocate(lnd%cell)
        if (allocated(lnd%id_map)) deallocate(lnd%id_map)
        if (allocated(lnd%ij_1d))  deallocate(lnd%ij_1d)

        allocate(lnd%vc(nx,ny,lnd%n_vc))
        allocate(lnd%cell(nx,ny))
        allocate(lnd%id_map(nx,ny))
        allocate(lnd%ij_1d(2,nx*ny))

        lnd%id_map = 0
        lnd%ij_1d  = 0

        return
    end subroutine lndvc_alloc

    subroutine lndvc_dealloc(lnd)
        implicit none
        type(lndvc_class), intent(inout) :: lnd

        if (allocated(lnd%vc))      deallocate(lnd%vc)
        if (allocated(lnd%cell))    deallocate(lnd%cell)
        if (allocated(lnd%id_map))  deallocate(lnd%id_map)
        if (allocated(lnd%ij_1d))   deallocate(lnd%ij_1d)
        if (allocated(lnd%z_edges)) deallocate(lnd%z_edges)
        if (allocated(lnd%leaf_1d)) deallocate(lnd%leaf_1d)
        if (allocated(lnd%z_lake))  deallocate(lnd%z_lake)

        return
    end subroutine lndvc_dealloc

    ! === Per-vc dispatch ======================================================

    subroutine lndvc_update_vc(vc)
        ! Update a single virtual cell: downscale forcing to its elevation, then
        ! dispatch to the single-column physics for its surface class.

        implicit none

        type(vc_t), intent(inout) :: vc

        call lndvc_downscale_forcing(vc%forc, vc%desc)

        select case(vc%desc%class)
            case(1)   ! land
                call lndvc_update_land(vc)
            case(2)   ! lake
                call lndvc_update_lake(vc)
            case(3)   ! ice
                call lndvc_update_ice(vc)
        end select

        return

    end subroutine lndvc_update_vc

    subroutine lndvc_init_land(vc, m_theta_sat, m_k_sat, m_psi_sat, m_Bi, m_lambda_s, m_lambda_dry, &
                               c13_c12_atm, c14_c_atm)
        ! One-time physical init of a vegetated land vc (port C decision B):
        ! seed the soil parameters from the coarse-cell mineral texture (the
        ! reference lnd-init block, src/lnd/lnd_model.f90 ~1206) then run the
        ! ported init_cell_veg state init. frac_surf is set from pft_frac (bare
        ! is the remainder) since the multi-class surface_frac_up does not apply
        ! to a single-class veg vc. Called once per land vc from cmn_to_lndvc.

        implicit none

        type(vc_t), intent(inout) :: vc
        real(wp), intent(in) :: m_theta_sat, m_k_sat, m_psi_sat, m_Bi, m_lambda_s, m_lambda_dry
        real(wp), intent(in) :: c13_c12_atm, c14_c_atm

        integer  :: k
        real(wp) :: sum_pft

        ! --- soil parameters from mineral texture -----------------------------
        do k = 1, nl
            vc%soil%theta_sat(k) = m_theta_sat
            vc%soil%psi_sat(k)   = m_psi_sat
            vc%soil%k_sat(k)     = m_k_sat
            vc%soil%k_exp(k)     = 2*nint(m_Bi)+3
            vc%soil%psi_exp(k)   = -nint(m_Bi)
            vc%soil%theta_field(k) = vc%soil%theta_sat(k) * (0.1_wp / (86400._wp * vc%soil%k_sat(k)))**(1._wp/vc%soil%k_exp(k))
            vc%soil%theta_wilt(k)  = vc%soil%theta_sat(k) * (hydro_par%p_psi_min / vc%soil%psi_sat(k))**(1._wp/vc%soil%psi_exp(k))
            vc%soil%lambda_s(k)   = m_lambda_s
            vc%soil%lambda_dry(k) = m_lambda_dry
        enddo
        if (soil_par%uniform_porosity) vc%soil%theta_sat(:) = soil_par%theta_sat_u
        if (soil_par%uniform_soil_par_therm) then
            vc%soil%lambda_s(:)   = soil_par%lambda_s_u
            vc%soil%lambda_dry(:) = soil_par%lambda_dry_u
        endif

        ! --- physical veg/soil state (ported init_cell_veg) -------------------
        call lndvc_init_cell_veg(c13_c12_atm, c14_c_atm, vc%desc%w, &
            vc%veg%pft_frac, vc%flx%t_skin, vc%soil%t_soil, vc%flx%w_can, vc%flx%s_can, &
            vc%snow%w_snow, vc%snow%w_snow_max, vc%snow%h_snow, vc%snow%mask_snow, &
            vc%soil%theta, vc%soil%theta_w, vc%soil%theta_i, vc%soil%w_w, vc%soil%w_i, vc%soil%theta_sat, &
            vc%soil%alt, vc%veg%gdd5, vc%veg%gdd, vc%veg%phen, vc%veg%phen_acc, vc%veg%lai_bal, &
            vc%veg%lai, vc%veg%sai, vc%soil%root_frac, vc%carb%litter_in_frac, &
            vc%veg%gamma_dist, vc%veg%gamma_fire, vc%veg%npp_ann, vc%veg%npp13_ann, vc%veg%npp14_ann, &
            vc%veg%leaf_c, vc%veg%root_c, vc%veg%stem_c, vc%veg%veg_h, vc%veg%veg_c, vc%veg%veg_c13, vc%veg%veg_c14, &
            vc%carb%f_peat, vc%carb%f_peat_pot, vc%soil%w_table_min, vc%soil%w_table_peat, vc%carb%dCpeat_dt, &
            vc%flx%z0m)

        ! --- sub-tile fractions within the vc (bare = remainder) --------------
        sum_pft = 0._wp
        do k = 1, npft
            vc%flx%frac_surf(k) = vc%veg%pft_frac(k)
            sum_pft = sum_pft + vc%veg%pft_frac(k)
        enddo
        vc%flx%frac_surf(i_bare) = max(0._wp, 1._wp - sum_pft)

        return

    end subroutine lndvc_init_land

    subroutine lndvc_update_land(vc)
        implicit none
        type(vc_t), intent(inout) :: vc
        ! TODO: ported single-column land physics (veg energy balance, soil
        !       thermal+permafrost, hydrology, snow, soil carbon).
        return
    end subroutine lndvc_update_land

    subroutine lndvc_update_lake(vc)
        ! Single-column lake-surface update (pattern A): unpack the vc's lake /
        ! snow / flux blocks into the ported single-tile lake chain, mirroring the
        ! reference lnd lake path (surface params -> lake thermal params -> energy
        ! balance -> lake temperature -> skin update -> surface hydrology). The
        ! downstream sublake soil column (sublake_par_thermal/sublake_temp) is
        ! deferred: it needs soil params not on the lake block and feeds only
        ! sublake-soil / lake-carbon state, not the surface coupling currency.
        ! t_skin_old is snapshotted inside ebal_lake; the melt-iso redistribution
        ! that lives in the reference lnd_model wrapper is deferred (like SMB iso).

        implicit none

        type(vc_t), intent(inout) :: vc

        integer, parameter :: ii = 0, jj = 0   ! debug-print indices only
        real(wp) :: z0m_lake
        real(wp) :: calving, runoff_sur              ! lake calving/runoff not yet aggregated
        real(wp) :: calving_iso(nwiso), runoff_sur_iso(nwiso)

        ! carry-over snapshot for the step
        vc%snow%w_snow_old = vc%snow%w_snow

        ! no canopy over lake: throughfall = precipitation (tag iso at VSMOW)
        vc%flx%rain_ground(1) = vc%forc%rain
        vc%flx%snow_ground(1) = vc%forc%snow
        if (l_wiso) then
            vc%flx%rain_ground_iso(1,i_o18) = Rstd(i_o18) * vc%forc%rain
            vc%flx%snow_ground_iso(1,i_o18) = Rstd(i_o18) * vc%forc%snow
        endif

        ! lake surface momentum roughness (fixed lake/ice value)
        z0m_lake = lnd_surf_par%z0m_lake_ice
        vc%flx%z0m(1) = z0m_lake

        ! --- surface parameters -----------------------------------------------
        call resist_aer_lake(vc%snow%h_snow, vc%forc%tatm, vc%flx%t_skin(1), vc%forc%wind, &
            z0m_lake, vc%flx%rough_m(1), vc%flx%rough_h(1), vc%flx%Ch(1), vc%flx%r_a(1), vc%flx%Ri(1))

        call snow_albedo_lake(vc%flx%t_skin(1), vc%forc%snow, vc%snow%w_snow, vc%snow%w_snow_max, &
            vc%forc%dust, vc%forc%coszm, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_nir_dif, &
            vc%snow%snow_grain, vc%snow%dust_con)

        call surface_albedo_lake(vc%snow%h_snow, vc%forc%coszm, vc%lake%f_lake_ice, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_nir_dif, &
            vc%snow%f_snow, vc%flx%alb_vis_dir(1), vc%flx%alb_vis_dif(1), vc%flx%alb_nir_dir(1), vc%flx%alb_nir_dif(1), &
            vc%flx%albedo(1))

        call resist_sur_lake(vc%flx%beta_s(1), vc%flx%r_s(1))

        ! --- lake thermal parameters ------------------------------------------
        call lake_par_thermal(vc%snow%h_snow, vc%lake%h_lake, vc%lake%t_lake(1:nl_l), vc%lake%f_i_lake, &
            vc%forc%wind, vc%desc%lat, &
            vc%lake%cap_lake, vc%lake%lambda_lake, vc%lake%lambda_int_lake)

        ! --- surface energy balance -------------------------------------------
        call ebal_lake(vc%snow%mask_snow, vc%snow%h_snow, vc%lake%lambda_lake, &
            vc%flx%t_skin(1), vc%flx%t_skin_old(1), vc%lake%t_lake, vc%forc%tatm, vc%forc%qatm, vc%forc%pressure, &
            vc%forc%swnet, vc%forc%lwdown, &
            vc%flx%beta_s(1), vc%flx%r_s(1), vc%flx%r_a(1), &
            vc%flx%flx_g(1), vc%flx%dflxg_dT(1), vc%flx%flx_melt(1), vc%flx%t_skin_amp(1), &
            vc%flx%num_lh(1), vc%flx%num_sh(1), vc%flx%num_sw(1), vc%flx%num_lw(1), &
            vc%flx%denom_lh(1), vc%flx%denom_sh(1), vc%flx%denom_lw(1), &
            vc%flx%f_sh(1), vc%flx%f_e(1), vc%flx%f_le(1), vc%flx%f_lw(1), vc%flx%qsat_e(1), vc%flx%dqsatdT_e(1), &
            vc%energy_cons_surf1, ii, jj)

        ! --- lake temperature (convection handled inside lake_temp) ------------
        call lake_temp(vc%snow%mask_snow, vc%snow%h_snow, vc%lake%h_lake, vc%lake%cap_lake, vc%lake%lambda_int_lake, &
            vc%flx%flx_g(1), vc%flx%dflxg_dT(1), vc%flx%flx_melt(1), &
            vc%lake%t_lake, vc%snow%w_snow, vc%lake%w_w_lake, vc%lake%w_i_lake, vc%lake%f_i_lake, vc%lake%f_lake_ice, &
            vc%snow%snowmelt, vc%lake%t_lake_old, vc%snow%w_snow_old, &
            vc%lake%h_lake_conv, vc%lake%h_lake_mix, &
            vc%lake%energy_cons_lake, ii, jj)

        ! --- skin temperature + surface fluxes --------------------------------
        call update_tskin_lake(vc%snow%mask_snow, vc%flx%t_skin_old(1), vc%flx%dflxg_dT(1), &
            vc%forc%tatm, vc%forc%qatm, vc%forc%swnet, vc%forc%lwdown, &
            vc%lake%t_lake, vc%lake%t_lake_old, vc%flx%flx_g(1), vc%flx%flx_melt(1), &
            vc%flx%t_skin(1), vc%flx%flx_sh(1), vc%flx%flx_lwu(1), vc%flx%flx_lh(1), vc%flx%evap_surface(1), vc%flx%et(1), &
            vc%flx%num_lh(1), vc%flx%num_sh(1), vc%flx%num_sw(1), vc%flx%num_lw(1), &
            vc%flx%denom_lh(1), vc%flx%denom_sh(1), vc%flx%denom_lw(1), &
            vc%flx%f_sh(1), vc%flx%f_e(1), vc%flx%f_le(1), vc%flx%f_lw(1), vc%flx%qsat_e(1), vc%flx%dqsatdT_e(1), &
            vc%energy_cons_surf2, ii, jj, &
            vc%snow%w_snow, vc%snow%w_snow_iso, vc%lake%w_w_lake(1), vc%lake%w_w_lake_iso(1,:), &
            vc%flx%evap_surface_iso(1,:), vc%flx%et_iso(1,:))

        ! --- surface hydrology (snow budget + lake water balance -> runoff) ----
        call surface_hydrology_lake(vc%snow%mask_snow, vc%flx%evap_surface(1), vc%flx%snow_ground(1), vc%flx%rain_ground(1), &
            vc%snow%snowmelt, vc%lake%cap_lake(1), vc%lake%t_lake(1), vc%snow%w_snow_old, vc%snow%w_snow, vc%snow%w_snow_max, &
            vc%snow%h_snow, calving, runoff_sur, vc%lake%lake_water_tendency, &
            vc%flx%evap_surface_iso(1,:), vc%flx%snow_ground_iso(1,:), vc%flx%rain_ground_iso(1,:), &
            vc%snow%snowmelt_iso, vc%snow%w_snow_iso, calving_iso, runoff_sur_iso)

        return
    end subroutine lndvc_update_lake

    subroutine lndvc_update_ice(vc)
        ! Single-column ice-surface update via the ported SEMI driver (pattern A):
        ! the vc's per-class blocks are unpacked into SEMI's scalar signature.
        ! Genuine carry-over state binds to block fields (snow/ice/flx); pure
        ! per-call diagnostics are call-local (kept lean, see design sec 12).
        ! SEMI does its own elevation downscaling from forc%z_sur_i to desc%z.

        implicit none

        type(vc_t), intent(inout) :: vc

        integer  :: ii, jj
        real(wp) :: dz, dT
        real(wp) :: f_ele, pressure

        ! downscaled/diagnostic scalars SEMI writes but we do not persist
        real(wp) :: tam, t2m, q2m, qsat, dqsatdT, r_a
        real(wp) :: u700, v700, wind, snow_rate, rain_rate, prc, f_wind
        real(wp) :: alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif
        real(wp) :: alb_bg, cld, swnet, swnet_min, swdown, lwdown
        real(wp) :: dflxg_dT, flx_melt, flx_lwu
        real(wp) :: num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw
        real(wp) :: f_sh, f_e, f_lh, f_lw

        ii = 0; jj = 0   ! SEMI uses these only for debug prints

        ! --- geometry-derived forcing (framework side) ------------------------
        ! barometric surface pressure at the vc elevation
        pressure = p0 * exp(-vc%desc%z / h_atm)
        vc%forc%pressure = pressure
        ! precipitation elevation-correction factor (Clausius-Clapeyron), cf. topo_factors
        if (prc_par%l_elevation_corr .and. vc%desc%z >= prc_par%z_sur_crit_fele) then
            dz = vc%desc%z - vc%forc%z_sur_i
            dT = -gamma*dz
            f_ele = exp(prc_par%dP_dT*dT)
        else
            f_ele = 1._wp
        end if
        if (vc%desc%z >= prc_par%z_sur_high_fele) f_ele = 0.1_wp*f_ele
        vc%desc%f_ele = f_ele
        ! pure ice vc: fully ice-covered patch, firn background albedo
        vc%forc%f_ice   = 1._wp
        vc%forc%alb_ice = surf_par%alb_firn

        ! --- surface energy + mass balance (ported SEMI, pattern A) ------------
        call semi(ii, jj, vc%forc%f_ice, vc%forc%alb_ice, &
            vc%desc%z, vc%forc%z_sur_i, vc%desc%z_sur_std, &
            vc%desc%dz_dx, vc%desc%dz_dy, vc%desc%dz_sur, vc%desc%f_ele, &
            vc%forc%tam_i, vc%forc%t2m_bias_i, vc%forc%dTvar, vc%forc%gam_i, vc%forc%tstd_i, vc%forc%ram_i, &
            pressure, vc%forc%u700_i, vc%forc%v700_i, vc%forc%wind_i, vc%forc%prc_i, vc%forc%prc_bias_i, &
            vc%forc%alb_vis_dir_i, vc%forc%alb_nir_dir_i, vc%forc%alb_vis_dif_i, vc%forc%alb_nir_dif_i, &
            vc%forc%swd_sur_vis_dir_i, vc%forc%swd_sur_nir_dir_i, vc%forc%swd_sur_vis_dif_i, vc%forc%swd_sur_nir_dif_i, &
            vc%forc%dswd_dalb_vis_dir_i, vc%forc%dswd_dalb_nir_dir_i, vc%forc%dswd_dalb_vis_dif_i, vc%forc%dswd_dalb_nir_dif_i, &
            vc%forc%dswd_dz_nir_dir_i, vc%forc%dswd_dz_nir_dif_i, vc%forc%dust_i, vc%forc%coszm_i, &
            vc%forc%swd_toa_i, vc%forc%swd_toa_min_i, vc%forc%cld_i, vc%forc%lwdown_i, vc%forc%gam_lw_i, &
            tam, t2m, vc%flx%t_skin(1), vc%flx%t_skin_old(1), vc%flx%t_skin_amp(1), &
            vc%ice%t_prof, vc%ice%t_prof_old, &
            q2m, qsat, dqsatdT, r_a, &
            u700, v700, wind, snow_rate, rain_rate, prc, f_wind, &
            vc%snow%mask_snow, vc%snow%f_snow, vc%snow%h_snow, vc%snow%w_snow, vc%snow%w_snow_old, vc%snow%w_snow_max, &
            vc%snow%snow_grain, vc%snow%dust_con, &
            alb_vis_dir, alb_nir_dir, alb_vis_dif, alb_nir_dif, &
            vc%snow%alb_snow_vis_dir, vc%snow%alb_snow_nir_dir, vc%snow%alb_snow_vis_dif, vc%snow%alb_snow_nir_dif, &
            vc%snow%dt_snowfree, alb_bg, vc%flx%albedo(1), cld, swnet, swnet_min, swdown, lwdown, &
            vc%flx%flx_g(1), dflxg_dT, flx_melt, vc%flx%flx_sh(1), flx_lwu, vc%flx%flx_lh(1), vc%flx%evap_surface(1), &
            num_lh, num_sh, num_sw, num_lw, denom_lh, denom_sh, denom_lw, &
            f_sh, f_e, f_lh, f_lw, &
            vc%snow%snowmelt, vc%snow%icemelt, vc%snow%refreezing, vc%snow%refreezing_sum, vc%ice%f_rfz_to_snow)

        ! --- SMB mass-budget diagnostic (kg/m2/s) -----------------------------
        vc%ice%melt   = vc%snow%snowmelt + vc%snow%icemelt
        vc%ice%runoff = vc%snow%snowmelt + vc%snow%icemelt + rain_rate - vc%snow%refreezing
        vc%ice%smb    = snow_rate - vc%flx%evap_surface(1) - vc%snow%snowmelt - vc%snow%icemelt + vc%snow%refreezing

        return
    end subroutine lndvc_update_ice

end module lndvc_model

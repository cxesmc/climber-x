module lndvc_decomp
    ! Decomposition service for the virtual-cell framework.
    !
    ! Turns each coarse coupler cell into a list of leaf virtual cells via a
    ! composable pipeline of refine stages (see docs/design/virtual-cells.md, sec 4.3):
    !
    !   coarse cell --[optional bilinear regrid]--> mid-res cell
    !               --[hypsometric split]--------> virtual cell {z,dz,w,slope,class}
    !
    ! Also owns the conservative state-remap that runs when the decomposition
    ! changes (ice advance/retreat, topography shift).
    !
    ! Physics is grid-agnostic: all horizontal coupling lives here, never in the
    ! column physics.
    !
    ! Phase 1 (identity baseline): a single elevation band per cell, with one
    ! leaf per present surface class at the cell-mean elevation. This reproduces
    ! today's surface-type tiling; the elevation dimension is trivial (1 band).

    use precision, only : wp
    use lndvc_def
    use lndvc_grid
    use const_m,    only : T0
    use timer,      only : nmon_year      ! months per year (cwd_mon accumulator size)
    use smb_grid_m, only : nl_smb => nl   ! snow+ice/firn layer count (=4; distinct from land nl=5)
    use wiso_params, only : nwiso         ! water-isotope tracer count (lake iso arrays)

    implicit none

    private
    public :: lndvc_decompose        ! build the leaf virtual-cell list for the domain
    public :: lndvc_remap_state       ! conservative state transfer when decomposition changes

contains

    subroutine lndvc_decompose(lnd, mask_lnd, f_veg, f_ice, f_lake, z_veg, z_ice)
        ! Build the leaf virtual-cell descriptors from coarse-cell geometry.
        ! Identity baseline: one leaf per present class at the cell-mean elevation.
        ! (Later: mid-res regrid + hypsometric split into multiple bands.)

        implicit none

        type(lndvc_class), intent(inout) :: lnd
        integer,  intent(in) :: mask_lnd(:,:)
        real(wp), intent(in) :: f_veg(:,:), f_ice(:,:), f_lake(:,:)
        real(wp), intent(in) :: z_veg(:,:), z_ice(:,:)

        integer :: i, j, k, nx, ny

        nx = size(f_veg,1)
        ny = size(f_veg,2)

        lnd%ncells = 0

        do j = 1, ny
        do i = 1, nx

            ! reset leaves for this cell
            do k = 1, lnd%n_vc
                lnd%vc(i,j,k)%desc%class = 0
                lnd%vc(i,j,k)%desc%w     = 0._wp
            end do
            lnd%id_map(i,j) = 0

            if (mask_lnd(i,j) /= 1) cycle

            lnd%ncells = lnd%ncells + 1
            lnd%ij_1d(1,lnd%ncells) = i
            lnd%ij_1d(2,lnd%ncells) = j
            lnd%id_map(i,j) = lnd%ncells

            k = 0
            if (f_veg(i,j)  > 0._wp) call set_leaf(lnd%vc(i,j,:), k, 1, z_veg(i,j), f_veg(i,j))
            if (f_lake(i,j) > 0._wp) call set_leaf(lnd%vc(i,j,:), k, 2, z_veg(i,j), f_lake(i,j))
            if (f_ice(i,j)  > 0._wp) call set_leaf(lnd%vc(i,j,:), k, 3, z_ice(i,j),  f_ice(i,j))

        end do
        end do

        return

    end subroutine lndvc_decompose

    subroutine set_leaf(vc, k, class, z, w)
        ! Populate leaf k+1 with a class/elevation/weight, allocate its class
        ! blocks, and advance k. (Inner arrays are sized during the physics port.)

        implicit none

        type(vc_t), intent(inout) :: vc(:)
        integer,    intent(inout) :: k
        integer,    intent(in)    :: class
        real(wp),   intent(in)    :: z, w

        k = k + 1
        if (k > size(vc)) return   ! guard: n_vc too small for present classes

        vc(k)%desc%class = class
        vc(k)%desc%z     = z
        vc(k)%desc%dz    = 0._wp
        vc(k)%desc%w     = w
        ! identity baseline: single band at cell-mean elevation, no sub-grid
        ! slope/variance/orography (populated once hypsometry + slopes are added)
        vc(k)%desc%dz_dx     = 0._wp
        vc(k)%desc%dz_dy     = 0._wp
        vc(k)%desc%grad      = 0._wp
        vc(k)%desc%z_sur_std = 0._wp
        vc(k)%desc%dz_sur    = 0._wp
        vc(k)%desc%f_ele     = 1._wp

        select case(class)
            case(1)   ! land
                if (.not. allocated(vc(k)%veg))  allocate(vc(k)%veg)
                if (.not. allocated(vc(k)%soil)) allocate(vc(k)%soil)
                if (.not. allocated(vc(k)%carb)) allocate(vc(k)%carb)
                if (.not. allocated(vc(k)%snow)) allocate(vc(k)%snow)
                call alloc_land_blocks(vc(k))
            case(2)   ! lake
                if (.not. allocated(vc(k)%lake)) allocate(vc(k)%lake)
                if (.not. allocated(vc(k)%snow)) allocate(vc(k)%snow)
                call alloc_lake_blocks(vc(k))
            case(3)   ! ice
                if (.not. allocated(vc(k)%ice))  allocate(vc(k)%ice)
                if (.not. allocated(vc(k)%snow)) allocate(vc(k)%snow)
                call alloc_ice_blocks(vc(k))
        end select

        return

    end subroutine set_leaf

    subroutine alloc_ice_blocks(vc)
        ! Size + initialize the inner arrays an ice virtual cell needs for the
        ! SEMI single-column SMB physics (pattern A). Only the fields SEMI reads
        ! or writes as carry-over state are allocated (kept lean, see
        ! docs/design/virtual-cells.md sec 12). Done once, guarded by allocation.

        implicit none

        type(vc_t), intent(inout) :: vc

        ! First-time setup only: t_prof allocation doubles as the once-guard, so
        ! prognostic carry-over state is not reset on every decompose.
        if (allocated(vc%ice%t_prof)) return

        ! snow+ice/firn thermal profile (0:nl, index 0 = skin layer)
        allocate(vc%ice%t_prof(0:nl_smb))
        allocate(vc%ice%t_prof_old(0:nl_smb))
        vc%ice%t_prof(:)     = T0
        vc%ice%t_prof_old(:) = T0
        vc%ice%smb    = 0._wp
        vc%ice%melt   = 0._wp
        vc%ice%runoff = 0._wp
        vc%ice%f_rfz_to_snow = 0._wp

        ! surface fluxes consumed by aggregation + t_skin carry-over (size 1)
        allocate(vc%flx%t_skin(1), vc%flx%t_skin_old(1), vc%flx%t_skin_amp(1))
        allocate(vc%flx%albedo(1))
        allocate(vc%flx%flx_sh(1), vc%flx%flx_lh(1), vc%flx%flx_g(1))
        allocate(vc%flx%evap_surface(1))
        vc%flx%t_skin(1)       = T0
        vc%flx%t_skin_old(1)   = T0
        vc%flx%t_skin_amp(1)   = 0._wp
        vc%flx%albedo(1)       = 0._wp
        vc%flx%flx_sh(1)       = 0._wp
        vc%flx%flx_lh(1)       = 0._wp
        vc%flx%flx_g(1)        = 0._wp
        vc%flx%evap_surface(1) = 0._wp

        ! snowpack carry-over scalars
        vc%snow%mask_snow      = 0
        vc%snow%f_snow         = 0._wp
        vc%snow%h_snow         = 0._wp
        vc%snow%w_snow         = 0._wp
        vc%snow%w_snow_old     = 0._wp
        vc%snow%w_snow_max     = 0._wp
        vc%snow%snow_grain     = 0._wp
        vc%snow%dust_con       = 0._wp
        vc%snow%refreezing     = 0._wp
        vc%snow%refreezing_sum = 0._wp
        vc%snow%dt_snowfree    = 0._wp
        vc%snow%alb_snow_vis_dir = 0._wp
        vc%snow%alb_snow_nir_dir = 0._wp
        vc%snow%alb_snow_vis_dif = 0._wp
        vc%snow%alb_snow_nir_dif = 0._wp

        return

    end subroutine alloc_ice_blocks

    subroutine alloc_lake_blocks(vc)
        ! Size + initialize the inner arrays a lake virtual cell needs for the
        ! single-column lake surface+thermal chain (pattern A). Array bounds
        ! mirror the reference lnd allocation (src/lnd/lnd_model.f90). Lean
        ! identity init (thermal profiles = T0, everything else 0); the physical
        ! lake spin-up (init_cell_lake: t_lake ramp, sublake theta from soil,
        ! cross-tile skin mean) is a follow-up. Done once, guarded by allocation.

        implicit none

        type(vc_t), intent(inout) :: vc

        ! t_lake allocation doubles as the once-guard
        if (allocated(vc%lake%t_lake)) return

        ! lake thermal column (0:nl_l)
        allocate(vc%lake%t_lake(0:nl_l), vc%lake%t_lake_old(0:nl_l))
        allocate(vc%lake%lambda_lake(0:nl_l), vc%lake%lambda_int_lake(0:nl_l), vc%lake%cap_lake(0:nl_l))
        allocate(vc%lake%w_w_lake(nl_l), vc%lake%w_i_lake(nl_l), vc%lake%f_i_lake(nl_l))
        allocate(vc%lake%w_w_lake_iso(nl_l,nwiso), vc%lake%w_i_lake_iso(nl_l,nwiso))
        vc%lake%t_lake(:)          = T0
        vc%lake%t_lake_old(:)      = T0
        vc%lake%lambda_lake(:)     = 0._wp
        vc%lake%lambda_int_lake(:) = 0._wp
        vc%lake%cap_lake(:)        = 0._wp
        vc%lake%w_w_lake(:)        = 0._wp
        vc%lake%w_i_lake(:)        = 0._wp
        vc%lake%f_i_lake(:)        = 0._wp
        vc%lake%w_w_lake_iso(:,:)  = 0._wp
        vc%lake%w_i_lake_iso(:,:)  = 0._wp

        ! sublake soil column (t_sublake 0:nl, rest 1:nl)
        allocate(vc%lake%t_sublake(0:nl))
        allocate(vc%lake%lambda_int_sublake(0:nl), vc%lake%cap_sublake(nl))
        allocate(vc%lake%theta_w_sublake(nl), vc%lake%theta_i_sublake(nl))
        allocate(vc%lake%w_w_sublake(nl), vc%lake%w_i_sublake(nl))
        allocate(vc%lake%t_sublake_cum(nl), vc%lake%theta_w_sublake_cum(nl), vc%lake%theta_i_sublake_cum(nl))
        vc%lake%t_sublake(:)           = T0
        vc%lake%lambda_int_sublake(:)  = 0._wp
        vc%lake%cap_sublake(:)         = 0._wp
        vc%lake%theta_w_sublake(:)     = 0._wp
        vc%lake%theta_i_sublake(:)     = 0._wp
        vc%lake%w_w_sublake(:)         = 0._wp
        vc%lake%w_i_sublake(:)         = 0._wp
        vc%lake%t_sublake_cum(:)       = 0._wp
        vc%lake%theta_w_sublake_cum(:) = 0._wp
        vc%lake%theta_i_sublake_cum(:) = 0._wp

        ! lake scalars
        vc%lake%h_lake              = 0._wp
        vc%lake%h_lake_conv         = 0._wp
        vc%lake%h_lake_mix          = 0._wp
        vc%lake%f_lake_ice          = 0._wp
        vc%lake%lake_water_tendency = 0._wp
        vc%lake%energy_cons_lake    = 0._wp

        ! full surface-flux block (single lake tile)
        call alloc_surface_flux(vc%flx, 1)

        ! snowpack carry-over scalars + iso arrays
        vc%snow%mask_snow      = 0
        vc%snow%f_snow         = 0._wp
        vc%snow%h_snow         = 0._wp
        vc%snow%w_snow         = 0._wp
        vc%snow%w_snow_old     = 0._wp
        vc%snow%w_snow_max     = 0._wp
        vc%snow%snowmelt       = 0._wp
        vc%snow%icemelt        = 0._wp
        vc%snow%icesub         = 0._wp
        vc%snow%refreezing     = 0._wp
        vc%snow%refreezing_sum = 0._wp
        vc%snow%dt_snowfree    = 0._wp
        vc%snow%snow_grain     = 0._wp
        vc%snow%dust_con       = 0._wp
        vc%snow%alb_snow_vis_dir = 0._wp
        vc%snow%alb_snow_nir_dir = 0._wp
        vc%snow%alb_snow_vis_dif = 0._wp
        vc%snow%alb_snow_nir_dif = 0._wp
        allocate(vc%snow%w_snow_iso(nwiso), vc%snow%w_snow_iso_old(nwiso))
        allocate(vc%snow%snowmelt_iso(nwiso), vc%snow%icemelt_iso(nwiso), vc%snow%icesub_iso(nwiso))
        vc%snow%w_snow_iso(:)     = 0._wp
        vc%snow%w_snow_iso_old(:) = 0._wp
        vc%snow%snowmelt_iso(:)   = 0._wp
        vc%snow%icemelt_iso(:)    = 0._wp
        vc%snow%icesub_iso(:)     = 0._wp

        return

    end subroutine alloc_lake_blocks

    subroutine alloc_land_blocks(vc)
        ! Size + neutral-initialize the inner arrays a vegetated land virtual
        ! cell needs for the port-C surface energy-balance chain. Array bounds
        ! mirror the reference lnd allocation (src/lnd/lnd_model.f90 ~1760):
        ! surface_flux over the nveg (bare+PFT) sub-tiles, soil column over nl,
        ! veg tiles over npft, one shared snowpack. Values are neutral here
        ! (T0/0); the physical state is seeded from the coarse cell in
        ! cmn_to_lndvc (port C decision: seed veg+soil state, decision B).
        ! Done once, guarded by allocation.

        implicit none

        type(vc_t), intent(inout) :: vc

        ! t_soil allocation doubles as the once-guard
        if (allocated(vc%soil%t_soil)) return

        ! --- shared surface-flux block over the nsurf tiles -------------------
        ! Sized nsurf (not nveg): the pattern-A veg routines (ebal_veg,
        ! update_tskin_veg, canopy_water) use whole-array flag_veg/flag_pft
        ! masks of length nsurf, so the flux arrays must conform. The veg tiles
        ! (1..nveg) are active; the lake/ice slots (i_lake,i_ice) stay zero.
        call alloc_surface_flux(vc%flx, nsurf)

        ! --- soil column (thermal) -------------------------------------------
        allocate(vc%soil%t_soil(0:nl), vc%soil%t_soil_old(0:nl), vc%soil%t_soil_max(0:nl))
        allocate(vc%soil%lambda_soil(0:nl), vc%soil%lambda_int_soil(0:nl), vc%soil%cap_soil(0:nl))
        allocate(vc%soil%theta_w(nl), vc%soil%theta_i(nl), vc%soil%theta(nl))
        allocate(vc%soil%w_w(nl), vc%soil%w_i(nl), vc%soil%w_w_old(nl), vc%soil%w_i_old(nl))
        allocate(vc%soil%w_w_phase(nl), vc%soil%w_i_phase(nl))
        allocate(vc%soil%lambda_s(nl), vc%soil%lambda_dry(nl), vc%soil%kappa_int(nl))
        allocate(vc%soil%theta_sat(nl), vc%soil%k_sat(nl), vc%soil%psi_sat(nl))
        allocate(vc%soil%theta_field(nl), vc%soil%theta_wilt(nl), vc%soil%psi(nl))
        allocate(vc%soil%k_exp(nl), vc%soil%psi_exp(nl))
        allocate(vc%soil%t_soil_cum(nl), vc%soil%theta_w_cum(nl), vc%soil%theta_i_cum(nl))
        ! permafrost
        allocate(vc%soil%frozen_years(nl), vc%soil%thaw_timer(nl))
        ! hydrology / water balance (single land tile)
        allocate(vc%soil%runoff(1), vc%soil%runoff_sur(1), vc%soil%calving(1))
        allocate(vc%soil%drainage(1), vc%soil%water_cons(1))
        ! water deficit + static wetland parameters
        allocate(vc%soil%cwd_mon(nmon_year))
        allocate(vc%soil%cti_cdf(15))
        ! rooting / wilting (nl,npft)
        allocate(vc%soil%wilt(nl,npft), vc%soil%root_frac(nl,npft))
        ! soil water isotopes
        allocate(vc%soil%w_w_iso(nl,nwiso), vc%soil%w_i_iso(nl,nwiso))
        allocate(vc%soil%w_w_iso_old(nl,nwiso), vc%soil%w_i_iso_old(nl,nwiso))
        allocate(vc%soil%infiltration_iso(nwiso))
        allocate(vc%soil%runoff_iso(nwiso), vc%soil%runoff_sur_iso(nwiso))
        allocate(vc%soil%drainage_iso(nwiso), vc%soil%calving_iso(nwiso))
        allocate(vc%soil%water_iso_cons(nwiso))

        vc%soil%t_soil(:)        = T0
        vc%soil%t_soil_old(:)    = T0
        vc%soil%t_soil_max(:)    = T0
        vc%soil%lambda_soil(:)   = 0._wp
        vc%soil%lambda_int_soil(:) = 0._wp
        vc%soil%cap_soil(:)      = 0._wp
        vc%soil%theta_w(:)       = 0._wp
        vc%soil%theta_i(:)       = 0._wp
        vc%soil%theta(:)         = 0._wp
        vc%soil%w_w(:)           = 0._wp
        vc%soil%w_i(:)           = 0._wp
        vc%soil%w_w_old(:)       = 0._wp
        vc%soil%w_i_old(:)       = 0._wp
        vc%soil%w_w_phase(:)     = 0._wp
        vc%soil%w_i_phase(:)     = 0._wp
        vc%soil%lambda_s(:)      = 0._wp
        vc%soil%lambda_dry(:)    = 0._wp
        vc%soil%kappa_int(:)     = 0._wp
        vc%soil%theta_sat(:)     = 0._wp
        vc%soil%k_sat(:)         = 0._wp
        vc%soil%psi_sat(:)       = 0._wp
        vc%soil%theta_field(:)   = 0._wp
        vc%soil%theta_wilt(:)    = 0._wp
        vc%soil%psi(:)           = 0._wp
        vc%soil%k_exp(:)         = 0
        vc%soil%psi_exp(:)       = 0
        vc%soil%t_soil_cum(:)    = 0._wp
        vc%soil%theta_w_cum(:)   = 0._wp
        vc%soil%theta_i_cum(:)   = 0._wp
        vc%soil%frozen_years(:)  = 0._wp
        vc%soil%thaw_timer(:)    = 0._wp
        vc%soil%runoff(:)        = 0._wp
        vc%soil%runoff_sur(:)    = 0._wp
        vc%soil%calving(:)       = 0._wp
        vc%soil%drainage(:)      = 0._wp
        vc%soil%water_cons(:)    = 0._wp
        vc%soil%wilt(:,:)        = 0._wp
        vc%soil%root_frac(:,:)   = 0._wp
        vc%soil%w_w_iso(:,:)     = 0._wp
        vc%soil%w_i_iso(:,:)     = 0._wp
        vc%soil%w_w_iso_old(:,:) = 0._wp
        vc%soil%w_i_iso_old(:,:) = 0._wp
        vc%soil%infiltration_iso(:) = 0._wp
        vc%soil%runoff_iso(:)    = 0._wp
        vc%soil%runoff_sur_iso(:) = 0._wp
        vc%soil%drainage_iso(:)  = 0._wp
        vc%soil%calving_iso(:)   = 0._wp
        vc%soil%water_iso_cons(:) = 0._wp
        ! soil scalars
        vc%soil%alt          = 0._wp
        vc%soil%infiltration = 0._wp
        vc%soil%w_table      = 0._wp
        vc%soil%w_table_peat = 0._wp
        vc%soil%w_table_cum  = 0._wp
        vc%soil%w_table_min  = 0._wp
        vc%soil%f_wet        = 0._wp
        vc%soil%f_wet_cum    = 0._wp
        vc%soil%f_wet_max    = 0._wp
        vc%soil%f_wetland    = 0._wp
        vc%soil%cti_lim      = 0._wp
        vc%soil%f_wet_mon    = 0._wp
        vc%soil%w_table_mon  = 0._wp
        vc%soil%f_wet_long   = 0._wp
        vc%soil%runoff_ann   = 0._wp
        vc%soil%pet          = 0._wp
        vc%soil%mcwd         = 0._wp
        vc%soil%mcwd_clim    = 0._wp
        vc%soil%cwd_mon(:)   = 0._wp
        ! static wetland parameters (seeded from topmodel/dyptop in lndvc_init_land)
        vc%soil%cti_mean     = 0._wp
        vc%soil%cti_cdf(:)   = 0._wp
        vc%soil%dyptop_k     = 0._wp
        vc%soil%dyptop_v     = 0._wp
        vc%soil%dyptop_xm    = 0._wp
        vc%soil%dyptop_fmax  = 0._wp

        ! --- vegetation tiles (npft) -----------------------------------------
        allocate(vc%veg%ci(npft), vc%veg%g_can(npft), vc%veg%gpp(npft), vc%veg%npp(npft))
        allocate(vc%veg%npp13(npft), vc%veg%npp14(npft), vc%veg%aresp(npft), vc%veg%discrimination(npft))
        allocate(vc%veg%lai(npft), vc%veg%sai(npft), vc%veg%phen(npft), vc%veg%phen_acc(npft))
        allocate(vc%veg%gdd(npft), vc%veg%gamma_leaf(npft), vc%veg%lambda(npft), vc%veg%lai_bal(npft))
        allocate(vc%veg%npp_cum(npft), vc%veg%npp13_cum(npft), vc%veg%npp14_cum(npft))
        allocate(vc%veg%npp_ann(npft), vc%veg%npp13_ann(npft), vc%veg%npp14_ann(npft))
        allocate(vc%veg%veg_c(npft), vc%veg%veg_h(npft), vc%veg%pft_frac(npft), vc%veg%seed_frac(npft))
        allocate(vc%veg%veg_c_below(nl), vc%veg%veg_c13_below(nl), vc%veg%veg_c14_below(nl))
        allocate(vc%veg%leaf_c(npft), vc%veg%stem_c(npft), vc%veg%root_c(npft), vc%veg%veg_c13(npft), vc%veg%veg_c14(npft))
        allocate(vc%veg%fire_c_flux_pft(npft), vc%veg%fire_c13_flux_pft(npft), vc%veg%fire_c14_flux_pft(npft))
        allocate(vc%veg%gamma_fire(npft), vc%veg%gamma_fire_cum(npft), vc%veg%gamma_luc(npft))
        allocate(vc%veg%gamma_ice(npft), vc%veg%gamma_dist(npft), vc%veg%gamma_dist_cum(npft))

        vc%veg%ci=0._wp; vc%veg%g_can=0._wp; vc%veg%gpp=0._wp; vc%veg%npp=0._wp
        vc%veg%npp13=0._wp; vc%veg%npp14=0._wp; vc%veg%aresp=0._wp; vc%veg%discrimination=0._wp
        vc%veg%lai=0._wp; vc%veg%sai=0._wp; vc%veg%phen=0._wp; vc%veg%phen_acc=0._wp
        vc%veg%gdd=0._wp; vc%veg%gamma_leaf=0._wp; vc%veg%lambda=0._wp; vc%veg%lai_bal=0._wp
        vc%veg%npp_cum=0._wp; vc%veg%npp13_cum=0._wp; vc%veg%npp14_cum=0._wp
        vc%veg%npp_ann=0._wp; vc%veg%npp13_ann=0._wp; vc%veg%npp14_ann=0._wp
        vc%veg%veg_c=0._wp; vc%veg%veg_h=0._wp; vc%veg%pft_frac=0._wp; vc%veg%seed_frac=0._wp
        vc%veg%veg_c_below=0._wp; vc%veg%veg_c13_below=0._wp; vc%veg%veg_c14_below=0._wp
        vc%veg%leaf_c=0._wp; vc%veg%stem_c=0._wp; vc%veg%root_c=0._wp; vc%veg%veg_c13=0._wp; vc%veg%veg_c14=0._wp
        vc%veg%fire_c_flux_pft=0._wp; vc%veg%fire_c13_flux_pft=0._wp; vc%veg%fire_c14_flux_pft=0._wp
        vc%veg%gamma_fire=0._wp; vc%veg%gamma_fire_cum=0._wp; vc%veg%gamma_luc=0._wp
        vc%veg%gamma_ice=0._wp; vc%veg%gamma_dist=0._wp; vc%veg%gamma_dist_cum=0._wp
        vc%veg%gdd5=0._wp; vc%veg%gdd5_temp=0._wp
        vc%veg%t2m_min_mon=0._wp; vc%veg%t2m_ann_mean=0._wp
        vc%veg%f_crop=0._wp; vc%veg%f_pasture=0._wp; vc%veg%df_crop=0._wp; vc%veg%df_pasture=0._wp
        vc%veg%fire_c_flux=0._wp; vc%veg%fire_c13_flux=0._wp; vc%veg%fire_c14_flux=0._wp
        vc%veg%carbon_bal_veg=0._wp; vc%veg%carbon13_bal_veg=0._wp; vc%veg%carbon14_bal_veg=0._wp
        allocate(vc%veg%disturbance(npft)); vc%veg%disturbance=0._wp

        ! --- soil carbon fields the thermal chain / init read ----------------
        allocate(vc%carb%soil_resp_l(nl,ncarb), vc%carb%litter_in_frac(nl))
        vc%carb%soil_resp_l(:,:)  = 0._wp
        vc%carb%litter_in_frac(:) = 0._wp
        vc%carb%f_peat     = 0._wp
        vc%carb%f_peat_pot = 0._wp
        vc%carb%dCpeat_dt  = 0._wp
        ! carbon pools read/written by the veg-dynamics chain (soil_par_update,
        ! dynveg_par, dyn_veg). Lean zero-init (port C.3, option A); physical
        ! carbon cold-start + soil_carbon updates are port C.4.
        allocate(vc%carb%litter_c(nl), vc%carb%fast_c(nl), vc%carb%slow_c(nl))
        allocate(vc%carb%cato_c(nl), vc%carb%frac_soc(nl))
        allocate(vc%carb%litterfall(nlc,ncarb), vc%carb%litterfall13(nlc,ncarb), vc%carb%litterfall14(nlc,ncarb))
        vc%carb%litter_c(:) = 0._wp; vc%carb%fast_c(:) = 0._wp; vc%carb%slow_c(:) = 0._wp
        vc%carb%cato_c(:)   = 0._wp; vc%carb%frac_soc(:) = 0._wp
        vc%carb%litterfall(:,:) = 0._wp; vc%carb%litterfall13(:,:) = 0._wp; vc%carb%litterfall14(:,:) = 0._wp
        vc%carb%litter_c_peat = 0._wp; vc%carb%acro_c = 0._wp

        ! --- shared snowpack (single land snow model) ------------------------
        vc%snow%mask_snow      = 0
        vc%snow%f_snow         = 0._wp
        vc%snow%h_snow         = 0._wp
        vc%snow%w_snow         = 0._wp
        vc%snow%w_snow_old     = 0._wp
        vc%snow%w_snow_max     = 0._wp
        vc%snow%snowmelt       = 0._wp
        vc%snow%icemelt        = 0._wp
        vc%snow%icesub         = 0._wp
        vc%snow%refreezing     = 0._wp
        vc%snow%refreezing_sum = 0._wp
        vc%snow%dt_snowfree    = 0._wp
        vc%snow%snow_grain     = 0._wp
        vc%snow%dust_con       = 0._wp
        vc%snow%alb_snow_vis_dir = 0._wp
        vc%snow%alb_snow_nir_dir = 0._wp
        vc%snow%alb_snow_vis_dif = 0._wp
        vc%snow%alb_snow_nir_dif = 0._wp
        allocate(vc%snow%w_snow_iso(nwiso), vc%snow%w_snow_iso_old(nwiso))
        allocate(vc%snow%snowmelt_iso(nwiso), vc%snow%icemelt_iso(nwiso), vc%snow%icesub_iso(nwiso))
        vc%snow%w_snow_iso(:)     = 0._wp
        vc%snow%w_snow_iso_old(:) = 0._wp
        vc%snow%snowmelt_iso(:)   = 0._wp
        vc%snow%icemelt_iso(:)    = 0._wp
        vc%snow%icesub_iso(:)     = 0._wp

        return

    end subroutine alloc_land_blocks

    subroutine alloc_surface_flux(flx, n)
        ! Allocate the full shared surface-flux block for n sub-tiles
        ! (init 0, skin temperatures T0). Lake uses n=1.

        implicit none

        type(surface_flux_t), intent(inout) :: flx
        integer, intent(in) :: n

        allocate(flx%rough_m(n), flx%rough_h(n), flx%Ch(n), flx%z0m(n), flx%Ri(n))
        allocate(flx%r_a(n), flx%r_s(n), flx%beta_s(n))
        allocate(flx%r_a_can(n), flx%r_s_can(n), flx%beta_s_can(n))
        allocate(flx%albedo(n), flx%alb_vis_dir(n), flx%alb_vis_dif(n), flx%alb_nir_dir(n), flx%alb_nir_dif(n))
        allocate(flx%flx_sh(n), flx%flx_lh(n), flx%flx_g(n), flx%dflxg_dT(n), flx%flx_melt(n), flx%flx_lwu(n), flx%lwnet(n))
        allocate(flx%t_skin(n), flx%t_skin_old(n), flx%t_skin_amp(n))
        allocate(flx%num_lh(n), flx%num_sh(n), flx%num_sw(n), flx%num_lw(n), flx%denom_lh(n), flx%denom_sh(n), flx%denom_lw(n))
        allocate(flx%f_sh(n), flx%f_e(n), flx%f_t(n), flx%f_le(n), flx%f_lt(n), flx%f_lw(n), flx%lh_ecan(n))
        allocate(flx%qsat_e(n), flx%dqsatdT_e(n), flx%qsat_t(n), flx%dqsatdT_t(n))
        allocate(flx%transpiration(n), flx%evap_surface(n), flx%et(n))
        allocate(flx%rain_ground(n), flx%evap_can(n), flx%snow_ground(n), flx%subl_can(n))
        allocate(flx%w_can(n), flx%w_can_old(n), flx%s_can(n), flx%s_can_old(n), flx%f_wat_can(n), flx%f_snow_can(n))
        allocate(flx%f_snow(n))
        allocate(flx%frac_surf(n))
        ! water isotopes (n,nwiso)
        allocate(flx%rain_iso(n,nwiso), flx%snow_iso(n,nwiso))
        allocate(flx%rain_ground_iso(n,nwiso), flx%snow_ground_iso(n,nwiso))
        allocate(flx%evap_can_iso(n,nwiso), flx%subl_can_iso(n,nwiso))
        allocate(flx%transpiration_iso(n,nwiso), flx%evap_surface_iso(n,nwiso), flx%et_iso(n,nwiso))
        allocate(flx%w_can_iso(n,nwiso), flx%w_can_iso_old(n,nwiso), flx%s_can_iso(n,nwiso), flx%s_can_iso_old(n,nwiso))

        flx%rough_m=0._wp; flx%rough_h=0._wp; flx%Ch=0._wp; flx%z0m=0._wp; flx%Ri=0._wp
        flx%r_a=0._wp; flx%r_s=0._wp; flx%beta_s=0._wp
        flx%r_a_can=0._wp; flx%r_s_can=0._wp; flx%beta_s_can=0._wp
        flx%albedo=0._wp; flx%alb_vis_dir=0._wp; flx%alb_vis_dif=0._wp; flx%alb_nir_dir=0._wp; flx%alb_nir_dif=0._wp
        flx%flx_sh=0._wp; flx%flx_lh=0._wp; flx%flx_g=0._wp; flx%dflxg_dT=0._wp; flx%flx_melt=0._wp; flx%flx_lwu=0._wp; flx%lwnet=0._wp
        flx%t_skin=T0; flx%t_skin_old=T0; flx%t_skin_amp=0._wp
        flx%num_lh=0._wp; flx%num_sh=0._wp; flx%num_sw=0._wp; flx%num_lw=0._wp; flx%denom_lh=0._wp; flx%denom_sh=0._wp; flx%denom_lw=0._wp
        flx%f_sh=0._wp; flx%f_e=0._wp; flx%f_t=0._wp; flx%f_le=0._wp; flx%f_lt=0._wp; flx%f_lw=0._wp; flx%lh_ecan=0._wp
        flx%qsat_e=0._wp; flx%dqsatdT_e=0._wp; flx%qsat_t=0._wp; flx%dqsatdT_t=0._wp
        flx%transpiration=0._wp; flx%evap_surface=0._wp; flx%et=0._wp
        flx%rain_ground=0._wp; flx%evap_can=0._wp; flx%snow_ground=0._wp; flx%subl_can=0._wp
        flx%w_can=0._wp; flx%w_can_old=0._wp; flx%s_can=0._wp; flx%s_can_old=0._wp; flx%f_wat_can=0._wp; flx%f_snow_can=0._wp
        flx%f_snow=0._wp
        flx%frac_surf=0._wp
        flx%rain_iso=0._wp; flx%snow_iso=0._wp
        flx%rain_ground_iso=0._wp; flx%snow_ground_iso=0._wp
        flx%evap_can_iso=0._wp; flx%subl_can_iso=0._wp
        flx%transpiration_iso=0._wp; flx%evap_surface_iso=0._wp; flx%et_iso=0._wp
        flx%w_can_iso=0._wp; flx%w_can_iso_old=0._wp; flx%s_can_iso=0._wp; flx%s_can_iso_old=0._wp

        return

    end subroutine alloc_surface_flux

    subroutine lndvc_remap_state(lnd)
        ! Conservatively transfer prognostic state (carbon, heat, snow, water)
        ! between leaf virtual cells when band membership or class changes.
        ! Replaces the ad-hoc "initialize newly vegetated/ice/lake cell" branches
        ! scattered through lnd_update today.

        implicit none

        type(lndvc_class), intent(inout) :: lnd

        ! TODO: detect decomposition change; conservative redistribution.

        return

    end subroutine lndvc_remap_state

end module lndvc_decomp

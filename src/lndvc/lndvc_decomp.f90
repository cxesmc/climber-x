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
    use smb_grid_m, only : nl_smb => nl   ! snow+ice/firn layer count (=4; distinct from land nl=5)

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

        select case(class)
            case(1)   ! land
                if (.not. allocated(vc(k)%veg))  allocate(vc(k)%veg)
                if (.not. allocated(vc(k)%soil)) allocate(vc(k)%soil)
                if (.not. allocated(vc(k)%carb)) allocate(vc(k)%carb)
                if (.not. allocated(vc(k)%snow)) allocate(vc(k)%snow)
            case(2)   ! lake
                if (.not. allocated(vc(k)%lake)) allocate(vc(k)%lake)
                if (.not. allocated(vc(k)%snow)) allocate(vc(k)%snow)
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

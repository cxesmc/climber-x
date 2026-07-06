module lndvc_downscale
    ! Forcing-downscaling service for the virtual-cell framework.
    !
    ! Shared infrastructure (promoted from smb/downscaling.f90) used by every
    ! surface class. Given coarse/mid-res forcing and the reference elevation it
    ! is valid at, plus a virtual cell's {z, slope}, returns vc-level forcing
    ! (T/q/LW/SW/precip/wind). See docs/design/virtual-cells.md, sec 5.
    !
    ! Invariant: interpolate forcing AND its reference elevation together, so the
    ! downscaling delta (z_vc - z_ref) is correct (cf. SMB z_sur_i vs z_sur).

    use precision, only : wp
    use lndvc_def

    implicit none

    private
    public :: lndvc_downscale_forcing    ! coarse/mid-res forcing + vc {z,slope} -> vc forcing

contains

    subroutine lndvc_downscale_forcing(smb, vc)
        ! Downscale forcing to a single virtual cell's elevation and slope.
        ! Temperature/humidity via lapse rates; longwave via LW lapse rate;
        ! shortwave via first-order (albedo, elevation) sensitivities;
        ! precipitation via elevation factor + orographic wind-slope factor.

        implicit none

        type(lndvc_smb_class), intent(inout) :: smb   ! placeholder carrier for forcing block
        type(lndvc_vc_class),  intent(in)    :: vc

        ! TODO: port temp/wind/prc/rad downscaling from smb/downscaling.f90 as
        !       single-column operations keyed on (vc%zsrf - z_ref).

        return

    end subroutine lndvc_downscale_forcing

end module lndvc_downscale

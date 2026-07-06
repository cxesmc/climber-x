module lndvc_downscale
    ! Forcing-downscaling service for the virtual-cell framework.
    !
    ! Shared infrastructure (promoted from smb/downscaling.f90) used by every
    ! surface class. Given coarse/mid-res forcing and the reference elevation it
    ! is valid at, plus a virtual cell's {z, slope}, returns vc-level forcing.
    ! See docs/design/virtual-cells.md, sec 5.
    !
    ! Invariant: interpolate forcing AND its reference elevation together, so the
    ! downscaling delta (vc%z - z_ref) is correct (cf. SMB z_sur_i vs z_sur).

    use precision, only : wp
    use lndvc_def, only : forcing_t, vc_desc_t

    implicit none

    private
    public :: lndvc_downscale_forcing    ! reference forcing + vc {z,slope} -> vc forcing

contains

    subroutine lndvc_downscale_forcing(forc, desc)
        ! Downscale forcing to a single virtual cell's elevation and slope.
        ! Temperature/humidity via lapse rates; longwave via LW lapse rate;
        ! shortwave via first-order (albedo, elevation) sensitivities;
        ! precipitation via elevation factor + orographic wind-slope factor.

        implicit none

        type(forcing_t), intent(inout) :: forc
        type(vc_desc_t), intent(in)    :: desc

        ! TODO: port temp/wind/prc/rad downscaling from smb/downscaling.f90 as
        !       single-column operations keyed on (desc%z - z_ref).

        return

    end subroutine lndvc_downscale_forcing

end module lndvc_downscale

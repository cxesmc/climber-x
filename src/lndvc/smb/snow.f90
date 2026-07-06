module smb_snow_m
    ! Snow-layer update for the SMB single-column physics.
    !
    ! Ported into the lndvc framework from src/smb/snow.f90 (module snow_mod)
    ! with the body unchanged; only the parameter/constant dependencies are
    ! re-homed onto the framework-owned smb_par_m so lndvc is self-contained.

    use precision, only : wp
    use smb_par_m, only : dt, snow_par, check_water

    implicit none

    private
    public :: snow_update

contains

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   Subroutine :  s n o w _ u p d a t e
    !   Purpose    :  update snow layer
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine snow_update(mask_snow,evp, &
                                w_snow,w_snow_old,w_snow_max, &
                                h_snow)

        implicit none

        integer, intent(inout) :: mask_snow
        real(wp), intent(in) :: evp
        real(wp), intent(inout) :: w_snow, w_snow_old, w_snow_max
        real(wp), intent(out) :: h_snow


        ! snow water equivalent evolution
        ! remove sublimation, snowfall has been added already and snowmelt already removed and refreezing added during soil temperature update
        if (mask_snow.eq.1) then
          w_snow = w_snow - evp*dt  ! kg/m2
        endif

        if (check_water .and. w_snow.lt.0._wp) print *,'WARNING w_snow < 0',w_snow
        w_snow = max(0._wp,w_snow)  ! not strictly conserving water here!

        ! limit w_snow
        if (w_snow.gt.snow_par%w_snow_max) w_snow = snow_par%w_snow_max

        ! update snow thickness
        h_snow = w_snow / snow_par%rho  ! m

        ! update snow mask
        if( w_snow .gt. snow_par%w_snow_crit ) then
          mask_snow = 1
        else
          mask_snow = 0
        endif

        ! save seasonal maximum snow swe
        if (w_snow.gt.w_snow_old .or. w_snow.eq.snow_par%w_snow_max) w_snow_max = w_snow

        return

    end subroutine snow_update

end module smb_snow_m

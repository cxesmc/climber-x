module bgc_model

  use precision, only : wp
  use bgc_def, only : bgc_class

  implicit none

  private
  public :: bgc_ini, bgc_update, bgc_end, bgc_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   subroutine :  b g c _ i n i
  !   purpose    :  initialize marine bio-geo-chemistry module
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_ini(bgc, ni, nj, lon, lat, layer_thk, layer_depth, level_depth, kbo, ocn_area, A_shelf_geo, l_daily_input_save)

    type(bgc_class), intent(out) :: bgc          !! state variables of bgc model, to initialize

    integer, intent(in) :: ni, nj              !! horizontal dimensions
    real(wp), intent(in) :: lon(:) 
    real(wp), intent(in) :: lat(:)  
    real(wp), intent(in) :: layer_thk(:)                  !! layer thickness [m] -- first index for surface layer
    real(wp), intent(in) :: layer_depth(:)                !! depth of layer centers [m] -- first index for surface layer
    real(wp), intent(in) :: level_depth(:)                !! depth of layer interfaces [m] -- first index for surface layer
    integer, intent(in) :: kbo(:, :)         !! index of bottom layer (1<= <= maxk if ocean, 0 otherwise)
    real(wp), intent(in) :: ocn_area(:,:)                 !! ocean grid cell area [m2]
    real(wp), intent(in) :: A_shelf_geo
    logical, intent(in) :: l_daily_input_save

    return

  end subroutine bgc_ini


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! subroutine :  b g c _ u p d a t e
  ! purpose  :  wrapper for column bgc update
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_update(bgc)

    type(bgc_class), intent(inout) :: bgc

  end subroutine bgc_update


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! subroutine :  b g c _ u p d a t e
  ! purpose  :  end bgc model
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_end(bgc)

    type(bgc_class), intent(inout) :: bgc

  end subroutine bgc_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! subroutine :  b g c _ w r i t e _ r e s t a r t
  ! purpose  :  write restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_write_restart(fnm,bgc)

    implicit none

    character (len=*) :: fnm
    type(bgc_class) :: bgc

   return

  end subroutine bgc_write_restart


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! subroutine :  b g c _ r e a d _ r e s t a r t
  ! purpose  :  read restart netcdf file 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bgc_read_restart(fnm,bgc)

    implicit none

    character (len=*) :: fnm
    type(bgc_class) :: bgc

   return

  end subroutine bgc_read_restart

end module


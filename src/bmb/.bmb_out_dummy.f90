module bmb_out

  use bmb_def, only : bmb_class

  implicit none

  private
  public :: bmb_diag, bmb_diag_init


contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ d i a g _ i n i t
  ! Purpose  :  Initialize netcdf output for bmb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_diag_init(bmb)

    implicit none

    type(bmb_class) :: bmb

   return

  end subroutine bmb_diag_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b m b _ d i a g
  !   Purpose    :  sea ice diagnostics
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_diag(bmb)

    implicit none

    type(bmb_class) :: bmb

   return

  end subroutine bmb_diag

end module bmb_out

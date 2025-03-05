module bmb_model

  use precision, only : wp, dp
  use coord, only : grid_class
  use bmb_def, only : bmb_class

  implicit none

  private
  public :: bmb_init, bmb_update, bmb_end, bmb_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b m b _ u p d a t e
  !   Purpose    :  update ice melt to ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_update(bmb)

  !$  use omp_lib

  implicit none

  type(bmb_class) :: bmb

  return

  end subroutine bmb_update
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  b m b _ i n i t
  !   Purpose    :  initialize ice melt to ocean 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_init(bmb,grid,cmn_grid)

    implicit none

    type(bmb_class) :: bmb
    type(grid_class), intent(in) :: grid 
    type(grid_class), intent(in) :: cmn_grid 

  return

  end subroutine bmb_init
      
     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ e n d 
  ! Purpose  :  end bmb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_end(bmb)

    implicit none

    type(bmb_class) :: bmb

    return

  end subroutine bmb_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  b m b _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine bmb_write_restart(fnm,bmb)

    implicit none

    character (len=*) :: fnm
    type(bmb_class) :: bmb

  end subroutine bmb_write_restart

end module bmb_model

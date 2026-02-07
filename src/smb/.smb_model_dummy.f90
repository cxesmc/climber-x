module smb_model

  use precision, only : wp
  use coord, only : grid_class
  use smb_def, only : smb_class, smb_in_class

  implicit none

  private
  public :: smb_init, smb_update, smb_end, smb_write_restart

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ u p d a t e
  !   Purpose    :  update surface energy and mass balance interface
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_update(smb_in,smb)

  !$  use omp_lib

  implicit none

  type(smb_in_class), intent(in) :: smb_in
  type(smb_class) :: smb

  return

  end subroutine smb_update
      

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  s m b _ i n i t
  !   Purpose    :  initialize surface energy and mass balance interface
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_init(smb_in,smb,i_domain,grid,cmn_grid)

    implicit none

    type(smb_in_class) :: smb_in
    type(smb_class) :: smb
    integer, intent(in) :: i_domain
    type(grid_class), intent(in) :: grid 
    type(grid_class), intent(in) :: cmn_grid

  return

  end subroutine smb_init
      
     
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ e n d 
  ! Purpose  :  end smb
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_end(smb_in, smb)

    implicit none

    type(smb_in_class) :: smb_in
    type(smb_class) :: smb

    return

  end subroutine smb_end


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Function :  s m b _ w r i t e _ r e s t a r t
  ! Purpose  :  Write restart file
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine smb_write_restart(fnm,smb)

    implicit none

    character (len=*) :: fnm
    type(smb_class) :: smb

  end subroutine smb_write_restart

end module smb_model

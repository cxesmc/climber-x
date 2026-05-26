!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : w i s o _ p a r a m s
!
!  Purpose : water isotope tracer parameters and utilities
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module wiso_params

  use precision, only : wp

  implicit none

  integer, parameter :: nwiso = 1
  integer, parameter :: i_o18 = 1

  logical :: l_wiso = .false.

  ! VSMOW standard ratios
  real(wp), parameter :: Rstd_o18 = 2005.20e-6_wp  ! 18O/16O

  ! VSMOW ratios array indexed by isotope species
  real(wp), parameter, dimension(nwiso) :: Rstd = (/ Rstd_o18 /)

contains

  pure function wiso_ratio(m_wiso, m_bulk) result(R)
    real(wp), intent(in) :: m_wiso, m_bulk
    real(wp) :: R
    if (m_bulk > 0._wp) then
      R = m_wiso / m_bulk
    else
      R = 0._wp
    endif
  end function

  pure function delta_from_ratio(R, Rstd_val) result(delta)
    real(wp), intent(in) :: R, Rstd_val
    real(wp) :: delta
    if (Rstd_val > 0._wp) then
      delta = (R / Rstd_val - 1._wp) * 1000._wp
    else
      delta = 0._wp
    endif
  end function

end module wiso_params

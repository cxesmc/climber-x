module smb_grid_m
    ! Vertical grid for the SMB single-column physics (snow + ice/firn layers).
    !
    ! Ported into the lndvc framework from src/smb/smb_grid.f90 (module smb_grid)
    ! with the body unchanged; module renamed to avoid clashing with the original
    ! during the transition.

    use precision, only : wp

    implicit none

    ! levels
    integer,  parameter :: nl = 4
    real(wp), dimension(0:nl) :: z = (/0._wp,0.2_wp,1._wp,5._wp,15._wp/)! z(0) is overwritten with 0.5*h_snow

    real(wp), dimension(0:nl) :: z_int, dz, rdz, rdz_pos, rdz_neg

contains

    subroutine smb_grid_init

    implicit none

    integer :: k


      ! vertical grid
      z(0) = 0._wp
      dz(0) = 0._wp
      ! vertical layers thickness
      dz(1) = 0.5_wp * ( z(1) + z(2) )
      do k=2,nl-1
       dz(k) = 0.5_wp * ( z(k+1) - z(k-1) )
      enddo
      dz(nl) = z(nl) - z(nl-1)

      ! depth of vertical layer interfaces
      z_int(0) = 0._wp ! snow - soil interface
      do k=1,nl-1
       z_int(k) = 0.5_wp * ( z(k) + z(k+1) )
      enddo
      z_int(nl) = z(nl) + 0.5_wp * dz(nl)

      ! reciprocals to speed up fortran
      rdz(1:nl) = 1._wp/dz(1:nl)
      do k=1,nl-1
       rdz_pos(k) = 1._wp/(z(k+1)-z(k))
      enddo
      rdz_pos(nl) = 0._wp
      do k=2,nl
       rdz_neg(k) = 1._wp/(z(k)-z(k-1))
      enddo

    return

    end subroutine smb_grid_init

end module smb_grid_m

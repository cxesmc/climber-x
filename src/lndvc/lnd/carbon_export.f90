!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : c a r b o n _ e x p o r t _ m o d
!
!  Purpose : river export of organic carbon 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Andrey Ganopolski
!
! This file is part of CLIMBER-X.
!
! CLIMBER-X is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! CLIMBER-X is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with CLIMBER-X.  If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module lndvc_carbon_export_mod

  use precision, only : wp
  use lnd_grid, only : nl, nlc, dz
  use lnd_params, only : kexport

  implicit none

  private
  public :: carbon_export

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  c a r b o n _ e x p o r t
  !   Purpose    :  land carbon export (POC and DOC) through rivers to ocean
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine carbon_export(runoff_ann, f_veg, f_peat, &
           litter_c, litter_c13, litter_c14, fast_c, fast_c13, fast_c14, slow_c, slow_c13, slow_c14, &
           litter_c_peat, litter_c13_peat, litter_c14_peat, acro_c, acro_c13, acro_c14, cato_c, cato_c13, cato_c14, &
           poc_export, poc13_export, poc14_export, &
           doc_export, doc13_export, doc14_export)

    implicit none

    real(wp), intent(in) :: runoff_ann
    real(wp), intent(in) :: f_veg, f_peat
    real(wp), dimension(:), intent(inout) :: litter_c, litter_c13, litter_c14
    real(wp), dimension(:), intent(inout) :: fast_c, fast_c13, fast_c14
    real(wp), dimension(:), intent(inout) :: slow_c, slow_c13, slow_c14
    real(wp), intent(inout) :: litter_c_peat, litter_c13_peat, litter_c14_peat
    real(wp), intent(inout) :: acro_c, acro_c13, acro_c14
    real(wp), dimension(:), intent(inout) :: cato_c, cato_c13, cato_c14
    real(wp), intent(out) :: poc_export, poc13_export, poc14_export 
    real(wp), intent(out) :: doc_export, doc13_export, doc14_export

    integer :: n
    real(wp) :: c_tot, c13_tot, c14_tot, export, d_c


    doc_export   = 0._wp
    doc13_export = 0._wp
    doc14_export = 0._wp
    poc_export   = 0._wp
    poc13_export = 0._wp
    poc14_export = 0._wp

    if (f_veg.gt.0._wp) then

      ! total soil carbon in vegetated grid cell (f_veg)
      c_tot   = ( (f_veg-f_peat)*sum((litter_c(1:nl)+fast_c(1:nl)+slow_c(1:nl))*dz(1:nl))       &
        + f_peat*(litter_c_peat+acro_c+sum(cato_c(1:nl)*dz(1:nl)))       )  ! kgC/m2
      c13_tot = ( (f_veg-f_peat)*sum((litter_c13(1:nl)+fast_c13(1:nl)+slow_c13(1:nl))*dz(1:nl)) &
        + f_peat*(litter_c13_peat+acro_c13+sum(cato_c13(1:nl)*dz(1:nl))) )  ! kgC/m2
      c14_tot = ( (f_veg-f_peat)*sum((litter_c14(1:nl)+fast_c14(1:nl)+slow_c14(1:nl))*dz(1:nl)) &
        + f_peat*(litter_c14_peat+acro_c14+sum(cato_c14(1:nl)*dz(1:nl))) )  ! kgC/m2

      ! carbon export, proportional to runoff and soil carbon content
      export = kexport * runoff_ann * c_tot  ! kgC/m2 = kg water/m2/yr * kgC/m2 * m2/kg water 

      if (export>0._wp) then

        ! remove exported carbon from soil carbon pools
        do n=1,nl
          d_c = (export * litter_c(n)*dz(n) / c_tot) / dz(n)
          if (litter_c(n)>0._wp) then
            litter_c13(n) = litter_c13(n) - d_c * litter_c13(n)/litter_c(n)
            litter_c14(n) = litter_c14(n) - d_c * litter_c14(n)/litter_c(n)
          endif
          litter_c(n) = litter_c(n) - d_c  ! kgC/m3
          d_c = (export * fast_c(n)*dz(n) / c_tot) / dz(n)
          if (fast_c(n)>0._wp) then
            fast_c13(n) = fast_c13(n) - d_c * fast_c13(n)/fast_c(n)
            fast_c14(n) = fast_c14(n) - d_c * fast_c14(n)/fast_c(n)
          endif
          fast_c(n)   = fast_c(n)   - d_c  ! kgC/m3
          d_c = (export * slow_c(n)*dz(n) / c_tot) / dz(n)
          if (slow_c(n)>0._wp) then
            slow_c13(n) = slow_c13(n) - d_c * slow_c13(n)/slow_c(n)
            slow_c14(n) = slow_c14(n) - d_c * slow_c14(n)/slow_c(n)
          endif
          slow_c(n)   = slow_c(n)   - d_c  ! kgC/m3
          d_c = (export * cato_c(n)*dz(n) / c_tot) / dz(n)
          if (cato_c(n)>0._wp) then
            cato_c13(n) = cato_c13(n) - d_c * cato_c13(n)/cato_c(n)
            cato_c14(n) = cato_c14(n) - d_c * cato_c14(n)/cato_c(n)
          endif
          cato_c(n)   = cato_c(n)   - d_c  ! kgC/m3
        enddo
        d_c = export * litter_c_peat / c_tot
        if (litter_c_peat>0._wp) then
          litter_c13_peat = litter_c13_peat - d_c * litter_c13_peat/litter_c_peat
          litter_c14_peat = litter_c14_peat - d_c * litter_c14_peat/litter_c_peat
        endif
        litter_c_peat   = litter_c_peat   - d_c  ! kgC/m2
        d_c = export * acro_c / c_tot
        if (acro_c>0._wp) then
          acro_c13 = acro_c13 - d_c * acro_c13/acro_c
          acro_c14 = acro_c14 - d_c * acro_c14/acro_c
        endif
        acro_c   = acro_c   - d_c  ! kgC/m2

        ! split export into dissolved and particulate carbon
        doc_export   = 0.5_wp*export / f_veg    ! kgC/m2 grid cell
        poc_export   = 0.5_wp*export / f_veg    ! kgC/m2 grid cell

        ! isotopes
        doc13_export = doc_export * c13_tot/c_tot
        doc14_export = doc_export * c14_tot/c_tot
        poc13_export = poc_export * c13_tot/c_tot
        poc14_export = poc_export * c14_tot/c_tot

      endif

    endif


    return

  end subroutine carbon_export

end module lndvc_carbon_export_mod


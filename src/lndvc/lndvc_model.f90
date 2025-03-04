module lndvc_model
    ! Interface to the land-virtual-cell (lndvc) model.
    ! 
    use ncio
    use precision, only : sp, dp, wp

    use lndvc_def
    use lndvc_grid 

    implicit none

    private
    public :: lndvc_init 
    public :: lndvc_aggregate_cell
    public :: lndvc_update_vc
    public :: lndvc_end

contains

    ! === Routines for whole domain control =============

    subroutine lndvc_init()

        implicit none



        return

    end subroutine lndvc_init

    subroutine lndvc_update()

        implicit none

        return

    end subroutine lndvc_update

    subroutine lndvc_end()

        implicit none

        call lndvc_dealloc()

        return

    end subroutine lndvc_end

    subroutine lndvc_alloc()

        implicit none



        return

    end subroutine lndvc_alloc

    subroutine lndvc_dealloc()

        implicit none



        return

    end subroutine lndvc_dealloc

    ! === Routines for vc column control ==============

    subroutine lndvc_set_vc(vc)
        ! Routine to set characteristics of a column
        ! of virtual cells. Set fractions etc. Needs
        ! high resolution topography as input probably. 

        implicit none

        type(lndvc_class), intent(INOUT) :: vc
        
        return

    end subroutine lndvc_set_vc

    subroutine lndvc_update_vc(veg,smb,vc)
        ! Update steps for a given virtual cell (vc)

        implicit none

        type(lndvc_veg_class), intent(INOUT) :: veg 
        type(lndvc_smb_class), intent(INOUT) :: smb 
        type(lndvc_vc_class),  intent(IN) :: vc

        ! Local variables
        ! ...

        select case(vc%surf)

            case(1) ! land

                ! Call all updates you need for a land point
                
                !call lndvc_update_veg(veg,vc)

                !call lndvc_update_soil(vc,soil,under_ice=.FALSE.)

                !call lndvc_update_carboncycle()

            case(2) ! lake

                ! Call all updates you need for a lake point

                !call lndvc_update_lake()

                !call lndvc_update_soil()

                !call lndvc_update_carboncycle()
                

            case(3) ! ice 

                ! Call all updates you need for a ice point

                call lndvc_update_smb(smb,vc)

                !call lndvc_update_soil(vc,soil,under_ice=.TRUE.)

        end select

        return

    end subroutine lndvc_update_vc

    subroutine lndvc_update_veg(veg,vc)

        implicit none

        type(lndvc_veg_class), intent(INOUT) :: veg 
        type(lndvc_vc_class),  intent(IN) :: vc



        return

    end subroutine lndvc_update_veg
    
    subroutine lndvc_update_smb(smb,vc)

        implicit none

        type(lndvc_smb_class), intent(INOUT) :: smb 
        type(lndvc_vc_class),  intent(IN) :: vc



        return

    end subroutine lndvc_update_smb
    
    subroutine lndvc_aggregate_cell(veg,smb,veg_vc,smb_vc,vc)
        ! Given a column of virtual cells, aggregate the values to 
        ! give the representative average cell values. 
        ! This routine can either be used to get the cell values
        ! for a large cell, or to interpolate down to a specific point.
        ! This will be determined by the weighting function obtained
        ! through analysis of vc(:).

        type(lndvc_veg_class), intent(OUT) :: veg 
        type(lndvc_smb_class), intent(OUT) :: smb 
        
        type(lndvc_veg_class), intent(IN) :: veg_vc(:)
        type(lndvc_smb_class), intent(IN) :: smb_vc(:)
        type(lndvc_vc_class),  intent(IN) :: vc(:)

        ! Local variables
        integer :: k, n
        real(wp), allocatable :: wt(:) 

        ! How many virtual cells are there?
        n = size(vc,1)

        allocate(wt(n))

        ! === Step 1 =======================
        ! Determine weights for averaging between all virtual cells

        wt = 0.0 

        do k = 1, n 

        end do 


        ! Normalize weights to sum to 1.0 
        if (sum(wt) .gt. 0.0) wt = wt / sum(wt)

        ! === Step 2 ======================
        ! Apply weighting function to all variables

        ! To do....
        ! call lndvc_aggregate_veg(veg,veg_vc,wt)
        ! call lndvc_aggregate_smb(smb,smb_vc,wt)
        

        return

    end subroutine lndvc_aggregate_cell

end module lndvc_model

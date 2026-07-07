module ice_model
    
    use precision, only : wp
    use ice_def, only : ice_class
    use coords, only : grid_class

    implicit none 

    private 
    public :: ice_update
    public :: ice_init_domains
    public :: ice_init 
    public :: ice_end 
    public :: ice_write_restart

contains 
    
    subroutine ice_update(ice,idx_dom,model,n_year_ice,time,time_out_ice)
        ! This routine serves as the climberx interface to an ice sheet
        ! model. climberx only uses fields that exist in ice,
        ! which have been transformed from each ice-sheet's unique
        ! field characteristics (naming, etc) into a standard 
        ! ice sheet object, which can then be transformed to 
        ! be used by other components. 

        ! This routine is designed to work on one domain at a time. 

        implicit none 

        type(ice_class),  intent(INOUT) :: ice   
        integer,          intent(IN)    :: idx_dom 
        character(len=*), intent(IN)    :: model
        integer,          intent(IN)    :: n_year_ice       ! number of years to integrate ice model    
        real(wp),         intent(IN)    :: time             ! current astronomical time in years
        logical,          intent(IN)    :: time_out_ice     ! Flag for whether to write 2D ice output
        
        return 

    end subroutine ice_update

    subroutine ice_init_domains(n_dom,model)
        ! Initialize the ice object with a predefined grid size 
        ! and the ice sheet object associated with it.

        implicit none 

        integer,          intent(IN)    :: n_dom 
        character(len=*), intent(IN)    :: model    
        
        return 

    end subroutine ice_init_domains

    subroutine ice_init(ice, idx_dom, model, time, l_restart, grid, cmn_grid, geo_grid, &
        z_bed_geo, z_bed_rel_geo, h_ice_geo, q_geo_geo, h_sed_geo)
        ! Initialize the ice object with a predefined grid size 
        ! and the ice sheet object associated with it.

        implicit none 

        type(ice_class),     intent(INOUT) :: ice
        integer,             intent(IN)    :: idx_dom 
        character(len=*),    intent(IN)    :: model     ! model name 
        real(wp),            intent(IN)    :: time 
        logical,             intent(IN)    :: l_restart 
        type(grid_class),    intent(IN)    :: grid      ! ice sheet grid
        type(grid_class),    intent(IN)    :: cmn_grid  ! coupler grid
        type(grid_class),    intent(IN)    :: geo_grid  ! geo grid
        real(wp),            intent(IN)    :: z_bed_geo(:,:)    ! bedrock elevation on geo grid [m]
        real(wp),            intent(IN)    :: z_bed_rel_geo(:,:)    ! icefree relaxed bedrock elevation on geo grid [m]
        real(wp),            intent(IN)    :: h_ice_geo(:,:)    ! ice thickness on geo grid [m]
        real(wp),            intent(IN)    :: q_geo_geo(:,:)    ! geothermal heat flux on geo grid [W/m2]
        real(wp),            intent(IN)    :: h_sed_geo(:,:)    ! sediment thickness on geo grid [m]
   
        return 

    end subroutine ice_init

    subroutine ice_end(ice,idx_dom,model,time)
        ! Finalize the ice object 

        implicit none 

        type(ice_class),  intent(INOUT) :: ice 
        integer,          intent(IN)    :: idx_dom
        character(len=*), intent(IN)    :: model  
        real(wp),         intent(IN)    :: time 

        return 

    end subroutine ice_end

    
    subroutine ice_write_restart(model,ice,idx_dom,time,restart_out_dir)
        ! write restart file

        implicit none 

        character(len=*), intent(IN) :: model    
        type(ice_class),  intent(IN) :: ice 
        integer,          intent(IN) :: idx_dom 
        real(wp),         intent(IN) :: time 
        character(len=*), intent(in) :: restart_out_dir
        
        return 

    end subroutine ice_write_restart

end module ice_model

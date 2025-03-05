!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : i c e _ m o d e l
!
!  Purpose : wrapper for ice sheet models
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Alexander Robinson and Matteo Willeit
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
module ice_model
    
    use precision, only : wp, sp, dp 
    use control, only : out_dir
    use coord, only : grid_class, grid_init
    use coord, only : map_scrip_class, map_scrip_init, map_scrip_field
    use ice_def, only : ice_class
    use yelmo, only : wp_yelmo, yelmo_class, yregions_class, yelmo_init_grid, yelmo_init, yelmo_init_state, &
                      yelmo_update, yelmo_end, yelmo_write_init, yelmo_write_reg_init, yelmo_write_reg_step, &
                      yelmo_write_step_model_metrics, yelmo_restart_write 
    use sicopolis, only : sico_class, sico_init, sico_update, sico_end, sico_write_restart
    use sico_out, only : sico_diag_init, sico_diag
    use ice_id_mod, only : set_ice_id
    use nml
    use ncio 

    implicit none 

    type(sico_class),  allocatable, target :: sico_doms(:)
    type(yelmo_class), allocatable, target :: ylmo_doms(:) 

    ! Save regional averaging calculations over [ndom,ntimes]
    ! to avoid writing netcdf output every year (which is slow)
    integer, parameter :: ylmo_reg_nt_tot = 100     ! [a] Write after nt_tot yrs
    integer,  allocatable :: ylmo_reg_nt_now(:) 
    real(wp), allocatable :: ylmo_reg_times(:,:)
    type(yregions_class), allocatable :: ylmo_regs(:,:) 
    
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
        integer,          intent(IN)    :: n_year_ice
        real(wp),         intent(IN)    :: time             ! current astronomical time in years
        logical,          intent(IN)    :: time_out_ice     ! Flag for whether to write 2D ice output
        
        ! Local variables 
        type(yelmo_class), pointer :: ylmo
        type(sico_class),  pointer :: sico 

        integer :: n, niter, k 
        character(len=1024) :: file1D 
        character(len=1024) :: file2D 
        character(len=1024) :: file_restart

        ! Update ice sheet 
        select case(trim(model))

            case("yelmo")
                ! Yelmo ice-sheet model interface 

                ! Get current domain from global domains
                ylmo => ylmo_doms(idx_dom) 

                ! First update ice sheet boundaries with ice data 
                call ice_to_yelmo(ylmo,ice) 

                ! Call ice sheet model 
                call yelmo_update(ylmo,real(time,wp_yelmo))

                ! Update ice with new ice-sheet information 
                call yelmo_to_ice(ice,ylmo)
                
                if (time_out_ice) then 
                    ! Write 2D output 

                    file2D = trim(out_dir)//"/ice_"//trim(ylmo%par%grid_name)//".nc"
                    call yelmo_write_step_2D(ylmo,file2D,time=real(time,sp))

                    !file_restart = trim(out_dir)//"/ice_restart_"//trim(ylmo%par%grid_name)//".nc"
                    !call yelmo_restart_write(ylmo,file_restart,time=real(time,sp)) 

                end if 

                ! Write 1D output (every year)
                !file1D = trim(out_dir)//"/ice_"//trim(ylmo%par%grid_name)//"_ts.nc"
                !call yelmo_write_reg_step(ylmo,file1D,time=real(time,sp))  

                ! Store regional calculations for this timestep in buffer
                ylmo_reg_nt_now(idx_dom) = ylmo_reg_nt_now(idx_dom)+1 
                ylmo_regs(idx_dom,ylmo_reg_nt_now)      = ylmo%reg 
                ylmo_reg_times(idx_dom,ylmo_reg_nt_now) = time 
                
                if (ylmo_reg_nt_now(idx_dom) .eq. ylmo_reg_nt_tot) then 
                    ! Buffer full, write output to file 

                    ! Write 1D output for each year in buffer
                    do k = 1, ylmo_reg_nt_tot
                        file1D = trim(out_dir)//"/ice_"//trim(ylmo%par%grid_name)//"_ts.nc"
                        call yelmo_write_reg_step(ylmo,file1D,time=real(ylmo_reg_times(idx_dom,k),wp_yelmo), &
                                                                            reg_now=ylmo_regs(idx_dom,k))
                    end do 

                    ! Reset total timesteps to zero 
                    ylmo_reg_nt_now(idx_dom) = 0 


                end if 

                nullify(ylmo) 

            case("sico")
                ! SICOPOLIS ice-sheet model interface 

                ! Get current domain from global domains
                sico => sico_doms(idx_dom)  ! fixme, causes segmentation fault
                
                ! First update ice sheet boundaries with ice data 
                call ice_to_sico(sico,ice) 
                
                niter = nint(n_year_ice/sico%timer%dtime_a)   ! integration steps/ice sheet call step
                ! Call ice sheet model 
                do n=1,niter
                  call sico_update(sico) 
                enddo
                call sico_diag(sico)

                ! Update ice with new ice-sheet information 
                call sico_to_ice(ice,sico)

                nullify(sico)
             
            case DEFAULT 

                write(*,*) "ice_init_domains:: Error: model not recognized."
                write(*,*) "model = ", trim(model) 
                stop 

        end select
        
        return 

    end subroutine ice_update

    subroutine ice_init_domains(n_dom,model)
        ! Initialize the ice object with a predefined grid size 
        ! and the ice sheet object associated with it.

        implicit none 

        integer,          intent(IN)    :: n_dom 
        character(len=*), intent(IN)    :: model    
        
        ! Allocate ice-sheet domains
        select case(trim(model))

            case("yelmo")
                ! Yelmo ice-sheet model interface 
                
                if (allocated(ylmo_doms)) deallocate(ylmo_doms)
                allocate(ylmo_doms(n_dom))

                ! Also allocate yreg object to store regional calcs 
                if (allocated(ylmo_regs)) deallocate(ylmo_regs) 
                if (allocated(ylmo_reg_nt_now)) deallocate(ylmo_reg_nt_now) 
                if (allocated(ylmo_reg_times))  deallocate(ylmo_reg_times) 
                allocate(ylmo_regs(n_dom,ylmo_reg_nt_tot))
                allocate(ylmo_reg_nt_now(n_dom))
                allocate(ylmo_reg_times(n_dom,ylmo_reg_nt_tot))

                ylmo_reg_nt_now = 0 
                ylmo_reg_times  = -9999.0 

            case("sico")
                ! SICOPOLIS ice-sheet model interface 

                if (allocated(sico_doms)) deallocate(sico_doms)
                allocate(sico_doms(n_dom))

            case DEFAULT 

                write(*,*) "ice_init_domains:: Error: model not recognized."
                write(*,*) "model = ", trim(model) 
                stop 

        end select
        
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
   
        ! Local variables 
        integer :: i, j, ii, jj
        real(wp) :: lon1, lon2, lat1, lat2, dlon_sur, dlat_sur
        real(wp), allocatable :: z_bed(:,:)    ! bedrock elevation on ice sheet grid [m]
        real(wp), allocatable :: z_bed_fil(:,:)    ! filtered bedrock elevation on ice sheet grid [m]
        real(wp), allocatable :: z_bed_rel(:,:)    ! icefree relaxed bedrock elevation on ice sheet grid [m]
        real(wp), allocatable :: h_ice(:,:)    ! ice thickness on ice sheet grid [m]
        real(wp), allocatable :: q_geo(:,:)    ! geothermal heat flux on ice sheet grid [W/m2]
        real(wp), allocatable :: h_sed(:,:)    ! sediment thickness on ice sheet grid [m]
        integer, allocatable :: id_mask(:,:)    ! ice id mask 
        type(map_scrip_class) :: maps_geo_to_ice

        type(yelmo_class), pointer :: ylmo
        type(sico_class),  pointer :: sico 

        ! Needed for Yelmo initialization 
        character(len=56)   :: domain 
        character(len=1024) :: par_path 
        character(len=1024) :: file_restart 
        character(len=1024) :: file1D 
        character(len=1024) :: file2D 

        ! Allocate ice object 
        call ice_alloc(ice,grid%G%nx,grid%G%ny)

        ! Store grid information 
        ice%grid = grid 

        ! derive correspondence between indexes on ice and coupler grids
        allocate(ice%grid_ice_to_cmn%i_lowres(ice%grid%G%nx,ice%grid%G%ny))
        allocate(ice%grid_ice_to_cmn%j_lowres(ice%grid%G%nx,ice%grid%G%ny))
        allocate(ice%grid_ice_to_cmn%ncells(cmn_grid%G%nx,cmn_grid%G%ny))
        dlon_sur = cmn_grid%lon(2,1)-cmn_grid%lon(1,1)
        dlat_sur = cmn_grid%lat(1,2)-cmn_grid%lat(1,1)
        ice%grid_ice_to_cmn%ncells = 0
        do i=1,cmn_grid%G%nx
          do j=1,cmn_grid%G%ny
            lon1 = cmn_grid%lon(i,j)-0.5_wp*dlon_sur
            if (i.eq.1) lon1 = lon1-1.e-3_wp
            lon2 = cmn_grid%lon(i,j)+0.5_wp*dlon_sur
            lat1 = cmn_grid%lat(i,j)-0.5_wp*dlat_sur
            if (j.eq.1) lat1 = lat1-1.e-3_wp
            lat2 = cmn_grid%lat(i,j)+0.5_wp*dlat_sur
            do ii=1,ice%grid%G%nx
              do jj=1,ice%grid%G%ny 
                if (ice%grid%lon(ii,jj).gt.lon1 .and. ice%grid%lon(ii,jj).le.lon2 &
                  .and. ice%grid%lat(ii,jj).gt.lat1 .and. ice%grid%lat(ii,jj).le.lat2) then
                  ice%grid_ice_to_cmn%i_lowres(ii,jj) = i
                  ice%grid_ice_to_cmn%j_lowres(ii,jj) = j
                  ice%grid_ice_to_cmn%ncells(i,j) = ice%grid_ice_to_cmn%ncells(i,j) + 1
                endif
              enddo
            enddo
          enddo
        enddo

        allocate(z_bed(grid%G%nx,grid%G%ny))
        allocate(z_bed_fil(grid%G%nx,grid%G%ny))
        allocate(z_bed_rel(grid%G%nx,grid%G%ny))
        allocate(h_ice(grid%G%nx,grid%G%ny))
        allocate(q_geo(grid%G%nx,grid%G%ny))
        allocate(h_sed(grid%G%nx,grid%G%ny))
        allocate(id_mask(grid%G%nx,grid%G%ny))
 
        ! set ice ID mask
        call set_ice_id(ice%grid, id_mask)

        call map_scrip_init(maps_geo_to_ice,geo_grid,ice%grid,method="con",fldr="maps",load=.TRUE.,clean=.FALSE.)
        call map_scrip_field(maps_geo_to_ice,"z_bed",z_bed_geo,z_bed,method="mean")
        call map_scrip_field(maps_geo_to_ice,"z_bed",z_bed_geo,z_bed_fil,method="mean", &
          filt_method="gaussian",filt_par=[100._dp,ice%grid%G%dx])
        call map_scrip_field(maps_geo_to_ice,"z_bed_rel",z_bed_rel_geo,z_bed_rel,method="mean")
        call map_scrip_field(maps_geo_to_ice,"h_ice",h_ice_geo,h_ice,method="mean")
        call map_scrip_field(maps_geo_to_ice,"q_geo",q_geo_geo,q_geo,method="mean", missing_value=-9999._dp, &
          filt_method="gaussian",filt_par=[100._dp,ice%grid%G%dx])
        call map_scrip_field(maps_geo_to_ice,"h_sed",h_sed_geo,h_sed,method="mean", missing_value=-9999._dp, &
          filt_method="none",filt_par=[100._dp,ice%grid%G%dx])

        where (h_ice<10._wp) h_ice = 0._wp

        ! Initialize ice sheet 
        select case(trim(model))

            case("yelmo")
                ! Yelmo ice-sheet model interface 

                ! Assign current domain to ylmo pointer
                ylmo => ylmo_doms(idx_dom)

                ! First define the Yelmo grid information from cmn%grid 
                ! Note: yelmo expects grid axes in [m], so use conversion factor.
                call yelmo_init_grid(grd=ylmo%grd,grid_name=grid%name, &
                        xc=real(grid%G%x*grid%xy_conv,wp_yelmo),yc=real(grid%G%y*grid%xy_conv,wp_yelmo), &
                        lon=real(grid%lon,wp_yelmo),lat=real(grid%lat,wp_yelmo),area=real(grid%area,wp_yelmo))

                ! Define the parameter file we expect for this simulation
                ! Note: currently the same parameters are used for all domains. 
                !par_path = trim(out_dir)//"/yelmo_"//trim(grid%name)//".nml" 
                par_path = trim(out_dir)//"/ice_yelmo_par.nml" 

                ! For now, set domain equal to the grid name
                ! (usually this should be "North", "Greenland", "Antarctica", etc.) 
                domain = trim(grid%name) 

                ! Now initialize remaining Yelmo components
                ! (set `grid_def="none"` to use the grid definition from above)
                call yelmo_init(ylmo,par_path,grid_def="none",time=real(time,wp_yelmo), &
                            load_topo=.FALSE.,domain=domain,grid_name=grid%name)

                ! Define topography based on input data
                ylmo%tpo%now%H_ice   = h_ice 
                ylmo%bnd%z_bed       = z_bed 
                ylmo%bnd%z_bed_sd    = 0.0_wp 
                
                ! Populate dummy boundary information
                ylmo%bnd%z_sl        = 0.0_wp 
                ylmo%bnd%Q_geo       = q_geo*1e3  ! [W/m2] => mW/m2]
                ylmo%bnd%T_srf       = 240.0_wp 
                ylmo%bnd%T_shlf      = 273.15_wp 
                ylmo%bnd%smb         = 0.0_wp 
                ylmo%bnd%bmb_shlf    = 0.0_wp 
                ylmo%bnd%H_sed       = H_sed
                
                ! Allow ice to grow over the whole domain by default
                ! (later this will be set by ice%mask_extent)
                ylmo%bnd%ice_allowed = .TRUE. 

                ! Set reference ice thickness and bed elevation to initial one
                ! (only useful for comparisons, probably not used in climber) 
                ylmo%bnd%H_ice_ref   = ylmo%tpo%now%H_ice 
                ylmo%bnd%z_bed_ref   = ylmo%bnd%z_bed 


                ! Finish state initialization of all Yelmo components
                call yelmo_init_state(ylmo,real(time,wp_yelmo),thrm_method="robin-cold")

                ! Done: yelmo domain is now ready to run. 

                ! Fill ice with ice-sheet information 
                call yelmo_to_ice(ice,ylmo) 


                ! Initialize 2D output file 
                file2D = trim(out_dir)//"/ice_"//trim(ylmo%par%grid_name)//".nc"

                call yelmo_write_init(ylmo,file2D,time_init=real(time,wp_yelmo),units="years")  
                call yelmo_write_step_2D(ylmo,file2D,time=real(time,sp))
                
                ! Initialize 1D output file 
                file1D = trim(out_dir)//"/ice_"//trim(ylmo%par%grid_name)//"_ts.nc"

                call yelmo_write_reg_init(ylmo,file1D,time_init=real(time,wp_yelmo),units="years", &
                                                            mask=ylmo%bnd%ice_allowed)
                call yelmo_write_reg_step(ylmo,file1D,time=real(time,wp_yelmo))  

                ! Disassociate the local pointer ylmo 
                nullify(ylmo) 

            case("sico")
                ! SICOPOLIS ice-sheet model interface 
                
                ! Assign pointer to current ice domain 
                sico => sico_doms(idx_dom)

                sico%grid%grid1 = grid 
                call sico_init(sico, grid, l_restart, &
                  real(z_bed,dp), real(z_bed_fil,dp), real(z_bed_rel,dp), real(h_ice,dp), real(q_geo,dp), real(h_sed,dp), id_mask)

                call sico_diag_init(sico)

                ! fill ice with ice-sheet information 
                call sico_to_ice(ice,sico)

                ! Disassociate the local pointer sico
                nullify(sico)
             
            case DEFAULT 

                write(*,*) "ice_init:: Error: model not recognized."
                write(*,*) "model = ", trim(model) 
                stop 

        end select

        deallocate(z_bed)
        deallocate(z_bed_fil)
        deallocate(z_bed_rel)
        deallocate(h_ice)
        deallocate(q_geo)
        deallocate(h_sed)
        deallocate(id_mask)

        return 

    end subroutine ice_init

    subroutine ice_end(ice,idx_dom,model,time)
        ! Finalize the ice object 

        implicit none 

        type(ice_class),  intent(INOUT) :: ice 
        integer,          intent(IN)    :: idx_dom
        character(len=*), intent(IN)    :: model  
        real(wp),         intent(IN)    :: time 

        ! Deallocate ice object
        call ice_dealloc(ice)

        ! Terminate ice sheet instance 
        select case(trim(model))

            case("yelmo")
                ! Yelmo ice-sheet model interface 

                call yelmo_end(ylmo_doms(idx_dom),time=real(time,wp_yelmo))

            case("sico")
                ! SICOPOLIS ice-sheet model interface 

                call sico_end(sico_doms(idx_dom)) 

            case DEFAULT 

                write(*,*) "ice_end:: Error: model not recognized."
                write(*,*) "model = ", trim(model) 
                stop 

        end select

        return 

    end subroutine ice_end

    subroutine yelmo_to_ice(ice,ylmo)
        ! Populate ice fields with Yelmo information 

        implicit none 

        type(ice_class),   intent(INOUT) :: ice 
        type(yelmo_class), intent(IN)    :: ylmo 

        ! mask of maximum allowed ice extent
        ! Note: Yelmo expects this as a boundary condition,
        ! so providing it here is a bit circular. But it is done 
        ! for consistency for now.
        ice%mask_extent = 0 
        where(ylmo%bnd%ice_allowed) ice%mask_extent =  1 

        ! ice thickness 
        ice%H_ice = ylmo%tpo%now%H_ice 

        ! surface elevation
        ice%z_sur = ylmo%tpo%now%z_srf 

        ! ice base elevation 
        ice%z_base = ylmo%tpo%now%z_base 

        ! bedrock elevation 
        ice%z_bed = ylmo%bnd%z_bed

        ! basal melt [m(ice)/yr] => [m(ice)/s]
        ! (negativ bc yelmo convention is basal mass balance with melt negative)
        ice%Q_b = -ylmo%tpo%now%bmb*ylmo%tpo%now%f_ice/ylmo%bnd%c%sec_year 

        ! calving \[m(ice)/yr] => \[m(ice)/s]
        ! (negative bc yelmo convention is 'calving mass balance' with loss negative)
        ! include here both explicit calving (cmb) and subgrid discharge (dmb)
        ice%calv = -(ylmo%tpo%now%cmb+ylmo%tpo%now%dmb)/ylmo%bnd%c%sec_year
        
        return 

    end subroutine yelmo_to_ice

    subroutine ice_to_yelmo(ylmo,ice)
        ! Populate Yelmo fields with ice information 

        implicit none 

        type(yelmo_class), intent(INOUT) :: ylmo 
        type(ice_class),   intent(IN)    :: ice 
        
        ! sea level, always zero in climber [m]
        ylmo%bnd%z_sl  = 0.0_wp 

        ! bedrock elevation [m]
        ylmo%bnd%z_bed = ice%z_bed 
        ylmo%bnd%z_bed_sd = ice%z_bed_std

        ! annual surface mass balance [m(ice equivalent)/s] => [m(ice equivalent)/yr]
        ylmo%bnd%smb = ice%smb*ylmo%bnd%c%sec_year 

        ! basal mass balance of floating ice [m(ice equivalent)/s] => [! m(ice equivalent)/yr]
        ylmo%bnd%bmb_shlf = ice%bmb_shlf*ylmo%bnd%c%sec_year 

        ! ice surface temperature [degC] => [K] 
        ylmo%bnd%T_srf = min(0.0_wp,ice%temp_s) + 273.15_wp

        ! ground temperature (shelf temperature?) [degC] => [K]
        ylmo%bnd%T_shlf = ice%temp_g + 273.15_wp

        ! geothermal heat flux [W/m2] => mW/m2]
        ylmo%bnd%Q_geo = ice%q_geo*1e3 

        ! sediment thickness 
        ylmo%bnd%H_sed = ice%H_sed

        return 

    end subroutine ice_to_yelmo

    subroutine sico_to_ice(ice,sico)
        ! Populate ice fields with SICOPOLIS information 

        implicit none 

        type(ice_class),  intent(INOUT) :: ice 
        type(sico_class), intent(IN)    :: sico 

        integer :: i, j

        ! Note: SICOPOLIS grid is (j,i) and indexes start from 0

        ! error flag
        ice%error = sico%error

        do i=1,ice%grid%G%nx
          do j=1,ice%grid%G%ny
            ! mask of maximum allowed ice extent
            ice%mask_extent(i,j) = sico%state%mask_maxextent(j-1,i-1) 
            ! ice thickness 
            ice%H_ice(i,j) = sico%state%H_c(j-1,i-1) + sico%state%H_t(j-1,i-1)     ! m
            ! surface elevation
            ice%z_sur(i,j) = sico%state%zs(j-1,i-1)     ! m
            ! ice base elevation
            ice%z_base(i,j) = sico%state%zb(j-1,i-1)    ! m
            ! bedrock elevation
            ice%z_bed(i,j) = sico%state%zl(j-1,i-1)    ! m
            ! basal melt 
            ice%Q_b(i,j) = sico%state%Q_b_apl(j-1,i-1)     ! m(ice)/s
            ! calving
            ice%calv(i,j) = sico%state%calving_apl(j-1,i-1) ! m(ice)/s
          enddo
        enddo

        return 

    end subroutine sico_to_ice

    subroutine ice_to_sico(sico,ice)
        ! Populate SICOPOLIS fields with ice information 

        implicit none 

        type(sico_class), intent(INOUT) :: sico 
        type(ice_class),  intent(IN)    :: ice 
        
        integer :: i, j, ii, jj, i_f, j_f, nx, ny, n_filter
        real(wp) :: dx, dy, sigma_ref, sigma_filter, dist, sum_weigh, weigh
        real(wp), dimension(:,:), allocatable :: zl

        ! Note: SICOPOLIS grid is (j,i) and indexes start from 0

        do i=1,ice%grid%G%nx
          do j=1,ice%grid%G%ny
            ! sea level
            sico%state%z_sl(j-1,i-1) = ice%z_sl(i,j)     ! m
            ! bedrock elevation
            sico%state%zl_neu(j-1,i-1) = ice%z_bed(i,j)     ! m
            ! bedrock elevation, filtered
            sico%state%zl_fil(j-1,i-1) = ice%z_bed_fil(i,j)     ! m
            ! standard deviation of surface elevation
            sico%state%zs_std(j-1,i-1) = ice%z_sur_std(i,j)     ! m
            ! standard deviation of bedrock topography
            sico%state%zl_std(j-1,i-1) = ice%z_bed_std(i,j)     ! m
            ! annual surface mass balance
            sico%state%as_perp(j-1,i-1) = ice%smb(i,j)      ! m(ice equivalent)/s
            ! annual accumulation (including liquid precipitation)
            sico%state%accum(j-1,i-1)  = ice%accum(i,j)     ! m(ice equivalent)/s
            ! annual runoff
            sico%state%runoff(j-1,i-1) = ice%runoff(i,j)    ! m(ice equivalent)/s
            ! basal melt of floating ice 
            sico%state%Q_bm_float(j-1,i-1) = -ice%bmb_shlf(i,j) ! m(ice equivalent)/s
            ! ice surface temperature
            sico%state%temp_s(j-1,i-1) = min(-0.001_wp,ice%temp_s(i,j)) ! degC
            ! ground temperature
            sico%state%temp_g(j-1,i-1) = ice%temp_g(i,j)    ! degC
            ! geothermal heat flux
            sico%state%q_geo(j-1,i-1) = ice%q_geo(i,j)      ! W/m2
            ! sediment thickness 
            sico%state%H_sed(j-1,i-1) = ice%H_sed(i,j)      ! m
            ! ocean temperature
            sico%state%t_ocn(j-1,i-1) = ice%t_ocn(i,j)    ! degC
            ! ocean salinity
            sico%state%s_ocn(j-1,i-1) = ice%s_ocn(i,j)    ! psu
          enddo
        enddo

      ! filter lithosphere topography for ice sheet model
      if (sico%par%l_smooth_zl > 0) then

        nx = sico%grid%grid1%G%nx
        ny = sico%grid%grid1%G%ny
        allocate(zl(0:ny-1,0:nx-1))
        zl = sico%state%zl_neu
!        sico%state%zl_neu = 0._wp
        dx = sico%grid%grid1%G%dx ! km
        dy = sico%grid%grid1%G%ny ! km

        if(sico%par%l_smooth_zl == 1) then
          sigma_filter = sico%par%sigma_filter_zl/dx   ! half span of filtered area, in grid points
          n_filter     = ceiling(2.0_wp*sigma_filter)
          do i=0,nx-1 
          do j=0,ny-1
            sum_weigh = 0.0_wp
            sico%state%zl_neu(j,i) = 0._wp
            do ii=-n_filter, n_filter
            do jj=-n_filter, n_filter
              i_f = i+ii
              j_f = j+jj
              if (i_f <  0) i_f = 0
              if (i_f > nx-1) i_f = nx-1
              if (j_f <  0) j_f = 0
              if (j_f > ny-1) j_f = ny-1
              dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
              weigh     = exp(-(dist/sigma_filter)**2)
              sum_weigh = sum_weigh + weigh
              sico%state%zl_neu(j,i) = sico%state%zl_neu(j,i) + weigh*zl(j_f,i_f)
            end do
            end do
            sico%state%zl_neu(j,i) = sico%state%zl_neu(j,i)/sum_weigh
          end do
          end do
        else if(sico%par%l_smooth_zl == 2) then
          sigma_ref = sico%par%sigma_filter_zl/dx
          do i=0,nx-1 
          do j=0,ny-1
            ! half span of filtered area, in grid points
            sigma_filter = min(max(sigma_ref, sigma_ref-sico%par%s_filter_zl* &
                                 (zl(j,i)-sico%par%zl_cont)), 500.0_wp/dx)
            n_filter     = ceiling(2.0_wp*sigma_filter)
            sum_weigh = 0.0_wp
            sico%state%zl_neu(j,i) = 0._wp
            do ii=-n_filter, n_filter
            do jj=-n_filter, n_filter
              i_f = i+ii
              j_f = j+jj
              if (i_f <  0) i_f = 0
              if (i_f > nx-1) i_f = nx-1
              if (j_f <  0) j_f = 0
              if (j_f > ny-1) j_f = ny-1
              dist      = sqrt(real(ii,wp)**2+real(jj,wp)**2)
              weigh     = exp(-(dist/sigma_filter)**2)
              sum_weigh = sum_weigh + weigh
              sico%state%zl_neu(j,i) = sico%state%zl_neu(j,i) + weigh*zl(j_f,i_f)
            end do
            end do
            sico%state%zl_neu(j,i) = sico%state%zl_neu(j,i)/sum_weigh
          end do
          end do
        end if
        deallocate(zl)
      endif

      ! rate of change of lithosphere elevation, needed? fixme
      sico%state%dzl_dtau = (sico%state%zl_neu - sico%state%zl) * sico%timer%dtime_inv

        return 

    end subroutine ice_to_sico
    
    subroutine ice_write_restart(model,ice,idx_dom,time,restart_out_dir)
        ! write restart file

        implicit none 

        character(len=*), intent(IN) :: model    
        type(ice_class),  intent(IN) :: ice 
        integer,          intent(IN) :: idx_dom 
        real(wp),         intent(IN) :: time 
        character(len=*), intent(in) :: restart_out_dir
        
        ! Local variables 
        character(len=1024) :: file_restart 

        ! Allocate ice-sheet domains
        select case(trim(model))

            case("yelmo")
                ! Yelmo ice-sheet model interface 

                file_restart = trim(restart_out_dir)//"/ice_yelmo_"//trim(ice%grid%name)//"_restart.nc"

                call yelmo_restart_write(ylmo_doms(idx_dom),file_restart,time=real(time,wp_yelmo))

            case("sico")
                ! SICOPOLIS ice-sheet model interface 

                file_restart = trim(restart_out_dir)//"/ice_sico_"//trim(ice%grid%name)//"_restart.nc"
                call sico_write_restart(file_restart,sico_doms(idx_dom))

            case DEFAULT 

                write(*,*) "ice_write_restart:: Error: model not recognized."
                write(*,*) "model = ", trim(model) 
                stop 

        end select
        
        return 

    end subroutine ice_write_restart

    subroutine ice_alloc(ice,nx,ny)
        ! Allocate ice fields 

        implicit none 

        type(ice_class), intent(INOUT) :: ice 
        integer, intent(IN) :: nx 
        integer, intent(IN) :: ny 
        
        ! First ensure all fields are deallocated
        call ice_dealloc(ice) 

        ! Allocate fields 
        allocate(ice%z_sur(nx,ny))
        allocate(ice%z_sur_std(nx,ny))
        allocate(ice%z_bed_std(nx,ny))
        allocate(ice%z_base(nx,ny))
        allocate(ice%z_bed(nx,ny))
        allocate(ice%z_bed_fil(nx,ny))
        allocate(ice%H_ice(nx,ny))
        allocate(ice%mask_ocn_lake(nx,ny))
        allocate(ice%z_sl(nx,ny))
        allocate(ice%mask_extent(nx,ny))
        allocate(ice%calv(nx,ny))
        allocate(ice%Q_b(nx,ny))
        allocate(ice%bmb_shlf(nx,ny))
        allocate(ice%smb(nx,ny))
        allocate(ice%accum(nx,ny))
        allocate(ice%runoff(nx,ny))
        allocate(ice%temp_s(nx,ny))
        allocate(ice%temp_g(nx,ny))
        allocate(ice%q_geo(nx,ny))
        allocate(ice%H_sed(nx,ny))
        allocate(ice%t_ocn(nx,ny))
        allocate(ice%s_ocn(nx,ny))
        
        ! Initialize to zero 
        ice%z_sur       = 0.0_wp
        ice%z_sur_std   = 0.0_wp
        ice%z_bed_std   = 0.0_wp
        ice%z_base      = 0.0_wp
        ice%z_bed       = 0.0_wp
        ice%z_bed_fil   = 0.0_wp
        ice%H_ice       = 0.0_wp
        ice%z_sl        = 0.0_wp
        ice%mask_extent = 0.0_wp
        ice%calv        = 0.0_wp
        ice%Q_b         = 0.0_wp
        ice%bmb_shlf    = 0.0_wp
        ice%smb         = 0.0_wp
        ice%accum       = 0.0_wp
        ice%runoff      = 0.0_wp
        ice%temp_s      = 0.0_wp
        ice%temp_g      = 0.0_wp
        ice%q_geo       = 0.0_wp
        ice%H_sed       = 0.0_wp
        ice%t_ocn       = 0.0_wp
        ice%s_ocn       = 0.0_wp

        return 

    end subroutine ice_alloc
    
    subroutine ice_dealloc(ice)
        ! Deallocate ice fields 

        implicit none 

        type(ice_class), intent(INOUT) :: ice 
        
        if (allocated(ice%z_sur))       deallocate(ice%z_sur)
        if (allocated(ice%z_sur_std))   deallocate(ice%z_sur_std)
        if (allocated(ice%z_bed_std))   deallocate(ice%z_bed_std)
        if (allocated(ice%z_base))      deallocate(ice%z_base)
        if (allocated(ice%z_bed))       deallocate(ice%z_bed)
        if (allocated(ice%z_bed_fil))   deallocate(ice%z_bed_fil)
        if (allocated(ice%H_ice))       deallocate(ice%H_ice)
        if (allocated(ice%mask_ocn_lake))deallocate(ice%mask_ocn_lake)
        if (allocated(ice%z_sl))        deallocate(ice%z_sl)
        if (allocated(ice%mask_extent)) deallocate(ice%mask_extent)
        if (allocated(ice%calv))        deallocate(ice%calv)
        if (allocated(ice%Q_b))         deallocate(ice%Q_b)
        if (allocated(ice%bmb_shlf))    deallocate(ice%bmb_shlf)
        if (allocated(ice%smb))         deallocate(ice%smb)
        if (allocated(ice%accum))       deallocate(ice%accum)
        if (allocated(ice%runoff))      deallocate(ice%runoff)
        if (allocated(ice%temp_s))      deallocate(ice%temp_s)
        if (allocated(ice%temp_g))      deallocate(ice%temp_g)
        if (allocated(ice%q_geo))       deallocate(ice%q_geo)
        if (allocated(ice%H_sed))       deallocate(ice%H_sed)
        if (allocated(ice%t_ocn))       deallocate(ice%t_ocn)
        if (allocated(ice%s_ocn))       deallocate(ice%s_ocn)
        
        return 

    end subroutine ice_dealloc


! === Yelmo 2D output routine =======

    subroutine yelmo_write_step_2D(ylmo,filename,time)

        implicit none 
        
        type(yelmo_class), intent(IN) :: ylmo
        character(len=*),  intent(IN) :: filename
        real(sp), intent(IN) :: time

        ! Local variables
        integer    :: ncid, n
        real(sp)   :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename,"time",ncid)
        call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
        if (abs(time-time_prev).gt.1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)

        ! Write model metrics (model speed, dt, eta)
        call yelmo_write_step_model_metrics(filename,ylmo,n,ncid)

        if (n .eq. 1) then 

            ! Write some constant fields 

            call nc_write(filename,"Q_geo",ylmo%bnd%Q_geo,units="mW/m^2",long_name="Geothermal heat flux", &
                      dim1="xc",dim2="yc",ncid=ncid)
        
        end if 

        ! == Boundaries == 

        call nc_write(filename,"z_bed",ylmo%bnd%z_bed,units="m",long_name="Bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"z_bed_sd",ylmo%bnd%z_bed_sd,units="m",long_name="Sub-grid standard deviation of bedrock elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid) 
        call nc_write(filename,"H_sed",ylmo%bnd%H_sed,units="m",long_name="Sediment thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"smb",ylmo%bnd%smb,units="m/a ice equiv.",long_name="Surface mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"bmb_shlf",ylmo%bnd%bmb_shlf,units="m/a ice equiv.",long_name="Shelf basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! == Topo == 

        call nc_write(filename,"H_ice",ylmo%tpo%now%H_ice,units="m",long_name="Ice thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"z_srf",ylmo%tpo%now%z_srf,units="m",long_name="Surface elevation", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mask_bed",ylmo%tpo%now%mask_bed,units="",long_name="Bed mask", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dist_grline",ylmo%tpo%now%dist_grline,units="m",long_name="Distance to grounding line", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_grnd",ylmo%tpo%now%f_grnd,units="1",long_name="Grounded ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_ice",ylmo%tpo%now%f_ice,units="1",long_name="Total ice fraction", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"cmb",ylmo%tpo%now%cmb,units="m/a ice equiv.",long_name="Calving mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dmb",ylmo%tpo%now%dmb,units="m/a ice equiv.",long_name="Discharge mass balance (subgrid)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_net",ylmo%tpo%now%mb_net,units="m",long_name="Net mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"mb_resid",ylmo%tpo%now%mb_resid,units="m",long_name="Residual mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"bmb",ylmo%tpo%now%bmb,units="m/a ice equiv.",long_name="Basal mass balance", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"dHidt",ylmo%tpo%now%dHidt,units="m/a",long_name="Ice thickness change", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"H_calv",ylmo%tpo%now%H_calv,units="m",long_name="Threshold ice thickness for calving", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)

        ! == Dynamics == 

        call nc_write(filename,"N_eff",ylmo%dyn%now%N_eff,units="Pa",long_name="Effective pressure", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"c_bed",ylmo%dyn%now%c_bed,units="Pa",long_name="Bed friction coefficient", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"beta",ylmo%dyn%now%beta,units="Pa a m^-1",long_name="Basal friction coefficient", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"visc_eff_int",ylmo%dyn%now%visc_eff_int,units="Pa a m",long_name="Depth-integrated effective viscosity (SSA)", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taud",ylmo%dyn%now%taud,units="Pa",long_name="Driving stress", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"taub",ylmo%dyn%now%taub,units="Pa",long_name="Basal dragging stress (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"enh_bar",ylmo%mat%now%enh_bar,units="1",long_name="Vertically averaged enhancement factor", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! call nc_write(filename,"ux_s",ylmo%dyn%now%ux_s,units="m/a",long_name="Surface velocity (x)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"uy_s",ylmo%dyn%now%uy_s,units="m/a",long_name="Surface velocity (y)", &
        !                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"uxy_s",ylmo%dyn%now%uxy_s,units="m/a",long_name="Surface velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"uxy_b",ylmo%dyn%now%uxy_b,units="m/a",long_name="Basal velocity (magnitude)", &
                       dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        ! call nc_write(filename,"uxy_bar",ylmo%dyn%now%uxy_bar,units="m/a",long_name="Vertically averaged velocity magnitude", &
        !               dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        
        ! call nc_write(filename,"uz",ylmo%dyn%now%uz,units="m/a",long_name="3D velocity (z)", &
        !                dim1="xc",dim2="yc",dim3="zeta_ac",dim4="time",start=[1,1,1,n],ncid=ncid)
        
        ! == Thermo == 
        
        call nc_write(filename,"Q_ice_b",ylmo%thrm%now%Q_ice_b,units="mW m^-2",long_name="Basal ice heat flow", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_rock",ylmo%thrm%now%Q_rock,units="mW m^-2",long_name="Bedrock surface heat flow", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"Q_b",ylmo%thrm%now%Q_b,units="mW m^-2",long_name="Basal frictional heating", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        call nc_write(filename,"T_ice",ylmo%thrm%now%T_ice,units="K",long_name="Ice temperature", &
                      dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        ! call nc_write(filename,"T_prime",ylmo%thrm%now%T_ice-ylmo%thrm%now%T_pmp,units="deg C",long_name="Homologous ice temperature", &
        !               dim1="xc",dim2="yc",dim3="zeta",dim4="time",start=[1,1,1,n],ncid=ncid)
        call nc_write(filename,"T_prime_b",ylmo%thrm%now%T_prime_b,units="K",long_name="Basal homologous ice temperature", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"f_pmp",ylmo%thrm%now%f_pmp,units="1",long_name="Fraction of grid point at pmp", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        call nc_write(filename,"H_w",ylmo%thrm%now%H_w,units="m water equiv.",long_name="Basal water layer thickness", &
                      dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)
        
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine yelmo_write_step_2D

end module ice_model

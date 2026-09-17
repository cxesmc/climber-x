module lndvc_params
    ! Namelist-loaded parameters for the lndvc output layer.
    !
    ! Reads nml/lndvc_par.nml (staged into out_dir at run-time by runme), the
    ! lndvc analog of nml/lnd_par.nml. Keeps output selectors + hires field
    ! toggles in one place so lndvc_out (and the future lndvc_restart) don't
    ! duplicate the namelist plumbing.

    use precision, only : wp
    use nml
    use control,   only : out_dir

    implicit none

    ! ---- output selectors --------------------------------------------------
    logical :: write_surf   = .true.
    logical :: write_hires  = .true.
    logical :: write_cons   = .false.

    logical :: l_monthly_output = .true.
    logical :: l_daily_output   = .false.

    ! ---- hires field subset ------------------------------------------------
    logical :: hires_write_h_snow  = .true.
    logical :: hires_write_w_snow  = .true.
    logical :: hires_write_t_skin  = .true.
    logical :: hires_write_smb     = .true.
    logical :: hires_write_melt    = .true.
    logical :: hires_write_f_ice   = .true.
    logical :: hires_write_veg_c   = .false.
    logical :: hires_write_soil_c  = .false.

    ! ---- restart control ---------------------------------------------------
    integer :: i_write_restart_lndvc = 1
    logical :: lndvc_restart         = .false.

    public

contains

    subroutine lndvc_par_init
        ! Read lndvc_par.nml into the module variables above. Mirrors the
        ! reference smb_par_init / lnd_params_init idiom: parameterless,
        ! filename from control::out_dir. Safe to call multiple times.
        implicit none

        character(len=256) :: filename

        filename = trim(out_dir) // "/lndvc_par.nml"

        call nml_read(filename,"lndvc_par","write_surf",            write_surf)
        call nml_read(filename,"lndvc_par","write_hires",           write_hires)
        call nml_read(filename,"lndvc_par","write_cons",            write_cons)
        call nml_read(filename,"lndvc_par","l_monthly_output",      l_monthly_output)
        call nml_read(filename,"lndvc_par","l_daily_output",        l_daily_output)

        call nml_read(filename,"lndvc_par","hires_write_h_snow",    hires_write_h_snow)
        call nml_read(filename,"lndvc_par","hires_write_w_snow",    hires_write_w_snow)
        call nml_read(filename,"lndvc_par","hires_write_t_skin",    hires_write_t_skin)
        call nml_read(filename,"lndvc_par","hires_write_smb",       hires_write_smb)
        call nml_read(filename,"lndvc_par","hires_write_melt",      hires_write_melt)
        call nml_read(filename,"lndvc_par","hires_write_f_ice",     hires_write_f_ice)
        call nml_read(filename,"lndvc_par","hires_write_veg_c",     hires_write_veg_c)
        call nml_read(filename,"lndvc_par","hires_write_soil_c",    hires_write_soil_c)

        call nml_read(filename,"lndvc_par","i_write_restart_lndvc", i_write_restart_lndvc)
        call nml_read(filename,"lndvc_par","lndvc_restart",         lndvc_restart)

        return
    end subroutine lndvc_par_init

end module lndvc_params

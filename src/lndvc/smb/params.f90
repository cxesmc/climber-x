module smb_par_m
    ! Parameters for the ported SMB single-column physics.
    !
    ! Brought into the lndvc framework so it is self-contained (does not depend
    ! on src/smb). Values mirror src/smb/smb_params (the "smb_par" namelist
    ! group) so behaviour matches the reference SMB model. Expanded as more
    ! physics files are ported.

    use precision, only : wp
    use nml

    implicit none

    ! time step [s] (set in smb_par_init; default 1 day)
    real(wp) :: dt  = 86400._wp
    real(wp) :: rdt = 1._wp/86400._wp

    ! water-budget debug switch (mirrors control::check_water)
    logical  :: check_water = .false.

    type snow_par_type
        logical  :: lsnow_dust
        logical  :: lsnow_aging
        integer  :: isnow_albedo
        real(wp) :: dalb_snow_vis
        real(wp) :: dalb_snow_nir
        real(wp) :: w_snow_crit = 5._wp        ! kg/m2, min SWE for explicit snow layer
        real(wp) :: snow_grain_fresh           ! um, fresh snow grain size
        real(wp) :: snow_grain_old = 1000._wp  ! um, old snow grain size
        real(wp) :: f_age_t                    ! temperature parameter for snow aging
        real(wp) :: dT_age
        real(wp) :: snow_0                      ! kg/m2/s, critical snowfall rate for aging
        real(wp) :: snow_1                      ! kg/m2/s, minimum snowfall rate for aging
        real(wp) :: k_sigma_orog
        real(wp) :: sigma_orog_crit
        real(wp) :: c_fsnow
        real(wp) :: c_fsnow_orog
        logical  :: l_fsnow_orog
        real(wp) :: rho
        real(wp) :: lambda
        real(wp) :: w_snow_max  = 1000._wp     ! kg/m2, maximum SWE
        real(wp) :: w_snow_dust
        real(wp) :: dust_con_scale
        integer  :: i_rfz
        real(wp) :: porosity
        real(wp) :: f_rfz_max                   ! max fraction of refreezing
        real(wp) :: f_rfz_to_snow_max           ! fraction of refreezing into snow
        real(wp) :: wsnow_crit_rfz
    end type
    type(snow_par_type) :: snow_par

contains

    subroutine smb_par_init(filename, dt_smb)
        ! Read the snow parameters from the "smb_par" namelist group so values
        ! match the reference SMB model. dt is set from the passed time step.

        implicit none

        character(len=*), intent(in)           :: filename
        real(wp),         intent(in), optional :: dt_smb

        if (present(dt_smb)) dt = dt_smb
        rdt = 1._wp/dt

        call nml_read(filename,"smb_par","rho_snow",          snow_par%rho)
        call nml_read(filename,"smb_par","lambda_snow",       snow_par%lambda)
        call nml_read(filename,"smb_par","lsnow_aging",       snow_par%lsnow_aging)
        call nml_read(filename,"smb_par","lsnow_dust",        snow_par%lsnow_dust)
        call nml_read(filename,"smb_par","f_age_t",           snow_par%f_age_t)
        call nml_read(filename,"smb_par","dT_age",            snow_par%dT_age)
        call nml_read(filename,"smb_par","snow_0",            snow_par%snow_0)
        call nml_read(filename,"smb_par","snow_1",            snow_par%snow_1)
        call nml_read(filename,"smb_par","snow_grain_fresh",  snow_par%snow_grain_fresh)
        call nml_read(filename,"smb_par","isnow_albedo",      snow_par%isnow_albedo)
        call nml_read(filename,"smb_par","dalb_snow_vis",     snow_par%dalb_snow_vis)
        call nml_read(filename,"smb_par","dalb_snow_nir",     snow_par%dalb_snow_nir)
        call nml_read(filename,"smb_par","k_sigma_orog",      snow_par%k_sigma_orog)
        call nml_read(filename,"smb_par","sigma_orog_crit",   snow_par%sigma_orog_crit)
        call nml_read(filename,"smb_par","c_fsnow",           snow_par%c_fsnow)
        call nml_read(filename,"smb_par","c_fsnow_orog",      snow_par%c_fsnow_orog)
        call nml_read(filename,"smb_par","l_fsnow_orog",      snow_par%l_fsnow_orog)
        call nml_read(filename,"smb_par","w_snow_max",        snow_par%w_snow_max)
        call nml_read(filename,"smb_par","w_snow_dust",       snow_par%w_snow_dust)
        call nml_read(filename,"smb_par","dust_con_scale",    snow_par%dust_con_scale)
        call nml_read(filename,"smb_par","i_rfz",             snow_par%i_rfz)
        call nml_read(filename,"smb_par","porosity",          snow_par%porosity)
        call nml_read(filename,"smb_par","f_rfz_max",         snow_par%f_rfz_max)
        call nml_read(filename,"smb_par","f_rfz_to_snow_max", snow_par%f_rfz_to_snow_max)
        call nml_read(filename,"smb_par","wsnow_crit_rfz",    snow_par%wsnow_crit_rfz)

        return

    end subroutine smb_par_init

end module smb_par_m

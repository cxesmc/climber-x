!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!
!  Module : t i m e r
!
!  Purpose : model timer 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2022 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit
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
module timer

  use precision, only : wp, dp
  use nml
  use control, only : flag_atm, flag_ocn, flag_lnd, flag_sic, flag_ice, flag_smb, flag_bmb, flag_bgc, flag_geo, ifake_geo
  use control, only : l_spinup_cc, nyears_spinup_bgc, year_start_offline, nyear_avg_offline
  use control, only : i_write_restart, n_year_write_restart, years_write_restart

  implicit none

! doy: day of year
! year: ordinary year
! year_clim:  climate year (reduced due to acceleration)
! year_smb:  smb year (reduced due to acceleration times smb call)
  integer :: soy, doy, mon, year, year_geo, year_clim, year_smb, nyears, nyears_clim, nyears_geo, &
  year_ini, year_now
  integer :: step, nstep, nstep_day, nstep_mon, nstep_year
  integer :: step_atm, step_lnd, step_ocn, step_sic, step_bgc, step_smb, step_bmb, step_ice

  integer, parameter :: nday_year = 360  !! days per year [day/year]
  integer, parameter :: nmon_year = 12   !! months per year [mon/year]
  integer, parameter :: nday_mon  = 30   !! days per month [day/mon]

  real(wp), parameter :: day_year = dble(nday_year) !! days per year [day/year] 
  real(wp), parameter :: day_mon  = dble(nday_mon)  !! months per year [mon/year] 
  real(wp), parameter :: mon_year = dble(nmon_year) !! days per month [day/mon]
  real(wp), parameter :: sec_year = 31556926._wp
  real(wp), parameter :: sec_year_inv = 1._wp/sec_year
  real(wp), parameter :: sec_day  = sec_year / day_year   ! 8.765813d4, actual seconds per day
  real(wp), parameter :: sec_mon  = sec_year / mon_year
  real(dp), parameter :: sec_year_dp = 31556926._dp
  real(dp), parameter :: sec_day_dp  = sec_year_dp / day_year   ! 8.765813d4, actual seconds per day

  real(wp), parameter :: dt_day_atm  = 1._wp                ! atmosphere time step in days 
  real(wp), parameter :: dt_day_lnd  = 1._wp                ! land time step in days 
  real(wp), parameter :: dt_day_sic  = 1._wp                ! sea ice time step in days 
  real(wp)            :: dt_day_ice  != real(n_year_ice,wp)*real(nday_year,wp)   ! ice sheet model time step in days

  real(wp) :: dt_day_fastest, dt_day_ocn, dt_day_bgc, dt_day_smb, dt_day_bmb
  real(wp) :: dt_atm, dt_lnd, dt_ocn, dt_bgc, dt_sic, dt_smb, dt_bmb, dt_geo, dt_fastest
  integer :: n_year_ice, n_year_smb, n_year_geo
  integer :: n_accel

  integer :: nstep_day_atm, nstep_mon_atm, nstep_year_atm
  integer :: nstep_day_ocn, nstep_mon_ocn, nstep_year_ocn
  integer :: nstep_day_bgc, nstep_mon_bgc, nstep_year_bgc
  integer :: nstep_day_sic, nstep_mon_sic, nstep_year_sic
  integer :: nstep_day_lnd, nstep_mon_lnd, nstep_year_lnd
  integer :: nstep_day_smb, nstep_mon_smb, nstep_year_smb
  integer :: nstep_day_bmb, nstep_mon_bmb, nstep_year_bmb
  logical :: time_soy, time_eoy, time_som, time_eom
  logical :: time_soy_atm, time_eoy_atm, time_som_atm, time_eom_atm, time_sod_atm, time_eod_atm
  logical :: time_soy_ocn, time_eoy_ocn, time_som_ocn, time_eom_ocn, time_sod_ocn, time_eod_ocn
  logical :: time_soy_bgc, time_eoy_bgc, time_som_bgc, time_eom_bgc, time_sod_bgc, time_eod_bgc
  logical :: time_soy_sic, time_eoy_sic, time_som_sic, time_eom_sic, time_sod_sic, time_eod_sic
  logical :: time_soy_lnd, time_eoy_lnd, time_som_lnd, time_eom_lnd, time_sod_lnd, time_eod_lnd
  logical :: time_soy_smb, time_eoy_smb, time_som_smb, time_eom_smb, time_sod_smb, time_eod_smb
  logical :: time_soy_bmb, time_eoy_bmb, time_som_bmb, time_eom_bmb, time_sod_bmb, time_eod_bmb
  logical :: time_soy_bnd
  logical :: time_feedback_save, time_feedback_analysis
  logical :: time_spinup_cc_1, time_spinup_cc_2, time_call_daily_input_save, time_use_daily_input_save
  logical :: time_write_restart

  integer :: year_out_start
  integer :: y_out_ts, y_out_ts_geo, y_out_ts_clim, y_out_ts_smb, ny_out_ts, ny_out_ts_geo, &
  ny_out_ts_accel, ny_out_ts_smb
  integer :: nyout_cmn, nyout_atm, nyout_lnd, nyout_ocn, nyout_sic, nyout_bgc, &
  nyout_smb, nyout_bmb, nyout_ice, nyout_geo
  logical :: time_out_ts, time_out_ts_geo, time_out_ts_clim, time_out_ts_smb, time_out_cmn, &
  time_out_atm, time_out_lnd, time_out_ocn, time_out_sic, time_out_bgc, time_out_smb, &
  time_out_bmb, time_out_ice, time_out_geo

  logical :: time_call_atm, time_call_lnd, time_call_ocn, time_call_sic, time_call_bgc, &
  time_call_smb, time_call_bmb, time_call_ice, time_call_geo, time_call_clim

contains

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t i m e r _ i n i t
  !   Purpose    :  initialize climber model timer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine timer_init(filename)

    implicit none

    character(len=*), intent(in) :: filename

    call nml_read(filename,"control","nyears",nyears)
    call nml_read(filename,"control","year_ini",year_ini)
    call nml_read(filename,"control","dt_day_ocn",dt_day_ocn)
    call nml_read(filename,"control","dt_day_bgc",dt_day_bgc)
    call nml_read(filename,"control","dt_day_smb",dt_day_smb)
    call nml_read(filename,"control","dt_day_bmb",dt_day_bmb)
    call nml_read(filename,"control","n_year_ice",n_year_ice)
    call nml_read(filename,"control","n_year_smb",n_year_smb)
    call nml_read(filename,"control","n_year_geo",n_year_geo)
    call nml_read(filename,"control","year_out_start",year_out_start)
    call nml_read(filename,"control","nyout_cmn",nyout_cmn)
    call nml_read(filename,"control","nyout_atm",nyout_atm)
    call nml_read(filename,"control","nyout_lnd",nyout_lnd)
    call nml_read(filename,"control","nyout_ocn",nyout_ocn)
    call nml_read(filename,"control","nyout_bgc",nyout_bgc)
    call nml_read(filename,"control","nyout_sic",nyout_sic)
    call nml_read(filename,"control","nyout_smb",nyout_smb)
    call nml_read(filename,"control","nyout_bmb",nyout_bmb)
    call nml_read(filename,"control","nyout_ice",nyout_ice)
    call nml_read(filename,"control","nyout_geo",nyout_geo)
    call nml_read(filename,"control","n_accel",n_accel)

    ! make sure that SMB is not called more frequently than the acceleration factor
    n_year_smb = max(n_year_smb,n_accel)
    ! fixme: it should be possible that n_year_smb is any integer >= 1 
    if(mod(n_year_smb,n_accel).ne.0) then
      print *,'ERROR, timer parameter n_year_smb has to be an integer multiple of n_accel'
      print *,'n_year_smb',n_year_smb
      print *,'n_accel',n_accel
      stop
    end if

    ! ice sheet time step in days
    dt_day_ice  = real(n_year_ice,wp)*real(nday_year,wp)   ! ice sheet model time step in days

    ! find fastest timestep
    dt_day_fastest = 9999._wp
    if (flag_atm ) dt_day_fastest = minval([dt_day_atm, dt_day_fastest])
    if (flag_lnd ) dt_day_fastest = minval([dt_day_lnd, dt_day_fastest])
    if (flag_ocn ) dt_day_fastest = minval([dt_day_ocn, dt_day_fastest])
    if (flag_bgc ) dt_day_fastest = minval([dt_day_bgc, dt_day_fastest])
    if (flag_sic ) dt_day_fastest = minval([dt_day_sic, dt_day_fastest])
    if (flag_smb ) dt_day_fastest = minval([dt_day_smb, dt_day_fastest])
    if (flag_bmb ) dt_day_fastest = minval([dt_day_bmb, dt_day_fastest])
    if (flag_ice ) dt_day_fastest = minval([dt_day_ice, dt_day_fastest])
    if (flag_geo ) dt_day_fastest = minval([real(nday_year*n_year_geo,wp), dt_day_fastest])
    dt_day_fastest = minval([day_year,dt_day_fastest])
    dt_fastest = dt_day_fastest * sec_day

    print *,'dt_day_fastest',dt_day_fastest

    ! total number of time steps climate
    nyears_clim=ceiling(nyears/real(n_accel))*n_accel ! enlarge nyears to multiple of n_accel
    ! total number of time steps geo
    nyears_geo=ceiling(nyears/real(n_year_geo))*n_year_geo ! enlarge nyears to multiple of n_year_geo
    if (flag_geo) then
      ! total number of time steps required
      nyears=max(nyears_clim, nyears_geo)
    endif
    write(6,'(a,i8)') 'timer: nyears=', nyears
    nstep = nint(nyears*day_year/dt_day_fastest)


    ! set how often time series output is written to file
    ny_out_ts = min(100,nyears)

    ! number of time steps per year
    nstep_day = nint(1._wp/dt_day_fastest)
    nstep_mon = nint(day_mon/dt_day_fastest)
    nstep_year = nint(day_year/dt_day_fastest)
    ny_out_ts_accel = ny_out_ts*n_accel
    ny_out_ts_smb = ny_out_ts*n_year_smb
    ny_out_ts_geo = ny_out_ts*n_year_geo

    ! integer multiples of fastest time step
    step_atm = nint(dt_day_atm/dt_day_fastest)
    step_lnd = nint(dt_day_lnd/dt_day_fastest)
    step_ocn = nint(dt_day_ocn/dt_day_fastest)
    step_bgc = nint(dt_day_bgc/dt_day_fastest)
    step_sic = nint(dt_day_sic/dt_day_fastest)
    step_smb = nint(dt_day_smb/dt_day_fastest)
    step_bmb = nint(dt_day_bmb/dt_day_fastest)
    step_ice = nint(dt_day_ice/dt_day_fastest)

    nstep_day_atm = max(1, nint(1._wp/dt_day_atm))
    nstep_day_ocn = max(1, nint(1._wp/dt_day_ocn))
    nstep_day_bgc = max(1, nint(1._wp/dt_day_bgc))
    nstep_day_sic = max(1, nint(1._wp/dt_day_sic))
    nstep_day_lnd = max(1, nint(1._wp/dt_day_lnd))
    nstep_day_smb = max(1, nint(1._wp/dt_day_smb))
    nstep_day_bmb = max(1, nint(1._wp/dt_day_bmb))

    nstep_mon_atm = nint(nday_mon/dt_day_atm)
    nstep_mon_ocn = nint(nday_mon/dt_day_ocn)
    nstep_mon_bgc = nint(nday_mon/dt_day_bgc)
    nstep_mon_sic = nint(nday_mon/dt_day_sic)
    nstep_mon_lnd = nint(nday_mon/dt_day_lnd)
    nstep_mon_smb = nint(nday_mon/dt_day_smb)
    nstep_mon_bmb = nint(nday_mon/dt_day_bmb)

    nstep_year_atm = nint(nday_year/dt_day_atm)
    nstep_year_ocn = nint(nday_year/dt_day_ocn)
    nstep_year_bgc = nint(nday_year/dt_day_bgc)
    nstep_year_sic = nint(nday_year/dt_day_sic)
    nstep_year_lnd = nint(nday_year/dt_day_lnd)
    nstep_year_smb = nint(nday_year/dt_day_smb)
    nstep_year_bmb = nint(nday_year/dt_day_bmb)

    ! timesteps in seconds
    dt_atm = dt_day_atm * sec_day
    dt_lnd = dt_day_lnd * sec_day
    dt_ocn = dt_day_ocn * sec_day
    dt_bgc = dt_day_bgc * sec_day
    dt_sic = dt_day_sic * sec_day
    dt_smb = dt_day_smb * sec_day
    dt_bmb = dt_day_bmb * sec_day
    dt_geo = n_year_geo * sec_year

    ! initialize
    doy = 0
    mon = 0
    year = 0

!   write(6,'(a)') "---------------------------------------------"
!   write(6,'(a)') "From subroutine timer_init"
!   write(6,*) n_accel
!   write(6,*) step_atm, step_lnd, step_ocn, step_sic, step_bgc, step_smb, step_bmb, step_ice
!   write(6,'(a)') "---------------------------------------------"

   return

  end subroutine timer_init


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   Subroutine :  t i m e r _ u p d a t e
  !   Purpose    :  update climber model timer
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine timer_update

    implicit none
    real(dp) :: time_now
    logical :: year_call_accel
    logical :: year_call_smb

    time_now = step * dt_fastest   ! cumulated time in sec

    ! step of year
    soy = mod(step,nstep_year)
    if(soy .eq. 0) soy = nstep_year
    ! day of year
    doy = mod(floor((time_now-0.01_dp)/sec_day_dp)+1,nday_year)
    if(doy .eq. 0) doy = nday_year
    mon  = (doy-1)/nday_mon + 1
    year = floor((time_now-0.5_dp*dt_fastest)/sec_year_dp) + 1

    year_clim = ceiling(dble(year)/dble(n_accel))
    year_smb = ceiling(dble(year)/dble(n_year_smb))
    year_geo = ceiling(dble(year)/dble(n_year_geo))

    !year = nint(time_now/sec_year_dp)

    year_now = year_ini+year

    !print *,''
    !print *,'step, time_now',step, time_now
    !print *,'year,mon,doy,soy',year,mon,doy,soy
    !print *,(time_now-1._dp)/sec_day_dp,floor((time_now-1._dp)/sec_day_dp)+1,nint(time_now/sec_day_dp)
    !print *,(time_now-0.5_dp*dt_fastest)/sec_year_dp, year

!    if (doy.ne.soy) stop 'stop in subroutine timer_update'

    if (l_spinup_cc) then
      time_call_daily_input_save = year.gt.(year_start_offline-nyear_avg_offline) .and. year.le.year_start_offline
      time_use_daily_input_save = year.gt.year_start_offline
      time_spinup_cc_1 = soy.eq.1 .and. year.eq.(year_start_offline+1)  ! time after which to switch to ocn/bgc only setup
      time_spinup_cc_2 = soy.eq.1 .and. year.eq.nyears_spinup_bgc  ! time to end compensation of sediment fluxes and sediment spinup
    else
      time_call_daily_input_save = .false.
      time_use_daily_input_save = .false.
      time_spinup_cc_1 = .false. 
      time_spinup_cc_2 = .false. 
    endif

    year_call_accel = (mod(year,n_accel) .eq. 0).or.(year.eq.1)
    year_call_smb   = (mod(year,n_year_smb).eq.0).or.(year.eq.1)

    if (flag_atm .and. year_call_accel) then
      time_call_atm = (mod(step,step_atm) .eq. 1) .or. (step_atm.eq.1)
      time_soy_atm = soy .eq. 1
      time_sod_atm = (mod(soy,nstep_day_atm) .eq. 1) .or. (nstep_day_atm.eq.1)
      time_eod_atm = (mod(soy,nstep_day_atm) .eq. 0)
      time_eom_atm = (mod(doy,nday_mon) .eq. 0) .and. time_eod_atm
      time_eoy_atm = (mod(doy,nday_year) .eq. 0) .and. time_eod_atm
      time_out_atm = (mod(year,nyout_atm) .eq. 0) .and. year_now.ge.year_out_start
      time_feedback_save     = year.eq.nyears/2  ! save fields 
      time_feedback_analysis = year.eq.nyears    ! do feedback analysis 
    else
      time_call_atm =.false. 
      time_soy_atm = .false.
      time_sod_atm = .false.
      time_eod_atm = .false.
      time_eom_atm = .false.
      time_eoy_atm = .false.
      time_out_atm = .false.
      time_feedback_save = .false. 
      time_feedback_analysis = .false.
    endif

    if (flag_ocn .and. year_call_accel) then
      !time_call_ocn = (mod(step,step_ocn) .eq. 0)
      !time_soy_ocn = step_ocn .eq. soy
      time_call_ocn = (mod(step,step_ocn) .eq. 1) .or. (step_ocn.eq.1)
      time_soy_ocn = soy .eq. 1
      time_sod_ocn = (mod(soy,nstep_day_ocn) .eq. 1) .or. (nstep_day_ocn.eq.1)
      time_eod_ocn = (mod(soy,nstep_day_ocn) .eq. 0)
      time_eom_ocn = (mod(doy,nday_mon) .eq. 0) .and. time_eod_ocn
      time_eoy_ocn = (mod(doy,nday_year) .eq. 0) .and. time_eod_ocn
      time_out_ocn = (mod(year,nyout_ocn) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_ocn =.false. 
      time_soy_ocn = .false.
      time_sod_ocn = .false.
      time_eod_ocn = .false.
      time_eom_ocn = .false.
      time_eoy_ocn = .false.
      time_out_ocn = .false.
    endif

    if (flag_bgc .and. year_call_accel) then
      time_call_bgc = (mod(step,step_bgc) .eq. 1) .or. (step_bgc.eq.1)
      time_soy_bgc = soy .eq. 1
      time_sod_bgc = (mod(soy,nstep_day_bgc) .eq. 1) .or. (nstep_day_bgc.eq.1)
      time_eod_bgc = (mod(soy,nstep_day_bgc) .eq. 0)
      time_eom_bgc = (mod(doy,nday_mon) .eq. 0) .and. time_eod_bgc
      time_eoy_bgc = (mod(doy,nday_year) .eq. 0) .and. time_eod_bgc
      !time_eom_bgc = (mod(soy,(nstep_mon_bgc-1)*step_bgc+1) .eq. 0)
      !time_eoy_bgc = (mod(soy,(nstep_year_bgc-1)*step_bgc+1) .eq. 0)
      time_out_bgc = (mod(year,nyout_bgc) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_bgc =.false. 
      time_soy_bgc = .false.
      time_sod_bgc = .false.
      time_eod_bgc = .false.
      time_eom_bgc = .false.
      time_eoy_bgc = .false.
      time_out_bgc = .false.
    endif

    if (flag_lnd .and. year_call_accel) then
      time_call_lnd = (mod(step,step_lnd) .eq. 1) .or. (step_lnd.eq.1)
      time_soy_lnd = soy .eq. 1
      time_sod_lnd = (mod(soy,nstep_day_lnd) .eq. 1) .or. (nstep_day_lnd.eq.1)
      time_eod_lnd = (mod(soy,nstep_day_lnd) .eq. 0)
      time_eom_lnd = (mod(doy,nday_mon) .eq. 0) .and. time_eod_lnd
      time_eoy_lnd = (mod(doy,nday_year) .eq. 0) .and. time_eod_lnd
      time_out_lnd = (mod(year,nyout_lnd) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_lnd =.false. 
      time_soy_lnd = .false.
      time_sod_lnd = .false.
      time_eod_lnd = .false.
      time_eom_lnd = .false.
      time_eoy_lnd = .false.
      time_out_lnd = .false.
    endif

    if (flag_sic .and. year_call_accel) then
      time_call_sic = (mod(step,step_sic) .eq. 1) .or. (step_sic.eq.1)
      time_soy_sic = soy .eq. 1
      time_sod_sic = (mod(soy,nstep_day_sic) .eq. 1) .or. (nstep_day_sic.eq.1)
      time_eod_sic = (mod(soy,nstep_day_sic) .eq. 0)
      time_eom_sic = (mod(doy,nday_mon) .eq. 0) .and. time_eod_sic
      time_eoy_sic = (mod(doy,nday_year) .eq. 0) .and. time_eod_sic
      time_out_sic = (mod(year,nyout_sic) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_sic =.false. 
      time_soy_sic = .false.
      time_sod_sic = .false.
      time_eod_sic = .false.
      time_eom_sic = .false.
      time_eoy_sic = .false.
      time_out_sic = .false.
    endif

    if (flag_smb .and. year_call_smb) then 
      time_call_smb = (mod(step,step_smb) .eq. 1) .or. (step_smb.eq.1)
      time_soy_smb = soy .eq. 1
      time_sod_smb = (mod(soy,nstep_day_smb) .eq. 1) .or. (nstep_day_smb.eq.1)
      time_eod_smb = (mod(soy,nstep_day_smb) .eq. 0)
      time_eom_smb = (mod(doy,nday_mon) .eq. 0) .and. time_eod_smb
      time_eoy_smb = (mod(doy,nday_year) .eq. 0) .and. time_eod_smb
      time_out_smb = (mod(year,nyout_smb) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_smb =.false. 
      time_soy_smb = .false.
      time_sod_smb = .false.
      time_eod_smb = .false.
      time_eom_smb = .false.
      time_eoy_smb = .false.
      time_out_smb = .false.
    endif
    
    if (flag_bmb .and. year_call_accel) then 
      time_call_bmb = (mod(step,step_bmb) .eq. 1) .or. (step_bmb.eq.1)
      time_soy_bmb = soy .eq. 1
      time_sod_bmb = (mod(soy,nstep_day_bmb) .eq. 1) .or. (nstep_day_bmb.eq.1)
      time_eod_bmb = (mod(soy,nstep_day_bmb) .eq. 0)
      time_eom_bmb = (mod(soy,(nstep_mon_bmb-1)*step_bmb+1) .eq. 0) .or. (nstep_mon_bmb.eq.1)
      time_eoy_bmb = (mod(soy,(nstep_year_bmb-1)*step_bmb+1) .eq. 0) 
      time_out_bmb = (mod(year,nyout_bmb) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_bmb =.false. 
      time_soy_bmb = .false.
      time_eom_bmb = .false.
      time_eoy_bmb = .false.
      time_out_bmb = .false.
    endif
    
    if (flag_ice) then
      time_call_ice = (mod(step+nstep_day-1,step_ice) .eq. 0)
      time_out_ice = (mod(year,nyout_ice) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_ice =.false. 
      time_out_ice = .false.
    endif

    if (flag_geo .or. ifake_geo.eq.1) then
      time_call_geo = soy.eq.1 .and. mod(year,n_year_geo).eq.0
      time_out_geo = (mod(year,nyout_geo) .eq. 0) .and. year_now.ge.year_out_start
    else
      time_call_geo =.false. 
      time_out_geo = .false.
    endif

    ! time flags for start of year
    time_soy = soy .eq. 1
    time_soy_bnd = soy .eq. 1

    ! time flags seasonal cycle
    if (year_call_accel) then
      time_call_clim = (mod(step, 1) .eq. 0)
      ! time flag for end of month
      if (nstep_mon.eq.0) then
        time_eom = .false.
      else
        time_eom     = (mod(soy,nstep_mon) .eq. 0)
      endif
      ! time flags for end of year
      time_eoy = (mod(soy,nstep_year) .eq. 0)
    else
      time_call_clim = .false.
      time_eom = .false.
      time_eoy = .false.
    end if

    ! time flags for writing output
    time_out_cmn = (mod(year,nyout_cmn) .eq. 0) .and. year_now.ge.year_out_start

    ! control for time series output for geo
    if (time_soy) then
      y_out_ts_geo = mod(year,ny_out_ts_geo)
      if (y_out_ts_geo==0) y_out_ts_geo = ny_out_ts_geo
      if (year==nyears_geo .or. y_out_ts_geo==ny_out_ts_geo) then
        time_out_ts_geo = .true.
      else
        time_out_ts_geo = .false.
      endif
      y_out_ts_geo = y_out_ts_geo/n_year_geo
    endif

    ! control for time series output for climate (except ice sheets and smb)
    if (time_soy) then
      y_out_ts_clim = mod(year,ny_out_ts_accel)
      if (y_out_ts_clim==0) y_out_ts_clim = ny_out_ts_accel
      if (year==nyears_clim .or. y_out_ts_clim==ny_out_ts_accel) then
        time_out_ts_clim = .true.
      else
        time_out_ts_clim = .false.
      endif
      y_out_ts_clim = y_out_ts_clim/n_accel
    endif

    ! control for time series output for smb
    if (time_soy) then
      y_out_ts_smb = mod(year,ny_out_ts_smb)
      if (y_out_ts_smb==0) y_out_ts_smb = ny_out_ts_smb
      if (year==nyears .or. y_out_ts_smb==ny_out_ts_smb) then
        time_out_ts_smb = .true.
      else
        time_out_ts_smb = .false.
      endif
      y_out_ts_smb = y_out_ts_smb/n_year_smb
      y_out_ts_smb = max(1,y_out_ts_smb)
    endif

    if(n_year_geo.gt.1.and.year.eq.1) then
      y_out_ts_geo = 1
    end if

   if(n_accel.gt.1.and.year.eq.1) then
      y_out_ts_clim = 1
      y_out_ts_smb = 1
    end if
    
    ! control for time series output for ice sheets
    if (time_soy) then
      y_out_ts = mod(year,ny_out_ts)
      if (y_out_ts==0) y_out_ts = ny_out_ts
      if (year==nyears .or. y_out_ts==ny_out_ts) then
        time_out_ts = .true.
      else
        time_out_ts = .false.
      endif
    endif

    ! control of writing of restart files
    time_write_restart = .false.
    if (i_write_restart.eq.0) then
      ! end of simulation only
      time_write_restart = (mod(soy,nstep_year).eq.0) .and. year.eq.nyears
    else if (i_write_restart.eq.1) then
      ! regular with frequency n_year_write_restart
      time_write_restart = (mod(soy,nstep_year).eq.0) .and. (mod(year,n_year_write_restart).eq.0 .or. year.eq.nyears)
    else if (i_write_restart.eq.2) then
      ! at specified times
      time_write_restart = (mod(soy,nstep_year).eq.0) .and. any(years_write_restart==year_now)
    endif

!    print *
!    write(6,'(4(a,i6))') "step = ", step, ", year = ", year, ", year_clim = ", year_clim, &
!    ", mon = ",mon
!    write(6,'(a,l1)') 'year_call_accel=', year_call_accel
!    write(6,'(a,i4)') "doy  = ", doy
!    write(6,'(a,i4,a,l1)') "soy  = ", soy, ', time_soy=', time_soy
!    write(6,'(2(a,l1),a,i8)') 'time_out_ts_clim=', time_out_ts_clim, &
!    ', time_out_ts=', time_out_ts, ', y_out_ts=', y_out_ts
!    write(6,'(5(a,l1))') 'time_call_atm=', time_call_atm, ', time_soy_atm=', time_soy_atm, &
!    ', time_eom_atm=', time_eom_atm, ', time_eoy_atm=', time_eoy_atm, &
!    ', time_out_atm=', time_out_atm
!    write(6,'(2(a,l1))') 'time_call_smb=', time_call_smb, ', time_soy_smb=', time_soy_smb
!    write(6,'(3(a,l1))') 'time_out_cmn=', time_out_cmn, ', time_eom=', time_eom, ', time_eoy=', time_eoy
!    write(6,'(2(a,l1))') 'time_call_ice=', time_call_ice, &
!   ', time_out_ice=', time_out_ice

   return

  end subroutine timer_update

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  m o n t h l y 2 d a i l y
  ! Purpose    :  Return indices of months and weights
  !               to be used for linear daily interpolation
  !               from monthly values
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine monthly2daily(m0,m1,wtm0,wtm1)

    implicit none

    integer, dimension(nday_year), intent(out) :: m0, m1
    real(wp), dimension(nday_year), intent(out) :: wtm0, wtm1

    integer :: k, m, n, mid, temp
    real(wp) :: wttot

    ! Get length of arrays to fill, midpoint of month (in days)
    ! and the total weight (total number of days in a month)
    n = size(m0)

    mid = nday_mon / 2

    ! ####################################################################
    ! Interpolate data in time: monthly => daily
    ! ####################################################################   
    do k = 1, n

      do m = 1, nmon_year+1
        if ( m*nday_mon-mid .gt. k ) exit
      end do
      m1(k) = m; m0(k) = m1(k)-1

!      if ( m1(k) .gt. 12 ) m1(k) =  1
!      if ( m0(k) .lt.  1 ) m0(k) = 12

      temp = mod(k-mid,nday_mon)
      if( temp .lt. 0 ) temp = nday_mon+k-mid
      wtm1(k) = real(abs(mod(temp,nday_mon)),wp)
      wtm0(k) = day_mon - wtm1(k)

      wttot = wtm0(k) + wtm1(k)
      wtm0(k) = wtm0(k)/wttot; wtm1(k) = wtm1(k)/wttot

    end do

    return

  end subroutine monthly2daily


  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  t i m e r _ p r i n t 
  ! Purpose    :  print timer information 
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine timer_print

    implicit none

    if (time_call_atm) print *,'atm: sod=',time_sod_atm,' eod=',time_eod_atm,' eom=',time_eom_atm,' eoy=',time_eoy_atm
    if (time_call_ocn) print *,'ocn: sod=',time_sod_ocn,' eod=',time_eod_ocn,' eom=',time_eom_ocn,' eoy=',time_eoy_ocn
    if (time_call_sic) print *,'sic: sod=',time_sod_sic,' eod=',time_eod_sic,' eom=',time_eom_sic,' eoy=',time_eoy_sic
    if (time_call_lnd) print *,'lnd: sod=',time_sod_lnd,' eod=',time_eod_lnd,' eom=',time_eom_lnd,' eoy=',time_eoy_lnd
    if (time_call_bgc) print *,'bgc: sod=',time_sod_bgc,' eod=',time_eod_bgc,' eom=',time_eom_bgc,' eoy=',time_eoy_bgc
    if (time_call_smb) print *,'smb: sod=',time_sod_smb,' eod=',time_eod_smb,' eom=',time_eom_smb,' eoy=',time_eoy_smb
    if (time_call_bmb) print *,'bmb: sod=',time_sod_bmb,' eod=',time_eod_bmb,' eom=',time_eom_bmb,' eoy=',time_eoy_bmb

    return

  end subroutine timer_print

end module timer


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module : i c e _ s y n _ t o p o
!
!  Purpose : Synthetic ice-sheet topography kernels: signed distance to
!            a target-extent mask and linear / perfect-plastic surface
!            elevation profiles built from it. Pure routines with no
!            model state; shared by the synthetic ice model (ice_syn)
!            and the simple SMB scheme (smb_simple).
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2017-2024 Potsdam Institute for Climate Impact Research,
!                         Matteo Willeit and Reinhard Calov
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
module ice_syn_topo

    use precision, only : wp
    implicit none
    private

    public :: compute_signed_distance
    public :: compute_z_syn_linear
    public :: compute_z_syn_plastic

contains

    !=================================================================
    ! Signed distance to mask boundary (m)
    !=================================================================

    !-----------------------------------------------------------------
    !> Per-cell signed distance to the nearest target-mask boundary, in
    !> metres. Positive inside the mask, negative outside.
    !>
    !> x(:,:), y(:,:) are the (unique) 2-D coordinates of each cell.
    !> `units` selects how they are interpreted (optional, default "m"):
    !>   "m"       Cartesian, already metres: sqrt(dx^2 + dy^2)
    !>   "km"      Cartesian in km; result converted to metres (x1000)
    !>   "degrees" x = lon, y = lat (deg); local flat-Earth metric
    !>             (R = 6.371e6 m, mid-latitude cos(phi) factor for dlon,
    !>             longitude wrap to [-180, 180)).
    !> The result is always in metres, regardless of `units`.
    !-----------------------------------------------------------------
    subroutine compute_signed_distance(d, mask_target, x, y, units)
        real(wp), intent(out) :: d(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(wp), intent(in)  :: x(:,:)
        real(wp), intent(in)  :: y(:,:)
        character(len=*), intent(in), optional :: units

        integer  :: nx, ny, i, j, ib, jb, k, nb, kk
        logical  :: m, is_bnd, latlon
        real(wp) :: xq, yq, xb, yb, dx, dy, dphi, dlam
        real(wp) :: phi_mid, dy_m, dx_m, dist, dmin, scale
        character(len=16) :: units_use
        integer, allocatable :: bi(:), bj(:)
        logical, allocatable :: bm(:)
        real(wp), parameter :: R_earth = 6.371e6_wp   ! m
        real(wp), parameter :: deg2rad = 3.141592653589793_wp / 180.0_wp

        nx = size(mask_target, 1)
        ny = size(mask_target, 2)

        units_use = "m"
        if (present(units)) units_use = units

        if (size(d, 1) /= nx .or. size(d, 2) /= ny) then
            error stop "compute_signed_distance: d shape /= mask_target"
        end if
        if (size(x, 1) /= nx .or. size(x, 2) /= ny .or. &
            size(y, 1) /= nx .or. size(y, 2) /= ny) then
            error stop "compute_signed_distance: x, y shape /= mask_target"
        end if
        if (nx < 2 .or. ny < 2) then
            error stop "compute_signed_distance: nx, ny must be >= 2"
        end if

        select case (trim(units_use))
        case ("degrees")
            latlon = .true.;  scale = 1.0_wp
        case ("km")
            latlon = .false.; scale = 1000.0_wp
        case ("m")
            latlon = .false.; scale = 1.0_wp
        case default
            error stop "compute_signed_distance: units must be 'm', 'km', or 'degrees'"
        end select

        ! Pass 1: count boundary cells (4-connectivity).
        nb = 0
        do j = 1, ny
            do i = 1, nx
                m      = mask_target(i, j)
                is_bnd = .false.
                if (i > 1) then
                    if (mask_target(i - 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. i < nx) then
                    if (mask_target(i + 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j > 1) then
                    if (mask_target(i, j - 1) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j < ny) then
                    if (mask_target(i, j + 1) .neqv. m) is_bnd = .true.
                end if
                if (is_bnd) nb = nb + 1
            end do
        end do

        if (nb == 0) then
            if (mask_target(1, 1)) then
                d = huge(1.0_wp)
            else
                d = -huge(1.0_wp)
            end if
            return
        end if

        allocate(bi(nb), bj(nb), bm(nb))

        ! Pass 2: fill boundary lists.
        k = 0
        do j = 1, ny
            do i = 1, nx
                m      = mask_target(i, j)
                is_bnd = .false.
                if (i > 1) then
                    if (mask_target(i - 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. i < nx) then
                    if (mask_target(i + 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j > 1) then
                    if (mask_target(i, j - 1) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j < ny) then
                    if (mask_target(i, j + 1) .neqv. m) is_bnd = .true.
                end if
                if (is_bnd) then
                    k = k + 1
                    bi(k) = i
                    bj(k) = j
                    bm(k) = m
                end if
            end do
        end do

        ! Pass 3: min distance to opposite-mask boundary, per query cell.
        do j = 1, ny
            do i = 1, nx
                m     = mask_target(i, j)
                xq    = x(i, j)
                yq    = y(i, j)
                dmin  = huge(1.0_wp)
                do kk = 1, nb
                    if (bm(kk) .eqv. m) cycle
                    ib    = bi(kk)
                    jb    = bj(kk)
                    xb    = x(ib, jb)
                    yb    = y(ib, jb)
                    if (latlon) then
                        dphi    = yq - yb
                        dlam    = xq - xb
                        dlam    = modulo(dlam + 180.0_wp, 360.0_wp) - 180.0_wp
                        phi_mid = 0.5_wp * (yq + yb)
                        dy_m    = R_earth * dphi * deg2rad
                        dx_m    = R_earth * cos(phi_mid * deg2rad) * dlam * deg2rad
                        dist    = sqrt(dy_m * dy_m + dx_m * dx_m)
                    else
                        dx   = xq - xb
                        dy   = yq - yb
                        dist = scale * sqrt(dx * dx + dy * dy)
                    end if
                    if (dist < dmin) dmin = dist
                end do
                if (m) then
                    d(i, j) =  dmin
                else
                    d(i, j) = -dmin
                end if
            end do
        end do

        deallocate(bi, bj, bm)
    end subroutine compute_signed_distance

    !=================================================================
    ! Synthetic-elevation profiles (linear and plastic)
    !=================================================================

    !-----------------------------------------------------------------
    !> Linear wedge profile both inside and outside the mask. d_m is the
    !> signed distance in metres and slope is in m/m.
    !>     z_syn_raw  = clamp(slope * d_m, -z_max_out, +z_max_in)
    !>     z_syn(i,j) = max(z_syn_raw, z_sur(i,j))  if  mask_target(i,j)
    !>                = min(z_syn_raw, z_sur(i,j))  otherwise
    !-----------------------------------------------------------------
    subroutine compute_z_syn_linear(z_syn, d_m, z_sur, mask_target, &
                                    slope, z_max_in, z_max_out)
        real(wp), intent(out) :: z_syn(:,:)
        real(wp), intent(in)  :: d_m(:,:)
        real(wp), intent(in)  :: z_sur(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(wp), intent(in)  :: slope
        real(wp), intent(in)  :: z_max_in
        real(wp), intent(in)  :: z_max_out

        integer  :: nx, ny, i, j
        real(wp) :: zr

        nx = size(z_syn, 1)
        ny = size(z_syn, 2)

        if (size(d_m, 1)         /= nx .or. size(d_m, 2)         /= ny .or. &
            size(z_sur, 1)       /= nx .or. size(z_sur, 2)       /= ny .or. &
            size(mask_target, 1) /= nx .or. size(mask_target, 2) /= ny) then
            error stop "compute_z_syn_linear: shape mismatch among inputs"
        end if
        if (z_max_in < 0.0_wp .or. z_max_out < 0.0_wp) then
            error stop "compute_z_syn_linear: z_max_in and z_max_out must be >= 0"
        end if

        do j = 1, ny
            do i = 1, nx
                zr = slope * d_m(i, j)
                if (zr >  z_max_in)  zr =  z_max_in
                if (zr < -z_max_out) zr = -z_max_out
                if (mask_target(i, j)) then
                    z_syn(i, j) = max(zr, z_sur(i, j))
                else
                    z_syn(i, j) = min(zr, z_sur(i, j))
                end if
            end do
        end do
    end subroutine compute_z_syn_linear

    !-----------------------------------------------------------------
    !> Perfect-plasticity (Nye/Vialov) profile inside; linear outside.
    !> d_m is the signed distance in metres and slope_out is in m/m.
    !>   inside  (d >= 0):  z_syn_raw = min(z_max_in,  C * sqrt(d_m))
    !>                       C = sqrt(2 * tau0 / (rho_ice*g))
    !>   outside (d <  0):  z_syn_raw = max(-z_max_out, slope_out*d_m)
    !-----------------------------------------------------------------
    subroutine compute_z_syn_plastic(z_syn, d_m, z_sur, mask_target,         &
                                     tau0, slope_out, z_max_in, z_max_out,   &
                                     rho_ice, g)
        real(wp), intent(out) :: z_syn(:,:)
        real(wp), intent(in)  :: d_m(:,:)
        real(wp), intent(in)  :: z_sur(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(wp), intent(in)  :: tau0
        real(wp), intent(in)  :: slope_out
        real(wp), intent(in)  :: z_max_in
        real(wp), intent(in)  :: z_max_out
        real(wp), optional, intent(in) :: rho_ice
        real(wp), optional, intent(in) :: g

        integer  :: nx, ny, i, j
        real(wp) :: zr, d, C, rho_use, g_use

        nx = size(z_syn, 1)
        ny = size(z_syn, 2)

        if (size(d_m, 1)         /= nx .or. size(d_m, 2)         /= ny .or. &
            size(z_sur, 1)       /= nx .or. size(z_sur, 2)       /= ny .or. &
            size(mask_target, 1) /= nx .or. size(mask_target, 2) /= ny) then
            error stop "compute_z_syn_plastic: shape mismatch among inputs"
        end if
        if (z_max_in < 0.0_wp .or. z_max_out < 0.0_wp) then
            error stop "compute_z_syn_plastic: z_max_in and z_max_out must be >= 0"
        end if
        if (tau0 <= 0.0_wp) then
            error stop "compute_z_syn_plastic: tau0 must be > 0"
        end if

        rho_use = 910.0_wp
        g_use   =   9.81_wp
        if (present(rho_ice)) rho_use = rho_ice
        if (present(g))       g_use   = g

        C = sqrt(2.0_wp * tau0 / (rho_use * g_use))

        do j = 1, ny
            do i = 1, nx
                d = d_m(i, j)
                if (d >= 0.0_wp) then
                    zr = C * sqrt(d)
                    if (zr > z_max_in) zr = z_max_in
                else
                    zr = slope_out * d
                    if (zr < -z_max_out) zr = -z_max_out
                end if
                if (mask_target(i, j)) then
                    z_syn(i, j) = max(zr, z_sur(i, j))
                else
                    z_syn(i, j) = min(zr, z_sur(i, j))
                end if
            end do
        end do
    end subroutine compute_z_syn_plastic

end module ice_syn_topo

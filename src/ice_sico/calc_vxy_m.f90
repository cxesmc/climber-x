!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Module :  c a l c _ v x y _ m
!
!> @file
!!
!! Computation of the horizontal velocity vx, vy.
!!
!! @section Copyright
!!
!! Copyright 2009-2017 Ralf Greve, Tatsuru Sato, Thomas Goelles, Jorge Bernales
!!
!! @section License
!!
!! This file is part of SICOPOLIS.
!!
!! SICOPOLIS is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! SICOPOLIS is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with SICOPOLIS.  If not, see <http://www.gnu.org/licenses/>.
!<
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy.
!<------------------------------------------------------------------------------
module calc_vxy_m

  use timer, only : sec_year, year

  use sico_types_m
  use sico_state
  use sico_grid_mod
  use sico_timer
  use sico_params, only : sico_par_class, RHO_SW, eps_wp
  use ice_material_properties_m, only : ratefac_c, ratefac_t, ratefac_c_t, creep, viscosity

  !$  use omp_lib

  implicit none

  private
  public :: calc_dzs_dxy_aux, calc_vxy_b_sia, calc_vxy_sia, calc_vxy_static, calc_vxy_ssa

contains

!-------------------------------------------------------------------------------
!> Computation of the auxiliary surface gradients dzs_dx_aux, dzs_dy_aux
!! (optional one-sided gradients at the grounding line).
!<------------------------------------------------------------------------------
subroutine calc_dzs_dxy_aux(st,grd)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd

integer :: i, j
real(wp)     :: inv_dx, inv_dy
real(wp)     :: rhosw_rho_ratio
real(wp)     :: H_mid, zl_mid, zs_mid

inv_dx          = 1.0_wp/grd%dxi
inv_dy          = 1.0_wp/grd%deta
rhosw_rho_ratio = RHO_SW/RHO

st%dzs_dx_aux = st%dzs_dxi
st%dzs_dy_aux = st%dzs_deta

do i=0, grd%IMAX-1
do j=1, grd%JMAX-1
   ! inner point on the staggered grid in x-direction

   if ( (st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j,i+1)) &
        .or. &
        (st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j,i+1)) &
      ) then
        ! one neighbour is floating ice and the other is grounded ice
        ! (grounding line)

      H_mid  = 0.5_wp*((st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j,i+1)+st%H_t(j,i+1)))
      zl_mid = 0.5_wp*(st%zl(j,i)+st%zl(j,i+1))
      zs_mid = 0.5_wp*(st%zs(j,i)+st%zs(j,i+1))

      if (H_mid < (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) then
         ! floating condition is satisfied

         if ( &
              (st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j,i+1)) &
              .and. &
              (i+2 <= grd%IMAX) &
            ) then

            if ((st%maske(j,i+2) == 3).or.(st%maske(j,i+2) == 2)) &
               st%dzs_dx_aux(j,i) = (0.5_wp*(st%zs(j,i+1)+st%zs(j,i+2))-zs_mid)*inv_dx
                                 ! one-sided gradient into floating ice

         else if ( &
              (st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j,i+1)) &
              .and. &
              (i-1 >= 0) &
            ) then

            if ((st%maske(j,i-1) == 3).or.(st%maske(j,i-1) == 2)) &
               st%dzs_dx_aux(j,i) = (zs_mid-0.5_wp*(st%zs(j,i)+st%zs(j,i-1)))*inv_dx
                                 ! one-sided gradient into floating ice

         end if

      else   ! H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio,
             ! floating condition is not satisfied

         if ( &
              (st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j,i+1)) &
              .and. &
              (i-1 >= 0) &
            ) then

            if ((st%maske(j,i-1) == 0).or.(st%maske(j,i-1) == 1)) &
               st%dzs_dx_aux(j,i) = (zs_mid-0.5_wp*(st%zs(j,i)+st%zs(j,i-1)))*inv_dx
                                 ! one-sided gradient into grounded ice

         else if ( &
              (st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j,i+1)) &
              .and. &
              (i+2 <= grd%IMAX) &
            ) then

            if ((st%maske(j,i+2) == 0).or.(st%maske(j,i+2) == 1)) &
               st%dzs_dx_aux(j,i) = (0.5_wp*(st%zs(j,i+1)+st%zs(j,i+2))-zs_mid)*inv_dx
                                 ! one-sided gradient into grounded ice

         end if

      end if


   end if

end do
end do

do i=1, grd%IMAX-1
do j=0, grd%JMAX-1
   ! inner point on the staggered grid in y-direction

   if ( (st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j+1,i)) &
        .or. &
        (st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j+1,i)) &
      ) then
        ! one neighbour is floating ice and the other is grounded ice
        ! (grounding line)

      H_mid  = 0.5_wp*((st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j+1,i)+st%H_t(j+1,i)))
      zl_mid = 0.5_wp*(st%zl(j,i)+st%zl(j+1,i))
      zs_mid = 0.5_wp*(st%zs(j,i)+st%zs(j+1,i))

      if (H_mid < (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) then
         ! floating condition is satisfied

         if ( &
              (st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j+1,i)) &
              .and. &
              (j+2 <= grd%JMAX) &
            ) then

            if ((st%maske(j+2,i) == 3).or.(st%maske(j+2,i) == 2)) &
               st%dzs_dy_aux(j,i) = (0.5_wp*(st%zs(j+1,i)+st%zs(j+2,i))-zs_mid)*inv_dy
                                 ! one-sided gradient into floating ice

         else if ( &
              (st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j+1,i)) &
              .and. &
              (j-1 >= 0) &
            ) then

            if ((st%maske(j-1,i) == 3).or.(st%maske(j-1,i) == 2)) &
               st%dzs_dy_aux(j,i) = (zs_mid-0.5_wp*(st%zs(j,i)+st%zs(j-1,i)))*inv_dy
                                 ! one-sided gradient into floating ice

         end if

      else   ! H_mid >= (z_sl-zl_mid)*rhosw_rho_ratio,
             ! floating condition is not satisfied

         if ( &
              (st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j+1,i)) &
              .and. &
              (j-1 >= 0) &
            ) then

            if ((st%maske(j-1,i) == 0).or.(st%maske(j-1,i) == 1)) &
               st%dzs_dy_aux(j,i) = (zs_mid-0.5_wp*(st%zs(j,i)+st%zs(j-1,i)))*inv_dy
                                 ! one-sided gradient into grounded ice

         else if ( &
              (st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j+1,i)) &
              .and. &
              (j+2 <= grd%JMAX) &
            ) then

            if ((st%maske(j+2,i) == 0).or.(st%maske(j+2,i) == 1)) &
               st%dzs_dy_aux(j,i) = (0.5_wp*(st%zs(j+1,i)+st%zs(j+2,i))-zs_mid)*inv_dy
                                 ! one-sided gradient into grounded ice

         end if

      end if

   end if

end do
end do

end subroutine calc_dzs_dxy_aux

!-------------------------------------------------------------------------------
!> Computation of the basal horizontal velocity vx_b, vy_b in the shallow ice
!! approximation.
!<------------------------------------------------------------------------------
subroutine calc_vxy_b_sia(st,grd,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(in) :: par

integer :: i, j
integer :: r_smw_1, s_smw_1, r_smw_2, s_smw_2
integer :: mask_sed
real(wp) :: gamma_slide_inv
real(wp) :: smw_coeff_1, smw_coeff_2
real(wp) :: cvxy1, cvxy1a, cvxy1b, ctau1, ctau1a, ctau1b
real(wp) :: f_pmp
real(wp) :: temp_diff
real(wp) :: Hw0_slide_inv, ratio_Hw_slide
real(wp) :: ramp_up_factor
real(wp) :: vh_max, vh_max_inv
real(wp) :: year_sec_inv
real(wp) :: w_sed

year_sec_inv  = 1.0_wp/sec_year

gamma_slide_inv = 1._wp/max(par%gamma_slide, eps)

!-------- Sliding-law coefficients --------

! different sliding laws for hard rock and soft sediment

do i=0, grd%IMAX
do j=0, grd%JMAX

  if (par%i_slide_rock_sed==1) then
    ! separate rock and sediment sliding using a critical sediment thickness threshold

   if (st%h_sed(j,i)>=par%h_sed_thresh) then
     mask_sed = 1
     st%f_sed(j,i) = 1._wp
   else
     mask_sed = 0
     st%f_sed(j,i) = 0._wp
   endif
   st%sub_melt_flag(j,i)   = (par%GAMMA_SLIDE >= eps)
   st%gamma_slide_inv(j,i) = gamma_slide_inv
   if (mask_sed==0) then   ! hard bedrock 
      st%c_slide(j,i)         = par%C_SLIDE*year_sec_inv
      st%p_weert(j,i)         = par%P_WEERT
      st%q_weert(j,i)         = par%Q_WEERT
      st%c_fric(j,i)          = par%c_fric_rock
      st%vb_t(j,i)            = par%vb_t_rock
      st%delta(j,i)           = par%delta_rock
    else if (mask_sed==1) then ! soft sediment
      st%c_slide(j,i)         = par%C_SLIDE_SEDI*year_sec_inv
      st%p_weert(j,i)         = par%P_WEERT_SEDI
      st%q_weert(j,i)         = par%Q_WEERT_SEDI
      st%c_fric(j,i)          = par%c_fric_sed
      st%vb_t(j,i)            = par%vb_t_sed
      st%delta(j,i)           = par%delta_sed
   end if

  else if (par%i_slide_rock_sed==2) then
    ! continuous transition between rock and sediment sliding
    ! NOTE: this assumes sediment sliding law to be the SAME as rock 

    w_sed = (st%h_sed(j,i)-par%h_sed_min)/(par%h_sed_max-par%h_sed_min)
    w_sed = max(0._wp,w_sed)
    w_sed = min(1._wp,w_sed)
    st%f_sed(j,i) = w_sed
    st%c_slide(j,i)         = ((1._wp-w_sed)*par%C_SLIDE+w_sed*par%C_SLIDE_SEDI)*year_sec_inv
    st%p_weert(j,i)         = par%P_WEERT
    st%q_weert(j,i)         = par%Q_WEERT
    st%c_fric(j,i)          = (1._wp-w_sed)*par%c_fric_rock+w_sed*par%c_fric_sed
    st%vb_t(j,i)            = (1._wp-w_sed)*par%vb_t_rock+w_sed*par%vb_t_sed
    st%delta(j,i)           = (1._wp-w_sed)*par%delta_rock+w_sed*par%delta_sed
    st%gamma_slide_inv(j,i) = gamma_slide_inv
    st%sub_melt_flag(j,i)   = (par%GAMMA_SLIDE >= eps)

  endif

  ! special treatment for marine sediments, enhance sliding
  if (st%zl0(j,i).lt.0._wp) then
    ! assume marine sediments below present-day sea level
    st%c_slide(j,i) = par%enh_slide_marine_sed*st%c_slide(j,i)
    st%c_fric(j,i)  = 1._wp/par%enh_slide_marine_sed*st%c_fric(j,i)
  endif

end do
end do

do i=0, grd%IMAX
do j=0, grd%JMAX
   st%p_weert_inv(j,i) = 1.0_wp/max(real(st%p_weert(j,i),wp), eps)
end do
end do

!  ------ Ramping up basal sliding

if (par%nyears_ramp_up_slide>0) then
  if (year<=par%nyears_ramp_up_slide) then

    ramp_up_factor = real(year,wp)/real(par%nyears_ramp_up_slide,wp)
    ramp_up_factor = max(min(ramp_up_factor, 1.0_wp), 0.0_wp) ! constrain to interval [0,1]
    ! make transition smooth (quintic function)
    ramp_up_factor = ramp_up_factor*ramp_up_factor*ramp_up_factor*(10.0_wp + ramp_up_factor*(-15.0_wp+6.0_wp*ramp_up_factor))

    ! scale sliding coefficient
    st%c_slide = st%c_slide * ramp_up_factor

    ! scale friction coefficient
    st%c_fric = st%c_fric * 1._wp/max(eps_wp,ramp_up_factor)

  end if
endif


!-------- Computation of basal stresses --------

!  ------ Basal pressure p_b, basal water pressure p_b_w,
!         reduced pressure p_b_red

do i=0, grd%IMAX
do j=0, grd%JMAX

   if ((st%maske(j,i) == 0).or.st%flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

      ! overburden pressure
      st%p_b(j,i)         = max(RHO*G*(st%H_c(j,i)+st%H_t(j,i)), 0.0_wp)

      ! basal water pressure and reduced basal pressure
      if (par%i_p_b_red.eq.1) then

        st%p_b_w(j,i)       = RHO_SW*G*max((st%z_sl(j,i)-st%zb(j,i)), 0.0_wp)
        st%p_b_red(j,i)     = max(st%p_b(j,i)-st%p_b_w(j,i), 0.0_wp)
        ! limit in order to avoid very small values, which may lead to huge sliding velocities in the SIA
        st%p_b_red_lim(j,i) = max(st%p_b_red(j,i), par%red_pres_limit_fact*st%p_b(j,i))   

      else if (par%i_p_b_red.eq.2) then
        ! reduced basal pressure = delta * overburden pressure (where temperate base)

        if (st%n_cts(j,i) == -1) then   ! cold ice base
          if (st%sub_melt_flag(j,i)) then
            temp_diff = max((st%temp_c_m(0,j,i)-st%temp_c(0,j,i)), 0.0_wp)
            f_pmp = exp(-st%gamma_slide_inv(j,i)*temp_diff)  ! 'virtual' fraction at pressure melting point
          else
            f_pmp = 0._wp
          end if
          st%p_b_red(j,i) = (1._wp-f_pmp)*st%p_b(j,i) + f_pmp*st%delta(j,i)*st%p_b(j,i)
        else
          st%p_b_red(j,i) = st%delta(j,i)*st%p_b(j,i)
        endif
        st%p_b_red_lim(j,i) = st%p_b_red(j,i)
        st%p_b_w(j,i)       = st%p_b(j,i) - st%p_b_red(j,i) 

      else if (par%i_p_b_red.eq.3) then
        ! reduced basal pressure = delta * overburden pressure (where temperate base) - water pressure below sea level 

        if (st%n_cts(j,i) == -1) then   ! cold ice base
          if (st%sub_melt_flag(j,i)) then
            temp_diff = max((st%temp_c_m(0,j,i)-st%temp_c(0,j,i)), 0.0_wp)
            f_pmp = exp(-st%gamma_slide_inv(j,i)*temp_diff)  ! 'virtual' fraction at pressure melting point
          else
            f_pmp = 0._wp
          end if
          st%p_b_red(j,i) = min((1._wp-f_pmp)*st%p_b(j,i) + f_pmp*st%delta(j,i)*st%p_b(j,i), st%p_b(j,i)- RHO_SW*G*max((st%z_sl(j,i)-st%zb(j,i)), 0.0_wp))
          st%p_b_red(j,i) = max(par%red_pres_limit_fact*st%p_b(j,i), st%p_b_red(j,i))
        else
          st%p_b_red(j,i) = min(st%delta(j,i)*st%p_b(j,i), st%p_b(j,i) - RHO_SW*G*max((st%z_sl(j,i)-st%zb(j,i)), 0.0_wp))
          st%p_b_red(j,i) = max(par%red_pres_limit_fact*st%p_b(j,i), st%p_b_red(j,i))
        endif
        st%p_b_red_lim(j,i) = st%p_b_red(j,i)
        st%p_b_w(j,i)       = st%p_b(j,i) - st%p_b_red(j,i) 

      else if (par%i_p_b_red.eq.4) then
        ! basal water pressure computed from water content in the sediment (rock) + water pressure below sea level

        st%p_b_w(j,i) = max((1._wp-st%delta(j,i))*st%p_b(j,i)*(st%H_w(j,i)/par%H_w_max), RHO_SW*G*max((st%z_sl(j,i)-st%zb(j,i)), 0.0_wp)) 
        st%p_b_red(j,i)     = max(par%red_pres_limit_fact*st%p_b(j,i), st%p_b(j,i)-st%p_b_w(j,i))
        st%p_b_red_lim(j,i) = st%p_b_red(j,i)
      endif

   else   ! maske(j,i) == 1, 2 or 3 away from the grounding line

      st%p_b(j,i)         = 0.0_wp
      st%p_b_w(j,i)       = 0.0_wp
      st%p_b_red(j,i)     = 0.0_wp
      st%p_b_red_lim(j,i) = 0.0_wp

   end if

end do
end do

!  ------ Absolute value of the basal shear stress, tau_b

do i=0, grd%IMAX
do j=0, grd%JMAX
      st%tau_b(j,i) = st%p_b(j,i)*sqrt(st%dzs_dxi_g(j,i)**2+st%dzs_deta_g(j,i)**2)
end do
end do

!-------- Computation of d_help_b (defined on the grid points (i,j)) --------

do i=0, grd%IMAX
do j=0, grd%JMAX

   if ((st%maske(j,i) == 0).or.st%flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

if (par%slide_law==1) then
  cvxy1 = st%c_slide(j,i) &
    * ((st%tau_b(j,i)+eps_wp)**(st%p_weert(j,i)-1)/(st%p_b(j,i)+eps_wp)**st%q_weert(j,i)) &
    * st%p_b(j,i)
  ctau1 = 1.0_wp/(st%c_slide(j,i)+eps_wp)**st%p_weert_inv(j,i) &
    * (st%p_b(j,i)+eps_wp)**(st%q_weert(j,i)*st%p_weert_inv(j,i))
else if (par%slide_law==2) then
  cvxy1 = st%c_slide(j,i) &
    * ((st%tau_b(j,i)+eps_wp)**(st%p_weert(j,i)-1)/(st%p_b_red_lim(j,i)+eps_wp)**st%q_weert(j,i)) &
    * st%p_b(j,i)
  ctau1 = 1.0_wp/(st%c_slide(j,i)+eps_wp)**st%p_weert_inv(j,i) &
    * (st%p_b_red_lim(j,i)+eps_wp)**(st%q_weert(j,i)*st%p_weert_inv(j,i))
else if (par%slide_law==3) then
  cvxy1 = st%c_slide(j,i) &
    * ((st%tau_b(j,i)+eps_wp)**(st%p_weert(j,i)-1)/(st%p_b_red_lim(j,i)+eps_wp)**st%q_weert(j,i)) &
    * st%p_b(j,i)
  ctau1 = 1.0_wp/(st%c_slide(j,i)+eps_wp)**st%p_weert_inv(j,i) &
    * (st%p_b_red(j,i)+eps_wp)**(st%q_weert(j,i)*st%p_weert_inv(j,i))
else if (par%slide_law==4) then
  ! regolarized Coulomb friction law, dummy
  cvxy1 = 0._wp     ! dummy value
  ctau1 = 0._wp     ! dummy value 
else if (par%slide_law==5) then
  ! pseudo-plastic power-law friction law, dummy
  cvxy1 = 0._wp     ! dummy value
  ctau1 = 0._wp     ! dummy value 
endif

!  ------ d_help_b, c_drag

      if (st%n_cts(j,i) == -1) then   ! cold ice base

         if (st%sub_melt_flag(j,i)) then
            temp_diff = max((st%temp_c_m(0,j,i)-st%temp_c(0,j,i)), 0.0_wp)
            cvxy1a    = exp(-st%gamma_slide_inv(j,i)*temp_diff)  ! sub-melt sliding
            ctau1a    = 1.0_wp/(cvxy1a+eps_wp)**st%p_weert_inv(j,i)
         else
            cvxy1a    = 0.0_wp   ! no sub-melt sliding
            ctau1a    = 1.0_wp/eps_wp**st%p_weert_inv(j,i)   ! dummy value
         end if

         st%d_help_b(j,i)   = cvxy1*cvxy1a
         st%c_drag(j,i)  = ctau1*ctau1a

      else if (st%n_cts(j,i) == 0) then   ! temperate ice base

         st%d_help_b(j,i)   = cvxy1   ! basal sliding
         st%c_drag(j,i)  = ctau1   ! (pressure-melting conditions)

      else   ! st%n_cts(j,i) == 1, temperate ice layer

         st%d_help_b(j,i)   = cvxy1   ! basal sliding
         st%c_drag(j,i)  = ctau1   ! (pressure-melting conditions)

      end if

!    ---- Contribution of the basal water layer

if (par%basal_hydrology==1) then

  Hw0_slide_inv = 1.0_wp/max(par%Hw0_slide, eps_wp)

  if (par%hydro_slide_sat_fct==0) then
    ! exponential saturation function
    ! by Kleiner and Humbert (2014, J. Glaciol. 60)

    ratio_Hw_slide = max(st%H_w(j,i)*Hw0_slide_inv, 0.0_wp)
    ! constrain to interval [0,infty)
    cvxy1b = 1.0_wp + par%c_Hw_slide*(1.0_wp-exp(-ratio_Hw_slide))

  else if (par%hydro_slide_sat_fct==1) then
    ! linear saturation function

    ratio_Hw_slide = max(min(st%H_w(j,i)*Hw0_slide_inv, 1.0_wp), 0.0_wp)
    ! constrain to interval [0,1]
    cvxy1b = 1.0_wp + par%c_Hw_slide*ratio_Hw_slide

  else if (par%hydro_slide_sat_fct==2) then
    ! cubic S-shape saturation function

    ratio_Hw_slide = max(min(st%H_w(j,i)*Hw0_slide_inv, 1.0_wp), 0.0_wp)
    ! constrain to interval [0,1]
    cvxy1b = 1.0_wp + par%c_Hw_slide*ratio_Hw_slide*ratio_Hw_slide*(3.0_wp-2.0_wp*ratio_Hw_slide)

  else if (par%hydro_slide_sat_fct==3) then
    ! quintic S-shape saturation function

    ratio_Hw_slide = max(min(st%H_w(j,i)*Hw0_slide_inv, 1.0_wp), 0.0_wp)
    ! constrain to interval [0,1]
    cvxy1b = 1.0_wp + par%c_Hw_slide*ratio_Hw_slide*ratio_Hw_slide*ratio_Hw_slide*(10.0_wp + ratio_Hw_slide *(-15.0_wp+6.0_wp*ratio_Hw_slide))

  endif

  ctau1b = 1.0_wp/(cvxy1b+eps_wp)**st%p_weert_inv(j,i)


  st%d_help_b(j,i) = st%d_help_b(j,i)  *cvxy1b
  st%c_drag(j,i)   = st%c_drag(j,i) *ctau1b

endif

   else   ! maske(j,i) == 1, 2 or 3 away from the grounding line

      st%d_help_b(j,i) = 0.0_wp
      st%c_drag(j,i)   = 0.0_wp

   end if

end do
end do

!-------- Computation of vx_b (defined at (i+1/2,j)) --------

do i=0, grd%IMAX-1
do j=1, grd%JMAX-1
   st%vx_b(j,i) = -0.5_wp*(st%d_help_b(j,i)+st%d_help_b(j,i+1))*st%dzs_dx_aux(j,i)
end do
end do

!-------- Computation of vy_b (defined at (i,j+1/2)) --------

do i=1, grd%IMAX-1
do j=0, grd%JMAX-1
   st%vy_b(j,i) = -0.5_wp*(st%d_help_b(j,i)+st%d_help_b(j+1,i))*st%dzs_dy_aux(j,i)
end do
end do

!-------- Computation of vx_b_g and vy_b_g (defined at (i,j)) --------

do i=0, grd%IMAX
do j=0, grd%JMAX
   st%vx_b_g(j,i) = -st%d_help_b(j,i)*st%dzs_dxi_g(j,i)
   st%vy_b_g(j,i) = -st%d_help_b(j,i)*st%dzs_deta_g(j,i)
end do
end do

!-------- Limitation of computed vx_b, vy_b, vx_b_g, vy_b_g to the interval
!         [-VH_MAX, VH_MAX] --------

vh_max     = max(par%vh_max, eps_wp)/sec_year
vh_max_inv = 1.0_wp/vh_max

call velocity_limiter_gradual(st%vx_b, vh_max, vh_max_inv)
call velocity_limiter_gradual(st%vy_b, vh_max, vh_max_inv)

call velocity_limiter_gradual(st%vx_b_g, vh_max, vh_max_inv)
call velocity_limiter_gradual(st%vy_b_g, vh_max, vh_max_inv)

!-------- Discard basal velocities for HYB_MODE==1 --------

if (par%dynamics==2 .and. par%hyb_mode==1) then

st%d_help_b = 0.0_wp
st%vx_b     = 0.0_wp
st%vy_b     = 0.0_wp
st%vx_b_g   = 0.0_wp
st%vy_b_g   = 0.0_wp

! c_slide and c_drag are not reset because they will be used in the
! computation of the SStA velocity components

endif

end subroutine calc_vxy_b_sia

!-------------------------------------------------------------------------------
!> Computation of the shear stresses txz, tyz, the effective shear stress
!! sigma, the depth-averaged fluidity flui_ave_sia, the horizontal
!! velocity vx, vy and the horizontal volume flux qx, qy in the shallow ice
!! approximation.
!<------------------------------------------------------------------------------
subroutine calc_vxy_sia(st,grd,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(in) :: par

integer :: i, j, kc, kt
real(wp), dimension(:,:,:), allocatable :: d_help_c
real(wp), dimension(:,:,:), allocatable :: d_help_t
real(wp), dimension(:,:,:), allocatable :: ctxyz1 
real(wp), dimension(:,:,:), allocatable :: ctxyz2
real(wp) :: flui_t(0:100), flui_c(0:100)
real(wp) :: cflui0(0:100), cflui1(0:100)
real(wp) :: cvxy2(0:100), cvxy3(0:100)
real(wp) :: cqxy0(0:100), cqxy1(0:100)
real(wp) :: vh_max, vh_max_inv
real(wp) :: flui_min, flui_max, flui_init
real(wp) :: ratio_sl_threshold
real(wp) :: ratef

allocate(d_help_c(0:grd%KCMAX,0:grd%JMAX,0:grd%IMAX))
allocate(d_help_t(0:grd%KTMAX,0:grd%JMAX,0:grd%IMAX))

flui_min  = 1.0_wp/par%visc_max
flui_max  = 1.0_wp/par%visc_min
flui_init = 1.0_wp/par%visc_init_ssa

vh_max     = max(par%vh_max, eps_wp)/sec_year
vh_max_inv = 1.0_wp/vh_max

!-------- Computation of stresses --------

!  ------ Term abbreviations
allocate(ctxyz1(0:grd%KCMAX,0:grd%JMAX,0:grd%IMAX))
allocate(ctxyz2(0:grd%KTMAX,0:grd%JMAX,0:grd%IMAX))

!$omp parallel do collapse(2) private(i,j,kc,kt)
do i=0, grd%IMAX
do j=0, grd%JMAX

   if ((st%maske(j,i) == 0).or.st%flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

      do kc=0, grd%KCMAX
         ctxyz1(kc,j,i) = RHO*G*st%H_c(j,i)*(1.0_wp-grd%eaz_c_quotient(kc))
      end do

      if (st%n_cts(j,i) == 1) then   ! temperate layer

         do kt=0, grd%KTMAX
            ctxyz2(kt,j,i) = RHO*G*st%H_t(j,i)*(1.0_wp-grd%zeta_t(kt))
         end do

      else   ! cold base (-1), temperate base (0)

         do kt=0, grd%KTMAX
            ctxyz2(kt,j,i) = 0.0_wp
         end do

      end if

   else   ! maske(j,i) == 1, 2 or 3 away from the grounding line

      do kc=0, grd%KCMAX
         ctxyz1(kc,j,i) = 0.0_wp
      end do

      do kt=0, grd%KTMAX
         ctxyz2(kt,j,i) = 0.0_wp
      end do

   end if

end do
end do
!$omp end parallel do

!$omp parallel do collapse(2) private(i,j,kc,kt,flui_t,cflui0,flui_c,cflui1,ratef,cvxy2,cvxy3)
do i=0, grd%IMAX
do j=0, grd%JMAX

  !  ------ Shear stress txz (defined at (i+1/2,j,kc/t))

  if (i.lt.grd%IMAX) then
   do kc=0, grd%KCMAX
      st%txz_c(kc,j,i) = -0.5_wp*(ctxyz1(kc,j,i)+ctxyz1(kc,j,i+1)) &
                      *st%dzs_dx_aux(j,i)
   end do

   do kt=0, grd%KTMAX
      st%txz_t(kt,j,i) = st%txz_c(0,j,i) &
                      -0.5_wp*(ctxyz2(kt,j,i)+ctxyz2(kt,j,i+1)) &
                      *st%dzs_dx_aux(j,i)
   end do
 endif

  !  ------ Shear stress tyz (defined at (i,j+1/2,kc/t))

  if (j.lt.grd%JMAX) then
   do kc=0, grd%KCMAX
      st%tyz_c(kc,j,i) = -0.5_wp*(ctxyz1(kc,j,i)+ctxyz1(kc,j+1,i)) &
                      *st%dzs_dy_aux(j,i)
   end do

   do kt=0, grd%KTMAX
      st%tyz_t(kt,j,i) = st%tyz_c(0,j,i) &
                      -0.5_wp*(ctxyz2(kt,j,i)+ctxyz2(kt,j+1,i)) &
                      *st%dzs_dy_aux(j,i)
   end do
 endif
  
 !  ------ Effective shear stress sigma (defined at (i,j,kc/t))

   do kc=0, grd%KCMAX
      st%sigma_c(kc,j,i) = ctxyz1(kc,j,i) &
           *sqrt(st%dzs_dxi_g(j,i)**2+st%dzs_deta_g(j,i)**2)
   end do

   do kt=0, grd%KTMAX
      st%sigma_t(kt,j,i) = st%sigma_c(0,j,i) &
         + ctxyz2(kt,j,i) &
           *sqrt(st%dzs_dxi_g(j,i)**2+st%dzs_deta_g(j,i)**2)
   end do

!-------- Computation of the depth-averaged fluidity
!                 (defined on the grid points (i,j)) --------

   if ((st%maske(j,i) == 0).or.st%flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Fluidity, abbreviations

      do kt=0, grd%KTMAX
         flui_t(kt) = 2.0_wp &
                      *st%enh_t(kt,j,i) &
                      *ratefac_t(st%omega_t(kt,j,i)) &
                      *creep(par%fin_visc,par%flow_law,st%sigma_t(kt,j,i))
         cflui0(kt) = st%H_t(j,i)*flui_t(kt)*grd%dzeta_t
      end do

      do kc=0, grd%KCMAX
         if (par%calcmod==0 .or. par%calcmod==1 .or. par%calcmod==-1) then
           ratef = ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i))
         else if (par%calcmod==2 .or. par%calcmod==3) then
           ratef = ratefac_c_t(st%temp_c(kc,j,i), st%omega_c(kc,j,i),st%temp_c_m(kc,j,i))
         endif
         flui_c(kc) = 2.0_wp &
                      *st%enh_c(kc,j,i) &
                      *ratef &
                      *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i))
         cflui1(kc) = grd%aqxy1(kc)*st%H_c(j,i)*flui_c(kc)
      end do

!  ------ Depth average

      st%flui_ave_sia(j,i) = 0.0_wp

      if (st%n_cts(j,i) == 1) then

         do kt=0, grd%KTMAX-1
            st%flui_ave_sia(j,i) = st%flui_ave_sia(j,i)+0.5_wp*(cflui0(kt+1)+cflui0(kt))
         end do

      end if

      do kc=0, grd%KCMAX-1
         st%flui_ave_sia(j,i) = st%flui_ave_sia(j,i)+0.5_wp*(cflui1(kc+1)+cflui1(kc))
      end do

      st%flui_ave_sia(j,i) = st%flui_ave_sia(j,i)/max((st%H_c(j,i)+st%H_t(j,i)), eps_wp)
      st%flui_ave_sia(j,i) = max(min(st%flui_ave_sia(j,i), flui_max), flui_min)

   else   ! maske(j,i) == 1, 2 or 3 away from the grounding line

      st%flui_ave_sia(j,i) = flui_init

   end if

!-------- Computation of d_help_c/t
!         (defined on the grid points (i,j,kc/t)) --------

   if ((st%maske(j,i) == 0).or.st%flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

      do kt=0, grd%KTMAX
         cvxy2(kt) = 2.0_wp*st%H_t(j,i) &
                     *st%enh_t(kt,j,i) &
                     *ratefac_t(st%omega_t(kt,j,i)) &
                     *creep(par%fin_visc,par%flow_law,st%sigma_t(kt,j,i)) &
                     *(ctxyz1(0,j,i)+ctxyz2(kt,j,i)) &
                     *grd%dzeta_t
      end do

      do kc=0, grd%KCMAX
         if (par%calcmod==0 .or. par%calcmod==1 .or. par%calcmod==-1) then
           ratef = ratefac_c(st%temp_c(kc,j,i), st%temp_c_m(kc,j,i))
         else if (par%calcmod==2 .or. par%calcmod==3) then
           ratef = ratefac_c_t(st%temp_c(kc,j,i), st%omega_c(kc,j,i),st%temp_c_m(kc,j,i))
         endif
         cvxy3(kc) = 2.0_wp*grd%avxy3(kc)*st%H_c(j,i) &
                     *st%enh_c(kc,j,i) &
                     *ratef &
                     *creep(par%fin_visc,par%flow_law,st%sigma_c(kc,j,i)) &
                     *ctxyz1(kc,j,i)
      end do

!  ------ d_help_c, d_help_t

      if (st%n_cts(j,i) == -1) then   ! cold ice base

         do kt=0, grd%KTMAX
            d_help_t(kt,j,i) = st%d_help_b(j,i)
         end do

         d_help_c(0,j,i) = d_help_t(grd%KTMAX,j,i)

         do kc=0, grd%KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_wp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      else if (st%n_cts(j,i) == 0) then   ! temperate ice base

         do kt=0, grd%KTMAX
            d_help_t(kt,j,i) = st%d_help_b(j,i)
         end do

         d_help_c(0,j,i) = d_help_t(grd%KTMAX,j,i)

         do kc=0, grd%KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_wp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      else   ! n_cts(j,i) == 1, temperate ice layer

         d_help_t(0,j,i) = st%d_help_b(j,i)

         do kt=0, grd%KTMAX-1
            d_help_t(kt+1,j,i) = d_help_t(kt,j,i) &
                                +0.5_wp*(cvxy2(kt+1)+cvxy2(kt))
         end do

         d_help_c(0,j,i) = d_help_t(grd%KTMAX,j,i)

         do kc=0, grd%KCMAX-1
            d_help_c(kc+1,j,i) = d_help_c(kc,j,i) &
                                +0.5_wp*(cvxy3(kc+1)+cvxy3(kc))
         end do

      end if

   else   ! maske(j,i) == 1, 2 or 3 away from the grounding line

      do kt=0, grd%KTMAX
         d_help_t(kt,j,i) = 0.0_wp
      end do

      do kc=0, grd%KCMAX
         d_help_c(kc,j,i) = 0.0_wp
      end do

   end if

end do
end do
!$omp end parallel do

!$omp parallel do collapse(2) private(i,j,kt,kc,cqxy0,cqxy1)
do i=0, grd%IMAX
do j=0, grd%JMAX

  !-------- Computation of vx_c/t (defined at (i+1/2,j,kc/t)) --------

  if (j>0 .and. i<grd%IMAX .and. j<grd%JMAX) then
   do kt=0, grd%KTMAX
      st%vx_t(kt,j,i) = -0.5_wp*(d_help_t(kt,j,i)+d_help_t(kt,j,i+1)) &
                     *st%dzs_dx_aux(j,i)
   end do

   do kc=0, grd%KCMAX
      st%vx_c(kc,j,i) = -0.5_wp*(d_help_c(kc,j,i)+d_help_c(kc,j,i+1)) &
                     *st%dzs_dx_aux(j,i)
   end do
 endif

  !-------- Computation of vy_c/t (defined at (i,j+1/2,kc/t)) --------

  if (i>0 .and. i<grd%IMAX .and. j<grd%JMAX) then
   do kt=0, grd%KTMAX
      st%vy_t(kt,j,i) = -0.5_wp*(d_help_t(kt,j,i)+d_help_t(kt,j+1,i)) &
                     *st%dzs_dy_aux(j,i)
   end do

   do kc=0, grd%KCMAX
      st%vy_c(kc,j,i) = -0.5_wp*(d_help_c(kc,j,i)+d_help_c(kc,j+1,i)) &
                     *st%dzs_dy_aux(j,i)
   end do
 endif

  !-------- Computation of the surface velocities vx_s_g and vy_s_g (defined at (i,j)) --------

   st%vx_s_g(j,i) = -d_help_c(grd%KCMAX,j,i)*st%dzs_dxi_g(j,i)
   st%vy_s_g(j,i) = -d_help_c(grd%KCMAX,j,i)*st%dzs_deta_g(j,i)

!-------- Limitation of computed vx_c/t, vy_c/t, vx_s_g, vy_s_g
!         to the interval [-VH_MAX, VH_MAX] --------

call velocity_limiter_gradual(st%vx_s_g(j,i), vh_max, vh_max_inv)
call velocity_limiter_gradual(st%vy_s_g(j,i), vh_max, vh_max_inv)
   do kt=0, grd%KTMAX
call velocity_limiter_gradual(st%vx_t(kt,j,i), vh_max, vh_max_inv)
call velocity_limiter_gradual(st%vy_t(kt,j,i), vh_max, vh_max_inv)
enddo
   do kc=0, grd%KCMAX
call velocity_limiter_gradual(st%vx_c(kc,j,i), vh_max, vh_max_inv)
call velocity_limiter_gradual(st%vy_c(kc,j,i), vh_max, vh_max_inv)
enddo

!-------- Computation of h_diff
!         (defined on the grid points (i,j)) --------

   if ((st%maske(j,i) == 0).or.st%flag_grounding_line_2(j,i)) then
                     ! grounded ice, or floating ice at the grounding line

!  ------ Abbreviations

      do kt=0, grd%KTMAX
         cqxy0(kt) = st%H_t(j,i)*d_help_t(kt,j,i)*grd%dzeta_t
      end do

      do kc=0, grd%KCMAX
         cqxy1(kc) = grd%aqxy1(kc)*st%H_c(j,i)*d_help_c(kc,j,i)
      end do

!  ------ h_diff

      st%h_diff(j,i) = 0.0_wp

      if (st%n_cts(j,i) == 1) then

         do kt=0, grd%KTMAX-1
            st%h_diff(j,i) = st%h_diff(j,i)+0.5_wp*(cqxy0(kt+1)+cqxy0(kt))
         end do

      end if

      do kc=0, grd%KCMAX-1
         st%h_diff(j,i) = st%h_diff(j,i)+0.5_wp*(cqxy1(kc+1)+cqxy1(kc))
      end do

!  ------ Limitation of h_diff

      if (st%h_diff(j,i) < par%hd_min) st%h_diff(j,i) = 0.0_wp
      if (st%h_diff(j,i) > par%hd_max) st%h_diff(j,i) = par%hd_max

   else   ! maske(j,i) == 1, 2 or 3 away from the grounding line

      st%h_diff(j,i) = 0.0_wp

   end if

end do
end do
!$omp end parallel do

!-------- Computation of the horizontal volume flux
!                            and the depth-averaged velocity --------

do i=0, grd%IMAX-1
do j=0, grd%JMAX

   st%qx(j,i) = -0.5_wp*(st%h_diff(j,i)+st%h_diff(j,i+1))*st%dzs_dx_aux(j,i)

   if ( (st%maske(j,i)==0).or.(st%maske(j,i+1)==0) ) then   ! at least one neighbour
                                                      ! point is grounded ice

      st%vx_m(j,i)       = st%qx(j,i) &
                        / ( 0.5_wp*(st%H_c(j,i)+st%H_t(j,i)+st%H_c(j,i+1)+st%H_t(j,i+1)) )

      call velocity_limiter_gradual(st%vx_m(j,i), vh_max, vh_max_inv)

      st%ratio_sl_x(j,i) = abs(st%vx_t(0,j,i)) / max(abs(st%vx_c(grd%KCMAX,j,i)), eps_wp)

   else 

      st%vx_m(j,i)       = 0.0_wp
      st%ratio_sl_x(j,i) = 0.0_wp

   end if

end do
end do

do i=0, grd%IMAX
do j=0, grd%JMAX-1

   st%qy(j,i) = -0.5_wp*(st%h_diff(j,i)+st%h_diff(j+1,i))*st%dzs_dy_aux(j,i)

   if ( (st%maske(j,i)==0).or.(st%maske(j+1,i)==0) ) then   ! at least one neighbour
                                                      ! point is grounded ice

      st%vy_m(j,i)       = st%qy(j,i) &
                        / ( 0.5_wp*(st%H_c(j,i)+st%H_t(j,i)+st%H_c(j+1,i)+st%H_t(j+1,i)) )

      call velocity_limiter_gradual(st%vy_m(j,i), vh_max, vh_max_inv)

      st%ratio_sl_y(j,i) = abs(st%vy_t(0,j,i)) / max(abs(st%vy_c(grd%KCMAX,j,i)), eps_wp)

   else 

      st%vy_m(j,i)       = 0.0_wp
      st%ratio_sl_y(j,i) = 0.0_wp

   end if

end do
end do

st%ratio_sl = 0.0_wp

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1

   if (st%maske(j,i) == 0) &   ! grounded ice
      st%ratio_sl(j,i) = 0.25_wp &
                        * (   st%ratio_sl_x(j,i-1) + st%ratio_sl_x(j,i) &
                            + st%ratio_sl_y(j-1,i) + st%ratio_sl_y(j,i) )
end do
end do

!-------- Detection of shelfy stream points --------

st%flag_shelfy_stream_x = .false.
st%flag_shelfy_stream_y = .false.
st%flag_shelfy_stream   = .false.

if (par%dynamics==0 .or. par%dynamics==1) then

  ratio_sl_threshold = 1.11e+11_wp   ! dummy value

else if (par%dynamics==2) then

  ratio_sl_threshold = par%ratio_sl_thresh 

  do i=0, grd%IMAX-1
    do j=0, grd%JMAX
      if (par%hyb_mode==0) then
        if (st%ratio_sl_x(j,i) > ratio_sl_threshold) st%flag_shelfy_stream_x(j,i) = .true.
      else if (par%hyb_mode==1) then
        if ( (st%maske(j,i)==0).or.(st%maske(j,i+1)==0) ) st%flag_shelfy_stream_x(j,i) = .true.
      endif
    end do
  end do

  do i=0, grd%IMAX
    do j=0, grd%JMAX-1
      if (par%hyb_mode==0) then
        if (st%ratio_sl_y(j,i) > ratio_sl_threshold) st%flag_shelfy_stream_y(j,i) = .true.
      else if (par%hyb_mode==1) then
        if ((st%maske(j,i)==0).or.(st%maske(j+1,i)==0)) st%flag_shelfy_stream_y(j,i) = .true.
      endif
    end do
  end do

  do i=1, grd%IMAX-1
    do j=1, grd%JMAX-1

      if ( (st%maske(j,i) == 0) &   ! grounded ice
        .and. &
        (     st%flag_shelfy_stream_x(j,i-1)   &   ! at least
        .or.st%flag_shelfy_stream_x(j,i)     &   ! one neighbour
        .or.st%flag_shelfy_stream_y(j-1,i)   &   ! on the staggered grid
        .or.st%flag_shelfy_stream_y(j,i)   ) &   ! is a shelfy stream point
        ) then

        st%flag_shelfy_stream(j,i) = .true.

      end if

    end do
  end do

else

  stop ' >>> calc_vxy_sia: DYNAMICS must be 0, 1 or 2!'

endif

!-------- Save mean (depth-averaged) horizontal velocities from SIA --------

st%vx_m_sia = st%vx_m
st%vy_m_sia = st%vy_m

!-------- Initialisation of the variable q_gl_g
!         (volume flux across the grounding line, to be
!         computed in the routine calc_vxy_ssa
!         if ice shelves are present)

st%q_gl_g = 0.0_wp

deallocate(d_help_c)
deallocate(d_help_t)
deallocate(ctxyz1)
deallocate(ctxyz2)

end subroutine calc_vxy_sia

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy, the horizontal volume flux
!> qx, qy etc. for static ice.
!<------------------------------------------------------------------------------
subroutine calc_vxy_static(st,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_par_class), intent(in) :: par

real(wp) :: flui_init


flui_init = 1.0_wp/par%visc_init_ssa 

st%c_slide = 0.0_wp
st%p_weert = 0
st%q_weert = 0
st%p_b_w   = 0.0_wp

st%c_drag   = 0.0_wp

st%vx_b   = 0.0_wp
st%vy_b   = 0.0_wp
st%vx_b_g = 0.0_wp
st%vy_b_g = 0.0_wp

!st%txz_c = 0.0_wp
!st%txz_t = 0.0_wp
!
!st%tyz_c = 0.0_wp
!st%tyz_t = 0.0_wp

st%sigma_c = 0.0_wp
st%sigma_t = 0.0_wp

st%flui_ave_sia = flui_init
st%de_ssa       = 0.0_wp
st%vis_int_g    = 0.0_wp

st%vx_c = 0.0_wp
st%vy_c = 0.0_wp

st%vx_t = 0.0_wp
st%vy_t = 0.0_wp

st%vx_s_g = 0.0_wp
st%vy_s_g = 0.0_wp

st%h_diff = 0.0_wp

st%qx = 0.0_wp
st%qy = 0.0_wp

st%vx_m = 0.0_wp
st%vy_m = 0.0_wp

st%ratio_sl_x = 0.0_wp
st%ratio_sl_y = 0.0_wp

st%flag_shelfy_stream_x = .false.
st%flag_shelfy_stream_y = .false.
st%flag_shelfy_stream   = .false.

st%vx_m_sia = 0.0_wp
st%vy_m_sia = 0.0_wp

st%q_gl_g = 0.0_wp

end subroutine calc_vxy_static

!-------------------------------------------------------------------------------
!> Computation of the horizontal velocity vx, vy, the horizontal volume flux
!! qx, qy and the flux across the grounding line q_gl_g in the shallow shelf
!! approximation (SSA) or the shelfy stream approximation (SStA).
!<------------------------------------------------------------------------------
subroutine calc_vxy_ssa(st,grd,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(in) :: par

integer :: i, j, m, kc, kt
real(wp) :: res_vxy_m_ssa_1, res_vxy_m_ssa_2, res_vxy_m_ssa
real(wp) :: vh_max, vh_max_inv
real(wp) :: ratio_sl_threshold, ratio_help, v_ref
real(wp) :: qx_gl_g, qy_gl_g

    !$ real(wp) :: time1,time2

!-------- Parameters for the relaxation scheme --------

vh_max     = max(par%vh_max, eps_wp)/sec_year
vh_max_inv = 1.0_wp/vh_max

!write(6,'(10x,a)') 'calc_vxy_ssa:'

!-------- Iterations --------
    
!$ if(par%l_write_timer) print *

res_vxy_m_ssa = 1.11e+11_wp   ! initial, very large value of the residual

m=0

do while ((m < par%n_iter_ssa).and.(res_vxy_m_ssa > par%tol_iter_ssa))

  m = m+1

!  write(6,'(13x,a,i0,a)', advance='no') 'Iter ', m, ': '

  !  ------ Save velocities from previous iteration

  st%vx_m_prev = st%vx_m_ssa
  st%vy_m_prev = st%vy_m_ssa

  !  ------ Depth-integrated viscosity vis_int_g

  if (m > 1) then

    !$ time1 = omp_get_wtime()
    call calc_vis_ssa(st,grd,par)
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'m'
    !$ if(par%l_write_timer) print *,'vis_ssa',time2-time1

  else   ! (m == 1, first iteration)

    if (par%iter_init_ssa==1) then
      ! constant viscosity times ice thickness
      st%vis_int_g = (st%H_c+st%H_t)*par%visc_init_ssa
    else if (par%iter_init_ssa==2) then
      ! previous depth-averaged viscosity times ice thickness
      st%vis_int_g = (st%H_c+st%H_t)*st%vis_ave_g
    else if (par%iter_init_ssa==3) then
      ! standard computation by subroutine calc_vis_ssa
      call calc_vis_ssa(st,grd,par)
    endif

  end if

  !  ------ Horizontal velocity vx_m_ssa, vy_m_ssa

  st%flag_calc_vxy_ssa_x = .false.   ! initialization
  st%flag_calc_vxy_ssa_y = .false.   ! initialization

    !$ time1 = omp_get_wtime()
  call calc_vxy_ssa_matrix(st,grd,par)
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'vxy_ssa_matrix',time2-time1

    !$ time1 = omp_get_wtime()
  call velocity_limiter_gradual(st%vx_m_ssa, vh_max, vh_max_inv)
  call velocity_limiter_gradual(st%vy_m_ssa, vh_max, vh_max_inv)

  !  ------ Relaxation scheme

  if (m > 1) then
    st%vx_m_ssa = par%relax_fact_ssa*st%vx_m_ssa + (1.0_wp-par%relax_fact_ssa)*st%vx_m_prev
    st%vy_m_ssa = par%relax_fact_ssa*st%vy_m_ssa + (1.0_wp-par%relax_fact_ssa)*st%vy_m_prev
  end if

  !  ------ Residual

  res_vxy_m_ssa_1 &
    = sqrt( sum((st%vx_m_ssa-st%vx_m_prev)*(st%vx_m_ssa-st%vx_m_prev), mask=st%flag_calc_vxy_ssa_x) &
    +sum((st%vy_m_ssa-st%vy_m_prev)*(st%vy_m_ssa-st%vy_m_prev), mask=st%flag_calc_vxy_ssa_y) )
  res_vxy_m_ssa_2 &
    = sqrt( sum((st%vx_m_ssa+st%vx_m_prev)*(st%vx_m_ssa+st%vx_m_prev), mask=st%flag_calc_vxy_ssa_x) &
    +sum((st%vy_m_ssa+st%vy_m_prev)*(st%vy_m_ssa+st%vy_m_prev), mask=st%flag_calc_vxy_ssa_y) )

  res_vxy_m_ssa = 2.0_wp*res_vxy_m_ssa_1/max(res_vxy_m_ssa_2, eps_wp)
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'residuals',time2-time1

!  write(6,'(a,es9.2)') 'res =', res_vxy_m_ssa

end do

!  ------ Depth-integrated viscosity vis_int_g

call calc_vis_ssa(st,grd,par)

!-------- 3-D velocities, basal velocities and volume flux --------

if (par%dynamics==0 .or. par%dynamics==1) then

  ratio_sl_threshold = 1.11e+11_wp   ! dummy value
  ratio_help         = 0.0_wp
  v_ref              = 1.11e+11_dp   ! dummy value
else if (par%dynamics==2) then

  ratio_sl_threshold = par%ratio_sl_thresh
  ratio_help = 1.0_wp/(1.0_wp-ratio_sl_threshold)
  v_ref = par%hyb_ref_speed/sec_year
endif

!  ------ x-component

do i=0, grd%IMAX-1
  do j=0, grd%JMAX

    st%weigh_ssta_sia_x(j,i) = 0.0_wp

    if (st%flag_shelfy_stream_x(j,i)) then   ! shelfy stream

      if (par%hyb_mode==0) then  ! Ralf's approach

        st%weigh_ssta_sia_x(j,i) = (st%ratio_sl_x(j,i)-ratio_sl_threshold)*ratio_help

        st%weigh_ssta_sia_x(j,i) = max(min(st%weigh_ssta_sia_x(j,i), 1.0_wp), 0.0_wp)
        ! constrain to interval [0,1]

        if (par%ssta_sia_weigh_fct==1) then
          ! make transition smooth (cubic function)

          st%weigh_ssta_sia_x(j,i) = st%weigh_ssta_sia_x(j,i)*st%weigh_ssta_sia_x(j,i)*(3.0_wp-2.0_wp*st%weigh_ssta_sia_x(j,i))

        else if (par%ssta_sia_weigh_fct==2) then
          ! make transition even smoother (quintic function)

          st%weigh_ssta_sia_x(j,i) = st%weigh_ssta_sia_x(j,i)*st%weigh_ssta_sia_x(j,i)*st%weigh_ssta_sia_x(j,i) &
            *(10.0_wp + st%weigh_ssta_sia_x(j,i)*(-15.0_wp+6.0_wp*st%weigh_ssta_sia_x(j,i)))

        endif

        do kt=0, grd%KTMAX
          st%vx_t(kt,j,i) = st%weigh_ssta_sia_x(j,i)*st%vx_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_x(j,i))*st%vx_t(kt,j,i)
        end do

        do kc=0, grd%KCMAX
          st%vx_c(kc,j,i) = st%weigh_ssta_sia_x(j,i)*st%vx_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_x(j,i))*st%vx_c(kc,j,i)
        end do

        st%vx_b(j,i) = st%vx_t(0,j,i)

        st%vx_m(j,i) = st%weigh_ssta_sia_x(j,i)*st%vx_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_x(j,i))*st%vx_m_sia(j,i)

      else if (par%hyb_mode==1) then ! Jorge's approach

        st%weigh_ssta_sia_x(j,i) = (2.0_wp/pi) * atan( (abs(st%vx_m_ssa(j,i))**2) &
          / (v_ref**2) )

        do kt=0, grd%KTMAX
          st%vx_t(kt,j,i) = st%vx_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_x(j,i))*st%vx_t(kt,j,i)
          st%vx_t(kt,j,i) = max(st%vx_t(kt,j,i), -vh_max)
          st%vx_t(kt,j,i) = min(st%vx_t(kt,j,i),  vh_max)
        end do

        do kc=0, grd%KCMAX
          st%vx_c(kc,j,i) = st%vx_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_x(j,i))*st%vx_c(kc,j,i)
          st%vx_c(kc,j,i) = max(st%vx_c(kc,j,i), -vh_max)
          st%vx_c(kc,j,i) = min(st%vx_c(kc,j,i),  vh_max)
        end do

        st%vx_b(j,i) = st%vx_t(0,j,i)

        st%vx_m(j,i) = st%vx_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_x(j,i))*st%vx_m_sia(j,i)
        st%vx_m(j,i) = max(st%vx_m(j,i), -vh_max)
        st%vx_m(j,i) = min(st%vx_m(j,i),  vh_max)

      endif

      st%qx(j,i)   = st%vx_m(j,i) * 0.5_wp * ( (st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j,i+1)+st%H_t(j,i+1)) )

    else if (st%flag_calc_vxy_ssa_x(j,i)) then   ! floating ice

      do kt=0, grd%KTMAX
        st%vx_t(kt,j,i) = st%vx_m_ssa(j,i)
      end do

      do kc=0, grd%KCMAX
        st%vx_c(kc,j,i) = st%vx_m_ssa(j,i)
      end do

      st%vx_b(j,i) = st%vx_m_ssa(j,i)

      st%vx_m(j,i) = st%vx_m_ssa(j,i)

      st%qx(j,i)   = st%vx_m(j,i) * 0.5_wp * ( (st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j,i+1)+st%H_t(j,i+1)) )

      !  else
      !     In all other cases, the depth-averaged velocities vx_m_ssa(j,i) computed
      !     by the SSA/SStA solver are discarded.

    end if

  end do
end do

!  ------ y-component

do i=0, grd%IMAX
  do j=0, grd%JMAX-1

    st%weigh_ssta_sia_y(j,i) = 0.0_wp

    if (st%flag_shelfy_stream_y(j,i)) then   ! shelfy stream

      if (par%hyb_mode==0) then  ! Ralf's approach

        st%weigh_ssta_sia_y(j,i) = (st%ratio_sl_y(j,i)-ratio_sl_threshold)*ratio_help

        st%weigh_ssta_sia_y(j,i) = max(min(st%weigh_ssta_sia_y(j,i), 1.0_wp), 0.0_wp)
        ! constrain to interval [0,1]

        if (par%ssta_sia_weigh_fct==1) then
          ! make transition smooth (cubic function)

          st%weigh_ssta_sia_y(j,i) = st%weigh_ssta_sia_y(j,i)*st%weigh_ssta_sia_y(j,i)*(3.0_wp-2.0_wp*st%weigh_ssta_sia_y(j,i))

        else if (par%ssta_sia_weigh_fct==2) then
          ! make transition even smoother (quintic function)

          st%weigh_ssta_sia_y(j,i) = st%weigh_ssta_sia_y(j,i)*st%weigh_ssta_sia_y(j,i)*st%weigh_ssta_sia_y(j,i) &
            *(10.0_wp + st%weigh_ssta_sia_y(j,i)*(-15.0_wp+6.0_wp*st%weigh_ssta_sia_y(j,i)))

        endif

        do kt=0, grd%KTMAX
          st%vy_t(kt,j,i) = st%weigh_ssta_sia_y(j,i)*st%vy_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_y(j,i))*st%vy_t(kt,j,i)
        end do

        do kc=0, grd%KCMAX
          st%vy_c(kc,j,i) = st%weigh_ssta_sia_y(j,i)*st%vy_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_y(j,i))*st%vy_c(kc,j,i)
        end do

        st%vy_b(j,i) = st%vy_t(0,j,i)

        st%vy_m(j,i) = st%weigh_ssta_sia_y(j,i)*st%vy_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_y(j,i))*st%vy_m_sia(j,i)

      else if (par%hyb_mode==1) then ! Jorge's approach

        st%weigh_ssta_sia_y(j,i) = (2.0_wp/pi) * atan( (abs(st%vy_m_ssa(j,i))**2) &
          / (v_ref**2) )

        do kt=0, grd%KTMAX
          st%vy_t(kt,j,i) = st%vy_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_y(j,i))*st%vy_t(kt,j,i)
          st%vy_t(kt,j,i) = max(st%vy_t(kt,j,i), -vh_max)
          st%vy_t(kt,j,i) = min(st%vy_t(kt,j,i),  vh_max)
        end do

        do kc=0, grd%KCMAX
          st%vy_c(kc,j,i) = st%vy_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_y(j,i))*st%vy_c(kc,j,i)
          st%vy_c(kc,j,i) = max(st%vy_c(kc,j,i), -vh_max)
          st%vy_c(kc,j,i) = min(st%vy_c(kc,j,i),  vh_max)
        end do

        st%vy_b(j,i) = st%vy_t(0,j,i)

        st%vy_m(j,i) = st%vy_m_ssa(j,i) + (1.0_wp-st%weigh_ssta_sia_y(j,i))*st%vy_m_sia(j,i)
        st%vy_m(j,i) = max(st%vy_m(j,i), -vh_max)
        st%vy_m(j,i) = min(st%vy_m(j,i),  vh_max)

      endif

      st%qy(j,i)   = st%vy_m(j,i) * 0.5_wp * ( (st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j+1,i)+st%H_t(j+1,i)) )

    else if (st%flag_calc_vxy_ssa_y(j,i)) then   ! floating ice

      do kt=0, grd%KTMAX
        st%vy_t(kt,j,i) = st%vy_m_ssa(j,i)
      end do

      do kc=0, grd%KCMAX
        st%vy_c(kc,j,i) = st%vy_m_ssa(j,i)
      end do

      st%vy_b(j,i) = st%vy_m_ssa(j,i)

      st%vy_m(j,i) = st%vy_m_ssa(j,i)

      st%qy(j,i)   = st%vy_m(j,i) * 0.5_wp * ( (st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j+1,i)+st%H_t(j+1,i)) )

      !  else
      !     In all other cases, the depth-averaged velocities vy_m_ssa(j,i) computed
      !     by the SSA/SStA solver are discarded.

    end if

  end do
end do

!-------- Surface and basal velocities vx_s_g vy_s_g, vx_b_g vy_b_g
!                                                (defined at (i,j)) --------

do i=1, grd%IMAX-1
  do j=1, grd%JMAX-1

    if (st%flag_shelfy_stream(j,i)) then   ! shelfy stream

      st%vx_s_g(j,i) = 0.5_wp*(st%vx_c(grd%KCMAX,j,i-1)+st%vx_c(grd%KCMAX,j,i))
      st%vx_b_g(j,i) = 0.5_wp*(st%vx_b(          j,i-1)+st%vx_b(          j,i))

      st%vy_s_g(j,i) = 0.5_wp*(st%vy_c(grd%KCMAX,j-1,i)+st%vy_c(grd%KCMAX,j,i))
      st%vy_b_g(j,i) = 0.5_wp*(st%vy_b(          j-1,i)+st%vy_b(          j,i))

    else if (st%maske(j,i)==3) then   ! floating ice

      st%vx_s_g(j,i) = 0.5_wp*(st%vx_m(j,i-1)+st%vx_m(j,i))
      st%vx_b_g(j,i) = st%vx_s_g(j,i)

      st%vy_s_g(j,i) = 0.5_wp*(st%vy_m(j-1,i)+st%vy_m(j,i))
      st%vy_b_g(j,i) = st%vy_s_g(j,i)

    end if

  end do
end do

!-------- Computation of the flux across the grounding line q_gl_g

do i=1, grd%IMAX-1
  do j=1, grd%JMAX-1

    if ( st%flag_grounding_line_1(j,i) ) then   ! grounding line

      qx_gl_g = 0.5_wp*(st%qx(j,i)+st%qx(j,i-1))
      qy_gl_g = 0.5_wp*(st%qy(j,i)+st%qy(j-1,i))

      st%q_gl_g(j,i) = sqrt(qx_gl_g*qx_gl_g+qy_gl_g*qy_gl_g)

    end if

  end do
end do

end subroutine calc_vxy_ssa

!-------------------------------------------------------------------------------
!> Solution of the system of linear equations for the horizontal velocities
!! vx_m, vy_m in the shallow shelf approximation.
!<------------------------------------------------------------------------------
subroutine calc_vxy_ssa_matrix(st,grd,par)

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(in) :: par

integer :: i, j, k, n
integer :: i1, j1
real(wp) :: inv_dxi, inv_deta, inv_dxi_deta, inv_dxi2, inv_deta2
real(wp) :: factor_rhs_1, factor_rhs_2, factor_rhs_3a, factor_rhs_3b, rhosw_rho_ratio
real(wp) :: H_mid, zl_mid
real(wp) :: temp_diff, vb
real(wp), parameter :: vb_min = 1.e-3/sec_year  ! m/s
character(len=256) :: ch_solver_set_option

    !$ real(wp) :: time1,time2

! Include header for lis solver fortran interface
#include "lisf.h"

LIS_INTEGER :: ierr
LIS_INTEGER :: nc, nr
LIS_INTEGER :: lin_iter
LIS_MATRIX  :: lgs_a
LIS_VECTOR  :: lgs_b, lgs_x
LIS_SOLVER  :: solver

LIS_INTEGER :: nmax 
LIS_INTEGER :: n_sprs

LIS_INTEGER, allocatable, dimension(:) :: lgs_a_ptr, lgs_a_index
LIS_SCALAR,  allocatable, dimension(:) :: lgs_a_value, lgs_b_value, lgs_x_value

nmax = 2*(grd%IMAX+1)*(grd%JMAX+1)
n_sprs = 20*(grd%IMAX+1)*(grd%JMAX+1)

!-------- Abbreviations --------

inv_dxi      = 1.0_wp/grd%dxi
inv_deta     = 1.0_wp/grd%deta
inv_dxi_deta = 1.0_wp/(grd%dxi*grd%deta)
inv_dxi2     = 1.0_wp/(grd%dxi*grd%dxi)
inv_deta2    = 1.0_wp/(grd%deta*grd%deta)

rhosw_rho_ratio = RHO_SW/RHO

factor_rhs_1  = RHO*G
factor_rhs_2  = 0.5_wp*RHO*G*(RHO_SW-RHO)/RHO_SW
factor_rhs_3a = 0.5_wp*RHO*G
factor_rhs_3b = 0.5_wp*RHO_SW*G

!-------- Depth-integrated viscosity on the staggered grid
!                                       [at (i+1/2,j+1/2)] --------
!$ time1 = omp_get_wtime()

st%vis_int_sgxy = 0.0_wp   ! initialisation

do i=0, grd%IMAX-1
do j=0, grd%JMAX-1

   k=0

   if ((st%maske(j,i)==0).or.(st%maske(j,i)==3)) then
      k = k+1                              ! floating or grounded ice
      st%vis_int_sgxy(j,i) = st%vis_int_sgxy(j,i) + st%vis_int_g(j,i)
   end if

   if ((st%maske(j,i+1)==0).or.(st%maske(j,i+1)==3)) then
      k = k+1                                  ! floating or grounded ice
      st%vis_int_sgxy(j,i) = st%vis_int_sgxy(j,i) + st%vis_int_g(j,i+1)
   end if

   if ((st%maske(j+1,i)==0).or.(st%maske(j+1,i)==3)) then
      k = k+1                                  ! floating or grounded ice
      st%vis_int_sgxy(j,i) = st%vis_int_sgxy(j,i) + st%vis_int_g(j+1,i)
   end if

   if ((st%maske(j+1,i+1)==0).or.(st%maske(j+1,i+1)==3)) then
      k = k+1                                      ! floating or grounded ice
      st%vis_int_sgxy(j,i) = st%vis_int_sgxy(j,i) + st%vis_int_g(j+1,i+1)
   end if

   if (k>0) st%vis_int_sgxy(j,i) = st%vis_int_sgxy(j,i)/real(k,wp)

end do
end do

!-------- Basal drag parameter (for shelfy stream) --------

st%beta_drag = 0.0_wp   ! initialisation

do i=1, grd%IMAX-1
do j=1, grd%JMAX-1

   if (st%flag_shelfy_stream(j,i)) then

     if (par%slide_law.eq.1 .or. par%slide_law.eq.2 .or. par%slide_law.eq.3) then

      st%beta_drag(j,i) = st%c_drag(j,i) &
                     / sqrt( (   (0.5_wp*(st%vx_m_ssa(j,i)+st%vx_m_ssa(j,i-1)))**2  &
                               + (0.5_wp*(st%vy_m_ssa(j,i)+st%vy_m_ssa(j-1,i)))**2 ) &
                               + eps_wp**2 ) &
                                     **(1.0_wp-st%p_weert_inv(j,i))

     else if (par%slide_law.eq.4) then
       ! regularized Coulomb friction law

       ! basal velocity
       vb = sqrt((0.5_wp*(st%vx_m_ssa(j,i)+st%vx_m_ssa(j,i-1)))**2 + (0.5_wp*(st%vy_m_ssa(j,i)+st%vy_m_ssa(j-1,i)))**2 + vb_min**2)

       st%beta_drag(j,i) = st%c_fric(j,i) * st%p_b_red_lim(j,i) * (vb / (vb + st%vb_t(j,i)))**par%q_coulomb * (1._wp/vb)

     else if (par%slide_law.eq.5) then
       ! pseudo-plastic power-law friction law

       ! basal velocity
       vb = sqrt((0.5_wp*(st%vx_m_ssa(j,i)+st%vx_m_ssa(j,i-1)))**2 + (0.5_wp*(st%vy_m_ssa(j,i)+st%vy_m_ssa(j-1,i)))**2 + vb_min**2)

       st%beta_drag(j,i) = st%c_fric(j,i) * st%p_b_red_lim(j,i) * (vb / st%vb_t(j,i))**par%q_coulomb * (1._wp/vb)

     endif

   end if

end do
end do
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'beta_drag',time2-time1

!-------- Assembly of the system of linear equations
!                         (matrix storage: compressed sparse row CSR) --------

allocate(lgs_a_value(n_sprs), lgs_a_index(n_sprs), lgs_a_ptr(nmax+1))
allocate(lgs_b_value(nmax), lgs_x_value(nmax))

lgs_a_value = 0.0_wp
lgs_a_index = 0
lgs_a_ptr   = 0

lgs_b_value = 0.0_wp
lgs_x_value = 0.0_wp

lgs_a_ptr(1) = 1

k = 0

!$ time1 = omp_get_wtime()
!!!$omp parallel do private(n,i,j,i1,j1,nr,nc,H_mid,zl_mid)
do n=1, nmax-1, 2

   i = grd%n2i((n+1)/2)
   j = grd%n2j((n+1)/2)

!  ------ Equations for vx_m_ssa (at (i+1/2,j))

   nr = n   ! row counter

   if ( (i /= grd%IMAX).and.(j /= 0).and.(j /= grd%JMAX) ) then
      ! inner point on the staggered grid in x-direction

      H_mid  = 0.5_wp*((st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j,i+1)+st%H_t(j,i+1)))
      zl_mid = 0.5_wp*(st%zl(j,i)+st%zl(j,i+1))

      if ( &
           ( (st%maske(j,i)==3).and.(st%maske(j,i+1)==3) ) &
           .or. &
           ( st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j,i+1) &
             .and.(H_mid < (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
           .or. &
           ( st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j,i+1) &
             .and.(H_mid < (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
         ) then
           ! both neighbours are floating ice
           !   or
           ! one neighbour is floating ice and the other is grounded ice
           ! (grounding line)
           ! and floating conditions are satisfied;
           ! discretization of the x-component of the PDE

         st%flag_calc_vxy_ssa_x(j,i) = .true.

         st%flag_shelfy_stream_x(j,i) = .false.
                                   ! make sure not to treat as shelfy stream

         nc = 2*grd%ij2n(j,i-1)-1
                  ! smallest nc (column counter), for vx_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_wp*inv_dxi2*st%vis_int_g(j,i)
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j-1,i)-1
                  ! next nc (column counter), for vx_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_deta2*st%vis_int_sgxy(j-1,i)
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j-1,i)
                  ! next nc (column counter), for vy_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i)-1
                  ! next nc (column counter), for vx_m(j,i)
         if (nc /= nr) then   ! (diagonal element)
            stop ' >>> calc_vxy_ssa_matrix: Check for diagonal element failed!'
         end if
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -4.0_wp*inv_dxi2 &
                                 *(st%vis_int_g(j,i+1)+st%vis_int_g(j,i)) &
                          -inv_deta2 &
                                 *(st%vis_int_sgxy(j,i)+st%vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i)
                  ! next nc (column counter), for vy_m(j,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j+1,i)-1
                  ! next nc (column counter), for vx_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_deta2*st%vis_int_sgxy(j,i)
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j-1,i+1)
                  ! next nc (column counter), for vy_m(j-1,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j,i+1)+st%vis_int_sgxy(j-1,i))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i+1)-1
                  ! next nc (column counter), for vx_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_wp*inv_dxi2*st%vis_int_g(j,i+1)
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i+1)
                  ! largest nc (column counter), for vy_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j,i+1)+st%vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         lgs_b_value(nr) = factor_rhs_1*H_mid*st%dzs_dx_aux(j,i)

         lgs_x_value(nr) = st%vx_m_ssa(j,i)
 
      else if (st%flag_shelfy_stream_x(j,i)) then
           ! shelfy stream (as determined by routine calc_vxy_sia)

         st%flag_calc_vxy_ssa_x(j,i) = .true.

#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)

         if ( &
              ( ( st%flag_grounded_front_b_1(j,i) &
                       .and.st%flag_grounded_front_b_2(j,i+1) ) &
                .or. &
                ( st%flag_grounded_front_b_2(j,i) &
                       .and.st%flag_grounded_front_b_1(j,i+1) ) ) &
              .and. &
              ( zl_mid < st%z_sl(j,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#elif (BC_SSA_LTGF==2)

         if ( &
              ( st%flag_grounded_front_b_1(j,i) &
                     .and.st%flag_grounded_front_b_2(j,i+1) ) &
              .or. &
              ( st%flag_grounded_front_b_2(j,i) &
                     .and.st%flag_grounded_front_b_1(j,i+1) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#endif

            if (st%flag_grounded_front_b_1(j,i)) then
               i1 = i     ! grounded ice marker
            else   ! flag_grounded_front_b_1(j,i+1)==.true.
               i1 = i+1   ! grounded ice marker
            end if

            if (.not.( st%flag_grounded_front_b_2(j,i1-1) &
                       .and. &
                       st%flag_grounded_front_b_2(j,i1+1) ) ) then
               ! discretization of the x-component of the BC

               nc = 2*grd%ij2n(j,i1-1)-1
                        ! smallest nc (column counter), for vx_m(j,i1-1)
               k  = k+1
               lgs_a_value(k) = -4.0_wp*inv_dxi*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j-1,i1)
                        ! next nc (column counter), for vy_m(j-1,i1)
               k  = k+1
               lgs_a_value(k) = -2.0_wp*inv_deta*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j,i1)-1
                        ! next nc (column counter), for vx_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 4.0_wp*inv_dxi*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j,i1)
                        ! largest nc (column counter), for vy_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 2.0_wp*inv_deta*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(st%H_c(j,i1)+st%H_t(j,i1))**2 &
                               - factor_rhs_3b &
                                    *(max((st%z_sl(j,i)-st%zb(j,i1)), 0.0_wp))**2

               lgs_x_value(nr) = st%vx_m_ssa(j,i)

            else   !      (flag_grounded_front_b_2(j,i1-1)==.true.)
                   ! .and.(flag_grounded_front_b_2(j,i1+1)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_wp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_wp

               lgs_x_value(nr) = 0.0_wp

            end if

#if (BC_SSA_LTGF==2)

         else if ( &
              ( st%flag_grounded_front_a_1(j,i) &
                     .and.st%flag_grounded_front_a_2(j,i+1) ) &
              .or. &
              ( st%flag_grounded_front_a_2(j,i) &
                     .and.st%flag_grounded_front_a_1(j,i+1) ) &
            ) then
            ! one neighbour is grounded ice and the other is ice-free land
            ! (land-terminating grounded front)

            if (st%flag_grounded_front_a_1(j,i)) then
               i1 = i     ! grounded ice marker
            else   ! flag_grounded_front_a_1(j,i+1)==.true.
               i1 = i+1   ! grounded ice marker
            end if

            if (.not.( st%flag_grounded_front_a_2(j,i1-1) &
                       .and. &
                       st%flag_grounded_front_a_2(j,i1+1) ) ) then
               ! discretization of the x-component of the BC

               nc = 2*grd%ij2n(j,i1-1)-1
                        ! smallest nc (column counter), for vx_m(j,i1-1)
               k  = k+1
               lgs_a_value(k) = -4.0_wp*inv_dxi*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*st%ij2n(j-1,i1)
                        ! next nc (column counter), for vy_m(j-1,i1)
               k  = k+1
               lgs_a_value(k) = -2.0_wp*inv_deta*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j,i1)-1
                        ! next nc (column counter), for vx_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 4.0_wp*inv_dxi*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j,i1)
                        ! largest nc (column counter), for vy_m(j,i1)
               k  = k+1
               lgs_a_value(k) = 2.0_wp*inv_deta*st%vis_int_g(j,i1)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(st%H_c(j,i1)+st%H_t(j,i1))**2

               lgs_x_value(nr) = st%vx_m_ssa(j,i)

            else   !      (flag_grounded_front_a_2(j,i1-1)==.true.)
                   ! .and.(flag_grounded_front_a_2(j,i1+1)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_wp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_wp

               lgs_x_value(nr) = 0.0_wp

            end if

#endif

         else
            ! inner shelfy stream
            !   or
            ! one neighbour is floating ice and the other is grounded ice
            ! (grounding line)
            ! and floating conditions are not satisfied
#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)
            !   or
            ! land-terminating grounded front
#endif

            nc = 2*grd%ij2n(j,i-1)-1
                     ! smallest nc (column counter), for vx_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_wp*inv_dxi2*st%vis_int_g(j,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j-1,i)-1
                     ! next nc (column counter), for vx_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*st%vis_int_sgxy(j-1,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j-1,i)
                     ! next nc (column counter), for vy_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j-1,i))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i)-1
                     ! next nc (column counter), for vx_m(j,i)
            if (nc /= nr) then   ! (diagonal element)
               stop ' >>> calc_vxy_ssa_matrix: Check for diagonal element failed!'
            end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -4.0_wp*inv_dxi2 &
                                    *(st%vis_int_g(j,i+1)+st%vis_int_g(j,i)) &
                             -inv_deta2 &
                                    *(st%vis_int_sgxy(j,i)+st%vis_int_sgxy(j-1,i)) &
                             -0.5_wp*(st%beta_drag(j,i+1)+st%beta_drag(j,i))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i)
                     ! next nc (column counter), for vy_m(j,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j+1,i)-1
                     ! next nc (column counter), for vx_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_deta2*st%vis_int_sgxy(j,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j-1,i+1)
                     ! next nc (column counter), for vy_m(j-1,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_wp*st%vis_int_g(j,i+1)+st%vis_int_sgxy(j-1,i))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i+1)-1
                     ! next nc (column counter), for vx_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_wp*inv_dxi2*st%vis_int_g(j,i+1)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i+1)
                     ! largest nc (column counter), for vy_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_wp*st%vis_int_g(j,i+1)+st%vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_1*H_mid*st%dzs_dx_aux(j,i)

            lgs_x_value(nr) = st%vx_m_ssa(j,i)

         end if

      else if ( &
                ( st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j,i+1) &
                  .and.(H_mid >= (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
                .or. &
                ( st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j,i+1) &
                  .and.(H_mid >= (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
              ) then 
              ! one neighbour is floating ice and the other is grounded ice
              ! (grounding line)
              ! and floating conditions are not satisfied;
              ! velocity taken from the SIA solution for grounded ice

         st%flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = st%vx_m_sia(j,i)

         lgs_x_value(nr) = st%vx_m_sia(j,i)

      else if ( &
                ( (st%maske(j,i)==3).and.(st%maske(j,i+1)==1) ) &
                .or. &
                ( (st%maske(j,i)==1).and.(st%maske(j,i+1)==3) ) &
              ) then
              ! one neighbour is floating ice and the other is ice-free land;
              ! velocity assumed to be zero

         st%flag_calc_vxy_ssa_x(j,i) = .true.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_wp

         lgs_x_value(nr) = 0.0_wp

      else if ( &
                ( st%flag_calving_front_1(j,i).and.st%flag_calving_front_2(j,i+1) ) &
                .or. &
                ( st%flag_calving_front_2(j,i).and.st%flag_calving_front_1(j,i+1) ) &
              ) then
              ! one neighbour is floating ice and the other is ocean
              ! (calving front)

         st%flag_calc_vxy_ssa_x(j,i) = .true.

         if (st%flag_calving_front_1(j,i)) then
            i1 = i     ! floating ice marker
         else   ! flag_calving_front_1(j,i+1)==.true.
            i1 = i+1   ! floating ice marker
         end if

         if (.not.( st%flag_calving_front_2(j,i1-1) &
                    .and. &
                    st%flag_calving_front_2(j,i1+1) ) ) then
            ! discretization of the x-component of the BC

            nc = 2*grd%ij2n(j,i1-1)-1
                     ! smallest nc (column counter), for vx_m(j,i1-1)
            k  = k+1
            lgs_a_value(k) = -4.0_wp*inv_dxi*st%vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j-1,i1)
                     ! next nc (column counter), for vy_m(j-1,i1)
            k  = k+1
            lgs_a_value(k) = -2.0_wp*inv_deta*st%vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i1)-1
                     ! next nc (column counter), for vx_m(j,i1)
            k  = k+1
            lgs_a_value(k) = 4.0_wp*inv_dxi*st%vis_int_g(j,i1)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i1)
                     ! largest nc (column counter), for vy_m(j,i1)
            k  = k+1
            lgs_a_value(k) = 2.0_wp*inv_deta*st%vis_int_g(j,i1)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_2 &
                                 *(st%H_c(j,i1)+st%H_t(j,i1))**2

            lgs_x_value(nr) = st%vx_m_ssa(j,i)

         else   !      (flag_calving_front_2(j,i1-1)==.true.)
                ! .and.(flag_calving_front_2(j,i1+1)==.true.);
                ! velocity assumed to be zero

            k  = k+1
            lgs_a_value(k) = 1.0_wp   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0_wp

            lgs_x_value(nr) = 0.0_wp

         end if

      else if ( (st%maske(j,i)==0).or.(st%maske(j,i+1)==0) ) then
           ! neither neighbour is floating ice, but at least one neighbour is
           ! grounded ice; velocity taken from the SIA solution for grounded ice

         st%flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = st%vx_m_sia(j,i)

         lgs_x_value(nr) = st%vx_m_sia(j,i)

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         st%flag_calc_vxy_ssa_x(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_wp

         lgs_x_value(nr) = 0.0_wp

      end if

   else   ! boundary condition, velocity assumed to be zero

      st%flag_calc_vxy_ssa_x(j,i) = .false.

      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
      lgs_a_value(k) = 1.0_wp   ! diagonal element only
      lgs_a_index(k) = nr

      lgs_b_value(nr) = 0.0_wp

      lgs_x_value(nr) = 0.0_wp

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row

!  ------ Equations for vy_m_ssa (at (i,j+1/2))

   nr = n+1   ! row counter

   if ( (j /= grd%JMAX).and.(i /= 0).and.(i /= grd%IMAX) ) then
      ! inner point on the staggered grid in y-direction

      H_mid  = 0.5_wp*((st%H_c(j,i)+st%H_t(j,i))+(st%H_c(j+1,i)+st%H_t(j+1,i)))
      zl_mid = 0.5_wp*(st%zl(j,i)+st%zl(j+1,i))
   
      if ( &
           ( (st%maske(j,i)==3).and.(st%maske(j+1,i)==3) ) &
           .or. &
           ( st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j+1,i) &
             .and.(H_mid < (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
           .or. &
           ( st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j+1,i) &
             .and.(H_mid < (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
         ) then
           ! both neighbours are floating ice
           !   or
           ! one neighbour is floating ice and the other is grounded ice
           ! (grounding line)
           ! and floating conditions are satisfied;
           ! discretization of the y-component of the PDE

         st%flag_calc_vxy_ssa_y(j,i) = .true.

         st%flag_shelfy_stream_y(j,i) = .false.
                                   ! make sure not to treat as shelfy stream

         nc = 2*grd%ij2n(j,i-1)-1
                  ! smallest nc (column counter), for vx_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i-1)
                  ! next nc (column counter), for vy_m(j,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi2*st%vis_int_sgxy(j,i-1)
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j+1,i-1)-1
                  ! next nc (column counter), for vx_m(j+1,i-1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j+1,i)+st%vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j-1,i)
                  ! next nc (column counter), for vy_m(j-1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_wp*inv_deta2*st%vis_int_g(j,i)
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i)-1
                  ! next nc (column counter), for vx_m(j,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i)
                  ! next nc (column counter), for vy_m(j,i)
         if (nc /= nr) then   ! (diagonal element)
            stop ' >>> calc_vxy_ssa_matrix: Check for diagonal element failed!'
         end if
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = -4.0_wp*inv_deta2 &
                                 *(st%vis_int_g(j+1,i)+st%vis_int_g(j,i)) &
                          -inv_dxi2 &
                                 *(st%vis_int_sgxy(j,i)+st%vis_int_sgxy(j,i-1))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j+1,i)-1
                  ! next nc (column counter), for vx_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi_deta &
                                 *(2.0_wp*st%vis_int_g(j+1,i)+st%vis_int_sgxy(j,i))
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j+1,i)
                  ! next nc (column counter), for vy_m(j+1,i)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 4.0_wp*inv_deta2*st%vis_int_g(j+1,i)
         lgs_a_index(k) = nc

         nc = 2*grd%ij2n(j,i+1)
                  ! largest nc (column counter), for vy_m(j,i+1)
         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = inv_dxi2*st%vis_int_sgxy(j,i)
         lgs_a_index(k) = nc

         lgs_b_value(nr) = factor_rhs_1*H_mid*st%dzs_dy_aux(j,i)

         lgs_x_value(nr) = st%vy_m_ssa(j,i)
 
      else if (st%flag_shelfy_stream_y(j,i)) then
           ! shelfy stream (as determined by routine calc_vxy_sia)

         st%flag_calc_vxy_ssa_y(j,i) = .true.

#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)

         if ( &
              ( ( st%flag_grounded_front_b_1(j,i) &
                       .and.st%flag_grounded_front_b_2(j+1,i) ) &
                .or. &
                ( st%flag_grounded_front_b_2(j,i) &
                       .and.st%flag_grounded_front_b_1(j+1,i) ) ) &
              .and. &
              ( zl_mid < st%z_sl(j,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#elif (BC_SSA_LTGF==2)

         if ( &
              ( st%flag_grounded_front_b_1(j,i) &
                     .and.st%flag_grounded_front_b_2(j+1,i) ) &
              .or. &
              ( st%flag_grounded_front_b_2(j,i) &
                     .and.st%flag_grounded_front_b_1(j+1,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ocean
            ! (ocean-terminating grounded front)

#endif

            if (st%flag_grounded_front_b_1(j,i)) then
               j1 = j     ! grounded ice marker
            else   ! flag_grounded_front_b_1(j+1,i)==.true.
               j1 = j+1   ! grounded ice marker
            end if

            if (.not.( st%flag_grounded_front_b_2(j1-1,i) &
                       .and. &
                       st%flag_grounded_front_b_2(j1+1,i) ) ) then
               ! discretization of the y-component of the BC

               nc = 2*grd%ij2n(j1,i-1)-1
                        ! smallest nc (column counter), for vx_m(j1,i-1)
               k  = k+1
               lgs_a_value(k) = -2.0_wp*inv_dxi*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j1-1,i)
                        ! next nc (column counter), for vy_m(j1-1,i)
               k  = k+1
               lgs_a_value(k) = -4.0_wp*inv_deta*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j1,i)-1
                        ! next nc (column counter), for vx_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 2.0_wp*inv_dxi*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j1,i)
                        ! largest nc (column counter), for vy_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 4.0_wp*inv_deta*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(st%H_c(j1,i)+st%H_t(j1,i))**2 &
                               - factor_rhs_3b &
                                    *(max((st%z_sl(j,i)-st%zb(j1,i)), 0.0_wp))**2

               lgs_x_value(nr) = st%vy_m_ssa(j,i)

            else   !      (flag_grounded_front_b_2(j1-1,i)==.true.)
                   ! .and.(flag_grounded_front_b_2(j1+1,i)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_wp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_wp

               lgs_x_value(nr) = 0.0_wp

            end if

#if (BC_SSA_LTGF==2)

         else if ( &
              ( st%flag_grounded_front_a_1(j,i) &
                     .and.st%flag_grounded_front_a_2(j+1,i) ) &
              .or. &
              ( st%flag_grounded_front_a_2(j,i) &
                     .and.st%flag_grounded_front_a_1(j+1,i) ) &
            ) then
            ! one neighbour is grounded ice and the other is ice-free land
            ! (land-terminating grounded front)

            if (st%flag_grounded_front_a_1(j,i)) then
               j1 = j     ! grounded ice marker
            else   ! flag_grounded_front_a_1(j+1,i)==.true.
               j1 = j+1   ! grounded ice marker
            end if

            if (.not.( st%flag_grounded_front_a_2(j1-1,i) &
                       .and. &
                       st%flag_grounded_front_a_2(j1+1,i) ) ) then
               ! discretization of the y-component of the BC

               nc = 2*grd%ij2n(j1,i-1)-1
                        ! smallest nc (column counter), for vx_m(j1,i-1)
               k  = k+1
               lgs_a_value(k) = -2.0_wp*inv_dxi*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j1-1,i)
                        ! next nc (column counter), for vy_m(j1-1,i)
               k  = k+1
               lgs_a_value(k) = -4.0_wp*inv_deta*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j1,i)-1
                        ! next nc (column counter), for vx_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 2.0_wp*inv_dxi*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               nc = 2*grd%ij2n(j1,i)
                        ! largest nc (column counter), for vy_m(j1,i)
               k  = k+1
               lgs_a_value(k) = 4.0_wp*inv_deta*st%vis_int_g(j1,i)
               lgs_a_index(k) = nc

               lgs_b_value(nr) = factor_rhs_3a &
                                    *(st%H_c(j1,i)+st%H_t(j1,i))**2

               lgs_x_value(nr) = st%vy_m_ssa(j,i)

            else   !      (flag_grounded_front_a_2(j1-1,i)==.true.)
                   ! .and.(flag_grounded_front_a_2(j1+1,i)==.true.);
                   ! velocity assumed to be zero

               k  = k+1
               lgs_a_value(k) = 1.0_wp   ! diagonal element only
               lgs_a_index(k) = nr

               lgs_b_value(nr) = 0.0_wp

               lgs_x_value(nr) = 0.0_wp

            end if

#endif

         else
            ! inner shelfy stream
            !   or
            ! one neighbour is floating ice and the other is grounded ice
            ! (grounding line)
            ! and floating conditions are not satisfied
#if (!defined(BC_SSA_LTGF) || BC_SSA_LTGF==1)
            !   or
            ! land-terminating grounded front
#endif

            nc = 2*grd%ij2n(j,i-1)-1
                     ! smallest nc (column counter), for vx_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j,i-1))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i-1)
                     ! next nc (column counter), for vy_m(j,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*st%vis_int_sgxy(j,i-1)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j+1,i-1)-1
                     ! next nc (column counter), for vx_m(j+1,i-1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                  *(2.0_wp*st%vis_int_g(j+1,i)+st%vis_int_sgxy(j,i-1))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j-1,i)
                     ! next nc (column counter), for vy_m(j-1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_wp*inv_deta2*st%vis_int_g(j,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i)-1
                     ! next nc (column counter), for vx_m(j,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -inv_dxi_deta &
                                    *(2.0_wp*st%vis_int_g(j,i)+st%vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i)
                     ! next nc (column counter), for vy_m(j,i)
            if (nc /= nr) then   ! (diagonal element)
               stop ' >>> calc_vxy_ssa_matrix: Check for diagonal element failed!'
            end if
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = -4.0_wp*inv_deta2 &
                                    *(st%vis_int_g(j+1,i)+st%vis_int_g(j,i)) &
                             -inv_dxi2 &
                                    *(st%vis_int_sgxy(j,i)+st%vis_int_sgxy(j,i-1)) &
                             -0.5_wp*(st%beta_drag(j+1,i)+st%beta_drag(j,i))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j+1,i)-1
                     ! next nc (column counter), for vx_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi_deta &
                                    *(2.0_wp*st%vis_int_g(j+1,i)+st%vis_int_sgxy(j,i))
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j+1,i)
                     ! next nc (column counter), for vy_m(j+1,i)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = 4.0_wp*inv_deta2*st%vis_int_g(j+1,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j,i+1)
                     ! largest nc (column counter), for vy_m(j,i+1)
            k  = k+1
            ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
            lgs_a_value(k) = inv_dxi2*st%vis_int_sgxy(j,i)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_1*H_mid*st%dzs_dy_aux(j,i)

            lgs_x_value(nr) = st%vy_m_ssa(j,i)

         end if

      else if ( &
                ( st%flag_grounding_line_1(j,i).and.st%flag_grounding_line_2(j+1,i) &
                  .and.(H_mid >= (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
                .or. &
                ( st%flag_grounding_line_2(j,i).and.st%flag_grounding_line_1(j+1,i) &
                  .and.(H_mid >= (st%z_sl(j,i)-zl_mid)*rhosw_rho_ratio) ) &
              ) then
              ! one neighbour is floating ice and the other is grounded ice
              ! (grounding line)
              ! and floating conditions are not satisfied;
              ! velocity taken from the SIA solution for grounded ice

         st%flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = st%vy_m_sia(j,i)

         lgs_x_value(nr) = st%vy_m_sia(j,i)

      else if ( &
                ( (st%maske(j,i)==3).and.(st%maske(j+1,i)==1) ) &
                .or. &
                ( (st%maske(j,i)==1).and.(st%maske(j+1,i)==3) ) &
              ) then
           ! one neighbour is floating ice and the other is ice-free land;
           ! velocity assumed to be zero

         st%flag_calc_vxy_ssa_y(j,i) = .true.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_wp

         lgs_x_value(nr) = 0.0_wp

      else if ( &
                ( st%flag_calving_front_1(j,i).and.st%flag_calving_front_2(j+1,i) ) &
                .or. &
                ( st%flag_calving_front_2(j,i).and.st%flag_calving_front_1(j+1,i) ) &
              ) then
              ! one neighbour is floating ice and the other is ocean
              ! (calving front)

         st%flag_calc_vxy_ssa_y(j,i) = .true.

         if (st%flag_calving_front_1(j,i)) then
            j1 = j     ! floating ice marker
         else   ! flag_calving_front_1(j+1,i)==.true.
            j1 = j+1   ! floating ice marker
         end if

         if (.not.( st%flag_calving_front_2(j1-1,i) &
                    .and. &
                    st%flag_calving_front_2(j1+1,i) ) ) then
            ! discretization of the y-component of the BC

            nc = 2*grd%ij2n(j1,i-1)-1
                     ! smallest nc (column counter), for vx_m(j1,i-1)
            k  = k+1
            lgs_a_value(k) = -2.0_wp*inv_dxi*st%vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j1-1,i)
                     ! next nc (column counter), for vy_m(j1-1,i)
            k  = k+1
            lgs_a_value(k) = -4.0_wp*inv_deta*st%vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j1,i)-1
                     ! next nc (column counter), for vx_m(j1,i)
            k  = k+1
            lgs_a_value(k) = 2.0_wp*inv_dxi*st%vis_int_g(j1,i)
            lgs_a_index(k) = nc

            nc = 2*grd%ij2n(j1,i)
                     ! largest nc (column counter), for vy_m(j1,i)
            k  = k+1
            lgs_a_value(k) = 4.0_wp*inv_deta*st%vis_int_g(j1,i)
            lgs_a_index(k) = nc

            lgs_b_value(nr) = factor_rhs_2 &
                                 *(st%H_c(j1,i)+st%H_t(j1,i))**2

            lgs_x_value(nr) = st%vy_m_ssa(j,i)

         else   !      (flag_calving_front_2(j1-1,i)==.true.)
                ! .and.(flag_calving_front_2(j1+1,i)==.true.);
                ! velocity assumed to be zero

            k  = k+1
            lgs_a_value(k) = 1.0_wp   ! diagonal element only
            lgs_a_index(k) = nr

            lgs_b_value(nr) = 0.0_wp

            lgs_x_value(nr) = 0.0_wp

         end if

      else if ( (st%maske(j,i)==0).or.(st%maske(j+1,i)==0) ) then
           ! neither neighbour is floating ice, but at least one neighbour is
           ! grounded ice; velocity taken from the SIA solution for grounded ice

         st%flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = st%vy_m_sia(j,i)

         lgs_x_value(nr) = st%vy_m_sia(j,i)

      else   ! neither neighbour is floating or grounded ice,
             ! velocity assumed to be zero

         st%flag_calc_vxy_ssa_y(j,i) = .false.

         k  = k+1
         ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
         lgs_a_value(k) = 1.0_wp   ! diagonal element only
         lgs_a_index(k) = nr

         lgs_b_value(nr) = 0.0_wp

         lgs_x_value(nr) = 0.0_wp

      end if

   else   ! boundary condition, velocity assumed to be zero

      st%flag_calc_vxy_ssa_y(j,i) = .false.

      k  = k+1
      ! if (k > n_sprs) stop ' >>> calc_vxy_ssa_matrix: n_sprs too small!'
      lgs_a_value(k) = 1.0_wp   ! diagonal element only
      lgs_a_index(k) = nr

      lgs_b_value(nr) = 0.0_wp

      lgs_x_value(nr) = 0.0_wp

   end if

   lgs_a_ptr(nr+1) = k+1   ! row is completed, store index to next row
    
end do
!!!$omp end parallel do
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'fill_vectors',time2-time1

!-------- Settings for Lis --------

call lis_initialize(ierr)
call CHKERR(ierr)

call lis_matrix_create(LIS_COMM_WORLD, lgs_a, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_b, ierr)
call lis_vector_create(LIS_COMM_WORLD, lgs_x, ierr)

call lis_matrix_set_size(lgs_a, 0, nmax, ierr)
call lis_vector_set_size(lgs_b, 0, nmax, ierr)
call lis_vector_set_size(lgs_x, 0, nmax, ierr)

!$ time1 = omp_get_wtime()
do nr=1, nmax

   do nc=lgs_a_ptr(nr), lgs_a_ptr(nr+1)-1
      call lis_matrix_set_value(LIS_INS_VALUE, nr, lgs_a_index(nc), &
                                               lgs_a_value(nc), lgs_a, ierr)
   end do

   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_b_value(nr), lgs_b, ierr)
   call lis_vector_set_value(LIS_INS_VALUE, nr, lgs_x_value(nr), lgs_x, ierr)

end do

call lis_matrix_set_type(lgs_a, LIS_MATRIX_CSR, ierr)

!1: LIS_INT i,k,n,nnz;
!2: LIS_INT *ptr,*index;
!3: LIS_SCALAR *value;
!4: LIS_MATRIX A;
!5: n = 4; nnz = 10; k = 0;
!6: lis_matrix_malloc_csr(n,nnz,&ptr,&index,&value);
!7: lis_matrix_create(0,&A);
!8: lis_matrix_set_size(A,0,n); /* or lis_matrix_set_size(A,n,0); */
!9:
!10: for(i=0;i<n;i++)
!11: {
!12: if( i>0 ) {index[k] = i-1; value[k] = 1; k++;}
!13: index[k] = i; value[k] = 2; k++;
!14: if( i<n-1 ) {index[k] = i+1; value[k] = 1; k++;}
!15: ptr[i+1] = k;
!16: }
!17: ptr[0] = 0;
!18: lis_matrix_set_csr(nnz,ptr,index,value,A);
!19: lis_matrix_assemble(A);
!lis_matrix_set_csr(LIS_INTEGER nnz, LIS_INTEGER ptr(), LIS_INTEGER index(), LIS_SCALAR value(), LIS_MATRIX A, LIS_INTEGER ierr)

! todo, could this be faster? not working yet
!call lis_matrix_malloc_csr(lgs_a, nmax, k, ierr)
!print *,ierr
!!!call lis_matrix_set_csr(k, lgs_a_ptr-1, lgs_a_index-1, lgs_a_value, lgs_a, ierr)
!call lis_matrix_set_csr(k, lgs_a_ptr, lgs_a_index, lgs_a_value, lgs_a, ierr)
!print *,ierr
!call lis_matrix_assemble(lgs_a, ierr)
!call lis_vector_set_values(LIS_INS_VALUE, nmax, [1:nmax], lgs_b_value(:), lgs_b, ierr)
!call lis_vector_set_values(LIS_INS_VALUE, nmax, [1:nmax], lgs_x_value(:), lgs_x, ierr)
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'lis_solver_Setup',time2-time1

!$ time1 = omp_get_wtime()
call lis_matrix_assemble(lgs_a, ierr)
!print *,ierr
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'lis_matrix_assemble',time2-time1


!-------- Solution of the system of linear equations with Lis --------

call lis_solver_create(solver, ierr)

    ch_solver_set_option = trim(par%lis_opts)

call lis_solver_set_option(trim(ch_solver_set_option), solver, ierr)
call CHKERR(ierr)

    !$ time1 = omp_get_wtime()
call lis_solve(lgs_a, lgs_b, lgs_x, solver, ierr)
call CHKERR(ierr)
    !$ time2 = omp_get_wtime()
    !$ if(par%l_write_timer) print *,'lis_solver',time2-time1

call lis_solver_get_iter(solver, lin_iter, ierr)

!write(6,'(a,i0,a)', advance='no') 'lin_iter = ', lin_iter, ', '

!call lis_solver_get_time(solver,solver_time,ierr)
!print *, 'calc_vxy_ssa_matrix: time (s) = ', solver_time

lgs_x_value = 0.0_wp
call lis_vector_gather(lgs_x, lgs_x_value, ierr)
call lis_matrix_destroy(lgs_a, ierr)
call lis_vector_destroy(lgs_b, ierr)
call lis_vector_destroy(lgs_x, ierr)
call lis_solver_destroy(solver, ierr)

call lis_finalize(ierr)      
call CHKERR(ierr)

do n=1, nmax-1, 2

   i = grd%n2i((n+1)/2)
   j = grd%n2j((n+1)/2)

   nr = n
   st%vx_m_ssa(j,i) = lgs_x_value(nr)

   nr = n+1
   st%vy_m_ssa(j,i) = lgs_x_value(nr)

end do

deallocate(lgs_a_value, lgs_a_index, lgs_a_ptr)
deallocate(lgs_b_value, lgs_x_value)

end subroutine calc_vxy_ssa_matrix


!-------------------------------------------------------------------------------
!> Computation of the depth-integrated viscosity vis_int_g in the
!! shallow shelf approximation.
!<------------------------------------------------------------------------------
subroutine calc_vis_ssa(st,grd,par)

!$ use omp_lib

implicit none

type(sico_state_class), intent(inout) :: st
type(sico_grid_class), intent(in) :: grd
type(sico_par_class), intent(in) :: par

integer :: i, j, kc, kt, m
integer :: m_smooth_abs
real(wp) :: visc_min, visc_max, visc_init
real(wp) :: dvx_dxi, dvx_deta, dvy_dxi, dvy_deta
real(wp) :: cvis0(0:100), cvis1(0:100)


visc_min  = par%visc_min 
visc_max  = par%visc_max 
visc_init = par%visc_init_ssa

!-------- Computation of the depth-integrated viscosity --------

!!$omp parallel do private(i,j,kc,kt,dvx_dxi, dvx_deta, dvy_dxi, dvy_deta,cvis0,cvis1)
do i=0, grd%IMAX
do j=0, grd%JMAX

   if ((st%maske(j,i)==0).and.(.not.st%flag_shelfy_stream(j,i))) then
                                                   ! grounded ice, but
                                                   ! not shelfy stream
      st%de_ssa(j,i) = 0.0_wp   ! dummy value

      st%vis_ave_g(j,i) = 1.0_wp/st%flui_ave_sia(j,i)
      st%vis_int_g(j,i) = (st%H_c(j,i)+st%H_t(j,i)) * st%vis_ave_g(j,i) 

   else if ((st%maske(j,i)==1).or.(st%maske(j,i)==2)) then
                                                   ! ice-free land or ocean
      st%de_ssa(j,i) = 0.0_wp   ! dummy value

      st%vis_ave_g(j,i) = visc_init   ! dummy value
      st%vis_int_g(j,i) = 0.0_wp   ! dummy value

   else   ! (maske(j,i)==3).or.(st%flag_shelfy_stream(j,i)),
          ! floating ice or shelfy stream; 
          ! must not be at the margin of the computational domain

!  ------ Effective strain rate

      dvx_dxi  = (st%vx_m_ssa(j,i)-st%vx_m_ssa(j,i-1))*grd%dxi_inv 
      dvy_deta = (st%vy_m_ssa(j,i)-st%vy_m_ssa(j-1,i))*grd%deta_inv

      dvx_deta = 0.25_wp*grd%deta_inv &
                 *(st%vx_m_ssa(j+1,i)+st%vx_m_ssa(j+1,i-1)-st%vx_m_ssa(j-1,i)-st%vx_m_ssa(j-1,i-1))
      dvy_dxi  = 0.25_wp*grd%dxi_inv &
                 *(st%vy_m_ssa(j,i+1)+st%vy_m_ssa(j-1,i+1)-st%vy_m_ssa(j,i-1)-st%vy_m_ssa(j-1,i-1))

      st%de_ssa(j,i) = sqrt( dvx_dxi*dvx_dxi &
                          + dvy_deta*dvy_deta &
                          + dvx_dxi*dvy_deta &
                          + 0.25_wp*(dvx_deta+dvy_dxi)*(dvx_deta+dvy_dxi) )

!  ------ Term abbreviations

      if (par%dynamics==2) then

        if (.not.st%flag_shelfy_stream(j,i)) then

          if (par%calcmod==1) then
            cvis0 = 0._wp ! CMW fix, needed to avoid uninitialized variable below
          endif

          do kc=0, grd%KCMAX
            cvis1(kc) = grd%aqxy1(kc)*st%H_c(j,i) &
              *viscosity(par%fin_visc,par%flow_law,st%de_ssa(j,i), &
              st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
              st%enh_c(kc,j,i), 0)
          end do
          ! Ice shelves (floating ice) are assumed to consist of cold ice only

        else   ! flag_shelfy_stream(j,i) == .true.

          if (par%calcmod==-1 .or. par%calcmod==0) then

            do kc=0, grd%KCMAX
              cvis1(kc) = grd%aqxy1(kc)*st%H_c(j,i) &
                *viscosity(par%fin_visc,par%flow_law,st%de_ssa(j,i), &
                st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
                st%enh_c(kc,j,i), 0)
            end do

          else if (par%calcmod==1) then

            do kt=0, grd%KTMAX
              cvis0(kt) = grd%dzeta_t*st%H_t(j,i) &
                *viscosity(par%fin_visc,par%flow_law,st%de_ssa(j,i), &
                st%temp_t_m(kt,j,i), st%temp_t_m(kt,j,i), st%omega_t(kt,j,i), &
                st%enh_t(kt,j,i), 1)
            end do

            do kc=0, grd%KCMAX
              cvis1(kc) = grd%aqxy1(kc)*st%H_c(j,i) &
                *viscosity(par%fin_visc,par%flow_law,st%de_ssa(j,i), &
                st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
                st%enh_c(kc,j,i), 0)
            end do

          else if (par%calcmod==2 .or. par%calcmod==3) then

            do kc=0, grd%KCMAX
              cvis1(kc) = grd%aqxy1(kc)*st%H_c(j,i) &
                *viscosity(par%fin_visc,par%flow_law,st%de_ssa(j,i), &
                st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), st%omega_c(kc,j,i), &
                st%enh_c(kc,j,i), 2)
            end do

          else
            stop ' >>> calc_vis_ssa: CALCMOD must be -1, 0, 1, 2 or 3!'
          endif

        end if

      else  ! par%dynamics.ne.2

        do kc=0, grd%KCMAX
          cvis1(kc) = grd%aqxy1(kc)*st%H_c(j,i) &
            *viscosity(par%fin_visc,par%flow_law,st%de_ssa(j,i), &
            st%temp_c(kc,j,i), st%temp_c_m(kc,j,i), 0.0_wp, &
            st%enh_c(kc,j,i), 0)
        end do
        ! Ice shelves (floating ice) are assumed to consist of cold ice only

      endif

!  ------ Depth-integrated viscosity

      st%vis_int_g(j,i) = 0.0_wp

      if (par%calcmod==1) then
        do kt=0, grd%KTMAX-1
           st%vis_int_g(j,i) = st%vis_int_g(j,i)+0.5_wp*(cvis0(kt+1)+cvis0(kt))
        end do
      endif

      do kc=0, grd%KCMAX-1
         st%vis_int_g(j,i) = st%vis_int_g(j,i)+0.5_wp*(cvis1(kc+1)+cvis1(kc))
      end do

!  ------ Depth-averaged viscosity

      st%vis_ave_g(j,i) = st%vis_int_g(j,i)/max((st%H_c(j,i)+st%H_t(j,i)), eps_wp)

      st%vis_ave_g(j,i) = max(min(st%vis_ave_g(j,i), visc_max), visc_min)

   end if
  
end do
end do
!!$omp end parallel do

!-------- Smoothing of the depth-averaged viscosity --------

if (par%n_visc_smooth /= 0) then

   m_smooth_abs = abs(par%n_visc_smooth)

   if (par%n_visc_smooth < 0) st%vis_ave_g = log(st%vis_ave_g)   ! logarithmic smoothing

   st%vis_ave_g_smooth = st%vis_ave_g

   do m=1, m_smooth_abs

      do i=1, grd%IMAX-1
      do j=1, grd%JMAX-1
         st%vis_ave_g_smooth(j,i) = (1.0_wp-4.0_wp*par%visc_smooth_diff)*st%vis_ave_g(j,i) &
                                    + par%visc_smooth_diff &
                                       *( (st%vis_ave_g(j,i+1)+st%vis_ave_g(j,i-1)) &
                                         +(st%vis_ave_g(j+1,i)+st%vis_ave_g(j-1,i)) )
      end do
      end do

      st%vis_ave_g = st%vis_ave_g_smooth

   end do

   if (par%n_visc_smooth < 0) st%vis_ave_g = exp(st%vis_ave_g)   ! logarithmic smoothing

end if

!-------- Final depth-integrated viscosity --------

st%vis_int_g = st%vis_ave_g*(st%H_c+st%H_t)


end subroutine calc_vis_ssa

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!> Gradual limitation of computed horizontal velocities to the interval
!! [-vel_max, vel_max].
!<------------------------------------------------------------------------------
elemental subroutine velocity_limiter_gradual(velocity, vel_max, vel_max_inv)

implicit none

real(wp), intent(in)    :: vel_max, vel_max_inv
real(wp), intent(inout) :: velocity

real(wp) :: vel_abs, vel_sign, vel_scaled, vel_scaled_lim

vel_abs = abs(velocity)

if (vel_abs >= 1.1_wp*vel_max) then

   vel_sign = sign(1.0_wp, velocity)
   velocity = vel_sign * vel_max

else if (vel_abs > 0.9_wp*vel_max) then

   ! gradual limitation between 0.9*vel_max and 1.1*vel_max

   vel_sign = sign(1.0_wp, velocity)

   vel_scaled     = (vel_abs-0.9_wp*vel_max)*(10.0_wp*vel_max_inv)
                       ! between 0 and 2
   vel_scaled_lim = vel_scaled &
                       *(1.0_wp-0.25_wp*vel_scaled*vel_scaled &
                                       *(1.0_wp-0.25_wp*vel_scaled))
                       ! between 0 and 1

   velocity = vel_sign * vel_max * (0.9_wp + 0.1_wp*vel_scaled_lim)

end if

end subroutine velocity_limiter_gradual

!-------------------------------------------------------------------------------

end module calc_vxy_m
!

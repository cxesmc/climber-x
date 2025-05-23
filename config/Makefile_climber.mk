# Driver related source files and directories

########################################################################
# gitversion recording
dir_version = $(srcdir)/main/
files_version = gitversion.f90
tmp_version = $(patsubst %.f90, %.o, $(files_version) )
obj_version = $(patsubst %, $(objdir)/%, $(tmp_version) )
########################################################################

########################################################################
# atmospheric model related source files
dir_atm = $(srcdir)/atm/
files_atm = atm_params.f90 atm_def.f90 atm_grid.f90 smooth_atm.f90 \
						adifa.f90 synop.f90 wvel.f90 clouds.f90 vesta.f90 crisa.f90 slp.f90 u2d.f90 u3d.f90 \
						time_step.f90 lwr.f90 swr.f90 feedbacks.f90 rad_kernels.f90 dust.f90 \
						atm_model.f90 atm_out.f90
tmp_atm = $(patsubst %.f90, %.o, $(files_atm) )
obj_atm = $(patsubst %, $(objdir)/%, $(tmp_atm) )
########################################################################

########################################################################
# ocean biogeochemistry model related source files
dir_bgc = $(srcdir)/bgc/
files_bgc = bgc_grid.f90 bgc_params.f90 bgc_def.f90 mo_biomod.f90 \
	          bgc_diag.f90 apply_net_fw.f90 chemcon.f90 \
                  mo_m4ago_types.f90 mo_m4ago_params.f90 \
                  mo_m4ago_physics.f90 mo_m4ago_core.f90 mo_m4ago_HAMOCCinit.f90 \
                  mo_m4ago_HAMOCCPrimPart.f90 mo_m4ago_climberx.f90 \
	          ocprod.f90 carchm.f90 cyano.f90 \
	          powadi.f90 dipowa.f90 powach.f90 sedshi.f90 \
	          mo_beleg_bgc.f90 tracer_cons.f90 flx_sed.f90 \
                  flx_sed_net.f90 flx_bur.f90 sed_grid_update_tracers.f90 \
                  bgc_model.f90 bgc_out.f90
tmp_bgc = $(patsubst %.f90, %.o, $(files_bgc) )
obj_bgc = $(patsubst %, $(objdir)/%, $(tmp_bgc))

# addition of M4AGO sinking scheme (git submodule in separate directory)
dir_m4ago = $(srcdir)/bgc/m4ago/src/
files_m4ago = mo_m4ago_types.f90 mo_m4ago_params.f90 mo_m4ago_physics.f90 mo_m4ago_core.f90 \
              mo_m4ago_HAMOCCinit.f90 mo_m4ago_HAMOCCPrimPart.f90 mo_m4ago_climberx.f90	      
tmp_m4ago = $(patsubst %.f90, %.o, $(files_m4ago) )
obj_m4ago = $(patsubst %, $(objdir)/%, $(tmp_m4ago))

# dummy for climate only setup
dir_bgc_dummy = $(srcdir)/bgc-dummy/
files_bgc_dummy = .bgc_model.f90 .bgc_out.f90 .bgc_params.f90 .bgc_def.f90
tmp_bgc_dummy = $(patsubst %.f90, %.o, $(files_bgc_dummy) )
obj_bgc_dummy = $(patsubst %, $(objdir)/%, $(tmp_bgc_dummy) )
########################################################################

########################################################################
# boundary source files
dir_bnd = $(srcdir)/bnd/
files_bnd = bnd.f90 insolation.f90 solar.f90 volc.f90 co2.f90 co2_rad.f90 ch4.f90 ch4_rad.f90 n2o.f90 so4.f90 o3.f90 cfc.f90 \
						luc.f90 dist.f90 sea_level.f90 d13c_atm.f90 D14c_atm.f90 \
						fake_atm.f90 fake_dust.f90 fake_lnd.f90 fake_ocn.f90 fake_sic.f90 fake_ice.f90 fake_geo.f90
tmp_bnd = $(patsubst %.f90, %.o, $(files_bnd) )
obj_bnd = $(patsubst %, $(objdir)/%, $(tmp_bnd) )
########################################################################

########################################################################
# atmospheric CH4 model related source files
dir_ch4 = $(srcdir)/ch4/
files_ch4 = ch4_model.f90 ch4_def.f90 ch4_out.f90
tmp_ch4 = $(patsubst %.f90, %.o, $(files_ch4) )
obj_ch4 = $(patsubst %, $(objdir)/%, $(tmp_ch4) )
########################################################################

########################################################################
# atmospheric CO2 model related source files
dir_co2 = $(srcdir)/co2/
files_co2 = co2_model.f90 co2_def.f90 co2_out.f90
tmp_co2 = $(patsubst %.f90, %.o, $(files_co2) )
obj_co2 = $(patsubst %, $(objdir)/%, $(tmp_co2) )
########################################################################

########################################################################
# geo source files
dir_geo = $(srcdir)/geo/
files_geo = geo_params.f90 geo_grid.f90 geo_def.f90 runoff_routing.f90 lakes.f90 fill_ocean.f90 topo_filter.f90 topo_fill.f90 \
						connect_ocn.f90 vilma.f90 gia.f90 q_geo.f90 sed.f90 \
						corals_topo.f90 fix_runoff.f90 hires_to_lowres.f90 coast_cells.f90 drainage_basins.f90 geo.f90 geo_out.f90
tmp_geo = $(patsubst %.f90, %.o, $(files_geo) )
obj_geo = $(patsubst %, $(objdir)/%, $(tmp_geo) )
########################################################################

########################################################################
# ice related source files (ice)
dir_ice = $(srcdir)/ice/
files_ice = ice_model.f90 ice_def.f90 ice_id.f90
tmp_ice = $(patsubst %.f90, %.o, $(files_ice) )
obj_ice = $(patsubst %, $(objdir)/%, $(tmp_ice))
# dummy
files_ice_dummy = .ice_model_dummy.f90 ice_def.f90 ice_id.f90
tmp_ice_dummy = $(patsubst %.f90, %.o, $(files_ice_dummy) )
obj_ice_dummy = $(patsubst %, $(objdir)/%, $(tmp_ice_dummy))
########################################################################

########################################################################
# ice related source files (SICOPOLIS)
dir_ice_sico = $(srcdir)/ice_sico/
files_ice_sico = sico_types_m.f90 compare_float_m.f90 sico_params.f90 ice_material_properties_m.f90 stereo_proj_m.f90 \
						     sico_timer.f90 sico_grid.f90 sico_state.f90 sico_check.f90 \
                 sico_maths_m.f90 init_temp_water_age_m.f90 \
                 calc_vxy_m.f90 calc_vz_m.f90 calc_dxyz_m.f90 topograd_m.f90 \
                 calc_thk_m.f90 enth_temp_omega_m.f90 calc_temp_m.f90 calc_temp_enth_m.f90 \
                 calc_enhance_m.f90 flag_update_gf_gl_cf_m.f90\
                 calc_temp_melt_bas_m.f90 calc_bas_melt_m.f90 calc_thk_water_bas_m.f90 hydro_m.f90 \
								 discharge_workers_m.f90 boundary_m.f90 \
                 sico_model.f90 sico_out.f90
tmp_ice_sico = $(patsubst %.f90, %.o, $(files_ice_sico) )
obj_ice_sico = $(patsubst %, $(objdir)/%, $(tmp_ice_sico))
########################################################################

########################################################################
# bmb related source files
dir_bmb = $(srcdir)/bmb/
files_bmb = bmb_model.f90 bmb_params.f90 bmb_grid.f90 bmb_def.f90 bmb_out.f90 
tmp_bmb = $(patsubst %.f90, %.o, $(files_bmb) )
obj_bmb = $(patsubst %, $(objdir)/%, $(tmp_bmb) )
# dummy
files_bmb_dummy = .bmb_model_dummy.f90 bmb_def.f90 .bmb_out_dummy.f90 
tmp_bmb_dummy = $(patsubst %.f90, %.o, $(files_bmb_dummy) )
obj_bmb_dummy = $(patsubst %, $(objdir)/%, $(tmp_bmb_dummy) )
########################################################################

########################################################################
# land model related source files
dir_lnd = $(srcdir)/lnd/
files_lnd = lnd_model.f90 lnd_grid.f90 lnd_params.f90 lnd_def.f90 \
	          soil_par.f90 ice_par.f90 lake_par.f90 lake_rho.f90 shelf_par.f90 surface_par_lnd.f90 veg_par.f90 \
	          soil_carbon_par.f90 shelf_carbon_par.f90 ice_carbon_par.f90 lake_carbon_par.f90 \
            photosynthesis.f90 ebal_veg.f90 ebal_ice.f90 ebal_lake.f90 surface_hydro.f90 \
            soil_temp.f90 ice_temp.f90 shelf_temp.f90 lake_temp.f90 lake_convection.f90 sublake_temp.f90 soil_hydro.f90 water_check.f90 \
	          dyn_veg.f90 soil_carbon.f90 peat_carbon.f90 shelf_carbon.f90 ice_carbon.f90 lake_carbon.f90 \
	          carbon_trans.f90 dust_emis.f90 weathering.f90 carbon_export.f90 init_cell.f90 end_cell.f90 carbon_inventory.f90 carbon_flx_atm_lnd.f90 \
            lnd_model.f90 lnd_out.f90
tmp_lnd = $(patsubst %.f90, %.o, $(files_lnd) )
obj_lnd = $(patsubst %, $(objdir)/%, $(tmp_lnd) )
########################################################################

########################################################################
# land virtual cell (lndvc) model related source files
dir_lndvc = $(srcdir)/lndvc
files_lndvc = lndvc_def.f90 lndvc_grid.f90 lndvc_model.f90 
tmp_lndvc = $(patsubst %.f90, %.o, $(files_lndvc) )
obj_lndvc = $(patsubst %, $(objdir)/%, $(tmp_lndvc) )
########################################################################

########################################################################
# main climber source files
dir_main = $(srcdir)/main/
files_main = constants.f90 control.f90 timer.f90 climber_grid.f90 
tmp_main = $(patsubst %.f90, %.o, $(files_main) )
obj_main = $(patsubst %, $(objdir)/%, $(tmp_main) ) $(objdir)/coupler.o $(objdir)/cmn_out.o $(objdir)/climber.o
obj_main_clim = $(patsubst %, $(objdir)/%, $(tmp_main) ) $(objdir)/coupler_nobgc.o $(objdir)/cmn_out_nobgc.o $(objdir)/climber_clim.o
obj_main_clim_bgc = $(patsubst %, $(objdir)/%, $(tmp_main) ) $(objdir)/coupler.o $(objdir)/cmn_out.o $(objdir)/climber_clim_bgc.o
obj_main_clim_ice = $(patsubst %, $(objdir)/%, $(tmp_main) ) $(objdir)/coupler_nobgc.o $(objdir)/cmn_out_nobgc.o $(objdir)/climber_clim_ice.o
########################################################################

########################################################################
# ocean model related source files
dir_ocn = $(srcdir)/ocn/
files_ocn = ocn_model.f90 ocn_params.f90 ocn_grid.f90 ocn_def.f90 ocn_out.f90 ocn_check.f90 \
						momentum.f90 jbar.f90 ubarsolv.f90 island.f90 matinv.f90 wind.f90 invert.f90 velc.f90 \
            transport_ocn.f90 advection.f90 diffusion.f90 convection.f90 krausturner.f90 eos.f90 \
            restore_salinity.f90 hosing.f90 noise.f90 flux_adj.f90 ocn_grid_update_state.f90 \
						cfc_flux.f90 free_surface.f90 bering.f90
tmp_ocn = $(patsubst %.f90, %.o, $(files_ocn) )
obj_ocn = $(patsubst %, $(objdir)/%, $(tmp_ocn) )
########################################################################

########################################################################
# sea ice model related source files
dir_sic = $(srcdir)/sic/
files_sic = sic_model.f90 sic_grid.f90 sic_params.f90 sic_def.f90 sic_out.f90 \
            sic_dyn.f90 transport_sic.f90 surface_par_sic.f90 ebal_sic.f90 ebal_ocn.f90
tmp_sic = $(patsubst %.f90, %.o, $(files_sic) )
obj_sic = $(patsubst %, $(objdir)/%, $(tmp_sic) )
########################################################################

########################################################################
# smb related source files
dir_smb = $(srcdir)/smb/
files_smb = smb_model.f90 smb_params.f90 smb_grid.f90 smb_def.f90 smb_out.f90 \
				 	  topo.f90 ice.f90 downscaling.f90 smb_surface_par.f90 smb_ebal.f90 smb_temp.f90 snow.f90 smb_simple.f90 smb_pdd.f90 semi.f90 \
				 	  fake_atm_hires.f90 bias_corr.f90
tmp_smb = $(patsubst %.f90, %.o, $(files_smb) )
obj_smb = $(patsubst %, $(objdir)/%, $(tmp_smb) )
# dummy
files_smb_dummy = .smb_model_dummy.f90 smb_def.f90 .smb_out_dummy.f90
tmp_smb_dummy = $(patsubst %.f90, %.o, $(files_smb_dummy) )
obj_smb_dummy = $(patsubst %, $(objdir)/%, $(tmp_smb_dummy) )
########################################################################

########################################################################
# utilities source files
dir_utils = $(srcdir)/utils/
files_utils = precision.f90 dim_name.f90 tridiag.f90 filter.f90 #hyster.f90
tmp_utils = $(patsubst %.f90, %.o, $(files_utils) )
obj_utils = $(patsubst %, $(objdir)/%, $(tmp_utils) )
########################################################################

#############################################################
##							
## Rules for individual libraries or modules
##
#############################################################

#########################
# git version recording
.PHONY: $(dir_version)gitversion.f90
$(dir_version)gitversion.f90:
	sed -e "s/HEAD/`git rev-parse HEAD`/" $(dir_version)gitversion.f90.in > $(dir_version)gitversion.f90.new
	cmp -s $(dir_version)gitversion.f90.new $(dir_version)gitversion.f90 \
		|| mv $(dir_version)gitversion.f90.new $(dir_version)gitversion.f90
	rm -f $(dir_version)gitversion.f90.new 2>/dev/null

$(objdir)/gitversion.o : $(dir_version)gitversion.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

#########################
# main climber rules ####
$(objdir)/climber.o : $(dir_main)climber.f90 \
	$(objdir)/gitversion.o \
	$(objdir)/control.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/bnd.o $(objdir)/geo.o \
	$(objdir)/atm_model.o $(objdir)/ocn_model.o $(objdir)/sic_model.o $(objdir)/lnd_model.o \
	$(objdir)/bgc_model.o $(objdir)/co2_model.o $(objdir)/ch4_model.o \
	$(objdir)/ice_model.o $(objdir)/smb_model.o $(objdir)/bmb_model.o \
	$(objdir)/atm_out.o $(objdir)/ocn_out.o $(objdir)/sic_out.o $(objdir)/lnd_out.o \
	$(objdir)/bgc_out.o $(objdir)/co2_out.o $(objdir)/ch4_out.o \
	$(objdir)/smb_out.o $(objdir)/bmb_out.o $(objdir)/geo_out.o \
	$(objdir)/coupler.o $(objdir)/cmn_out.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/climber_clim.o : $(dir_main)climber.f90 \
	$(objdir)/gitversion.o \
	$(objdir)/control.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/bnd.o $(objdir)/geo.o \
	$(objdir)/atm_model.o $(objdir)/ocn_model.o $(objdir)/sic_model.o $(objdir)/lnd_model.o \
	$(objdir)/.bgc_model.o $(objdir)/co2_model.o $(objdir)/ch4_model.o \
	$(objdir)/.ice_model_dummy.o $(objdir)/.smb_model_dummy.o $(objdir)/.bmb_model_dummy.o \
	$(objdir)/atm_out.o $(objdir)/ocn_out.o $(objdir)/sic_out.o $(objdir)/lnd_out.o \
	$(objdir)/.bgc_out.o $(objdir)/co2_out.o $(objdir)/ch4_out.o \
	$(objdir)/.smb_out_dummy.o $(objdir)/.bmb_out_dummy.o $(objdir)/geo_out.o \
	$(objdir)/coupler_nobgc.o $(objdir)/cmn_out_nobgc.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/climber_clim_bgc.o : $(dir_main)climber.f90 \
	$(objdir)/gitversion.o \
	$(objdir)/control.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/bnd.o $(objdir)/geo.o \
	$(objdir)/atm_model.o $(objdir)/ocn_model.o $(objdir)/sic_model.o $(objdir)/lnd_model.o \
	$(objdir)/bgc_model.o $(objdir)/co2_model.o $(objdir)/ch4_model.o \
	$(objdir)/.ice_model_dummy.o $(objdir)/.smb_model_dummy.o $(objdir)/.bmb_model_dummy.o \
	$(objdir)/atm_out.o $(objdir)/ocn_out.o $(objdir)/sic_out.o $(objdir)/lnd_out.o \
	$(objdir)/bgc_out.o $(objdir)/co2_out.o $(objdir)/ch4_out.o \
	$(objdir)/.smb_out_dummy.o $(objdir)/.bmb_out_dummy.o $(objdir)/geo_out.o \
	$(objdir)/coupler.o $(objdir)/cmn_out.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/climber_clim_ice.o : $(dir_main)climber.f90 \
	$(objdir)/gitversion.o \
	$(objdir)/control.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/bnd.o $(objdir)/geo.o \
	$(objdir)/atm_model.o $(objdir)/ocn_model.o $(objdir)/sic_model.o $(objdir)/lnd_model.o \
	$(objdir)/.bgc_model.o $(objdir)/co2_model.o $(objdir)/ch4_model.o \
	$(objdir)/ice_model.o $(objdir)/smb_model.o $(objdir)/bmb_model.o \
	$(objdir)/atm_out.o $(objdir)/ocn_out.o $(objdir)/sic_out.o $(objdir)/lnd_out.o \
	$(objdir)/.bgc_out.o $(objdir)/co2_out.o $(objdir)/ch4_out.o \
	$(objdir)/smb_out.o $(objdir)/bmb_out.o $(objdir)/geo_out.o \
	$(objdir)/coupler_nobgc.o $(objdir)/cmn_out_nobgc.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/coupler.o : $(dir_main)coupler.f90 \
	$(objdir)/control.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/bnd.o \
	$(objdir)/atm_def.o $(objdir)/ocn_def.o $(objdir)/sic_def.o $(objdir)/lnd_def.o \
	$(objdir)/bgc_params.o $(objdir)/bgc_def.o $(objdir)/co2_def.o $(objdir)/ch4_def.o \
	$(objdir)/ice_def.o $(objdir)/smb_def.o $(objdir)/bmb_def.o $(objdir)/geo_def.o 
	$(FC) $(LDFLAGS) -c -o $@ $<	

$(objdir)/coupler_nobgc.o : $(dir_main)coupler.f90 \
	$(objdir)/control.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/bnd.o \
	$(objdir)/atm_def.o $(objdir)/ocn_def.o $(objdir)/sic_def.o $(objdir)/lnd_def.o \
	$(objdir)/.bgc_params.o $(objdir)/.bgc_def.o $(objdir)/co2_def.o $(objdir)/ch4_def.o \
	$(objdir)/ice_def.o $(objdir)/smb_def.o $(objdir)/bmb_def.o $(objdir)/geo_def.o 
	$(FC) $(LDFLAGS) -c -o $@ $<	

$(objdir)/cmn_out.o : $(dir_main)cmn_out.f90 $(objdir)/coupler.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/cmn_out_nobgc.o : $(dir_main)cmn_out.f90 $(objdir)/coupler_nobgc.o $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/timer.o : $(dir_main)timer.f90 $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/climber_grid.o : $(dir_main)climber_grid.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/control.o : $(dir_main)control.f90 $(objdir)/constants.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/constants.o : $(dir_main)constants.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

######################
# utilities rules ####
$(objdir)/precision.o : $(dir_utils)precision.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/%.o : $(dir_utils)%.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

#######################
# atmosphere rules ####
$(objdir)/atm_model.o : $(dir_atm)atm_model.f90 $(objdir)/atm_grid.o $(objdir)/atm_params.o $(objdir)/atm_def.o $(objdir)/adifa.o \
	$(objdir)/synop.o $(objdir)/smooth_atm.o $(objdir)/wvel.o $(objdir)/clouds.o $(objdir)/vesta.o \
	$(objdir)/crisa.o $(objdir)/slp.o $(objdir)/u2d.o $(objdir)/u3d.o $(objdir)/time_step.o \
	$(objdir)/lwr.o $(objdir)/swr.o $(objdir)/feedbacks.o $(objdir)/rad_kernels.o $(objdir)/dust.o
	$(FC) $(LDFLAGS) -c -o $@ $<
						
$(objdir)/smooth_atm.o : $(dir_atm)smooth_atm.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/atm_grid.o : $(dir_atm)atm_grid.f90 $(objdir)/climber_grid.o $(objdir)/constants.o \
					   $(objdir)/control.o $(objdir)/smooth_atm.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/atm_params.o : $(dir_atm)atm_params.f90 $(objdir)/constants.o $(objdir)/control.o $(objdir)/timer.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/atm_def.o : $(dir_atm)atm_def.f90 $(objdir)/atm_params.o $(objdir)/atm_grid.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/adifa.o : $(dir_atm)adifa.f90 $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/synop.o : $(dir_atm)synop.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/wvel.o : $(dir_atm)wvel.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/clouds.o : $(dir_atm)clouds.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/vesta.o : $(dir_atm)vesta.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/crisa.o : $(dir_atm)crisa.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o $(objdir)/smooth_atm.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/slp.o : $(dir_atm)slp.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o $(objdir)/smooth_atm.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/u2d.o : $(dir_atm)u2d.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/u3d.o : $(dir_atm)u3d.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/time_step.o : $(dir_atm)time_step.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o $(objdir)/smooth_atm.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lwr.o : $(dir_atm)lwr.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/swr.o : $(dir_atm)swr.f90 $(objdir)/constants.o $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/feedbacks.o : $(dir_atm)feedbacks.f90 $(objdir)/constants.o $(objdir)/control.o $(objdir)/timer.o $(objdir)/atm_grid.o $(objdir)/atm_params.o $(objdir)/lwr.o $(objdir)/swr.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/rad_kernels.o : $(dir_atm)rad_kernels.f90 $(objdir)/constants.o $(objdir)/control.o $(objdir)/atm_grid.o $(objdir)/atm_params.o $(objdir)/lwr.o $(objdir)/swr.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/dust.o : $(dir_atm)dust.f90 $(objdir)/atm_grid.o $(objdir)/atm_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/atm_out.o : $(dir_atm)atm_out.f90 $(objdir)/atm_grid.o $(objdir)/atm_params.o $(objdir)/atm_def.o $(objdir)/dim_name.o
	$(FC) $(LDFLAGS) -c -o $@ $<

##################
# ocean rules ####
$(objdir)/ocn_model.o : $(dir_ocn)ocn_model.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o $(objdir)/ocn_def.o $(objdir)/ocn_check.o \
	$(objdir)/momentum.o $(objdir)/transport_ocn.o $(objdir)/eos.o $(objdir)/restore_salinity.o \
	$(objdir)/hosing.o $(objdir)/noise.o $(objdir)/flux_adj.o $(objdir)/ocn_grid_update_state.o $(objdir)/cfc_flux.o \
	$(objdir)/free_surface.o $(objdir)/bering.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ocn_params.o : $(dir_ocn)ocn_params.f90 $(objdir)/timer.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ocn_grid.o : $(dir_ocn)ocn_grid.f90 $(objdir)/ocn_params.o $(objdir)/control.o $(objdir)/climber_grid.o $(objdir)/constants.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ocn_def.o : $(dir_ocn)ocn_def.f90 $(objdir)/precision.o $(dir_ocn)ocn_grid.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ocn_out.o : $(dir_ocn)ocn_out.f90 $(objdir)/ocn_def.o $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ocn_check.o : $(dir_ocn)ocn_check.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/momentum.o : $(dir_ocn)momentum.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o $(objdir)/wind.o $(objdir)/jbar.o \
	$(objdir)/ubarsolv.o $(objdir)/island.o $(objdir)/matinv.o $(objdir)/velc.o $(objdir)/invert.o  $(objdir)/eos.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/convection.o : $(dir_ocn)convection.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o $(objdir)/eos.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/eos.o : $(dir_ocn)eos.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/island.o : $(dir_ocn)island.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/matinv.o : $(dir_ocn)matinv.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/wind.o : $(dir_ocn)wind.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/jbar.o : $(dir_ocn)jbar.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ubarsolv.o : $(dir_ocn)ubarsolv.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/invert.o : $(dir_ocn)invert.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/krausturner.o : $(dir_ocn)krausturner.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o $(objdir)/eos.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/transport_ocn.o : $(dir_ocn)transport_ocn.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o \
	$(objdir)/advection.o $(objdir)/diffusion.o $(objdir)/eos.o $(objdir)/convection.o $(objdir)/krausturner.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/advection.o : $(dir_ocn)advection.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/diffusion.o : $(dir_ocn)diffusion.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/velc.o : $(dir_ocn)velc.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/restore_salinity.o : $(dir_ocn)restore_salinity.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/hosing.o : $(dir_ocn)hosing.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/noise.o : $(dir_ocn)noise.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/flux_adj.o : $(dir_ocn)flux_adj.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ocn_grid_update_state.o : $(dir_ocn)ocn_grid_update_state.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/cfc_flux.o : $(dir_ocn)cfc_flux.f90 $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/free_surface.o : $(dir_ocn)free_surface.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bering.o : $(dir_ocn)bering.f90 $(objdir)/ocn_params.o $(objdir)/ocn_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

####################
# sea ice rules ####
$(objdir)/sic_model.o : $(dir_sic)sic_model.f90 $(objdir)/sic_grid.o $(objdir)/sic_params.o $(objdir)/sic_def.o \
	$(objdir)/sic_dyn.o $(objdir)/transport_sic.o $(objdir)/surface_par_sic.o \
	$(objdir)/ebal_sic.o $(objdir)/ebal_ocn.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sic_grid.o : $(dir_sic)sic_grid.f90 $(objdir)/constants.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sic_params.o : $(dir_sic)sic_params.f90 $(objdir)/constants.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sic_def.o : $(dir_sic)sic_def.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sic_dyn.o : $(dir_sic)sic_dyn.f90 $(objdir)/sic_grid.o $(objdir)/sic_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/transport_sic.o : $(dir_sic)transport_sic.f90 $(objdir)/sic_grid.o $(objdir)/sic_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/surface_par_sic.o : $(dir_sic)surface_par_sic.f90 $(objdir)/sic_grid.o $(objdir)/sic_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ebal_sic.o : $(dir_sic)ebal_sic.f90 $(objdir)/sic_grid.o $(objdir)/sic_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ebal_ocn.o : $(dir_sic)ebal_ocn.f90 $(objdir)/sic_grid.o $(objdir)/sic_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sic_out.o : $(dir_sic)sic_out.f90 $(objdir)/sic_grid.o $(objdir)/sic_params.o $(objdir)/sic_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

#################
# land rules ####
$(objdir)/lnd_model.o : $(dir_lnd)lnd_model.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/lnd_def.o \
	$(objdir)/surface_par_lnd.o $(objdir)/ebal_veg.o $(objdir)/ebal_ice.o $(objdir)/ebal_lake.o \
	$(objdir)/soil_par.o $(objdir)/ice_par.o $(objdir)/shelf_par.o $(objdir)/lake_par.o $(objdir)/lake_rho.o \
	$(objdir)/soil_temp.o $(objdir)/ice_temp.o $(objdir)/shelf_temp.o $(objdir)/lake_temp.o $(objdir)/lake_convection.o $(objdir)/sublake_temp.o \
	$(objdir)/surface_hydro.o $(objdir)/soil_hydro.o $(objdir)/water_check.o \
	$(objdir)/veg_par.o $(objdir)/photosynthesis.o $(objdir)/dyn_veg.o \
	$(objdir)/soil_carbon.o $(objdir)/peat_carbon.o $(objdir)/shelf_carbon.o $(objdir)/ice_carbon.o $(objdir)/lake_carbon.o \
	$(objdir)/soil_carbon_par.o $(objdir)/shelf_carbon_par.o $(objdir)/ice_carbon_par.o $(objdir)/lake_carbon_par.o \
	$(objdir)/carbon_trans.o $(objdir)/dust_emis.o $(objdir)/weathering.o $(objdir)/carbon_export.o $(objdir)/init_cell.o $(objdir)/end_cell.o \
	$(objdir)/carbon_inventory.o $(objdir)/carbon_flx_atm_lnd.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lnd_grid.o : $(dir_lnd)lnd_grid.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lnd_params.o : $(dir_lnd)lnd_params.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lnd_def.o : $(dir_lnd)lnd_def.f90 $(objdir)/precision.o $(objdir)/lnd_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lnd_out.o : $(dir_lnd)lnd_out.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/lnd_def.o $(objdir)/veg_par.o
	$(FC) $(LDFLAGS) -c -o $@ $<
	
$(objdir)/veg_par.o : $(dir_lnd)veg_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/soil_carbon_par.o : $(dir_lnd)soil_carbon_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/shelf_carbon_par.o : $(dir_lnd)shelf_carbon_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice_carbon_par.o : $(dir_lnd)ice_carbon_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lake_carbon_par.o : $(dir_lnd)lake_carbon_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/soil_par.o : $(dir_lnd)soil_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice_par.o : $(dir_lnd)ice_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lake_par.o : $(dir_lnd)lake_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/lake_rho.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lake_rho.o : $(dir_lnd)lake_rho.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/shelf_par.o : $(dir_lnd)shelf_par.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/surface_par_lnd.o : $(dir_lnd)surface_par_lnd.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/veg_par.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/photosynthesis.o : $(dir_lnd)photosynthesis.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/veg_par.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ebal_veg.o : $(dir_lnd)ebal_veg.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/surface_par_lnd.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ebal_ice.o : $(dir_lnd)ebal_ice.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/surface_par_lnd.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ebal_lake.o : $(dir_lnd)ebal_lake.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/surface_par_lnd.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/surface_hydro.o : $(dir_lnd)surface_hydro.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/soil_temp.o : $(dir_lnd)soil_temp.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/tridiag.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice_temp.o : $(dir_lnd)ice_temp.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/shelf_temp.o : $(dir_lnd)shelf_temp.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lake_temp.o : $(dir_lnd)lake_temp.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/lake_convection.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lake_convection.o : $(dir_lnd)lake_convection.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/lake_rho.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sublake_temp.o : $(dir_lnd)sublake_temp.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/soil_hydro.o : $(dir_lnd)soil_hydro.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/veg_par.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/water_check.o : $(dir_lnd)water_check.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/dyn_veg.o : $(dir_lnd)dyn_veg.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/veg_par.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/soil_carbon.o : $(dir_lnd)soil_carbon.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/peat_carbon.o : $(dir_lnd)peat_carbon.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/shelf_carbon.o : $(dir_lnd)shelf_carbon.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice_carbon.o : $(dir_lnd)ice_carbon.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lake_carbon.o : $(dir_lnd)lake_carbon.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/carbon_trans.o : $(dir_lnd)carbon_trans.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/dust_emis.o : $(dir_lnd)dust_emis.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/weathering.o : $(dir_lnd)weathering.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/carbon_export.o : $(dir_lnd)carbon_export.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/init_cell.o : $(dir_lnd)init_cell.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/veg_par.o $(objdir)/surface_par_lnd.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/end_cell.o : $(dir_lnd)end_cell.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/veg_par.o $(objdir)/surface_par_lnd.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/carbon_inventory.o : $(dir_lnd)carbon_inventory.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/timer.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/carbon_flx_atm_lnd.o : $(dir_lnd)carbon_flx_atm_lnd.f90 $(objdir)/lnd_grid.o $(objdir)/lnd_params.o $(objdir)/timer.o
	$(FC) $(LDFLAGS) -c -o $@ $<

#######################
# lndvc rules #########

$(objdir)/lndvc_def.o : $(dir_lndvc)/lndvc_def.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lndvc_grid.o : $(dir_lndvc)/lndvc_grid.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lndvc_model.o : $(dir_lndvc)/lndvc_model.f90 $(objdir)/lndvc_def.o $(objdir)/lndvc_grid.o \
						$(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

##################################
# ocean biogeochemistry rules ####
$(objdir)/bgc_model.o : $(dir_bgc)bgc_model.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o \
	$(objdir)/mo_biomod.o $(objdir)/bgc_diag.o $(objdir)/apply_net_fw.o $(objdir)/chemcon.o \
	$(objdir)/mo_m4ago_types.o $(objdir)/mo_m4ago_params.o $(objdir)/mo_m4ago_physics.o \
	$(objdir)/mo_m4ago_core.o $(objdir)/mo_m4ago_HAMOCCinit.o $(objdir)/mo_m4ago_HAMOCCPrimPart.o \
	$(objdir)/mo_m4ago_climberx.o \
	$(objdir)/ocprod.o $(objdir)/carchm.o $(objdir)/cyano.o $(objdir)/powadi.o \
	$(objdir)/dipowa.o $(objdir)/powach.o $(objdir)/sedshi.o $(objdir)/mo_beleg_bgc.o $(objdir)/tracer_cons.o \
	$(objdir)/flx_sed.o $(objdir)/flx_sed_net.o $(objdir)/flx_bur.o $(objdir)/sed_grid_update_tracers.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bgc_grid.o : $(dir_bgc)/bgc_grid.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_m4ago_types.o : $(dir_m4ago)/mo_m4ago_types.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_m4ago_physics.o : $(dir_m4ago)/mo_m4ago_physics.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_m4ago_core.o : $(dir_m4ago)/mo_m4ago_core.f90 ${objdir}/mo_m4ago_types.o ${objdir}/mo_m4ago_physics.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_m4ago_params.o : $(dir_m4ago)/mo_m4ago_params.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_m4ago_HAMOCCinit.o : $(dir_m4ago)/mo_m4ago_HAMOCCinit.F90 ${objdir}/mo_m4ago_params.o ${objdir}/mo_m4ago_core.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bgc_params.o : $(dir_bgc)/bgc_params.f90 $(objdir)/bgc_grid.o ${objdir}/mo_m4ago_HAMOCCinit.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bgc_def.o : $(dir_bgc)/bgc_def.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_diag.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_biomod.o : $(dir_bgc)/mo_biomod.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bgc_diag.o : $(dir_bgc)/bgc_diag.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/apply_net_fw.o : $(dir_bgc)/apply_net_fw.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/chemcon.o : $(dir_bgc)/chemcon.f90 $(objdir)/bgc_def.o $(objdir)/bgc_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_m4ago_HAMOCCPrimPart.o : $(dir_m4ago)/mo_m4ago_HAMOCCPrimPart.F90 ${objdir}/mo_m4ago_types.o ${objdir}/mo_m4ago_HAMOCCinit.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_m4ago_climberx.o : $(dir_m4ago)/mo_m4ago_climberx.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o \
	${objdir}/mo_m4ago_types.o $(objdir)/mo_m4ago_params.o $(objdir)/mo_m4ago_core.o  $(objdir)/mo_m4ago_physics.o \
	${objdir}/mo_m4ago_HAMOCCinit.o $(objdir)/mo_m4ago_HAMOCCPrimPart.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ocprod.o : $(dir_bgc)/ocprod.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o ${objdir}/bgc_diag.o \
	$(objdir)/mo_m4ago_physics.o $(objdir)/mo_m4ago_HAMOCCinit.o ${objdir}/mo_m4ago_climberx.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/carchm.o : $(dir_bgc)/carchm.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o $(objdir)/bgc_diag.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/cyano.o : $(dir_bgc)/cyano.f90 $(objdir)/bgc_diag.o $(objdir)/bgc_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/powadi.o : $(dir_bgc)/powadi.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/dipowa.o : $(dir_bgc)/dipowa.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/powach.o : $(dir_bgc)/powach.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o $(objdir)/bgc_diag.o \
	$(objdir)/powadi.o $(objdir)/dipowa.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sedshi.o : $(dir_bgc)/sedshi.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mo_beleg_bgc.o : $(dir_bgc)/mo_beleg_bgc.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/tracer_cons.o : $(dir_bgc)/tracer_cons.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/flx_sed.o : $(dir_bgc)/flx_sed.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o $(objdir)/bgc_diag.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/flx_sed_net.o : $(dir_bgc)/flx_sed_net.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o $(objdir)/bgc_diag.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/flx_bur.o : $(dir_bgc)/flx_bur.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o $(objdir)/bgc_diag.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sed_grid_update_tracers.o : $(dir_bgc)/sed_grid_update_tracers.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bgc_out.o : $(dir_bgc)bgc_out.f90 $(objdir)/bgc_grid.o $(objdir)/bgc_params.o $(objdir)/bgc_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.bgc_model.o : $(dir_bgc_dummy).bgc_model.f90 $(objdir)/.bgc_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.bgc_out.o : $(dir_bgc_dummy).bgc_out.f90 $(objdir)/.bgc_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.bgc_params.o : $(dir_bgc_dummy).bgc_params.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.bgc_def.o : $(dir_bgc_dummy).bgc_def.f90 $(objdir)/.bgc_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

################
# geo rules ###
$(objdir)/geo.o : $(dir_geo)geo.f90 $(objdir)/climber_grid.o $(objdir)/control.o \
	$(objdir)/geo_params.o $(objdir)/geo_grid.o $(objdir)/geo_def.o \
	$(objdir)/runoff_routing.o $(objdir)/lakes.o $(objdir)/fill_ocean.o $(objdir)/topo_filter.o $(objdir)/topo_fill.o $(objdir)/connect_ocn.o \
	$(objdir)/q_geo.o $(objdir)/sed.o $(objdir)/vilma.o $(objdir)/gia.o $(objdir)/fix_runoff.o \
	$(objdir)/hires_to_lowres.o $(objdir)/corals_topo.o $(objdir)/coast_cells.o $(objdir)/drainage_basins.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/vilma.o : $(dir_geo)vilma.F90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.vilma_dummy.o : $(dir_geo).vilma_dummy.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/gia.o : $(dir_geo)gia.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/q_geo.o : $(dir_geo)q_geo.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sed.o : $(dir_geo)sed.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/corals_topo.o : $(dir_geo)corals_topo.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fix_runoff.o : $(dir_geo)fix_runoff.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/hires_to_lowres.o : $(dir_geo)hires_to_lowres.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/coast_cells.o : $(dir_geo)coast_cells.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/drainage_basins.o : $(dir_geo)drainage_basins.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/runoff_routing.o : $(dir_geo)runoff_routing.f90 $(objdir)/lakes.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/lakes.o : $(dir_geo)lakes.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fill_ocean.o : $(dir_geo)fill_ocean.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/topo_filter.o : $(dir_geo)topo_filter.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/topo_fill.o : $(dir_geo)topo_fill.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/connect_ocn.o : $(dir_geo)connect_ocn.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/geo_params.o : $(dir_geo)geo_params.f90 $(objdir)/control.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/geo_grid.o : $(dir_geo)geo_grid.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/geo_def.o : $(dir_geo)geo_def.f90 $(objdir)/precision.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/geo_out.o : $(dir_geo)geo_out.f90 $(objdir)/geo.o $(objdir)/control.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

################
# co2 rules ####
$(objdir)/co2_model.o : $(dir_co2)co2_model.f90 $(objdir)/co2_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/co2_def.o : $(dir_co2)co2_def.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/co2_out.o : $(dir_co2)co2_out.f90 $(objdir)/co2_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

################
# ch4 rules ####
$(objdir)/ch4_model.o : $(dir_ch4)ch4_model.f90 $(objdir)/ch4_def.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ch4_def.o : $(dir_ch4)ch4_def.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ch4_out.o : $(dir_ch4)ch4_out.f90 $(objdir)/ch4_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

######################
# ice model rules ####
$(objdir)/ice_model.o : $(dir_ice)ice_model.f90 $(objdir)/ice_def.o $(objdir)/ice_id.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.ice_model_dummy.o : $(dir_ice).ice_model_dummy.f90 $(objdir)/ice_def.o $(objdir)/ice_id.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice_def.o : $(dir_ice)ice_def.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice_id.o : $(dir_ice)ice_id.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

############################
# SICOPOLIS model rules ####
$(objdir)/sico_model.o : $(dir_ice_sico)sico_model.f90 $(objdir)/control.o $(objdir)/timer.o \
	$(objdir)/sico_types_m.o $(objdir)/sico_params.o $(objdir)/sico_check.o \
	$(objdir)/sico_timer.o $(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/ice_material_properties_m.o \
	$(objdir)/init_temp_water_age_m.o $(objdir)/calc_vxy_m.o $(objdir)/calc_vz_m.o $(objdir)/calc_dxyz_m.o \
	$(objdir)/topograd_m.o $(objdir)/calc_thk_m.o $(objdir)/enth_temp_omega_m.o $(objdir)/calc_temp_m.o \
	$(objdir)/calc_temp_enth_m.o $(objdir)/calc_enhance_m.o $(objdir)/flag_update_gf_gl_cf_m.o $(objdir)/calc_temp_melt_bas_m.o \
	$(objdir)/calc_bas_melt_m.o $(objdir)/calc_thk_water_bas_m.o $(objdir)/boundary_m.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_types_m.o : $(dir_ice_sico)sico_types_m.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/compare_float_m.o : $(dir_ice_sico)compare_float_m.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_params.o : $(dir_ice_sico)sico_params.f90 $(objdir)/sico_types_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_timer.o : $(dir_ice_sico)sico_timer.f90 $(objdir)/sico_types_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_state.o : $(dir_ice_sico)sico_state.f90 $(objdir)/sico_types_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_maths_m.o : $(dir_ice_sico)sico_maths_m.f90 $(objdir)/sico_types_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice_material_properties_m.o : $(dir_ice_sico)ice_material_properties_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/init_temp_water_age_m.o : $(dir_ice_sico)init_temp_water_age_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o \
	$(objdir)/sico_state.o $(objdir)/sico_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/stereo_proj_m.o : $(dir_ice_sico)stereo_proj_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_grid.o : $(dir_ice_sico)sico_grid.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o $(objdir)/sico_timer.o \
	$(objdir)/ice_material_properties_m.o $(objdir)/stereo_proj_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_check.o : $(dir_ice_sico)sico_check.f90 $(objdir)/sico_state.o $(objdir)/sico_grid.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_vxy_m.o : $(dir_ice_sico)calc_vxy_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_timer.o $(objdir)/ice_material_properties_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_vz_m.o : $(dir_ice_sico)calc_vz_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_grid.o $(objdir)/sico_state.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_dxyz_m.o : $(dir_ice_sico)calc_dxyz_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/topograd_m.o : $(dir_ice_sico)topograd_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_timer.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_thk_m.o : $(dir_ice_sico)calc_thk_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_timer.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_params.o $(objdir)/sico_maths_m.o $(objdir)/topograd_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/enth_temp_omega_m.o : $(dir_ice_sico)enth_temp_omega_m.f90 $(objdir)/sico_types_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_temp_m.o : $(dir_ice_sico)calc_temp_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_timer.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_params.o $(objdir)/sico_maths_m.o $(objdir)/ice_material_properties_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_temp_enth_m.o : $(dir_ice_sico)calc_temp_enth_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_timer.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_params.o $(objdir)/sico_maths_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_enhance_m.o : $(dir_ice_sico)calc_enhance_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_timer.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/flag_update_gf_gl_cf_m.o : $(dir_ice_sico)flag_update_gf_gl_cf_m.f90 $(objdir)/sico_types_m.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_temp_melt_bas_m.o : $(dir_ice_sico)calc_temp_melt_bas_m.f90 $(objdir)/sico_types_m.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_bas_melt_m.o : $(dir_ice_sico)calc_bas_melt_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_timer.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_params.o $(objdir)/ice_material_properties_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/calc_thk_water_bas_m.o : $(dir_ice_sico)calc_thk_water_bas_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o \
	$(objdir)/sico_state.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/hydro_m.o : $(dir_ice_sico)hydro_m.f90 $(objdir)/sico_types_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/discharge_workers_m.o : $(dir_ice_sico)discharge_workers_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/boundary_m.o : $(dir_ice_sico)boundary_m.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/discharge_workers_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sico_out.o : $(dir_ice_sico)sico_out.f90 $(objdir)/sico_types_m.o $(objdir)/sico_params.o \
	$(objdir)/sico_grid.o $(objdir)/sico_state.o $(objdir)/sico_model.o $(objdir)/enth_temp_omega_m.o
	$(FC) $(LDFLAGS) -c -o $@ $<

################
# smb rules ####
$(objdir)/smb_model.o : $(dir_smb)smb_model.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o $(objdir)/smb_def.o \
	$(objdir)/topo.o $(objdir)/ice.o $(objdir)/smb_simple.o $(objdir)/smb_pdd.o $(objdir)/semi.o \
	$(objdir)/fake_atm_hires.o $(objdir)/bias_corr.o $(objdir)/filter.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_grid.o : $(dir_smb)smb_grid.f90 $(objdir)/smb_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_params.o : $(dir_smb)smb_params.f90 $(objdir)/timer.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_def.o : $(dir_smb)smb_def.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_surface_par.o : $(dir_smb)smb_surface_par.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/topo.o : $(dir_smb)topo.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ice.o : $(dir_smb)ice.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/downscaling.o : $(dir_smb)downscaling.f90 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_ebal.o : $(dir_smb)smb_ebal.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_temp.o : $(dir_smb)smb_temp.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/snow.o : $(dir_smb)snow.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_simple.o : $(dir_smb)smb_simple.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_pdd.o : $(dir_smb)smb_pdd.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o $(objdir)/downscaling.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/semi.o : $(dir_smb)semi.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o $(objdir)/downscaling.o \
	$(objdir)/downscaling.o $(objdir)/smb_ebal.o $(objdir)/smb_temp.o $(objdir)/snow.o $(objdir)/smb_surface_par.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_atm_hires.o : $(dir_smb)fake_atm_hires.f90 $(objdir)/smb_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bias_corr.o : $(dir_smb)bias_corr.f90 $(objdir)/timer.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/smb_out.o : $(dir_smb)smb_out.f90 $(objdir)/smb_grid.o $(objdir)/smb_params.o $(objdir)/smb_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.smb_model_dummy.o : $(dir_smb).smb_model_dummy.f90 $(objdir)/smb_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.smb_out_dummy.o : $(dir_smb).smb_out_dummy.f90 $(objdir)/smb_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

################
# bmb rules ####
$(objdir)/bmb_model.o : $(dir_bmb)bmb_model.f90 $(objdir)/bmb_grid.o $(objdir)/bmb_params.o $(objdir)/bmb_def.o 
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bmb_grid.o : $(dir_bmb)bmb_grid.f90 $(objdir)/bmb_params.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bmb_params.o : $(dir_bmb)bmb_params.f90 $(objdir)/timer.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bmb_def.o : $(dir_bmb)bmb_def.f90 $(objdir)/precision.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/bmb_out.o : $(dir_bmb)bmb_out.f90 $(objdir)/bmb_grid.o $(objdir)/bmb_params.o $(objdir)/bmb_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.bmb_model_dummy.o : $(dir_bmb).bmb_model_dummy.f90 $(objdir)/bmb_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/.bmb_out_dummy.o : $(dir_bmb).bmb_out_dummy.f90 $(objdir)/bmb_def.o
	$(FC) $(LDFLAGS) -c -o $@ $<

#####################
# boundary rules ####
$(objdir)/bnd.o : $(dir_bnd)bnd.f90 $(objdir)/insolation.o $(objdir)/solar.o $(objdir)/volc.o \
	$(objdir)/co2.o $(objdir)/co2_rad.o $(objdir)/ch4.o $(objdir)/ch4_rad.o $(objdir)/n2o.o $(objdir)/so4.o $(objdir)/o3.o \
	$(objdir)/cfc.o $(objdir)/luc.o $(objdir)/dist.o $(objdir)/sea_level.o $(objdir)/d13c_atm.o $(objdir)/D14c_atm.o \
	$(objdir)/fake_atm.o $(objdir)/fake_dust.o $(objdir)/fake_lnd.o $(objdir)/fake_ocn.o $(objdir)/fake_sic.o $(objdir)/fake_ice.o $(objdir)/fake_geo.o \
	$(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/constants.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/insolation.o : $(dir_bnd)insolation.f90 $(objdir)/constants.o $(objdir)/timer.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/solar.o : $(dir_bnd)solar.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/volc.o : $(dir_bnd)volc.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/co2.o : $(dir_bnd)co2.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/co2_rad.o : $(dir_bnd)co2_rad.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/sea_level.o : $(dir_bnd)sea_level.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/d13c_atm.o : $(dir_bnd)d13c_atm.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/D14c_atm.o : $(dir_bnd)D14c_atm.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ch4.o : $(dir_bnd)ch4.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/ch4_rad.o : $(dir_bnd)ch4_rad.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/n2o.o : $(dir_bnd)n2o.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/so4.o : $(dir_bnd)so4.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/o3.o : $(dir_bnd)o3.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/cfc.o : $(dir_bnd)cfc.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/luc.o : $(dir_bnd)luc.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/dist.o : $(dir_bnd)dist.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_atm.o : $(dir_bnd)fake_atm.f90 $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_dust.o : $(dir_bnd)fake_dust.f90 $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_lnd.o : $(dir_bnd)fake_lnd.f90 $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_ocn.o : $(dir_bnd)fake_ocn.f90 $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_sic.o : $(dir_bnd)fake_sic.f90 $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_ice.o : $(dir_bnd)fake_ice.f90 $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/fake_geo.o : $(dir_bnd)fake_geo.f90 $(objdir)/timer.o $(objdir)/climber_grid.o $(objdir)/control.o
	$(FC) $(LDFLAGS) -c -o $@ $<

########################################################################
# Combinations of all object files to build given configurations

# climber-clim: minimal configuration with ocn,atm,lnd,sic #

climber_clim_obj = $(obj_utils) $(obj_bnd) $(obj_geo) \
			  $(obj_atm) $(obj_ocn) $(obj_sic) $(obj_lnd) \
				$(obj_bgc_dummy) $(obj_co2) $(obj_ch4) \
			  $(obj_ice_dummy) $(obj_smb_dummy) $(obj_bmb_dummy) $(obj_main_clim) 

# climber-clim-bgc: clim plus with bgc #

climber_clim_bgc_obj = $(obj_utils) $(obj_bnd) $(obj_geo)\
			  $(obj_atm) $(obj_ocn) $(obj_sic) $(obj_lnd) \
				$(obj_bgc) $(obj_m4ago) $(obj_co2) $(obj_ch4) \
			  $(obj_ice_dummy) $(obj_smb_dummy) $(obj_bmb_dummy) $(obj_main_clim_bgc)

# climber-clim-ice: clim plus with ice #

climber_clim_ice_obj = $(obj_utils) $(obj_bnd) $(obj_geo) \
			  $(obj_atm) $(obj_ocn) $(obj_sic) $(obj_lnd) \
				$(obj_bgc_dummy) $(obj_co2) $(obj_ch4) \
			  $(obj_ice) $(obj_ice_sico) $(obj_smb) $(obj_bmb) $(obj_main_clim_ice)

# climber-clim-bgc-ice: clim plus with bgc and ice #

climber_clim_bgc_ice_obj = $(obj_utils) $(obj_bnd) $(obj_geo) \
			  $(obj_atm) $(obj_ocn) $(obj_sic) $(obj_lnd) \
				$(obj_bgc) $(obj_m4ago) $(obj_co2) $(obj_ch4) \
			  $(obj_ice) $(obj_ice_sico) $(obj_smb) $(obj_bmb) $(obj_main)

########################################################################



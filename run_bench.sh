#!/bin/bash

############################################################################
# script to run a complete set of CLIMBER-X simulations for model benchmark
############################################################################

# Specify desidered name of output directory
outdir=bench_v1.4.0

# A complete model benchmark involves the following steps:
# step=1    : preindustrial equilibrium spinup, climate only
#           + do manual copy of restart files into restart/restart_pi/
# step=2    : preindustrial equilibrium spinup for biogeochemistry
#           + do manual copy of bgc and lnd restart files into restart/restart_pi/
#           + all climate-only simulations 
#           + do manual copy of LGM restart files into restart/restart_lgm/
#           + compute 2m temperature bias from the hist_biascorr run using the matlab script in matlab/t2m_bias/t2mb.m and copy to input/ directory
# step=3    : all simulations with carbon cycle and interactive ice sheets 
step=3

# Specify carbon cycle setup: 'close' or 'open'
cc_setup='close'


############################################################################

echo outdir $outdir

echo step $step

echo cc_setup $cc_setup


# closed carbon cycle
cc_string_close_bgc='bgc.l_sediments=F' 
cc_string_close_lnd='lnd.l_river_export=F'
# open with conserved Phosphorus
cc_string_open_bgc='bgc.l_sediments=T bgc.l_conserve_phos=T bgc.l_conserve_sil=T' 
cc_string_open_lnd='lnd.l_river_export=F'
# open 
cc_string_open2_bgc='bgc.l_sediments=T bgc.l_conserve_phos=F bgc.l_conserve_sil=T' 
cc_string_open2_lnd='lnd.l_river_export=T'

if [ $cc_setup='close' ]
then
  cc_string_bgc=$cc_string_close_bgc
  cc_string_lnd=$cc_string_close_lnd
elif [ $cc_setup='open' ]
then
  cc_string_bgc=$cc_string_open_bgc
  cc_string_lnd=$cc_string_open_lnd
elif [ $cc_setup='open2' ]
then
  cc_string_bgc=$cc_string_open2_bgc
  cc_string_lnd=$cc_string_open2_lnd
fi

echo cc_string_bgc $cc_string_bgc
echo cc_string_lnd $cc_string_lnd

####################
# step 1 
if [ $step -eq 1 ]
then
# preindustrial equilibrium spinup, climate only
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/preind -p ctl.nyears=10000 ctl.i_write_restart=1
fi

####################
# step 2
if [ $step -eq 2 ]
then

# preindustrial equilibrium with ocean biogeochemistry, closed carbon cycle setup
./runme -rs -q medium -w 60:00:00 --omp 32 -o output/$outdir/preind_cc_close -p ctl.nyears=10000 ctl.flag_bgc=T ctl.bgc_restart=F ctl.i_write_restart=1 bgc.l_spinup_bgc=F bgc.l_spinup_sed=F bgc.l_sediments=F lnd.l_river_export=F
# preindustrial equilibrium with ocean biogeochemistry, open carbon cycle setup
./runme -rs -q long -w 220:00:00 --omp 32 -o output/$outdir/preind_cc_open -p ctl.nyears=100000 ctl.flag_bgc=T ctl.l_spinup_cc=T ctl.nyears_spinup_bgc=5000 ctl.year_start_offline=1000000 ctl.bgc_restart=F ctl.l_weathering=T ctl.i_write_restart=1 bgc.l_spinup_bgc=T bgc.l_spinup_sed=T ${cc_string_open_bgc} ${cc_string_open_lnd}

# historical with annual atm output to compute T2m anomalies needed for anomaly approach in SEMI
./runme -rs -q short -w 10:00:00 --omp 32 -o output/$outdir/hist_biascorr -p ctl.nyears=220 ctl.year_ini=-200 ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ctl.year_out_start=-19 ctl.nyout_atm=1

# 1/4xCO2
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/025xCO2 -p ctl.nyears=10000 ctl.co2_const=70  atm.l_daily_output=T ocn.l_daily_output=T sic.l_daily_output=T
# 1/2xCO2
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/05xCO2  -p ctl.nyears=10000 ctl.co2_const=140 atm.l_daily_output=T ocn.l_daily_output=T sic.l_daily_output=T
# 180 ppm
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/180ppm  -p ctl.nyears=10000 ctl.co2_const=180
# 2xCO2
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/2xCO2   -p ctl.nyears=10000 ctl.co2_const=560 atm.l_daily_output=T ocn.l_daily_output=T sic.l_daily_output=T
# 4xCO2
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/4xCO2   -p ctl.nyears=10000 ctl.co2_const=1120 atm.l_daily_output=T ocn.l_daily_output=T sic.l_daily_output=T

# with fixed vegetation for vegetation feedback determination
# 1/2xCO2
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/05xCO2_fixveg -p ctl.nyears=5000 ctl.co2_const=140 lnd.l_dynveg=F lnd.l_fixlai=T lnd.l_co2_fert_lim=T
# 180 ppm
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/180ppm_fixveg -p ctl.nyears=5000 ctl.co2_const=180 lnd.l_dynveg=F lnd.l_fixlai=T lnd.l_co2_fert_lim=T
# 2xCO2
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/2xCO2_fixveg  -p ctl.nyears=5000 ctl.co2_const=560 lnd.l_dynveg=F lnd.l_fixlai=T lnd.l_co2_fert_lim=T

# feedback analysis
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/fb/CO2_1x_2x         -p ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=280
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/fb/CO2_1x_2x_fixveg  -p ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=280 lnd.l_dynveg=F lnd.l_fixlai=T lnd.l_co2_fert_lim=T
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/fb/CO2_05x_1x        -p ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=140
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/fb/CO2_05x_1x_fixveg -p ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=140 lnd.l_dynveg=F lnd.l_fixlai=T lnd.l_co2_fert_lim=T
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/fb/CO2_2x_4x         -p ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=560
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/fb/CO2_2x_4x_fixveg  -p ctl.l_feedbacks=T ctl.nyears=10000 ctl.co2_const=560 lnd.l_dynveg=F lnd.l_fixlai=T lnd.l_co2_fert_lim=T

# abrupt 0p5xCO2
./runme -rs -q short -w 1:00:00 --omp 32 -o output/$outdir/abrupt0p5xCO2 -p ctl.nyears=150 ctl.co2_const=140  ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# abrupt 2xCO2
./runme -rs -q short -w 1:00:00 --omp 32 -o output/$outdir/abrupt2xCO2   -p ctl.nyears=150 ctl.co2_const=560  ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# abrupt 4xCO2
./runme -rs -q short -w 1:00:00 --omp 32 -o output/$outdir/abrupt4xCO2   -p ctl.nyears=150 ctl.co2_const=1120 ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# abrupt 4% increase in solar constant
./runme -rs -q short -w 1:00:00 --omp 32 -o output/$outdir/abruptsolp4p  -p ctl.nyears=150 ctl.sol_const=1416.5 ctl.nyout_atm=150 ctl.nyout_ocn=150 ctl.nyout_sic=150 ctl.nyout_lnd=150

# aquaplanet
./runme -rs -q short -w 1:00:00 --omp 32 -o output/$outdir/aqua -p ctl.nyears=10 ctl.l_aquaplanet=T ctl.flag_ocn=F ctl.flag_sic=F ctl.flag_lnd=F ctl.flag_bgc=F ctl.nyout_atm=10

# mid-Holocene
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/midHolo        -p ctl.nyears=5000 ctl.iorbit=1 ctl.ecc_const=0.018682 ctl.obl_const=24.105 ctl.per_const=0.87 ctl.co2_const=264.4 ctl.ch4_const=597 ctl.n2o_const=262 
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/midHolo_fixveg -p ctl.nyears=5000 ctl.iorbit=1 ctl.ecc_const=0.018682 ctl.obl_const=24.105 ctl.per_const=0.87 ctl.co2_const=264.4 ctl.ch4_const=597 ctl.n2o_const=262 lnd.l_dynveg=F lnd.l_co2_fert_lim=T lnd.l_fixlai=T

# last interglacial
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/lig127k        -p ctl.nyears=5000 ctl.iorbit=1 ctl.ecc_const=0.039378 ctl.obl_const=24.04 ctl.per_const=275.41 ctl.co2_const=275 ctl.ch4_const=685 ctl.n2o_const=255 
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/lig127k_fixveg -p ctl.nyears=5000 ctl.iorbit=1 ctl.ecc_const=0.039378 ctl.obl_const=24.04 ctl.per_const=275.41 ctl.co2_const=275 ctl.ch4_const=685 ctl.n2o_const=255 lnd.l_dynveg=F lnd.l_co2_fert_lim=T lnd.l_fixlai=T

# last interglacial (125 ka) with bgc
./runme -rs -q medium -w 48:00:00 --omp 32 -o output/$outdir/lig125k_cc_close -p ctl.year_ini=-125000 ctl.nyears=15000  ctl.flag_bgc=T ctl.co2_const=278 ctl.ch4_const=640 ctl.n2o_const=260 ctl.bgc_restart=F ctl.i_write_restart=1  ${cc_string_close_bgc} ${cc_string_close_lnd}
./runme -rs -q medium             --omp 32 -o output/$outdir/lig125k_cc_open  -p ctl.year_ini=-125000 ctl.nyears=100000 ctl.flag_bgc=T ctl.co2_const=278 ctl.ch4_const=640 ctl.n2o_const=260 ctl.l_spinup_cc=T ctl.nyears_spinup_bgc=5000 ctl.year_start_offline=7000 ctl.bgc_restart=F ctl.l_weathering=T ctl.i_write_restart=1 bgc.l_spinup_bgc=T bgc.l_spinup_sed=T ${cc_string_open_bgc} ${cc_string_open_lnd}

# LGM
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/lgm           -p ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc         ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 lnd.lithology_uhh_file=input/Lithology_lgm_UHH.nc
#./runme -rs -q short -w 24:00:00 --omp 32 -o output/bench_v1.3/lgm_deglac -p ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_tarasov_deglac_21ka.nc ctl.fake_geo_ref_file=input/geo_ice_tarasov_deglac_0ka.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_deglac_21ka.nc ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 lnd.lithology_uhh_file=input/Lithology_lgm_UHH.nc
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/lgm_fixveg    -p ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc         ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 lnd.l_dynveg=F lnd.l_co2_fert_lim=T lnd.l_fixlai=T
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/lgm_peltier   -p ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_peltier_lgm.nc         ctl.fake_geo_ref_file=input/geo_ice_peltier_deglac_0ka.nc ctl.fake_ice_const_file=input/geo_ice_peltier_lgm.nc ctl.co2_const=190 ctl.ch4_const=375 ctl.n2o_const=200 

# hysteresis
./runme -rs -q medium -w 60:00:00 --omp 32 -o output/$outdir/hyst2050/up   -p ctl.nyears=30000 ocn.l_hosing=T ocn.hosing_ini=-0.3 ocn.hosing_trend=0.02  ocn.lat_min_hosing=20 ocn.lat_max_hosing=50
./runme -rs -q medium -w 60:00:00 --omp 32 -o output/$outdir/hyst2050/down -p ctl.nyears=30000 ocn.l_hosing=T ocn.hosing_ini=0.3  ocn.hosing_trend=-0.02 ocn.lat_min_hosing=20 ocn.lat_max_hosing=50
./runme -rs -q medium -w 60:00:00 --omp 32 -o output/$outdir/hyst5070/up   -p ctl.nyears=30000 ocn.l_hosing=T ocn.hosing_ini=-0.3 ocn.hosing_trend=0.02  ocn.lat_min_hosing=50 ocn.lat_max_hosing=70
./runme -rs -q medium -w 60:00:00 --omp 32 -o output/$outdir/hyst5070/down -p ctl.nyears=30000 ocn.l_hosing=T ocn.hosing_ini=0.3  ocn.hosing_trend=-0.02 ocn.lat_min_hosing=50 ocn.lat_max_hosing=70

# PaleoDEM
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/pliomip -p ctl.nyears=5000 ctl.fake_geo_const_file=input/geo_Pliomip2.nc         ctl.fake_ice_const_file=input/geo_Pliomip2.nc ctl.co2_const=400 geo.geo_ref_file=input/geo_Pliomip2.nc
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/65Ma    -p ctl.nyears=5000 ctl.fake_geo_const_file=input/topog_065Ma_climberX.nc ctl.fake_ice_const_file=input/topog_065Ma_climberX.nc ctl.atm_restart=F ctl.ocn_restart=F ctl.sic_restart=F ctl.lnd_restart=F ctl.dt_day_ocn=1 ocn.i_isl=0 ocn.i_init=1 geo.l_close_panama=F geo.geo_ref_file=input/topog_065Ma_climberX.nc
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/250Ma   -p ctl.nyears=5000 ctl.fake_geo_const_file=input/topog_250Ma_climberX.nc ctl.fake_ice_const_file=input/topog_250Ma_climberX.nc ctl.atm_restart=F ctl.ocn_restart=F ctl.sic_restart=F ctl.lnd_restart=F ctl.dt_day_ocn=1 ocn.i_isl=0 ocn.i_init=1 geo.l_close_panama=F geo.geo_ref_file=input/topog_250Ma_climberX.nc 
fi

####################
# step 3
if [ $step -eq 3 ]
then

# historical 
./runme -rs -q short -w 6:00:00 --omp 32 -o output/$outdir/hist       -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ctl.flag_ch4=T ctl.ich4_rad=2 ctl.ich4_emis=1 ctl.id13c=1 ctl.iD14c=1 ctl.year_out_start=-19 ctl.nyout_atm=1 ctl.nyout_ocn=1 ctl.nyout_sic=1 ctl.nyout_lnd=1 ctl.nyout_bgc=1 ocn.l_cfc=T ${cc_string_bgc} ${cc_string_lnd} 
./runme -rs -q short -w 6:00:00 --omp 32 -o output/$outdir/hist-nat   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=0 ctl.ich4=0 ctl.in2o=0 ctl.icfc=0 ctl.io3=1 ctl.iso4=0 ctl.iluc=1 ctl.isol=1 ctl.ivolc=1 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 6:00:00 --omp 32 -o output/$outdir/hist-ghg   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=0 ctl.iluc=1 ctl.isol=0 ctl.ivolc=0 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 6:00:00 --omp 32 -o output/$outdir/hist-aer   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=0 ctl.ich4=0 ctl.in2o=0 ctl.icfc=0 ctl.io3=1 ctl.iso4=1 ctl.iluc=1 ctl.isol=0 ctl.ivolc=0 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 6:00:00 --omp 32 -o output/$outdir/hist-noluc -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=1 ctl.isol=1 ctl.ivolc=1 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 6:00:00 --omp 32 -o output/$outdir/hist-esm   -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.iorbit=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ${cc_string_bgc} ${cc_string_lnd}

# historical with SEMI
./runme -rs -q short -w 10:00:00 --omp 32 -o output/$outdir/hist_smb_grl -p ctl.nyears=215 ctl.year_ini=-200 ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ctl.year_out_start=-19 ctl.n_year_smb=1 ctl.flag_smb=T ctl.nyout_smb=1 ctl.ice_domain_name=GRL-8KM smb.l_monthly_output=T smb.i_alb_ice=0

# SSP scenarios
./runme -rs -q short -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp119 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp119.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp119.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp119.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp119.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp126 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp126.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp126.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp126.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp126.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp245 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp245.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp245.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp245.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp245.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp370 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp370.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp370.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp370.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp370.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp460 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp460.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp460.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp460.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp460.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp585 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp585.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp585.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp585.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp585.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}

# SSP scenarios, emission driven
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/ssp_emis/ssp119 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp119.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp119_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp119.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp119.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1  ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/ssp_emis/ssp126 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp126.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp126_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp126.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp126.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1  ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/ssp_emis/ssp245 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp245.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp245_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp245.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp245.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1  ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/ssp_emis/ssp370 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp370.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp370_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp370.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp370.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1  ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/ssp_emis/ssp460 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp460.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp460_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp460.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp460.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1  ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/ssp_emis/ssp585 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp585.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp585_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp585.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp585.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1  ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ${cc_string_bgc} ${cc_string_lnd}

# 1%/yr CO2 increase
./runme -rs -q short -w 1:00:00  --omp 32 -o output/$outdir/1pctCO2/cpl      -p ctl.nyears=140 ctl.ico2=2 ctl.ico2_rad=0 ctl.co2_const=280 ctl.flag_bgc=T ctl.nyout_atm=140 ctl.nyout_ocn=140 ctl.nyout_sic=140 ctl.nyout_lnd=140 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 10:00:00 --omp 32 -o output/$outdir/1pctCO2/cpl_out1 -p ctl.nyears=140 ctl.ico2=2 ctl.ico2_rad=0 ctl.co2_const=280 ctl.flag_bgc=T ctl.nyout_atm=140 ctl.nyout_ocn=1   ctl.nyout_sic=140 ctl.nyout_lnd=140 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 1:00:00  --omp 32 -o output/$outdir/1pctCO2/bgc      -p ctl.nyears=140 ctl.ico2=2 ctl.ico2_rad=1 ctl.co2_const=280 ctl.flag_bgc=T ctl.nyout_atm=140 ctl.nyout_ocn=140 ctl.nyout_sic=140 ctl.nyout_lnd=140 ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q short -w 1:00:00  --omp 32 -o output/$outdir/1pctCO2/rad      -p ctl.nyears=140 ctl.ico2=0 ctl.ico2_rad=3 ctl.co2_const=280 ctl.flag_bgc=T ctl.nyout_atm=140 ctl.nyout_ocn=140 ctl.nyout_sic=140 ctl.nyout_lnd=140 ${cc_string_bgc} ${cc_string_lnd}

# ZECMIP
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/zecmip -p ctl.nyears=10000 ctl.flag_co2=T ctl.ico2_emis=4 ctl.co2_pulse=1000 ctl.flag_bgc=T

# preindustrial with interactive CO2, to check for drift
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/preind_co2_close -p ctl.nyears=10000 ctl.flag_co2=T ctl.flag_bgc=T ${cc_string_close_bgc} ${cc_string_close_lnd}

# LGM with prognostic CO2
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/lgm_co2_close_Cicebur   -p ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.flag_co2=T ctl.flag_bgc=T ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.ch4_const=375 ctl.n2o_const=200 lnd.k_ice=10  ${cc_string_close_bgc} ${cc_string_close_lnd}
./runme -rs -q short -w 24:00:00 --omp 32 -o output/$outdir/lgm_co2_close_Cicenobur -p ctl.nyears=10000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.flag_co2=T ctl.flag_bgc=T ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.ch4_const=375 ctl.n2o_const=200 lnd.k_ice=1e6 ${cc_string_close_bgc} ${cc_string_close_lnd}

# deglaciation
#./runme -rs -q medium --omp 32 -o output/$outdir/deglac_tarasov ctl.nyears=21000 ctl.year_ini=-21000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.iorbit=2 ctl.fake_geo_var_file=input/geo_ice_tarasov_lgc.nc ctl.fake_ice_var_file=input/geo_ice_tarasov_lgc.nc ctl.n_year_geo=100 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_geo=100 ctl.restart_in_dir=restart_lgm
#./runme -rs -q medium --omp 32 -o output/$outdir/deglac_peltier ctl.nyears=21000 ctl.year_ini=-21000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.iorbit=2 ctl.fake_geo_var_file=input/geo_ice_peltier_deglac.nc ctl.fake_geo_ref_file=input/geo_ice_peltier_deglac_0ka.nc ctl.fake_ice_var_file=input/geo_ice_peltier_deglac.nc ctl.n_year_geo=100 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_geo=100 ctl.restart_in_dir=restart_lgm_peltier
#./runme -rs -q medium --omp 32 -o output/$outdir/deglac_gowan   ctl.nyears=21000 ctl.year_ini=-21000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.iorbit=2 ctl.fake_geo_var_file=input/geo_ice_gowan.nc ctl.fake_geo_ref_file=input/geo_ice_gowan_0ka.nc ctl.fake_ice_var_file=input/geo_ice_gowan.nc ctl.n_year_geo=100 ctl.nyout_atm=100 nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_geo=100 ctl.restart_in_dir=restart_lgm

# last glacial cycle with prescribed ice sheets
./runme -rs -q long -w 240:00:00 --omp 32 -o output/$outdir/lgc_tarasov -p ctl.nyears=125000 ctl.year_ini=-125000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.orbit=2 ctl.fake_geo_var_file=input/geo_ice_tarasov_lgc.nc ctl.fake_ice_var_file=input/geo_ice_tarasov_lgc.nc ctl.n_year_geo=100 ocn.l_noise_fw=T ocn.noise_amp_fw=0 
./runme -rs -q long -w 240:00:00 --omp 32 -o output/$outdir/lgc_tarasov_noise -p ctl.nyears=125000 ctl.year_ini=-125000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.orbit=2 ctl.fake_geo_var_file=input/geo_ice_tarasov_lgc.nc ctl.fake_ice_var_file=input/geo_ice_tarasov_lgc.nc ctl.n_year_geo=100 ocn.l_noise_fw=T ocn.noise_amp_fw=0.5 

# last glacial cycle with interactive ice sheets
#./runme -rs -q long   -w 500:00:00 --omp 32 -o output/$outdir/lgc_ice_sico        -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=1  ctl.ice_model_name=sico 
#./runme -rs -q long   -w 500:00:00 --omp 32 -o output/$outdir/lgc_ice_sico_acc2   -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=2  ctl.ice_model_name=sico 
#./runme -rs -q medium -w 168:00:00 --omp 32 -o output/$outdir/lgc_ice_sico_acc5   -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=5  ctl.ice_model_name=sico 
#./runme -rs -q medium -w 168:00:00 --omp 32 -o output/$outdir/lgc_ice_sico_acc10  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=10 ctl.ice_model_name=sico 
#./runme -rs -q long   -w 500:00:00 --omp 32 -o output/$outdir/lgc_ice_yelmo       -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=1  ctl.ice_model_name=yelmo 
#./runme -rs -q long   -w 500:00:00 --omp 32 -o output/$outdir/lgc_ice_yelmo_acc2  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=2  ctl.ice_model_name=yelmo 
#./runme -rs -q medium -w 168:00:00 --omp 32 -o output/$outdir/lgc_ice_yelmo_acc5  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=5  ctl.ice_model_name=yelmo 
#./runme -rs -q medium -w 168:00:00 --omp 32 -o output/$outdir/lgc_ice_yelmo_acc10 -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=10 ctl.ice_model_name=yelmo 

#
#./runme -rs -q standby --omp 32 -o output/$outdir/deglac_ice_acc10_itempinit1 -p ctl.year_ini=-30000 ctl.nyears=30000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=10 ctl.ifake_ice=2 ctl.ifake_geo=2 ices.i_temp_init=1
#./runme -rs -q standby --omp 32 -o output/$outdir/deglac_ice_acc10_itempinit4 -p ctl.year_ini=-30000 ctl.nyears=30000 ctl.iorbit=2 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.ice_domain_name=NH-32KM ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.n_accel=10 ctl.ifake_ice=2 ctl.ifake_geo=2 ices.i_temp_init=4
#
## preind with interactive ice sheets
#./runme -rs -q standby --omp 32 -o output/$outdir/preind_ice_nh -p ctl.nyears=100000 ctl.n_accel=10 ctl.n_year_smb=10 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM
## MIS11 with interactive ice sheets
#./runme -rs -q standby --omp 32 -o output/$outdir/mis11_ice_nh  -p ctl.year_ini=-398000 ctl.nyears=100000 ctl.n_accel=10 ctl.n_year_smb=10 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM

fi

####################
# step 5
if [ $step -eq 5 ]
then

# preind with prognostic CO2 and open carbon cycle, to check for drift
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/preind_co2_open -p ctl.nyears=10000 ctl.flag_co2=T ctl.flag_bgc=T ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ctl.l_c14=F ${cc_string_open_bgc} ${cc_string_open_lnd}

# historical 
./runme -rs -q standby -w 10:00:00 --omp 32 -o output/$outdir/hist_cc_open     -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.ico2=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ctl.flag_ch4=T ctl.ich4_rad=2 ctl.ich4_emis=1 ctl.id13c=1 ctl.iD14c=1 ctl.year_out_start=-19 ctl.nyout_atm=1 ctl.nyout_ocn=1 ctl.nyout_sic=1 ctl.nyout_lnd=1 ctl.nyout_bgc=1 ctl.restart_in_dir=restart_pi_cc_open ctl.l_weathering=T ocn.l_cfc=T ${cc_string_open_bgc} ${cc_string_open_lnd} 
./runme -rs -q standby -w 10:00:00 --omp 32 -o output/$outdir/hist-esm_cc_open -p ctl.nyears=1015 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.iorbit=1 ctl.ich4=1 ctl.in2o=1 ctl.icfc=1 ctl.io3=2 ctl.iso4=1 ctl.iluc=2 ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}

# ZECMIP
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/zecmip_cc_open   -p ctl.nyears=10000 ctl.flag_co2=T ctl.ico2_emis=4 ctl.co2_pulse=1000 ctl.flag_bgc=T ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}

# SSP scenarios, emission driven
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open/ssp119 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp119.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp119_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp119.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp119.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open/ssp126 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp126.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp126_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp126.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp126.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open/ssp245 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp245.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp245_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp245.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp245.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open/ssp370 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp370.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp370_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp370.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp370.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open/ssp460 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp460.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp460_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp460.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp460.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
./runme -rs -q standby -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open/ssp585 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp585.nc ctl.flag_ch4=T ctl.ich4_emis=1 ctl.ch4_emis_file=input/ch4_emis_historical_ssp585_ext.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp585.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp585.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}

# SSP scenarios, emission driven
#./runme -rs -q medium -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp119 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp119.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp119.nc  ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp119.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp119.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
#./runme -rs -q medium -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp126 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp126.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp126.nc  ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp126.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp126.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
#./runme -rs -q medium -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp245 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp245.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp245.nc  ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp245.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp245.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
#./runme -rs -q medium -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp370 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp370.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp370.nc  ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp370.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp370.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
#./runme -rs -q medium -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp460 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp460.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp460.nc  ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp460.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp460.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
#./runme -rs -q medium -w 48:00:00 --omp 32 -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp585 -p ctl.nyears=6000 ctl.year_ini=-1000 ctl.flag_bgc=T ctl.iorbit=1 ctl.flag_co2=T ctl.ico2_emis=1 ctl.id13C_emis=1 ctl.co2_emis_file=input/co2_emis_hist_ssp585.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp585.nc  ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp585.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp585.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}

fi

####################
# step 6
if [ $step -eq 6 ]
then

co2_volc_lgc=0.05617464

# control to check for drift
# close carbon cycle
./runme -rs -q standby --omp 32 -o output/$outdir/lig125k_co2_close_control       -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=25000 co2_const=278 flag_co2=T flag_bgc=T ch4_const=640 n2o_const=260 restart_in_dir=restart_lig_cc_close l_c14=F" ${cc_string_close_bgc} ${cc_string_close_lnd}
# open carbon cycle
./runme -rs -q standby --omp 32 -o output/$outdir/lig125k_co2_open_volLIG_control -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=25000 co2_const=278 flag_co2=T flag_bgc=T ch4_const=640 n2o_const=260 ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
./runme -rs -q standby --omp 32 -o output/$outdir/lig125k_co2_open_volLGC_control -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=25000 co2_const=278 flag_co2=T flag_bgc=T ch4_const=640 n2o_const=260 ico2_degas=0 co2_degas_const=$co2_volc_lgc l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# inception with interactive ice sheets and one-way coupled CO2 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_close             -p ctl.year_ini=-125000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close ${cc_string_close_bgc} ${cc_string_close_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_close_nolnd       -p ctl.year_ini=-125000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLIG       -p ctl.year_ini=-125000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLIG_nolnd -p ctl.year_ini=-125000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_lnd_co2=F ${cc_string_open_bgc} ${cc_string_open_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLGC       -p ctl.year_ini=-125000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=0 ctl.co2_degas_const=$co2_volc_lgc ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLGC_nolnd -p ctl.year_ini=-125000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=0 ctl.co2_degas_const=$co2_volc_lgc ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_lnd_co2=F ${cc_string_open_bgc} ${cc_string_open_lnd} 

./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_close             -p ctl.year_ini=-120000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=280 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_pi ${cc_string_close_bgc} ${cc_string_close_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_close_nolnd       -p ctl.year_ini=-120000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=280 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_pi ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volPI        -p ctl.year_ini=-120000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=280 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volPI_nolnd  -p ctl.year_ini=-120000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=280 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ctl.l_lnd_co2=F ${cc_string_open_bgc} ${cc_string_open_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volLGC       -p ctl.year_ini=-120000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=280 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=0 ctl.co2_degas_const=$co2_volc_lgc ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd} 
./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volLGC_nolnd -p ctl.year_ini=-120000 ctl.nyears=50000 ctl.iorbit=2 ctl.co2_const=280 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=0 ctl.co2_degas_const=$co2_volc_lgc ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open ctl.l_lnd_co2=F ${cc_string_open_bgc} ${cc_string_open_lnd} 

# inception with fully interactive ice sheets and CO2 
#./runme -rs -q standby --omp 32 -o output/$outdir/inception_co2/incep_ice_co2_close -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close ctl.l_c14=F ${cc_string_close_bgc} ${cc_string_close_lnd}
#./runme -rs -q standby --omp 32 -o output/$outdir/inception_co2/incep_ice_co2_open  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_c14=F ${cc_string_open_bgc} ${cc_string_open_lnd}

# inception with fixed ice sheets and one-way coupled CO2 
#./runme -rs -q standby --omp 32 -o output/$outdir/inception_co2/incep_co2_rad_close       -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close ${cc_string_close_bgc} ${cc_string_close_lnd}
#./runme -rs -q standby --omp 32 -o output/$outdir/inception_co2/incep_co2_rad_close_nolnd -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd}
#./runme -rs -q medium  --omp 32 -o output/$outdir/inception_co2/incep_co2_rad_open        -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
#./runme -rs -q standby --omp 32 -o output/$outdir/inception_co2/incep_co2_rad_open_nolnd  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_lnd_co2=F ${cc_string_open_bgc} ${cc_string_open_lnd}

## inception with fixed ice sheets and interactive CO2 
#./runme -rs -q standby --omp 32 -o output/$outdir/inception_co2/incep_co2_close -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close ctl.l_c14=F ${cc_string_close_bgc} ${cc_string_close_lnd}
#./runme -rs -q standby --omp 32 -o output/$outdir/inception_co2/incep_co2_open  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_c14=F ${cc_string_open_bgc} ${cc_string_open_lnd}


## lgm with prognostic CO2 and open carbon cycle
#./runme -rs -q standby -w 72:00:00 --omp 32 -o output/$outdir/lgm_co2_open -p ctl.nyears=15000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.flag_co2=T ctl.flag_bgc=T ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.ch4_const=375 ctl.n2o_const=200 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open lnd.k_ice=10 ${cc_string_open_bgc} ${cc_string_open_lnd} lnd.lithology_uhh_file=input/Lithology_lgm_UHH.nc

##
##### deglaciation with interactive CO2
### ./runme -rs -q medium --omp 32 -o output/$outdir/deglac_co2_close         -p ctl.nyears=21000 ctl.year_ini=-21000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.flag_co2=T ctl.flag_bgc=T ctl.co2_restart=T ctl.ico2_rad=2 ctl.ich4=1 ctl.in2o=1 ctl.iorbit=2 ctl.restart_in_dir=restart_lgm_cc_close ${cc_string_close_bgc} ${cc_string_close_lnd}
### ./runme -rs -q medium --omp 32 -o output/$outdir/deglac_co2_close_control -p ctl.nyears=21000 ctl.year_ini=-21000 ctl.flag_co2=T ctl.flag_bgc=T ctl.co2_restart=T ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ch4_const=375 ctl.n2o_const=200 ctl.restart_in_dir=restart_lgm_cc_close ${cc_string_close_bgc} ${cc_string_close_lnd}
####./runme -rs -q medium --omp 32 -o output/$outdir/deglac_co2_open          -p ctl.nyears=21000 ctl.year_ini=-21000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.flag_co2=T ctl.flag_bgc=T ctl.co2_restart=T ctl.ico2_rad=2 ctl.ich4=1 in2o=1 ctl.iorbit=2ico2_degas=1 ctl.restart_in_dir=restart_lgm_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
####./runme -rs -q medium --omp 32 -o output/$outdir/deglac_co2_open_control  -p ctl.nyears=21000 ctl.year_ini=-21000 ctl.flag_co2=T ctl.flag_bgc=T ctl.co2_restart=T ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ch4_const=375 ctl.n2o_const=200 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lgm_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd}
####./runme -rs -q medium --omp 32 -o output/$outdir/deglac_co2_open2         -p ctl.nyears=21000 ctl.year_ini=-21000 ctl.ifake_ice=1 ctl.ifake_geo=1 ctl.flag_co2=T ctl.flag_bgc=T ctl.co2_restart=T ctl.ico2_rad=2 ctl.ich4=1 in2o=1 ctl.iorbit=2ico2_degas=1 ctl.restart_in_dir=restart_lgm_cc_open2 ${cc_string_open2_bgc} ${cc_string_open2_lnd}
####./runme -rs -q medium --omp 32 -o output/$outdir/deglac_co2_open2_control -p ctl.nyears=21000 ctl.year_ini=-21000 ctl.flag_co2=T ctl.flag_bgc=T ctl.co2_restart=T ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ch4_const=375 ctl.n2o_const=200 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lgm_cc_open2 ${cc_string_open2_bgc} ${cc_string_open2_lnd}
##
fi

####################
# step 7
if [ $step -eq 7 ]
then

./runme -rs -q standby --omp 32 -o output/$outdir/incep_ice_co2_rad_open_seddiss -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ${cc_string_open_bgc} ${cc_string_open_lnd} bgc.l_sed_dry_dissolve=T

#./runme -rs -q standby --omp 32 -o output/$outdir/incep_ice_co2_open_dT2  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_c14=F ${cc_string_open_bgc} ${cc_string_open_lnd} smb.t2m_bias_corr_uniform=-2
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_ice_co2_open_fw01 -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_c14=F ${cc_string_open_bgc} ${cc_string_open_lnd} ocn.l_hosing=T ocn.hosing_ini=0.1 ocn.year_hosing_ini=6000
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_ice_co2_open_fw02 -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_c14=F ${cc_string_open_bgc} ${cc_string_open_lnd} ocn.l_hosing=T ocn.hosing_ini=0.2 ocn.year_hosing_ini=6000

#./runme -rs -q standby --omp 32 -o output/$outdir/incep_ice_co2_rad_close_nolnd          -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd}

#./runme -rs -q standby --omp 32 -o output/$outdir/incep_co2_close_nolnd                  -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close ctl.l_c14=F ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd}

#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_ice_co2_rad_close       -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close_fixRrain ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_ice_co2_rad_close_nolnd -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close_fixRrain l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2
#
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_co2_close       -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close_fixRrain ctl.l_c14=F  ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_co2_close_nolnd -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close_fixRrain ctl.l_c14=F ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2

#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_corals_ice_co2_rad_close       -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close_fixRrain_corals ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2 bgc.l_corals=T
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_corals_ice_co2_rad_close_nolnd -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.restart_in_dir=restart_lig_cc_close_fixRrain_corals ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2 bgc.l_corals=T
#
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_corals_co2_close       -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close_fixRrain_corals ctl.l_c14=F ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2 bgc.l_corals=T
#./runme -rs -q standby --omp 32 -o output/$outdir/incep_fixRrain_corals_co2_close_nolnd -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.restart_in_dir=restart_lig_cc_close_fixRrain_corals ctl.l_c14=F ctl.l_lnd_co2=F ${cc_string_close_bgc} ${cc_string_close_lnd} bgc.calmax=0.07 bgc.i_delcar=2 bgc.l_corals=T

#./runme -rs -q standby --omp 32 -o output/$outdir/incep_ice_co2_rad_open_nolnd          -p ctl.year_ini=-125000 ctl.nyears=125000 ctl.iorbit=2 ctl.co2_const=278 ctl.flag_co2=T ctl.ico2_rad=2 ctl.flag_bgc=T ctl.ich4=1 ctl.in2o=1 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.ice_domain_name=NH-32KM ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_lig_cc_open ctl.l_lnd_co2=F ${cc_string_open_bgc} ${cc_string_open_lnd}

fi

####################
# step 8
if [ $step -eq 8 ]
then

./runme -rs -q standby -w 72:00:00 --omp 32 -o output/$outdir/lgm_co2_open2 ctl.nyears=15000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.flag_co2=T ctl.flag_bgc=T ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.ch4_const=375 ctl.n2o_const=200 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open lnd.k_ice=10 ${cc_string_open_bgc} ${cc_string_open_lnd}
#./runme -rs -q standby -w 72:00:00 --omp 32 -o output/$outdir/lgm_co2_open ctl.nyears=15000 ctl.iorbit=1 ctl.ecc_const=0.018994 ctl.obl_const=22.949 ctl.per_const=114.42 ctl.flag_co2=T ctl.flag_bgc=T ctl.fake_geo_const_file=input/geo_ice_tarasov_lgm.nc ctl.fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ctl.ico2_rad=1 ctl.co2_rad_const=190 ctl.ch4_const=375 ctl.n2o_const=200 ctl.ico2_degas=1 ctl.l_weathering=T ctl.restart_in_dir=restart_pi_cc_open lnd.k_ice=10 ${cc_string_open_bgc} ${cc_string_open_lnd} lnd.lithology_uhh_file=input/Lithology_lgm_UHH.nc

fi

####################
# step 9
if [ $step -eq 9 ]
then

drag_topo_fac=3.5
restart_in_dir=restart_pi_v22

# SSP scenarios
./runme -rs -q standby -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp119 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp119.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp119.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp119.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp119.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q standby -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp126 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp126.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp126.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp126.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp126.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q standby -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp245 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp245.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp245.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp245.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp245.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q standby -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp370 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp370.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp370.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp370.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp370.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q standby -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp460 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp460.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp460.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp460.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp460.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac ${cc_string_bgc} ${cc_string_lnd}
./runme -rs -q standby -w 20:00:00 --omp 32 -o output/$outdir/ssp_conc/ssp585 -p ctl.nyears=3000 ctl.year_ini=-1000 ctl.iorbit=1 ctl.ico2=1 ctl.co2_file=input/co2_Koehler2017_ssp585.nc ctl.ich4=1 ctl.ch4_file=input/ch4_Koehler2017_ssp585.nc ctl.in2o=1 ctl.n2o_file=input/n2o_Koehler2017_ssp585.nc ctl.icfc=1 ctl.cfc_file=input/CFCs_historical_ssp585.nc ctl.io3=2 ctl.iso4=1 ctl.so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc ctl.iluc=2 ctl.luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc ctl.isol=1 ctl.ivolc=1 ctl.nyout_atm=100 ctl.nyout_ocn=100 ctl.nyout_sic=100 ctl.nyout_lnd=100 ctl.nyout_bgc=100 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac ${cc_string_bgc} ${cc_string_lnd}

./runme -rs -q standby -w 24:00:00 --omp 32 -o output/$outdir/pliomip_400     -p ctl.year_ini=0 ctl.fake_geo_const_file=input/geo_Pliomip2.nc ctl.fake_ice_const_file=input/geo_Pliomip2.nc ctl.co2_const=400 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac
./runme -rs -q standby -w 24:00:00 --omp 32 -o output/$outdir/pliomip_500     -p ctl.year_ini=0 ctl.fake_geo_const_file=input/geo_Pliomip2.nc ctl.fake_ice_const_file=input/geo_Pliomip2.nc ctl.co2_const=500 ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac

jobrun ./runme -rs -q standby --omp 32 -a -o output/$outdir/DOevents/     -p ctl.year_ini=-70000 ctl.nyears=50000 ctl.iorbit=2 ctl.ice_domain_name=NH-32KM ctl.co2_const=180,190,200,210,220,230,240 ctl.flag_geo=T ctl.flag_ice=T ctl.flag_smb=T ctl.flag_bmb=T ctl.restart_in_dir=$restart_in_dir ocn.drag_topo_fac=$drag_topo_fac

fi


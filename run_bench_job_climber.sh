#!/bin/bash

############################################################################
# script to run a complete set of CLIMBER-X simulations for model benchmark
############################################################################

# Specify desidered name of output directory
outdir=bench_v1.3

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
cc_string_close_bgc='bgc_par=l_sediments=F' 
cc_string_close_lnd='lnd_par=l_river_export=F'
# open with conserved Phosphorus
cc_string_open_bgc='bgc_par=l_sediments=T l_conserve_phos=T l_conserve_sil=T' 
cc_string_open_lnd='lnd_par=l_river_export=F'
# open 
cc_string_open2_bgc='bgc_par=l_sediments=T l_conserve_phos=F l_conserve_sil=T' 
cc_string_open2_lnd='lnd_par=l_river_export=T'

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
./job_climber -s -f -o output/$outdir/preind -c short -j parallel -n 32 \&control="nyears=10000 i_write_restart=1" 
fi

####################
# step 2
if [ $step -eq 2 ]
then

# preindustrial equilibrium with ocean biogeochemistry, closed carbon cycle setup
./job_climber -s -f -o output/$outdir/preind_cc_close -c medium -w 60 -j parallel -n 32 \&control="nyears=10000 flag_bgc=T bgc_restart=F i_write_restart=1" bgc_par="l_spinup_bgc=F l_spinup_sed=F l_sediments=F" lnd_par="l_river_export=F"
# preindustrial equilibrium with ocean biogeochemistry, open carbon cycle setup
./job_climber -s -f -o output/$outdir/preind_cc_open -c long -w 220 -j parallel -n 32 \&control="nyears=100000 flag_bgc=T l_spinup_cc=T nyears_spinup_bgc=5000 year_start_offline=1000000 bgc_restart=F l_weathering=T i_write_restart=1" bgc_par="l_spinup_bgc=T l_spinup_sed=T" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# historical with annual atm output to compute T2m anomalies needed for anomaly approach in SEMI
./job_climber -s -f -o output/$outdir/hist_biascorr -c short -w 10 -j parallel -n 32 \&control="nyears=220 year_ini=-200 iorbit=1 ico2=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=1 iluc=2 isol=1 ivolc=1 year_out_start=-19 nyout_atm=1"

# 1/4xCO2
./job_climber -s -f -o output/$outdir/025xCO2 -c short -w 24 -j parallel -n 32 \&control="nyears=10000 co2_const=70" atm_par="l_daily_output=T" ocn_par="l_daily_output=T" sic_par="l_daily_output=T"
# 1/2xCO2
./job_climber -s -f -o output/$outdir/05xCO2  -c short -w 24 -j parallel -n 32 \&control="nyears=10000 co2_const=140" atm_par="l_daily_output=T" ocn_par="l_daily_output=T" sic_par="l_daily_output=T"
# 180 ppm
./job_climber -s -f -o output/$outdir/180ppm  -c short -w 24 -j parallel -n 32 \&control="nyears=10000 co2_const=180"
# 2xCO2
./job_climber -s -f -o output/$outdir/2xCO2   -c short -w 24 -j parallel -n 32 \&control="nyears=10000 co2_const=560" atm_par="l_daily_output=T" ocn_par="l_daily_output=T" sic_par="l_daily_output=T"
# 4xCO2
./job_climber -s -f -o output/$outdir/4xCO2   -c short -w 24 -j parallel -n 32 \&control="nyears=10000 co2_const=1120" atm_par="l_daily_output=T" ocn_par="l_daily_output=T" sic_par="l_daily_output=T"

# with fixed vegetation for vegetation feedback determination
# 1/2xCO2
./job_climber -s -f -o output/$outdir/05xCO2_fixveg -c short -w 24 -j parallel -n 32 \&control="nyears=5000 co2_const=140" lnd_par="l_dynveg=F l_fixlai=T l_co2_fert_lim=T"
# 180 ppm
./job_climber -s -f -o output/$outdir/180ppm_fixveg -c short -w 24 -j parallel -n 32 \&control="nyears=5000 co2_const=180" lnd_par="l_dynveg=F l_fixlai=T l_co2_fert_lim=T"
# 2xCO2
./job_climber -s -f -o output/$outdir/2xCO2_fixveg  -c short -w 24 -j parallel -n 32 \&control="nyears=5000 co2_const=560" lnd_par="l_dynveg=F l_fixlai=T l_co2_fert_lim=T"

# feedback analysis
./job_climber -s -f -o output/$outdir/fb/CO2_1x_2x         -c short -w 24 -j parallel -n 32 \&control="l_feedbacks=T nyears=10000 co2_const=280"
./job_climber -s -f -o output/$outdir/fb/CO2_1x_2x_fixveg  -c short -w 24 -j parallel -n 32 \&control="l_feedbacks=T nyears=10000 co2_const=280" lnd_par="l_dynveg=F l_fixlai=T l_co2_fert_lim=T"
./job_climber -s -f -o output/$outdir/fb/CO2_05x_1x        -c short -w 24 -j parallel -n 32 \&control="l_feedbacks=T nyears=10000 co2_const=140"
./job_climber -s -f -o output/$outdir/fb/CO2_05x_1x_fixveg -c short -w 24 -j parallel -n 32 \&control="l_feedbacks=T nyears=10000 co2_const=140" lnd_par="l_dynveg=F l_fixlai=T l_co2_fert_lim=T"
./job_climber -s -f -o output/$outdir/fb/CO2_2x_4x         -c short -w 24 -j parallel -n 32 \&control="l_feedbacks=T nyears=10000 co2_const=560"
./job_climber -s -f -o output/$outdir/fb/CO2_2x_4x_fixveg  -c short -w 24 -j parallel -n 32 \&control="l_feedbacks=T nyears=10000 co2_const=560" lnd_par="l_dynveg=F l_fixlai=T l_co2_fert_lim=T"

# abrupt 0p5xCO2
./job_climber -s -f -o output/$outdir/abrupt0p5xCO2 -c short -w 1 -j parallel -n 32 \&control="nyears=150 co2_const=140 nyout_atm=150 nyout_ocn=150 nyout_sic=150 nyout_lnd=150"

# abrupt 2xCO2
./job_climber -s -f -o output/$outdir/abrupt2xCO2 -c short -w 1 -j parallel -n 32 \&control="nyears=150 co2_const=560 nyout_atm=150 nyout_ocn=150 nyout_sic=150 nyout_lnd=150"

# abrupt 4xCO2
./job_climber -s -f -o output/$outdir/abrupt4xCO2 -c short -w 1 -j parallel -n 32 \&control="nyears=150 co2_const=1120 nyout_atm=150 nyout_ocn=150 nyout_sic=150 nyout_lnd=150"

# abrupt 4% increase in solar constant
./job_climber -s -f -o output/$outdir/abruptsolp4p -c short -w 1 -j parallel -n 32 \&control="nyears=150 sol_const=1416.5 nyout_atm=150 nyout_ocn=150 nyout_sic=150 nyout_lnd=150"

# aquaplanet
./job_climber -s -f -o output/$outdir/aqua -c short -w 1 -j parallel -n 32 \&control="nyears=10 l_aquaplanet=T flag_ocn=F flag_sic=F flag_lnd=F flag_bgc=F nyout_atm=10"

# mid-Holocene
./job_climber -s -f -o output/$outdir/midHolo -c short -w 24 -j parallel -n 32 \&control="nyears=5000 iorbit=1 ecc_const=0.018682 obl_const=24.105 per_const=0.87 co2_const=264.4 ch4_const=597 n2o_const=262" 
./job_climber -s -f -o output/$outdir/midHolo_fixveg -c short -w 24 -j parallel -n 32 \&control="nyears=5000 iorbit=1 ecc_const=0.018682 obl_const=24.105 per_const=0.87 co2_const=264.4 ch4_const=597 n2o_const=262" lnd_par="l_dynveg=F l_co2_fert_lim=T l_fixlai=T"

# last interglacial
./job_climber -s -f -o output/$outdir/lig127k -c short -w 24 -j parallel -n 32 \&control="nyears=5000 iorbit=1 ecc_const=0.039378 obl_const=24.04 per_const=275.41 co2_const=275 ch4_const=685 n2o_const=255" 
./job_climber -s -f -o output/$outdir/lig127k_fixveg -c short -w 24 -j parallel -n 32 \&control="nyears=5000 iorbit=1 ecc_const=0.039378 obl_const=24.04 per_const=275.41 co2_const=275 ch4_const=685 n2o_const=255" lnd_par="l_dynveg=F l_co2_fert_lim=T l_fixlai=T"

# last interglacial (125 ka) with bgc
./job_climber -s -f -o output/$outdir/lig125k_cc_close -c medium -w 48 -j parallel -n 32 \&control="year_ini=-125000 nyears=15000 flag_bgc=T co2_const=278 ch4_const=640 n2o_const=260 bgc_restart=F i_write_restart=1"  "${cc_string_close_bgc}" "${cc_string_close_lnd}"
./job_climber -s -f -o output/$outdir/lig125k_cc_open -c medium -j parallel -n 32 \&control="year_ini=-125000 nyears=100000 flag_bgc=T co2_const=278 ch4_const=640 n2o_const=260 l_spinup_cc=T nyears_spinup_bgc=5000 year_start_offline=7000 bgc_restart=F l_weathering=T i_write_restart=1" bgc_par="l_spinup_bgc=T l_spinup_sed=T" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 

# LGM
./job_climber -s -f -o output/$outdir/lgm -c short -w 24 -j parallel -n 32 \&control="nyears=10000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42  fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc co2_const=190 ch4_const=375 n2o_const=200" lnd_par="lithology_uhh_file=input/Lithology_lgm_UHH.nc"
#./job_climber -s -f -o output/bench_v1.3/lgm_deglac -c short -w 24 -j parallel -n 32 \&control="nyears=10000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42  fake_geo_const_file=input/geo_ice_tarasov_deglac_21ka.nc fake_geo_ref_file=input/geo_ice_tarasov_deglac_0ka.nc fake_ice_const_file=input/geo_ice_tarasov_deglac_21ka.nc co2_const=190 ch4_const=375 n2o_const=200" lnd_par="lithology_uhh_file=input/Lithology_lgm_UHH.nc"
./job_climber -s -f -o output/$outdir/lgm_fixveg -c short -w 24 -j parallel -n 32 \&control="nyears=10000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42  fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc co2_const=190 ch4_const=375 n2o_const=200" lnd_par="l_dynveg=F l_co2_fert_lim=T l_fixlai=T"
./job_climber -s -f -o output/$outdir/lgm_peltier -c short -w 24 -j parallel -n 32 \&control="nyears=10000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42  fake_geo_const_file=input/geo_ice_peltier_lgm.nc fake_geo_ref_file=input/geo_ice_peltier_deglac_0ka.nc fake_ice_const_file=input/geo_ice_peltier_lgm.nc co2_const=190 ch4_const=375 n2o_const=200" 

# hysteresis
./job_climber -s -f -o output/$outdir/hyst2050/up   -c medium -w 60 -j parallel -n 32 \&control="nyears=30000" ocn_par="l_hosing=T hosing_ini=-0.3 hosing_trend=0.02 lat_min_hosing=20 lat_max_hosing=50"
./job_climber -s -f -o output/$outdir/hyst2050/down -c medium -w 60 -j parallel -n 32 \&control="nyears=30000" ocn_par="l_hosing=T hosing_ini=0.3 hosing_trend=-0.02 lat_min_hosing=20 lat_max_hosing=50"
./job_climber -s -f -o output/$outdir/hyst5070/up   -c medium -w 60 -j parallel -n 32 \&control="nyears=30000" ocn_par="l_hosing=T hosing_ini=-0.3 hosing_trend=0.02 lat_min_hosing=50 lat_max_hosing=70"
./job_climber -s -f -o output/$outdir/hyst5070/down -c medium -w 60 -j parallel -n 32 \&control="nyears=30000" ocn_par="l_hosing=T hosing_ini=0.3 hosing_trend=-0.02 lat_min_hosing=50 lat_max_hosing=70"

# PaleoDEM
./job_climber -s -f -o output/$outdir/pliomip -c short -w 24 -j parallel -n 32 \&control="nyears=5000 fake_geo_const_file=input/geo_Pliomip2.nc fake_ice_const_file=input/geo_Pliomip2.nc co2_const=400" geo_par="geo_ref_file=input/geo_Pliomip2.nc"
./job_climber -s -f -o output/$outdir/65Ma -c short -w 24 -j parallel -n 32 \&control="nyears=5000 fake_geo_const_file=input/topog_065Ma_climberX.nc fake_ice_const_file=input/topog_065Ma_climberX.nc atm_restart=F ocn_restart=F sic_restart=F lnd_restart=F dt_day_ocn=1" ocn_par="i_isl=0 i_init=1" geo_par="l_close_panama=F geo_ref_file=input/topog_065Ma_climberX.nc" 
./job_climber -s -f -o output/$outdir/250Ma -c short -w 24 -j parallel -n 32 \&control="nyears=5000 fake_geo_const_file=input/topog_250Ma_climberX.nc fake_ice_const_file=input/topog_250Ma_climberX.nc atm_restart=F ocn_restart=F sic_restart=F lnd_restart=F dt_day_ocn=1" ocn_par="i_isl=0 i_init=1" geo_par="l_close_panama=F geo_ref_file=input/topog_250Ma_climberX.nc" 
fi

####################
# step 3
if [ $step -eq 3 ]
then

# historical 
./job_climber -s -f -o output/$outdir/hist       -c short -w 10 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=1 iluc=2 isol=1 ivolc=1 flag_ch4=T ich4_rad=2 ich4_emis=1 id13c=1 iD14c=1 year_out_start=-19 nyout_atm=1 nyout_ocn=1 nyout_sic=1 nyout_lnd=1 nyout_bgc=1" ocn_par="l_cfc=T" "${cc_string_bgc}" "${cc_string_lnd}" 
./job_climber -s -f -o output/$outdir/hist-nat   -c short -w 6 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T iorbit=1 ico2=0 ich4=0 in2o=0 icfc=0 io3=1 iso4=0 iluc=1 isol=1 ivolc=1" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/hist-ghg   -c short -w 6 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=0 iluc=1 isol=0 ivolc=0" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/hist-aer   -c short -w 6 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T iorbit=1 ico2=0 ich4=0 in2o=0 icfc=0 io3=1 iso4=1 iluc=1 isol=0 ivolc=0" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/hist-noluc -c short -w 6 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=1 iluc=1 isol=1 ivolc=1" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/hist-esm   -c short -w 6 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T flag_co2=T ico2_emis=1 id13C_emis=1 iorbit=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=1 iluc=2 isol=1 ivolc=1" "${cc_string_bgc}" "${cc_string_lnd}"

# historical with SEMI
./job_climber -s -f -o output/$outdir/hist_smb_grl -c short -w 10 -j parallel -n 32 \&control="nyears=215 year_ini=-200 iorbit=1 ico2=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=1 iluc=2 isol=1 ivolc=1 year_out_start=-19 n_year_smb=1 flag_smb=T nyout_smb=1 ice_domain_name=GRL-8KM" smb_par="l_monthly_output=T i_alb_ice=0"

# SSP scenarios
./job_climber -s -f -o output/$outdir/ssp_conc/ssp119 -c short -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp119.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp119.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp119.nc icfc=1 cfc_file=input/CFCs_historical_ssp119.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp126 -c short -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp126.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp126.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp126.nc icfc=1 cfc_file=input/CFCs_historical_ssp126.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp245 -c short -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp245.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp245.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp245.nc icfc=1 cfc_file=input/CFCs_historical_ssp245.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp370 -c short -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp370.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp370.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp370.nc icfc=1 cfc_file=input/CFCs_historical_ssp370.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp460 -c short -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp460.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp460.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp460.nc icfc=1 cfc_file=input/CFCs_historical_ssp460.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp585 -c short -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp585.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp585.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp585.nc icfc=1 cfc_file=input/CFCs_historical_ssp585.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"

# SSP scenarios, emission driven
./job_climber -s -f -o output/$outdir/ssp_emis/ssp119 -c short -w 24 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp119.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp119_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp119.nc icfc=1 cfc_file=input/CFCs_historical_ssp119.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc isol=1 ivolc=1  nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis/ssp126 -c short -w 24 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp126.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp126_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp126.nc icfc=1 cfc_file=input/CFCs_historical_ssp126.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc isol=1 ivolc=1  nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis/ssp245 -c short -w 24 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp245.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp245_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp245.nc icfc=1 cfc_file=input/CFCs_historical_ssp245.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc isol=1 ivolc=1  nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis/ssp370 -c short -w 24 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp370.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp370_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp370.nc icfc=1 cfc_file=input/CFCs_historical_ssp370.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc isol=1 ivolc=1  nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis/ssp460 -c short -w 24 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp460.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp460_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp460.nc icfc=1 cfc_file=input/CFCs_historical_ssp460.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc isol=1 ivolc=1  nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis/ssp585 -c short -w 24 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp585.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp585_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp585.nc icfc=1 cfc_file=input/CFCs_historical_ssp585.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc isol=1 ivolc=1  nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100" "${cc_string_bgc}" "${cc_string_lnd}"

# 1%/yr CO2 increase
./job_climber -s -f -o output/$outdir/1pctCO2/cpl -c short -w 1 -j parallel -n 32 \&control="nyears=140 ico2=2 ico2_rad=0 co2_const=280 flag_bgc=T nyout_atm=140 nyout_ocn=140 nyout_sic=140 nyout_lnd=140" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/1pctCO2/cpl_out1 -c short -w 10 -j parallel -n 32 \&control="nyears=140 ico2=2 ico2_rad=0 co2_const=280 flag_bgc=T nyout_atm=140 nyout_ocn=1 nyout_sic=140 nyout_lnd=140" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/1pctCO2/bgc -c short -w 1 -j parallel -n 32 \&control="nyears=140 ico2=2 ico2_rad=1 co2_const=280 flag_bgc=T nyout_atm=140 nyout_ocn=140 nyout_sic=140 nyout_lnd=140" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/1pctCO2/rad -c short -w 1 -j parallel -n 32 \&control="nyears=140 ico2=0 ico2_rad=3 co2_const=280 flag_bgc=T nyout_atm=140 nyout_ocn=140 nyout_sic=140 nyout_lnd=140" "${cc_string_bgc}" "${cc_string_lnd}"

# ZECMIP
./job_climber -s -f -o output/$outdir/zecmip -c short -w 24 -j parallel -n 32 \&control="nyears=10000 flag_co2=T ico2_emis=4 co2_pulse=1000 flag_bgc=T"

# preindustrial with interactive CO2, to check for drift
./job_climber -s -f -o output/$outdir/preind_co2_close -c short -w 24 -j parallel -n 32 \&control="nyears=10000 flag_co2=T flag_bgc=T" "${cc_string_close_bgc}" "${cc_string_close_lnd}"

# LGM with prognostic CO2
./job_climber -s -f -o output/$outdir/lgm_co2_close_Cicebur -c short -w 24 -j parallel -n 32 \&control="nyears=10000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 flag_co2=T flag_bgc=T fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ico2_rad=1 co2_rad_const=190 ch4_const=375 n2o_const=200" lnd_par="k_ice=10" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
./job_climber -s -f -o output/$outdir/lgm_co2_close_Cicenobur -c short -w 24 -j parallel -n 32 \&control="nyears=10000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 flag_co2=T flag_bgc=T fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ico2_rad=1 co2_rad_const=190 ch4_const=375 n2o_const=200" lnd_par="k_ice=1e6" "${cc_string_close_bgc}" "${cc_string_close_lnd}"

# deglaciation
#./job_climber -s -f -o output/$outdir/deglac_tarasov -c medium -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 ifake_ice=1 ifake_geo=1 ico2=1 ich4=1 in2o=1 iorbit=2 fake_geo_var_file=input/geo_ice_tarasov_lgc.nc fake_ice_var_file=input/geo_ice_tarasov_lgc.nc n_year_geo=100 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_geo=100 restart_in_dir=restart_lgm"
#./job_climber -s -f -o output/$outdir/deglac_peltier -c medium -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 ifake_ice=1 ifake_geo=1 ico2=1 ich4=1 in2o=1 iorbit=2 fake_geo_var_file=input/geo_ice_peltier_deglac.nc fake_geo_ref_file=input/geo_ice_peltier_deglac_0ka.nc fake_ice_var_file=input/geo_ice_peltier_deglac.nc n_year_geo=100 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_geo=100 restart_in_dir=restart_lgm_peltier"
#./job_climber -s -f -o output/$outdir/deglac_gowan -c medium -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 ifake_ice=1 ifake_geo=1 ico2=1 ich4=1 in2o=1 iorbit=2 fake_geo_var_file=input/geo_ice_gowan.nc fake_geo_ref_file=input/geo_ice_gowan_0ka.nc fake_ice_var_file=input/geo_ice_gowan.nc n_year_geo=100 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_geo=100 restart_in_dir=restart_lgm"

# last glacial cycle with prescribed ice sheets
./job_climber -s -f -o output/$outdir/lgc_tarasov -c long -w 240 -j parallel -n 32 \&control="nyears=125000 year_ini=-125000 ifake_ice=1 ifake_geo=1 ico2=1 ich4=1 in2o=1 iorbit=2 fake_geo_var_file=input/geo_ice_tarasov_lgc.nc fake_ice_var_file=input/geo_ice_tarasov_lgc.nc n_year_geo=100" #ocn_par="l_noise_fw=T noise_amp_fw=0,0.5 scale_runoff_ice=0,1"

# last glacial cycle with interactive ice sheets
#./job_climber -s -f -o output/$outdir/lgc_ice_sico -c long -w 500 -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=1 ice_model_name=sico" 
#./job_climber -s -f -o output/$outdir/lgc_ice_sico_acc2 -c long -w 500 -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=2 ice_model_name=sico" 
#./job_climber -s -f -o output/$outdir/lgc_ice_sico_acc5 -c medium -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=5 ice_model_name=sico" 
#./job_climber -s -f -o output/$outdir/lgc_ice_sico_acc10 -c medium -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=10 ice_model_name=sico" 
#./job_climber -s -f -o output/$outdir/lgc_ice_yelmo -c long -w 500 -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=1 ice_model_name=yelmo" 
#./job_climber -s -f -o output/$outdir/lgc_ice_yelmo_acc2 -c long -w 500 -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=2 ice_model_name=yelmo" 
#./job_climber -s -f -o output/$outdir/lgc_ice_yelmo_acc5 -c medium -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=5 ice_model_name=yelmo" 
#./job_climber -s -f -o output/$outdir/lgc_ice_yelmo_acc10 -c medium -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=10 ice_model_name=yelmo" 

#
#./job_climber -s -f -o output/$outdir/deglac_ice_acc10_itempinit1 -c standby -j parallel -n 32 \&control="year_ini=-30000 nyears=30000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=10 ifake_ice=2 ifake_geo=2" ice_sico_par="i_temp_init=1"
#./job_climber -s -f -o output/$outdir/deglac_ice_acc10_itempinit4 -c standby -j parallel -n 32 \&control="year_ini=-30000 nyears=30000 iorbit=2 ico2=1 ich4=1 in2o=1 ice_domain_name=NH-32KM flag_geo=T flag_ice=T flag_smb=T flag_bmb=T n_accel=10 ifake_ice=2 ifake_geo=2" ice_sico_par="i_temp_init=4"
#
## preind with interactive ice sheets
#./job_climber -s -f -o output/$outdir/preind_ice_nh -c standby -j parallel -n 32 \&control="nyears=100000 n_accel=10 n_year_smb=10 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM"
## MIS11 with interactive ice sheets
#./job_climber -s -f -o output/$outdir/mis11_ice_nh -c standby -j parallel -n 32 \&control="year_ini=-398000 nyears=100000 n_accel=10 n_year_smb=10 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM"

fi

####################
# step 5
if [ $step -eq 5 ]
then

# preind with prognostic CO2 and open carbon cycle, to check for drift
./job_climber -s -f -o output/$outdir/preind_co2_open -c standby -w 48 -j parallel -n 32 \&control="nyears=10000 flag_co2=T flag_bgc=T ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# historical 
./job_climber -s -f -o output/$outdir/hist_cc_open       -c standby -w 10 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T iorbit=1 ico2=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=1 iluc=2 isol=1 ivolc=1 flag_ch4=T ich4_rad=2 ich4_emis=1 id13c=1 iD14c=1 year_out_start=-19 nyout_atm=1 nyout_ocn=1 nyout_sic=1 nyout_lnd=1 nyout_bgc=1 restart_in_dir=restart_pi_cc_open l_weathering=T" ocn_par="l_cfc=T" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 
./job_climber -s -f -o output/$outdir/hist-esm_cc_open   -c standby -w 10 -j parallel -n 32 \&control="nyears=1015 year_ini=-1000 flag_bgc=T flag_co2=T ico2_emis=1 id13C_emis=1 iorbit=1 ich4=1 in2o=1 icfc=1 io3=2 iso4=1 iluc=2 isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# ZECMIP
./job_climber -s -f -o output/$outdir/zecmip_cc_open -c standby -w 48 -j parallel -n 32 \&control="nyears=10000 flag_co2=T ico2_emis=4 co2_pulse=1000 flag_bgc=T ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# SSP scenarios, emission driven
./job_climber -s -f -o output/$outdir/ssp_emis_cc_open/ssp119 -c standby -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp119.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp119_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp119.nc icfc=1 cfc_file=input/CFCs_historical_ssp119.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis_cc_open/ssp126 -c standby -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp126.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp126_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp126.nc icfc=1 cfc_file=input/CFCs_historical_ssp126.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis_cc_open/ssp245 -c standby -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp245.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp245_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp245.nc icfc=1 cfc_file=input/CFCs_historical_ssp245.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis_cc_open/ssp370 -c standby -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp370.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp370_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp370.nc icfc=1 cfc_file=input/CFCs_historical_ssp370.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis_cc_open/ssp460 -c standby -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp460.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp460_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp460.nc icfc=1 cfc_file=input/CFCs_historical_ssp460.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
./job_climber -s -f -o output/$outdir/ssp_emis_cc_open/ssp585 -c standby -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp585.nc flag_ch4=T ich4_emis=1 ch4_emis_file=input/ch4_emis_historical_ssp585_ext.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp585.nc icfc=1 cfc_file=input/CFCs_historical_ssp585.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# SSP scenarios, emission driven
#./job_climber -s -f -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp119 -c medium -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp119.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp119.nc  in2o=1 n2o_file=input/n2o_Koehler2017_ssp119.nc icfc=1 cfc_file=input/CFCs_historical_ssp119.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
#./job_climber -s -f -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp126 -c medium -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp126.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp126.nc  in2o=1 n2o_file=input/n2o_Koehler2017_ssp126.nc icfc=1 cfc_file=input/CFCs_historical_ssp126.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
#./job_climber -s -f -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp245 -c medium -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp245.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp245.nc  in2o=1 n2o_file=input/n2o_Koehler2017_ssp245.nc icfc=1 cfc_file=input/CFCs_historical_ssp245.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
#./job_climber -s -f -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp370 -c medium -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp370.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp370.nc  in2o=1 n2o_file=input/n2o_Koehler2017_ssp370.nc icfc=1 cfc_file=input/CFCs_historical_ssp370.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
#./job_climber -s -f -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp460 -c medium -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp460.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp460.nc  in2o=1 n2o_file=input/n2o_Koehler2017_ssp460.nc icfc=1 cfc_file=input/CFCs_historical_ssp460.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
#./job_climber -s -f -o output/$outdir/ssp_emis_cc_open_sspCH4/ssp585 -c medium -w 48 -j parallel -n 32 \&control="nyears=6000 year_ini=-1000 flag_bgc=T iorbit=1 flag_co2=T ico2_emis=1 id13C_emis=1 co2_emis_file=input/co2_emis_hist_ssp585.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp585.nc  in2o=1 n2o_file=input/n2o_Koehler2017_ssp585.nc icfc=1 cfc_file=input/CFCs_historical_ssp585.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc isol=1 ivolc=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

fi

####################
# step 6
if [ $step -eq 6 ]
then

co2_volc_lgc=0.05617464

# control to check for drift
# close carbon cycle
./job_climber -s -f -o output/$outdir/lig125k_co2_close_control -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=25000 co2_const=278 flag_co2=T flag_bgc=T ch4_const=640 n2o_const=260 restart_in_dir=restart_lig_cc_close l_c14=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
# open carbon cycle
./job_climber -s -f -o output/$outdir/lig125k_co2_open_volLIG_control -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=25000 co2_const=278 flag_co2=T flag_bgc=T ch4_const=640 n2o_const=260 ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
./job_climber -s -f -o output/$outdir/lig125k_co2_open_volLGC_control -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=25000 co2_const=278 flag_co2=T flag_bgc=T ch4_const=640 n2o_const=260 ico2_degas=0 co2_degas_const=$co2_volc_lgc l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# inception with interactive ice sheets and one-way coupled CO2 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=50000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close" "${cc_string_close_bgc}" "${cc_string_close_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=50000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLIG -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=50000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLIG_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=50000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_lnd_co2=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLGC -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=50000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=0 co2_degas_const=$co2_volc_lgc l_weathering=T restart_in_dir=restart_lig_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini125k_ice_co2_rad_open_volLGC_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=50000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=0 co2_degas_const=$co2_volc_lgc l_weathering=T restart_in_dir=restart_lig_cc_open l_lnd_co2=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 

./job_climber -s -f -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_close -c standby -j parallel -n 32 \&control="year_ini=-120000 nyears=50000 iorbit=2 co2_const=280 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_pi" "${cc_string_close_bgc}" "${cc_string_close_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-120000 nyears=50000 iorbit=2 co2_const=280 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_pi l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volPI -c standby -j parallel -n 32 \&control="year_ini=-120000 nyears=50000 iorbit=2 co2_const=280 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volPI_nolnd -c standby -j parallel -n 32 \&control="year_ini=-120000 nyears=50000 iorbit=2 co2_const=280 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open l_lnd_co2=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volLGC -c standby -j parallel -n 32 \&control="year_ini=-120000 nyears=50000 iorbit=2 co2_const=280 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=0 co2_degas_const=$co2_volc_lgc l_weathering=T restart_in_dir=restart_pi_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 
./job_climber -s -f -o output/$outdir/incep_co2/incep_ini120kPI_ice_co2_rad_open_volLGC_nolnd -c standby -j parallel -n 32 \&control="year_ini=-120000 nyears=50000 iorbit=2 co2_const=280 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=0 co2_degas_const=$co2_volc_lgc l_weathering=T restart_in_dir=restart_pi_cc_open l_lnd_co2=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}" 

# inception with fully interactive ice sheets and CO2 
#./job_climber -s -f -o output/$outdir/inception_co2/incep_ice_co2_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close l_c14=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
#./job_climber -s -f -o output/$outdir/inception_co2/incep_ice_co2_open -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

# inception with fixed ice sheets and one-way coupled CO2 
#./job_climber -s -f -o output/$outdir/inception_co2/incep_co2_rad_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
#./job_climber -s -f -o output/$outdir/inception_co2/incep_co2_rad_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
#./job_climber -s -f -o output/$outdir/inception_co2/incep_co2_rad_open -c medium -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
#./job_climber -s -f -o output/$outdir/inception_co2/incep_co2_rad_open_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_lnd_co2=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

## inception with fixed ice sheets and interactive CO2 
#./job_climber -s -f -o output/$outdir/inception_co2/incep_co2_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close l_c14=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
#./job_climber -s -f -o output/$outdir/inception_co2/incep_co2_open -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"


## lgm with prognostic CO2 and open carbon cycle
#./job_climber -s -f -o output/$outdir/lgm_co2_open -c standby -w 72 -j parallel -n 32 \&control="nyears=15000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 flag_co2=T flag_bgc=T fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ico2_rad=1 co2_rad_const=190 ch4_const=375 n2o_const=200 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" lnd_par="k_ice=10" "${cc_string_open_bgc}" "${cc_string_open_lnd}" lnd_par="lithology_uhh_file=input/Lithology_lgm_UHH.nc"

##
##### deglaciation with interactive CO2
###./job_climber -s -f -o output/$outdir/deglac_co2_close -c medium -w 168 -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 ifake_ice=1 ifake_geo=1 flag_co2=T flag_bgc=T co2_restart=T ico2_rad=2 ich4=1 in2o=1 iorbit=2 restart_in_dir=restart_lgm_cc_close" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
###./job_climber -s -f -o output/$outdir/deglac_co2_close_control -c medium -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 flag_co2=T flag_bgc=T co2_restart=T ico2_rad=1 co2_rad_const=190 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ch4_const=375 n2o_const=200 restart_in_dir=restart_lgm_cc_close" "${cc_string_close_bgc}" "${cc_string_close_lnd}"
####./job_climber -s -f -o output/$outdir/deglac_co2_open -c medium -w 168 -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 ifake_ice=1 ifake_geo=1 flag_co2=T flag_bgc=T co2_restart=T ico2_rad=2 ich4=1 in2o=1 iorbit=2ico2_degas=1 restart_in_dir=restart_lgm_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
####./job_climber -s -f -o output/$outdir/deglac_co2_open_control -c medium -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 flag_co2=T flag_bgc=T co2_restart=T ico2_rad=1 co2_rad_const=190 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ch4_const=375 n2o_const=200 ico2_degas=1 l_weathering=T restart_in_dir=restart_lgm_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
####./job_climber -s -f -o output/$outdir/deglac_co2_open2 -c medium -w 168 -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 ifake_ice=1 ifake_geo=1 flag_co2=T flag_bgc=T co2_restart=T ico2_rad=2 ich4=1 in2o=1 iorbit=2ico2_degas=1 restart_in_dir=restart_lgm_cc_open2" "${cc_string_open2_bgc}" "${cc_string_open2_lnd}"
####./job_climber -s -f -o output/$outdir/deglac_co2_open2_control -c medium -j parallel -n 32 \&control="nyears=21000 year_ini=-21000 flag_co2=T flag_bgc=T co2_restart=T ico2_rad=1 co2_rad_const=190 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ch4_const=375 n2o_const=200 ico2_degas=1 l_weathering=T restart_in_dir=restart_lgm_cc_open2" "${cc_string_open2_bgc}" "${cc_string_open2_lnd}"
##
fi

####################
# step 7
if [ $step -eq 7 ]
then

./job_climber -s -f -o output/$outdir/incep_ice_co2_rad_open_seddiss -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open" "${cc_string_open_bgc}" "${cc_string_open_lnd}" bgc_par="l_sed_dry_dissolve=T"

#./job_climber -s -f -o output/$outdir/incep_ice_co2_open_dT2 -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}" smb_par="t2m_bias_corr_uniform=-2"
#./job_climber -s -f -o output/$outdir/incep_ice_co2_open_fw01 -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}" ocn_par="l_hosing=T hosing_ini=0.1 year_hosing_ini=6000"
#./job_climber -s -f -o output/$outdir/incep_ice_co2_open_fw02 -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_c14=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}" ocn_par="l_hosing=T hosing_ini=0.2 year_hosing_ini=6000"

#./job_climber -s -f -o output/$outdir/incep_ice_co2_rad_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}"

#./job_climber -s -f -o output/$outdir/incep_co2_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close l_c14=F l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}"

#./job_climber -s -f -o output/$outdir/incep_fixRrain_ice_co2_rad_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close_fixRrain" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2"
#./job_climber -s -f -o output/$outdir/incep_fixRrain_ice_co2_rad_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close_fixRrain l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2"
#
#./job_climber -s -f -o output/$outdir/incep_fixRrain_co2_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close_fixRrain l_c14=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2"
#./job_climber -s -f -o output/$outdir/incep_fixRrain_co2_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close_fixRrain l_c14=F l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2"

#./job_climber -s -f -o output/$outdir/incep_fixRrain_corals_ice_co2_rad_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close_fixRrain_corals" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2 l_corals=T"
#./job_climber -s -f -o output/$outdir/incep_fixRrain_corals_ice_co2_rad_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM restart_in_dir=restart_lig_cc_close_fixRrain_corals l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2 l_corals=T"
#
#./job_climber -s -f -o output/$outdir/incep_fixRrain_corals_co2_close -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close_fixRrain_corals l_c14=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2 l_corals=T"
#./job_climber -s -f -o output/$outdir/incep_fixRrain_corals_co2_close_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T flag_bgc=T ich4=1 in2o=1 restart_in_dir=restart_lig_cc_close_fixRrain_corals l_c14=F l_lnd_co2=F" "${cc_string_close_bgc}" "${cc_string_close_lnd}" bgc_par="calmax=0.07 i_delcar=2 l_corals=T"

#./job_climber -s -f -o output/$outdir/incep_ice_co2_rad_open_nolnd -c standby -j parallel -n 32 \&control="year_ini=-125000 nyears=125000 iorbit=2 co2_const=278 flag_co2=T ico2_rad=2 flag_bgc=T ich4=1 in2o=1 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T ice_domain_name=NH-32KM ico2_degas=1 l_weathering=T restart_in_dir=restart_lig_cc_open l_lnd_co2=F" "${cc_string_open_bgc}" "${cc_string_open_lnd}"

fi

####################
# step 8
if [ $step -eq 8 ]
then

./job_climber -s -f -o output/$outdir/lgm_co2_open2 -c standby -w 72 -j parallel -n 32 \&control="nyears=15000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 flag_co2=T flag_bgc=T fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ico2_rad=1 co2_rad_const=190 ch4_const=375 n2o_const=200 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" lnd_par="k_ice=10" "${cc_string_open_bgc}" "${cc_string_open_lnd}"
#./job_climber -s -f -o output/$outdir/lgm_co2_open -c standby -w 72 -j parallel -n 32 \&control="nyears=15000 iorbit=1 ecc_const=0.018994 obl_const=22.949 per_const=114.42 flag_co2=T flag_bgc=T fake_geo_const_file=input/geo_ice_tarasov_lgm.nc fake_ice_const_file=input/geo_ice_tarasov_lgm.nc ico2_rad=1 co2_rad_const=190 ch4_const=375 n2o_const=200 ico2_degas=1 l_weathering=T restart_in_dir=restart_pi_cc_open" lnd_par="k_ice=10" "${cc_string_open_bgc}" "${cc_string_open_lnd}" lnd_par="lithology_uhh_file=input/Lithology_lgm_UHH.nc"

fi

####################
# step 9
if [ $step -eq 9 ]
then

drag_topo_fac=3.5
restart_in_dir=restart_pi_v22

# SSP scenarios
./job_climber -s -f -o output/$outdir/ssp_conc/ssp119 -c standby -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp119.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp119.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp119.nc icfc=1 cfc_file=input/CFCs_historical_ssp119.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp119_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP19_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp126 -c standby -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp126.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp126.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp126.nc icfc=1 cfc_file=input/CFCs_historical_ssp126.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp126_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP1_RCP26_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp245 -c standby -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp245.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp245.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp245.nc icfc=1 cfc_file=input/CFCs_historical_ssp245.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp245_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP2_RCP45_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp370 -c standby -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp370.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp370.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp370.nc icfc=1 cfc_file=input/CFCs_historical_ssp370.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp370_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP3_RCP70_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp460 -c standby -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp460.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp460.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp460.nc icfc=1 cfc_file=input/CFCs_historical_ssp460.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp460_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP4_RCP60_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac" "${cc_string_bgc}" "${cc_string_lnd}"
./job_climber -s -f -o output/$outdir/ssp_conc/ssp585 -c standby -w 20 -j parallel -n 32 \&control="nyears=3000 year_ini=-1000 iorbit=1 ico2=1 co2_file=input/co2_Koehler2017_ssp585.nc ich4=1 ch4_file=input/ch4_Koehler2017_ssp585.nc in2o=1 n2o_file=input/n2o_Koehler2017_ssp585.nc icfc=1 cfc_file=input/CFCs_historical_ssp585.nc io3=2 iso4=1 so4_file=input/so4_historical_ssp585_ext_anom_CMIP6_ensmean_ann_5x5.nc iluc=2 luc_file=input/LUH2_historical_SSP5_RCP85_850_2100_5x5.nc isol=1 ivolc=1 nyout_atm=100 nyout_ocn=100 nyout_sic=100 nyout_lnd=100 nyout_bgc=100 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac" "${cc_string_bgc}" "${cc_string_lnd}"

./job_climber -s -f -o output/$outdir/pliomip_400 -c standby -w 24 -j parallel -n 32 \&control="year_ini=0 fake_geo_const_file=input/geo_Pliomip2.nc fake_ice_const_file=input/geo_Pliomip2.nc co2_const=400 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac"
./job_climber -s -f -o output/$outdir/pliomip_500 -c standby -w 24 -j parallel -n 32 \&control="year_ini=0 fake_geo_const_file=input/geo_Pliomip2.nc fake_ice_const_file=input/geo_Pliomip2.nc co2_const=500 restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac"

./job_climber -s -f -a output/$outdir/DOevents/ -c standby -j parallel -n 32 \&control="year_ini=-70000 nyears=50000 iorbit=2 ice_domain_name=NH-32KM co2_const=180,190,200,210,220,230,240 flag_geo=T flag_ice=T flag_smb=T flag_bmb=T restart_in_dir=$restart_in_dir" ocn_par="drag_topo_fac=$drag_topo_fac"

fi


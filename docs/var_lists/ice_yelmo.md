| Variable | Dimensions | Units | Long Name |
|-----------|------------------------------|--------|------------|
| age_iso | age_iso | kyr | N/A |
| month | month | month | N/A |
| pc_steps | pc_steps | 1 | N/A |
| pd_age_iso | pd_age_iso | kyr | N/A |
| time | time | years | N/A |
| xc | xc | km | N/A |
| yc | yc | km | N/A |
| zeta | zeta | 1 | N/A |
| zeta_ac | zeta_ac | 1 | N/A |
| zeta_rock | zeta_rock | 1 | N/A |
| area | xc, yc | km^2 | N/A |
| basins | xc, yc | (0 - 8) | Hydrological basins |
| beta | xc, yc, time | Pa a m^-1 | Basal friction coefficient |
| bmb | xc, yc, time | m/a ice equiv. | Basal mass balance |
| bmb_shlf | xc, yc, time | m/a ice equiv. | Shelf basal mass balance |
| cmb | xc, yc, time | m/a ice equiv. | Calving mass balance |
| dHidt | xc, yc, time | m/a | Ice thickness change |
| dist_grline | xc, yc, time | m | Distance to grounding line |
| dmb | xc, yc, time | m/a ice equiv. | Discharge mass balance (subgrid) |
| dt_avg | time | yr | Average timestep |
| enh_bar | xc, yc, time | 1 | Vertically averaged enhancement factor |
| eta_avg | time | m a**-1 | Average eta (maximum PC truncation error) |
| f_grnd | xc, yc, time | 1 | Grounded ice fraction |
| f_ice | xc, yc, time | 1 | Total ice fraction |
| f_pmp | xc, yc, time | 1 | Fraction of grid point at pmp |
| H_calv | xc, yc, time | m | Threshold ice thickness for calving |
| H_ice | xc, yc, time | m | Ice thickness |
| H_sed | xc, yc | m | Stdev(z_bed) |
| H_w | xc, yc, time | m water equiv. | Basal water layer thickness |
| ice_allowed | xc, yc | N/A | Ice-allowed mask |
| mask_bed | xc, yc, time | N/A | Bed mask |
| mb_net | xc, yc, time | m | Net mass balance |
| mb_resid | xc, yc, time | m | Residual mass balance |
| N_eff | xc, yc, time | Pa | Effective pressure |
| pc_tau_max | xc, yc, time | m a**-1 | Maximum truncation error over last N timestep (magnitude) |
| Q_b | xc, yc, time | mW m^-2 | Basal frictional heating |
| Q_geo | xc, yc | mW/m^2 | Geothermal heat flux |
| Q_ice_b | xc, yc, time | mW m^-2 | Basal ice heat flow |
| Q_rock | xc, yc, time | mW m^-2 | Bedrock surface heat flow |
| regions | xc, yc | (0 - 8) | Domain regions |
| smb | xc, yc, time | m/a ice equiv. | Surface mass balance |
| speed | time | kyr/hr | Model speed (Yelmo only) |
| ssa_iter_avg | time | N/A | Average Picard iterations for SSA convergence |
| T_ice | xc, yc, zeta, time | K | Ice temperature |
| T_prime_b | xc, yc, time | K | Basal homologous ice temperature |
| taub | xc, yc, time | Pa | Basal dragging stress (magnitude) |
| taud | xc, yc, time | Pa | Driving stress |
| uxy_b | xc, yc, time | m/a | Basal velocity (magnitude) |
| uxy_s | xc, yc, time | m/a | Surface velocity (magnitude) |
| visc_eff_int | xc, yc, time | Pa a m | Depth-integrated effective viscosity (SSA) |
| x2D | xc, yc | km | N/A |
| y2D | xc, yc | km | N/A |
| z_bed | xc, yc, time | m | Bedrock elevation |
| z_bed_sd | xc, yc | m | Stdev(z_bed) |
| z_srf | xc, yc, time | m | Surface elevation |

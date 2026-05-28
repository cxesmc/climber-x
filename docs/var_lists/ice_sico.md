| Variable | Dimensions | Units | Long Name |
|-----------|------------------------------|--------|------------|
| kc | kc | 1 | N/A |
| kr | kr | 1 | N/A |
| kt | kt | 1 | N/A |
| time | time | years BP | N/A |
| x | x | km | N/A |
| y | y | km | N/A |
| accum | x, y, time | m/yr | accumulation flux |
| accum_apl | x, y, time | m/yr | applied accumulation flux |
| am_perp | x, y, time | m/yr | Volume flux across the z=zm interface |
| as_perp | x, y, time | m/yr | surface mass balance |
| as_perp_apl | x, y, time | m/yr | applied surface mass balance |
| beta_drag | x, y, time | / | basal drag parameter for shelfy stream |
| c_fric | x, y, time | / | basal friction coefficient |
| c_slide | x, y, time | m/(a*Pa^(p-q)) | basal sliding coefficient |
| calving | x, y, time | m/yr | calving flux |
| calving_apl | x, y, time | m/yr | applied calving flux |
| cos_grad_tc | x, y, time | / | Cosine of angle between surface gradient and cst dist gradient |
| cst_dist | x, y, time | km | coastal distance |
| delta | x, y, time | / | Fraction of overburden pressure at saturation |
| dH_c_dt | x, y, time | m/yr | Rate of change of the thickness of the upper (kc) ice layer |
| dH_dt | x, y, time | m/yr | Rate of change of the ice thickness |
| dH_t_dt | x, y, time | m/yr | Rate of change of the thickness of the lower (kt) ice layer |
| dis_perp | x, y, time | m/yr | applied calving flux from discharge parameterisation |
| dzb_dt | x, y, time | m/yr | Rate of change of the topography of the ice base |
| dzl_dt | x, y, time | m/yr | Rate of change of the topography of the lithosphere surface |
| dzm_dt | x, y, time | m/yr | Rate of change of the topography of the z=zm interface |
| dzs_dt | x, y, time | m/yr | Rate of change of the topography of the free surface |
| eta | y | m | y-coordinate |
| f_sed | x, y, time | / | sediment fraction |
| flag_shelfy_stream | x, y, time | / | Shelfy stream flag |
| flag_shelfy_stream_x | x, y, time | / | Shelfy stream flag in x-direction, at (i+1/2,j) |
| flag_shelfy_stream_y | x, y, time | / | Shelfy stream flag in y-direction, at (i,j+1/2) |
| H | x, y, time | m | Thickness of ice |
| H_calv | x, y, time | m | Calvin Threshold |
| H_cold | x, y, time | m | Thickness of the cold ice layer |
| H_eff | x, y, time | m | Effective thickness at shelf ice front |
| H_flow | x, y, time | m | Thickness of ice after dynamics (without source term) |
| H_sed | x, y, time | m | sediment thickness |
| H_temp | x, y, time | m | thickness of the temperate ice layer |
| H_w | x, y, time | m | Thickness of the water column under the ice base |
| id_mask | x, y, time | / | Mask indicating ice IDs |
| kc_cts | x, y, time | / | Grid index of the CTS position |
| lambda | x, y | degE | geographic longitude |
| mask_ablation_type | x, y, time | / | Mask indicating ablation type |
| mask_mar | x, y, time | 1 | margina ring mask |
| maske | x, y, time | / | ice-land-sea mask |
| maske_old | x, y, time | / | ice-land-sea mask (old) |
| mb_source | x, y, time | m/yr | total mass balance |
| mb_source_apl | x, y, time | m/yr | applied total mass balance |
| n_cts | x, y, time | / | mask for polythermal conditions |
| p_b | x, y, time | Pa | Basal pressure |
| p_b_red_lim | x, y, time | Pa | Reduced basal pressure |
| p_b_w | x, y, time | Pa | Basal water pressure |
| phi | x, y | degN | geographic latitude |
| Q_b | x, y, time | m/yr | Basal melting rate + Water drainage from the temperate layer |
| Q_b_apl | x, y, time | m/yr | applied Basal melting rate + Water drainage from the temperate layer |
| Q_bm | x, y, time | m/yr | Basal melting rate |
| q_geo | x, y, time | W/m2 | geothermal heat flux |
| q_gl_g | x, y, time | m2/yr | Horizontal volume flux across the grounding line |
| Q_tld | x, y, time | m/yr | Water drainage from the temperate layer |
| qx | x, y, time | m2/yr | Horizontal volume flux qx |
| qy | x, y, time | m2/yr | Horizontal volume flux qy |
| ratio_sl_x | x, y, time | / | Ratio of basal to surface velocity (slip ratio) in x-direction at (i+1/2,j) |
| ratio_sl_y | x, y, time | / | Ratio of basal to surface velocity (slip ratio) in y-direction at (i+1/2,j)? |
| runoff | x, y, time | m/yr | runoff flux |
| runoff_apl | x, y, time | m/yr | applied runoff flux |
| sigma_level_c | kc | 1 | sigma level in cold ice |
| sigma_level_r | kr | 1 | sigma level in litosphere |
| sigma_level_t | kt | 1 | sigma level in temperate ice |
| tau_b_drag | x, y, time | Pa | Basal drag |
| tau_b_driving | x, y, time | Pa | Driving stress |
| temp_b | x, y, time | degC | Temperature at the ice base |
| temp_g | x, y, time | degC | ground temperature |
| temp_s | x, y, time | degC | surface temperature |
| temph_b | x, y, time | degC | Temperature at the ice base relative to the pressure melting point |
| vb_t | x, y, time | m/a | threshold basal velocity for regularized Coulomb law |
| vh_b | x, y, time | m/yr | Horizontal velocity vh at the ice base |
| vh_m | x, y, time | m/yr | Vertical mean of horizontal velocity vh |
| vh_m_sia | x, y, time | m/yr | Vertical mean of SIA horizontal velocity |
| vh_m_ssa | x, y, time | m/yr | Vertical mean of SSA horizontal velocity |
| vh_s | x, y, time | m/yr | Horizontal velocity vh at the ice surface |
| vis_int_g | x, y, time | Pa s m | Depth-integrated viscosity |
| vx_b_g | x, y, time | m/yr | Horizontal velocity vx at the ice base |
| vx_m_g | x, y, time | m/yr | Vertical mean of horizontal velocity vx |
| vx_m_sia | x, y, time | m/yr | x depth-averaged horizontal velocity from SIA |
| vx_m_ssa | x, y, time | m/yr | x depth-averaged horizontal velocity from SSA |
| vx_s_g | x, y, time | m/yr | Horizontal velocity vx at the ice surface |
| vy_b_g | x, y, time | m/yr | Horizontal velocity vy at the ice base |
| vy_m_g | x, y, time | m/yr | Vertical mean of horizontal velocity vy |
| vy_m_sia | x, y, time | m/yr | y depth-averaged horizontal velocity from SIA |
| vy_m_ssa | x, y, time | m/yr | y depth-averaged horizontal velocity from SSA |
| vy_s_g | x, y, time | m/yr | Horizontal velocity vy at the ice surface |
| vz_b | x, y, time | m/yr | Vertical velocity vz at the ice base |
| vz_s | x, y, time | m/yr | Vertical velocity vz at the ice surface |
| weigh_ssta_sia_x | x, y, time | / | weight of SStA vs SIA in x-direction at (i+1/2,j) |
| weigh_ssta_sia_y | x, y, time | / | weight of SStA vs SIA in y-direction at (i,j+1/2) |
| xi | x | m | x-coordinate |
| z_sl | x, y, time | m | sea level |
| zb | x, y, time | m | Topography of the ice base |
| zl | x, y, time | m | Topography of the lithosphere surface |
| zl0 | x, y, time | m | Topography of isostatically relaxed lithosphere surface |
| zl_fil | x, y, time | m | Topography of the filtered lithosphere surface |
| zl_std | x, y, time | m | sub-grid standard deviation of topography of the lithosphere surface |
| zm | x, y, time | m | Topography of z=zm interface |
| zs | x, y, time | m | Topography of the free surface |

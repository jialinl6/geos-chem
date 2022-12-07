import Doc_string_extensions

"""

"""
module state_diag_mod

mutable struct DgnState
  # Standard Simulation Diagnostic Arrays

  # Restart file fields

  species_rst::Array{AbstractFloat,4}
  archive_species_rst::Bool

  #  Boundary condition fields

  species_bc::Array{AbstractFloat,4}
  map_species_bc::DgnMap
  archive_species_bc::Bool

  #  Concentrations

  species_conc::Array{AbstractFloat,4}
  map_species_conc::DgnMap
  archive_species_conc::Bool

  #  ML diagnostics
  conc_before_chem::Array{AbstractFloat,4}
  map_conc_before_chem::DgnMap
  archive_conc_before_chem::Bool

  conc_after_chem::Array{AbstractFloat,4}
  map_conc_after_chem::DgnMap
  archive_conc_after_chem::Bool

  #ifdef ADJOINT
  # Adjoint variables for diagnostic output
  species_adj::Array{AbstractFloat,4}
  map_species_adj::DgnMap
  archive_species_adj::Bool

  # Concentrations
  scale_ics_adj::Array{AbstractFloat,4}
  map_scale_ics_adj::DgnMap
  archive_scale_ics_adj::Bool
  #endif

  # Budget diagnostics

  budget_emis_dry_dep_full::Array{AbstractFloat,3}
  map_budget_emis_dry_dep_full::DgnMap
  archive_budget_emis_dry_dep_full::Bool

  budget_emis_dry_dep_trop::Array{AbstractFloat,3}
  map_budget_emis_dry_dep_trop::DgnMap
  archive_budget_emis_dry_dep_trop::Bool

  budget_emis_dry_dep_pbl::Array{AbstractFloat,3}
  map_budget_emis_dry_dep_pbl::DgnMap
  archive_budget_emis_dry_dep_pbl::Bool

  budget_transport_full::Array{AbstractFloat,3}
  map_budget_transport_full::DgnMap
  archive_budget_transport_full::Bool

  budget_transport_trop::Array{AbstractFloat,3}
  map_budget_transport_trop::DgnMap
  archive_budget_transport_trop::Bool

  budget_transport_pbl::Array{AbstractFloat,3}
  map_budget_transport_pbl::DgnMap
  archive_budget_transport_pbl::Bool

  budget_mixing_full::Array{AbstractFloat,3}
  map_budget_mixing_full::DgnMap
  archive_budget_mixing_full::Bool

  budget_mixing_trop::Array{AbstractFloat,3}
  map_budget_mixing_trop::DgnMap
  archive_budget_mixing_trop::Bool

  budget_mixing_pbl::Array{AbstractFloat,3}
  map_budget_mixing_pbl::DgnMap
  archive_budget_mixing_pbl::Bool

  budget_convection_full::Array{AbstractFloat,3}
  map_budget_convection_full::DgnMap
  archive_budget_convection_full::Bool

  budget_convection_trop::Array{AbstractFloat,3}
  map_budget_convection_trop::DgnMap
  archive_budget_convection_trop::Bool

  budget_convection_pbl::Array{AbstractFloat,3}
  map_budget_convection_pbl::DgnMap
  archive_budget_convection_pbl::Bool

  budget_chemistry_full::Array{AbstractFloat,3}
  map_budget_chemistry_full::DgnMap
  archive_budget_chemistry_full::Bool

  budget_chemistry_trop::Array{AbstractFloat,3}
  map_budget_chemistry_trop::DgnMap
  archive_budget_chemistry_trop::Bool

  budget_chemistry_pbl::Array{AbstractFloat,3}
  map_budget_chemistry_pbl::DgnMap
  archive_budget_chemistry_pbl::Bool

  budget_wet_dep_full::Array{AbstractFloat,3}
  map_budget_wet_dep_full::DgnMap
  archive_budget_wet_dep_full::Bool

  budget_wet_dep_trop::Array{AbstractFloat,3}
  map_budget_wet_dep_trop::DgnMap
  archive_budget_wet_dep_trop::Bool

  budget_wet_dep_pbl::Array{AbstractFloat,3}
  map_budget_wet_dep_pbl::DgnMap
  archive_budget_wet_dep_pbl::Bool

  budget_column_mass::Array{AbstractFloat,4}
  archive_budget_emis_dry_dep::Bool
  archive_budget_transport::Bool
  archive_budget_mixing::Bool
  archive_budget_convection::Bool
  archive_budget_chemistry::Bool
  archive_budget_wet_dep::Bool
  archive_budget::Bool

  # Dry deposition

  dry_dep_chm::Array{AbstractFloat,3}
  map_dry_dep_chm::DgnMap
  archive_dry_dep_chm::Bool

  dry_dep_mix::Array{AbstractFloat,3}
  map_dry_dep_mix::DgnMap
  archive_dry_dep_mix::Bool

  dry_dep::Array{AbstractFloat,3}
  map_dry_dep::DgnMap
  archive_dry_dep::Bool

  dry_dep_vel::Array{AbstractFloat,3}
  map_dry_dep_vel::DgnMap
  archive_dry_dep_vel::Bool

  # Photolysis

  jval::Array{AbstractFloat,4}
  map_jval::DgnMap
  archive_jval::Bool

  jval_o3o1d::Array{AbstractFloat,3}
  archive_jval_o3o1d::Bool

  jval_o3o3p::Array{AbstractFloat,3}
  archive_jval_o3o3p::Bool

  jnoon::Array{AbstractFloat,4}
  map_jnoon::DgnMap
  archive_jnoon::Bool

  jnoon_frac::Matrix{AbstractFloat}
  archive_jnoon_frac::Bool

  uv_flux_diffuse::Array{AbstractFloat,4}
  map_uv_flux_diffuse::DgnMap
  archive_uv_flux_diffuse::Bool

  uv_flux_direct::Array{AbstractFloat,4}
  map_uv_flux_direct::DgnMap
  archive_uv_flux_direct::Bool

  uv_flux_net::Array{AbstractFloat,4}
  map_uv_flux_net::DgnMap
  archive_uv_flux_net::Bool

  # Chemistry

  rxn_rate::Array{AbstractFloat,4}
  map_rxn_rate::DgnMap
  archive_rxn_rate::Bool

  oh_reactivity::Array{AbstractFloat,3}
  archive_oh_reactivity::Bool

  oh_conc_after_chem::Array{AbstractFloat,3}
  archive_oh_conc_after_chem::Bool

  ho2_conc_after_chem::Array{AbstractFloat,3}
  archive_ho2_conc_after_chem::Bool

  o1d_conc_after_chem::Array{AbstractFloat,3}
  archive_o1_dconc_after_chem::Bool

  o3_pconc_after_chem::Array{AbstractFloat,3}
  archive_o3_pconc_after_chem::Bool

  ch4_pseudo_flux::Matrix{AbstractFloat}
  archive_ch4_pseudo_flux::Bool

  loss::Array{AbstractFloat,4}
  map_loss::DgnMap
  archive_loss::Bool

  prod::Array{AbstractFloat,4}
  map_prod::DgnMap
  archive_prod::Bool

  #ifdef MODEL_GEOS
  nox_tau::Array{AbstractFloat,3}
  archive_nox_tau::Bool

  trop_nox_tau::Matrix{AbstractFloat}
  archive_trop_nox_tau::Bool
  #endif

  # Aerosol characteristics

  aer_hyg_growth::Array{AbstractFloat,4}
  map_aer_hyg_growth::DgnMap
  archive_aer_hyg_growth::Bool

  aer_aq_vol::Array{AbstractFloat,3}
  archive_aer_aq_vol::Bool

  aer_surf_area_hyg::Array{AbstractFloat,4}
  map_aer_surf_area_hyg::DgnMap
  archive_aer_surf_area_hyg::Bool

  aer_surf_area_dust::Array{AbstractFloat,3}
  archive_aer_surf_area_dust::Bool

  aer_surf_area_sla::Array{AbstractFloat,3}
  archive_aer_surf_area_sla::Bool

  aer_surf_area_psc::Array{AbstractFloat,3}
  archive_aer_surf_area_psc::Bool

  aer_num_den_sla::Array{AbstractFloat,3}
  archive_aer_num_den_sla::Bool

  aer_num_den_psc::Array{AbstractFloat,3}
  archive_aer_num_den_psc::Bool

  # Aerosol optical depths

  aod_dust::Array{AbstractFloat,3}
  archive_aod_dust::Bool
  archive_aod::Bool

  aod_dust_wl1::Array{AbstractFloat,4}
  map_aod_dust_wl1::DgnMap
  archive_aod_dust_wl1::Bool

  aod_dust_wl2::Array{AbstractFloat,4}
  map_aod_dust_wl2::DgnMap
  archive_aod_dust_wl2::Bool

  aod_dust_wl3::Array{AbstractFloat,4}
  map_aod_dust_wl3::DgnMap
  archive_aod_dust_wl3::Bool

  aod_hyg_wl1::Array{AbstractFloat,4}
  map_aod_hyg_wl1::DgnMap
  archive_aod_hyg_wl1::Bool

  aod_hyg_wl2::Array{AbstractFloat,4}
  map_aod_hyg_wl2::DgnMap
  archive_aod_hyg_wl2::Bool

  aod_hyg_wl3::Array{AbstractFloat,4}
  map_aod_hyg_wl3::DgnMap
  archive_aod_hyg_wl3::Bool

  aod_soa_from_aq_isop_wl1::Array{AbstractFloat,3}
  archive_aod_soa_from_aq_isop_wl1::Bool

  aod_soa_from_aq_isop_wl2::Array{AbstractFloat,3}
  archive_aod_soa_from_aq_isop_wl2::Bool

  aod_soa_from_aq_isop_wl3::Array{AbstractFloat,3}
  archive_aod_soa_from_aq_isop_wl3::Bool

  aod_sla_wl1::Array{AbstractFloat,3}
  archive_aod_sla_wl1::Bool
  archive_aod_strat::Bool

  aod_sla_wl2::Array{AbstractFloat,3}
  archive_aod_sla_wl2::Bool

  aod_sla_wl3::Array{AbstractFloat,3}
  archive_aod_sla_wl3::Bool

  aod_psc_wl1::Array{AbstractFloat,3}
  archive_aod_psc_wl1::Bool

  aod_psc_wl2::Array{AbstractFloat,3}
  archive_aod_psc_wl2::Bool

  aod_psc_wl3::Array{AbstractFloat,3}
  archive_aod_psc_wl3::Bool

  # Aerosol mass and PM2.5

  aer_mass_asoa::Array{AbstractFloat,3}
  archive_aer_mass_asoa::Bool
  archive_aer_mass::Bool

  aer_mass_bc::Array{AbstractFloat,3}
  archive_aer_mass_bc::Bool

  aer_mass_hms::Array{AbstractFloat,3}
  archive_aer_mass_hms::Bool

  aer_mass_indiol::Array{AbstractFloat,3}
  archive_aer_mass_indiol::Bool

  aer_mass_isn1oa::Array{AbstractFloat,3}
  archive_aer_mass_lvocoa::Bool

  aer_mass_lvocoa::Array{AbstractFloat,3}
  archive_aer_mass_isn1oa::Bool

  aer_mass_nh4::Array{AbstractFloat,3}
  archive_aer_mass_nh4::Bool

  aer_mass_nit::Array{AbstractFloat,3}
  archive_aer_mass_nit::Bool

  aer_mass_opoa::Array{AbstractFloat,3}
  archive_aer_mass_opoa::Bool

  aer_mass_poa::Array{AbstractFloat,3}
  archive_aer_mass_poa::Bool

  aer_mass_sal::Array{AbstractFloat,3}
  archive_aer_mass_sal::Bool

  aer_mass_so4::Array{AbstractFloat,3}
  archive_aer_mass_so4::Bool

  aer_mass_soagx::Array{AbstractFloat,3}
  archive_aer_mass_soagx::Bool

  aer_mass_soaie::Array{AbstractFloat,3}
  archive_aer_mass_soaie::Bool

  aer_mass_tsoa::Array{AbstractFloat,3}
  archive_aer_mass_tsoa::Bool

  beta_no::Array{AbstractFloat,3}
  archive_beta_no::Bool

  pm25::Array{AbstractFloat,3}
  archive_pm25::Bool

  #zhaisx
  pm10::Array{AbstractFloat,3}
  archive_pm10::Bool

  total_oa::Array{AbstractFloat,3}
  archive_total_oa::Bool

  total_oc::Array{AbstractFloat,3}
  archive_total_oc::Bool

  total_biogenic_oa::Array{AbstractFloat,3}
  archive_total_biogenic_oa::Bool

  # Advection

  adv_flux_zonal::Array{AbstractFloat,4}
  map_adv_flux_zonal::DgnMap
  archive_adv_flux_zonal::Bool

  adv_flux_merid::Array{AbstractFloat,4}
  map_adv_flux_merid::DgnMap
  archive_adv_flux_merid::Bool

  adv_flux_vert::Array{AbstractFloat,4}
  map_adv_flux_vert::DgnMap
  archive_adv_flux_vert::Bool

  # Mixing

  pbl_mix_frac::Array{AbstractFloat,3}
  archive_pbl_mix_frac::Bool

  pbl_flux::Array{AbstractFloat,4}
  map_pbl_flux::DgnMap
  archive_pbl_flux::Bool

  # Convection and Wet_dep

  cloud_conv_flux::Array{AbstractFloat,4}
  map_cloud_conv_flux::DgnMap
  archive_cloud_conv_flux::Bool

  wet_loss_conv_frac::Array{AbstractFloat,4}
  map_wet_loss_conv_frac::DgnMap
  archive_wet_loss_conv_frac::Bool

  wet_loss_conv::Array{AbstractFloat,4}
  map_wet_loss_conv::DgnMap
  archive_wet_loss_conv::Bool

  wet_loss_ls::Array{AbstractFloat,4}
  map_wet_loss_ls::DgnMap
  archive_wet_loss_ls::Bool

  # These are obsolete diagnostics
  #REAL(f4),  POINTER :: Precip_fracLS    (:,:,:  )
  #REAL(f4),  POINTER :: Rain_fracLS      (:,:,:,:)
  #REAL(f4),  POINTER :: Wash_fracLS      (:,:,:,:)
  #Archive_precip_fracLS::Bool
  #Archive_rain_fracLS::Bool
  #Archive_wash_fracLS::Bool

  # Carbon aerosols

  prod_bcpi_from_bcpo::Array{AbstractFloat,3}
  archive_prod_bcpi_from_bcpo::Bool

  prod_ocpi_from_ocpo::Array{AbstractFloat,3}
  archive_prod_ocpi_from_ocpo::Bool

  #  Sulfur aerosols prod & loss
  prod_so2_from_dm_sand_oh::Array{AbstractFloat,3}
  archive_prod_so2_from_dm_sand_oh::Bool

  prod_so2_from_dm_sand_no3::Array{AbstractFloat,3}
  archive_prod_so2_from_dm_sand_no3::Bool

  prod_so2_from_dms::Array{AbstractFloat,3}
  archive_prod_so2_from_dms::Bool

  prod_msa_from_dms::Array{AbstractFloat,3}
  archive_prod_msa_from_dms::Bool

  prod_nit_from_hno3_uptake_on_dust::Array{AbstractFloat,3}
  archive_prod_nit_from_hno3_uptake_on_dust::Bool

  prod_so4_from_gas_phase::Array{AbstractFloat,3}
  archive_prod_so4_from_gas_phase::Bool

  prod_so4_from_h2o2_in_cloud::Array{AbstractFloat,3}
  archive_prod_so4_from_h2o2_in_cloud::Bool

  prod_so4_from_o3_in_cloud::Array{AbstractFloat,3}
  archive_prod_so4_from_o3_in_cloud::Bool

  prod_so4_from_o2_in_cloud_metal::Array{AbstractFloat,3}
  archive_prod_so4_from_o2_in_cloud_metal::Bool

  prod_so4_from_o3_in_sea_salt::Array{AbstractFloat,3}
  archive_prod_so4_from_o3_in_sea_salt::Bool

  prod_so4_from_oxidation_on_dust::Array{AbstractFloat,3}
  archive_prod_so4_from_oxidation_on_dust::Bool

  prod_so4_from_uptake_of_h2so4g::Array{AbstractFloat,3}
  archive_prod_so4_from_uptake_of_h2so4g::Bool

  prod_so4_from_hobr_in_cloud::Array{AbstractFloat,3}
  archive_prod_so4_from_hobr_in_cloud::Bool

  prod_so4_from_sro3::Array{AbstractFloat,3}
  archive_prod_so4_from_sro3::Bool

  prod_so4_from_srhobr::Array{AbstractFloat,3}
  archive_prod_so4_from_srhobr::Bool

  prod_so4_from_o3s::Array{AbstractFloat,3}
  archive_prod_so4_from_o3s::Bool

  loss_hno3_on_sea_salt::Array{AbstractFloat,3}
  archive_loss_hno3_on_sea_salt::Bool

  prod_hms_from_so2_and_hcho_in_cloud::Array{AbstractFloat,3}
  archive_prod_hms_from_so2_and_hcho_in_cloud::Bool

  prod_so2_and_hcho_from_hms_in_cloud::Array{AbstractFloat,3}
  archive_prod_so2_and_hcho_from_hms_in_cloud::Bool

  prod_so4_from_hms_in_cloud::Array{AbstractFloat,3}
  archive_prod_so4_from_hms_in_cloud::Bool

  # O3 and HNO3 at a given height above the surface

  dry_dep_ra_alt1::Matrix{AbstractFloat}
  archive_dry_dep_ra_alt1::Bool

  dry_dep_vel_for_alt1::Array{AbstractFloat,3}
  archive_dry_dep_vel_for_alt1::Bool

  species_conc_alt1::Array{AbstractFloat,3}
  archive_species_conc_alt1::Bool
  archive_conc_above_sfc::Bool

  # Time spent in the troposphere

  frac_of_time_in_trop::Array{AbstractFloat,3}
  archive_frac_of_time_in_trop::Bool

  # KPP solver diagnostics

  kpp_int_counts::Array{AbstractFloat,3}
  archive_kpp_int_counts::Bool

  kpp_jac_counts::Array{AbstractFloat,3}
  archive_kpp_jac_counts::Bool

  kpp_tot_steps::Array{AbstractFloat,3}
  archive_kpp_tot_steps::Bool

  kpp_acc_steps::Array{AbstractFloat,3}
  archive_kpp_acc_steps::Bool

  kpp_rej_steps::Array{AbstractFloat,3}
  archive_kpp_rej_steps::Bool

  kpp_lu_decomps::Array{AbstractFloat,3}
  archive_kpp_lu_decomps::Bool

  kpp_substs::Array{AbstractFloat,3}
  archive_kpp_substs::Bool

  kpp_sm_decomps::Array{AbstractFloat,3}
  archive_kpp_sm_decomps::Bool

  archive_kpp_diags::Bool

  # Chemistry metrics (e.g. mean oh, MCF lifetime, CH4 lifetime)

  air_mass_column_full::Matrix{AbstractFloat}
  archive_air_mass_column_full::Bool
  archive_metrics::Bool

  air_mass_column_trop::Matrix{AbstractFloat}
  archive_air_mass_column_trop::Bool

  ch4_emission::Matrix{AbstractFloat}
  archive_ch4_emission::Bool

  ch4_mass_column_full::Matrix{AbstractFloat}
  archive_ch4_mass_column_full::Bool

  ch4_mass_column_trop::Matrix{AbstractFloat}
  archive_ch4_mass_column_trop::Bool

  loss_oh_by_ch4_column_trop::Matrix{AbstractFloat}
  archive_loss_oh_by_ch4_column_trop::Bool

  loss_oh_by_mcf_column_trop::Matrix{AbstractFloat}
  archive_loss_oh_by_mcf_column_trop::Bool

  oh_wgt_by_air_mass_column_full::Matrix{AbstractFloat}
  archive_oh_wgt_by_air_mass_column_full::Bool

  oh_wgt_by_air_mass_column_trop::Matrix{AbstractFloat}
  archive_oh_wgt_by_air_mass_column_trop::Bool

  #----------------------------------------------------------------------
  # Specialty Simulation Diagnostic Arrays
  #----------------------------------------------------------------------

  # Transport_tracers simulation

  pb_from_rn_decay::Array{AbstractFloat,3}
  archive_pb_from_rn_decay::Bool

  rad_decay::Array{AbstractFloat,4}
  map_rad_decay::DgnMap
  archive_rad_decay::Bool

  # CO2 specialty simulation

  prod_co2_from_co::Array{AbstractFloat,3}
  archive_prod_co2_from_co::Bool

  # CH4 specialty simulation

  loss_ch4_by_clin_trop::Array{AbstractFloat,3}
  archive_loss_ch4_by_clin_trop::Bool

  loss_ch4_by_oh_in_trop::Array{AbstractFloat,3}
  archive_loss_ch4_by_oh_in_trop::Bool

  loss_ch4_in_strat::Array{AbstractFloat,3}
  archive_loss_ch4_in_strat::Bool

  #  Tagged CO simulation
  prod_co_from_ch4::Array{AbstractFloat,3}
  archive_prod_co_from_ch4::Bool

  prod_co_from_nmvoc::Array{AbstractFloat,3}
  archive_prod_co_from_nmvoc::Bool

  # Persistent Organic Pollutants (POPS) simulation

  loss_poppocpo_by_gas_phase::Array{AbstractFloat,3}
  archive_loss_poppocpo_by_gas_phase::Bool

  prod_poppocpo_from_gas_phase::Array{AbstractFloat,3}
  archive_prod_poppocpo_from_gas_phase::Bool

  loss_poppbcpo_by_gas_phase::Array{AbstractFloat,3}
  archive_loss_poppbcpo_by_gas_phase::Bool

  prod_poppbcpo_from_gas_phase::Array{AbstractFloat,3}
  archive_prod_poppbcpo_from_gas_phase::Bool

  prod_popg_from_oh::Array{AbstractFloat,3}
  archive_prod_popg_from_oh::Bool

  prod_poppocpo_from_o3::Array{AbstractFloat,3}
  archive_prod_poppocpo_from_o3::Bool

  prod_poppocpi_from_o3::Array{AbstractFloat,3}
  archive_prod_poppocpi_from_o3::Bool

  prod_poppbcpi_from_o3::Array{AbstractFloat,3}
  archive_prod_poppbcpi_from_o3::Bool

  prod_poppbcpo_from_o3::Array{AbstractFloat,3}
  archive_prod_poppbcpo_from_o3::Bool

  prod_poppocpo_from_no3::Array{AbstractFloat,3}
  archive_prod_poppocpo_from_no3::Bool

  prod_poppocpi_from_no3::Array{AbstractFloat,3}
  archive_prod_poppocpi_from_no3::Bool

  prod_poppbcpi_from_no3::Array{AbstractFloat,3}
  archive_prod_poppbcpi_from_no3::Bool

  prod_poppbcpo_from_no3::Array{AbstractFloat,3}
  archive_prod_poppbcpo_from_no3::Bool

  # Hg specialty simulation
  #  -- emissions quantities (e.g. for HEMCO manual diagnostics)
  emis_hg0_anthro::Matrix{AbstractFloat}
  archive_emis_hg0_anthro::Bool

  emis_hg0_biomass::Matrix{AbstractFloat}
  archive_emis_hg0_biomass::Bool

  emis_hg0_geogenic::Matrix{AbstractFloat}
  archive_emis_hg0_geogenic::Bool

  emis_hg0_land::Matrix{AbstractFloat}
  archive_emis_hg0_land::Bool

  emis_hg0_ocean::Matrix{AbstractFloat}
  archive_emis_hg0_ocean::Bool

  emis_hg0_snow::Matrix{AbstractFloat}
  archive_emis_hg0_snow::Bool

  emis_hg0_soil::Matrix{AbstractFloat}
  archive_emis_hg0_soil::Bool

  emis_hg0_vegetation::Matrix{AbstractFloat}
  archive_emis_hg0_vegetation::Bool

  emis_hg2hg_panthro::Matrix{AbstractFloat}
  archive_emis_hg2hg_panthro::Bool

  emis_hg2_snow_to_ocean::Matrix{AbstractFloat}
  archive_emis_hg2_snow_to_ocean::Bool

  emis_hg2_rivers::Matrix{AbstractFloat}
  archive_emis_hg2_rivers::Bool

  flux_hg2hgp_from_air_to_snow::Matrix{AbstractFloat}
  archive_flux_hg2hgp_from_air_to_snow::Bool
  #
  #  -- oceanic quantities
  flux_hg0_from_ocean_to_air::Matrix{AbstractFloat}
  archive_flux_hg0_from_air_to_ocean::Bool

  flux_hg0_from_air_to_ocean::Matrix{AbstractFloat}
  archive_flux_hg0_from_ocean_to_air::Bool

  flux_hg2hgp_from_air_to_ocean::Matrix{AbstractFloat}
  archive_flux_hg2hgp_from_air_to_ocean::Bool

  flux_hg2_to_deep_ocean::Matrix{AbstractFloat}
  archive_flux_hg2_to_deep_ocean::Bool

  flux_oc_to_deep_ocean::Matrix{AbstractFloat}
  archive_flux_oc_to_deep_ocean::Bool

  mass_hg0_in_ocean::Matrix{AbstractFloat}
  archive_mass_hg0_in_ocean::Bool

  mass_hg2_in_ocean::Matrix{AbstractFloat}
  archive_mass_hg2_in_ocean::Bool

  mass_hgp_in_ocean::Matrix{AbstractFloat}
  archive_mass_hgp_in_ocean::Bool

  mass_hg_total_in_ocean::Matrix{AbstractFloat}
  archive_mass_hg_total_in_ocean::Bool
  #
  #  -- chemistry quantities
  conc_br::Array{AbstractFloat,3}
  archive_conc_br::Bool

  conc_bro::Array{AbstractFloat,3}
  archive_conc_bro::Bool

  loss_hg2_by_sea_salt::Array{AbstractFloat,3}
  archive_loss_hg2_by_sea_salt::Bool

  loss_rate_hg2_by_sea_salt::Matrix{AbstractFloat}
  archive_loss_rate_hg2_by_sea_salt::Bool

  polar_conc_br::Array{AbstractFloat,3}
  archive_polar_conc_br::Bool

  polar_conc_bro::Array{AbstractFloat,3}
  archive_polar_conc_bro::Bool

  polar_conc_o3::Array{AbstractFloat,3}
  archive_polar_conc_o3::Bool

  prod_hg2_from_br::Array{AbstractFloat,3}
  archive_prod_hg2_from_br::Bool

  prod_hg2_from_bry::Array{AbstractFloat,3}
  archive_prod_hg2_from_bry::Bool

  prod_hg2_from_cly::Array{AbstractFloat,3}
  archive_prod_hg2_from_cly::Bool

  prod_hg2_from_hg0::Array{AbstractFloat,3}
  archive_prod_hg2_from_hg0::Bool

  prod_hg2_from_hgbr_plus_br2::Array{AbstractFloat,3}
  archive_prod_hg2_from_hgbr_plus_br2::Bool

  prod_hg2_from_hgbr_plus_brbro::Array{AbstractFloat,3}
  archive_prod_hg2_from_hgbr_plus_brbro::Bool

  prod_hg2_from_hgbr_plus_brclo::Array{AbstractFloat,3}
  archive_prod_hg2_from_hgbr_plus_brclo::Bool

  prod_hg2_from_hgbr_plus_brho2::Array{AbstractFloat,3}
  archive_prod_hg2_from_hgbr_plus_brho2::Bool

  prod_hg2_from_hgbr_plus_brno2::Array{AbstractFloat,3}
  archive_prod_hg2_from_hgbr_plus_brno2::Bool

  prod_hg2_from_hgbr_plus_broh::Array{AbstractFloat,3}
  archive_prod_hg2_from_hgbr_plus_broh::Bool

  prod_hg2_from_oh::Array{AbstractFloat,3}
  archive_prod_hg2_from_oh::Bool

  prod_hg2_from_o3::Array{AbstractFloat,3}
  archive_prod_hg2_from_o3::Bool

  particulate_bound_hg::Array{AbstractFloat,3}
  archive_particulate_bound_hg::Bool

  reactive_gaseous_hg::Array{AbstractFloat,3}
  archive_reactive_gaseous_hg::Bool

  # From Viral Shah (MSL, 7.1.21)
  hgbr_after_chem::Array{AbstractFloat,3}
  archive_hgbr_after_chem::Bool

  hgcl_after_chem::Array{AbstractFloat,3}
  archive_hgcl_after_chem::Bool

  hgoh_after_chem::Array{AbstractFloat,3}
  archive_hgoh_after_chem::Bool

  hgbro_after_chem::Array{AbstractFloat,3}
  archive_hgbro_after_chem::Bool

  hgclo_after_chem::Array{AbstractFloat,3}
  archive_hgclo_after_chem::Bool

  hgoho_after_chem::Array{AbstractFloat,3}
  archive_hgoho_after_chem::Bool

  hg2g_to_hg2p::Array{AbstractFloat,3}
  archive_hg2g_to_hg2p::Bool

  hg2p_to_hg2g::Array{AbstractFloat,3}
  archive_hg2p_to_hg2g::Bool

  hg2_gas_to_hg2_strp::Array{AbstractFloat,3}
  archive_hg2_gas_to_hg2_strp::Bool

  hg2_gas_to_ssa::Array{AbstractFloat,3}
  archive_hg2_gas_to_ssa::Bool

  # Simulation with RRTMG

  n_rad_out::Integer
  rad_out_ind::Vector{Integer}
  rad_out_name::Vector{String}

  rad_all_sky_lw_surf::Array{AbstractFloat,3}
  archive_rad_all_sky_lw_surf::Bool

  rad_all_sky_lw_toa::Array{AbstractFloat,3}
  archive_rad_all_sky_lw_toa::Bool

  rad_all_sky_sw_surf::Array{AbstractFloat,3}
  archive_rad_all_sky_sw_surf::Bool

  rad_all_sky_sw_toa::Array{AbstractFloat,3}
  archive_rad_all_sky_sw_toa::Bool

  rad_clr_sky_lw_surf::Array{AbstractFloat,3}
  archive_rad_clr_sky_lw_surf::Bool

  rad_clr_sky_lw_toa::Array{AbstractFloat,3}
  archive_rad_clr_skyLWTOA::Bool

  rad_clr_sky_sw_surf::Array{AbstractFloat,3}
  archive_rad_clr_sky_sw_surf::Bool

  rad_clr_sky_sw_toa::Array{AbstractFloat,3}
  archive_rad_clr_sky_sw_toa::Bool

  rad_aod_wl1::Array{AbstractFloat,3}
  archive_rad_aod_wl1::Bool

  rad_aod_wl2::Array{AbstractFloat,3}
  archive_rad_aod_wl2::Bool

  rad_aod_wl3::Array{AbstractFloat,3}
  archive_rad_aod_wl3::Bool

  rad_ssa_wl1::Array{AbstractFloat,3}
  archive_rad_ssa_wl1::Bool

  rad_ssa_wl2::Array{AbstractFloat,3}
  archive_rad_ssa_wl2::Bool

  rad_ssa_wl3::Array{AbstractFloat,3}
  archive_rad_ssa_wl3::Bool

  rad_asym_wl1::Array{AbstractFloat,3}
  archive_rad_asym_wl1::Bool

  rad_asym_wl2::Array{AbstractFloat,3}
  archive_rad_asym_wl2::Bool

  rad_asym_wl3::Array{AbstractFloat,3}
  archive_rad_asym_wl3::Bool

  archive_rad_optics::Bool

  #----------------------------------------------------------------------
  # Variables for the Obs_pack diagnostic
  # NOTE: Obs_pack archives point data, so don't register these
  # as the Obs_pack file format won't be COARDS-compliant#
  #----------------------------------------------------------------------

  # Obs_pack File variables
  do_obs_pack::Bool
  obs_pack_f_id::Integer
  obs_pack_in_file::String
  obs_pack_out_file::String

  # Obs_pack Inputs
  obs_pack_n_obs::Integer
  obs_pack_id::Vector{String}
  obs_pack_n_samples::Vector{Integer}
  obs_pack_strategy::Vector{Integer}
  obs_pack_latitude::Vector{AbstractFloat}
  obs_pack_longitude::Vector{AbstractFloat}
  obs_pack_altitude::Vector{AbstractFloat}

  # Obs_pack time and averaging interval variables
  obs_pack_ival_length::AbstractFloat
  obs_pack_ival_start::Vector{AbstractFloat}
  obs_pack_ival_center::Vector{AbstractFloat}
  obs_pack_ival_end::Vector{AbstractFloat}

  # Obs_pack outputs (add more if necessary)
  obs_pack_p::Vector{AbstractFloat}
  obs_pack_u::Vector{AbstractFloat}
  obs_pack_v::Vector{AbstractFloat}
  obs_pack_blh::Vector{AbstractFloat}
  obs_pack_q::Vector{AbstractFloat}
  obs_pack_t::Vector{AbstractFloat}

  # Obs_pack species and metadata variables
  obs_pack_n_species::Integer
  obs_pack_species::Matrix{AbstractFloat}
  obs_pack_species_ind::Vector{Integer}
  obs_pack_species_name::Vector{String}
  obs_pack_species_lname::Vector{String}

  #ifdef MODEL_GEOS
  #----------------------------------------------------------------------
  # The following diagnostics are only used when
  # GEOS-Chem is interfaced into the NASA-GEOS ESM
  #----------------------------------------------------------------------

  monin_obukhov::Matrix{AbstractFloat}
  archive_monin_obukhov::Bool

  bry::Array{AbstractFloat,3}
  archive_bry::Bool

  noy::Array{AbstractFloat,3}
  archive_noy::Bool

  cly::Array{AbstractFloat,3}
  archive_cly::Bool

  organic_cl::Array{AbstractFloat,3}
  archive_organic_cl::Bool

  o3_mass::Array{AbstractFloat,3}
  archive_o3_mass::Bool

  gccto3::Matrix{AbstractFloat}
  archive_gccto3::Bool

  gcctto3::Matrix{AbstractFloat}
  archive_gcctto3::Bool

  o3_mass::Array{AbstractFloat,3}
  archive_o3_mass::Bool

  chem_top::Matrix{AbstractFloat}
  archive_chem_top::Bool

  chem_tropp::Matrix{AbstractFloat}
  archive_chem_tropp::Bool

  convcld_top::Matrix{AbstractFloat}
  archive_convcld_top::Bool

  extra_lnlevs::Matrix{AbstractFloat}
  archive_extra_lnlevs::Bool

  extra_lniter::Matrix{AbstsractFloat}
  archive_extra_lniter::Bool

  lightning_potential::Matrix{AbstractFloat}
  archive_lightning_potential::Bool

  # Chemistry diagnostics

  kpp_error::Array{AbstractFloat,3}
  archive_kpp_error::Bool

  o3_conc_after_chem::Array{AbstractFloat,3}
  archive_o3_conc_after_chem::Bool

  ro2_conc_after_chem::Array{AbstractFloat,3}
  archive_ro2_conc_after_chem::Bool

  # PM2.5 diagnostics

  "PM25 nitrates"
  pm25ni::Array{AbstractFloat,3}
  archive_pm25ni::Bool

  "PM25 sulfates"
  pm25su::Array{AbstractFloat,3}
  archive_pm25su::Bool

  "PM25 OC"
  pm25oc::Array{AbstractFloat,3}
  archive_pm25oc::Bool

  "PM25 BC"
  pm25bc::Array{AbstractFloat,3}
  archive_pm25bc::Bool

  "PM25 dust"
  pm25du::Array{AbstractFloat,3}
  archive_pm25du::Bool

  "PM25 sea salt"
  pm25ss::Array{AbstractFloat,3}
  archive_pm25ss::Bool

  "PM25 SOA"
  pm25soa::Array{AbstractFloat,3}
  archive_pm25soa::Bool

  # Species diagnostics
  trop_col::Array{AbstractFloat,3}
  map_trop_col::DgnMap
  archive_trop_col::Bool

  tot_col::Array{AbstractFloat,3}
  map_tot_col::DgnMap
  archive_tot_col::Bool
  #endif

  #ifdef MODEL_WRF
  #----------------------------------------------------------------------
  # The following diagnostics are only used when
  # GEOS-Chem is interfaced into WRF (as WRF-GC)
  #----------------------------------------------------------------------
  kpp_error::Array{AbstractFloat,3}
  archive_kpp_error::Bool
  #endif

  #----------------------------------------------------------------------
  # Registry of variables contained within State_diag
  #----------------------------------------------------------------------
  "Name of this state"
  state::String = "DIAG"
  "Registry object"
  registry::MetaRegItem
  "Lookup table"
  reg_dict::dictionary_t
end

end

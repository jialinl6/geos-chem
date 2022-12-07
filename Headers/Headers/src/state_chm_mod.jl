import DocStringExtensions

"""
Contains the derived type used to define the Chemistry State object for GEOS-Chem.

This module also contains the routines that allocate and deallocate memory to the Chemistry State object. The chemistry state object is not defined in this module. It must be be declared as variable in the top-level driver routine, and then passed to lower-level routines as an argument.
"""
module state_chm_mod

export chm_state

"""
Derived type for Chemistry State

$(TYPEDFIELDS)
"""
mutable struct ChmState
	# Count of each type of species
	"species (all)"
	n_species::Integer
	"advected species"
	n_advect::Integer
	"of Aerosol Species"
	n_aero_spc::Integer
	"of Aerosol Types"
	n_aero_type::Integer
	"dryalt species"
	n_dry_alt::Integer
	"drydep species"
	n_dry_dep::Integer
	"gas phase species"
	n_gas_spc::Integer
	"hygroscopic growth"
	n_hyg_grth::Integer
	"KPP variable species"
	n_kpp_var::Integer
	"KPP fixed species"
	n_kpp_fix::Integer
	"KPP chem species"
	n_kpp_spc::Integer
	"of loss species"
	n_loss::Integer
	"of omitted species"
	n_omitted::Integer
	"photolysis species"
	n_photol::Integer
	"of prod species"
	n_prod::Integer
	"of radionuclides"
	n_rad_nucl::Integer
	"wetdep species"
	n_wet_dep::Integer

	# Mapping vectors to subset types of species
	"Advected species IDs"
	map_advect::Vector{Integer}
	"Aerosol species IDs"
	map_aero::Vector{Integer}
	"Dryalt species IDs"
	map_dry_alt::Vector{Integer}
	"Drydep species IDs"
	map_dry_dep::Vector{Integer}
	"Gas species IDs"
	map_gas_spc::Vector{Integer}
	"HygGrth species IDs"
	map_hyg_grth::Vector{Integer}
	"Kpp variable spc IDs"
	map_kpp_var::Vector{Integer}
	"KPP fixed species IDs"
	map_kpp_fix::Vector{Integer}
	"KPP chem species IDs"
	Map_kpp_spc::Vector{Integer}
	"Loss diag species"
	map_loss::Vector{Integer}
	"ID's and names"
	name_loss::String
	"Photolysis species IDs"
	map_photol::Vector{Integer}
	"Prod diag species"
	map_prod::Vector{Integer}
	"ID and names"
	name_prod::String
	"Radionuclide IDs"
	map_rad_nucl::Vector{Integer}
	"Wetdep species IDs"
	map_wet_dep::Vector{Integer}
	"Wavelength bins in fjx"
	map_wl::Vector{Integer}

	# Physical properties & indices for each species
	"GC Species database"
	spc_data::Vector{SpcPtr}
	"Species dictionary"
	spc_dict::dictionary_t

	# Chemical Species
	"Vector for species concentrations [kg/kg dry air]"
	species::Vector{SpcConc}
	"Species units"
	spc_units::String
	# TODO:
	#ifdef ADJOINT
	"Species adjoint variables"
	species_adj::Array{AbstractFloat, 4}
	"cost function volume mask"
	cost_func_mask::Array{AbstractFloat, 3}
	#endif

	# Boundary conditions
	"Boundary conditions [kg/kg dry air]"
	boundary_cond::Array{AbstractFloat, 4}
	
	# RRTMG state variables
	rrtmg_i_seed::Integer
	rrtmg_i_cld::Integer

	# Aerosol quantities
	"Aerosol Area [cm2/cm3]"
	aero_area::Array{AbstractFloat, 4}
	"Aerosol Radius [cm]"
	aero_radi::Array{AbstractFloat, 4}
	"Aerosol Area [cm2/cm3]"
	wet_aero_area::Array{AbstractFloat, 4}
	"Aerosol Radius [cm]"
	wet_aero_radi::Array{AbstractFloat, 4}
	"Aerosol water [cm3/cm3]"
	aero_h2o::Array{AbstractFloat, 4}
	"N2O5 aerosol uptake"
	gamma_n2o5::Array{AbstractFloat, 4}
	"Sea-salt alkalinity[-]"
	ss_alk::Array{AbstractFloat, 4}
	"H2O2, SO2 [v/v]"
	h2o2_after_chem::Array{AbstractFloat, 3}
	" after sulfate chem"
	so2_after_chem::Array{AbstractFloat, 3}
	"OM:OC Ratio [unitless]"
	omoc::Array{AbstractFloat, 2}
	"OM:OC Ratio (OCFPOA)"
	omoc_poa::Array{AbstractFloat, 2}
	"OM:OC Ratio (OCFOPOA)"
	omoc_opoa::Array{AbstractFloat, 2}
	"Fine Cl- Area [cm2/cm3]"
	acl_area::Array{AbstractFloat, 3}
	"Fine Cl- Radius [cm]"
	acl_radi::Array{AbstractFloat, 3}
	qlxph_cloud::Array{AbstractFloat, 3}
	"Soil dust [kg/m3]"
	soil_dust::Array{AbstractFloat, 4}
	"Sesquiterpenes mass [kg/box]"
	orvc_sesq::Array{AbstractFloat, 3}

	# Fields for nitrogen deposition
	"Dry deposited N"
	dry_dep_nitrogen::Matrix{AbstractFloat}
	"Wet deposited N"
	wet_dep_nitrogen::Matrix{AbstractFloat}

	# Cloud quantities
	"Cloud pH [-]"
	ph_cloud::Array{AbstractFloat, 3}
	"Cloud presence [-]"
	is_cloud::Array{AbstractFloat, 3}

	# Fields for KPP solver
	"H-value for Rosenbrock solver"
	kpph_value::Array{AbstractFloat, 3}
	
	# Fields for UCX mechanism
	"PSC type (see Kirner et al. 2011, GMD)"
	state_psc::Array{AbstractFloat, 3}
	"Strat. liquid aerosol reaction cofactors"
	kheti_sla::Array{AbstractFloat, 4}

	# For isoprene SOA via ISORROPIA
	"ISORROPIA aero pH"
	isorrop_aero_ph::Array{AbstractFloat, 4}
	"H+ conc [M]"
	isorrop_hplus::Array{AbstractFloat, 4}
	"ISORROPIA aero H2O"
	isorrop_aero_h2o::Array{AbstractFloat, 4}
	"Sulfate conc [M]"
	isorrop_sulfate::Array{AbstractFloat, 3}
	"Nitrate conc [M]"
	isorrop_nitrate::Array{AbstractFloat, 4}
	"Chloride conc [M]"
	isorrop_chloride::Array{AbstractFloat, 4}
	"Bisulfate conc [M]"
	isorrop_bisulfate::Array{AbstractFloat, 3}

	# For the tagged Hg simulation
	"Hg(0)  ocean mass [kg]"
	ocean_hg0::Matrix{AbstractFloat}
	"Hg(II) ocean mass [kg]"
	ocean_hg2::Matrix{AbstractFloat}
	"HgP    ocean mass [kg]"
	ocean_hgp::Matrix{AbstractFloat}
	"Reducible Hg snowpack on ocean [kg]"
	snow_hg_ocean::Matrix{AbstractFloat}
	"Reducible Hg snowpack on land [kg]"
	snow_hg_land::Matrix{AbstractFloat}
	"Non-reducible Hg snowpack on ocean"
	snow_hg_ocean_stored::Matrix{AbstractFloat}
	"Non-reducible Hg snowpack on land"
	snow_hg_land_stored::Matrix{AbstractFloat}

	# For HOBr + S(IV) heterogeneous chemistry
	"Cloud bisulfite/SO2 ratio"
	hso3_aq::Array{AbstractFloat}
	"Cloud sulfite/SO2 ratio"
	so3_aq::Array{AbstractFloat}
	"Correction factor for HOBr removal by SO2 [unitless]"
	f_update_hobr::Array{AbstractFloat}
	"Correction factor for HOCl removal by SO2 [unitless]"
	f_update_hocl::Array{AbstractFloat}
	
	# Fields for dry deposition
	"Ocn surf iodide [nM]"
	iodide::Matrix{AbatractFloat}
	"Ocn surf salinity [PSU]"
	salinity::Matrix{AbatractFloat}
	"Drydep freq [s-1]"
	dry_dep_freq::Array{AbatractFloat, 3}
	"Dry deposition velocities [m/s] - use REAL8 in drydep"
	dry_dep_vel::Array{AbatractFloat, 3}
	# TODO:
	#if defined( MODEL_GEOS )
	"2m aerodynamic resistance"
	dry_dep_ra2m::Matrix{AbstractFloat}
	"10m aerodynamic resistance"
	dry_dep_ra10m::Matrix{AbstractFloat}
	# #endif
	"OH J-value"
	j_oh::Matrix{AbstractFloat}
	"NO2 J-value"
	j_no2::Matrix{AbstractFloat}

	# Fields for non-local PBL mixing
	surface_flux::Array{AbstractFloat, 3}

	# Fields for Linoz stratospheric ozone algorithm
	"TLSTT (I,J,L,LINOZ_NFIELDS)"
	tlstt::Array{AbstractFloat, 4}
	
	# Fields for Gan Luo et al Wetdep scheme (GMD-12-3439-2019)
	k_rate::Array{AbstractFloat}
	qq3d::Array{AbstractFloat}
	"Rain pH [-]"
	ph_rain::Array{AbstractFloat}
	"Rain pH*QQ3D [-]"
	qq_ph_rain::Array{AbstractFloat}
	"Rain QQ3D [-]"
	qq_rain::Array{AbstractFloat}

	# Fields for setting mean surface CH4 from HEMCO
	sfc_ch4::Matrix{AbstractFloat}

	# Fields for TOMS overhead ozone column data
	"Daily overhead ozone"
	to3_daily::Matrix{AbstractFloat}
	toms1::Matrix{AbstractFloat}
	toms2::Matrix{AbstractFloat}

	# Switches to enable SO2 cloud chemistry and seasalt chemistry in sulfate_mod (TRUE) or in the KPP mechanism (FALSE).
	do_sulfate_mod_cld::Bool
	do_sulfate_mod_seasalt::Bool

	# Registry of variables contained within State_Chm
	"Name of this state"
	state::String = "CHEM"
	"Registry object"
	# I think this would be better as an array than as a linked list
	registry::Ptr{MetaRegItem}
	"Registry lookup table"
	reg_dict::dictionary_t
end

end

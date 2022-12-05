"""
Contains routines for the TPCORE transport scheme, as implemented in the GMI model (cf. John Tannahill), based on Lin Rood 1995. The Harvard Atmospheric Chemistry Modeling Group has added modifications to implement the Philip-Cameron Smith pressure fixer for mass conservation. Mass flux diagnostics have also been added.

### References
- Lin, S.-J., and R. B. Rood, 1996: _Multidimensional flux form semi-Lagrangian transport schemes_, **Mon. Wea. Rev.**, **124**, 2046-2070.
- Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: _A class of the van Leer-type transport schemes and its applications to the moisture transport in a General Circulation Model_, **Mon. Wea. Rev.**, **122**, 1575-1593.

### Selecting E/W, N/S and vertical advection options
The flags IORD, JORD, KORD select which transport schemes are used in the E/W, N/S, and vertical directions, respectively. Here is a list of the possible values that IORD, JORD, KORD may be set to (original notes from S-J Lin):

- 1st order upstream scheme (too diffusive, not a real option; it can be used for debugging purposes; this is THE only known "linear" monotonic advection scheme.).
- 2nd order van Leer (full monotonicity constraint; see Lin et al 1994, MWR)
- monotonic PPM* (Collela & Woodward 1984)
- semi-monotonic PPM (same as 3, but overshoots are allowed)
- positive-definite PPM (constraint on the subgrid distribution is only strong enough to prevent generation of negative values; both overshoots & undershoots are possible).
- un-constrained PPM (nearly diffusion free; faster but positivity of the subgrid distribution is not quaranteed. Use this option only when the fields and winds are very smooth.
- Huynh/Van Leer/Lin full monotonicity constraint. Only KORD can be set to 7 to enable the use of Huynh's 2nd monotonicity constraint for piece-wise parabolic distribution.

Recommended values:

- IORD=JORD=3 for high horizontal resolution.
- KORD=3 or 7

The implicit numerical diffusion decreases as _ORD increases. DO NOT use option 4 or 5 for non-positive definite scalars (such as Ertel Potential Vorticity).

In GEOS-Chem we have been using IORD=3, JORD=3, KORD=7. We have tested the OpenMP parallelization with these options. GEOS-Chem users who wish to use different (I,J,K)ORD options should consider doing single-procsessor vs. multi-processor tests to test the implementation of the parallelization.

### GEOS-4 and GEOS-5 Hybrid Grid Definition

For GEOS-4 and GEOS-5 met fields, the pressure at the bottom edge of grid box (I,J,L) is defined as follows:

``P_{edge}(i, j, l) = A_{k}(L) + [ B_{k}(L) * P_{surface}(I,J) ]``

where

- ``P_{surface}(I,J)`` is the "true" surface pressure at lon,lat (I,J)
- ``A_{k}(L)`` has the same units as surface pressure [hPa]
- ``B_{k}(L)`` is a unitless constant given at level edges

``A_{k}(L)`` and ``B_{k}(L)`` are supplied to us by GMAO.

## REMARKS:

Ak(L) and Bk(L) are defined at layer edges.
```txt
         +- Model top P=ak(1) ------ ak(1),    bk(1)
 delp(1) |  q(i, j, 1)
         +-------------------------- ak(2),    bk(2)
				 
				 
         +-------------------------- ak(k),    bk(k)
 delp(k) |  q(i, j, k)
         +-------------------------- ak(k+1),  bk(k+1)
				 
				 
         +-------------------------- ak(km),   bk(km)
delp(km) |  q(i, j, km)
         +- Earth's surface P=Psfc - ak(km+1), bk(km+1)
```
"""
module tpcore_fvdas_mod

"""
Subroutine Init_Tpcore allocates and initializes all module variables,

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function init_tpcore!(
	im,
	jm,
	km,
	jfirst,
	jlast,
	ng,
	mg,
	dt,
	ae,
	clat,
	rc
)::Nothing
	# ! !USES:
	# !
	# USE PhysConstants
	# USE ErrCode_Mod
	# !
	# ! !INPUT PARAMETERS:
	# !
	# INTEGER,        INTENT(IN)  :: IM         ! Global E-W dimension
	# INTEGER,        INTENT(IN)  :: JM         ! Global N-S dimension
	# INTEGER,        INTENT(IN)  :: KM         ! Vertical dimension
	# INTEGER,        INTENT(IN)  :: NG         ! large ghost width
	# INTEGER,        INTENT(IN)  :: MG         ! small ghost width
	# REAL(fp),       INTENT(IN)  :: dt         ! Time step in seconds
	# REAL(fp),       INTENT(IN)  :: ae         ! Earth"s radius (m)
	# REAL(fp),       INTENT(IN)  :: clat(JM)   ! latitude in radian
	# !
	# ! !OUTPUT PARAMETERS:
	# !
	# INTEGER,        INTENT(OUT) :: JFIRST     ! Local first index for N-S axis
	# INTEGER,        INTENT(OUT) :: JLAST      ! Local last  index for N-S axis
	# INTEGER,        INTENT(OUT) :: RC         ! Success or failure

	elat = zeros(AbstractFloat, jm + 1) # cell edge latitude in radian
	sine = zeros(AbstractFloat, jm + 1)
	sine_25 = zeros(AbstractFloat, jm + 1)

	# Initialize
	rc = gc_success
	errmsg = ""
	thisloc = " -> at Init_Tpcore (in module GeosCore/tpcore_fvas_mod.F90)"

	# NOTE: since we are not using MPI parallelization, we can set JFIRST and JLAST to the global grid limits in latitude. (bmy, 12/3/08)
	jfirst = 1
	jlast = jm

	if jlast - jfirst < 2
		println("Minimum size of subdomain is 3")
	end

	# Allocate arrays
	
	# TODO: Translate to Julia
	
	# ALLOCATE( cosp( JM ), STAT=RC )
	# CALL GC_CheckVar( "tpcore_fvdas_mod.F90:cosp",  0, RC )
	# IF ( RC /= GC_SUCCESS ) RETURN

	# ALLOCATE( cose( JM ), STAT=RC )
	# CALL GC_CheckVar( "tpcore_fvdas_mod.F90:cose",  0, RC )
	# IF ( RC /= GC_SUCCESS ) RETURN

	# ALLOCATE( gw( JM ), STAT=RC )
	# CALL GC_CheckVar( "tpcore_fvdas_mod.F90:gw",    0, RC )
	# IF ( RC /= GC_SUCCESS ) RETURN

	# ALLOCATE( dtdx5( JM ), STAT=RC )
	# CALL GC_CheckVar( "tpcore_fvdas_mod.F90:dtdx5", 0, RC )
	# IF ( RC /= GC_SUCCESS ) RETURN

	# ALLOCATE( dtdy5( JM ), STAT=RC )
	# CALL GC_CheckVar( "tpcore_fvdas_mod.F90:dtdy5", 0, RC )
	# IF ( RC /= GC_SUCCESS ) RETURN

	# ALLOCATE( DLAT( JM ), STAT=RC )                 ! For PJC pressure-fixer
	# CALL GC_CheckVar( "tpcore_fvdas_mod.F90:dlat",  0, RC )
	# IF ( RC /= GC_SUCCESS ) RETURN

	# Define quantities
	
	dlon = 2. * pi / im

	# S. Pole
	elat[1] = -0.5 * π
	sine[1] = -1.0
	sine_25[1] = -1.0
	cose[1] = 0.0

	for j = 2:jm
		elat[j] = 0.5 * (clat[j - 1] + clat[j])
		sine[j] = sin(elat[j])
		sine_25[j] = sin(clat[j])
		cose[j] = cos(elat[j])
	end

	# N. Pole
	elat[jm + 1] = 0.5 * π
	sine[jm + 1] = 1.0
	sine_25[jm + 1] = 1.0

	# Polar cap (S. Pole)
	dlat[1] = 2.0 * (elat[2] - elat[1])
	for j = 2:(jm - 1)
		dlat[j] = elat[j + 1] - elat[j]
	end

	# Polar cap (N. Pole)
	dlat[jm] = 2.0 * (elat[jm + 1] - elat[jm])

	for j = 1:jm
		gw[j] = sine[j + 1] - sine[j]
		cosp[j] = gw[j] / dlat[j]

		dtdx5[j] = 0.5 * dt / (dlon * ae * cosp[j])
		dtdy5[j] = 0.5 * dt / (ae * dlat[j])
	end

	# Echo info to stdout
	println(repeat("=", 79))
	println("TPCORE_FVDAS (based on GMI) Tracer Transport Module successfully initialized")
	println(repeat("=", 79))
end

tpcore_fvdas_first = true
tpcore_fvdas_ilmt = missing::Union{Missing, Integer}
tpcore_fvdas_jlmt = missing::Union{Missing, Integer}
tpcore_fvdas_klmt = missing::Union{Missing, Integer}

"""
Takes horizontal winds on sigma (or hybrid sigma-p) surfaces and calculates mass fluxes, and then updates the 3D mixing ratio fields one time step (tdt).  The basic scheme is a Multi-Dimensional Flux Form Semi-Lagrangian (FFSL) based on the van Leer or PPM (see Lin and Rood, 1995).
	
## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function tpcore_fvdas!(
	dt::AbstractFloat,
	ae::AbstractFloat,
	im::Integer,
	jm::Integer,
	km::Integer,
	jfirst::Integer,
	jlast::Integer,
	ng::Integer,
	mg::Integer,
	nq::Integer,
	ak::Array{AbstractFloat},
	bk::Array{AbstractFloat},
	u::Array{AbstractFloat, 3},
	v::Array{AbstractFloat, 3},
	ps1::Matrix{AbstractFloat},
	ps2::Matrix{AbstractFloat},
	ps::Matrix{AbstractFloat},
	iord::Integer,
	jord::Integer,
	kord::Integer,
	n_adj::Integer,
	xmass::Array{AbstractFloat, 3},
	ymass::Array{AbstractFloat, 3},
	fill::Bool,
	area_m2::Array{AbstractFloat},
	state_chm::ChmState,
	state_diag::DgnState
)::Nothing
	# !
	# ! !USES:
	# !
	# 	! Include files w/ physical constants and met values
	# 	USE PhysConstants
	# 	USE ErrCode_Mod
	# 	USE State_Chm_Mod,  ONLY : ChmState
	# 	USE State_Diag_Mod, ONLY : DgnState
	# 	USE error_mod
	# !
	# ! !INPUT PARAMETERS:
	# !
	# 	! Transport time step [s]
	# 	REAL(fp),  INTENT(IN)  :: dt

	# 	! Earth"s radius [m]
	# 	REAL(fp),  INTENT(IN)  :: ae

	# 	! Global E-W, N-S, and vertical dimensions
	# 	INTEGER, INTENT(IN)    :: IM
	# 	INTEGER, INTENT(IN)    :: JM
	# 	INTEGER, INTENT(IN)    :: KM

	# 	! Latitude indices for local first box and local last box
	# 	! (NOTE: for global grids these are 1 and JM, respectively)
	# 	INTEGER, INTENT(IN)    :: JFIRST
	# 	INTEGER, INTENT(IN)    :: JLAST

	# 	! Primary ghost region
	# 	! (NOTE: only required for MPI parallelization; use 0 otherwise)
	# 	INTEGER, INTENT(IN)    :: ng

	# 	! Secondary ghost region
	# 	! (NOTE: only required for MPI parallelization; use 0 otherwise)
	# 	INTEGER, INTENT(IN)    :: mg

	# 	! Ghosted latitudes (3 required by PPM)
	# 	! (NOTE: only required for MPI parallelization; use 0 otherwise)
	# 	INTEGER, INTENT(IN)    :: nq

	# 	! Flags to denote E-W, N-S, and vertical transport schemes
	# 	INTEGER, INTENT(IN)    :: iord
	# 	INTEGER, INTENT(IN)    :: jord
	# 	INTEGER, INTENT(IN)    :: kord

	# 	! Number of adjustments to air_mass_flux (0 = no adjustment)
	# 	INTEGER, INTENT(IN)    :: n_adj

	# 	! Ak and Bk coordinates to specify the hybrid grid
	# 	! (see the REMARKS section below)
	# 	REAL(fp),  INTENT(IN)  :: ak(KM+1)
	# 	REAL(fp),  INTENT(IN)  :: bk(KM+1)

	# 	! u-wind (m/s) at mid-time-level (t=t+dt/2)
	# 	REAL(fp),  INTENT(IN)  :: u(:,:,:)

	# 	! E/W and N/S mass fluxes [kg/s]
	# 	! (These are computed by the pressure fixer, and passed into TPCORE)
	# 	REAL(fp),  INTENT(IN)  :: XMASS(:,:,:)
	# 	REAL(fp),  INTENT(IN)  :: YMASS(:,:,:)

	# 	! Grid box surface area for mass flux diag [m2]
	# 	REAL(fp),  INTENT(IN)  :: AREA_M2(JM)

	# 	LOGICAL, INTENT(IN)    :: FILL    ! Fill negatives ?
	# !
	# ! !INPUT/OUTPUT PARAMETERS:
	# !
	# 	! V-wind (m/s) at mid-time-level (t=t+dt/2)
	# 	REAL(fp),  INTENT(INOUT) :: v(:,:,:)

	# 	! surface pressure at current time
	# 	REAL(fp),  INTENT(INOUT) :: ps1(IM, JFIRST:JLAST)

	# 	! surface pressure at future time=t+dt
	# 	REAL(fp),  INTENT(INOUT) :: ps2(IM, JFIRST:JLAST)

	# 	! Diagnostics state object
	# 	TYPE(ChmState), INTENT(INOUT) :: State_Chm
	# 	TYPE(DgnState), INTENT(INOUT) :: State_Diag
	# !
	# ! !OUTPUT PARAMETERS:
	# !
	# 	! "Predicted" surface pressure [hPa]
	# 	REAL(fp),  INTENT(OUT)   :: ps(IM,JFIRST:JLAST)

	# !
	# ! !DEFINED PARAMETERS:
	
	advec_consrv_opt = 2 # 2 = floating pressure
	cross = true
	
	jn = zeros(Integer, km)
	js = zeros(Integer, km)
	
	dap = zeros(AbstractFloat, km)
	dbk = zeros(AbstractFloat, km)
	
	cx = zeros(AbstractFloat, im, (jfirst - ng):(jlast + ng), km) # E-W CFL # on C-grid
	cy = zeros(AbstractFloat, im, jfirst:(jlast + mg), km) # N-S CFL # on C-grid
	delp1 = zeros(AbstractFloat, im, jm, km)
	delp2 = zeros(AbstractFloat, im, jm, km)
	delpm = zeros(AbstractFloat, im, jm, km)
	pu = zeros(AbstractFloat, im, jm, km)
	dpi = zeros(AbstractFloat, im, jm, km)
	geofac = zeros(AbstractFloat, jm) # geometrical factor for meridional advection; geofac uses correct spherical geometry, and replaces RGW_25. (ccc, 4/1/09)
	
	dps_ctm = zeros(AbstractFloat, im,jm)
	ua = zeros(AbstractFloat, im, jm, km)
	va = zeros(AbstractFloat, im, jm, km)
	wz = zeros(AbstractFloat, im, jm, km)
	dq1 = zeros(AbstractFloat, im, (jfirst - ng):(jlast + ng), km)

	# qqu, qqv, adx and ady are now 2d arrays for parallelization purposes. (ccc, 4/1/08)
	qqu = zeros(AbstractFloat, im, jm)
	qqv = zeros(AbstractFloat, im, jm)
	adx = zeros(AbstractFloat, im, jm)
	ady = zeros(AbstractFloat, im, jm)

	# fx, fy, fz and qtp are now 4D arrays for parallelization purposes. (ccc, 4/1/09)
	fx = zeros(AbstractFloat, im, jm, km, nq)
	fy = zeros(AbstractFloat, im, jm + 1, km, nq) # one more for edges
	fz = zeros(AbstractFloat, im, jm, km, nq)


	# 	# Add pointer to avoid array temporary in call to FZPPM (bmy, 6/5/13)
	# 	REAL(fp),  POINTER :: q_ptr(:,:,:)

	# Add definition of j1p and j2p for enlarge polar cap. (ccc, 11/20/08)
	j1p = 3
	j2p = jm - j1p + 1

	# MACROS
	#ifdef TOMAS
		# For TOMAS microphysics: zero out UA and VA.
		# Segregate this block from the code with an #ifdef block. We can"t bring this into the standard GEOS-Chem yet, since that will make it hard to compare benchmark results to prior versions.  When we do bring this change into the standard code, we will have to benchmark it. (sfarina, bmy, 5/30/13)
		for ik = 1:km, ij = 1:jm, il = 1:im
			va[il, ij, ik] = 0.0
			ua[il, ij, ik] = 0.0
		end
	#endif

	# Average surf. pressures in the polar cap. (ccc, 11/20/08)
	# TODO:
	average_press_poles!(area_m2, ps1, 1, im, 1, jm, 1, im, 1, jm)
	# TODO:
	average_press_poles!(area_m2, ps2, 1, im, 1, jm, 1, im, 1, jm)
	
	# Calculation of some geographic factors. (ccc, 11/20/08)
	rj2m1 = jm - 1
	dp = π / rj2m1

	for ij = 1:jm
		geofac[ij] = dp / (2.0 * area_m2[ij] / (sum(area_m2) * im) * im)
	end

	geofac_pc = dp / (2.0 * (sum(area_m2[1:2]) / (sum(area_m2) * im)) * im)

	if tpcore_fvdas_first
		tpcore_fvdas_first = false
		# TODO:
		set_lmts!(tpcore_fvdas_ilmt, tpcore_fvdas_jlmt, tpcore_fvdas_klmt, im, jm, iord, jord, kord)
	end

	# Pressure calculations. (ccc, 11/20/08)
	for ik = 1:km
		dap[ik] = ak[ik + 1] - ak[ik]
		dbk[ik] = bk[ik + 1] - bk[ik]
	end

	# NOTE: Translate parallel for below to parallel loop in Julia
	# !$OMP PARALLEL DO        &
	# !$OMP DEFAULT( SHARED  ) &
	# !$OMP PRIVATE( IK, IQ, q_ptr )
	for ik = 1:km
		# TODO:
		set_press_terms!(dap[ik], dbk[ik], ps1, ps2, delp1[:, :, ik], delpm[:, :, ik], pu[:, :, ik], 1, jm, 1, im, 1, jm, j1p, j2p, 1, im, 1, jm)
				
		# ...intent(in)  dap - difference in ai across layer (mb)
		# ...intent(in)  dbk - difference in bi across layer (mb)
		# ...intent(in)  pres1 - surface pressure at t1 (mb)
		# ...intent(in)  pres2 - surface pressure at t1+tdt (mb)
		# ...intent(out) delp1 - pressure thickness at t1 (mb)
		# ...intent(out) delpm - pressure thickness at t1+tdt/2 (mb)
		# ...intent(out) pu - pressure at edges of box for "u" (mb)

		if j1p != 1 + 1
			for iq = 1:nq

				# TODO: translate to Julia
				# TODO: Translate % to Julia
				# q_ptr => State_Chm%Species(iq)%Conc(:,:,km:1:-1)
				
				average_const_poles!(dap[ik], dbk[ik], area_m2, ps1, q_ptr[:, :, ik], 1, jm, im, 1, im, 1, jm, 1, im, 1, jm)
				
				# TODO: translate to Julia
				# q_ptr => NULL()
			end
		end
		
		# TODO:
		calc_courant!(cose, delpm[:, :, ik], pu[:, :, ik], xmass[:, :, ik], ymass[:, :, ik], cx[:, :, ik], cy[:, :, ik], j1p, j2p, 1, jm, 1, im, 1, jm, 1, im, 1, jm)

		# TODO:
		calc_divergence!(true, geofac_pc, geofac, dpi[:, :, ik], xmass[:, :, ik], ymass[:, :, ik], j1p, j2p, 1, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm)
		
		# TODO:
		set_cross_terms!(cx[:, :, ik], cy[:, :, ik], ua[:, :, ik], va[:, :, ik], j1p, j2p, 1, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm, cross)
	end

	dps_ctm[:, :] .= sum(dpi[:, :, :])

	# TODO:
	calc_vert_mass_flux!(dbk, dps_ctm, dpi, wz, 1, im, 1, jm, 1, km)

	# .sds2.. have all mass flux here: east-west(xmass),
	#         north-south(ymass), vertical(wz)
	# .sds2.. save omega (vertical flux) as diagnostic
	
	# TODO:
	set_jn_js!(jn, js, cx, 1, im, 1, jm, 1, jm, j1p, j2p, 1, im, 1, jm, 1, km)
	
	if advec_consrv_opt == 0
		# NOTE: Translate parallel for below to parallel loop in Julia
		# !$OMP PARALLEL DO           &
		# !$OMP DEFAULT( SHARED     ) &
		# !$OMP PRIVATE( IK, IJ, IL )
		for ik = 1:km, ij = 1:jm, il = 1:im
			delp2[il, ij, ik] = dap[ik] + (dbk[ik] * (ps1[il, ij] + dps_ctm[il, ij]))
		end
	elseif advec_consrv_opt == 1 || advec_consrv_opt == 2
		# NOTE: Translate parallel for below to parallel loop in Julia
		# !$OMP PARALLEL DO           &
		# !$OMP DEFAULT( SHARED     ) &
		# !$OMP PRIVATE( IK, IJ, IL )
		for ik = 1:km, ij = 1:jm, il = 1:im
			delp2[il, ij, ik] = dap[ik] + (dbk[ik] * ps2[il, ij])
		end
	end

	# Calculate surf. pressure at t+dt. (ccc, 11/20/08)
	ps = ak[1] + sum(delp2)

	# > For time optimization : we parallelize over tracers and we loop over the levels outside horizontal transport functions. (ccc, 4/1/09)
	# > Also zeroed PRIVATE variables within the loop, and set jn(ik) and js(ik) to PRIVATE loop variables.  This seems to avoid small diffs in output.
	# > -- Bob Yantosca (04 Jan 2022)

	# NOTE: Translate parallel for below to parallel loop in Julia
	# !$OMP PARALLEL DO                                                     &
	# !$OMP DEFAULT( SHARED                                               ) &
	# !$OMP PRIVATE( iq, dq1, ik, adx, ady, q_ptr, qqu, qqv, north, south )
	for iq = 1:nq
		# TODO: figure out how to translate this to Julia
		# q_ptr => state_chm%species(iq)%conc(:,:,km:1:-1)

		# Zero 3-D arrays for each species
		dq1 = 0.0

		for ik = 1:km
			# Zero PRIVATE variables for safety"s sake
			adx = 0.0
			ady = 0.0
			qqu = 0.0
			qqv = 0.0

			# Northernmost and southernmost latitude indices by level
			north = jn[ik]
			south = js[ik]

			# .sds.. convert to "mass"
			dq1[:, :, ik] .= q_ptr[:, :, ik] * delp1[:, :, ik]

			# TODO:
			calc_advec_cross_terms!(north, south, q_ptr[:, :, ik], qqu, qqv, ua[:, :, ik], va[:, :, ik], j1p, j2p, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm, cross)
			
			# .sds.. notes on arrays
			#   q   (in)  - species mixing ratio
			#   qqu (out) - concentration contribution from E-W
			#                advection cross terms(mixing ratio)
			#   qqv (out) - concentration contribution from N-S
			#                advection cross terms(mixing ratio)
			#   ua  (in)  - average of Courant numbers from il and il+1
			#   va  (in)  - average of Courant numbers from ij and ij+1

			# Add advective form E-W operator for E-W cross terms.

			# TODO:
			xadv_dao2!(2, north, south, adx, qqv, ua[:, :, ik], 1, im, 1, jm, 1, jm, j1p, j2p, 1, im, 1, jm)

			# .sds notes on output arrays
			#   adx (out)- cross term due to E-W advection (mixing ratio)
			#   qqv (in) - concentration contribution from N-S
			#              advection (mixing ratio)
			#   ua  (in) - average of Courant numbers from il and il+1
			# .sds

			# Add advective form N-S operator for N-S cross terms.

			# TODO:
			yadv_dao2(2, ady, qqu, va[:, :, ik], 1, im, 1, jm, j1p, j2p, 1, im, 1, jm, 1, im, 1, jm)

			# .sds notes on output arrays
			#   ady (out)- cross term due to N-S advection (mixing ratio)
			#   qqu (in) - concentration contribution from N-S advection
			#              (mixing ratio)
			#   va  (in) - average of Courant numbers from il and il+1
			# .sds
			# 
			# .bmy notes: use a polar cap of 2 boxes (i.e. the "2" as
			#  the first argument to YADV_DAO2.  The older TPCORE only had
			#  a polar cap of 1 box (just the Pole itself).  Claire figured
			#  this out.  (bmy, 12/11/08)
			# ... update constituent array qq1 by adding in cross terms
			#            - use in fzppm
				
			q_ptr[:, :, ik] = q_ptr[:, :, ik] + ady + adx

			# TODO:
			xtp!(tpcore_fvdas_ilmt, north, south, pu[:, :, ik],  cx[:, :, ik], dq1[:, :, ik], qqv, xmass[:, :, ik], fx[:, :, ik, iq], j1p, j2p, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm, iord)

			# .sds notes on output arrays
			#   pu  (in)    - pressure at edges in "u" (mb)
			#   crx (in)    - Courant number in E-W direction
			#   dq1 (inout) - species density (mb) - updated with the E-W flux
			#                 fx in Xtp)
			#   qqv (inout) - continue oncentration contribution from N-S advection
			#                 (mixing ratio)
			#   xmass(in)   - horizontal mass flux in E-W direction (mb)
			#   fx  (out)   - species E-W mass flux
			# .sds
			
			# TODO:
			ytp!(tpcore_fvdas_jlmt, geofac_pc, geofac, cy[:, :, ik], dq1[:, :, ik], qqu, qqv, ymass[:, :, ik], fy[:, :, ik, iq], j1p, j2p, 1, im, 1, jm, im, 1, im, 1, jm, 1, im, 1, jm, jord)

			# .sds notes on output arrays
			#   cy (in)     - Courant number in N-S direction
			#   dq1 (inout) - species density (mb) - updated with the N-S flux
			#                 (fy in Ytp)
			#   qqu (in)    - concentration contribution from E-W advection
			#                 (mixing ratio)
			#   qqv (inout) - concentration contribution from N-S advection
			#                 (mixing ratio)
			#   ymass(in)   - horizontal mass flux in E-W direction (mb)
			#   fy  (out)   - species N-S mass flux (need to mult by geofac)
			# .sds
		end # ik
		
		# TODO:
		fzppm!(tpcore_fvdas_klmt, delp1, wz, dq1, q_ptr, fz[:, :, :, iq], j1p, 1, jm, 1, im, 1, jm, im, km, 1, im, 1, jm, 1, km)

		# .sds notes on output arrays
		#    wz  (in) : vertical mass flux
		#    dq1 (inout) : species density (mb)
		#    q (in) : species concentration (mixing ratio)
		# .sds

		if fill
			# TODO:
			qckxyz!(dq1, j1p, j2p, 1, jm, 1, im, 1, jm, 1, im, 1, jm, 1, km)
		end

		q_ptr[:, :, :] .= dq1 / delp2

		if j1p != 2
			q_ptr[:, 2, :] .= q_ptr[:, 1, :]
			q_ptr[:, jm - 1, :]  .= q_ptr[:, jm, :]
		end

		# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
		# Set tracer concentration to a small positive number if concentration is negative. Negative concentration may occur at the poles. This is an issue that should be looked into in the future. (ewl, 6/30/15)
			
		# TODO: Tranlate this to Julia
		# where ( q_ptr < 0.0_fp )
		# 	q_ptr = 1.0e-26_fp
		# end where

		# TODO: Tranlate this to Julia
		# q_ptr => NULL()
	end
	
	# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
	# HISTORY (aka netCDF diagnostics)
	# E/W flux of advected species [kg/s] (ccarouge 12/2/08)
	# The unit conversion is:
	# Mass/time = (P diff in grid box) * 100/1 * 1/g * (area of grid box area_m2) * 1/s
	# k/g = hPa/1 * Pa/hPa * (s^2)/m * (m^2)/1 * 1/ΔT
	
	# TODO: Translate % to Julia
	# if state_diag%archive_advfluxzonal
	# 	# Zero netCDF diagnostic array
	# 	state_diag%advfluxzonal = 0.0

	# 	# Calculate fluxes for diag. (ccc, 11/20/08)
	# 	js2g0  = max( j1p, jfirst ) # No ghosting
	# 	jn2g0  = min( j2p, jlast  ) # No ghosting

	# 	# Loop over diagnostic slots
	# 	# NOTE: Translate parallel for below to parallel loop in Julia
	# 	# !$OMP PARALLEL DO                           &
	# 	# !$OMP DEFAULT( SHARED                     ) &
	# 	# !$OMP PRIVATE( S, IQ, K, J, I, Kflip )
	# 	# TODO: Translate % to Julia
	# 	# for s = 1:state_diag%map_advfluxzonal%nslots
	# 	# 	# Get the advectId from the slotId
	# 	# 	# iq = state_diag%map_advfluxzonal%slot2id(s)
	# 	# 	# Loop over grid boxes
	# 	# 	for k = 1:km, j = js2g0:jn2g0, i = 1:im
	# 	# 		# Units: [kg/s]
	# 	# 		# But consider changing to area-independent units [kg/m2/s]
	# 	# 		kflip = km - k + 1 # flip vert
				
	# 	# 		# state_diag%advfluxzonal(i,j,kflip,s) = fx(i,j,k,iq) * area_m2(j) * g0_100 / dt
	# 	# 	end
	# 	# end
	# end

	# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
	# HISTORY (aka netCDF diagnostics)
	# N/S flux of tracer [kg/s] (bdf, bmy, 9/28/04, ccarouge 12/12/08)
	# NOTE, the unit conversion is the same as desciribed above for the E-W diagnostics.  The geometrical factor was already applied to fy in Ytp. (ccc, 4/1/09)
	
	
	# TODO: Translate % to Julia
	# if state_diag%archive_advfluxmerid
	# 	# Zero netCDF diagnostic array
	# 	state_diag%advfluxmerid = 0.0

	# 	# NOTE: Translate parallel for below to parallel loop in Julia
	# 	# !$OMP PARALLEL DO                           &
	# 	# !$OMP DEFAULT( SHARED                     ) &
	# 	# !$OMP PRIVATE( S, IQ, K, J, I, Kflip )
	# 	for s = 1:state_diag%map_advfluxmerid%nslots
	# 		# Get the advectId from the slotId
	# 		iq = state_diag%map_advfluxmerid%slot2id(s)

	# 		# Loop over grid boxes
	# 		for k = 1:km, j = 1:jm, i = 1:im
	# 			# Compute mass flux [kg/s]
	# 			# Units: [kg/s]
	# 			# But consider changing to area-independent units [kg/m2/s]
	# 			kflip = km - k + 1  ! flip vert
	# 			state_diag%advfluxmerid(i,j,kflip,s) = fy(i,j,k,iq) * area_m2(j) * g0_100 / dt
	# 		end
	# 	end
	# end

	# ! MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
	# ! HISTORY (aka netCDF diagnostics)
	# Up/down flux of tracer [kg/s] (bmy, bdf, 9/28/04, ccarouge 12/2/08)
	# The vertical transport done in qmap.  We need to find the difference in order to to interpret transport.
	# Break up diagnostic into up & down fluxes using the surface boundary conditions.  Start from top down (really surface up for flipped TPCORE)
	# By construction, MASSFLUP is flux into the bottom of the box. The flux at the bottom of KM (the surface box) is not zero by design. (phs, 3/4/08)
	
	# TODO: Translate % to Julia
	# if State_Diag%Archive_AdvFluxVert
	# 	# Zero netCDF diagnostic array
	# 	State_Diag%AdvFluxVert .= 0.0

	# 	# NOTE: Translate parallel for below to parallel loop in Julia
	# 	# !$OMP PARALLEL DO                           &
	# 	# !$OMP DEFAULT( SHARED                     ) &
	# 	# !$OMP PRIVATE( S, IQ, K, J, I, Kflip )
	# 	for s = 1:State_Diag%Map_AdvFluxVert%nSlots
	# 		# Get the advectId from the modelId
	# 		iq = State_Diag%Map_AdvFluxVert%slot2Id(S)

	# 		# Loop over grid boxes
	# 		for k = 1:km, j = 1:jm, i = 1:im
	# 			# Ilya Stanevic (stanevic@atmosp.physics.utoronto.ca) says that using FZ for the ND26 diagnostic should be correct. He writes:
	# 			# > To be safe you can use FZ variable from Fzppm. That is the real vertical species mass. And it is zero at the surface.
	# 			# > To be clear, Fz is present only in the tpcore for low resolution GLOBAL model (4x5, 2x2.5). Nested model has different way to calculate vertical advection and there is no such thing as FZ. Therefore, we have to calculate the species mass difference in the box before and after vertical advection in order to get vertical mass flux.
	# 			# > -- Bob Yantosca (28 Mar 2017)
	# 			# Units: [kg/s]
	# 			# But consider changing to area-independent units [kg/m2/s]
	# 			kflip = km - k + 1  !flip vert
	# 			state_diag%advfluxvert(i,j,kflip,s) = fz(i,j,k,iq) * area_m2(j) * g0_100 / dt
	# 		end
	# 	end
	# end
end

"""
Averages the species concentrations at the Poles when the Polar cap is enlarged. It makes the last two latitudes equal.

## Arguments
- `dap::AbstractFloat` - IN - Pressure difference across layer from (ai * pt) term [hPa]
- `dbk::Integer` - IN - Difference in bi across layer - the dSigma term
- `rel_area::Array{AbstractFloat}` - IN - (ju1:j2) - 
- `pctm1::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - 
- `const1::Matrix{AbstractFloat}` - INOUT - (i1:i2, ju1:j2) - 
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function average_const_poles!(
	dap::AbstractFloat,
	dbk::Integer,
	rel_area::Array{AbstractFloat},
	pctm1::Matrix{AbstractFloat},
	const1::Matrix{AbstractFloat},
	ju1_gl::Integer,
	j2_gl::Integer,
	i2_gl::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer
)::Nothing
	# pressure thickness at North Pole, the psudo-density in a hydrostatic system at t1 (mb)
	delp1n = zeros(AbstractFloat, i1:i2, (j2 - 1):j2)
	# pressure thickness at South Pole, the psudo-density in a hydrostatic system at t1 (mb)
	delp1s = zeros(AbstractFloat, i1:i2,  ju1:(ju1 + 1))

	if ju1 == ju1_gl
		delp1s[i1:i2, ju1:(ju1+1)] = dap + (dbk * pctm1[i1:i2, ju1:(ju1+1)])

		sum1 = 0.0
		sum2 = 0.0
		for il = i1:i2
			sum1 = sum1 + sum(const1[il, ju1:(ju1 + 1)] * delp1s[il, ju1:(ju1 + 1)] * rel_area[ju1:(ju1 + 1)]) / (sum(rel_area) * i2_gl)
			sum2 = sum2 + sum(delp1s[il, ju1:(ju1 + 1)] * rel_area[ju1:(ju1 + 1)]) / (sum(rel_area) * i2_gl)
		end

		meanq = sum1 / sum2

		const1[:, ju1:(ju1 + 1)] .= meanq
	end
	
	if j2 == j2_gl
		delp1n[i1:i2, (j2 - 1):j2] .= dap + (dbk * pctm1[i1:i2, (j2 - 1):j2])

		sum1 = 0.0
		sum2 = 0.0
		for il = i1:i2
			sum1 = sum1 + sum(const1[il, (j2 - 1):j2] * delp1n[il, (j2 - 1):j2] * rel_area[(j2 - 1):j2]) / (sum(rel_area) * i2_gl)
			sum2 = sum2 + sum(delp1n[il, (j2 - 1):j2] * rel_area[(j2 - 1):j2]) / (sum(rel_area) * i2_gl)
		end

		meanq = sum1 / sum2

		const1[:, (j2 - 1):j2] .= meanq
	end
end

"""
Sets the cross terms for E-W horizontal advection.

## Arguments
- `cry::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Courant number in E-W direction
- `crx::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Courant number in N-S direction
- `ua::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Average of Courant numbers from il and il+1
- `va::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Average of Courant numbers from ij and ij+1
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `cross::Bool` - IN - Logical switch.  If CROSS=T then cross-terms will be computed.

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function set_cross_terms!(
	crx::Matrix{AbstractFloat},
	cry::Matrix{AbstractFloat},
	ua::Matrix{AbstractFloat},
	va::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	cross::Bool
)::Nothing
	if !cross
		ua[:, :] .= 0.0
		va[:, :] .= 0.0
	else
		for ij = j1p:j2p
			for il = i1:(i2 - 1)
				ua[il, ij] = 0.5 * (crx[il, ij] + crx[il + 1, ij])
			end
			ua[i2, ij] = 0.5 * (crx[i2, ij] + crx[1, ij])
		end

		for ij = (ju1 + 1):(j2 - 1), il = i1:i2
			va[il, ij] = 0.5 * (cry[il, ij] + cry[il, ij + 1])
		end
		
		# TODO:
		do_cross_terms_pole_i2d2!(cry, va, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
	end
end

"""
Calculates the vertical mass flux.

## Arguments
- `dbk::Array{AbstractFloat}` - IN - (k1:k2) - Difference in bi across layer - the dSigma term
- `dps_ctm::Matrix{AbstractFloat}` - IN - (i1:i2, ju1:j2) -  CTM surface pressure tendency; sum over vertical of dpi calculated from original mass fluxes [hPa]
- `dpi::Array{AbstractFloat, 3}` - IN - (i1:i2, ju1:j2, k1:k2) - Divergence at a grid point; used to calculate vertical motion [mb]
- `wz::Array{AbstractFloat, 3}` - OUT - (i1:i2, ju1:j2, k1:k2) - Large scale mass flux (per time step tdt) in the vertical direction as diagnosed from the hydrostatic relationship [hPa]
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `k1::Integer` - IN - Local min altitude index
- `k2::Integer` - IN - Local max altitude index

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function calc_vert_mass_flux!(
	dbk::Array{AbstractFloat},
	dps_ctm::Matrix{AbstractFloat},
	dpi::Array{AbstractFloat, 3},
	wz::Array{AbstractFloat, 3},
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	k1::Integer,
	k2::Integer
)::Nothing
	# Compute vertical mass flux from mass conservation.

	# NOTE: Translate parallel for below to parallel loop in Julia
	# !$OMP PARALLEL DO       &
	# !$OMP DEFAULT( SHARED ) &
	# !$OMP PRIVATE( IJ, IL )
	for ij = ju1:j2, il = i1:i2
		wz[il, ij, k1] = dpi[il, ij, k1] - (dbk[k1] * dps_ctm[il, ij])

		wz[il, ij, k2] = 0.0
	end

	for ik = (k1 + 1):(k2 - 1)
		# NOTE: Translate parallel for below to parallel loop in Julia
		# !$OMP PARALLEL DO       &
		# !$OMP DEFAULT( SHARED ) &
		# !$OMP PRIVATE( IJ, IL )
		for ij = ju1:j2, il = i1:i2
			wz[il, ij, ik] = wz[il, ij, ik - 1] + dpi[il, ij, ik] - (dbk[ik] * dps_ctm[il, ij])
		end
	end
end

"""
Determines Jn and Js, by looking where Courant number is > 1.

## Arguments
- `jn::Array{AbstractFloat}` - OUT - (k1:k2) - Northward of latitude index = jn; Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `js::Array{AbstractFloat}` - OUT - (k1:k2) - Southward of latitude index = js; Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `crx::Array{AbstractFloat, 3}` - IN - (ilo:ihi, julo:jhi, k1:k2) - Courant number in E-W direction
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `k1::Integer` - IN - Local min altitude index
- `k2::Integer` - IN - Local max altitude index

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Remarks
We cannot parallelize this function because there is a CYCLE statement within the outer loop.

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function set_jn_js!(
	jn::Array{AbstractFloat},
	js::Array{AbstractFloat},
	crx::Array{AbstractFloat, 3},
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	j2p::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	k1::Integer,
	k2::Integer
)::Nothing
	js0 = (j2_gl + 1) / 2
	jn0 = j2_gl - js0 + 1

	jst = max(ju1, j1p)
	jend = min(j2, js0)

	for ik = k1:k2
		js[ik] = j1p

		for ij = jend:-1:jst, il = i1:i2
			if abs(crx[il, ij, ik]) > 1.0
				js[ik] = ij
				
				@goto ikloop1
			end
		end
		@label ikloop1
	end

	jst = max(ju1, jn0)
	jend = min(j2, j2p)

	for ik = k1:k2
		jn[ik] = j2p

		for ij = jst:jend, il = i1:i2
			if abs(crx[il, ij, ik]) > 1.0
				jn[ik] = ij
				
				@goto ikloop2
			end
		end
		@label ikloop2
	end
end

"""
Calculates the advective cross terms.

## Arguments
- `jn::Integer` - IN - Northward of latitude index = jn, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `js::Integer` - IN - Southward of latitude index = js, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `qq1::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Species concentration (mixing ratio)
- `qqu::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Concentration contribution from E-W advection [mixing ratio]
- `qqv::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - concentration contribution from N-S advection [mixing ratio]
- `ua::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Average of Courant numbers from il and il+1
- `va::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Average of Courant numbers from ij and ij+1
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `cross::Bool` - IN - Logical switch: If CROSS=T then cross-terms are being computed

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function calc_advec_cross_terms!(
	jn::Integer,
	js::Integer,
	qq1::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	ua::Matrix{AbstractFloat},
	va::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	cross::Bool
)::Nothing
	qtmp = zeros(AbstractFloat, (-i2 / 3):(i2 + i2 / 3), julo:jhi)

	for ij = julo:jhi
		for i = 1:i2
			qtmp[i, ij] = qq1[i, ij]
		end

		for il = (-i2 / 3):0
			qtmp[il, ij] = qq1[i2 + il, ij]
		end

		for il = (i2 + 1):(i2 + i2 / 3)
			qtmp[il, ij] = qq1[il - i2, ij]
		end
	end

	if !cross
		qqv[:, :] .= qq1[:, :]
		qqu[:, :] .= qq1[:, :]
	else
		qqu[:, :] .= 0.0
		qqv[:, :] .= 0.0

		for ij = j1p:j2p
			if ij <= js || ij >= jn
				# In Polar area, so need to deal with large courant numbers.
				for il = i1:i2
					# !c?
					iu = ua[il, ij]
					riu = iu
					ru = ua[il, ij] - riu
					iu = il - iu

					if ua[il, ij] >= 0.0
						qqu[il, ij] = qtmp[iu, ij] + ru * (qtmp[iu - 1, ij] - qtmp[iu, ij])
					else
						qqu[il, ij] = qtmp[iu, ij] + ru * (qtmp[iu, ij] - qtmp[iu + 1, ij])
					end

					qqu[il, ij] = qqu[il, ij] - qtmp[il, ij]
				end
			else # js < ij < jn
				# Do interior area (use PPM).

				for il = i1:i2
					ril = il
					iu = ril - ua[il, ij]

					qqu[il, ij] = ua[il, ij] * (qtmp[iu, ij] - qtmp[iu + 1, ij])
				end
			end

			for il = i1:i2
				# !c?
				rij = ij
				jv = rij - va[il, ij]

				qqv[il, ij] = va[il, ij] * (qtmp[il, jv] - qtmp[il, jv + 1])
			end
		end

		for ij = ju1:j2, il = i1:i2
			qqu[il, ij] = qtmp[il, ij] + (0.5 * qqu[il, ij])
			qqv[il, ij] = qtmp[il, ij] + (0.5 * qqv[il, ij])
		end
	end
end

"""
Checks for "filling".
	
## Arguments
- `dq1::Array{AbstractFloat, 3}` - INOUT - (ilo:ihi, julo:jhi, k1:k2) - Species density [hPa]
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `k1::Integer` - IN - Local min altitude index
- `k2::Integer` - IN - Local max altitude index

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function qckxyz!(
	dq1::Array{AbstractFloat, 3},
	j1p::Integer,
	j2p::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	k1::Integer,
	k2::Integer
)::Nothing
	fill_diag = false
	
	ip = 0
	
	# Top layer.

	k1p1 = k1 + 1

	# NOTE: Translate parallel for below to parallel loop in Julia
	# !$OMP PARALLEL DO          &
	# !$OMP DEFAULT( SHARED )    &
	# !$OMP PRIVATE( IJ, IL, IP )
	for ij = j1p:j2p, il = i1:i2
		if dq1[il, ij, k1] < 0.0
			ip = ip + 1

			dq1[il, ij, k1p1] = dq1[il, ij, k1p1] + dq1[il, ij, k1]
			dq1[il, ij, k1] = 0.0
		end
	end
	
	for ik = (k1 + 1):(k2 - 1)
		# NOTE: Translate parallel for below to parallel loop in Julia
		# !$OMP PARALLEL DO                         &
		# !$OMP DEFAULT( SHARED )                   &
		# !$OMP PRIVATE( IJ, IL, IP, QUP, QLY, DUP )
		for ij = j1p:j2p, il = i1:i2
			if dq1[il, ij, ik] < 0.0
				ip = ip + 1
				
				# From above.

				qup = dq1[il, ij, ik - 1]
				qly = -dq1[il, ij, ik]
				dup = min(qly, qup)

				dq1[il, ij, ik - 1] = qup - dup
				dq1[il, ij, ik] = dup - qly
				
				# From below.

				dq1[il, ij, ik + 1] = dq1[il, ij, ik + 1] + dq1[il, ij, ik]
				dq1[il, ij, ik] = 0.0
			end
		end
	end

	# Bottom layer.

	sum  = 0.0
	k2m1 = k2 - 1

	# NOTE: Translate parallel for below to parallel sum in Julia
	# NOTE: Sum seems to be not used in the loop below!
	# !$OMP PARALLEL DO                          &
	# !$OMP DEFAULT( SHARED )                    &
	# !$OMP PRIVATE( IJ, IL, IP, QUP, QLY, DUP ) &
	# !$OMP REDUCTION( +:SUM )
	for ij = j1p:j2p, il = i1:i2
		if dq1[il, ij, k2] < 0.0
			ip = ip + 1
			
			# From above.

			qup = dq1[il, ij, k2m1]
			qly = -dq1[il, ij, k2]
			dup = min[qly, qup]

			dq1[il, ij, k2m1] = qup - dup
			
			# From "below" the surface.

			sum = sum + qly - dup

			dq1[il, ij, k2] = 0.0
		end
	end

	# We don"t want to replace zero values by 1e-30. (ccc, 11/20/08)
	#   where ((dq1(i1:i2,j1p:j2p,:) < 1.0d-30)) &
	#           dq1(i1:i2,j1p:j2p,:) = 1.0d-30
end

"""
Sets ilmt, jlmt, klmt. (See REMARKS section of routine Tpcore_FvDas for more info)

## Arguments
- `ilmt::Integer` - OUT - Controls various options in E-W advection
- `jlmt::Integer` - OUT - Controls various options in N-S advection
- `klmt::Integer` - OUT - Controls various options in vertical advection
- `i2_gl::Integer` - IN - Global max latitude index
- `j2_gl::Integer` - IN - Global max longitude index
- `iord::Integer` - IN - Flag to denote E-W transport schemes
- `jord::Integer` - IN - Flag to denote N-S transport schemes
- `kord::Integer` - IN - Flag to denote vertical transport schemes

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function set_lmts!(
	ilmt::Integer,
	jlmt::Integer,
	klmt::Integer,
	i2_gl::Integer,
	j2_gl::Integer,
	iord::Integer,
	jord::Integer,
	kord::Integer
)::Nothing
	j2_glm1 = j2_gl - 1

	# !c?
	if iord <= 0
		if i2_gl >= 144
			ilmt = 0
		elseif i2_gl >= 72
			ilmt = 1
		else
			ilmt = 2
		end
	else
		ilmt = iord - 3
	end

	# !c?
	if jord <= 0
		if j2_glm1 >= 90
			jlmt = 0
		elseif j2_glm1 >= 45
			jlmt = 1
		else
			jlmt = 2
		end
	else
		jlmt = jord - 3
	end

	klmt = max(kord - 3, 0)
end

"""
Sets the pressure terms: delp1, delpm, pu.

## Arguments
- `dap::AbstractFloat` - Pressure difference across layer from (ai * pt) term [hPa]
- `dbk::AbstractFloat` - Difference in bi across layer - the dSigma term
- `pres1::Matrix{AbstractFloat}` - (ilo:ihi, julo:jhi) - Surface pressure at t1 [hPa]
- `pres2::Matrix{AbstractFloat}` - (ilo:ihi, julo:jhi) - Surface pressure at t1+tdt [hPa]
- `delp1::Matrix{AbstractFloat}` - (ilo:ihi, julo:jhi) - Pressure thickness, the pseudo-density in a hydrostatic system at t1 [hPa]
- `delpm::Matrix{AbstractFloat}` - (ilo:ihi, julo:jhi) - Pressure thickness, the pseudo-density in a hydrostatic system at t1+tdt/2 (approximate) [hPa]
- `pu::Matrix{AbstractFloat}` - (ilo:ihi, julo:jhi) - Pressure at edges in "u" [hPa]
- `ju1_gl::Integer` - Global min latitude index
- `j2_gl::Integer` - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function set_press_terms!(
	dap::AbstractFloat,
	dbk::AbstractFloat,
	pres1::Matrix{AbstractFloat},
	pres2::Matrix{AbstractFloat},
	delp1::Matrix{AbstractFloat},
	delpm::Matrix{AbstractFloat},
	pu::Matrix{AbstractFloat},
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	j1p::Integer,
	j2p::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	delp1[:, :] .= dap + (dbk * pres1[:, :])

	delpm[:, :] .= dap + (dbk * 0.5 * (pres1[:, :] + pres2[:, :]))
	
	for ij = j1p:j2p
		pu[1, ij] = 0.5 * (delpm[1, ij] + delpm[i2, ij])
		for il = (i1 + 1):i2
			pu[il, ij] = 0.5 * (delpm[il, ij] + delpm[il - 1, ij])
		end
	end
end

"""
Calculates courant numbers from the horizontal mass fluxes.

## Arguments
- `cose::Array{AbstractFloat}` - IN - (ju1_gl:j2_gl) - Cosine of grid box edges
- `delpm::Matrix{AbstractFloat}` - IN - (ju1_gl:j2_gl) - Pressure thickness, the pseudo-density in a hydrostatic system at t1+tdt/2 (approximate) (mb)
- `pu::Matrix{AbstractFloat}` - IN - (ju1_gl:j2_gl) - pressure at edges in "u"  (mb)
- `xmass::Matrix{AbstractFloat}` - IN - (ju1_gl:j2_gl) - horizontal mass flux in E-W directions [hPa]
- `ymass::Matrix{AbstractFloat}` - IN - (ju1_gl:j2_gl) - horizontal mass flux in N-S directions [hPa]
- `crx::Matrix{AbstractFloat}` - OUT - (ju1_gl:j2_gl) - Courant numbers in E-W directions
- `cry::Matrix{AbstractFloat}` - OUT - (ju1_gl:j2_gl) - Courant numbers in N-S directions
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `ju1_gl::Integer` - IN - Local min latitude index
- `j2_gl::Integer` - IN - Local max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO).
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function calc_courant!(
	cose::Array{AbstractFloat},
	delpm::Matrix{AbstractFloat},
	pu::Matrix{AbstractFloat},
	xmass::Matrix{AbstractFloat},
	ymass::Matrix{AbstractFloat},
	crx::Matrix{AbstractFloat},
	cry::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	crx[:, :] .= 0.0
	cry[:, :] .= 0.0
	
	# Calculate E-W and N-S horizontal mass fluxes.

	for ij = j1p:j2p
		crx[:, ij] = xmass[:, ij] / pu[:, ij]

		cry[:, ij] .= ymass[:, ij] / ((0.5 * cose[ij]) * (delpm[:, ij] + delpm[:, ij - 1]))
	end

	cry[:, j2p + 1] .= ymass[:, j2p + 1] / ((0.5 * cose[j2p + 1]) * (delpm[:, j2p + 1] + delp[:, j2p]))
end

"""
Calculates the divergence.
	
## Arguments
- `do_reduction::Bool` - IN - Set to F if called on root core or T if called by secondary cores
	- (NOTE: This is only for MPI parallelization, for OPENMP it should be F)
- `geofac_pc::AbstractFloat` - IN - Special geometrical factor (geofac) for Polar cap
- `geofac::Array{AbstractFloat}` - IN - Geometrical factor for meridional advection; geofac uses correct spherical geometry, and replaces acosp as the meridional geometrical factor in TPCORE
- `dpi::Matrix{AbstractFloat}` - OUT - Divergence at a grid point; used to calculate vertical motion [hPa]
- `xmass::Matrix{AbstractFloat}` - IN - Horizontal mass flux in E/W directions [hPa]
- `ymass::Matrix{AbstractFloat}` - IN - Horizontal mass flux in N/S directions [hPa]
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function calc_divergence!(
	do_reduction::Bool,
	geofac_pc::AbstractFloat,
	geofac::Array{AbstractFloat},
	dpi::Matrix{AbstractFloat},
	xmass::Matrix{AbstractFloat},
	ymass::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	# Calculate N-S divergence.

	for ij = j1p:j2p
		dpi[:, ij] .= (ymass[:, ij] - ymass[:, ij + 1]) * geofac[ij]

		# Calculate E-W divergence.
		
		for il = i1:(i2 - 1)
			dpi[il, ij] = dpi[il, ij] + xmass[il, ij] - xmass[il + 1, ij]
		end

		dpi[i2, ij] = dpi[i2, ij] + xmass[i2, ij] - xmass[1, ij]
	end

	# TODO:
	do_divergence_pole_sum!(do_reduction, geofac_pc, dpi, ymass, i1_gl, i2_gl, j1p, j2p, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
	
	if j1p != ju1_gl + 1
		# Polar cap enlarged:  copy dpi to polar ring.
		
		dpi[:, ju1 + 1] .= dpi[:, ju1]
		dpi[:, j2 - 1] .= dpi[:, j2]
	end
end

"""
Sets the divergence at the Poles.

## Arguments
- `do_reduction::Bool` - IN - Set to T if called on root core or F if called by secondary cores
	- NOTE: This seems not to be used here....)
- `geofac_pc::AbstractFloat` - IN - Special geometrical factor (geofac) for Polar cap
- `dpi::Matrix{AbstractFloat}` - OUT - (i1:i2, ju1:j2) - Divergence at a grid point; used to calculate vertical motion [hPa]
- `ymass::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Horizontal mass flux in N-S direction [hPa]
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
	
## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_divergence_pole_sum!(
	do_reduction::Bool,
	geofac_pc::AbstractFloat,
	dpi::Matrix{AbstractFloat},
	ymass::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	j1p::Integer,
	j2p::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	ri2 = i2_gl
	
	if ju1 == ju1_gl
		sumsp = 0.0
		for il = i1:i2
			sumsp = sumsp + ymass[il, j1p]
		end

		mean_sp = -sumsp / ri2 * geofac_pc

		for il = i1:i2
			dpi[il, ju1] = mean_sp
		end
	end
		
	if j2 == j2_gl
		sumnp = 0.0

		for il = i1:i2
			sumnp = sumnp + ymass[il, j2p + 1]
		end

		mean_np = sumnp / ri2 * geofac_pc

		for il = i1:i2
			dpi[il, j2] = mean_np
		end
	end
end

"""
Sets "va" at the Poles.

## Arguments
- `cry::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Courant number in N-S direction
- `va::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Average of Courant numbers from ij and ij+1
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude indices at the edge of the South polar cap; j1p=ju1_gl+1 for a polar cap of 1 latitude band
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_cross_terms_pole_i2d2!(
	cry::Matrix{AbstractFloat},
	va::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	i2d2 = i2_gl / 2
	
	if j1p == ju1_gl + 1
		# Polar Cap NOT Enlarged: Get cross terms for N-S horizontal advection.
		
		if ju1 == ju1_gl
			for il = i1:i2d2
				va[il, ju1] = 0.5 * (cry[il, ju1 + 1] - cry[il + i2d2, ju1 + 1])
				va[il + i2d2, ju1] = -va[il, ju1]
			end
		end

		if j2 == j2_gl
			for il = i1:i2d2
				va[il, j2] = 0.5 * (cry[il, j2] - cry[il + i2d2, j2 - 1])
				va[il + i2d2, j2] = -va[il, j2]
			end
		end
	end
end

"""
The advective form E-W operator for computing the adx (E-W) cross term.

## Arguments
- `iad::Integer` - IN
	- if iad = 1, use 1st order accurate scheme
	- if iad = 2, use 2nd order accurate scheme
- `jn::Integer` - IN - Northward of latitude index = jn, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `js::Integer` - IN - southward of latitude index = js, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `adx::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Cross term due to E-W advection [mixing ratio]
- `qqv::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Concentration contribution from N-S advection [mixing ratio]
- `ua::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Average of Courant numbers from il and il+1
- `ihi::Integer` - IN - Local max longitude index
- `ilo::Integer` - IN - Local min longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `ju1_gl::Integer` - IN - Local min latitude index
- `j2_gl::Integer` - IN - Local max latitude index
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function xadv_dao2!(
	iad::Integer,
	jn::Integer,
	js::Integer,
	adx::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	ua::Matrix{AbstractFloat},
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	j2p::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	qtmp = zeros(AbstractFloat, (i2 / 3):(i2 + i2 / 3), julo:jhi)

	# Zero output array
	adx = 0
	for ij = julo:jhi
		for il = 1:i2
			qtmp[il, ij] = qqv[il, ij]
		end

		for il = (-i2 / 3):0
			qtmp[il, ij] = qqv[i2 + il, ij]
		end

		for il = (i2 + 1):(i2 + i2 / 3)
			qtmp[il, ij] = qqv[il - i2, ij]
		end
	end
			
	if iad == 1
		# 1st order.
		
		for ij = j1p:j2p
			if ij <= js || ij >= jn
				# In Polar area.

				for il = i1:i2
					iu = ua[il, ij]
					riu = iu
					ru = ua[il, ij] - riu
					iu = il - iu

					if ua[il, ij] >= 0.0
						rdiff = qtmp[iu - 1, ij] - qtmp[iu, ij]
					else
						rdiff = qtmp[iu, ij] - qtmp[iu + 1, ij]
					end

					adx[il, ij] = (qtmp[iu, ij] - qtmp[il, ij]) + (ru * rdiff)
				end
			else # js < ij < jn
				# Eulerian upwind.
				
				for il = i1:i2
					ril = il
					iu = ril - ua[il, ij]
					
					adx[il, ij] = ua[il, ij] * (qtmp[iu, ij] - qtmp[iu + 1, ij])
				end
			end
		end
	elseif iad == 2
		for ij = j1p:j2p
			if ij <= js || ij >= jn
				# In Polar area.
				
				for il = i1:i2
					iu = round(ua[il, ij])
					riu = iu
					ru = riu - ua[il, ij]
					iu = il - iu

					a1 = 0.5 * (qtmp[iu + 1, ij] + qtmp[iu - 1, ij]) - qtmp[iu, ij]

					b1 = 0.5 * (qtmp[iu + 1, ij] - qtmp[iu - 1, ij])

					c1 = qtmp[iu, ij] - qtmp[il, ij]

					adx[il, ij] = (ru * ((a1 * ru) + b1)) + c1
				end
			else # js < ij < jn
				# Eulerian upwind.

				for il = i1:i2
					iu = round(ua[il, ij])
					riu = iu
					ru = riu - ua[il, ij]
					iu = il - iu

					a1 = 0.5 * (qtmp[iu + 1, ij] + qtmp[iu - 1, ij]) - qtmp[iu, ij]

					b1 = 0.5 * (qtmp[iu + 1, ij] - qtmp[iu - 1, ij])

					c1 = qtmp[iu, ij] - qtmp[il, ij]

					adx[il, ij] = (ru * ((a1 * ru) + b1)) + c1
				end
			end
		end
	end
	
	if ju1 == ju1_gl
		adx[i1:i2, ju1] .= 0.0

		if j1p != ju1_gl + 1
			adx[i1:i2, ju1 + 1] .= 0.0
		end
	end

	if j2 == j2_gl
		adx[i1:i2, j2] .= 0.0

		if j1p != ju1_gl + 1
			adx[i1:i2, j2 - 1] .= 0.0
		end
	end
end

"""
The advective form N-S operator for computing the ady (N-S) cross term.

## Arguments
- `iad::Integer` - IN
	- If iad = 1, use 1st order accurate scheme
	- If iad = 2, use 2nd order accurate scheme
- `ady::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Cross term due to N-S advection (mixing ratio)
- `qqu::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Concentration contribution from E-W advection [mixing ratio]
- `va::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Average of Courant numbers from ij and ij+1
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function yadv_dao2!(
	iad::Integer,
	ady::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	va::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	j2p::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	# We may need a small ghost zone depending on the polar cap used
	qquwk = zeros(AbstractFloat, ilo:ihi, (julo - 2):(jhi + 2))

	# Zero output array
	ady = 0

	# Make work array
	for ij = julo:jhi
		qquwk[:, ij] .= qqu[:, ij]
	end

	# This routine creates a ghost zone in latitude in case of not enlarged polar cap (ccc, 11/20/08)
	# TODO:
	do_yadv_pole_i2d2!(qqu, qquwk, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
	
	if iad == 1
		# 1st order.
		for ij = (j1p - 1):(j2p + 1), il = i1:i2
			# !c?
			rij = ij
			jv = rij - va[il, ij]

			ady[il, ij] = va[il, ij] * (qquwk[il, jv] - qquwk[il, jv + 1])
		end
	elseif iad == 2
		for ij = (j1p - 1):(j2p + 1), il = i1:i2
			# c?
			jv  = round(va[il, ij])
			rjv = jv
			rv  = rjv - va[il, ij]
			jv  = ij - jv

			a1 = 0.5 * (qquwk[il, jv + 1] + qquwk[il, jv - 1]) - qquwk[il, jv]

			b1 = 0.5 * (qquwk[il, jv + 1] - qquwk[il, jv - 1])

			c1 = qquwk[il, jv] - qquwk[il, ij]

			ady[il, ij] = (rv * ((a1 * rv) + b1)) + c1
		end
	end
	
	# TODO:
	do_yadv_pole_sum!( ady, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
end

"""
Sets "qquwk" at the Poles.

## Arguments
- `qqu::Integer` - IN - (ilo:ihi, julo:jhi) - concentration contribution from E-W advection [mixing ratio]
- `qquwk::Integer` - IN - (ilo:ihi, julo-2:jhi+2) - qqu working array [mixing ratio]
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude indices at the edges of the South polar cap. j1p=ju1_gl+1 for a polar cap of 1 latitude band
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
# 05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_yadv_pole_i2d2!(
	qqu::Matrix{AbstractFloat},
	qquwk::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	i2d2 = i2_gl / 2
	
	if j1p == ju1_gl + 1
		# Polar Cap NOT Enlarged.

		if ju1 == ju1_gl
			for il = i1:i2d2, inb = 1:2
				qquwk[il, ju1 - inb] = qqu[il + i2d2, ju1 + inb]
				qquwk[il + i2d2, ju1 - inb] = qqu[il, ju1 + inb]
			end
		end

		if j2 == j2_gl
			for il = i1:i2d2, inb = 1:2
				qquwk[il, j2 + inb] = qqu[il + i2d2, j2 - inb]
				qquwk[il + i2d2, j2 + inb] = qqu[il, j2 - inb]
			end
		end
	end
end

"""
Sets the cross term due to N-S advection at the Poles.

## Arguments
- `ady::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Cross term due to N-S advection (mixing ratio)
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude index at the edge of the South polar cap; j1p=ju1_gl+1; for a polar cap of 1 latitude band
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_yadv_pole_sum!(
	ady::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	# Test if we are using extended polar caps (i.e. the S pole and next N latitude and N. Pole and next S latitude).  Do this outside the loops. (bmy, 12/11/08)
	is_ext_polar_cap = j1p == ju1_gl + 2

	# South Pole

	sumsp = 0.0
	sumnp = 0.0

	if is_ext_polar_cap
		# For a 2-latitude polar cap (S. Pole + next Northward latitude)
		for il = i1:i2
			sumsp = sumsp + ady[il, ju1 + 1]
			sumnp = sumnp + ady[il, j2 - 1]
		end
	else
		# For a 1-latitude polar cap (S. Pole only)
		for il = i1:i2
			sumsp = sumsp + ady[il, ju1]
			sumnp = sumnp + ady[il, j2]
		end
	end

	sumsp = sumsp / i2_gl
	sumnp = sumnp / i2_gl
				
	if is_ext_polar_cap 
		# For a 2-latitude polar cap (S. Pole + next Northward latitude)
		for il = i1:i2
			ady[il, ju1 + 1] = sumsp
			ady[il, ju1] = sumsp
			ady[il, j2 - 1] = sumnp
			ady[il, j2] = sumnp
		end
	else
		# For a 1-latitude polar cap (S. Pole only)
		for il = i1:i2
			ady[il, ju1] = sumsp
			ady[il, j2] = sumnp
		end
	end
end

"""
Does horizontal advection in the E-W direction.

## Arguments
- `ilmt::Integer` - IN - Controls various options in E-W advection
- `jn::Integer` - IN - Northward of latitude index = jn, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `js::Integer` - IN - Southward of latitude index = js, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `pu::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - pressure at edges in "u" [hPa]
- `crx::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Courant number in E-W direction
- `dq1::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Species density [hPa]
- `qqv::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Concentration contribution from N-S advection [mixing ratio]
- `xmass::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Horizontal mass flux in E-W direction [hPa]
- `fx::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - E-W flux [mixing ratio]
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `iord::Integer` - IN - Option for E-W transport scheme.  See module header for more info.

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function xtp!(
	ilmt::Integer,
	jn::Integer,
	js::Integer,
	pu::Matrix{AbstractFloat},
	crx::Matrix{AbstractFloat},
	dq1::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	xmass::Matrix{AbstractFloat},
	fx::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	iord::Integer
)::Nothing
	isav = zeros(Integer, i1:i2)
	dcx = zeros(AbstractFloat, (-i2 / 3):(i2 + i2 / 3), julo:jhi)
	qtmp = zeros(AbstractFloat, (-i2 / 3):(i2 + i2 / 3), julo:jhi)
	fx[:, :] .= 0.0

	imp = i2 + 1

	# NOTE: these loops do not parallelize well (bmy, 12/5/08)

	# Populate qtmp
	for il = i1:i2
		qtmp[il, :] .= qqv[il, :]
	end

	for il = (-i2 / 3):0
		qtmp[il, :] = qqv[i2 + il, :]
	end

	for il = (i2 + 1):(i2 + i2 / 3)
		qtmp[il, :] = qqv[il - i2, :]
	end

	if iord != 1
		qtmp[i1 - 1, :] .= qqv[i2, :]
		qtmp[i1 - 2, :] .= qqv[i2 - 1, :]
		qtmp[i2 + 1, :] .= qqv[i1, :]
		qtmp[i2 + 2, :] .= qqv[i1 + 1, :]
		
		# TODO:
		xmist!(dcx, qtmp, j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
	end
	
	jvan = max(1, j2_gl / 18)
	
	for ij = j1p:j2p
		if ij > js && ij < jn
			# Do horizontal Eulerian advection in the E-W direction.
			
			if iord == 1 || ij == j1p || ij == j2p
				for il = i1:i2
					ril = il
					iu = ril - crx[il, ij]

					fx[il, ij] = qtmp[iu, ij]
				end
			else
				if iord == 2 || ij <= j1p + jvan || ij >= j2p - jvan
					for il = i1:i2
						ril = il
						iu = ril - crx[il, ij]

						fx[il, ij] = qtmp[iu, ij] + (dcx[iu, ij] * (sign(crx[il, ij]) - crx[il, ij]))
					end
				else
					# TODO:
					fxppm!(ij, ilmt, crx, dcx, fx, qtmp, -i2/3, i2+i2/3, julo, jhi, i1, i2)
					# qtmp (inout) - can be updated
				end
			end
			
			for il = i1:i2
				fx[il, ij] = fx[il, ij] * xmass[il, ij]
			end
		else
			# Do horizontal Conservative (flux-form) Semi-Lagrangian advection in the E-W direction (van Leer at high latitudes).
			if iord == 1 || ij == j1p  ij == j2p
				for il = i1:i2
					ic = crx[il, ij]
					isav[il] = il - ic
					ril = il
					iu = ril - crx[il, ij]
					ric = ic
					rc = crx[il, ij] - ric

					fx[il, ij] = rc * qtmp[iu, ij]
				end
			else
				for il = i1:i2
					ic = crx[il, ij]
					isav[il] = il - ic
					ril = il
					iu = ril - crx[il, ij]
					ric = ic
					rc = crx[il, ij] - ric

					fx[il, ij] = rc * (qtmp[iu, ij] + (dcx[iu, ij] * (sign(rc) - rc)))
				end
			end

			for il = i1:i2
				if crx[il, ij] > 1.0
					for ix = isav[il]:(il - 1)
						fx[il, ij] = fx[il, ij] + qtmp[ix, ij]
					end
				elseif crx[il, ij] < -1.0
					for ix = il:(isav[il] - 1)
						fx[il, ij] = fx[il, ij] - qtmp[ix, ij]
					end
				end
			end

			for il = i1:i2
				fx[il, ij] = pu[il, ij] * fx[il, ij]
			end
		end
	end
			
	# NOTE: This loop does not parallelize well (bmy, 12/5/08)
	for ij = j1p:j2p
		for il = i1:(i2 - 1)
			dq1[il, ij] = dq1[il, ij] + (fx[il, ij] - fx[il + 1, ij])
		end
		dq1[i2, ij] = dq1[i2, ij] + (fx[i2, ij] - fx[i1, ij])
	end
end

"""
Computes the linear tracer slope in the E-W direction. It uses the Lin et. al. 1994 algorithm.

## Arguments
- `dcx::Matrix{AbstractFloat}` - (-i2/3:i2+i2/3, julo:jhi) - Slope of concentration distribution in E-W direction [mixing ratio]
- `qqv::Matrix{AbstractFloat}` - (-i2/3:i2+i2/3, julo:jhi) - Concentration contribution from N-S advection [mixing ratio]
- `j1p::Integer` - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i2_gl::Integer` - Global max longitude index
- `ju1_gl::Integer` - Global min latitude index
- `j2_gl::Integer` - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function xmist!(
	dcx::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	r24 = 1.0 / 24.0

	for ij = (j1p + 1):(j2p - 1), il = i1:i2
		tmp = ((8.0 * (qqv[il + 1, ij] - qqv[il - 1, ij])) + qqv[il - 2, ij] - qqv[il + 2, ij]) * r24

		pmax = max(qqv[il - 1, ij], qqv[il, ij], qqv[il + 1, ij]) - qqv[il, ij]

		pmin = qqv[il, ij] - min(qqv[il - 1, ij], qqv[il, ij], qqv[il + 1, ij])

		dcx[il, ij] = sign(tmp) * min(abs(tmp), pmax, pmin)
	end

	# Populate ghost zones of dcx (ccc, 11/20/08)
	for ij = julo:jhi
		for il = (-i2 / 3):0
			dcx[il, ij] = dcx[i2 + il, ij]
		end

		for il = (i2 + 1):(i2 + i2 / 3)
			dcx[il, ij] = dcx[il - i2, ij]
		end
	end
end

"""
The 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the E-W direction.

## Arguments
- `ij::Integer` - IN - Latitude index
- `ilmt::Integer` - IN - Controls various options in E-W advection
- `crx::Matrix{AbstractFloat}` - IN - (i1:i2, julo:jhi)
- `dcx::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi)
- `fx::Matrix{AbstractFloat}` - OUT - (i1:i2, julo:jhi)
- `qqv::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi)
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Remarks
This routine is called from w/in a OpenMP parallel loop fro

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function fxppm!(
	ij::Integer,
	ilmt::Integer,
	crx::Matrix{AbstractFloat},
	dcx::Matrix{AbstractFloat},
	fx::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer
)::Nothing
	# Zero arrays (bmy, 12/5/08)
	a6 = zeros(AbstractFloat, ilo:ihi)
	al = zeros(AbstractFloat, ilo:ihi)
	ar = zeros(AbstractFloat, ilo:ihi)
	a61 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	al1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	ar1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	dcxi1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)
	qqvi1 = zeros(AbstractFloat, (ihi-1) - (ilo+1) + 1)

	r13 = 1.0 / 3.0
	r23 = 2.0 / 3.0

	for il = (ilo + 1):ihi
		rval = 0.5 * (qqv[il - 1, ij] + qqv[il, ij]) + (dcx[il - 1, ij] - dcx[il, ij]) * r13
		al[il] = rval
		ar[il - 1] = rval
	end

	for il = (ilo + 1):(ihi - 1)
		a6[il] = 3.0 * (qqv[il, ij] + qqv[il, ij] - (al[il] + ar[il]))
	end
	
	if ilmt <= 2
		a61[:] .= 0.0
		al1[:] .= 0.0
		ar1[:] .= 0.0

		dcxi1[:] .= 0.0
		qqvi1[:] .= 0.0

		lenx = 0
		for il = (ilo + 1):(ihi - 1)
			lenx = lenx + 1

			a61[lenx] = a6[il]
			al1[lenx] = al[il]
			ar1[lenx] = ar[il]

			dcxi1[lenx] = dcx[il, ij]
			qqvi1[lenx] = qqv[il, ij]
		end
		
		# TODO:
		lmtppm!(lenx, ilmt, a61, al1, ar1, dcxi1, qqvi1)

		lenx = 0
		for il = (ilo + 1):(ihi - 1)
			lenx = lenx + 1

			a6[il] = a61[lenx]
			al[il] = al1[lenx]
			ar[il] = ar1[lenx]

			dcx[il, ij] = dcxi1[lenx]
			qqv[il, ij] = qqvi1[lenx]
		end

		# Populate ghost zones of qqv and dcx with new values (ccc, 11/20/08)
		for il = (-i2 / 3):0
			dcx[il, ij] = dcx[i2 + il, ij]
			qqv[il, ij] = qqv[i2 + il, ij]
		end

		for il = (i2 + 1):(i2 + i2 / 3)
			dcx[il, ij] = dcx[il - i2, ij]
			qqv[il, ij] = qqv[il - i2, ij]
		end
	end

	for il = (i1 + 1):i2
		if crx[il, ij] > 0.0
			ilm1 = il - 1
			fx[il, ij] = ar[ilm1] + 0.5 * crx[il, ij] * (al[ilm1] - ar[ilm1] + (a6[ilm1] * (1.0 - (r23 * crx[il, ij]))))
		else
			fx[il, ij] = al[il] - 0.5 * crx[il, ij] * (ar[il] - al[il] + (a6[il] * (1.0 + (r23 * crx[il, ij]))))
		end
	end

	# First box case (ccc, 11/20/08)
	if crx[i1, ij] > 0.0
		ilm1 = i2
		fx[i1, ij] = ar[ilm1] + 0.5 * crx[i1, ij] * (al[ilm1] - ar[ilm1] + (a6[ilm1] * (1.0 - (r23 * crx[i1, ij]))))
	else
		fx[i1, ij] = al[i1] - 0.5 * crx[i1, ij] * (ar[i1] - al[i1] + (a6[i1] * (1.0 + (r23 * crx[i1, ij]))))
	end
end

"""
Enforces the full monotonic, semi-monotonic, or the positive-definite constraint to the sub-grid parabolic distribution of the Piecewise Parabolic Method (PPM).

## Arguments
- `lenx::Integer` - IN - Vector length
- `lmt::Integer` - IN
	- If 0 => full monotonicity;
	- If 1 => semi-monotonic constraint (no undershoots);
	- If 2 => positive-definite constraint
- `a6::Array{AbstractFloat}` - (lenx) - INOUT - Curvature of the test parabola
- `al::Array{AbstractFloat}` - (lenx) - INOUT - Left edge value of the test parabola
- `ar::Array{AbstractFloat}` - (lenx) - INOUT - Right edge value of the test parabola
- `dc::Array{AbstractFloat}` - (lenx) - INOUT - 0.5 * mismatch
- `qa::Array{AbstractFloat}` - (lenx) - INOUT - Cell-averaged value

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function lmtppm!(
	lenx::Integer,
	lmt::Integer,
	a6::Array{AbstractFloat},
	al::Array{AbstractFloat},
	ar::Array{AbstractFloat},
	dc::Array{AbstractFloat},
	qa::Array{AbstractFloat}
)::Nothing
	r12 = 1.0 / 12.0
	
	if lmt == 0
		# Full constraint.
		
		for il = 1:lenx
			if dc[il] == 0.0
					a6[il] = 0.0
					al[il] = qa[il]
					ar[il] = qa[il]
			else
				da1 = ar[il] - al[il]
				da2 = da1 * da1
				a6da = a6[il] * da1

				if a6da < -da2
					a6[il] = 3.0 * (al[il] - qa[il])
					ar[il] = al[il] - a6[il]
				elseif a6da > da2
					a6[il] = 3.0 * (ar[il] - qa[il])
					al[il] = ar[il] - a6[il]
				end
			end
		end
	elseif lmt == 1
		# Semi-monotonic constraint.
		
		for il = 1:lenx
			if abs(ar[il] - al[il]) < -a6[il]
				if qa[il] < ar[il] && qa[il] < al[il]
					a6[il] = 0.0
					al[il] = qa[il]
					ar[il] = qa[il]
				elseif ar[il] > al[il]
					a6[il] = 3.0 * (al[il] - qa[il])
					ar[il] = al[il] - a6[il]
				else
					a6[il] = 3.0 * (ar[il] - qa[il])
					al[il] = ar[il] - a6[il]
				end
			end
		end
	elseif lmt == 2
		# positive-definite constraint
		
		for il = 1:lenx
			if abs(ar[il] - al[il]) < -a6[il]

				ftmp = ar[il] - al[il]

				fmin = qa[il] + 0.25 * (ftmp * ftmp) / a6[il] + a6[il] * r12

				if fmin < 0.0
					if qa[il] < ar[il] && qa[il] < al[il]
						a6[il] = 0.0
						al[il] = qa[il]
						ar[il] = qa[il]
					elseif ar[il] > al[il]
						a6[il] = 3.0 * (al[il] - qa[il])
						ar[il] = al[il] - a6[il]
					else
						a6[il] = 3.0 * (ar[il] - qa[il])
						al[il] = ar[il] - a6[il]
					end
				end
			end
		end
	end
end

"""
Does horizontal advection in the N-S direction.

## Arguments
- `jlmt::Integer` - IN - Controls various options in N-S advection
- `geofac_pc::AbstractFloat` - IN - special geometrical factor (geofac) for Polar cap
- `geofac::Matrix{AbstractFloat}` - IN - (ju1_gl:j2_gl) - geometrical factor for meridional advection; geofac uses correct spherical geometry, and replaces acosp as the  meridional geometrical factor in tpcore
- `cry::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Courant number in N-S direction
- `dq1::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Species density [hPa]
- `qqu::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Concentration contribution from E-W advection [mixing ratio]
- `qqv::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Concentration contribution from N-S advection [mixing ratio]
- `ymass::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Horizontal mass flux in N-S direction [hPa]
- `fy::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi+1) - N-S flux [mixing ratio]
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilong::Integer` - IN - ???
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `jord::Integer` - IN - N-S transport scheme (see module header for more info)

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function ytp!(
	jlmt::Integer,
	geofac_pc::AbstractFloat,
	geofac::Matrix{AbstractFloat},
	cry::Matrix{AbstractFloat},
	dq1::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	qqv::Matrix{AbstractFloat},
	ymass::Matrix{AbstractFloat},
	fy::Matrix{AbstractFloat},
	j1p::Integer,
	j2p::Integer,
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilong::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer,
	jord::Integer
)::Nothing
	dcy = zeros(AbstractFloat, ilo:ihi, julo:jhi)
	fy[:, :] .= 0.0

	rj1p = j1p
	
	if jord == 1
		for ij = j1p:(j2p + 1), il = i1:i2
			# c?
			jv = rj1p - cry[il, ij]
			qqv[il, ij] = qqu[il, jv]
		end
	else
			# TODO:
			ymist!(4, dcy, qqu, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
		if jord <= 0 || jord >= 3
			# TODO:
			fyppm!(jlmt, cry, dcy, qqu, qqv, j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
		else
			for ij = j1p:(j2p + 1), il = i1:i2
				# c?
				jv = rj1p - cry[il, ij]
				qqv[il, ij] = qqu[il, jv] + ((sign(cry[il, ij]) - cry[il, ij]) * dcy[il, jv])
			end
		end
	end

	for ij = j1p:(j2p + 1)
		qqv[i1:i2, ij] = qqv[i1:i2, ij] * ymass[i1:i2, ij]
	end

	# .sds.. save N-S species flux as diagnostic
	for ij = i1:i2
		fy[ij, j1p:(j2p + 1)] = qqv[ij, j1p:(j2p + 1)] * geofac[j1p:(j2p + 1)]
	end

	# ... meridional flux update
	for ij = j1p:j2p
		dq1[i1:i2, ij] = dq1[i1:i2, ij] + (qqv[i1:i2, ij] - qqv[i1:i2, ij + 1]) * geofac[ij]
	end
	
	# TODO:
	do_ytp_pole_sum!(geofac_pc, dq1, qqv, fy, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
end

"""
Computes the linear tracer slope in the N-S direction.  It uses the Lin et. al. 1994 algorithm.

## Arguments
- `id::Integer` - IN - The "order" of the accuracy in the computed linear "slope" (or mismatch, Lin et al. 1994); it is either 2 or 4.
- `dcy::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Slope of concentration distribution in N-S direction [mixing ratio]
- `qqu::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Concentration contribution from E-W advection (mixing ratio)
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude index at the edge of the South polar cap; j1p=ju1_gl+1 for a polar cap of 1 latitude band
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function ymist!(
	id::Integer,
	dcy::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	# I suppose the values for these indexes are 0. It should work as the pole values are re-calculated in the pole functions. (ccc)
	qtmp = zeros(AbstractFloat, ilo:ihi, (julo - 2):(jhi + 2))

	r24  = 1.0 / 24.0

	# Populate qtmp
	qtmp = 0.0
	for ij = ju1:j2
			qtmp[:, ij] .= qqu[:, ij]
	end
	
	if id == 2
		for ij = (ju1 - 1):(j2 - 1), il = i1:i2
			tmp  = 0.25 * (qtmp[il, ij + 2] - qtmp[il, ij])

			pmax = max(qtmp[il, ij], qtmp[il, ij + 1], qtmp[il, ij + 2]) - qtmp[il, ij + 1]

			pmin = qtmp[il, ij + 1] - min(qtmp[il, ij], qtmp[il, ij + 1], qtmp[il, ij + 2])

			dcy[il, ij + 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	else
		# TODO:
		do_ymist_pole1_i2d2!(dcy, qtmp, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

		for ij = (ju1 - 2):(j2 - 2), il = i1:i2
			tmp = ((8.0 * (qtmp[il, ij + 3] - qtmp[il, ij + 1])) + qtmp[il, ij] - qtmp[il, ij + 4]) * r24

			pmax = max(qtmp[il, ij + 1], qtmp[il, ij + 2], qtmp[il, ij + 3]) - qtmp[il, ij + 2]

			pmin = qtmp[il, ij + 2] - min(qtmp[il, ij + 1], qtmp[il, ij + 2], qtmp[il, ij + 3])

			dcy[il, ij + 2] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	end
	
	# TODO:
	do_ymist_pole2_i2d2!(dcy, qtmp, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
end

"""
Sets "dcy" at the Poles.

## Arguments
- `dcy::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Slope of concentration distribution in N-S direction [mixing ratio]
- `qqu::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo-2:jhi+2) - Concentration contribution from E-W advection [mixing ratio]
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_ymist_pole1_i2d2!(
	dcy::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	i2d2 = i2_gl / 2

	r24  = 1.0 / 24.0
	
	if ju1 == ju1_gl
		for il = i1:i2d2
			tmp = ((8.0 * (qqu[il, ju1 + 2] - qqu[il, ju1])) + qqu[il + i2d2, ju1 + 1] - qqu[il, ju1 + 3]) * r24

			pmax = max(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2]) - qqu[il, ju1 + 1]

			pmin = qqu(il,ju1+1) - min(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2])

			dcy[il, ju1 + 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end

		for il = (i1 + i2d2):i2
			tmp = ((8.0 * (qqu[il, ju1 + 2] - qqu[il, ju1])) + qqu[il - i2d2, ju1 + 1] - qqu[il, ju1 + 3]) * r24

			pmax = max(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2]) - qqu[il, ju1 + 1]

			pmin = qqu[il, ju1 + 1] - min(qqu[il, ju1], qqu[il, ju1 + 1], qqu[il, ju1 + 2])

			dcy[il, ju1 + 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	end

	if j2 == j2_gl
		for il = i1:i2d2
			tmp = ((8.0 * (qqu[il, j2] - qqu[il, j2 - 2])) + qqu[il, j2 - 3] - qqu[il+ i2d2, j2 - 1]) * r24

			pmax = max(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2]) - qqu[il, j2 - 1]

			pmin = qqu[il, j2 - 1] - min(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2])

			dcy[il, j2 - 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end

		for il = (i1 + i2d2):i2
			tmp = ((8.0 * (qqu[il, j2] - qqu[il, j2 - 2])) + qqu[il, j2 - 3] - qqu[il - i2d2, j2 - 1]) * r24

			pmax = max(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2]) - qqu[il, j2 - 1]

			pmin = qqu[il, j2 - 1] - min(qqu[il, j2 - 2], qqu[il, j2 - 1], qqu[il, j2])

			dcy[il, j2 - 1] = sign(tmp) * min(abs(tmp), pmin, pmax)
		end
	end
end

"""
Sets "dcy" at the Poles.

## Arguments
- `dcy::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Slope of concentration distribution in N-S direction [mixing ratio]
- `qqu::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo-2:jhi+2) - Concentration contribution from E-W advection [mixing ratio]
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude index at the edge of the South polar cap;  j1p=ju1_gl+1 for a polar cap of 1 latitude band
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos-chem for complete history.
"""
function do_ymist_pole2_i2d2!(
	dcy::Matrix{AbstractFloat},
	qqu::Matrix{AbstractFloat},
	i1_gl::Integer,
	i2_gl::Integer,
	ju1_gl::Integer,
	j2_gl::Integer,
	j1p::Integer,
	ilo::Integer,
	ihi::Integer,
	julo::Integer,
	jhi::Integer,
	i1::Integer,
	i2::Integer,
	ju1::Integer,
	j2::Integer
)::Nothing
	i2d2 = i2_gl / 2
	
	if ju1 == ju1_gl
		if j1p != ju1_gl + 1
			dcy[i1:i2, ju1] .= 0.0
		else
			# Determine slope in South Polar cap for scalars.
			for il = i1:i2d2
				tmp = 0.25 * (qqu[il, ju1 + 1] - qqu[il + i2d2, ju1 + 1])

				pmax = max(qqu[il, ju1 + 1], qqu[il, ju1], qqu[il + i2d2, ju1 + 1]) - qqu[il, ju1]

				pmin = qqu[il, ju1] - min(qqu[il, ju1 + 1], qqu[il, ju1], qqu[il + i2d2, ju1 + 1])
							
				dcy[il, ju1] = sign(tmp) * min(abs(tmp), pmax, pmin)
			end

			for il = (i1 + i2d2):i2
				dcy[il, ju1] = -dcy[il - i2d2, ju1]
			end
		end
	end
		
	if j2 == j2_gl
		if j1p != ju1_gl + 1
			dcy[i1:i2, j2] .= 0.0
		else
			# Determine slope in North Polar cap for scalars.

			for il = i1:i2d2
				tmp  = 0.25 * (qqu[il + i2d2, j2 - 1] - qqu[il, j2 - 1])

				pmax = max(qqu[il + i2d2, j2 - 1], qqu[il, j2], qqu[il, j2 - 1]) - qqu[il, j2]

				pmin = qqu[il, j2] - min(qqu[il + i2d2, j2 - 1], qqu[il, j2], qqu[il, j2 - 1])

				dcy[il, j2] = sign(tmp) * min(abs(tmp), pmax, pmin)
			end

			for il = (i1 + i2d2):i2
				dcy[il, j2] = -dcy[il - i2d2, j2]
			end
		end
	end
end

"""
The 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the N - S direction.

## Arguments
- `jlmt::Integer` - IN - Controls various options in N-S advection
- `cry::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Courant number in N-S direction
- `dcy::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Slope of concentration distribution in N-S direction [mixing ratio]
- `qqu::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Concentration contribution from E-W advection [mixing ratio]
- `qqv::Matrix{AbstractFloat}` - OUT - (ilo:ihi, julo:jhi) - Concentration contribution from N-S advection [mixing ratio]
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+1; j2p=j2_gl-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; j1p=ju1_gl+2; j2p=j2_gl-2 for a polar cap of 2 latitude band
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilong::Integer` - IN - ilong ??
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history.
"""
function fyppm!(
  jlmt::Integer,
  cry::Matrix{AbstractFloat},
  dcy::Matrix{AbstractFloat},
  qqu::Matrix{AbstractFloat},
  qqv::Matrix{AbstractFloat},
  j1p::Integer,
  j2p::Integer,
  i1_gl::Integer,
  i2_gl::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  ilong::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer
)::Nothing
  a61 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  al1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  ar1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  dcy1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  qqu1 = zeros(AbstractFloat, ilong * ((jhi - 1) - (julo + 1) + 1))
  a6 = zeros(AbstractFloat, ilo:ihi, julo:jhi)
  al = zeros(AbstractFloat, ilo:ihi, julo:jhi)
  ar = zeros(AbstractFloat, ilo:ihi, julo:jhi)

  # NOTE: The code was writtein with I1:I2 as the first dimension of AL, AR, A6, AL1, A61, AR1. However, the limits should really should be ILO:IHI. In practice, however, for a global grid (and OpenMP parallelization) ILO=I1 and IHI=I2. Nevertheless, we will change the limits to ILO:IHI to be consistent and to avoid future problems. (bmy, 12/5/08)

  r13 = 1.0 / 3.0
  r23 = 2.0 / 3.0

  for ij = (julo + 1):jhi, il = ilo:ihi
    al[il, ij] = 0.5 * (qqu[il, ij - 1] + qqu[il, ij]) + (dcy[il, ij - 1] - dcy[il, ij]) * r13
    ar[il, ij - 1] = al[il, ij]
  end

  # TODO:
  do_fyppm_pole_i2d2!(al, ar, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

  for ij = (julo + 1):(jhi - 1), il = ilo:ihi
    a6[il, ij] = 3.0 * (qqu[il, ij] + qqu[il, ij] - (al[il, ij] + ar[il, ij]))
  end

  if jlmt <= 2
    lenx = 0

    for ij = (julo + 1):(jhi - 1), il = ilo:ihi
      lenx = lenx + 1

      a61[lenx] = a6[il, ij]
      al1[lenx] = al[il, ij]
      ar1[lenx] = ar[il, ij]
      dcy1[lenx] = dcy[il, ij]
      qqu1[lenx] = qqu[il, ij]
    end

    # TODO:
    lmtppm!(lenx, jlmt, a61, al1, ar1, dcy1, qqu1)

    lenx = 0

    for ij = (julo + 1):(jhi - 1), il = ilo:ihi
      lenx = lenx + 1

      a6[il, ij] = a61[lenx]
      al[il, ij] = al1[lenx]
      ar[il, ij] = ar1[lenx]
    end
  end

  for ij = j1p:j2p + 1
    ijm1 = ij - 1

    for il = ilo:ihi
      if cry[il, ij] > 0.0
        qqv[il, ij] = ar[il, ijm1] + 0.5 * cry[il, ij] * (al[il, ijm1] - ar[il, ijm1] + (a6[il, ijm1] * (1.0 - (r23 * cry[il, ij]))))
      else
        qqv[il, ij] = al[il, ij] - 0.5 * cry[il, ij] * (ar[il, ij] - al[il, ij] + (a6[il, ij] * (1.0 + (r23 * cry[il, ij]))))
      end
    end
  end
end

"""
Sets "al" & "ar" at the Poles.

## Arguments
- `al::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Left edge value of the test parabola
- `ar::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Right edge value of the test parabola
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history
"""
function do_fyppm_pole_i2d2!(
  al::Matrix{AbstractFloat},
  ar::Matrix{AbstractFloat},
  i1_gl::Integer,
  i2_gl::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer
)::Nothing
  i2d2 = i2_gl / 2

  for il = i1:i2d2
    al[il, ju1] = al[il + i2d2, ju1 + 1]
    al[il + i2d2, ju1] = al[il, ju1 + 1]
    ar[il, j2] = ar[il + i2d2, j2 - 1]
    ar[il + i2d2, j2] = ar[il, j2 - 1]
  end
end

"""
Sets "dq1" at the Poles.

## Arguments
- `geofac_pc::AbstractFloat` - IN - Special geometrical factor (geofac) for Polar cap
- `dq1::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi) - Species density [hPa]
- `qqv::Matrix{AbstractFloat}` - IN - (ilo:ihi, julo:jhi) - Concentration contribution from N-S advection [mixing ratio]
- `fy::Matrix{AbstractFloat}` - INOUT - (ilo:ihi, julo:jhi+1) - N-S mass flux [mixing ratio]
- `i1_gl::Integer` - IN - Global min longitude index
- `i2_gl::Integer` - IN - Global max longitude index
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `j1p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
- `j2p::Integer` - IN - Global latitude indices at the edges of the S/N polar caps; J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude band
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history.
"""
function do_ytp_pole_sum!(
  geofac_pc::AbstractFloat,
  dq1::Matrix{AbstractFloat},
  qqv::Matrix{AbstractFloat},
  fy::Matrix{AbstractFloat},
  i1_gl::Integer,
  i2_gl::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  j1p::Integer,
  j2p::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer
)::Nothing
  ri2 = i2_gl

  dqik = zeros(AbstractFloat, 2)

  # ... Integrate N - S flux around polar cap lat circle for each level

  sumsp = 0.0
  sumnp = 0.0
  for il = i1:i2
    sumsp = sumsp + qqv[il, j1p]
    sumnp = sumnp + qqv[il, j2p + 1]
  end

  # ... wrap in E - W direction
  if i1 == i1_gl
    dqik[1] = dq1[i1, ju1]
    dqik[2] = dq1[i1, j2]
  end

  # ... normalize and set inside polar cap

  dq_sp = dqik[1] - (sumsp / ri2 * geofac_pc)
  dq_np = dqik[2] + (sumnp / ri2 * geofac_pc)

  for il = i1:i2
    dq1[il, ju1] = dq_sp
    dq1[il, j2] = dq_np
    # ... save polar flux
    fy[il, ju1] =  - (sumsp / ri2 * geofac_pc)
    fy[il, j2 + 1] = (sumnp / ri2 * geofac_pc)
  end

  if j1p != ju1_gl + 1
    for il = i1:i2
      dq1[il, ju1 + 1] = dq_sp
      dq1[il, j2 - 1] = dq_np
      # ... save polar flux
      fy[il, ju1 + 1] =  - (sumsp / ri2 * geofac_pc)
      fy[il, j2] = sumnp / ri2 * geofac_pc
    end
  end
end

"""
The 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the vertical direction.

Modified by S. - J. Lin, 12/14/98, to allow the use of the KORD=7 (klmt=4) option. KORD=7 enforces the 2nd monotonicity constraint of Huynh (1996). Note that in Huynh"s original scheme, two constraints are necessary for the preservation of monotonicity. To use Huynh"s algorithm, it was modified as follows. The original PPM is still used to obtain the first guesses for the cell edges, and as such Huynh"s 1st constraint is no longer needed. Huynh"s median function is also replaced by a simpler yet functionally equivalent in - line algorithm.

## Arguments
- `klmt::Integer` - IN - Controls various options in vertical advection
- `delp1::Array{AbstractFloat, 3}` - IN - (ilo:ihi, julo:jhi, k1:k2) - Pressure thickness, the pseudo-density in a hydrostatic system at t1 [hPa]
- `wz::Array{AbstractFloat, 3}` - IN - (i1:i2, ju1:j2, k1:k2) - Large scale mass flux (per time step tdt) in the vertical direction as diagnosed from the hydrostatic relationship [hPa]
- `dq1::Array{AbstractFloat, 3}` - INOUT - (ilo:ihi, julo:jhi, k1:k2) - Species density [hPa]
- `qq1::Array{AbstractFloat, 3}` - IN - (:, :, :) - Species concentration [mixing ratio]
- `fz::Array{AbstractFloat, 3}` - (ilo:ihi, julo:jhi,  k1:k2) - OUT - Vertical flux [mixing ratio]
- `j1p::Integer` - IN - Global latitude index at the edges of the South polar cap. j1p=ju1_gl+1 for a polar cap of 1 latitude band. J1P=JU1_GL+2 for a polar cap of 2 latitude bands.
- `ju1_gl::Integer` - IN - Global min latitude index
- `j2_gl::Integer` - IN - Global max latitude index
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index
- `ilong::Integer` - IN - Dimensions in longitude ???
- `ivert::Integer` - IN - Dimensions in altitude ???
- `i1::Integer` - IN - Local min longitude index
- `i2::Integer` - IN - Local max longitude index
- `ju1::Integer` - IN - Local min latitude index
- `j2::Integer` - IN - Local max latitude index
- `k1::Integer` - IN - Local min altitude index
- `k2::Integer` - IN - Local max altitude index

## Author
Original code from Shian - Jiann Lin, DAO.
John Tannahill, LLNL (jrt@llnl.gov).

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history
"""
function fzppm!(
  klmt::Integer,
  delp1::Array{AbstractFloat,3},
  wz::Array{AbstractFloat,3},
  dq1::Array{AbstractFloat,3},
  qq1::Array{AbstractFloat,3},
  fz::Array{AbstractFloat,3},
  j1p::Integer,
  ju1_gl::Integer,
  j2_gl::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer,
  ilong::Integer,
  ivert::Integer,
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer,
  k1::Integer,
  k2::Integer
)::Nothing
  a61 = zeros(AbstractFloat, ilong * (ivert - 4))
  al1 = zeros(AbstractFloat, ilong * (ivert - 4))
  ar1 = zeros(AbstractFloat, ilong * (ivert - 4))
  dca1 = zeros(AbstractFloat, ilong * (ivert - 4))
  qq1a1 = zeros(AbstractFloat, ilong * (ivert - 4))

  a6 = zeros(AbstractFloat, i1:i2, k1:k2)
  al = zeros(AbstractFloat, i1:i2, k1:k2)
  ar = zeros(AbstractFloat, i1:i2, k1:k2)
  dca = zeros(AbstractFloat, i1:i2, k1:k2)
  dlp1a = zeros(AbstractFloat, i1:i2, k1:k2)
  qq1a = zeros(AbstractFloat, i1:i2, k1:k2)
  wza = zeros(AbstractFloat, i1:i2, k1:k2)

  dc = zeros(AbstractFloat, i1:i2, ju1:j2, k1:k2)

  # Work array
  dp = zeros(AbstractFloat, i1:i2, ju1:j2, k1:k2)

  # .sds... diagnostic vertical flux for species - set top to 0.0
  fz[:, :, :] .= zeros(AbstractFloat, ilo:ihi, julo:jhi, k1:k2)

  k1p1 = k1 + 1
  k1p2 = k1 + 2

  k2m1 = k2 - 1
  k2m2 = k2 - 2

  r13 = 1.0 / 3.0
  r23 = 2.0 / 3.0

  # Compute dc for PPM.

  for ik = k1:k2m1
    dpi[:, :, ik] .= qq1[i1:i2, ju1:j2, ik + 1] - qq1[i1:i2, ju1:j2, ik]
  end

  for ik = k1p1:k2m1, ij = ju1:j2, il = i1:i2
    c0 = delp1[il, ij, ik] / (delp1[il, ij, ik - 1] + delp1[il, ij, ik] + delp1[il, ij, ik + 1])
    c1 = (delp1[il, ij, ik - 1] + (0.5 * delp1[il, ij, ik])) / (delp1[il, ij, ik + 1] + delp1[il, ij, ik])
    c2 = (delp1[il, ij, ik + 1] + (0.5 * delp1[il, ij, ik])) / (delp1[il, ij, ik - 1] + delp1[il, ij, ik])

    tmp = c0 * ((c1 * dpi[il, ij, ik]) + (c2 * dpi[il, ij, ik - 1]))

    qmin = qq1[il, ij, ik] - min(qq1[il, ij, ik - 1], qq1[il, ij, ik], qq1[il, ij, ik + 1])
    qmax = max(qq1[il, ij, ik - 1], qq1[il, ij, ik], qq1[il, ij, ik + 1]) - qq1[il, ij, ik]

    dc[il, ij, ik] = min(abs(tmp), qmax, qmin) * sign(tmp)
  end

  # c?
  # Loop over latitudes (to save memory).

  for ij = ju1:j2 
    if (ij == ju1_gl + 1 || ij == j2_gl - 1) && (j1p != ju1_gl + 1)
      @goto ijloop
    end

    for ik = k1:k2, il = i1:i2
      dca[il, ik] = dc[il, ij, ik] # the monotone slope
      wza[il, ik] = wz[il, ij, ik]
      dlp1a[il, ik] = delp1[il, ij, ik]
      qq1a[il, ik] = qq1[il, ij, ik]
    end

    # Compute first guesses at cell interfaces. First guesses are required to be continuous. Three - cell parabolic subgrid distribution at model top; two - cell parabolic with zero gradient subgrid distribution at the surface.

    # First guess top edge value.

    for il = i1:i2
      # Three - cell PPM; compute a, b, & c of q = aP^2 + bP + c using cell averages and dlp1a.

      fac1 = dpi[il, ij, k1p1] - dpi[il, ij, k1] * (dlp1a[il, k1p1] + dlp1a[il, k1p2]) / (dlp1a[il, k1] + dlp1a[il, k1p1])

      fac2 = (dlp1a[il, k1p1] + dlp1a[il, k1p2]) * (dlp1a[il, k1] + dlp1a[il, k1p1] + dlp1a[il, k1p2])

      aa = 3.0 * fac1 / fac2

      bb = 2.0 * dpi[il, ij, k1] / (dlp1a[il, k1] + dlp1a[il, k1p1]) - r23 * aa * (2.0 * dlp1a[il, k1] + dlp1a[il, k1p1])

      al[il, k1] = qq1a[il, k1] - dlp1a[il, k1] * (r13 * aa * dlp1a[il, k1] + 0.5 * bb)

      al[il, k1p1] = dlp1a[il, k1] * (aa * dlp1a[il, k1] + bb) + al[il, k1]

      # Check if change sign.

      if qq1a[il, k1] * al[il, k1] <= 0.0
        al[il, k1] = 0.0
        dca[il, k1] = 0.0
      else
        dca[il, k1] = qq1a[il, k1] - al[il, k1]
      end
    end

    # Bottom.
    for il = i1:i2
      # 2 - cell PPM with zero gradient right at the surface.

      fac1 = dpi[il, ij, k2m1] * (dlp1a[il, k2] * dlp1a[il, k2]) / ((dlp1a[il, k2] + dlp1a[il, k2m1]) * (2.0 * dlp1a[il, k2] + dlp1a[il, k2m1]))

      ar[il, k2] = qq1a[il, k2] + fac1
      al[il, k2] = qq1a[il, k2] - (fac1 + fac1)

      if qq1a[il, k2] * ar[il, k2] <= 0.0
        ar[il, k2] = 0.0
      end

      dca[il, k2] = ar[il, k2] - qq1a[il, k2]
    end

    # 4th order interpolation in the interior.

    for ik = k1p2:k2m1, il = i1:i2
      c1 = (dpi[il, ij, ik - 1] * dlp1a[il, ik - 1]) / (dlp1a[il, ik - 1] + dlp1a[il, ik])
      c2 = 2.0 / (dlp1a[il, ik - 2] + dlp1a[il, ik - 1] + dlp1a[il, ik] + dlp1a[il, ik + 1])

      a1 = (dlp1a[il, ik - 2] + dlp1a[il, ik - 1]) / (2.0 * dlp1a[il, ik - 1] + dlp1a[il, ik])
      a2 = (dlp1a[il, ik] + dlp1a[il, ik + 1]) / (2.0 * dlp1a[il, ik] + dlp1a[il, ik - 1])

      al[il, ik] = qq1a[il, ik - 1] + c1 + c2 * (dlp1a[il, ik] * (c1 * (a1 - a2) + a2 * dca[il, ik - 1]) - dlp1a[il, ik - 1] * a1 * dca[il, ik])
    end

    for ik = k1:k2m1, il = i1:i2
      ar[il, ik] = al[il, ik + 1]
    end

    # Top & Bottom 2 layers always monotonic.

    lenx = i2 - i1 + 1

    for ik = k1:k1p1
      for il = i1:i2
        a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (al[il, ik] + ar[il, ik]))
      end

      # TODO:
      lmtppm!(lenx, 0, a6[i1, ik], al[i1, ik], ar[i1, ik], dca[i1, ik], qq1a[i1, ik])
    end

    for ik = k2m1:k2
      for il = i1:i2
        a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (al[il, ik] + ar[il, ik]))
      end

      lmtppm!(lenx, 0, a6[i1, ik], al[i1, ik], ar[i1, ik], dca[i1, ik], qq1a[i1, ik])
    end

    # Interior depending on klmt.

    if klmt == 4
      # KORD=7, Huynh"s 2nd constraint.

      for ik = k1p1:k2m1, il = i1:i2
        dca[il, ik] = dpi[il, ij, ik] - dpi[il, ij, ik - 1]
      end

      for ik = k1p2:k2m2, il = i1:i2
        # Right edges.

        qmp = qq1a[il, ik] + (2.0 * dpi[il, ij, ik - 1])
        lac = qq1a[il, ik] + (1.5 * dca[il, ik - 1]) + (0.5 * dpi[il, ij, ik - 1])
        qmin = min(qq1a[il, ik], qmp, lac)
        qmax = max(qq1a[il, ik], qmp, lac)

        ar[il, ik] = min(max(ar[il, ik], qmin), qmax)

        # Left edges.

        qmp = qq1a[il, ik] - (2.0 * dpi[il, ij, ik])
        lac = qq1a[il, ik] + (1.5 * dca[il, ik + 1]) - (0.5 * dpi[il, ij, ik])
        qmin = min(qq1a[il, ik], qmp, lac)
        qmax = max(qq1a[il, ik], qmp, lac)

        al[il, ik] = min(max(al[il, ik], qmin), qmax)

        # Recompute a6.

        a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (ar[il, ik] + al[il, ik]))
      end
    elseif klmt <= 2
      lenx = 0

      for ik = k1p2:k2m2, il = i1:i2
        lenx = lenx + 1

        al1[lenx] = al[il, ik]
        ar1[lenx] = ar[il, ik]
        dca1[lenx] = dca[il, ik]
        qq1a1[lenx] = qq1a[il, ik]

        a61[lenx] = 3.0 * (qq1a1[lenx] + qq1a1[lenx] - (al1[lenx] + ar1[lenx]))
      end

      # TODO:
      lmtppm!(lenx, klmt, a61, al1, ar1, dca1, qq1a1)

      lenx = 0
      for ik = k1p2:k2m2, il = i1:i2
        lenx = lenx + 1

        a6[il, ik] = a61[lenx]
        al[il, ik] = al1[lenx]
        ar[il, ik] = ar1[lenx]
        dca[il, ik] = dca1[lenx]
        qq1a[il, ik] = qq1a1[lenx]
      end
    end

    for ik = k1:k2m1, il = i1:i2
      if wza[il, ik] > 0.0
        cm = wza[il, ik] / dlp1a[il, ik]

        dca[il, ik + 1] = ar[il, ik] + 0.5 * cm * (al[il, ik] - ar[il, ik] + a6[il, ik] * (1.0 - r23 * cm))
      else
        cp = wza[il, ik] / dlp1a[il, ik + 1]

        dca[il, ik + 1] = al[il, ik + 1] + 0.5 * cp * (al[il, ik + 1] - ar[il, ik + 1] - a6[il, ik + 1] * [1.0 + r23 * cp])

      end
    end

    for ik = k1:k2m1, il = i1:i2
      dca[il, ik + 1] = wza[il, ik] * dca[il, ik + 1]
      # .sds.. save vertical flux for species as diagnostic
      fz[il, ij, ik + 1] = dca[il, ik + 1]
    end

    for il = i1:i2
      dq1[il, ij, k1] = dq1[il, ij, k1] - dca[il, k1p1]
      dq1[il, ij, k2] = dq1[il, ij, k2] + dca[il, k2]
    end

    for ik = k1p1:k2m1, il = i1:i2
      dq1[il, ij, ik] = dq1[il, ij, ik] + dca[il, ik] - dca[il, ik + 1]
    end
		
		@label ijloop
  end
end

"""
Averages pressure at the Poles when the Polar cap is enlarged. It makes the last two latitudes equal.

## Arguments
- `area_1d::Array{AbstractFloat}` - (ju1:j2) - IN - Surface area of grid box
- `press::Matrix{AbstractFloat}` - (ilo:ihi, julo:jhi) - INOUT - Surface pressure [hPa]
- `i1::Integer` - IN - Local min longitude
- `i2::Integer` - IN - Local max longitude
- `ju1::Integer` - IN - Local min latitude
- `j2::Integer` - IN - Local max latitude
- `ilo::Integer` - IN - Local min longitude index
- `ihi::Integer` - IN - Local max longitude index
- `julo::Integer` - IN - Local min latitude index
- `jhi::Integer` - IN - Local max latitude index

## Author
Philip Cameron - Smith and John Tannahill, GMI project @ LLNL (2003).
Implemented into GEOS - Chem by Claire Carouge (ccarouge@seas.harvard.edu).

## Remarks
function from pjc_pfix. Call this one once everything is working fine.

## Revision History
05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S - J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. See https://github.com/geoschem/geos - chem for complete history.
"""
function average_press_poles!(
  area_1d::Array{AbstractFloat},
  press::Matrix{AbstractFloat},
  i1::Integer,
  i2::Integer,
  ju1::Integer,
  j2::Integer,
  ilo::Integer,
  ihi::Integer,
  julo::Integer,
  jhi::Integer
)::Nothing
  # Compute the sum of surface area

  rel_area = zeros(AbstractFloat, ju1:j2)
  # TODO: looks like a map
  sum_area = sum[area_1d] * i2
  # calculate rel_area for each lat. (ccc, 11/20/08)
  for j = ju1:j2
    rel_area[j] = area_1d[j] / sum_area
  end

  # South Pole

  # Surface area of the S. Polar cap
  sum_area = sum(rel_area[ju1:(ju1 + 1)]) * i2

  # TODO: This looks like a reduce
  # Zero
  meanp = 0.0
  # Sum pressure * surface area over the S. Polar cap
  for j = ju1:(ju1 + 1), i = i1:i2
    meanp = meanp + (rel_area[j] * press[i, j])
  end

  # Normalize pressure in all grid boxes w/in the S. Polar cap
  press[:, ju1:(ju1 + 1)] .= meanp / sum_area

  # North Pole

  # Surface area of the N. Polar cap
  sum_area = sum(rel_area[(j2 - 1):j2]) * i2

  # TODO: This looks like a reduce
  # Zero
  meanp = 0.0
  # ! Sum pressure * surface area over the N. Polar cap
  for j = (j2 - 1):j2, i = i1:i2
    meanp = meanp + (rel_area[j] * press[i, j])
  end

  # ! Normalize pressure in all grid boxes w/in the N. Polar cap
  press[:, (j2 - 1):j2] .= meanp / sum_area
end

end # module tpcore_fvdas_mod

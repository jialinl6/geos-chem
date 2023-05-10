# $Id: tpcore_window_mod.f90, v 1.2 2009/10/15 17:46:23 bmy Exp $ 
module TPCORE_WINdoW_MOD
# 
#  
# Module TPCORE_WINdoW_MOD contains routines for the GEOS - 4/fvDAS 
# transport scheme.  Original code from S - J Lin and Kevin Yeh. 
# (bdf, bmy, 5/7/03, 10/29/04) 
# 
# The Harvard Atmospheric Chemistry Modeling Group has modified the original 
# code in order to implement the Philip - Cameron Smith pressure fixer for mass 
# conservation, and also to save out mass fluxes.  These changes are denoted 
#  in the code by comment tag lines # .  Also, all modifications to the 
# original code are written in ALL CAPITALS. 
# 
# Module Routines: 
# ============================================================================ 
# (1 ) INIT_TPCORE    : Initialization routine for TPCORE 
# (2 ) EXIT_TPCORE    : Cleanup and exit routine for TPCORE 
# (3 ) TPCORE_FVDAS   : Driver routine for GEOS - 4/TPCORE transport scheme 
# (4 ) AIR_MASS_FLUX  : TPCORE internal routine 
# (5 ) TP2G           : TPCORE internal routine 
# (6 ) TP2D           : TPCORE internal routine 
# (7 ) XTP            : TPCORE internal routine 
# (8 ) XMIST          : TPCORE internal routine 
# (9 ) FXPPM          : TPCORE internal routine 
# (10) LMPPM          : TPCORE internal routine 
# (11) HUYNH          : TPCORE internal routine 
# (12) YTP            : TPCORE internal routine 
# (13) YMIST          : TPCORE internal routine 
# (14) FYPPM          : TPCORE internal routine 
# (15) XPAVG          : TPCORE internal routine 
# (16) QMAP           : TPCORE internal routine 
# (17) MAP1_PPM       : TPCORE internal routine 
# (18) PPM2M          : TPCORE internal routine 
# (19) STEEPZ         : TPCORE internal routine 
# (20) KMPPM          : TPCORE internal routine 
# (21) FCT_X          : TPCORE internal routine 
# (22) FILLZ          : TPCORE internal routine 
# (23) PFIX           : TPCORE internal routine 
# (24) GMEAN          : TPCORE internal routine 
# (25) ADJ_FX         : TPCORE internal routine 
# 
# GEOS - CHEM modules referenced by "tpcore_fvdas_mod.f90" 
# ============================================================================ 
# none 
# 
# NOTES: 
# (1 ) Renamed this module from "transport_fv.F90" to "tpcore_fvdas_mod.f90" 
# to be more consistent with GEOS - CHEM naming convention. 
# (2 ) Renamed routine TPCORE to TPCORE_FVDAS to avoid conflict with 
# existing routine TPCORE from S - J Lin's version 7.1.0m. 
# (3 ) Added code for PJC pressure fixer.  Also now declare everything 
# PRIVATE except for INIT_TPCORE, TPCORE_FVDAS, and EXIT_TPCORE. 
# (bdf, bmy, 5/7/03) 
# (4 ) Added modifications to save mass fluxes in ND24, ND25, ND26 
# diagnostics.  Also now make places in the code which have been 
# modified by Harvard more clear to discern. (bdf, bmy, 9/28/04) 
# (5 ) Bug fix: Need to multiply ND25 N/S transport fluxes by the array 
# RGW_25 which accounts for the latitude factor (bdf, bmy, 10/29/04) 
# (6 ) Bug fix: In INIT_GEOS5_WINdoW, need to dimension COSE with JM + 1 
# instead of JM.  (Xiaoguang Gu, bmy, 1/20/09) 
# 09 Sep 2010 - C. Carouge  - Modify declarations of MASSFLEW, MASSFLNS and 
# MASSFLUP to save memory space. 
# 04 Nov 2015 - M. Sulprizio- Rename from tpcore_geos5_window_mod.F90 to 
# tpcore_window_mod.F90 for use with all nested 
# grids 
# 11 Jan 2016 - E. Lundgren - Add diagnostics for output to netcdf. 
# Block off both bpch and netcdf diagnostic 
# code within pre - processor statements. 
# 19 Jan 2016 - E. Lundgren - combine bpch and netcdf diagnostic code 
#  
# 
# The original module documentation header is listed here: 
# 
# TransPort module for NASA Goddard Chemistry Transport Model 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Last modified: February 27, 2002 
# 
# Purpose: perform the transport of  3 - D mixing ratio fields using 
# externally specified winds and surface pressure on the 
# hybrid Eta - coordinate. 
#          One # call to tpcore updates the 3 - D mixing ratio 
# fields for one time step (DT). 
# 
# Schemes: Multi - dimensional Flux Form Semi - Lagrangian (FFSL) schemes 
# (Lin and Rood 1996, MWR) with many unpublished modifications 
# 
# Messaging passing library based on "Pilgrim" developed by W. Sawyer 
# 
# Suggested compiler options: 
# SGI Origin: f90 -c -r8 -64 -O3 -mips4 -mp 
# loader: f90 -64 -mp 
# Linux Lahey/Fujitsu lf95 -c -CcdRR8 --tpp 
# 
# Send comments/suggestions to the algorithm developers: 
# 
# S.-J. Lin 
# Code 910.3, NASA/GSFC, Greenbelt, MD 20771 
# E - mail: slin@dao.gsfc.nasa.gov 
# 
# Kevin Yeh 
# Code 910.3, NASA/GSFC, Greenbelt, MD 20771 
# E - mail: kyeh@dao.gsfc.nasa.gov 
# 
# The algorithm is primarily based on the following papers: 
# 
# 1.0 Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi- 
# Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046 - 2070.0 
# 
# 2.0 Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of 
# the van Leer - type transport schemes and its applications to the moist- 
# ure transport in a General Circulation Model. Mon. Wea. Rev., 122, 
# 1575 - 1593.0 
#  
# 
#   # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#   # 
# Add MODULE PRIVATE declarations.  For safety's sake, declare all 
# routines and variables PRIVATE except for INIT_TPCORE and 
# TPCORE_FVDAS, which need to be seen outside.  (bdf, bmy, 5/7/03) 
#   # 
# export tpcore_window, init_window, exit_tpcore_window 
#   # 

# @static if @isdefined(SPMD) 
# define PRT_PREFIX  if gid==0) 
# @static if @isdefined(PILGRIM) 
# use decompmodule, only: decomptype 
# use ghostmodule, only: ghosttype 
# use parutilitiesmodule, only: gid, parpatterntype, parbegintransfer, parendtransfer 
# type(parpatterntype) :: pattern2dmg, pattern2dng 
# type(decomptype) :: decomp2d 
# type(ghosttype)  :: ghost2dmg, ghost2dng 
# else 
# use mod_comm, only: gid, mp_barrier, mp_send3d_ns, mp_recv3d_ns, mp_allgather1d 
# end 
# else 
# define PRT_PREFIX 
# end 

# real , ALLOCATABLE, save :: dtdx5(:) 
# real , ALLOCATABLE, save :: dtdy5(:) 
# real , ALLOCATABLE, save :: cosp(:) 
# real , ALLOCATABLE, save :: cose(:) 
# real , ALLOCATABLE, save ::  gw(:) 
# real , ALLOCATABLE, save :: rgw(:) 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added DLAT as allocatable array for PJC pressure fixer 
# (bdf, bmy, 5/7/03) 
#  # 
# REAL, ALLOCATABLE, SAVE :: DLAT(:) 
#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added RGW_25 as allocatable array for ND25 N/S mass flux diagnostic. 
# This accounts for the latitude factor. (bdf, bmy, 10/29/04) 
#  # 
# REAL, ALLOCATABLE, SAVE :: RGW_25(:) 
#  REAL, ALLOCATABLE, SAVE :: SINE_25(:) # (dan 0803) 
#  # 

# CONTAINS 

# ------------------------------------------------------------------------- 
# function init_window(state_grid, im, jm, km, jfirst, jlast, ng, mg, dt, ae, clat) 
# ------------------------------------------------------------------------- 

# @static if @isdefined(SPMD) 
# @static if @isdefined(PILGRIM) 
# use decompmodule, only : decompcreate 
# use ghostmodule, only : ghostcreate 
# use parutilitiesmodule, only : gid, gsize, commglobal, parpatterncreate 
# else 
# use mod_comm, only : gid, y_decomp 
# end 
# end 
# USE State_Grid_Mod, ONLY : GrdState 

# ------- 
# Input 
# ------- 

#  TYPE(GrdState), INTENT(IN) :: State_Grid # Grid State object 
#  integer, intent(in):: im # Global E - W dimension 
#  integer, intent(in):: jm # Global N - S dimension 
#  integer, intent(in):: km # Vertical dimension 
#  integer, intent(out):: jfirst # Local first index for N - S 
#  integer, intent(out):: jlast # Local last index for N - S 
#  integer, intent(in):: ng # large ghost width 
#  integer, intent(in):: mg # small ghost width 

#  real, intent(in):: dt # Time step in seconds 
#  real, intent(in):: ae # Earth's radius (m) 
#  real, intent(in):: clat(0:jm + 1) # latitude in radian  (dan) 

# ----- 
# Local 
# ----- 

#  elat::AbstractFloat = zeros(jm + 1) # cell edge latitude in radian 
# sine::AbstractFloat = zeros(jm + 1) 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Comment out this declaration of DLAT.  DLAT has now been declared 
# as ALLOCATABLE for use in TPCORE_FVDAS. (bdf, bmy, 5/7/03) 
#  # 
#  #dlat::AbstractFloat = zeros(jm) # delta - latitude in radian 
#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Add SINE_25 array as a local variable.  This is used to initialize 
# the RGW_25 array, which is necessary for the ND25 diagnostic. 
# (bdf, bmy, 10/29/04) 
#  # 
# REAL SINE_25(JM + 1) 
#  # 

# dlon::AbstractFloat 
# pi::AbstractFloat 
# i, j::Integer = 0 

# @static if @isdefined(SPMD) 
# @static if @isdefined(PILGRIM) 
# integer, allocatable :: xdist(:), ydist(:) 

# allocate(xdist(1)) 
# allocate(ydist(gsize)) 
# 
# Decomposition 
# 
# xdist(1) = im 
# call newdecomp(jm, gsize, ydist) 

# jfirst = 1 
# for i ∈ 1:gid 
# jfirst = jfirst + ydist(i) 
# end 
# jlast = jfirst + ydist(gid + 1) - 1 

#   # call decompcreate(1, gsize, xdist, ydist, decomp2d ) # 2D region with 1D lat decomposition 
# call ghostcreate(decomp2d, gid, im, 1, im, .false., jm, jfirst - mg, jlast + mg, .false.,  ghost2dmg ) 
# call ghostcreate(decomp2d, gid, im, 1, im, .false., jm, jfirst - ng, jlast + ng, .false.,  ghost2dng ) 
# call parpatterncreate(commglobal, ghost2dmg, pattern2dmg) 
# call parpatterncreate(commglobal, ghost2dng, pattern2dng) 
# else 
# 
# Default decomposition 
# 
# call y_decomp(jm, km, jfirst, jlast, 1, km, gid) 
# end 
# else 
# jfirst = 1 
# jlast  = jm 
# end 

# if  jlast - jfirst < 2  
# println("Minimum size of subdomain is 3") 
# end 

# ---------------- 
# Allocate arrays 
# ---------------- 

# allocate ( cosp(jm) ) 
# allocate ( cose(jm + 1) ) 
# allocate (   gw(jm) ) 
# allocate (  rgw(jm) ) 
# allocate ( dtdx5(jm) ) 
# allocate ( dtdy5(jm) ) 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# We must allocate DLAT here for PJC pressure fixer (bdf, bmy, 5/7/03) 
#  # 
# allocate ( DLAT(JM) ) 
#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# We must allocate RGW_25 here for ND25 N/S transport fluxes diagnostic. 
# This accounts for the latitude factor.  (bdf, bmy, 10/29/04) 
#  # 
# allocate ( RGW_25(JM) ) 
#  # 
#   # allocate ( SINE_25(JM + 1) ) # (dan 0803) 

# pi = 4.0 * atan(1.0) 

# if TRIM(State_GridGridRes) == "0.25x0.3125"  
#     dlon = 2.0*pi / float(1152) # (dan) 
# elseif TRIM(State_GridGridRes) == "0.5x0.625"  
#     dlon = 2.0*pi / float(576) # (dan) 
# end 

# dan for window 
#     #elat(1) = -0.5*pi # S. Pole 
# sine(1) = -1.0 
#     # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#     # 
# Initialize SINE_25 array (bmy, bdf, 10/29/04) 
#     # 
# SINE_25(1) = -1.0 
#     # 
# cose(1) =  0.0 



#  for j ∈ 1:jm + 1 # (dan) 
# elat(j) = 0.5*(clat(j - 1) + clat(j)) 
# sine(j) = sin(elat(j)) 
#     # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#     # 
# Initialize SINE_25 array (bmy, bdf, 10/29/04) 
#     # 
# SINE_25(J) = SIN( CLAT(J) ) 
#     # 
# cose(j) = cos(elat(j)) 
# end 

# dan for window 
#     #elat(jm + 1) = 0.5*pi # N. Pole 
# sine(jm + 1) = 1.0 
#     # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#     # 
# Initialize SINE_25 array (bmy, bdf, 10/29/04) 
#     # 
# SINE_25(JM + 1) = 1.0 
#     # 



#     #dlat(1) = 2.0*(elat(2) - elat(1)) # Polar cap  (dan) 

# for j ∈ 1:jm 
# dlat(j) = elat(j + 1) - elat(j) 
# end 
#     #dlat(jm) = 2.0*(elat(jm + 1) - elat(jm)) # Polar cap (dan) 

# for j ∈ 1:jm 
# gw(j) = sine(j + 1) - sine(j) 
# cosp(j) = gw(j) / dlat(j) 
# rgw(j) =  1.0 / gw(j) 
#     # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#     # 
# Initialize RGW_25 for the ND25 N/S transport fluxes diagnostic. 
# RGW_25 takes into account the latitude factor. (bdf, bmy, 10/29/04) 
#     # 
# RGW_25(J) = 1.0 / ( SINE_25(J + 1) - SINE_25(J) ) 
#     # 
# dtdx5(j) = 0.5 * dt / (dlon*ae*cosp(j)) 
# dtdy5(j) = 0.5 * dt / (ae*dlat(j)) 
# end 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
# Comment out the === lines (bmy, 7/20/04) 
#  ## Now use REPEAT cmd (bmy, 4/29/03) 
# PRT_PREFIX write( 6, "(a)" ) REPEAT( "=", 79 ) 
#  # 
# prt_prefix println("NASA - GSFC Tracer Transport Module successfully initialized") 
#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
# Comment out the === lines (bmy, 7/20/04) 
# PRT_PREFIX write( 6, "(a)" ) REPEAT( "=", 79 ) 
#  # 

# end 

# ------------------------------------------------------------------------- 
# function exit_tpcore_window() 
# ------------------------------------------------------------------------- 

# @static if @isdefined(SPMD) && defined(PILGRIM) 
# use decompmodule, only : decompfree 
# use ghostmodule, only : ghostfree 
# use parutilitiesmodule, only : commglobal, parpatternfree 

# call parpatternfree(commglobal, pattern2dmg) 
# call parpatternfree(commglobal, pattern2dng) 
# call ghostfree(ghost2dmg) 
# call ghostfree(ghost2dng) 
# call decompfree(decomp2d) 
# end 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Comment out original code below (bdf, bmy, 5/9/03) 
#  # 
# deallocate ( cosp ) 
# deallocate ( cose ) 
# deallocate (   gw ) 
# deallocate (  rgw ) 
# deallocate ( dtdx5 ) 
# deallocate ( dtdy5 ) 
#  # 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Now deallocate arrays only if they have been allocated (bmy, 5/9/03) 
# Also deallocate RGW_25 array (bdf, bmy, 10/29/04) 
#  # 
# if ( ALLOCATED( COSP   ) ) DEALLOCATE( COSP   ) 
# if ( ALLOCATED( COSE   ) ) DEALLOCATE( COSE   ) 
# if ( ALLOCATED( GW     ) ) DEALLOCATE( GW     ) 
# if ( ALLOCATED( RGW    ) ) DEALLOCATE( RGW    ) 
#  if ( ALLOCATED( RGW_25 ) ) DEALLOCATE( RGW_25 ) # (bdf, bmy, 10/29/04) 
# if ( ALLOCATED( DTDX5  ) ) DEALLOCATE( DTDX5  ) 
# if ( ALLOCATED( DTDY5  ) ) DEALLOCATE( DTDY5  ) 
# if ( ALLOCATED( DLAT   ) ) DEALLOCATE( DLAT   ) 
#  if ( ALLOCATED( SINE_25) ) DEALLOCATE( SINE_25) # (dan 0803) 
#  # 


# end 


# ---------------------------------------------------------------------------- 
# function TPCORE_WINdoW( dt,    ae, im, jm, km, jfirst, jlast, ng, mg, nq,    ak, bk, u, v, ps1, ps2, ps,  q, iord, jord, kord, n_adj, # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added XMASS, YMASS arguments to arg list of TPCORE_FVDAS for PJC/LLNL 
# pressure fixer (bdf, bmy, 5/7/03) 
#  # 
#                            XMASS,    YMASS, # 
# AREA_M2, State_Chm, State_Diag   ) 
# ---------------------------------------------------------------------------- 

# Uses: 
#     # USE PhysConstants # Physical constants g0_100 and AIRMW 
# USE State_Chm_Mod,  ONLY : ChmState 
# USE State_Diag_Mod, ONLY : DgnState 



# Input: 
#  integer, intent(in):: im # Global E - W dimension 
#  integer, intent(in):: jm # Global N - S dimension 
#  integer, intent(in):: km # Vertical dimension 
#  integer, intent(in):: jfirst # Local first index for N - S 
#  integer, intent(in):: jlast # Local last index for N - S 
#  integer, intent(in):: ng # Primary ghost region 
#  integer, intent(in):: mg # Secondary ghost region 
#  integer, intent(in):: nq # Ghosted latitudes (3 required by PPM) 
#  integer, intent(in):: iord # E - W transport scheme 
#  integer, intent(in):: jord # N - S transport scheme 
#  integer, intent(in):: kord # Vertical mapping scheme 
#  integer, intent(in):: n_adj # Number of adjustemnt to air_mass_flux 
# 0 --> no adjustment 

# Recommended values : iord=jord=4, kord=7 
# _ord: 
# --------------------------------------------------------------------------- 
# 1: 1st order upstream scheme 
# 2: 2nd order van Leer (full monotonicity constraint; 
# see Lin et al 1994, MWR) 
# 3: Standard monotonic PPM* (Collela & Woodward 1984) 
# 4: New & Improved monotonic PPM 
# 5: positive - definite PPM (constraint on the subgrid distribution is 
# only strong enough to prevent generation of negative values; 
# both overshoots & undershootes are possible). 
# 6: un - constrained PPM (nearly diffusion free; faster but 
# positivity of the subgrid distribution is not quaranteed. 
# 7: Huynh/Van Leer/Lin full monotonicity constraint 
# --------------------------------------------------------------------------- 
# Only kord can be set to 7 to enable the use of Huynh's 2nd monotonicity 
# constraint for piece - wise parabolic distribution. 
# *PPM: Piece - wise Parabolic Method 

#  real, intent(in):: ak(km + 1) # See below 
#  real, intent(in):: bk(km + 1) # See below 
#  real, intent(in):: u(:,:,:) # u - wind (m/s) at mid - time-level (t=t + dt/2) 
#  real, intent(inout):: v(:,:,:) # v - wind (m/s) at mid - time-level (t=t + dt/2) 

# ------------------------------------------------------ 
# The hybrid ETA - coordinate: 
# pressure at layer edges are defined as follows: 
# 
# p(i, j, k) = ak(k) + bk(k)*ps(i, j) 
# ------------------------------------------------------ 
# ak and bk are defined at layer edges. 

# ///////////////////////////////// 
# / \ ------ Model top P=ak(1) --------- ak(1), bk(1) 
# | 
# delp(1)    |  ........... q(i, j, 1) ............ 
# | 
# \ / ---------------------------------  ak(2), bk(2) 
# 
# 
# 
# / \ ---------------------------------  ak(k), bk(k) 
# | 
# delp(k)    |  ........... q(i, j, k) ............ 
# | 
# \ / ---------------------------------  ak(k + 1), bk(k + 1) 
# 
# 
# 
# / \ ---------------------------------  ak(km), bk(km) 
# | 
# delp(km)   |  ........... q(i, j, km) ......... 
# | 
# \ / -----Earth's surface P=Psfc ------ ak(km + 1), bk(km + 1) 
# ////////////////////////////////// 


# Note: surface pressure can be of any unit (e.g., pascal or mb) as long as it is 
# consistent with the definition of (ak, bk) defined above 
# Winds (u, v), ps, and q are assumed to be defined at the same points. 
# The latitudes are given by clat, input to the initialization routine: init_tpcore. 

#  real, intent(in):: ps1(im, jfirst:jlast) # surface pressure at current time 
#  real, intent(in):: ps2(im, jfirst:jlast) # surface pressure at future time=t + dt 
#  real, intent(in):: dt # Transport time step in seconds 
#  real, intent(in):: ae # Earth's radius (m) 

#  real, intent(inout):: q(:,:,:,:) # Tracer "mixing ratios" 
# q could easily be re - dimensioned 

#  real, intent(out):: ps(im, jfirst:jlast) # "predicted" surface pressure 

#  delp::AbstractFloat = zeros(im, jfirst:jlast, km) # Predicted thickness at future time (t=t + dt) 
#  pe::AbstractFloat = zeros(im, km + 1, jfirst:jlast) # Pressure at layer edges (predicted) 

#  fx::AbstractFloat = zeros(im, jfirst:jlast, km) # E - W air mass flux 
#  va::AbstractFloat = zeros(im, jfirst:jlast, km) # N - S CFL at cell center (scalar points) 

# Added XMASS, YMASS for the PJC pressure - fixer (bdf, bmy, 5/7/03) 
# REAL,    INTENT(IN)    :: XMASS(:,:,:), YMASS(:,:,:) 

#  REAL,    INTENT(IN)    :: AREA_M2(JM) # box area for mass flux diag 

# Chemistry and diagnostic state objects 
# TYPE(ChmState), INTENT(INOUT) :: State_Chm 
# TYPE(DgnState), INTENT(INOUT) :: State_Diag 

#  # 

# ----------------------- 
# Ghosted local arrays: 
# ----------------------- 

#  cx::AbstractFloat = zeros(im, jfirst - ng:jlast + ng, km) # E - W CFL number on C - grid 
#  delp1::AbstractFloat = zeros(im, jfirst - mg:jlast + mg, km) # Pressure thickness at current time (t) 
#  fy::AbstractFloat = zeros(im, jfirst:jlast + mg, km) # N - S air mass flux 
#  cy::AbstractFloat = zeros(im, jfirst:jlast + mg, km) # N - S CFL number on C - grid 
#  psg::AbstractFloat = zeros(im, jfirst - mg:jlast + mg, 2) # Was psm and psn 
#  q2::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # local 2D q array 
#  ffsl::Bool = zeros(jfirst - ng:jlast + ng, km) # Flag to compute Integer fluxes 

# Local variables: 
# i, j, k, iq::Integer = 0 
#  iord_bg::Integer # E - W scheme for background mass flux 
#  jord_bg::Integer # N - S scheme for background mass flux 
#  js1gd::Integer # Southern latitude border (1 on SP PE) 
#  jn1gd::Integer # Northern latitude border (JM on NP PE) 
#  nx::Integer # Internal E - W OpenMP decomposition 
#  iv::Integer # Monotonicity constraints for top and bottom 

# Binary diagnostics are retired but netcdf needs implementation. 
# Comment out for now (ewl, 3/10/22) 
# 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
# 
# Define DTC, QTEMP, TRACE_DIFF for ND26 diagnostic (bdf, bmy, 9/28/04) 
# 
# DTC::AbstractFloat = zeros(IM, JM, KM, NQ) # up/down flux temp array 
# QTEMP::AbstractFloat = zeros(IM, JM, KM, NQ) # up/down flux array 
# REAL TRACE_DIFF # up/down flux variable 
# 

#  const nx = 1 # Try 2 or 4 if large number of OMP threads 
# is to be used 

#  js1gd = max(1,  jfirst - ng) # NG latitudes on S (starting at 1) 
#  jn1gd = min(jm, jlast + ng) # NG latitudes on N (ending at jm) 

# Eloise Marais questioned why iord_bg and jord_bg are set to 1 instead of 
# being set to the iord and jord values from geoschem_config.yml.  I looked back 
# into the Git history and couldn't find where this change was made, so 
# this issue has probably persisted for quite some time.  We will set these 
# values to iord and jord in an attempt to avoid negative values. 
# See Github issue: https://github.com/geoschem/geos - chem/issues/840.0 
# -- Bob Yantosca (31 Aug 2021) 
# iord_bg = iord 
# jord_bg = jord 

#  iv = 1 # Enforce strong constraint at top & bottom 
# May want to change to iv=0 if diffusion is a problem 

#  # iv =0 # (dan.iv 0803) 
#  # iv=-1 # (dan 0803) 


# Ensure inputs are single - valued at poles: 


# for j ∈ jfirst:jlast 
# for i ∈ 1:im 

# psg(i, j, 1) = ps1(i, j) 
# psg(i, j, 2) = ps2(i, j) 
# end 
# end 

# if  jfirst == 1  
# call xpavg(psg(:, 1, 1), im) 
# call xpavg(psg(:, 1, 2), im) 
# end 

# if  jlast == jm  
# call xpavg(psg(:, jm, 1), im) 
# call xpavg(psg(:, jm, 2), im) 
# end 

# @static if @isdefined(SPMD) 
# Ghost v, psm and psn north/south --> now in one array psg 
# @static if @isdefined(PILGRIM) 
# call parbegintransfer(pattern2dmg, km, v, v) 
# call parbegintransfer(pattern2dmg, 2, psg, psg) 
# else 
# call mp_send3d_ns(im, jm, jfirst, jlast, 1, km, mg, mg, v, 1) 
# call mp_send3d_ns(im, jm, jfirst, jlast, 1, 2, mg, mg, psg, 2) 
# end 
# end 

# Average q at both poles 
# for iq ∈ 1:nq 
##$omp parallel do  ##$omp shared(im)  ##$omp private(k) 
# for k ∈ 1:km 
# if  jfirst == 1  
# call xpavg(q(:, 1, k, iq), im) 
# end 
# if  jlast == jm  
# call xpavg(q(:, jm, k, iq), im) 
# end 
# end 
# end 

# @static if @isdefined(SPMD) 
# @static if @isdefined(PILGRIM) 
# call parendtransfer(pattern2dmg, km, v, v) 
# call parendtransfer(pattern2dmg, 2, psg, psg) 
# else 
# call mp_barrier() 
# call mp_recv3d_ns(im, jm, jfirst, jlast, 1, km, mg, mg, v, 1) 
# call mp_recv3d_ns(im, jm, jfirst, jlast, 1, 2, mg, mg, psg, 2) 
# call mp_barrier() 
# end 
# end 

# ---------------------------------------------- 
# Compute background air mass fluxes 
# ---------------------------------------------- 

#   # call air_mass_flux(im, jm, km, jfirst, jlast, iord_bg, jord_bg,   ak, bk, psg, ps,  u, v, cx, cy, va, fx, fy, ng, mg, ffsl, delp1, delp, pe,  dt, ae, n_adj, # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added XMASS, YMASS to the arg list of AIR_MASS_FLUX 
# for the PJC/LLNL pressure - fixer (bdf, bmy, 5/7/03) 
#  # 
# XMASS, YMASS ) 
#  # 


# --------------------------------------------------- 
# Do tracer transport 
# --------------------------------------------------- 

# @static if @isdefined(SPMD) 
# non - blocking-send for tracer # 1 
# @static if @isdefined(PILGRIM) 
# call parbegintransfer(pattern2dng, km, q(:,:,:, 1), q(:,:,:, 1)) 
# else 
# call mp_send3d_ns(im, jm, jfirst, jlast, 1, km, ng, ng, q(1, jfirst - ng, 1, 1), 1) 
# end 
# end 


# Multi_Tracer: 
# for iq ∈ 1:nq 

# @static if @isdefined(SPMD) 
# Receive current tracer 
# @static if @isdefined(PILGRIM) 
# call parendtransfer(pattern2dng, km, q(:,:,:, iq), q(:,:,:, iq)) 
# if iq < nq 
# call parbegintransfer(pattern2dng, km, q(:,:,:, iq + 1), q(:,:,:, iq + 1)) 
# end 
# else 
# call mp_barrier() 
# call mp_recv3d_ns(im, jm, jfirst, jlast, 1, km, ng, ng, q(1, jfirst - ng, 1, iq), iq) 
# call mp_barrier() 
# if  iq < nq  
# non - blocking send for next tracer 
# call mp_send3d_ns(im, jm, jfirst, jlast, 1, km, ng, ng, q(1, jfirst - ng, 1, iq + 1), iq + 1) 
# end 
# end 
# end 

#$omp parallel do  #$omp default( shared ) # $omp private( i, j, k, q2 ) 

# Vertical_OMP: 

# for k ∈ 1:km 


# q2(:,:) = 0.0 

# Copying q to 2d work array for transport. This allows q to be dimensioned 
# differently from the # calling routine. 

# for j ∈ js1gd:jn1gd 
# for i ∈ 1:im 
# q2(i, j) = q(i, j, k, iq) 
# end 
# end 

# --- Previous to (ccc, 9/9/10) 
#     # call tp2g( q2(1, jfirst - ng),    va(1, jfirst, k),  #               cx(1, jfirst - ng, k),  cy(1, jfirst, k),  #               im,  jm,  iv,   iord,     jord,  #               ng,  mg,  fx(1, jfirst, k), fy(1, jfirst, k),  #               ffsl(jfirst - ng, k),    jfirst,   jlast,  #               delp1(1, jfirst - mg, k),    delp(1, jfirst, k),  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
# 
# Now pass MASSFLEW, MASSFLNS, AREA_M2, TCVV, ND24, ND25, DT as 
# arguments to routine TP2G for GEOS - CHEM mass flux diagnostics 
# (bdf, bmy, 9/28/04) 
# 
#               MASSFLEW(1, 1, K, IQ), MASSFLNS(1, 1, K, IQ), # AREA_M2, TCVV(IQ), ND24, ND25, DT ) 
# 

#      # call tp2g( q2(1, jfirst - ng),    va(1, jfirst, k), cx(1, jfirst - ng, k),  cy(1, jfirst, k), im,  jm,  iv,   iord,     jord, ng,  mg,  fx(1, jfirst, k), fy(1, jfirst, k), ffsl(jfirst - ng, k),    jfirst,   jlast, delp1(1, jfirst - mg, k),    delp(1, jfirst, k))#, # Binary diagnostics are retired but netcdf needs implementation. 
# Comment out for now (ewl, 3/10/22) 
# AREA_M2, ND24, ND25, DT ) 

# ------------------------------------------------------------------------------ 
# Prior to 4/1/15: 
# Preserve original code here.  Lin Zhang submitted the fix below. 
# for j ∈ jfirst:jlast 
#    do j=max(jfirst, jord + 1), min(jlast, jm - jord + 1) # Lin_20140518 
# for i ∈ 1:im 
# q(i, j, k, iq) = q2(i, j) 
# end 
# end 
# ------------------------------------------------------------------------------ 
# NOTE: This fix was submitted by Lin Zhang.  Not sure if it supersedes 
# the previous code but we'll put it here for now. (bmy, 4/1/15) 
#     for j ∈ (jfirst + 2):(jlast - 2) # (lzh, 05/10/2014) 
# for i ∈ 3:im - 2 
# q(i, j, k, iq) = q2(i, j) 
# end 
# end 

# end Vertical_OMP 
# end Multi_Tracer 

# end 
# end 

# --------------------------------------------------------------- 
# Perform Remapping back to the hybrid sigma - pressure coordinate 
# Mass will be conserved if predicted ps2 == psn (data/model) 
# --------------------------------------------------------------- 

# Binary diagnostics are retired but netcdf needs implementation. 
# Comment out for now (ewl, 3/10/22) 
#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Save tracer values before vertical transport (bdf, bmy, 9/28/04) 
#  # 
# if ND26 > 0  
# QTEMP = Q 
# end 
#  # 

# call qmap(pe, q, im, jm, km, nx, jfirst, jlast, ng, nq, ps, ak, bk, kord, iv) 

# Binary diagnostics are retired but netcdf needs implementation. 
# Comment out for now (ewl, 3/10/22) 
# 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
# 
# Implement ND26 diag: Up/down flux of tracer [kg/s] (bmy, bdf, 9/28/04) 
# 
# The vertical transport done in qmap.  We need to find the difference 
# in order to to interpret transport. 
# 
# Break up diagnostic into up & down fluxes using the surface boundary 
# conditions.  Start from top down (really surface up for flipped TPCORE) 
# 
# if ND26 > 0  
# 
# ----------------- 
# start with top 
# ----------------- 
# K = 1 
# 
##$OMP PARALLEL do  ##$OMP DEFAULT( SHARED )  ##$OMP PRIVATE( I, J, IQ ) 
# for IQ ∈ 1:NQ 
# for I ∈ 1:IM 
# for J ∈ 1:JM 
# 
#       DTC(I, J, K, IQ) = ( Q(I, J, K, IQ)     * DELP1(I, J, K)   -  #                         QTEMP(I, J, K, IQ) * DELP(I, J, K)  ) * # AREA_M2(J) * g0_100 
# 
# top layer should have no residual.  the small residual is from 
# a non - pressure fixed flux diag.  The z direction may be off by 
# a few percent. 
# MASSFLUP(I, J, K, IQ) = MASSFLUP(I, J, K, IQ) + DTC(I, J, K, IQ)/dt 
# end 
# end 
# end 
##$OMP end PARALLEL do 
# 
# ---------------------------------------------------- 
# get the other fluxes using a mass balance equation 
# ---------------------------------------------------- 
# for K ∈ 2:KM 
##$OMP PARALLEL do  ##$OMP DEFAULT( SHARED )  ##$OMP PRIVATE( I, J, IQ, TRACE_DIFF ) 
# for IQ ∈ 1:NQ 
# for I ∈ 1:IM 
# for J ∈ 1:JM 
# 
#       TRACE_DIFF         = ( Q(I, J, K, IQ)     * DELP1(I, J, K)  -  #                              QTEMP(I, J, K, IQ) * DELP(I, J, K) ) * # AREA_M2(J) * g0_100 
# 
# DTC(I, J, K, IQ)      = DTC(I, J, K - 1, IQ) + TRACE_DIFF 
# 
# end 
# end 
# end 
##$OMP end PARALLEL do 
# end 
# 
# end 

# end 


# function air_mass_flux(im, jm, km, jfirst, jlast, iord, jord, ak, bk, psg, ps, u, v, cx, cy, va, fx, fy, ng,  mg,  ffsl, delp1,  delp, pe, dt, ae,  n_adj, # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added XMASS, YMASS to the arg list of AIR_MASS_FLUX 
# for the PJC/LLNL pressure - fixer (bdf, bmy, 5/7/03) 
#  # 
# XMASS, YMASS ) 
#  # 

# ------------------------------------------------------ 
# The hybrid ETA - coordinate: 
# pressure at layer edges are defined as follows: 
# 
# p(i, j, k) = ak(k) + bk(k)*ps(i, j)          (1) 
# ------------------------------------------------------ 

# Input from Data/Model: 
# (u, v) is the time mean wind at Time=t + dt/2 
# delp1 is the layer thickness at Time=t 

# Output: 
# delp is the predicted thickness at Time=t + dt 
# (fx, fy): background air mass flxues 
# (cx, cy): CFL number 



# integer, intent(in):: im 
# integer, intent(in):: jm 
# integer, intent(in):: km 
# integer, intent(in):: jfirst 
# integer, intent(in):: jlast 
# integer, intent(in):: iord 
# integer, intent(in):: jord 
# integer, intent(in):: ng 
# integer, intent(in):: mg 
# integer, intent(in):: n_adj 

# real, intent(in):: dt 
# real, intent(in):: ae 
# real, intent(in):: ak(km + 1) 
# real, intent(in):: bk(km + 1) 
#  real, intent(in):: psg(im, jfirst - mg:jlast + mg, 2) # Was ps1 and ps2 
# real, intent(in):: u(im, jfirst:jlast, km) 
# real, intent(in):: v(im, jfirst - mg:jlast + mg, km) 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added XMASS, YMASS for PJC/LLNL pressure fixer (bdf, bmy, 5/7/03) 
#  # 
# REAL, INTENT(IN) :: XMASS(IM, JM, KM), YMASS(IM, JM, KM) 
#  # 

# Output: 
# logical, intent(out):: ffsl(jfirst - ng:jlast + ng, km) 
# real, intent(out):: cx(im, jfirst - ng:jlast + ng, km) 
# real, intent(out):: delp(im, jfirst:jlast, km) 

# real, intent(out):: ps(im, jfirst:jlast) 
# real, intent(out):: fx(im, jfirst:jlast, km) 
# real, intent(out):: cy(im, jfirst:jlast + mg, km) 
# real, intent(out):: fy(im, jfirst:jlast + mg, km) 
# real, intent(out):: va(im, jfirst:jlast, km) 

# real, intent(out):: delp1(im, jfirst - mg:jlast + mg, km) 

# real, intent(out):: pe(im, km + 1, jfirst:jlast) 

# Local: 
# yms::AbstractFloat = zeros(im, jfirst:jlast + mg, km) 

# tiny::AbstractFloat 
# const tiny = 1.0e-10 
# dak, dbk::AbstractFloat = 0 
# dtoa, vt::AbstractFloat = 0 
# i, j, k::Integer = 0 

# js2g0::Integer 
# jn2g0::Integer 
# jn1g1::Integer 
# js2gd, jn2gd::Integer = 0 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Declare extra variables PJC/LLNL pressure fixer (bdf, bmy, 5/7/03) 
#  # 
# REAL :: DELPM(IM, JM, KM), FACTY, UT 
#  # 

#  js2g0  = max(2, jfirst) # No ghosting 
#  jn2g0  = min(jm - 1, jlast) # No ghosting 
#  jn1g1  = min(jm, jlast + 1) # Ghost 1 on N 
#  js2gd = max(2,  jfirst - ng) # NG latitudes on S (starting at 1) 
#  jn2gd = min(jm - 1, jlast + ng) # NG latitudes on N (ending at jm - 1) 

# dtoa = .5*dt/ae 

# cx(:,:,:)=0D0 
# cy(:,:,:)=0D0 
# fx(:,:,:)=0D0 
# fy(:,:,:)=0D0 
# delp(:,:,:)=0D0 
# ps(:,:)=0D0 
# va(:,:,:)=0D0 
# delp1(:,:,:)=0D0 
# pe(:,:,:)=0D0 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Define DELPM for PJC pressure fixer (bdf, bmy, 5/7/03) 
#  # 
# for K ∈ 1:KM 
# for J ∈ 1:JM 
# for I ∈ 1:IM 
# DELPM(I, J, K) = ( AK(K + 1) - AK(K) ) + ( BK(K + 1) - BK(K) ) * ( 0.5d0 * ( PSG(I, J, 1) + PSG(I, J, 2 ) + 2d0 * AK(1) ) ) 
# end 
# end 
# end 
#  # 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added for PJC/LLNL pressure - fixer (bdf, bmy, 5/7/03) 
# Note that DTDY5 is the same everywhere except at the poles, so 
# we can just pick a value roughly close to the equator 
#  # 
# FACTY = DTDY5(JM/2) 
#  # 

# $omp parallel do private(i, j, k, vt, UT ) 

# for k ∈ 1:km 

# for j ∈ js2g0:jn1g1 
# for i ∈ 1:im 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Change calculation of VT for PJC pressure fixer (bdf, bmy, 5/7/03) 
#  # 
# VT = YMASS(I, J, K) / FACTY / COSE(J) / DELPM(I, J, K) + V(I, J - 1, K) * ( 1d0 - DELPM(I, J - 1, K) / DELPM(I, J, K) ) 
#  # 

# if  vt > 0.0  
# cy(i, j, k) = dtdy5(j - 1)*vt 
# else 
# cy(i, j, k) = dtdy5(j)*vt 
# end 
# yms(i, j, k) = dtoa*vt*cose(j) 
# end 
# end 

# for j ∈ js2g0:jn2g0 
# for i ∈ 1:im 
# if  cy(i, j, k)*cy(i, j + 1, k) > 0.0  
# if  cy(i, j, k) > 0.0  
# va(i, j, k) = cy(i, j, k) 
# else 
# va(i, j, k) = cy(i, j + 1, k) 
# end 
# else 
# va(i, j, k) = 0.0 
# end 
# end 
# end 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Removed this section for PJC pressure fixer (bdf, bmy, 5/7/03) 
# for j ∈ js2g0:jn2g0 
# cx(1, j, k) = dtdx5(j)*(u(1, j, k)+u(im, j, k)) 
# for i ∈ 2:im 
# cx(i, j, k) = dtdx5(j)*(u(i, j, k)+u(i - 1, j, k)) 
# end 
# end 
#  # 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Added this section for PJC pressure fixer (bdf, bmy, 5/7/03) 
#  # 
# for J ∈ JS2G0:JN2G0 
# UT        = XMASS(1, J, K) / DTDX5(J) / DELPM(1, J, K) + U(IM, J, K) * ( 1d0 - DELPM(IM, J, K) / DELPM(1, J, K) ) 
# CX(1, J, K) = DTDX5(J) * UT 

# for I ∈ 2:IM 
# UT        = XMASS(I, J, K) / DTDX5(J) / DELPM(I, J, K) + U(I - 1, J, K) * ( 1d0 - DELPM(I - 1, J, K) / DELPM(I, J, K) ) 
# CX(I, J, K) = DTDX5(J) * UT 
# end 
# end 
#  # 

# end 


# @static if @isdefined(SPMD) 
# No buffer version (km # calls to mpi_sendrecv) 
# @static if @isdefined(PILGRIM) 
# call parbegintransfer(pattern2dng, km, cx, cx) 
# call parendtransfer(pattern2dng, km, cx, cx) 
# else 
# call mp_send3d_ns(im, jm, jfirst, jlast, 1, km, ng, ng, cx, 3) 
# call mp_barrier() 
# call mp_recv3d_ns(im, jm, jfirst, jlast, 1, km, ng, ng, cx, 3) 
# call mp_barrier() 
# end 
# end 

# --------------------------------------------------- 
# Compute background mass - flux (fx, fy) and (cx, cy) 
# --------------------------------------------------- 

#$omp parallel do  #$omp shared(im, jm, iord, jord, mg, jfirst, jlast) # $omp private(i, j, k, dak, dbk) 

# for k ∈ 1:km 

#      for j ∈ js2gd:jn2gd # ffsl needed on N*ng S*ng 
# ffsl(j, k) = .false. 
# for i ∈ 1:im 
# if  abs(cx(i, j, k)) > 1.0  
# ffsl(j, k) = .true. 
# go to 2222 
# end 
# end 
# 2222  continue 
# end 

# dak = ak(k + 1) - ak(k) 
# dbk = bk(k + 1) - bk(k) 

# do j=max(1, jfirst - mg), min(jm, jlast + mg) 
# for i ∈ 1:im 
# delp1(i, j, k) = dak + dbk*psg(i, j, 1) 
# end 
# end 

# call tp2d(va(1, jfirst, k), delp1(1, jfirst - mg, k), cx(1, jfirst - mg, k), cy(1, jfirst, k), im, jm, iord, jord, mg,  mg, fx(1, jfirst, k), fy(1, jfirst, k),  ffsl(jfirst - mg, k), cx(1, jfirst, k), yms(1, jfirst, k), 0, jfirst, jlast) 

#  # 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
#  # 
# Fix mass fluxes in regions not over the courant limit (bdf, bmy, 5/7/03) 
#  # 
# for J ∈ 4:JM - 4 
# FX(:, J, K) = XMASS(:, J, K) 
# FY(:, J, K) = YMASS(:, J, K) * DLAT(J) 
# end 
#  # 

# for j ∈ js2g0:jn2g0 
# for i ∈ 1:im - 1 
# delp(i, j, k) = delp1(i, j, k) + fx(i, j, k) - fx(i + 1, j, k) + (fy(i, j, k)-fy(i, j + 1, k))*rgw(j) 
# end 
# delp(im, j, k) = delp1(im, j, k) + fx(im, j, k) - fx(1, j, k) + (fy(im, j, k)-fy(im, j + 1, k))*rgw(j) 
# end 

# if  jfirst ==  1  
# for i ∈ 1:im 
# delp(i, 1, k) = delp1(i, 1, k) - fy(i, 2, k)*rgw(1) 
# end 
# call xpavg(delp(:, 1, k), im) 
# end 

# if  jlast == jm  
# for i ∈ 1:im 
# delp(i, jm, k) = delp1(i, jm, k) + fy(i, jm, k)*rgw(jm) 
# end 
# call xpavg(delp(:, jm, k), im) 
# end 

# if  n_adj == 0  
# for j ∈ js2g0:jn2g0 
# if  ffsl(j, k)  
# for i ∈ 1:im 
# fx(i, j, k) = fx(i, j, k)/sign(max(abs(cx(i, j, k)), tiny), cx(i, j, k)) 
# end 
# end 
# end 
# end 

# end 

# -------------- 
# Compute ps: 
# -------------- 

# $omp parallel do private(i, j, k) 

# for j ∈ jfirst:jlast 
# for i ∈ 1:im 
# pe(i, 1, j) = ak(1) 
# end 

# for k ∈ 1:km 
# for i ∈ 1:im 
# pe(i, k + 1, j) = pe(i, k, j) + delp(i, j, k) 
# end 
# end 

# for i ∈ 1:im 
# ps(i, j) = pe(i, km + 1, j) 
# end 
# end 

# -------------------------------------------------------------- 
# Apply mass_flux adjuster to nudge predicted ps towards "data" 
# -------------------------------------------------------------- 

# if  n_adj > 0  
# call adj_fx(im, jm, km, jfirst, jlast, ak, bk, ffsl, ps, psg(:,:, 2), pe, delp, fx, cx, fy, ng, mg, tiny, n_adj) 
# end 

# end 

# function tp2g(h,  va, crx, cry, im, jm, iv, iord, jord, ng, mg, xfx, yfx, ffsl, jfirst, jlast, dp, dpp                ) 

# Uses: 
#     USE PhysConstants # Physical constants g0_100 



# INPUT PARAMETERS: 
#    integer, intent(in):: im, jm # Dimensions 
#    integer, intent(in):: jfirst, jlast # Latitude strip 
#    integer, intent(in):: iv # iv=-1 --> vector 
#    integer, intent(in):: iord, jord # Interpolation order in x, y 
#    integer, intent(in):: ng # Max. NS dependencies 
#    integer, intent(in):: mg # Secondary ghosting zones 
#    logical, intent(in):: ffsl(jfirst - ng:jlast + ng) # Use flux - form semi - Lagrangian trans.? 
#    real, intent(in):: va(im, jfirst:jlast) # CFL in y at cell center 
# real, intent(in):: dp(im, jfirst - mg:jlast + mg) 
# real, intent(in):: dpp(im, jfirst:jlast) 

#    real, intent(in):: crx(im, jfirst - ng:jlast + ng) # ( N*NG S*NG ) 
#    real, intent(in):: cry(im, jfirst:jlast + mg) # ( N like FY ) 

#    real, intent(in):: xfx(im, jfirst:jlast) # x - mass flux 
#    real, intent(in):: yfx(im, jfirst:jlast + mg) # y - mass flux 

# real, intent(inout) :: h(im, jfirst - ng:jlast + ng) 

# Binary diagnostics are retired but netcdf needs implementation. 
# Comment out for now (ewl, 3/10/22) 
#   REAL,    INTENT(IN)    :: AREA_M2(JM) # Grid bos surface area [m2] 
#   INTEGER, INTENT(IN)    :: ND24 # flux diag 
#   INTEGER, INTENT(IN)    :: ND25 # flux diag 
#   REAL,    INTENT(IN)    :: DT # time step for flux diagnostic 

# Local 
#    fx::AbstractFloat = zeros(im, jfirst:jlast) # tracer flux in x ( unghosted ) 
#    fy::AbstractFloat = zeros(im, jfirst:jlast + mg) # tracer flux in y ( N, see tp2c ) 

# i, j, js2g0, jn2g0::Integer = 0 
# sum1, DTC::AbstractFloat = 0 

#    js2g0  = max(2, jfirst) # No ghosting 
#    jn2g0  = min(jm - 1, jlast) # No ghosting 

# call tp2d(va, h(1, jfirst - ng), crx(1, jfirst - ng), cry, im, jm, iord, jord, ng, mg, fx, fy, ffsl(jfirst - ng), xfx, yfx, 1, jfirst, jlast) 

# ------------------------------------------------------------------------------ 
# Prior to 4/1/15: 
# Don't treat edges (Lin Zhang, 4/1/15) 
# for j ∈ js2g0:jn2g0 
# for i ∈ 1:im - 1 
# h(i, j) = h(i, j)*dp(i, j) + fx(i, j)-fx(i + 1, j)+(fy(i, j)-fy(i, j + 1))*rgw(j) 
# end 
# end 
# 
# for j ∈ js2g0:jn2g0 
# h(im, j) = h(im, j)*dp(im, j) + fx(im, j)-fx(1, j)+(fy(im, j)-fy(im, j + 1))*rgw(j) 
# end 
# ------------------------------------------------------------------------------ 
# for j ∈ js2g0:jn2g0 
# for i ∈ 2:im - 1 
# h(i, j) = h(i, j)*dp(i, j) + fx(i, j)-fx(i + 1, j)+(fy(i, j)-fy(i, j + 1))*rgw(j) 
# end 
# end 

# Poles 
# if  jfirst == 1  
# for i ∈ 1:im 
# h(i, 1) = h(i, 1)*dp(i, 1) - fy(i, 2)*rgw(1) 
# end 
# call xpavg(h(:, 1), im) 
# end 

# if  jlast == jm  
# for i ∈ 1:im 
# h(i, jm) = h(i, jm)*dp(i, jm) + fy(i, jm)*rgw(jm) 
# end 
# call xpavg(h(:, jm), im) 
# end 

# Binary diagnostics are retired but netcdf needs implementation. 
# Comment out for now (ewl, 3/10/22) 
# 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
# 
# Implement ND24 diag: E/W flux of tracer [kg/s] (bmy, bdf, 9/28/04) 
# 
# (1) H is in units of mixing ratio (input as Q) 
# (2) Unit conversion needs multiply from mixing 
# (airmass/tracer mass)/timestep to get into kg/s 
# (3) DP is current pressure thickness 
# 
# if ND24 > 0  
# for J ∈ JS2G0:JN2G0 
# 
# for I ∈ 1:IM - 1 
# 
# DTC = FX(I, J) * AREA_M2(J) * g0_100 / DT 
# 
# end 
# 
# DTC = FX(IM, J) * AREA_M2(J) * g0_100 / DT 
# 
# end 
# end 
# 
# 
# MODIFICATION by Harvard Atmospheric Chemistry Modeling Group 
# 
# Implement ND25 diag: N/S flux of tracer [kg/s] (bdf, bmy, 9/28/04) 
# Now multiply fluxes by latitude factor RGW_25 (bdf, bmy, 10/29/04) 
# 
# if ND25 > 0  
# for J ∈ JS2G0:JN2G0 
# for I ∈ 1:IM 
# 
# DTC = FY(I, J) * RGW_25(J) * AREA_M2(J) * g0_100 / DT 
# 
# end 
# end 
# 
# South Pole 
# if JFIRST == 1  
# for I ∈ 1:IM 
# 
# DTC = -FY(I, 2) * RGW_25(1) * AREA_M2(1) * g0_100 / DT 
# 
# end 
# end 
# 
# North Pole 
# if JLAST == JM  
# for I ∈ 1:IM 
# 
# DTC = FY(I, JM) * RGW_25(JM) * AREA_M2(JM) * g0_100 / DT 
# 
# end 
# end 
# end 
# 

# ------------------------------------------------------------------- 
# Apply a simple nearest neighbor flux correction to reduce negatives 
# ------------------------------------------------------------------- 
# if  iv != -1  
# call fct_x(h, im, jm, jfirst, jlast, ng, i) 
# end 

# for j ∈ jfirst:jlast 
# for i ∈ 1:im 
# h(i, j) = h(i, j) / dpp(i, j) 
# end 
# end 

# end 

# function tp2d(va, q, crx, cry, im, jm, iord, jord, ng, mg, fx, fy, ffsl, xfx, yfx, id, jfirst, jlast) 



# INPUT PARAMETERS: 
#  integer, intent(in):: im, jm # Dimensions 
#  integer, intent(in):: jfirst, jlast # Latitude strip 
#  iord, jord::Integer = 0 # Interpolation order in x, y 
#  ng::Integer # Max. NS dependencies 
#  mg::Integer # 
#  id::Integer # density (0)  (mfx = C) 
# mixing ratio (1) (mfx = mass flux) 
#  ffsl::Bool = zeros(jfirst - ng:jlast + ng) # Use flux - form semi - Lagrangian trans.? 
# ghosted N*ng S*ng 
#  va::AbstractFloat = zeros(im, jfirst:jlast) # Courant  (unghosted) 
#  q::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # transported scalar ( N*NG S*NG ) 
#  crx::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # Ask S.-J. ( N*NG S*NG ) 
#  cry::AbstractFloat = zeros(im, jfirst:jlast + mg) # Ask S.-J. ( N like FY ) 
#  xfx::AbstractFloat = zeros(im, jfirst:jlast) # Ask S.-J. ( unghosted like FX ) 
#  yfx::AbstractFloat = zeros(im, jfirst:jlast + mg) # Ask S.-J. ( N like FY ) 

# OUTPUT PARAMETERS: 
#  fx::AbstractFloat = zeros(im, jfirst:jlast) # Flux in x ( unghosted ) 
#  fy::AbstractFloat = zeros(im, jfirst:jlast + mg) # Flux in y ( N, see tp2c ) 

# Local: 
# i, j, iad, jp, js2g0, js2gng, jn2g0, jn2gng::Integer = 0 
# adx::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) 
# wk1::AbstractFloat = zeros(im) 
# dm::AbstractFloat = zeros(-im/3:im + im/3) 
# qtmp::AbstractFloat = zeros(-im/3:im + im/3) 
# al::AbstractFloat = zeros(-im/3:im + im/3) 
# ar::AbstractFloat = zeros(-im/3:im + im/3) 
# a6::AbstractFloat = zeros(-im/3:im + im/3) 

# Number of ghost latitudes 
#   js2g0  = max(2, jfirst) # No ghosting 
#   jn2g0  = min(jm - 1, jlast) # No ghosting 
#   js2gng = max(2, jfirst - ng) # Number needed on S 
#   jn2gng = min(jm - 1, jlast + ng) # Number needed on N 

# iad = 1 

#   for j ∈ js2gng:jn2gng # adx needed on N*ng S*ng 

# call xtp(im,  ffsl(j), wk1, q(1, j), crx(1, j), iad, crx(1, j), cosp(j), 0, dm, qtmp, al, ar, a6) 

# for i ∈ 1:im - 1 
# adx(i, j) = q(i, j) + 0.5 * (wk1(i)-wk1(i + 1) + q(i, j)*(crx(i + 1, j)-crx(i, j))) 
# end 
# adx(im, j) = q(im, j) + 0.5 * (wk1(im)-wk1(1) + q(im, j)*(crx(1, j)-crx(im, j))) 
# end 

# if  jfirst == 1  
# for i ∈ 1:im 
# adx(i, 1) = q(i, 1) 
# end 
# end 
# if  jlast == jm  
# for i ∈ 1:im 
# adx(i, jm) = q(i, jm) 
# end 
# end 

# call ytp(im, jm, fy, adx, cry, yfx, ng, mg, jord, 0, jfirst, jlast) 

# for j ∈ js2g0:jn2g0 
# for i ∈ 1:im 
# jp = j - va(i, j) 
# wk1(i) = q(i, j) +0.5*va(i, j)*(q(i, jp)-q(i, jp + 1)) 
# end 

# call xtp(im,  ffsl(j), fx(1, j), wk1, crx(1, j), iord, xfx(1, j), cosp(j), id, dm, qtmp, al, ar, a6) 
# end 
# end 


# function xtp(im, ffsl,  fx,  q,  c,  iord,  mfx, cosa, id, dm, qtmp, al, ar, a6) 



# INPUT PARAMETERS: 
#    id::Integer # ID = 0: density (mfx = C) 
# ID = 1: mixing ratio (mfx is mass flux) 

#    im::Integer # Total longitudes 
#    c::AbstractFloat = zeros(im) # Courant numbers 
# q::AbstractFloat = zeros(im) 
# mfx::AbstractFloat = zeros(im) 
# logical ffsl 
# iord::Integer 
# cosa::AbstractFloat 

# INPUT/OUTPUT PARAMETERS: 
#    qtmp::AbstractFloat = zeros(-im/3:im + im/3) # Input work arrays: 
# dm::AbstractFloat = zeros(-im/3:im + im/3) 
# al::AbstractFloat = zeros(-im/3:im + im/3) 
# ar::AbstractFloat = zeros(-im/3:im + im/3) 
# a6::AbstractFloat = zeros(-im/3:im + im/3) 

# OUTPUT PARAMETERS: 
# fx::AbstractFloat = zeros(im) 

# Local: 
#    cos_upw::AbstractFloat # critical cosine for upwind 
#    cos_van::AbstractFloat # critical cosine for van Leer 
#    cos_ppm::AbstractFloat # critical cosine for ppm 

#    const cos_upw = 0.05 # roughly at 87 deg. 
#    const cos_van = 0.25 # roughly at 75 deg. 
# const cos_ppm = 0.25 

# i, imp::Integer = 0 
# qmax, qmin::AbstractFloat = 0 
# rut, tmp::AbstractFloat = 0 
# iu, itmp, ist::Integer = 0 
# integer isave(im) 
# iuw, iue::Integer = 0 

# imp = im + 1 

# for i ∈ 1:im 
# qtmp(i) = q(i) 
# end 

# if  ffsl  

# Figure out ghost zone for the western edge: 
# iuw =  -c(1) 
# iuw = min(0, iuw) 

# for i ∈ iuw:0 
# qtmp(i) = q(im + i) 
# end 

# Figure out ghost zone for the eastern edge: 
# iue = im - c(im) 
# iue = max(imp, iue) 

# for i ∈ imp:iue 
# qtmp(i) = q(i - im) 
# end 

# if  iord == 1 .or. cosa < cos_upw 
# for i ∈ 1:im 
# iu = c(i) 
# if c(i) <= 0.0 
# itmp = i - iu 
# isave(i) = itmp - 1 
# else 
# itmp = i - iu - 1 
# isave(i) = itmp + 1 
# end 
# fx(i) = (c(i)-iu) * qtmp(itmp) 
# end 
# else 

# for i ∈ 1:im 
# 2nd order slope 
# tmp = 0.25*(qtmp(i + 1) - qtmp(i - 1)) 
# qmax = max(qtmp(i - 1), qtmp(i), qtmp(i + 1)) - qtmp(i) 
# qmin = qtmp(i) - min(qtmp(i - 1), qtmp(i), qtmp(i + 1)) 
# dm(i) = sign(min(abs(tmp), qmax, qmin), tmp) 
# end 

# for i ∈ iuw:0 
# dm(i) = dm(im + i) 
# end 

# for i ∈ imp:iue 
# dm(i) = dm(i - im) 
# end 

# if iord >= 3 && cosa > cos_ppm 
# call fxppm(im, c, mfx, qtmp, dm, fx, iord, al, ar, a6, iuw, iue, ffsl, isave) 
# else 
# for i ∈ 1:im 
# iu  = c(i) 
# rut = c(i) - iu 
# if c(i) .le. 0.0 
# itmp = i - iu 
# isave(i) = itmp - 1 
# fx(i) = rut*(qtmp(itmp)-dm(itmp)*(1.0 + rut)) 
# else 
# itmp = i - iu - 1 
# isave(i) = itmp + 1 
# fx(i) = rut*(qtmp(itmp)+dm(itmp)*(1.0 - rut)) 
# end 
# end 
# end 

# end 

# for i ∈ 1:im 
# if c(i) >= 1.0 
# do ist = isave(i), i - 1 
# fx(i) = fx(i) + qtmp(ist) 
# end 
# elseif c(i) <= -1.0 
# for ist ∈ i:isave(i) 
# fx(i) = fx(i) - qtmp(ist) 
# end 
# end 
# end 

# if id != 0 
# for i ∈ 1:im 
# fx(i) =  fx(i)*mfx(i) 
# end 
# end 

# else 
# Regular PPM (Eulerian without FFSL extension) 

# qtmp(imp) = q(1) 
# qtmp(  0) = q(im) 

# if iord == 1 .or. cosa < cos_upw 
# for i ∈ 1:im 
# iu = float(i) - c(i) 
# fx(i) = mfx(i)*qtmp(iu) 
# end 
# else 

# qtmp(-1)    = q(im - 1) 
# qtmp(imp + 1) = q(2) 

# if iord > 0 .or. cosa < cos_van 
# call xmist(im, qtmp, dm, 2) 
# else 
# call xmist(im, qtmp, dm, iord) 
# end 

# dm(0) = dm(im) 

# if  abs(iord) ==2 .or. cosa < cos_van  
# for i ∈ 1:im 
# iu = float(i) - c(i) 
# fx(i) =  mfx(i)*(qtmp(iu)+dm(iu)*(sign(1.0, c(i))-c(i))) 
# end 
# else 
# call fxppm(im, c, mfx, qtmp, dm, fx, iord, al, ar, a6, iuw, iue, ffsl, isave) 
# end 
# end 

# end 
# end 

# function xmist(im,  q,  dm,  id) 



# INPUT PARAMETERS: 
#  integer, intent(in):: im # Total number of longitudes 
#  integer, intent(in):: id # ID = 0: density (mfx = C) 
# ID = 1: mixing ratio (mfx is mass flux) 
#  real, intent(in):: q(-im/3:im + im/3) # scalar field 

# OUTPUT PARAMETERS: 
#  real, intent(out):: dm(-im/3:im + im/3)   # 

# Local 
# r24::AbstractFloat 
# parameter( r24 = 1.0/24.0) 
# i::Integer 
# qmin, qmax::AbstractFloat = 0 

# if id <= 2 
# for i ∈ 1:im 
# dm(i) = r24*(8.0*(q(i + 1) - q(i - 1)) + q(i - 2) - q(i + 2)) 
# end 
# else 
# for i ∈ 1:im 
# dm(i) = 0.25*(q(i + 1) - q(i - 1)) 
# end 
# end 

# if  id < 0 ) return 

# Apply monotonicity constraint (Lin et al. 1994, MWR) 
# for i ∈ 1:im 
# qmax = max( q(i - 1), q(i), q(i + 1) ) - q(i) 
# qmin = q(i) - min( q(i - 1), q(i), q(i + 1) ) 
# dm(i) = sign( min(abs(dm(i)), qmax, qmin), dm(i) ) 
# end 
# end 

# function fxppm(im, c, mfx,  p, dm, fx, iord, al, ar, a6, iuw, iue, ffsl, isave) 


# INPUT PARAMETERS: 
# integer, intent(in):: im, iord 
# integer, intent(in):: iuw, iue 
# logical, intent(in):: ffsl 
# real, intent(in):: c(im) 
# real, intent(in):: p(-im/3:im + im/3) 
# real, intent(in):: dm(-im/3:im + im/3) 
# real, intent(in):: mfx(im) 

# INPUT/OUTPUT PARAMETERS: 
# integer, intent(inout):: isave(im) 

# real, intent(out):: fx(im) 
# real, intent(out):: al(-im/3:im + im/3) 
# real, intent(out):: ar(-im/3:im + im/3) 
# real, intent(out):: a6(-im/3:im + im/3) 

# LOCAL VARIABLES: 
# r3, r23::AbstractFloat = 0 
# parameter ( r3 = 1.0/3.0, r23 = 2.0/3.0 ) 

# i, lmt::Integer = 0 
# iu, itmp::Integer = 0 
# ru::AbstractFloat 

# for i ∈ 1:im 
# al(i) = 0.5*(p(i - 1)+p(i)) + (dm(i - 1) - dm(i))*r3 
# end 

# for i ∈ 1:im - 1 
# ar(i) = al(i + 1) 
# end 
# ar(im) = al(1) 

# if iord == 7 
# call huynh(im, ar(1), al(1), p(1), a6(1), dm(1)) 
# else 
# if iord == 3 .or. iord == 5 
# for i ∈ 1:im 
# a6(i) = 3.0*(p(i)+p(i)  - (al(i)+ar(i))) 
# end 
# end 
# lmt = iord - 3 
# call lmppm( dm(1), a6(1), ar(1), al(1), p(1), im, lmt ) 
# end 

# if  ffsl  

# for i ∈ iuw:0 
# al(i) = al(im + i) 
# ar(i) = ar(im + i) 
# a6(i) = a6(im + i) 
# end 

# do i=im + 1, iue 
# al(i) = al(i - im) 
# ar(i) = ar(i - im) 
# a6(i) = a6(i - im) 
# end 

# for i ∈ 1:im 
# iu = c(i) 
# ru = c(i) - iu 
# if c(i) > 0.0 
# itmp = i - iu - 1 
# isave(i) = itmp + 1 
# fx(i) = ru*(ar(itmp)+0.5*ru*(al(itmp)-ar(itmp) + a6(itmp)*(1.0 - r23*ru)) ) 
# else 
# itmp = i - iu 
# isave(i) = itmp - 1 
# fx(i) = ru*(al(itmp)-0.5*ru*(ar(itmp)-al(itmp) + a6(itmp)*(1.0 + r23*ru)) ) 
# end 
# end 

# else 
# al(0) = al(im) 
# ar(0) = ar(im) 
# a6(0) = a6(im) 
# for i ∈ 1:im 
# if c(i) > 0.0 
# fx(i) = ar(i - 1) + 0.5*c(i)*(al(i - 1) - ar(i - 1) + a6(i - 1)*(1.0 - r23*c(i)) ) 
# else 
# fx(i) = al(i) - 0.5*c(i)*(ar(i) - al(i) + a6(i)*(1.0 + r23*c(i))) 
# end 
# fx(i) = mfx(i) * fx(i) 
# end 
# end 
# end 

# function lmppm(dm, a6, ar, al, p, im, lmt) 



# INPUT PARAMETERS: 
#  integer, intent(in):: im # Total longitudes 
#  integer, intent(in):: lmt # LMT = 0: full monotonicity 
# LMT = 1: Improved and simplified full monotonic constraint 
# LMT = 2: positive - definite constraint 
# LMT = 3: Quasi - monotone constraint 
# real, intent(in):: p(im) 
# real, intent(in):: dm(im) 

# real, intent(inout):: a6(im) 
# real, intent(inout):: ar(im) 
# real, intent(inout):: al(im) 

# LOCAL VARIABLES: 
# r12::AbstractFloat 
# parameter ( r12 = 1.0/12.0 ) 

# da1, da2, fmin, a6da::AbstractFloat = 0 
# dr, dl::AbstractFloat = 0 

# i::Integer 

# LMT = 0: full monotonicity 
# LMT = 1: Improved and simplified full monotonic constraint 
# LMT = 2: positive - definite constraint 
# LMT = 3: Quasi - monotone constraint 

# if  lmt == 0  

# Full constraint 
# for i ∈ 1:im 
# if dm(i) == 0.0 
# ar(i) = p(i) 
# al(i) = p(i) 
# a6(i) = 0.0 
# else 
# da1  = ar(i) - al(i) 
# da2  = da1^2 
# a6da = a6(i)*da1 
# if a6da < -da2 
# a6(i) = 3.0*(al(i)-p(i)) 
# ar(i) = al(i) - a6(i) 
# elseif a6da > da2 
# a6(i) = 3.0*(ar(i)-p(i)) 
# al(i) = ar(i) - a6(i) 
# end 
# end 
# end 

# elseif lmt == 1  

# Improved (Lin 200?) full constraint 
# for i ∈ 1:im 
# da1 = dm(i) + dm(i) 
# dl = sign(min(abs(da1), abs(al(i)-p(i))), da1) 
# dr = sign(min(abs(da1), abs(ar(i)-p(i))), da1) 
# ar(i) = p(i) + dr 
# al(i) = p(i) - dl 
# a6(i) = 3.0*(dl - dr) 
# end 

# elseif lmt == 2  
# Positive definite only constraint 
# do 250 i=1, im 
# if abs(ar(i)-al(i)) .ge. -a6(i)) go to 250 
# fmin = p(i) + 0.25*(ar(i)-al(i))^2/a6(i) + a6(i)*r12 
# if fmin.ge.0.0) go to 250 
# if p(i) < ar(i) && p(i) < al(i) 
# ar(i) = p(i) 
# al(i) = p(i) 
# a6(i) = 0.0 
# elseif ar(i) > al(i) 
# a6(i) = 3.0*(al(i)-p(i)) 
# ar(i) = al(i) - a6(i) 
# else 
# a6(i) = 3.0*(ar(i)-p(i)) 
# al(i) = ar(i) - a6(i) 
# end 
# 250   continue 

# elseif lmt == 3 
# Quasi - monotone constraint 
# for i ∈ 1:im 
# da1 = 4.0*dm(i) 
# dl = sign(min(abs(da1), abs(al(i)-p(i))), da1) 
# dr = sign(min(abs(da1), abs(ar(i)-p(i))), da1) 
# ar(i) = p(i) + dr 
# al(i) = p(i) - dl 
# a6(i) = 3.0*(dl - dr) 
# end 
# end 
# end 


# function huynh(im, ar, al, p, d2, d1) 



# INPUT PARAMETERS: 
# im::Integer 
# p::AbstractFloat = zeros(im) 

# OUTPUT PARAMETERS: 
# ar::AbstractFloat = zeros(im) 
# al::AbstractFloat = zeros(im) 
# d2::AbstractFloat = zeros(im) 
# d1::AbstractFloat = zeros(im) 

# LOCAL VARIABLES: 
# i::Integer 
# pmp::AbstractFloat 
# lac::AbstractFloat 
# pmin::AbstractFloat 
# pmax::AbstractFloat 

# Compute d1 and d2 
# d1(1) = p(1) - p(im) 
# for i ∈ 2:im 
# d1(i) = p(i) - p(i - 1) 
# end 

# for i ∈ 1:im - 1 
# d2(i) = d1(i + 1) - d1(i) 
# end 
# d2(im) = d1(1) - d1(im) 

# Constraint for AR 
# i = 1 
# pmp   = p(1) + 2.0 * d1(1) 
# lac   = p(1) + 0.5 * (d1(1)+d2(im)) + d2(im) 
# pmin  = min(p(1), pmp, lac) 
# pmax  = max(p(1), pmp, lac) 
# ar(1) = min(pmax, max(ar(1), pmin)) 

# for i ∈ 2:im 
# pmp   = p(i) + 2.0*d1(i) 
# lac   = p(i) + 0.5*(d1(i)+d2(i - 1)) + d2(i - 1) 
# pmin  = min(p(i), pmp, lac) 
# pmax  = max(p(i), pmp, lac) 
# ar(i) = min(pmax, max(ar(i), pmin)) 
# end 

# Constraint for AL 
# for i ∈ 1:im - 1 
# pmp   = p(i) - 2.0*d1(i + 1) 
# lac   = p(i) + 0.5*(d2(i + 1)-d1(i + 1)) + d2(i + 1) 
# pmin  = min(p(i), pmp, lac) 
# pmax  = max(p(i), pmp, lac) 
# al(i) = min(pmax, max(al(i), pmin)) 
# end 

# i=im 
# i = im 
# pmp    = p(im) - 2.0*d1(1) 
# lac    = p(im) + 0.5*(d2(1)-d1(1)) + d2(1) 
# pmin   = min(p(im), pmp, lac) 
# pmax   = max(p(im), pmp, lac) 
# al(im) = min(pmax, max(al(im), pmin)) 

# compute A6 (d2) 
# for i ∈ 1:im 
# d2(i) = 3.0*(p(i)+p(i)  - (al(i)+ar(i))) 
# end 
# end 


# function ytp(im, jm, fy, q, c, yfx, ng, mg, jord, iv, jfirst, jlast) 



# INPUT PARAMETERS: 
#  im, jm::Integer = 0 # Dimensions 
#  jfirst, jlast::Integer = 0 # Latitude strip 
#  ng::Integer # Max. NS dependencies 
#  mg::Integer # 
#  jord::Integer # order of subgrid dist 
#  iv::Integer # Scalar=0, Vector=1 
#  q::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # advected scalar N*jord S*jord 
#  c::AbstractFloat = zeros(im, jfirst:jlast + mg) # Courant   N (like FY) 
#  yfx::AbstractFloat = zeros(im, jfirst:jlast + mg) # Backgrond mass flux 

# OUTPUT PARAMETERS: 
#  fy::AbstractFloat = zeros(im, jfirst:jlast + mg) # Flux      N (see tp2c) 

# LOCAL VARIABLES: 
# i, j, jt::Integer = 0 
# js2g0, jn1g1::Integer = 0 

# work arrays (should pass in eventually for performance enhancement): 
# dm::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) 

#     ar::AbstractFloat = zeros(im, jfirst - 1:jlast + 1) # AR needs to be ghosted on NS 
#     al::AbstractFloat = zeros(im, jfirst - 1:jlast + 2) # AL needs to be ghosted on N2S 
#     a6::AbstractFloat = zeros(im, jfirst - 1:jlast + 1) # A6 needs to be ghosted on NS 


#    js2g0  = max(2, jfirst) # No ghosting 
#    jn1g1  = min(jm, jlast + 1) # Ghost N*1 

# if jord == 1 
# for j ∈ js2g0:jn1g1 
# for i ∈ 1:im 
# jt = float(j) - c(i, j) 
# fy(i, j) = q(i, jt) 
# end 
# end 
# else 

# 
# YMIST requires q on NS;  Only # call to YMIST here 
# 
# call ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast) 

# if  abs(jord) .ge. 3  

# call fyppm(c, q, dm, fy, im, jm, ng, mg, jord, iv, jfirst, jlast) 

# else 
# 
# JORD can either have the value 2 or -2 at this point 
# 
# for j ∈ js2g0:jn1g1 
# for i ∈ 1:im 
# jt = float(j) - c(i, j) 
# fy(i, j) = q(i, jt) + (sign(1.0, c(i, j))-c(i, j))*dm(i, jt) 
# end 
# end 
# end 
# end 

# for j ∈ js2g0:jn1g1 
# for i ∈ 1:im 
# fy(i, j) = fy(i, j)*yfx(i, j) 
# end 
# end 
# end 

# function ymist(im, jm, q, dm, ng, jord, iv, jfirst, jlast) 



# INPUT PARAMETERS: 
#  im, jm::Integer = 0 # Dimensions 
#  jfirst, jlast::Integer = 0 # Latitude strip 
#  ng::Integer # NS dependencies 
#  jord::Integer # order of subgrid distribution 
#  iv::Integer # Scalar (==0) Vector (==1) 
#  q::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # transported scalar  N*ng S*ng 

# OUTPUT PARAMETERS: 
#  dm::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # Slope only N*(ng - 1) S*(ng - 1) used 

# Local variables 

# i, j, jm1, im2, js2gng1, jn2gng1::Integer = 0 
# qmax, qmin, tmp::AbstractFloat = 0 

#     js2gng1 = max(2,   jfirst - ng + 1) # Number needed on S 
#     jn2gng1 = min(jm - 1, jlast + ng - 1) # Number needed on N 

# jm1 = jm - 1 
# im2 = im / 2 

# for j ∈ js2gng1:jn2gng1 
# for i ∈ 1:im 
# dm(i, j) = 0.25*(q(i, j + 1) - q(i, j - 1)) 
# end 
# end 

# if  iv == 0  

# if  jfirst == 1  
# S pole 
# for i ∈ 1:im2 
# tmp = 0.25*(q(i, 2)-q(i + im2, 2)) 
# qmax = max(q(i, 2), q(i, 1), q(i + im2, 2)) - q(i, 1) 
# qmin = q(i, 1) - min(q(i, 2), q(i, 1), q(i + im2, 2)) 
# dm(i, 1) = sign(min(abs(tmp), qmax, qmin), tmp) 
# end 

# do i=im2 + 1, im 
# dm(i, 1) =  - dm(i - im2, 1) 
# end 
# end 

# if  jlast == jm  
# N pole 
# for i ∈ 1:im2 
# tmp = 0.25*(q(i + im2, jm1)-q(i, jm1)) 
# qmax = max(q(i + im2, jm1), q(i, jm), q(i, jm1)) - q(i, jm) 
# qmin = q(i, jm) - min(q(i + im2, jm1), q(i, jm), q(i, jm1)) 
# dm(i, jm) = sign(min(abs(tmp), qmax, qmin), tmp) 
# end 

# do i=im2 + 1, im 
# dm(i, jm) =  - dm(i - im2, jm) 
# end 
# end 

# else 

# if  jfirst == 1  
# South 
# for i ∈ 1:im2 
# tmp  = 0.25*(q(i, 2)+q(i + im2, 2)) 
# qmax = max(q(i, 2), q(i, 1), -q(i + im2, 2)) - q(i, 1) 
# qmin = q(i, 1) - min(q(i, 2), q(i, 1),-q(i + im2, 2)) 
# dm(i, 1) = sign(min(abs(tmp), qmax, qmin), tmp) 
# end 

# do i=im2 + 1, im 
# dm(i, 1) = dm(i - im2, 1) 
# end 
# end 

# if  jlast == jm  
# North 
# for i ∈ 1:im2 
# tmp  = -0.25*(q(i + im2, jm1)+q(i, jm1)) 
# qmax = max(-q(i + im2, jm1), q(i, jm), q(i, jm1)) - q(i, jm) 
# qmin = q(i, jm) - min(-q(i + im2, jm1), q(i, jm), q(i, jm1)) 
# dm(i, jm) = sign(min(abs(tmp), qmax, qmin), tmp) 
# end 

# do i=im2 + 1, im 
# dm(i, jm) = dm(i - im2, jm) 
# end 
# end 

# end 

# if  jord > 0  
# 
# Applies monotonic slope constraint (off if jord less than zero) 
# 
# for j ∈ js2gng1:jn2gng1 
# for i ∈ 1:im 
# qmax = max(q(i, j - 1), q(i, j), q(i, j + 1)) - q(i, j) 
# qmin = q(i, j) - min(q(i, j - 1), q(i, j), q(i, j + 1)) 
# dm(i, j) = sign(min(abs(dm(i, j)), qmin, qmax), dm(i, j)) 
# end 
# end 
# end 
# end 

# function fyppm(c,  q,  dm, flux, im, jm, ng, mg, jord, iv, jfirst, jlast) 



# INPUT PARAMETERS: 
#  im, jm::Integer = 0 # Dimensions 
#  jfirst, jlast::Integer = 0 # Latitude strip 
#  ng::Integer # Max. NS dependencies 
#  mg::Integer # 
#  jord::Integer # Approximation order 
#  iv::Integer # Scalar=0, Vector=1 
#  q::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # mean value needed only N*2 S*2 
#  dm::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # Slope     needed only N*2 S*2 
#  c::AbstractFloat = zeros(im, jfirst:jlast + mg) # Courant   N (like FLUX) 

# INPUT/OUTPUT PARAMETERS: 
#  ar::AbstractFloat = zeros(im, jfirst - 1:jlast + 1) # AR needs to be ghosted on NS 
#  al::AbstractFloat = zeros(im, jfirst - 1:jlast + 2) # AL needs to be ghosted on N2S 
#  a6::AbstractFloat = zeros(im, jfirst - 1:jlast + 1) # A6 needs to be ghosted on NS 

# OUTPUT PARAMETERS: 
#  flux::AbstractFloat = zeros(im, jfirst:jlast + mg) # Flux      N (see tp2c) 

# Local 
# r3, r23::AbstractFloat = 0 
# parameter ( r3 = 1.0/3.0, r23 = 2.0/3.0 ) 
# i, j, imh, jm1, lmt::Integer = 0 
# js1g1, js2g0, js2g1, jn1g2, jn1g1, jn2g1::Integer = 0 

# --------------------------------------------------------------------- 
# Initialize local variables (bmy, 7/10/17) 
# ar = 0.0 
# al = 0.0 
# a6 = 0.0 
# --------------------------------------------------------------------- 

# imh = im / 2 
# jm1 = jm - 1 

#       js1g1  = max(1, jfirst - 1) # Ghost S*1 
#       js2g0  = max(2, jfirst) # No ghosting 
#       js2g1  = max(2, jfirst - 1) # Ghost S*1 
#       jn1g1  = min(jm, jlast + 1) # Ghost N*1 
#       jn1g2  = min(jm, jlast + 2) # Ghost N*2 
#       jn2g1  = min(jm - 1, jlast + 1) # Ghost N*1 

#       for j ∈ js2g1:jn1g2 # AL needed N2S 
#         for i ∈ 1:im # P, dm ghosted N2S2 (at least) 
# al(i, j) = 0.5*(q(i, j - 1)+q(i, j)) + r3*(dm(i, j - 1) - dm(i, j)) 
# end 
# end 

#       for j ∈ js1g1:jn2g1 # AR needed NS 
# for i ∈ 1:im 
#           ar(i, j) = al(i, j + 1) # AL ghosted N2S 
# end 
# end 

# Poles: 

# if  iv == 0  

# if  jfirst == 1  
# for i ∈ 1:imh 
# al(i,    1) = al(i + imh, 2) 
# al(i + imh, 1) = al(i,    2) 
# end 
# end 

# if  jlast == jm  
# for i ∈ 1:imh 
# ar(i,    jm) = ar(i + imh, jm1) 
# ar(i + imh, jm) = ar(i,    jm1) 
# end 
# end 

# else 

# if  jfirst == 1  
# for i ∈ 1:imh 
# al(i,    1) = -al(i + imh, 2) 
# al(i + imh, 1) = -al(i,    2) 
# end 
# end 

# if  jlast == jm  
# for i ∈ 1:imh 
# ar(i,    jm) = -ar(i + imh, jm1) 
# ar(i + imh, jm) = -ar(i,    jm1) 
# end 
# end 

# end 

# if  jord == 3 .or. jord == 5  
#       for j ∈ js1g1:jn1g1 # A6 needed NS 
# for i ∈ 1:im 
# a6(i, j) = 3.0*(q(i, j)+q(i, j) - (al(i, j)+ar(i, j))) 
# end 
# end 
# end 

# lmt = jord - 3 

# call lmppm(dm(1, js1g1), a6(1, js1g1), ar(1, js1g1), al(1, js1g1),  q(1, js1g1), im*(jn1g1 - js1g1 + 1), lmt) 

#       for j ∈ js2g0:jn1g1 # flux needed N 
# for i ∈ 1:im 
# if c(i, j) > 0.0 
# flux(i, j) = ar(i, j - 1) + 0.5*c(i, j)*(al(i, j - 1) - ar(i, j - 1) + a6(i, j - 1)*(1.0 - r23*c(i, j)) ) 
# else 
# flux(i, j) = al(i, j) - 0.5*c(i, j)*(ar(i, j) - al(i, j) + a6(i, j)*(1.0 + r23*c(i, j))) 
# end 
# end 
# end 
# end 

# function xpavg(p, im) 



# INPUT PARAMETERS: 
# im::Integer 

# INPUT/OUTPUT PARAMETERS: 
# p::AbstractFloat = zeros(im) 

# i::Integer 
# sum1::AbstractFloat 

# p(1:im) = sum(p(1:im))/im 

# sum1 = 0.0 
# for i ∈ 1:im 
# sum1 = sum1 + p(i) 
# end 
# sum1 = sum1 / im 

# for i ∈ 1:im 
# p(i) = sum1 
# end 
# end 

# function qmap(pe,  q, im, jm, km, nx, jfirst, jlast, ng, nq, ps,  ak, bk, kord, iv) 



# INPUT 
#    im, jm, km::Integer = 0 # x, y, z dimensions 
#    nq::Integer # number of tracers 
#    nx::Integer # number of SMP "decomposition" in x 
#    iv::Integer # monotonicity at top and bottom 
# iv=0 : weak constraint 
# iv=1 : strong constraint 
# iv=-1: for vector 
#    jfirst, jlast::Integer = 0 # starting & ending latitude index 
#    ng::Integer # width of ghost regions 
# real, intent(in):: ak(km + 1) 
# real, intent(in):: bk(km + 1) 
# real, intent(in):: pe(im, km + 1, jfirst:jlast) 

# INPUT/OUTPUT 
#    q::AbstractFloat = zeros(im, jfirst - ng:jlast + ng, km, nq) # tracers including specific humidity 
#    ps::AbstractFloat = zeros(im, jfirst:jlast) # surface pressure 

# Local arrays: 
# pe2::AbstractFloat = zeros(im, km + 1) 

# temp::AbstractFloat 
# i, j, k, iq::Integer = 0 
# ixj, jp, it, i1, i2::Integer = 0 
# kord::Integer 


# it = im / nx 
# jp = nx * ( jlast - jfirst + 1 ) 



#$omp parallel do  #$omp shared(im, km, jfirst, jlast, ng, iv, kord) # $omp private(i, j, k, iq, i1, i2, ixj, pe2) 

# do 2000 j=jfirst, jlast 
# do 2000 ixj=1, jp 

# j  = jfirst + (ixj - 1) / nx 
# i1 = 1 + it * mod(ixj - 1, nx) 
# i2 = i1 + it - 1 


# k=1 
# for i ∈ i1:i2 
# pe2(i, 1) = ak(1) 
# end 

# for k ∈ 2:km 
# for i ∈ i1:i2 
# pe2(i, k) = ak(k) + bk(k)*ps(i, j) 
# end 
# end 

# k=km + 1 
# for i ∈ i1:i2 
# pe2(i, km + 1) = ps(i, j) 
# end 

# temp = sum(q) 
# for iq ∈ 1:nq 
# call map1_ppm ( km,   pe(1, 1, j),   q(1, jfirst - ng, 1, iq), km,   pe2,         q(1, jfirst - ng, 1, iq), im, i1, i2,  j, jfirst, jlast, ng, iv, kord) 
# end 
# 2000  continue 

# end 

# function map1_ppm( km,   pe1,   q1, kn,   pe2,   q2, im, i1, i2, j, jfirst, jlast, ng, iv, kord) 



# INPUT PARAMETERS: 
#  i1::Integer # Starting longitude 
#  i2::Integer # Finishing longitude 
#  im::Integer # E - W dimension 
#  iv::Integer # Mode: 0 ==  constituents  1 == ??? 
#  kord::Integer # Method order 
#  j::Integer # Current latitude 
#  jfirst::Integer # Starting latitude 
#  jlast::Integer # Finishing latitude 
#  ng::Integer # Width of ghost regions 
#  km::Integer # Original vertical dimension 
#  kn::Integer # Target vertical dimension 

#  pe1::AbstractFloat = zeros(im, km + 1) # pressure at layer edges 
# (from model top to bottom surface) 
# in the original vertical coordinate 
#  pe2::AbstractFloat = zeros(im, kn + 1) # pressure at layer edges 
# (from model top to bottom surface) 
# in the new vertical coordinate 
#  q1::AbstractFloat = zeros(im, jfirst - ng:jlast + ng, km) # Field input 

# OUTPUT PARAMETERS: 
#  q2::AbstractFloat = zeros(im, jfirst - ng:jlast + ng, kn) # Field output 

# LOCAL VARIABLES: 

# dp1::AbstractFloat = zeros(i1:i2, km) 
# q4::AbstractFloat = zeros(4, i1:i2, km) 
# i, k, l, ll, k0::Integer = 0 
# real  pl, pr, qsum, delp, esl 
# real       r3, r23 
# temp::AbstractFloat 
# parameter (r3 = 1.0/3.0, r23 = 2.0/3.0) 

# Initialize local arrays (bmy, 7/10/17) 
# dp1 = 0.0 
# q4  = 0.0 

# for k ∈ 1:km 
# for i ∈ i1:i2 
# dp1(i, k) = pe1(i, k + 1) - pe1(i, k) 
# q4(1, i, k) = q1(i, j, k) 
# end 
# end 

# temp = sum(q4) 
# Compute vertical subgrid distribution 
# call ppm2m( q4, dp1, km, i1, i2, iv, kord ) 

# temp = sum(q2) 
# Mapping 
# do 1000 i=i1, i2 
# k0 = 1 
# do 555 k=1, kn 
# do 100 l=k0, km 
# locate the top edge: pe2(i, k) 
# if pe2(i, k) .ge. pe1(i, l) && pe2(i, k) .le. pe1(i, l + 1) 
# pl = (pe2(i, k)-pe1(i, l)) / dp1(i, l) 
# if pe2(i, k + 1) .le. pe1(i, l + 1) 
# entire new grid is within the original grid 
# pr = (pe2(i, k + 1)-pe1(i, l)) / dp1(i, l) 
# q2(i, j, k) = q4(2, i, l) + 0.5*(q4(4, i, l)+q4(3, i, l)-q4(2, i, l)) *(pr + pl)-q4(4, i, l)*r3*(pr*(pr + pl)+pl^2) 
# k0 = l 
# goto 555 
# else 
# Fractional area... 
# qsum = (pe1(i, l + 1)-pe2(i, k))*(q4(2, i, l)+0.5*(q4(4, i, l)+ q4(3, i, l)-q4(2, i, l))*(1.0 + pl)-q4(4, i, l)* (r3*(1.0 + pl*(1.0 + pl)))) 
# do ll=l + 1, km 
# locate the bottom edge: pe2(i, k + 1) 
# if pe2(i, k + 1) > pe1(i, ll + 1)  
# Whole layer.. 
# qsum = qsum + dp1(i, ll)*q4(1, i, ll) 
# else 
# delp = pe2(i, k + 1)-pe1(i, ll) 
# esl = delp / dp1(i, ll) 
# qsum = qsum + delp*(q4(2, i, ll)+0.5*esl* (q4(3, i, ll)-q4(2, i, ll)+q4(4, i, ll)*(1.0 - r23*esl))) 
# k0 = ll 
# goto 123 
# end 
# end 
# goto 123 
# end 
# end 
# 100   continue 
# 123   q2(i, j, k) = qsum / ( pe2(i, k + 1) - pe2(i, k) ) 
# 555   continue 
# 1000  continue 

# end 

function ppm2m!(
    a4,
    delp,
    km,
    i1,
    i2,
    iv,
    kord
)::Nothing
    # ! INPUT PARAMETERS: 
    #  integer, intent(in):: iv      ! iv =-1: winds 
    #                                ! iv = 0: positive definite scalars 
    #                                ! iv = 1: others 
    #  integer, intent(in):: i1      ! Starting longitude 
    #  integer, intent(in):: i2      ! Finishing longitude 
    #  integer, intent(in):: km      ! vertical dimension 
    #  integer, intent(in):: kord    ! Order (or more accurately method no.): 
    #                                ! 
    #  real, intent(in):: delp(i1:i2,km)     ! layer pressure thickness 

    #  real, intent(inout):: a4(4,i1:i2,km)  ! Interpolated values 

    # ! local arrays. 
    #  real dc(i1:i2,km) 
    #  real h2(i1:i2,km) 
    #  real delq(i1:i2,km) 
    #  real df2(i1:i2,km) 
    #  real d4(i1:i2,km) 

    #  real fac 
    #  real a1, a2, c1, c2, c3, d1, d2 
    #  real qmax, qmin, cmax, cmin 
    #  real qm, dq, tmp 
    #  real qmp, pmp 
    #  real lac 
    #  integer lmt 
    #  integer i, k, km1 
    #  integer it 

    #    km1 = km - 1 
    #    it = i2 - i1 + 1 

    #     do k=2,km 
    #        do i=i1,i2 
    #           delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1) 
    #             d4(i,k  ) = delp(i,k-1) + delp(i,k) 
    #        enddo 
    #     enddo 

    #     do k=2,km1 
    #        do i=i1,i2 
    #           c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1) 
    #           c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k) 
    #           tmp = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /       & 
    #                                   (d4(i,k)+delp(i,k+1)) 
    #           qmax = max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) - a4(1,i,k) 
    #           qmin = a4(1,i,k) - min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1)) 
    #            dc(i,k) = sign(min(abs(tmp),qmax,qmin), tmp) 
    #           df2(i,k) = tmp 
    #        enddo 
    #     enddo 

    # !------------------------------------------------------------ 
    # ! 4th order interpolation of the provisional cell edge value 
    # !------------------------------------------------------------ 

    #     do k=3,km1 
    #       do i=i1,i2 
    #       c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k) 
    #       a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1)) 
    #       a2 = d4(i,k+1) / (d4(i,k) + delp(i,k)) 
    #       a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *      & 
    #                 ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -            & 
    #                                 delp(i,k-1)*a1*dc(i,k  ) ) 
    #       enddo 
    #     enddo 

    #     if(kord>3) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4) 

    # ! Area preserving cubic with 2nd deriv. = 0 at the boundaries 
    # ! Top 
    #     do i=i1,i2 
    #       d1 = delp(i,1) 
    #       d2 = delp(i,2) 
    #       qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2) 
    #       dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2) 
    #       c1 = (a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) ) 
    #       c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2) 
    #       a4(2,i,2) = qm - c1*d1*d2*(d2+3.*d1) 
    #       a4(2,i,1) = d1*(8.*c1*d1**2-c3) + a4(2,i,2) 
    #       dc(i,1) =  a4(1,i,1) - a4(2,i,1) 
    # ! No over- and undershoot condition 
    #       cmax = max(a4(1,i,1), a4(1,i,2)) 
    #       cmin = min(a4(1,i,1), a4(1,i,2)) 
    #       a4(2,i,2) = max(cmin,a4(2,i,2)) 
    #       a4(2,i,2) = min(cmax,a4(2,i,2)) 
    #     enddo 

    #     if( iv == 0 ) then 
    #         do i=i1,i2 
    #             a4(2,i,1) = max(0.,a4(2,i,1)) 
    #         enddo 
    #     elseif ( iv ==  1 ) then 
    # ! Monotone tracers: 
    #         do i=i1,i2 
    #            dc(i,1) = 0. 
    #            a4(2,i,1) = a4(1,i,1) 
    #            a4(2,i,2) = a4(1,i,1) 
    #         enddo 
    #     elseif ( iv == -1 ) then 
    # ! Winds: 
    #         do i=i1,i2 
    #             if( a4(1,i,1)*a4(2,i,1) <=  0. ) then 
    #                 a4(2,i,1) = 0. 
    #             else 
    #                 a4(2,i,1) = sign(min(abs(a4(1,i,1)),      & 
    #                                      abs(a4(2,i,1))),     & 
    #                                          a4(1,i,1)  ) 
    #             endif 
    #         enddo 
    #     endif 

    # ! Bottom 
    # ! Area preserving cubic with 2nd deriv. = 0 at the surface 
    #     do i=i1,i2 
    #        d1 = delp(i,km) 
    #        d2 = delp(i,km1) 
    #        qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2) 
    #        dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2) 
    #        c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1))) 
    #        c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1**2) 
    #        a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1) 
    #        a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km) 
    #        dc(i,km) = a4(3,i,km) -  a4(1,i,km) 
    # ! No over- and under-shoot condition 
    #        cmax = max(a4(1,i,km), a4(1,i,km1)) 
    #        cmin = min(a4(1,i,km), a4(1,i,km1)) 
    #        a4(2,i,km) = max(cmin,a4(2,i,km)) 
    #        a4(2,i,km) = min(cmax,a4(2,i,km)) 
    #     enddo 

    # ! Enforce constraint at the surface 

    #     if ( iv == 0 ) then 
    # ! Positive definite scalars: 
    #          do i=i1,i2 
    #             a4(3,i,km) = max(0., a4(3,i,km)) 
    #          enddo 
    #     elseif ( iv ==  1 ) then 
    # ! Monotone tracers: 
    #          do i=i1,i2 
    #             dc(i,km) = 0. 
    #             a4(2,i,km) = a4(1,i,km) 
    #             a4(3,i,km) = a4(1,i,km) 
    #          enddo 
    #     elseif ( iv == -1 ) then 
    # ! Winds: 
    #          do i=i1,i2 
    #             if( a4(1,i,km)*a4(3,i,km) <=  0. ) then 
    #                 a4(3,i,km) = 0. 
    #             else 
    #                 a4(3,i,km) = sign( min(abs(a4(1,i,km)),      & 
    #                                        abs(a4(3,i,km))),     & 
    #                                            a4(1,i,km)  ) 
    #             endif 
    #          enddo 
    #     endif 

    #     do k=1,km1 
    #        do i=i1,i2 
    #           a4(3,i,k) = a4(2,i,k+1) 
    #        enddo 
    #     enddo 

    # ! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 ) 

    # ! Top 2 and bottom 2 layers always use monotonic mapping 
    #     do k=1,2 
    #        do i=i1,i2 
    #           a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))) 
    #        enddo 
    #        call kmppm(dc(i1,k), a4(1,i1,k), it, 0) 
    #     enddo 

    if kord >= 7
        # Huynh's 2nd constraint 

        for k ∈ 2:km1
            for i ∈ i1:i2
                # Method#1 
                # h2(i,k) = delq(i,k) - delq(i,k-1) 

                # Method#2 
                # h2(i, k) = 
                #     2.0 * (dc(i, k + 1) / delp(i, k + 1) - dc(i, k - 1) / delp(i, k - 1)) / 
                #     (delp(i, k) + 0.5 * (delp(i, k - 1) + delp(i, k + 1))) * 
                #     delp(i, k)^2 

                # Method#3 
                h2(i, k) = dc(i, k + 1) - dc(i, k - 1)
            end
        end

        if kord == 7
            fac = 1.5 # original quasi-monotone 
        else
            fac = 0.125 # full monotone 
        end

        for k ∈ 3:km-2
            for i ∈ i1:i2
                # Right edges 
                # qmp   = a4(1,i,k) + 2.0*delq(i,k-1) 
                # lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1) 
                pmp = 2.0 * dc(i, k)
                qmp = a4(1, i, k) + pmp
                lac = a4(1, i, k) + fac * h2(i, k - 1) + dc(i, k)
                qmin = min(a4(1, i, k), qmp, lac)
                qmax = max(a4(1, i, k), qmp, lac)
                a4(3, i, k) = min(max(a4(3, i, k), qmin), qmax)
                # Left  edges 
                # qmp = a4(1, i, k) - 2.0 * delq(i, k) 
                # lac = a4(1, i, k) + fac * h2(i, k + 1) - 0.5 * delq(i, k) 
                qmp = a4(1, i, k) - pmp
                lac = a4(1, i, k) + fac * h2(i, k + 1) - dc(i, k)
                qmin = min(a4(1, i, k), qmp, lac)
                qmax = max(a4(1, i, k), qmp, lac)
                a4(2, i, k) = min(max(a4(2, i, k), qmin), qmax)
                # Recompute A6 
                a4(4, i, k) = 3.0 * (2.0 * a4(1, i, k) - (a4(2, i, k) + a4(3, i, k)))
            end
            # Additional constraint to prevent negatives when kord=7 
            if iv != -1 && kord == 7
                # call kmppm(dc(i1,k), a4(1,i1,k), it, 2) 
            end
        end
    else
        lmt = kord - 3
        lmt = max(0, lmt)
        #          if (iv == 0) lmt = min(2, lmt) 

        for k ∈ 3:km-2
            if kord != 4
                for i ∈ i1:i2
                    a4(4, i, k) = 3.0 * (2.0 * a4(1, i, k) - (a4(2, i, k) + a4(3, i, k)))
                end
            end
            # call kmppm(dc(i1,k), a4(1,i1,k), it, lmt) 
        end
    end

    for k ∈ km1:km
        for i ∈ i1:i2
            a4(4, i, k) = 3.0 * (2.0 * a4(1, i, k) - (a4(2, i, k) + a4(3, i, k)))
        end
        #    call kmppm(dc(i1,k), a4(1,i1,k), it, 0) 
    end
end

"""

## Arguments
| Name     | Intent  | Dimensions           | Description
|:---------|:--------|:---------------------|:------------
"""
function steepz!(
    i1,
    i2,
    km,
    a4,
    df2,
    dm,
    dq,
    dp,
    d4
)::Nothing
    # INPUT PARAMETERS: 
    #  km::Integer # Total levels 
    #  i1::Integer # Starting longitude 
    #  i2::Integer # Finishing longitude 
    #  dp::AbstractFloat = zeros(i1:i2, km) # grid size 
    #  dq::AbstractFloat = zeros(i1:i2, km) # backward diff of q 
    #  d4::AbstractFloat = zeros(i1:i2, km) # backward sum:  dp(k)+ dp(k - 1) 
    #  df2::AbstractFloat = zeros(i1:i2, km) # first guess mismatch 
    #  dm::AbstractFloat = zeros(i1:i2, km) # monotonic mismatch 

    # INPUT/OUTPUT PARAMETERS: 
    #  a4::AbstractFloat = zeros(4, i1:i2, km) # first guess/steepened 

    # LOCAL VARIABLES: 
    # i, k::Integer = 0 
    # alfa::AbstractFloat = zeros(i1:i2, km) 
    # f::AbstractFloat = zeros(i1:i2, km) 
    # rat::AbstractFloat = zeros(i1:i2, km) 
    # dg2::AbstractFloat 

    # Compute ratio of dq/dp 
    for k ∈ 2:km
        for i ∈ i1:i2
            rat[i, k] = dq[i, k-1] / d4[i, k]
        end
    end

    # Compute F 
    for k ∈ 2:km-1
        for i ∈ i1:i2
            f[i, k] =
                (rat[i, k+1] - rat[i, k]) /
                (dp[i, k-1] + dp[i, k] + dp[i, k+1])
        end
    end

    for k ∈ 3:km-2
        for i ∈ i1:i2
            if f[i, k+1] * f[i, k-1] < 0.0 && df2[i, k] != 0.0
                dg2 =
                    (f[i, k+1] - f[i, k-1]) *
                    ((dp[i, k+1] - dp[i, k-1])^2 + d4[i, k] * d4[i, k+1])

                alfa[i, k] = max(0.0, min(0.5, -0.1875 * dg2 / df2[i, k]))
            else
                alfa[i, k] = 0.0
            end
        end
    end

    for k ∈ 4:km-2
        for i ∈ i1:i2
            a4[2, i, k] =
                (1.0 - alfa[i, k-1] - alfa[i, k]) * a4[2, i, k] +
                alfa[i, k-1] * (a4[1, i, k] - dm[i, k]) +
                alfa[i, k] * (a4[1, i, k-1] + dm[i, k-1])
        end
    end

end

"""

## Arguments
| Name     | Intent  | Dimensions           | Description
|:---------|:--------|:---------------------|:------------
"""
function kmppm!(
    dm,
    a4,
    itot,
    lmt
)::Nothing
    # INPUT PARAMETERS: 
    #       dm::AbstractFloat = zeros(*) # ?????? 
    #       itot::Integer # Total Longitudes 
    #       lmt::Integer # 0: Standard PPM constraint 
    # 1: Improved full monotonicity constraint (Lin) 
    # 2: Positive definite constraint 
    # 3: do nothing (return immediately) 

    # INPUT/OUTPUT PARAMETERS: 
    # real   a4(4,*) 
    # AA <-- a4(1, i) 
    # AL <-- a4(2, i) 
    # AR <-- a4(3, i) 
    # A6 <-- a4(4, i) 

    # LOCAL VARIABLES: 
    # r12::AbstractFloat 
    # parameter (r12 = 1.0/12.0) 

    # qmp::AbstractFloat 

    # i::Integer 
    # da1, da2, a6da::AbstractFloat = 0 
    # fmin::AbstractFloat 

    if lmt == 3
        return
    end

    if lmt == 0

        # Standard PPM constraint 
        for i ∈ 1:itot
            if dm[i] == 0.0
                a4[2, i] = a4[1, i]
                a4[3, i] = a4[1, i]
                a4[4, i] = 0.0
            else
                da1 = a4[3, i] - a4[2, i]
                da2 = da1^2
                a6da = a4[4, i] * da1

                if a6da < -da2
                    a4[4, i] = 3.0 * (a4[2, i] - a4[1, i])
                    a4[3, i] = a4[2, i] - a4[4, i]
                elseif a6da > da2
                    a4[4, i] = 3.0 * (a4[3, i] - a4[1, i])
                    a4[2, i] = a4[3, i] - a4[4, i]
                end
            end
        end

    elseif lmt == 1

        # Improved full monotonicity constraint (Lin) 
        # Note: no need to provide first guess of A6 <-- a4(4, i) 
        for i ∈ 1:itot
            qmp = 2.0 * dm(i)
            a4[2, i] = a4[1, i] - sign(qmp) * min(abs(qmp), abs(a4[2, i] - a4[1, i]))
            a4[3, i] = a4[1, i] + sign(qmp) * min(abs(qmp), abs(a4[3, i] - a4[1, i]))
            a4[4, i] = 3.0 * (2.0 * a4[1, i] - (a4[2, i] + a4[3, i]))
        end

    elseif lmt == 2

        # Positive definite constraint 
        for i ∈ 1:itot
            if abs(a4[3, i] - a4[2, i]) < -a4[4, i]
                fmin = a4[1, i] + 0.25 * (a4[3, i] - a4[2, i])^2 / a4[4, i] + a4[4, i] * r12
                if fmin < 0.0
                    if a4[1, i] < a4[3, i] && a4[1, i] < a4[2, i]
                        a4[3, i] = a4[1, i]
                        a4[2, i] = a4[1, i]
                        a4[4, i] = 0.0
                    elseif a4[3, i] > a4[2, i]
                        a4[4, i] = 3.0 * (a4[2, i] - a4[1, i])
                        a4[3, i] = a4[2, i] - a4[4, i]
                    else
                        a4[4, i] = 3.0 * (a4[3, i] - a4[1, i])
                        a4[2, i] = a4[3, i] - a4[4, i]
                    end
                end
            end
        end

    end

end

"""

## Arguments
| Name     | Intent  | Dimensions           | Description
|:---------|:--------|:---------------------|:------------
"""
function fct_x!(
    q,
    im,
    jm,
    jfirst,
    jlast,
    ng,
    ipx
)::Nothing
    # INPUT PARAMETERS: 
    #  im::Integer # Longitudes 
    #  jm::Integer # Total latitudes 
    #  jfirst::Integer # Starting latitude 
    #  jlast::Integer # Finishing latitude 
    # ng::Integer 

    const tiny = 1.0e-40 # A small number to pump up value 

    # INPUT/OUTPUT PARAMETERS: 
    #  q::AbstractFloat = zeros(im, jfirst - ng:jlast + ng) # Field to adjust 

    # OUTPUT PARAMETERS: 
    #  ipx::Integer # Flag:  0 if Q not change, 1 if changed 

    # LOCAL VARIABLES: 
    # d0, d1, d2::AbstractFloat = 0 
    # qtmp::AbstractFloat = zeros(jfirst:jlast, im) 

    # i, j, jm1, ip2::Integer = 0 
    # j1, j2::Integer = 0 

    # j1 = max( jfirst,   2 ) 
    # j2 = min( jlast, jm - 1 ) 
    # jm1 = jm - 1 
    # ipx = 0 

    # Copy & swap direction for vectorization. 
    for i ∈ 1:im
        for j ∈ j1:j2
            qtmp[j, i] = q[i, j]
        end
    end

    for i ∈ 2:im-1
        for j ∈ j1:j2
            if qtmp[j, i] < 0.0
                ipx = 1
                # west 
                d0 = max(0.0, qtmp[j, i-1])
                d1 = min(-qtmp[j, i], d0)
                qtmp[j, i-1] = qtmp[j, i-1] - d1
                qtmp[j, i] = qtmp[j, i] + d1
                # east 
                d0 = max(0.0, qtmp[j, i+1])
                d2 = min(-qtmp[j, i], d0)
                qtmp[j, i+1] = qtmp[j, i+1] - d2
                qtmp[j, i] = qtmp[j, i] + d2 + tiny
            end
        end
    end

    i = 1

    for j ∈ j1:j2
        if qtmp[j, i] < 0.0
            ipx = 1
            # west 
            d0 = max(0.0, qtmp[j, im])
            d1 = min(-qtmp[j, i], d0)
            qtmp[j, im] = qtmp[j, im] - d1
            qtmp[j, i] = qtmp[j, i] + d1
            # east 
            d0 = max(0.0, qtmp[j, i+1])
            d2 = min(-qtmp[j, i], d0)
            qtmp[j, i+1] = qtmp[j, i+1] - d2
            qtmp[j, i] = qtmp[j, i] + d2 + tiny
        end
    end

    i = im

    for j ∈ j1:j2
        if qtmp[j, i] < 0.0
            ipx = 1
            # west 
            d0 = max(0.0, qtmp[j, i-1])
            d1 = min(-qtmp[j, i], d0)
            qtmp[j, i-1] = qtmp[j, i-1] - d1
            qtmp[j, i] = qtmp[j, i] + d1
            # east 
            d0 = max(0.0, qtmp[j, 1])
            d2 = min(-qtmp[j, i], d0)
            qtmp[j, 1] = qtmp[j, 1] - d2

            qtmp[j, i] = qtmp[j, i] + d2 + tiny
        end
    end


    if ipx != 0
        # Final pass 
        for i ∈ 1:im-1
            for j ∈ j1:j2
                if qtmp[j, i] < 0.0
                    # Take mass from east (essentially adjusting fx(i + 1, j)) 
                    qtmp[j, i+1] = qtmp[j, i+1] + qtmp[j, i]
                    qtmp[j, i] = 0.0
                end
            end
        end

        # Final sweep 
        for i ∈ im:-1:2
            for j ∈ j1:j2
                if qtmp[j, i] < 0.0
                    # Take mass from west (essentially adjusting fx(i, j)) 
                    qtmp[j, i-1] = qtmp[j, i-1] + qtmp[j, i]
                    qtmp[j, i] = 0.0
                end
            end
            # Note: qtmp(j, 1) could still be negative 
        end

        for j ∈ j1:j2
            for i ∈ 1:im
                q[i, j] = qtmp[j, i]
                q[i, j] = max(0.0, qtmp[j, i]) # (dan 0803) 
            end
        end
    end

    # Check Poles. 
    if jfirst == 1
        ip2 = 0
        # SP 
        if q[1, 1] < 0.0
            # call pfix(q(1, 2), q(1, 1), im, ipx) 
        else
            # Check j=2 
            for i ∈ 1:im
                if q[i, 2] < 0.0
                    ip2 = 1
                    break
                end
            end

            if ip2 != 0
                # call pfix(q(1, 2), q(1, 1), im, ipx) 
            end
        end
    end

    if jlast == jm
        ip2 = 0
        # NP 
        if q[1, jm] < 0.0
            # call pfix(q(1, jm1), q(1, jm), im, ipx) 
        else
            # Check j=jm1 
            for i ∈ 1:im
                if q[i, jm1] < 0.0
                    ip2 = 1
                    break
                end
            end

            if ip2 != 0
                # call pfix(q(1, jm1), q(1, jm), im, ipx) 
            end
        end
    end
end

"""

## Arguments
| Name     | Intent  | Dimensions           | Description
|:---------|:--------|:---------------------|:------------
"""
function fillz!(
    im,
    i1,
    i2,
    km,
    nq,
    q,
    dp
)::Nothing
    # INPUT PARAMETERS: 
    #    integer, intent(in) :: im # No. of longitudes 
    #    integer, intent(in) :: km # No. of levels 
    #    integer, intent(in) :: i1 # Starting longitude 
    #    integer, intent(in) :: i2 # Finishing longitude 
    #    integer, intent(in) :: nq # Total number of tracers 
    #    real, intent(in) ::  dp(im, km) # pressure thickness 

    # INPUT/OUTPUT PARAMETERS: 
    #    real, intent(inout) :: q(im, km, nq) # tracer mixing ratio 

    # LOCAL VARIABLES: 
    i::Integer = k::Integer = ic::Integer = 0
    qup::AbstractFloat = qly::AbstractFloat = dup::AbstractFloat = 0

    for ic ∈ 1:nq
        # Top layer 
        for i ∈ i1:i2
            if q[i, 1, ic] < 0.0
                q[i, 2, ic] = q[i, 2, ic] + q[i, 1, ic] * dp[i, 1] / dp[i, 2]
                q[i, 1, ic] = 0.0
            end
        end

        # Interior 
        for k ∈ 2:km-1
            for i ∈ i1:i2
                if q[i, k, ic] < 0.0
                    # Borrow from above 
                    qup = q[i, k-1, ic] * dp[i, k-1]
                    qly = -q[i, k, ic] * dp[i, k]
                    dup = min(0.5 * qly, qup) # borrow no more than 50 
                    q[i, k-1, ic] = q[i, k-1, ic] - dup / dp[i, k-1]
                    # Borrow from below: q(i, k, ic) is still negative at this stage 
                    q[i, k+1, ic] = q[i, k+1, ic] + (dup - qly) / dp[i, k+1]
                    q[i, k, ic] = 0.0
                end
            end
        end

        # Bottom layer 
        k = km
        for i ∈ i1:i2
            if q[i, k, ic] < 0.0
                # Borrow from above 
                qup = q[i, k-1, ic] * dp[i, k-1]
                qly = -q[i, k, ic] * dp[i, k]
                dup = min(qly, qup)
                q[i, k-1, ic] = q[i, k-1, ic] - dup / dp[i, k-1]
                q[i, k, ic] = 0.0
            end
        end
    end
end

"""

## Arguments
| Name     | Intent  | Dimensions           | Description
|:---------|:--------|:---------------------|:------------
"""
function pfix!(
    q::Vector{AbstractFloat},
    qp::Vector{AbstractFloat},
    im::Integer,
    ipx::Integer
)::Nothing
    summ::AbstractFloat = 0.0
    sump::AbstractFloat = 0.0
    for i ∈ 1:im
        summ += q[i]
        sump += qp[i]
    end

    pmean::AbstractFloat = (sump * gw[1] + summ * gw[2]) / (im * (gw[1] + gw[2]))

    for i ∈ 1:im
        q[i] = pmean
        qp[i] = pmean
    end

    if qp[1] < 0.0
        ipx = 1
    end
end

"""

## Arguments
| Name     | Intent  | Dimensions           | Description
|:---------|:--------|:---------------------|:------------
| `im`     | `inout` |                      |
| `jm`     | `inout` |                      |
| `jfirst` | `inout` |                      |
| `jlast`  | `inout` |                      |
| `q`      | `in`    | `[im, jfirst:jlast]` |
| `qmean`  | `inout` |                      |
"""
function gmean!(
    im::Integer,
    jm::Integer,
    jfirst::Integer,
    jlast::Integer,
    q::Matrix{AbstractFloat},
    qmean::AbstractFloat,
)::Nothing
    @static if @isdefined(SPMD)
        @static if @isdefined(PILGRIM)
            # use parutilitiesmodule, only : parcollective, commglobal, sumop 
        else
            # use mod_comm, only: mp_allgather1d, gid 
        end
    end

    @static if @isdefined(SPMD)
        gsum::AbstractFloat = zeros(jm)
    end

    xsum = zeros(AbstractFloat, jfirst:jlast)

    for j ∈ jfirst:jlast
        xsum[j] = 0.0
        for i ∈ 1:im
            xsum[j] += q[i, j]
        end
        xsum[j] *= gw[j]
    end

    @static if @isdefined(SPMD)
        gsum = 0.0

        @static if @isdefined(PILGRIM)
            for j ∈ jfirst:jlast
                gsum[j] = xsum[j]
            end
            # TODO:
            # call parcollective(commglobal, sumop, jm, gsum) 
        else
            # TODO:
            # call mp_allgather1d(jm, jfirst, jlast, xsum(jfirst), gsum) 
        end

        if gid == 0
            qmean = 0.0
            for j ∈ 1:jm
                qmean += gsum(j)
            end
            qmean /= (2 * im)
        end
    else
        qmean = 0.0
        for j ∈ 1:jm
            qmean += xsum(j)
        end
        qmean /= (2 * im)
    end

end

"""

## Arguments
| Name     | Intent  | Dimensions                     | Description
|:---------|:--------|:-------------------------------|:------------
| `im`     | `in`    |                                |
| `jm`     | `in`    |                                |
| `km`     | `in`    |                                |
| `jfirst` | `in`    |                                |
| `jlast`  | `in`    |                                |
| `ak`     | `in`    | `[km+1]`                       |
| `bk`     | `in`    | `[km+1]`                       |
| `ffsl`   | `in`    | `[jfirst-ng:jlast+ng, km]`     |
| `ps0`    | `inout` | `[im, jfirst:jlast]`           |
| `ps2`    | `in`    | `[im, jfirst-mg:jlast+mg]`     |
| `pe`     | `inout` | `[im, km+1, jfirst:jlast]`     |
| `delp`   | `inout` | `[im, jfirst:jlast, km]`       |
| `fx3`    | `inout` | `[im, jfirst:jlast, km]`       |
| `cx`     | `in`    | `[im, jfirst-ng:jlast+ng, km]` |
| `fy3`    | `inout` | `[im, jfirst:jlast+mg, km]`    |
| `ng`     | `in`    |                                |
| `mg`     | `in`    |                                |
| `tiny`   | `in`    |                                |
| `n_adj`  | `in`    |                                |
"""
function adj_fx!(
    im::Integer,
    jm::Integer,
    km::Integer,
    jfirst::Integer,
    jlast::Integer,
    ak::Vector{AbstractFloat},
    bk::Vector{AbstractFloat},
    ffsl::Matrix{AbstractFloat},
    ps0::Matrix{AbstractFloat},
    ps2::Matrix{AbstractFloat},
    pe::Array{AbstractFloat,3},
    delp::Array{AbstractFloat,3},
    fx3::Array{AbstractFloat,3},
    cx::Array{AbstractFloat,3},
    fy3::Array{AbstractFloat,3},
    ng::Integer,
    mg::Integer,
    tiny::AbstractFloat,
    n_adj::Integer,
)::Nothing
    ps = zeros(AbstractFloat, im, jfirst-mg:jlast+mg)
    fy = zeros(AbstractFloat, im, jfirst:jlast+mg)
    fx = zeros(AbstractFloat, im + 1)
    dps = zeros(AbstractFloat, 0:im)
    dpy = zeros(AbstractFloat, im, jfirst-mg:jlast+mg)

    er0 = er1 = tmpf = dh = dbk = zero(AbstractFloat)
    const fac = 0.25

    js2g0 = max(2, jfirst)::Integer # No ghosting 
    jn2g0 = min(jm - 1, jlast)::Integer # No ghosting 

    for j ∈ jfirst:jlast
        for i ∈ 1:im
            ps[i, j] = ps0[i, j]
        end
    end

    # fx_iteration: 
    for it ∈ 1:n_adj

        @static if @isdefined(SPMD)
            @static if @isdefined(PILGRIM)
                # TODO:
                # call parbegintransfer(pattern2dmg, ps, ps) 
                # call parendtransfer(pattern2dmg, ps, ps) 
            else
                # TODO:
                # call mp_send3d_ns(im, jm, jfirst, jlast, 1, 1, mg, mg, ps, 2) 
                # call mp_barrier() 
                # call mp_recv3d_ns(im, jm, jfirst, jlast, 1, 1, mg, mg, ps, 2) 
                # call mp_barrier() 
            end
        end

        # adjust fx 
        for k ∈ 3:km
            dbk = bk[k+1] - bk[k]

            if dbk > 0.001
                for j ∈ js2g0:jn2g0
                    for i ∈ 1:im
                        dps[i] = (ps[i, j] - ps2[i, j]) * dbk
                    end

                    dps[0] = dps[im]

                    for i ∈ 1:im
                        fx[i] = fac * (dps[i-1] - dps[i])
                        tmpf = fx3[i, j, k] + fx[i]

                        if tmpf * fx3[i, j, k] > 0.0
                            fx3[i, j, k] = tmpf
                        else
                            fx[i] = fx3[i, j, k]
                            fx3[i, j, k] = sign(fx3[i, j, k]) * min(abs(tmpf), abs(fx3[i, j, k]))
                            fx[i] = fx3[i, j, k] - fx[i]
                        end
                    end

                    fx[im+1] = fx[1]

                    # update delp 
                    for i ∈ 1:im
                        delp[i, j, k] = delp[i, j, k] + fx[i] - fx[i+1]
                    end
                end # j - loop 

                # adjust fy 

                # Need ps at jlast + 1 
                for j ∈ max(jfirst - 1, 1):min(jm, jlast + 1)
                    for i ∈ 1:im
                        dpy[i, j] = (ps[i, j] - ps2[i, j]) * dbk * gw[j]
                    end
                end

                for j ∈ js2g0:min(jm, jlast + 1)
                    for i ∈ 1:im
                        fy[i, j] = fac * (dpy[i, j-1] - dpy[i, j])
                        tmpf = fy3[i, j, k] + fy[i, j]

                        if tmpf * fy3[i, j, k] > 0.0
                            fy3[i, j, k] = tmpf
                        else
                            fy[i, j] = fy3[i, j, k]
                            fy3[i, j, k] = sign(fy3[i, j, k]) * min(abs(tmpf), abs(fy3[i, j, k]))
                            fy[i, j] = fy3[i, j, k] - fy[i, j]
                        end
                    end
                end

                # update delp 
                for j ∈ js2g0:jn2g0
                    for i ∈ 1:im
                        delp[i, j, k] = delp[i, j, k] + (fy[i, j] - fy[i, j+1]) * rgw[j]
                    end
                end

                # Poles
                if jfirst == 1
                    for i ∈ 1:im
                        delp[i, 1, k] = delp[i, 1, k] - fy(i, 2) * rgw[1]
                    end
                    # Left over comment from fotran -> julia migration
                    # call xpavg(delp[:, 1, k], im) 
                end

                if jlast == jm
                    for i ∈ 1:im
                        delp[i, jm, k] = delp[i, jm, k] + fy[i, jm] * rgw[jm]
                    end
                    # Left over comment from fotran -> julia migration
                    # call xpavg(delp[:, jm, k], im)
                end
            end
        end # k - loop 

        # Update pe and ps

        for j ∈ jfirst:jlast
            for i ∈ 1:im
                pe[i, 1, j] = ak[1]
            end

            for k ∈ 1:km
                for i ∈ 1:im
                    pe[i, k+1, j] = pe[i, k, j] + delp[i, j, k]
                end
            end

            for i ∈ 1:im
                ps[i, j] = pe[i, km+1, j]
            end
        end

    end # fx_iteration

    for j ∈ js2g0:jn2g0
        for k ∈ km:-1:3
            dbk = bk[k+1] - bk[k]
            if dbk > 0.001
                for i ∈ 1:im
                    dps[i] = (ps[i, j] - ps2[i, j]) * dbk
                end

                dps[0] = dps[im]

                i = 1
                er0 = dps[i-1]
                er1 = dps[i]

                if er0 * er1 < 0.0
                    if er1 > 0.0
                        dh = min(-er0, er1)
                        fx3[i, j, k] = fx3[i, j, k] - dh
                        delp[im, j, k] = delp[im, j, k] + dh
                        delp[i, j, k] = delp[i, j, k] - dh
                    else
                        dh = min(er0, -er1)
                        fx3[i, j, k] = fx3[i, j, k] + dh
                        delp[im, j, k] = delp[im, j, k] - dh
                        delp[i, j, k] = delp[i, j, k] + dh
                    end
                end

                for i ∈ 2:im
                    er0 = dps[i-1]
                    er1 = dps[i]

                    if er0 * er1 < 0.0
                        if er1 > 0.0
                            dh = min(-er0, er1)
                            fx3[i, j, k] = fx3[i, j, k] - dh
                            delp[i-1, j, k] = delp[i-1, j, k] + dh
                            delp[i, j, k] = delp[i, j, k] - dh
                        else
                            dh = min(er0, -er1)
                            fx3[i, j, k] = fx3[i, j, k] + dh
                            delp[i-1, j, k] = delp[i-1, j, k] - dh
                            delp[i, j, k] = delp[i, j, k] + dh
                        end
                    end
                end
            end
        end # k - loop 

        for i ∈ 1:im
            pe[i, 1, j] = ak[1]
        end

        for k ∈ 1:km
            for i ∈ 1:im
                pe[i, k+1, j] = pe[i, k, j] + delp[i, j, k]
            end
        end

        for i ∈ 1:im
            ps[i, j] = pe[i, km+1, j]
        end

        for k ∈ 1:km
            if ffsl[j, k]
                for i ∈ 1:im
                    fx3[i, j, k] = fx3[i, j, k] / (sign(cx[i, j, k]) * max(abs(cx[i, j, k]), tiny))
                end
            end
        end
    end

    # Copy adjusted surface pressure 
    for j ∈ jfirst:jlast
        for i ∈ 1:im
            ps0[i, j] = ps[i, j]
        end
    end

end

end

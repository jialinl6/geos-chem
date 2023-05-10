# """
# ------------------------------------------------------------------------------
#                   GEOS-Chem Global Chemical Transport Model                  !
# ------------------------------------------------------------------------------
# BOP

#  !MODULE: tpcore_fvdas_mod.F90

#  !DESCRIPTION: \subsection*{Overview}
#   Module Tpcore\_Fvdas\_Mod contains routines for the TPCORE
#   transport scheme, as implemented in the GMI model (cf. John Tannahill),
#   based on Lin \ Rood 1995.  The Harvard Atmospheric Chemistry Modeling Group
#   has added modifications to implement the Philip-Cameron Smith pressure
#   fixer for mass conservation.  Mass flux diagnostics have also been added.

# \subsection*{References}
#   \begin{enumerate}
#   \item Lin, S.-J., and R. B. Rood, 1996: \emph{Multidimensional flux
#          form semi-Lagrangian transport schemes},
#          \underline{ Mon. Wea. Rev.}, \textbf{124}, 2046-2070.
#   \item Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994:
#          \emph{A class of the van Leer-type transport schemes and its
#          applications to the moisture transport in a General Circulation
#          Model}, \underline{Mon. Wea. Rev.}, \textbf{122}, 1575-1593.
#   \end{enumerate}

# \subsection*{Selecting E/W, N/S and vertical advection options}

#   The flags IORD, JORD, KORD select which transport schemes are used in the
#   E/W, N/S, and vertical directions, respectively.  Here is a list of the
#   possible values that IORD, JORD, KORD may be set to (original notes from
#   S-J Lin):

#   \begin{enumerate}
#   \item 1st order upstream scheme (too diffusive, not a real option;
#          it can be used for debugging purposes; this is THE only known
#          "linear" monotonic advection scheme.).
#   \item 2nd order van Leer (full monotonicity constraint;
#          see Lin et al 1994, MWR)
#   \item monotonic PPM* (Collela \& Woodward 1984)
#   \item semi-monotonic PPM (same as 3, but overshoots are allowed)
#   \item positive-definite PPM (constraint on the subgrid distribution is
#          only strong enough to prevent generation of negative values;
#          both overshoots \& undershoots are possible).
#   \item un-constrained PPM (nearly diffusion free; faster but
#          positivity of the subgrid distribution is not quaranteed. Use
#          this option only when the fields and winds are very smooth.
#   \item Huynh/Van Leer/Lin full monotonicity constraint.  Only KORD can be
#          set to 7 to enable the use of Huynh's 2nd monotonicity constraint
#          for piece-wise parabolic distribution.
#   \end {enumerate}

#   Recommended values:

#   \begin{itemize}
#   \item IORD=JORD=3 for high horizontal resolution.
#   \item KORD=3 or 7
#   \end{itemize}

#   The implicit numerical diffusion decreases as \_ORD increases.
#   DO NOT use option 4 or 5 for non-positive definite scalars
#   (such as Ertel Potential Vorticity).
# \\
# \\
#  In GEOS-Chem we have been using IORD=3, JORD=3, KORD=7.  We have tested
#  the OpenMP parallelization with these options.  GEOS-Chem users who wish to
#  use different (I,J,K)ORD options should consider doing single-procsessor
#  vs. multi-processor tests to test the implementation of the parallelization.

# \subsection*{GEOS-4 and GEOS-5 Hybrid Grid Definition}

#   For GEOS-4 and GEOS-5 met fields, the pressure at the bottom edge of
#   grid box (I,J,L) is defined as follows:

#      \$\$P_{edge}(I,J,L) = A_{k}(L) + [ B_{k}(L) * P_{surface}(I,J) ]\$\$

#   where

#   \begin{itemize}
#   \item \$P_{surface}\$(I,J) is the "true" surface pressure at lon,lat (I,J)
#   \item \$A_{k}\$(L) has the same units as surface pressure [hPa]
#   \item \$B_{k}\$(L) is a unitless constant given at level edges
#   \end{itemize}

#   \$A_{k}(L)\$ and \$B_{k}(L)\$ are supplied to us by GMAO.
# \\
# \\
#  !REMARKS:
#   Ak(L) and Bk(L) are defined at layer edges.


#                   /////////////////////////////////
#               / \ ------ Model top P=ak(1) --------- ak(1), bk(1)
#                |
#     delp(1)    |  ........... q(i,j,1) ............
#                |
#               \ / ---------------------------------  ak(2), bk(2)



#               / \ ---------------------------------  ak(k), bk(k)
#                |
#     delp(k)    |  ........... q(i,j,k) ............
#                |
#               \ / ---------------------------------  ak(k+1), bk(k+1)



#               / \ ---------------------------------  ak(km), bk(km)
#                |
#     delp(km)   |  ........... q(i,j,km) .........
#                |
#               \ / -----Earth's surface P=Psfc ------ ak(km+1), bk(km+1)
#                  //////////////////////////////////

#  Note: surface pressure can be of any unit (e.g., pascal or mb) as
#  long as it is consistent with the definition of (ak, bk) defined above.

# Winds (u,v), ps, and q are assumed to be defined at the same points.

# The latitudes are given to the initialization routine: init_tpcore.

# ## Author
# Original code from Shian-Jiann Lin, GMAO \\
# Modified for GMI model by John Tannahill, LLNL (jrt@llnl.gov) \\
# Implemented into GEOS-Chem by Claire Carouge (ccarouge@seas.harvard.edu) \\
# ProTeX documentation added by Bob Yantosca (yantosca@seas.harvard.edu) \\
# OpenMP parallelization added by Bob Yantosca (yantosca@seas.harvard.edu)

# ## Revision History
# 05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from the GMI model. This eliminates the polar overshoot in the stratosphere. \\
# See https://github.com/geoschem/geos-chem for complete history
# """
module tpcore_fvdas_mod

export init_tpcore!, exit_tpcore!, tpcore_fvdas!

# TODO: On integratiion, these should persist between loads of the module
global tpcore_fvdas_mod_dtdx5::Vector{Float64}
global tpcore_fvdas_mod_dtdy5::Vector{Float64}
global tpcore_fvdas_mod_cosp::Vector{Float64}
global tpcore_fvdas_mod_cose::Vector{Float64}
global tpcore_fvdas_mod_gw::Vector{Float64}
global tpcore_fvdas_mod_dlat::Vector{Float64}

mutable struct StateChm
end
mutable struct StateDiag
end

"""
Allocates and initializes all module variables.

## Arguments
- `im` (`in`) - Global E-W dimension
- `jm` (`in`) - Global N-S dimension
- `km` (`in`) - Vertical dimension
- `jfirst` (`out`) - Local first index for N-S axis
- `jlast` (`out`) - Local last index for N-S axis
- `ng` (`in`) - large ghost width
- `mg` (`in`) - small ghost width
- `dt` (`in`) - Time step in seconds
- `ae` (`in`) - Earth's radius (m)
- `clat` (`in`; `[jm]`) - latitude in radian
- `rc` (`out`) - Success or failure

## Revision History
05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function init_tpcore!(
    im::Int64,
    jm::Int64,
    km::Int64,
    jfirst::Ref{Int64},
    jlast::Ref{Int64},
    ng::Int64,
    mg::Int64,
    dt::Float64,
    ae::Float64,
    clat::Array{Float64,1},
    rc::Ref{Int64}
)::Nothing
    # cell edge latitude in radian
    elat::Vector{Float64} = zeros(jm + 1)
    sine::Vector{Float64} = zeros(jm + 1)
    sine_25::Vector{Float64} = zeros(jm + 1)
    dlon::Float64 = 0.0
    i::Int64, j::Int64 = 0, 0

    rc[] = GC_SUCCESS
    errmsg = ""
    thisloc = " -> at Init_Tpcore (in module GeosCore/tpcore_fvas_mod.F90)"

    # Since we are not using MPI parallelization, we can set JFIRST and JLAST to the global grid limits in latitude. (bmy, 12/3/08)
    jfirst[] = 1
    jlast[] = jm

    if jlast[] - jfirst[] < 2
        println("Minimum size of subdomain is 3")
    end

    global tpcore_fvdas_mod_dtdx5 = zeros(jm)
    global tpcore_fvdas_mod_dtdy5 = zeros(jm)
    global tpcore_fvdas_mod_cosp = zeros(jm)
    global tpcore_fvdas_mod_cose = zeros(jm)
    global tpcore_fvdas_mod_gw = zeros(jm)
    # For PJC pressure-fixer
    global tpcore_fvdas_mod_dlat = zeros(jm)

    dlon = 2.0 * pi / im

    elat[1] = -0.5 * pi
    sine[1] = -1.0
    sine_25[1] = -1.0
    tpcore_fvdas_mod_cose[1] = 0.0

    for j in 2:jm
        elat[j] = 0.5 * (clat[j-1] + clat[j])
        sine[j] = sin(elat[j])
        sine_25[j] = sin(clat[j])
        tpcore_fvdas_mod_cose[j] = cos(elat[j])
    end

    # N. Pole
    elat[jm+1] = 0.5 * π
    sine[jm+1] = 1.0
    sine_25[jm+1] = 1.0

    # Polar cap (S. Pole)
    tpcore_fvdas_mod_dlat[1] = 2.0 * (elat[2] - elat[1])
    for j ∈ 2:jm-1
        tpcore_fvdas_mod_dlat[j] = elat[j+1] - elat[j]
    end

    # Polar cap (N. Pole)
    tpcore_fvdas_mod_dlat[jm] = 2.0 * (elat[jm+1] - elat[jm])

    for j in 1:jm
        tpcore_fvdas_mod_gw[j] = sine[j+1] - sine[j]
        tpcore_fvdas_mod_cosp[j] = tpcore_fvdas_mod_gw[j] / tpcore_fvdas_mod_dlat[j]
        tpcore_fvdas_mod_dtdx5[j] = 0.5 * dt / (dlon * ae * tpcore_fvdas_mod_cosp[j])
        tpcore_fvdas_mod_dtdy5[j] = 0.5 * dt / (ae * tpcore_fvdas_mod_dlat[j])
    end

    # Echo info to stdout
    println("=" * 79)
    println("TPCORE_FVDAS (based on GMI) Tracer Transport Module successfully initialized")
    println("=" * 79)
end

"""
Deallocates all module variables.

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function exit_tpcore!()::Nothing
    # ! Deallocate arrays only if they are allocated
    # IF ( ALLOCATED( tpcore_fvdas_mod_COSP   ) ) DEALLOCATE( tpcore_fvdas_mod_COSP   )
    # IF ( ALLOCATED( tpcore_fvdas_mod_COSE   ) ) DEALLOCATE( tpcore_fvdas_mod_COSE   )
    # IF ( ALLOCATED( tpcore_fvdas_mod_GW     ) ) DEALLOCATE( tpcore_fvdas_mod_GW     )
    # IF ( ALLOCATED( tpcore_fvdas_mod_DTDX5  ) ) DEALLOCATE( tpcore_fvdas_mod_DTDX5  )
    # IF ( ALLOCATED( tpcore_fvdas_mod_DTDY5  ) ) DEALLOCATE( tpcore_fvdas_mod_DTDY5  )
    # IF ( ALLOCATED( tpcore_fvdas_mod_DLAT   ) ) DEALLOCATE( tpcore_fvdas_mod_DLAT   )
end

#   END SUBROUTINE Exit_Tpcore

tpcore_fvdas_save::Bool = true
# controls various options in E-W; N-S; vertcal advection
tpcore_fvdas_ilmt::Int, tpcore_fvdas_jlmt::Int, tpcore_fvdas_klmt::Int = 0, 0, 0

# 2=floating pressure
const TPCORE_FVDAS_ADVEC_CONSRV_OPT = 2
const TPCORE_FVDAS_CROSS = true

"""
Takes horizontal winds on sigma (or hybrid sigma-p) surfaces and calculates mass fluxes, and then updates the 3D mixing ratio fields one time step (tdt).  The basic scheme is a Multi-Dimensional Flux Form Semi-Lagrangian (FFSL) based on the van Leer or PPM (see Lin and Rood, 1995).

## Arguments
- `dt` (`in`) - time step [s]
- `ae` (`in`) - Earth's radius [m]
- `im` (`in`) - Global E-W dimension
- `jm` (`in`) - Global N-S dimension
- `km` (`in`) - Global vertical dimension
- `jfirst` (`in`) - Latitude index for local first box (NOTE: for global grids these are 1 and JM, respectively)
- `jlast` (`in`) - Latitude index for local last box (NOTE: for global grids these are 1 and JM, respectively)
- `ng` (`in`) - Primary ghost region
- `mg` (`in`) - Secondary ghost region
- `nq` (`in`) - Ghosted latitudes (3 required by PPM)
- `ak` (`in`; `[km+1]`) - Ak coordinates to specify the hybrid grid (see the REMARKS section below)
- `bk` (`in`; `[km+1]`) - Bk coordinates to specify the hybrid grid (see the REMARKS section below)
- `u` (`in`; `[:,:,:]`) - u-wind (m/s) at mid-time-level (t=t+dt/2)
- `v` (`in`; `[:,:,:]`) - v-wind (m/s) at mid-time-level (t=t+dt/2)
- `ps1` (`in`; `[im,jfirst:jlast]`) - surface pressure at current time
- `ps2` (`in`; `[im,jfirst:jlast]`) - surface pressure at future time=t+dt
- `ps` (`out`; `[im,jfirst:jlast]`) - "Predicted" surface pressure [hPa]
- `iord` (`in`) - Flags to denote E-W transport scheme
- `jord` (`in`) - Flags to denote N-S transport scheme
- `kord` (`in`) - Flags to denote vertical transport scheme
- `n_adj` (`in`) - Number of adjustments to air_mass_flux (0 = no adjustment)
- `xmass` (`in`; `[:,:,:]`) - E/W mass fluxes [kg/s] (These are computed by the pressure fixer, and passed into TPCORE)
- `ymass` (`in`; `[:,:,:]`) - N/S mass fluxes [kg/s] (These are computed by the pressure fixer, and passed into TPCORE)
- `fill` (`in`) - Fill negatives ?
- `area_m2` (`in`; `[jm]`) - Grid box surface area for mass flux diag [m2]
- `state_chm` (`inout`) - Diagnostics state object
- `state_diag` (`inout`) - Diagnostics state object

## Author
Original code from Shian-Jiann Lin, DAO) \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function tpcore_fvdas!(
    dt::Float64,
    ae::Float64,
    im::Int64,
    jm::Int64,
    km::Int64,
    jfirst::Int64,
    jlast::Int64,
    ng::Int64,
    mg::Int64,
    nq::Int64,
    ak::Vector{Float64},
    bk::Vector{Float64},
    u::Array{Float64,3},
    v::Array{Float64,3},
    ps1::Matrix{Float64},
    ps2::Matrix{Float64},
    ps::Matrix{Float64},
    iord::Int64,
    jord::Int64,
    kord::Int64,
    n_adj::Int64,
    xmass::Array{Float64,3},
    ymass::Array{Float64,3},
    fill::Bool,
    area_m2::Vector{Float64},
    state_chm::StateChm,
    state_diag::StateDiag
)::Nothing
    rj2m1::Int64
    j1p::Int64, j2p::Int64
    jn::Vector{Int64} = zeros(km)
    js::Vector{Int64} = zeros(km)
    il::Int64, ij::Int64, ik::Int64, iq::Int64, k::Int64, j::Int64, i::Int64, Kflip::Int64
    num::Int64, k2m1::Int64, s::Int64
    north::Int64, south::Int64

    dap::Vector{Float64} = zeros(km)
    dbk::Vector{Float64} = zeros(km)

    # E-W CFL # on C-grid
    cx::Array{Float64,3} = zeros(im, jfirst-ng:jlast+ng, km)
    # N-S CFL # on C-grid
    cy::Array{Float64,3} = zeros(im, jfirst:jlast+mg, km)

    delp1::Array{Float64,3} = zeros(im, jm, km)
    delp2::Array{Float64,3} = zeros(im, jm, km)
    delpm::Array{Float64,3} = zeros(im, jm, km)
    pu::Array{Float64,3} = zeros(im, jm, km)
    dpi::Array{Float64,3} = zeros(im, jm, km)

    # geometrical factor for meridional
    geofac::Vector{Float64} = zeros(jm)
    # geometrical gactor for poles.
    geofac_pc::Float64 = 0.0

    dp::Float64 = 0.0
    dps_ctm::Array{Float64,2} = zeros(im, jm)
    ua::Array{Float64,3} = zeros(im, jm, km)
    va::Array{Float64,3} = zeros(im, jm, km)
    wz::Array{Float64,3} = zeros(im, jm, km)
    dq1::Array{Float64,3} = zeros(im, jfirst-ng:jlast+ng, km)

    # qqu, qqv, adx and ady are now 2d arrays for parallelization purposes. (ccc, 4/1/08)
    qqu::Array{Float64,2} = zeros(im, jm)
    qqv::Array{Float64,2} = zeros(im, jm)
    adx::Array{Float64,2} = zeros(im, jm)
    ady::Array{Float64,2} = zeros(im, jm)

    # fx, fy, fz and qtp are now 4D arrays for parallelization purposes. (ccc, 4/1/09)
    fx::Array{Float64,4} = zeros(im, jm, km, nq)
    # one more for edges
    fy::Array{Float64,4} = zeros(im, jm + 1, km, nq)
    fz::Array{Float64,4} = zeros(im, jm, km, nq)

    js2g0::Int, jn2g0::Int = 0, 0

    # Add pointer to avoid array temporary in call to FZPPM (bmy, 6/5/13)
    q_ptr::Array{Float64,3} = zeros()

    # Add definition of j1p and j2p for enlarge polar cap. (ccc, 11/20/08)
    j1p = 3
    j2p = jm - j1p + 1

    # Average surf. pressures in the polar cap. (ccc, 11/20/08)
    average_press_poles!(area_m2, ps1, 1, im, 1, jm, 1, im, 1, jm)
    average_press_poles!(area_m2, ps2, 1, im, 1, jm, 1, im, 1, jm)

    # Calculation of some geographic factors. (ccc, 11/20/08)
    rj2m1 = jm - 1
    dp = π / rj2m1

    geofac[1:jm] .= dp ./ (2.0 .* area_m2[1:jm] ./ (sum(area_m2) * im) .* im)
    geofac_pc = dp / (2.0 * (sum(area_m2[1:2]) / (sum(area_m2) * im)) * im)

    if first
        first = false
        set_lmts!(tpcore_fvdas_ilmt, tpcore_fvdas_jlmt, tpcore_fvdas_klmt, im, jm, iord, jord, kord)
    end

    for ik ∈ 1:km
        dap[ik] = ak[ik+1] - ak[ik]
        dbk[ik] = bk[ik+1] - bk[ik]
    end

    for ik ∈ 1:km
        set_cross_terms!(dap[ik], dbk[ik], ps1, ps2, delp1[:, :, ik], delpm[:, :, ik], pu[:, :, ik], 1, jm, 1, im, 1, jm, j1p, j2p, 1, im, 1, jm)

        if j1p != 1 + 1
            for iq ∈ 1:nq
                q_ptr = state_chm.species[iq].conc[:, :, km:-1:1]
                average_const_poles!(dap[ik], dbk[ik], area_m2, ps1, q_ptr, 1, jm, im, 1, im, 1, jm, 1, im, 1, jm)
            end
        end

        calc_courant!(tpcore_fvdas_mod_cose, delpm[:, :, ik], pu[:, :, ik], xmass[:, :, ik], ymass[:, :, ik], cx[:, :, ik], cy[:, :, ik], j1p, j2p, 1, jm, 1, im, 1, jm, 1, im, 1, jm)

        calc_divergence!(true, geofac_pc, geofac, dpi[:, :, ik], xmass[:, :, ik], ymass[:, :, ik], j1p, j2p, 1, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm)

        set_cross_terms!(cx(:, :, ik), cy(:, :, ik), ua(:, :, ik), va(:, :, ik), j1p, j2p, 1, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm, CROSS)
    end

    dps_ctm .= sum(dpi, dims=3)

    calc_vert_mass_flux!(dbk, dps_ctm, dpi, wz, 1, im, 1, jm, 1, km)

    set_jn_js!(jn, js, cx, 1, im, 1, jm, 1, jm, j1p, j2p, 1, im, 1, jm, 1, km)

    if advec_consrv_opt == 0
        for ik ∈ 1:km
            for ij ∈ 1:jm
                for il ∈ 1:im
                    delp2[il, ij, ik] = dap[ik] + dbk[ik] * (ps1[il, ij] + dps_ctm[il, ij])
                end
            end
        end
    elseif advec_consrv_opt == 1 || advec_consrv_opt == 2
        for ik ∈ 1:km
            for ij ∈ 1:jm
                for il ∈ 1:im
                    delp2[il, ij, ik] = dap[ik] + dbk[ik] * ps2[il, ij]
                end
            end
        end
    end

    # calculate surf. pressure at t+dt. (ccc, 11/20/08)
    ps .= ak[1] .+ sum(delp2, dims=3)

    # For time optimization : we parallelize over tracers and we loop over the levels outside horizontal transport subroutines. (ccc, 4/1/09)
    # Also zeroed PRIVATE variables within the loop, and set jn(ik) and js(ik) to PRIVATE loop variables.  This seems to avoid small diffs in output. -- Bob Yantosca (04 Jan 2022)
    for iq ∈ 1:nq
        q_ptr = state_chm.species[iq].conc[:, :, km:-1:1]

        for ik ∈ 1:km
            adx = 0.0
            ady = 0.0
            qqu = 0.0
            qqv = 0.0

            # Northernmost and southernmost latitude indices by level
            north = jn[ik]
            south = js[ik]

            # .sds. convert to "mass"
            dq1[:, :, ik] .= q_ptr[:, :, ik] .* delp1[:, :, ik]

            calc_advec_cross_terms!(north, south, q_ptr[:, :, ik], qqu, qqv, ua[:, :, ik], va[:, :, ik], j1p, j2p, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm, CROSS)

            # Add advective form E-W operator for E-W cross terms.
            xadv_dao2!(2, north, south, adx, qqv, ua(:, :, ik), 1, im, 1, jm, 1, jm, j1p, j2p, 1, im, 1, jm)

            # Add advective form N-S operator for N-S cross terms.
            # use a polar cap of 2 boxes (i.e. the "2" as the first argument to YADV_DAO2.  The older TPCORE only had a polar cap of 1 box (just the Pole itself).  Claire figured this out.  (bmy, 12/11/08)
            # TODO: Possible mistranslation
            yadv_dao2!(2, ady, qqu, va[:, :, ik], j1p, j2p, 1, im, 1, jm, 1, jm, 1, im)

            # update constituent array qq1 by adding in cross terms
            q_ptr[:, :, ik] .= q_ptr[:, :, ik] .+ ady .+ adx

            xtp!(tpcore_fvdas_ilmt, north, south, pu[:, :, ik], cx[:, :, ik], dq1[:, :, ik], qqv, xmass[:, :, ik], fx(:, :, ik, iq), j1p, j2p, im, 1, jm, 1, im, 1, jm, 1, im, 1, jm, iord)

            # TODO: Possible mistranslation
            ytp!(tpcore_fvdas_jlmt, geofac_pc, geofac, cy[:, :, ik], dq1[:, :, ik], qqu, qqv, ymass[:, :, ik], fy[:, :, ik, iq], j1p, j2p, 1, im, 1, jm, im, 1, im, 1, jm, 1, im, 1, jm, jord)
        end

        fzppm!(tpcore_fvdas_klmt, delp1, wz, dq1, q_ptr, fz[:, :, :, iq], j1p, 1, jm, 1, im, 1, jm, im, km, 1, im, 1, jm, 1, km)

        if fill
            qckxyz!(dq1, j1p, j2p, 1, jm, 1, im, 1, jm, 1, im, 1, jm, 1, km)
        end

        q_ptr .= dq1 ./ delp2

        if j1p != 2
            q_ptr[:, 2, :] .= q_ptr[:, 1, :]
            q_ptr[:, jm-1, :] .= q_ptr[:, jm, :]
        end

        # MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
        # Set tracer concentration to a small positive number if concentration is negative. Negative concentration may occur at the poles. This is an issue that should be looked into in the future. (ewl, 6/30/15)
        for i ∈ 1:im
            for j ∈ 1:jm
                for k ∈ 1:km
                    if q_ptr[i, j, k] < 0.0
                        q_ptr[i, j, k] = 1.0e-26
                    end
                end
            end
        end

        # TODO: Reconsider translation
        # q_ptr => NULL()
    end

    # MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
    # HISTORY (aka netCDF diagnostics)
    # E/W flux of advected species [kg/s]  (ccarouge 12/2/08)
    # The unit conversion is:
    #
    # Mass    P diff     100      1       area of     kg tracer      1
    # ----- = in grid *  ---  *  ---   *  grid box * ------------ * ---
    # time    box         1       g       AREA_M2    kg moist air    s
    #
    #  kg      hPa     Pa     s^2    m^2        1
    # ----  = ----- * ----- * ---- * ---- * --------
    #  s        1      hPa     m      1       DeltaT

    if state_diag.archive_advfluxzonal
        # Zero netCDF diagnostic array
        state_diag.advfluxzonal .= 0.0

        # Calculate fluxes for diag. (ccc, 11/20/08)
        # No ghosting
        js2g0 = max(j1p, jfirst)
        # No ghosting
        jn2g0 = min(j2p, jlast)

        # Loop over diagnostic slots
        for s ∈ 1:state_diag.map_advfluxzonal.nslots
            # Get the advectId from the slotId
            iq = state_diag.map_advfluxzonal.slot2id[s]

            for k ∈ 1:km
                for j ∈ js2g0:jn2g0
                    for i ∈ 1:im
                        # Units: [kg/s]
                        # But consider changing to area-independent units [kg/m2/s]
                        kflip = km - k + 1 # flip vert
                        state_diag.advfluxzonal[i, j, kflip, s] = fx[i, j, k, iq] * area_m2[j] * g0_100 / dt
                    end
                end
            end
        end
    end

    # MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
    # HISTORY (aka netCDF diagnostics)
    # N/S flux of tracer [kg/s]
    # (bdf, bmy, 9/28/04, ccarouge 12/12/08)
    # NOTE, the unit conversion is the same as desciribed above for the E-W diagnostics.  The geometrical factor was already applied to fy in Ytp. (ccc, 4/1/09)
    if state_diag.archive_advfluxmerid
        # Zero netCDF diagnostic array
        state_diag.advfluxmerid .= 0.0

        # Loop over diagnostic slots
        for s ∈ 1:state_diag.map_advfluxmerid.nslots
            # Get the advectId from the slotId
            iq = state_diag.map_advfluxmerid.slot2id[s]

            for k ∈ 1:km
                for j ∈ 1:jm
                    for i ∈ is2g0:in2g0
                        # Compute mass flux [kg/s]
                        # Units: [kg/s]
                        # But consider changing to area-independent units [kg/m2/s]

                        # flip vert
                        kflip = km - k + 1
                        state_diag.advfluxmerid[i, j, kflip, s] = fy[i, j, k, iq] * area_m2[j] * g0_100 / dt
                    end
                end
            end
        end
    end

    # MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
    # HISTORY (aka netCDF diagnostics)
    # Up/down flux of tracer [kg/s]
    # (bmy, bdf, 9/28/04, ccarouge 12/2/08)
    # The vertical transport done in qmap.  We need to find the difference in order to to interpret transport.
    # Break up diagnostic into up & down fluxes using the surface boundary conditions.  Start from top down (really surface up for flipped TPCORE)
    # By construction, MASSFLUP is flux into the bottom of the box. The flux at the bottom of KM (the surface box) is not zero by design. (phs, 3/4/08)
    if state_diag.archive_advfluxvert
        # Zero netCDF diagnostic array
        state_diag.advfluxvert .= 0.0

        # Loop over diagnostic slots
        for s ∈ 1:state_diag.map_advfluxvert.nslots
            # Get the advectId from the slotId
            iq = state_diag.map_advfluxvert.slot2id[s]

            # Loop over grid boxes
            for k ∈ 1:km
                for j ∈ 1:jm
                    for i ∈ 1:im
                        # Ilya Stanevic (stanevic@atmosp.physics.utoronto.ca) says that using FZ for the ND26 diagnostic should be correct.  He writes:
                        # > To be safe you can use FZ variable from Fzppm. That is the real vertical species mass. And it is zero at the surface.
                        # > To be clear, Fz is present only in the tpcore for low resolution GLOBAL model (4x5, 2x2.5). Nested model has different way to calculate vertical advection and there is no such thing as FZ. Therefore, we have to calculate the species mass difference in the box before and after vertical advection in order to get vertical mass flux.
                        # > -- Bob Yantosca (28 Mar 2017)

                        # Units: [kg/s]
                        # But consider changing to area-independent units [kg/m2/s]

                        # flip vert
                        kflip = km - k + 1
                        state_diag.advfluxvert[i, j, kflip, s] = fz[i, j, k, iq] * area_m2[j] * g0_100 / dt
                    end
                end
            end
        end
    end
end

"""
Averages the species concentrations at the Poles when the Polar cap is enlarged.  It makes the last two latitudes equal.

## Arguments
- `dap` (`in`) - Pressure difference across layer from (ai * pt) term [hPa]
- `dbk` (`in`) - Difference in bi across layer - the dSigma term
- `rel_area` (`in`, `[ju1:j2]`) - Relative surface area of grid box [fraction]
- `pctm1` (`in`, `[ilo:ihi, julo:jhi]`) - CTM surface pressure at t1 [hPa]
- `const1` (`inout`, `[i1:i2, ju1:j2]`) - Species concentration, known at zone center [mixing ratio]
- `ju1_gl` (`in`) - Global latitude indices of the South Pole
- `j2_gl` (`in`) - Global latitude indices of the North Pole
- `i2_gl` (`in`) - Global max longitude index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO) \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function average_const_poles!(
    dap::Float64,
    dbk::Float64,
    rel_area::Vector{Float64},
    pctm1::Vector{Float64},
    const1::Vector{Float64},
    ju1_gl::Int64,
    j2_gl::Int64,
    i2_gl::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64
)::Nothing
    ik::Int64, il::Int64 = 0, 0
    meanq = 0.0

    # pressure thickness at North Pole, the psudo-density in a hydrostatic system at t1 (mb)
    delp1n::Matrix{Float64} = zeros(i2:i1, j2-1:j2)
    # pressure thickness at South Pole, the psudo-density in a hydrostatic system at t1 (mb)
    delp1s::Matrix{Float64} = zeros(i2:i1, ju1:ju1+1)

    if ju1 == ju1_gl
        delp1s[i1:i2, ju1:ju1+1] .= dap .+ dbk .* pctm1[i1:i2, ju1:ju1+1]

        sum1::Float64 = 0.0
        sum2::Float64 = 0.0
        for il in i1:i2
            sum1 += sum(const1[il, ju1:ju1+1] .* delp1s[il, ju1:ju1+1] .* rel_area[ju1:ju1+1]) / (sum(rel_area) * i2_gl)
            sum2 += sum(delp1s[il, ju1:ju1+1] .* rel_area[ju1:ju1+1]) / (sum(rel_area) * i2_gl)
        end

        meanq::Float64 = sum1 / sum2

        const1[:, ju1:ju1+1] .= meanq
    end
end

"""
Sets the cross terms for E-W horizontal advection.
    
## Arguments
- `crx` (`in`; `[ilo:ihi, `julo:jhi`]`) - Courant number in E-W direction
- `cry` (`in`; `[ilo:ihi, `julo:jhi`]`) - Courant number in N-S direction
- `ua` (`out`; `[ilo:ihi, `julo:jhi`]`) - Average of Courant numbers from il and il+1
- `va` (`out`; `[ilo:ihi, `julo:jhi`]`) - Average of Courant numbers from ij and ij+1
- `j1p` (`in`) - Global latitude indices at the edges of the S/N polar caps
- `j2p` (`in`) - Global latitude indices at the edges of the S/N polar caps
- `i1_gl` (`in`) - Global min longitude (I) indices
- `i2_gl` (`in`) - Global max longitude (I) indices
- `ju1_gl` (`in`) - Global min latitude (J) indices
- `j2_gl` (`in`) - Global max latitude (J) indices
- `ilo` (`in`) - Local min longitude (I) indices
- `ihi` (`in`) - Local max longitude (I) indices
- `julo` (`in`) - Local min latitude (J) indices
- `jhi` (`in`) - Local max latitude (J) indices
- `i1` (`in`) - Local min longitude (I) indices
- `i2` (`in`) - Local max longitude (I) indices
- `ju1` (`in`) - Local min latitude (J) indices
- `j2` (`in`) - Local max latitude (J) indices
- `cross` (`in`) - Logical switch. If CROSS=T then cross-terms will be computed.

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function set_cross_terms!(
    crx::Matrix{Float64},
    cry::Matrix{Float64},
    ua::Matrix{Float64},
    va::Matrix{Float64},
    j1p::Int,
    j2p::Int,
    i1_gl::Int,
    i2_gl::Int,
    ju1_gl::Int,
    j2_gl::Int,
    ilo::Int,
    ihi::Int,
    julo::Int,
    jhi::Int,
    i1::Int,
    i2::Int,
    ju1::Int,
    j2::Int,
    cross::Bool
)::Nothing
    # Grid box indices for lon & lat.
    il::Int64, ij::Int64 = 0, 0

    if !cross
        ua .= 0.0
        va .= 0.0
    else
        for ij in j1p:j2p
            for il in i1:i2-1
                ua[il, ij] = 0.5 * (crx[il, ij] + crx[il+1, ij])
            end
            ua[i2, ij] = 0.5 * (crx[i2, ij] + crx[1, ij])
        end

        for ij in ju1+1:j2-1
            for il in i1:i2
                va[il, ij] = 0.5 * (cry[il, ij] + cry[il, ij+1])
            end
        end

        do_cross_terms_pole_i2d2!(cry, va, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
    end
end

"""
Calculates the vertical mass flux.

## Arguments
- `dbk` (`in`; `[k1:k2]`) - Difference in bi across layer - the dSigma term
- `dps_ctm` (`in`; `[i1:i2, ju1:j2]`) - CTM surface pressure tendency; sum over vertical of dpi calculated from original mass fluxes [hPa]
- `dpi` (`in`; `[i1:i2, ju1:j2, k1:k2]`) - Divergence at a grid point; used to calculate vertical motion [mb]
- `wz` (`out`; `[i1:i2, ju1:j2, k1:k2]`) - Large scale mass flux (per time step tdt) in the vertical direction as diagnosed from the hydrostatic relationship [hPa]
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index
- `k1` (`in`) - Local min altitude (K) index
- `k2` (`in`) - Local max altitude (K) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function calc_vert_mass_flux!(
    dbk::Vector{Float64},
    dps_ctm::Matrix{Float64},
    dpi::Array{Float64,3},
    wz::Array{Float64,3},
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    k1::Int64,
    k2::Int64
)::Nothing
    ik::Int64, ij::Int64, il::Int64 = 0, 0, 0

    # Compute vertical mass flux from mass conservation.

    for ij ∈ ju1:j2
        for il ∈ i1:i2
            wz[il, ij, k1] = dpi[il, ij, k1] - (dbk[k1] * dps_ctm[il, ij])
            wz[il, ij, k2] = 0.0
        end
    end

    for ik ∈ k1+1:k2-1
        for ij ∈ ju1:j2
            for il ∈ i1:i2
                wz[il, ij, ik] = wz[il, ij, ik-1] + dpi[il, ij, ik] - (dbk[ik] * dps_ctm[il, ij])
            end
        end
    end
end

"""
Determines Jn and Js, by looking where Courant number is > 1.

## Arguments
- `jn` (`out`, `[k1:k2]`) - Northward of latitude index = jn; Courant numbers could be > 1
- `js` (`out`, `[k1:k2]`) - Southward of latitude index = js; Courant numbers could be > 1
- `crx` (`in`, `[ilo:ihi, julo:jhi, k1:k2]`) - Courant number in E-W direction
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude index at the edge of the S polar cap
- `j2p` (`in`) - Global latitude index at the edge of the N polar cap
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index
- `k1` (`in`) - Local min altitude (K) index
- `k2` (`in`) - Local max altitude (K) index

## Author
Original code from Shian-Jiann Lin, DAO) \\
John Tannahill, LLNL (jrt@llnl.gov)

## Remarks
We cannot parallelize this subroutine because there is a CYCLE statement within the outer loop.

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function set_jn_js!(
    jn::Vector{Int64},
    js::Vector{Int64},
    crx::Array{Float64,3},
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    j2p::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    k1::Int64,
    k2::Int64
)::Nothing
    il::Int64, ij::Int64, ik::Int64 = 0, 0, 0

    js0::Int64 = (j2_gl + 1) / 2
    jn0::Int64 = j2_gl - js0 + 1

    jst::Int64 = max(ju1, j1p)
    jend::Int64 = min(j2, js0)

    for ik ∈ k1:k2
        js[ik] = j1p
        for ij ∈ jend:-1:jst
            for il ∈ i1:i2
                if abs(crx[il, ij, ik]) > 1.0
                    js[ik] = ij
                    @goto ikloop1
                end
            end
        end
        @label ikloop1
    end

    jst = max(ju1, jn0)
    jend = min(j2, j2p)

    for ik ∈ k1:k2
        jn[ik] = j2p
        for ij ∈ jst:jend
            for il ∈ i1:i2
                if abs(crx[il, ij, ik]) > 1.0
                    jn[ik] = ij
                    @goto ikloop2
                end
            end
        end
        @label ikloop2
    end
end

"""
Calculates the advective cross terms.

## Arguments
- `jn` (`in`) - Northward of latitude index = jn, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `js` (`in`) - Southward of latitude index = js, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `qq1` (`in`) - Species concentration (mixing ratio)
- `qqu` (`in`) - Concentration contribution from E-W advection [mixing ratio]
- `qqv` (`in`) - Concentration contribution from N-S advection [mixing ratio]
- `ua` (`in`) - Average of Courant numbers from il and il+1
- `va` (`in`) - Average of Courant numbers from ij and ij+1
- `j1p` (`in`) - Global latitude indices at the edges of the S polar caps J1P=JU1_GL+1 for a polar cap of 1 latitude band J1P=JU1_GL+2 for a polar cap of 2 latitude bands
- `j2p` (`in`) - Global latitude indices at the edges of the N polar caps J2P=J2_GL-1 for a polar cap of 1 latitude band J2P=J2_GL-2 for a polar cap of 2 latitude bands
- `i2_gl` (`in`) - Global min longitude (I) indices
- `ju1_gl` (`in`) - Global min latitude (J) indices
- `j2_gl` (`in`) - Global max latitude (J) indices
- `ilo` (`in`) - Local min longitude (I) indices
- `ihi` (`in`) - Local max longitude (I) indices
- `julo` (`in`) - Local min latitude (J) indices
- `jhi` (`in`) - Local max latitude (J) indices
- `i1` (`in`) - Local min longitude (I) indices
- `i2` (`in`) - Local max longitude (I) indices
- `ju1` (`in`) - Local min latitude (J) indices
- `j2` (`in`) - Local max latitude (J) indices
- `cross` (`in`) - Logical switch: If CROSS=T then cross-terms are being computed

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function calc_advec_cross_terms!(
    jn::Int64,
    js::Int64,
    qq1::Matrix{Float64},
    qqu::Matrix{Float64},
    qqv::Matrix{Float64},
    ua::Matrix{Float64},
    va::Matrix{Float64},
    j1p::Int64,
    j2p::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    cross::Bool
)::Nothing
    i::Int64, imp::Int64, il::Int64, ij::Int64, iu::Int64 = 0, 0, 0, 0, 0
    jv::Int64, iuw::Int64, iue::Int64 = 0, 0, 0
    ril::Float64, rij::Float64, riu::Float64 = 0.0, 0.0, 0.0
    ru::Float64 = 0.0
    qtmp::Matrix{Float64} = zeros(-i2/3:i2+i2/3, julo:jhi)

    for ij ∈ julo:jhi
        qtmp[1:i2, ij] .= qq1[1:i2, ij]
        qtmp[-i2/3:0, ij] = qq1[(-i2/3:0).+i2, ij]
        qtmp[i2+1:i2+i2/3, ij] = qq1[(i2+1:i2+i2/3).-i2, ij]
    end

    if !CROSS
        qqu[:, :] .= qqv[:, :] .= qq1[:, :]
    else
        qqu[:, :] .= qqv[:, :] .= 0

        for ij ∈ j1p:j2p
            if ij <= js || ij >= jn
                # In Polar area, so need to deal with large courant numbers.
                for il ∈ i1:i2
                    iu = ua[il:ij]
                    riu = iu
                    ru = ua[il, ij] - riu
                    iu = il - iu

                    if ua[il, ij] >= 0
                        qqu[il, ij] =
                            qtmp[iu, ij] +
                            ru * (qtmp[iu-1, ij] - qtmp[iu, ij])
                    else
                        qqu[il, ij] =
                            qtmp[iu, ij] +
                            ru * (qtmp[iu, ij] - qtmp[iu+1, ij])
                    end

                    qqu[il, ij] -= qtmp[il, ij]
                end
            else # js < ij < jn
                # Do interior area (use PPM).

                for il ∈ i2:i2
                    ril = il
                    iu = ril - ua[il, ij]

                    qqu[il, ij] =
                        ua[il, ij] *
                        (qtmp[iu, ij] - qtmp[iu+1, ij])
                end
            end

            for il ∈ i1:i2
                rij = ij
                jv = rij - va[il, ij]

                qqv[il, ij] =
                    va[il, ij] *
                    (qtmp[il, jv] - qtmp[il, jv+1])
            end
        end

        for ij ∈ ju1:j2
            for il ∈ i1:ij
                qqu[il, ij] = qtmp[il, ij] + (0.5 * qqu[il, ij])
                qqv[il, ij] = qtmp[il, ij] + (0.5 * qqv[il, ij])
            end
        end
    end
end

const QCKXYZ_FILL_DIAG = false

"""
Check for "filling".
    
## Arguments
- `dq1` (`inout`; `[ilo:ihi]`)- Species density [hPa]
- `j1p` (`in`) - Global latitude indices at the edges of the S polar caps; J1P=JU1_GL+1 for a polar cap of 1 latitude band; J1P=JU1_GL+2 for a polar cap of 2 latitude bands
- `j2p` (`in`) - Global latitude indices at the edges of the N polar caps; J2P=J2_GL-1 for a polar cap of 1 latitude band; J2P=J2_GL-2 for a polar cap of 2 latitude bands
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index
- `k1` (`in`) - Local min altitude (K) index
- `k2` (`in`) - Local max altitude (K) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function qckxyz!(
    dq1::Array{Float64,3},
    j1p::Int64,
    j2p::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    k1::Int64,
    k2::Int64
)::Nothing
    il::Int64, ij::Int64, ik::Int64 = 0, 0, 0
    dup::Float64, qup::Float64 = 0.0, 0.0
    qly::Float64 = 0.0

    ip::Int64 = 0

    # Top layer.

    k1p1::Int64 = k1 + 1

    for ij ∈ j1p:j2p
        for il ∈ i1:i2
            if dq1[il, ij, k1] < 0.0
                ip += 1
                dq1[il, ij, k1p1] += dq1[il, ij, k1]
                dq1[il, ij, k1] = 0.0
            end
        end
    end

    for ik ∈ k1+1:k2-1
        for ij ∈ j1p:j2p
            for il ∈ i1:i2
                if dq1[il, ij, ik] < 0.0
                    ip += 1
                    qup = dq1[il, ij, ik-1]
                    qly = -dq1[il, ij, ik]
                    dup = min(qly, qup)
                    dq1[il, ij, ik-1] = qup - dup
                    dq1[il, ij, ik] = dup - qly
                    dq1[il, ij, ik+1] += dq1[il, ij, ik]
                    dq1[il, ij, ik] = 0.0
                end
            end
        end
    end

    sum::Float64 = 0.0

    k2m1::Int64 = k2 - 1

    # NOTE: Sum seems to be not used in the loop below!
    for ij ∈ j1p:j2p
        for il ∈ i1:i2
            if dq1[il, ij, k2] < 0.0
                ip += 1

                # From above.
                qup = dq1[il, ij, k2m1]
                qly = -dq1[il, ij, k2]
                dup = min(qly, qup)
                dq1[il, ij, k2m1] = qup - dup

                # From "below" the surface.
                sum += qly - dup
                dq1[il, ij, k2] = 0.0
            end
        end
    end

    # We don't want to replace zero values by 1e-30. (ccc, 11/20/08)
    # for k ∈ dq1[i1:i2, j1p:j2p, :]
    #     if k < 1.0
    #         dq1[i1:i2, j1p:j2p, k] = 1.0
    #     end
    # end
end

"""
Sets ILMT, JLMT, KLMT.
    
## Arguments
- `ilmt` (`out`) - Controls various options in E-W advection
- `jlmt` (`out`) - Controls various options in N-S advection
- `klmt` (`out`) - Controls various options in vertical advection
- `i2_gl` (`in`) - Global maximum longitude (I) indices
- `j2_gl` (`in`) - Global maximum longitude (J) indices
- `iord` (`out`) - Flags to denote E-W transport schemes (See REMARKS section of routine Tpcore_FvDas for more info)
- `jord` (`out`) - Flags to denote N-S transport schemes (See REMARKS section of routine Tpcore_FvDas for more info)
- `kord` (`out`) - Flags to denote vertical transport schemes (See REMARKS section of routine Tpcore_FvDas for more info)

## Author
Original code from Shian-Jiann Lin, DAO) \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function set_lmts!(
    ilmt::Ref{Int64},
    jlmt::Ref{Int64},
    klmt::Ref{Int64},
    i2_gl::Int64,
    j2_gl::Int64,
    iord::Ref{Int64},
    jord::Ref{Int64},
    kord::Ref{Int64}
)::Nothing
    j2_glm1::Int64 = j2_gl - 1

    if iord[] <= 0
        if i2_gl >= 144
            ilmt[] = 0
        elseif i2_gl >= 72
            ilmt[] = 1
        else
            ilmt[] = 2
        end
    else
        ilmt[] = iord[] - 3
    end

    if jord[] <= 0
        if j2_glm1 >= 90
            jlmt[] = 0
        elseif j2_glm1 >= 45
            jlmt[] = 1
        else
            jlmt[] = 2
        end
    else
        jlmt[] = jord[] - 3
    end

    klmt[] = max(kord[] - 3, 0)
end

"""
Sets the pressure terms: DELP1, DELPM, PU.

## Arguments
- `dap` (`in`) - Pressure difference across layer from (ai * pt) term [hPa]
- `dbk` (`in`) - Difference in bi across layer - the dSigma term
- `pres1` (`in`; `[ilo:ihi, julo:jhi]`) - Surface pressure at t1 [hPa]
- `pres2` (`in`; `[ilo:ihi, julo:jhi]`) - Surface pressure at t1+tdt [hPa]
- `delp1` (`out`; `[ilo:ihi, julo:jhi]`) - Pressure thickness, the pseudo-density in a hydrostatic system at t1 [hPa]
- `delpm` (`out`; `[ilo:ihi, julo:jhi]`) - Pressure thickness, the pseudo-density in a hydrostatic system at t1+tdt/2 (approximate) [hPa]
- `pu` (`out`; `[ilo:ihi, julo:jhi]`) - Pressure at edges in "u" [hPa]
- `ju1_gl` (`in`) - Global latitude indices at the edges of the S/N polar caps
- `j2_gl` (`in`) - Global latitude indices at the edges of the S/N polar caps
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `j1p` (`in`) - Global latitude indices at the edges of the S polar caps; J1P=JU1_GL+1 for a polar cap of 1 latitude band
- `j2p` (`in`) - Global latitude indices at the edges of the N polar caps; J2P=J2_GL-2 for a polar cap of 2 latitude bands
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO) \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function set_press_terms!(
    dap::Float64,
    dbk::Float64,
    pres1::Matrix{Float64},
    pres2::Matrix{Float64},
    delp1::Matrix{Float64},
    delpm::Matrix{Float64},
    pu::Matrix{Float64},
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    j1p::Int64,
    j2p::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64, ij::Int64 = 0, 0

    delp1 .= dap .+ (dbk .* pres1)
    delpm .= dap .+ (dbk .* 0.5 .* (pres1 .+ pres2))

    for ij ∈ j1p:j2p
        pu[1, ij] = 0.5 .* (delpm[1, ij] + delpm[i2, ij])
        for il ∈ i1+1:i2
            pu[il, ij] = 0.5 .* (delpm[il, ij] + delpm[il-1, ij])
        end
    end
end

"""
Calculates courant numbers from the horizontal mass fluxes.

## Arguments
- `tpcore_fvdas_mod_cose` (`in`; `[ju1_gl:j2_gl]`) - Cosine of grid box edges
- `delpm` (`in`; `[ilo:ihi, julo:jhi]`) - Pressure thickness, the pseudo-density in a hydrostatic system at t1+tdt/2 (approximate) (mb)
- `pu` (`in`; `[ilo:ihi, julo:jhi]`) - Pressure at edges in "u"  (mb)
- `xmass` (`in`; `[ilo:ihi, julo:jhi]`) - Horizontal mass flux in E-W direction [hPa]
- `ymass` (`in`; `[ilo:ihi, julo:jhi]`) - Horizontal mass flux in N-S direction [hPa]
- `crx` (`out`; `[ilo:ihi, julo:jhi]`) - Courant number in E-W direction
- `cry` (`out`; `[ilo:ihi, julo:jhi]`) - Courant number in N-S direction
- `j1p` (`in`) - Global latitude index at the edge of the S polar cap
- `j2p` (`in`) - Global latitude index at the edge of the N polar cap
- `ju1_gl` (`in`) - Global minimum latitude index
- `j2_gl` (`in`) - Global maximum latitude index
- `ilo` (`in`) - Local minimum longitude index
- `ihi` (`in`) - Local maximum longitude index
- `julo` (`in`) - Local minimum latitude index
- `jhi` (`in`) - Local maximum latitude index
- `i1` (`in`) - Local minimum longitude index
- `i2` (`in`) - Local maximum longitude index
- `ju1` (`in`) - Local minimum latitude index
- `j2` (`in`) - Local maximum latitude index

## Author
Original code from Shian-Jiann Lin, DAO) \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function calc_courant!(
    tpcore_fvdas_mod_cose::Vector{Float64},
    delpm::Matrix{Float64},
    pu::Matrix{Float64},
    xmass::Matrix{Float64},
    ymass::Matrix{Float64},
    crx::Matrix{Float64},
    cry::Matrix{Float64},
    j1p::Int64,
    j2p::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    ij::Int64 = 0

    crx .= 0.0
    cry .= 0.0

    for ij = j1p:j2p
        crx[:, ij] .= xmass[:, ij] ./ pu[:, ij]
        cry[:, ij] .= ymass[:, ij] ./ ((0.5 * tpcore_fvdas_mod_cose[ij]) .* (delpm[:, ij] .+ delpm[:, ij-1]))
    end

    cry[:, j2p+1] .= ymass[:, j2p+1] ./ ((0.5 * tpcore_fvdas_mod_cose[j2p+1]) .* (delpm[:, j2p+1] .+ delpm[:, j2p]))
end

"""
Calculate the divergence.

## Arguments
- `do_reduction` (`in`) - Set to F if called on root core or T if called by secondary cores (NOTE: This is only for MPI parallelization, for OPENMP it should be F)
- `geofac_pc` (`in`) - Special geometrical factor (geofac) for Polar cap
- `geofac` (`in`; `[ju1_gl:j2_gl]`) - Geometrical factor for meridional advection; geofac uses correct spherical geometry, and replaces tpcore_fvdas_mod_acosp as the meridional geometrical factor in TPCORE
- `dpi` (`out`; `[i1:i2, ju1:j2]`) - Divergence at a grid point; used to calculate vertical motion [hPa]
- `xmass` (`in`; `[ilo:ihi, julo:jhi]`) - Horizontal mass flux in E/W and N/S directions [hPa]
- `ymass` (`in`; `[ilo:ihi, julo:jhi]`) - Horizontal mass flux in E/W and N/S directions [hPa]
- `j1p` (`in`) - Global latitude indices at the edges of the S polar caps J1P=JU1_GL+1 for a polar cap of 1 latitude band J1P=JU1_GL-1 for a polar cap of 1 latitude bands
- `j2p` (`in`) - Global latitude indices at the edges of the N polar caps J1P=JU1_GL+2 for a polar cap of 1 latitude band J1P=JU1_GL-2 for a polar cap of 2 latitude bands
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function calc_divergence!(
    do_reduction::Bool,
    geofac_pc::Float64,
    geofac::Vector{Float64},
    dpi::Matrix{Float64},
    xmass::Matrix{Float64},
    ymass::Matrix{Float64},
    j1p::Int64,
    j2p::Int64,
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64, ij::Int64 = 0, 0

    # Calculate N-S divergence.

    for ij ∈ j1p:j2p
        dpi[:, ij] .= (ymass[:, ij] .- ymass[:, ij+1]) .* geofac[ij]

        # Calculate E-W divergence.

        for il ∈ i1:i2-1
            dpi[il, ij] = dpi[il, ij] + xmass[il, ij] - xmass[il+1, ij]
        end

        dpi[i2, ij] = dpi[i2, ij] + xmass[i2, ij] - xmass[1, ij]
    end

    do_divergence_pole_sum!(do_reduction, geofac_pc, dpi, ymass, i1_gl, i2_gl, j1p, j2p, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

    if j1p != ju1_gl + 1
        # Polar cap enlarged: copy dpi to polar ring.

        dpi[:, ju1+1] .= dpi[:, ju1]
        dpi[:, j2-1] .= dpi[:, j2]
    end
end

"""
Sets the divergence at the Poles.

## Arguments
- `do_reduction` (`in`) - Set to T if called on root core or F if called by secondary cores
- `geofac_pc` (`in`) - Special geometrical factor (geofac) for Polar cap
- `dpi` (`inout`; `[i1:i2, ju1:j2]`) - Divergence at a grid point; used to calculate vertical motion [hPa]
- `ymass` (`in`; `[ilo:ihi, julo:jhi]`) - Horizontal mass flux in N-S direction [hPa]
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `j1p` (`in`) - Global latitude index at the edge of the S polar cap
- `j2p` (`in`) - Global latitude index at the edge of the N polar cap
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_divergence_pole_sum!(
    do_reduction::Bool,
    geofac_pc::Float64,
    dpi::Matrix{Float64},
    ymass::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    j1p::Int64,
    j2p::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64

    mean_np::Float64 = 0.0
    mean_sp::Float64 = 0.0
    sumnp::Float64 = 0.0
    sumsp::Float64 = 0.0

    ri2::Float64 = i2_gl

    if ju1 == ju1_gl
        sumsp = 0.0

        for il ∈ i1:i2
            sumsp += ymass[il, j1p]
        end

        mean_sp = -sumsp / ri2 * geofac_pc

        for il = i1:i2
            dpi[il, ju1] = mean_sp
        end
    end

    if j2 == j2_gl
        sumnp = 0.0

        for il ∈ i1:i2
            sumnp += ymass[il, j2p+1]
        end

        mean_np = sumnp / ri2 * geofac_pc

        for il ∈ i1:i2
            dpi[il, j2] = mean_np
        end
    end
end

"""
Sets "va" at the Poles.

## Arguments
- `cry` (`in`; `[ilo:ihi, julo:jhi]`) - Courant number in N-S direction
- `va` (`out`; `[ilo:ihi, julo:jhi]`) - Average of Courant numbers from ij and ij+1
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude index at the edge of the South polar cap
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_cross_terms_pole_i2d2!(
    cry::Matrix{Float64},
    va::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64 = 0

    i2d2::Int64 = i2_gl / 2

    if j1p == ju1_gl + 1
        # Polar Cap NOT Enlarged: Get cross terms for N-S horizontal advection.
        if ju1 == ju1_gl
            for il ∈ i1:i2d2
                va[il, ju1] = 0.5 * (cry[il, ju1+1] - cry[il+i2d2, ju1+1])
                va[il+i2d2, ju1] = -va[il, ju1]s
            end
        end

        if j2 == j2_gl
            for il ∈ i1:i2d2
                va[il, j2] = 0.5 * (cry[il, j2] - cry[il+i2d2, j2-1])
                va[il+i2d2, j2] = -va[il, j2]
            end
        end
    end
end

"""
The advective form E-W operator for computing the adx (E-W) cross term.

## Arguments
- `iad` (`in`) - if iad = 1, use 1st order accurate scheme; if iad = 2, use 2nd order accurate scheme
- `jn` (`in`) - Northward of latitude index = jn, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `js` (`in`) - southward of latitude index = js, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `adx` (`out`; `[ilo:ihi, `julo:jhi`]`) - Cross term due to E-W advection [mixing ratio]
- `qqv` (`in`; `[ilo:ihi, `julo:jhi`]`) - Concentration contribution from N-S advection [mixing ratio]
- `ua` (`in`; `[ilo:ihi, `julo:jhi`]`) - Average of Courant numbers from il and il+1
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude index at the edge of the S polar cap
- `j2p` (`in`) - Global latitude index at the edge of the N polar cap
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function xadv_dao2!(
    iad::Int64,
    jn::Int64,
    js::Int64,
    adx::Matrix{Float64},
    qqv::Matrix{Float64},
    ua::Matrix{Float64},
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    j2p::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64, ij::Int64, iu::Int64 = 0, 0, 0
    imp::Int64, iue::Int64, iuw::Int64 = 0, 0, 0
    a1::Float64, b1::Float64, c1::Float64 = 0.0, 0.0, 0.0
    rdiff::Float64 = 0.0
    ril::Float64, riu::Float64 = 0.0, 0.0
    ru::Float64 = 0.0

    qtmp::Matrix{Float64} = zeros(-i2/3:i2+i2/3, julo:jhi)

    adx .= 0.0

    for ij ∈ julo:jhi
        for il = 1:i2
            qtmp[il, ij] = qqv[il, ij]
        end
        for il = -i2/3:0
            qtmp[il, ij] = qqv[i2+il, ij]
        end
        for il = i2+1:i2+i2/3
            qtmp[il, ij] = qqv[il-i2, ij]
        end
    end

    if iad == 1
        for ij ∈ j1p:j2p
            if ij <= js || ij >= jn
                for il ∈ i1:i2
                    iu = ua[il, ij]
                    riu = iu
                    ru = ua[il, ij] - riu
                    iu = il - iu
                    if ua[il, ij] >= 0.0
                        rdiff = qtmp[iu-1, ij] - qtmp[iu, ij]
                    else
                        rdiff = qtmp[iu, ij] - qtmp[iu+1, ij]
                    end
                    adx[il, ij] = (qtmp[iu, ij] - qtmp[il, ij]) + (ru * rdiff)
                end
            else # js < ij < jn
                for il ∈ i1:i2
                    ril = il
                    iu = ril - ua[il, ij]
                    adx[il, ij] = ua[il, ij] * (qtmp[iu, ij] - qtmp[iu+1, ij])
                end
            end
        end
    elseif iad == 2
        for ij ∈ j1p:j2p
            if ij <= js || ij >= jn
                for il ∈ i1:i2
                    iu = Nint(ua[il, ij])
                    riu = iu
                    ru = riu - ua[il, ij]
                    iu = il - iu
                    a1 = 0.5 * (qtmp[iu+1, ij] + qtmp[iu-1, ij]) - qtmp[iu, ij]
                    b1 = 0.5 * (qtmp[iu+1, ij] - qtmp[iu-1, ij])
                    c1 = qtmp[iu, ij] - qtmp[il, ij]
                    adx[il, ij] = (ru * ((a1 * ru) + b1)) + c1
                end
            else # js < ij < jn
                for il ∈ i1:i2
                    iu = Nint(ua[il, ij])
                    riu = iu
                    ru = riu - ua[il, ij]
                    iu = il - iu
                    a1 = 0.5 * (qtmp[iu+1, ij] + qtmp[iu-1, ij]) - qtmp[iu, ij]
                    b1 = 0.5 * (qtmp[iu+1, ij] - qtmp[iu-1, ij])
                    c1 = qtmp[iu, ij] - qtmp[il, ij]
                    adx[il, ij] = (ru * ((a1 * ru) + b1)) + c1
                end
            end
        end
    end

    if ju1 == ju1_gl
        adx[i1:i2, ju1] .= 0.0
        if j1p != ju1_gl + 1
            adx[i1:i2, ju1+1] .= 0.0
        end
    end

    if j2 == j2_gl
        adx[i1:i2, j2] .= 0.0
        if j1p != ju1_gl + 1
            adx[i1:i2, j2-1] .= 0.0
        end
    end
end

"""
The advective form N-S operator for computing the ady (N-S) cross term.

## Arguments
- `iad` (`in`) - If iad = 1, use 1st order accurate scheme; If iad = 2, use 2nd order accurate scheme
- `ady` (`out`; `[ilo:ihi, julo:jhi]`) - Cross term due to N-S advection (mixing ratio)
- `qqu` (`in`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from E-W advection [mixing ratio]
- `va` (`in`; `[ilo:ihi, julo:jhi]`) - Average of Courant numbers from ij and ij+1
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude index at the edge of the S polar cap
- `j2p` (`in`) - Global latitude index at the edge of the N polar cap
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function yadv_dao2!(
    iad::Int64,
    ady::Matrix{Float64},
    qqu::Matrix{Float64},
    va::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    j2p::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64, ij::Int64 = 0, 0
    jv::Int64 = 0
    a1::Float64, b1::Float64, c1::Float64 = 0.0, 0.0, 0.0
    rij::Float64, rjv::Float64 = 0.0, 0.0
    rv::Float64 = 0.0

    # We may need a small ghost zone depending on the polar cap used
    qquwk::Matrix{Float64} = zeros(ilo:ihi, julo-2:jhi+2)

    # Zero output array
    ady .= 0.0

    # Make work array
    qquwk[:, julo:jhi] .= qqu[:, julo:jhi]

    # This routine creates a ghost zone in latitude in case of not enlarged polar cap (ccc, 11/20/08)
    do_yadv_pole_i2d2!(qqu, qquwk, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2,)

    if iad == 1
        # 1st order.
        for ij in j1p-1:j2p+1
            for il in i1:i2
                #c?
                rij = ij
                jv = rij - va[il, ij]

                ady[il, ij] = va[il, ij] * (qquwk[il, jv] - qquwk[il, jv+1])
            end
        end
    elseif iad == 2
        for ij in j1p-1:j2p+1
            for il in i1:i2
                #c?
                jv = round(Int, va[il, ij])
                rjv = jv
                rv = rjv - va[il, ij]
                jv = ij - jv

                a1 = 0.5 * (qquwk[il, jv+1] + qquwk[il, jv-1]) - qquwk[il, jv]

                b1 = 0.5 * (qquwk[il, jv+1] - qquwk[il, jv-1])

                c1 = qquwk[il, jv] - qquwk[il, ij]

                ady[il, ij] = (rv * ((a1 * rv) + b1)) + c1
            end
        end
    end

    do_yadv_pole_sum!(ady, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2,)
end

"""
Sets "qquwk" at the Poles.

## Arguments
- `qqu` (`in`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from E-W advection [mixing ratio]
- `qquwk` (`inout`; `[ilo:ihi, julo-2:jhi+2]`) - qqu working array [mixing ratio]
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (J) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude indices at the edges of the South polar cap; J1P=JU1_GL+1 for a polar cap of 1 latitude band; J1P=JU1_GL+2 for a polar cap of 2 latitude bands

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_yadv_pole_i2d2!(
    qqu::Matrix{Float64},
    qquwk::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    qqu_copy::Matrix{Float64} = copy(qqu)

    il::Int64, ij::Int64 = 0, 0
    inb::Int64 = 0

    i2d2::Int64 = i2_gl / 2

    if j1p == ju1_gl + 1
        # Polar Cap NOT Enlarged.
        if ju1 == ju1_gl
            for il ∈ i1:i2d2
                for inb ∈ 1:2
                    qquwk[il, ju1-inb] = qqu_copy[il+i2d2, ju1+inb]
                    qquwk[il+i2d2, ju1-inb] = qqu_copy[il, ju1+inb]
                end
            end
        end

        if j2 == j2_gl
            for il ∈ i1:i2d2
                for inb ∈ 1:2
                    qquwk[il, j2+inb] = qqu_copy[il+i2d2, j2-inb]
                    qquwk[il+i2d2, j2+inb] = qqu_copy[il, j2-inb]
                end
            end
        end
    end
end

"""
Sets the cross term due to N-S advection at the Poles.

## Arguments
- `ady` (`inout`; `[ilo:ihi, julo:jhi]`) - Cross term due to N-S advection (mixing ratio)
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude index at the edge of the South polar cap
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_yadv_pole_sum!(
    ady::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64 = 0

    # Test if we are using extended polar caps (i.e. the S pole and next N latitude and N. Pole and next S latitude). Do this outside the loops. (bmy, 12/11/08)
    IS_EXT_POLAR_CAP::Bool = (j1p == ju1_gl + 2)

    # South Pole

    sumsp::Float64 = 0.0
    sumnp::Float64 = 0.0

    if IS_EXT_POLAR_CAP
        # For a 2-latitude polar cap (S. Pole + next Northward latitude)
        for il = i1:i2
            sumsp += ady[il, ju1+1]
            sumnp += ady[il, j2-1]
        end
    else
        # For a 1-latitude polar cap (S. Pole only)
        for il = i1:i2
            sumsp += ady[il, ju1]
            sumnp += ady[il, j2]
        end
    end

    sumsp /= i2_gl
    sumnp /= i2_gl

    if IS_EXT_POLAR_CAP
        # For a 2-latitude polar cap (S. Pole + next Northward latitude)
        for il = i1:i2
            ady[il, ju1+1] = sumsp
            ady[il, ju1] = sumsp
            ady[il, j2-1] = sumnp
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
- `ilmt` (`in`) - Controls various options in E-W advection
- `jn` (`in`) - Northward of latitude index = jn, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `js` (`in`) - Southward of latitude index = js, Courant numbers could be > 1, so use the flux-form semi-Lagrangian scheme
- `pu` (`in`; `[ilo:ihi, julo:jhi]`) - pressure at edges in "u" [hPa]
- `crx` (`in`; `[ilo:ihi, julo:jhi]`) - Courant number in E-W direction
- `dq1` (`inout`; `[ilo:ihi, julo:jhi]`) - Species density [hPa]
- `qqv` (`inout`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from N-S advection [mixing ratio]
- `xmass` (`in`; `[ilo:ihi, julo:jhi]`) - Horizontal mass flux in E-W direction [hPa]
- `fx` (`out`; `[ilo:ihi, julo:jhi]`) - E-W flux [mixing ratio]
- `j1p` (`in`) - Global latitude indices at the edges of the S/N polar caps
- `j2p` (`in`) - Global latitude indices at the edges of the S/N polar caps

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function xtp!(
    ilmt::Int64,
    jn::Int64,
    js::Int64,
    pu::Matrix{Float64},
    crx::Matrix{Float64},
    dq1::Matrix{Float64},
    qqv::Matrix{Float64},
    xmass::Matrix{Float64},
    fx::Matrix{Float64},
    j1p::Int64,
    j2p::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    iord::Int64
)::Nothing
    pu_copy::Matrix{Float64} = copy(pu)
    crx_copy::Matrix{Float64} = copy(crx)
    xmass_copy::Matrix{Float64} = copy(xmass)

    il::Int64, ij::Int64, ic::Int64 = 0, 0, 0
    iu::Int64, ix::Int64, iuw::Int64, iue::Int64, imp::Int64 = 0, 0, 0, 0, 0
    rc::Float64 = 0.0
    ric::Float64, ril::Float64 = 0.0, 0.0
    isav::Vector{Int64} = zeros(i1:i2)
    dcx::Matrix{Float64} = zeros(-i2/3:i2+i2/3, julo:jhi)
    qtmp::Matrix{Float64} = zeros(-i2/3:i2+i2/3, julo:jhi)

    fx .= 0.0

    imp = i2 + 1

    # NOTE: these loops do not parallelize well (bmy, 12/5/08)

    # Populate qtmp
    for il ∈ i1:i2
        qtmp[il, :] .= qqv[il, :]
    end

    for il ∈ -i2/3:0
        qtmp[il, :] .= qqv[i2+il, :]
    end

    for il ∈ i2+1:i2+i2/3
        qtmp[il, :] .= qqv[il-i2, :]
    end

    if iord != 1
        qtmp[i1-1, :] .= qqv[i2, :]
        qtmp[i1-2, :] .= qqv[i2-1, :]
        qtmp[i2+1, :] .= qqv[i1, :]
        qtmp[i2+2, :] .= qqv[i1+1, :]

        xmist!(dcx, qtmp, j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
    end

    jvan::Int64 = max(1, j2_gl / 18)

    for ij ∈ j1p:j2p
        if ij > js && ij < jn
            # Do horizontal Eulerian advection in the E-W direction.

            if iord == 1 || ij == j1p || ij == j2p
                for il ∈ i1:i2
                    ril = il
                    iu = ril - crx_copy[il, ij]
                    fx[il, ij] = qtmp[iu, ij]
                end
            else
                if iord == 2 || ij <= j1p + jvan || ij >= j2p - jvan
                    for il ∈ i1:i2
                        ril = il
                        iu = ril - crx_copy[il, ij]
                        fx[il, ij] = qtmp[iu, ij] + dcx[iu, ij] * (sign(crx_copy[il, ij]) - crx_copy[il, ij])
                    end
                else
                    fxppm!(ij, ilmt, crx_copy, dcx, fx, qtmp, -i2 / 3, i2 + i2 / 3, julo, jhi, i1, i2)
                end
            end

            for il ∈ i1:i2
                fx[il, ij] *= xmass_copy[il, ij]
            end
        else
            # Do horizontal Conservative (flux-form) Semi-Lagrangian advection in the E-W direction (van Leer at high latitudes).

            if iord == 1 || ij == j1p || ij == j2p
                for il ∈ i1:i2
                    ic = crx_copy[il, ij]
                    isav[il] = il - ic
                    ril = il
                    iu = ril - crx_copy[il, ij]
                    ric = ic
                    rc = crx_copy[il, ij] - ric
                    fx[il, ij] = rc * qtmp[iu, ij]
                end
            else
                for il ∈ i1:i2
                    ic = crx_copy[il, ij]
                    isav[il] = il - ic
                    ril = il
                    iu = ril - crx_copy[il, ij]
                    ric = ic
                    rc = crx_copy[il, ij] - ric
                    fx[il, ij] = rc * (qtmp[iu, ij] + dcx[iu, ij] * (sign(rc) - rc))
                end
            end

            for il ∈ i1:i2
                if crx_copy[il, ij] > 1.0
                    for ix ∈ isav[il]:il-1
                        fx[il, ij] += qtmp[ix, ij]
                    end
                elseif crx_copy[il, ij] < -1.0
                    for ix ∈ il:isav[il]-1
                        fx[il, ij] -= qtmp[ix, ij]
                    end
                end
            end

            for il ∈ i1:i2
                fx[il, ij] *= pu_copy[il, ij]
            end
        end
    end

    # NOTE: This loop does not parallelize well (bmy, 12/5/08)
    for ij ∈ j1p:j2p
        for il ∈ i1:i2-1
            dq1[il, ij] += fx[il, ij] - fx[il+1, ij]
        end
        dq1[i2, ij] += fx[i2, ij] - fx[i1, ij]
    end
end

"""
Computes the linear tracer slope in the E-W direction. It uses the Lin et. al. 1994 algorithm.

## Arguments
- `dcx` (`out`; `[-i2/3:i2/3, julo:jhi]`) - Slope of concentration distribution in E-W direction [mixing ratio]
- `qqv` (`in`; `[-i2/3:i2/3, julo:jhi]`) - Concentration contribution from N-S advection [mixing ratio]
- `j1p` (`in`) - Global latitude indices at the edges of the S polar caps; J1P=JU1_GL+1 for a polar cap of 1 latitude band; J1P=JU1_GL+2 for a polar cap of 2 latitude bands
- `j2p` (`in`) - Global latitude indices at the edges of the N polar caps; J2P=J2_GL-1 for a polar cap of 1 latitude band; J2P=J2_GL-2 for a polar cap of 2 latitude bands
- `i2_gl` (`in`) - Global min longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function xmist!(
    dcx::Matrix{Float64},
    qqv::Matrix{Float64},
    j1p::Int64,
    j2p::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    qqv_copy::Matrix{Float64} = copy(qqv)

    il::Int64, ij::Int64 = 0, 0
    pmax::Float64, pmin::Float64 = 0.0, 0.0
    tmp::Float64 = 0.0

    r24::Float64 = 1.0 / 24.0

    for ij ∈ j1p+1:j2p-1
        for il ∈ i1:i2
            tmp = (
                (8.0 * (qqv_copy[il+1, ij] - qqv_copy[il-1, ij])) +
                qqv_copy[il-2, ij] - qqv_copy[il+2, ij]
            ) * r24

            pmax = max(
                qqv_copy[il-1, ij], qqv_copy[il, ij], qqv_copy[il+1, ij]
            ) - qqv_copy[il, ij]

            pmin = qqv_copy[il, ij] - min(
                qqv_copy[il-1, ij], qqv_copy[il, ij], qqv_copy[il+1, ij]
            )

            dcx[il, ij] = sign(tmp) * min(abs(tmp), pmax, pmin)
        end
    end

    # Populate ghost zones of dcx (ccc, 11/20/08)

    for ij ∈ julo:jhi
        for il ∈ -i2_gl/3:0
            dcx[il, ij] = dcx[i2_gl+il, ij]
        end

        for il ∈ i2_gl+1:i2_gl+i2_gl/3
            dcx[il, ij] = dcx[il-i2_gl, ij]
        end
    end
end

"""
The 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the E-W direction.

## Arguments
- `ij` (`in`) - Latitude (IJ) and altitude (IK) indices
- `ilmt` (`in`) - Controls various options in E-W advection
- `crx` (`in`; `[i1:i2, julo:jhi]`) - Courant number in E-W direction
- `dcx` (`inout`; `[ilo:ihi, julo:jhi]`) - Slope of concentration distribution in E-W direction (mixing ratio)
- `fx` (`out`; `[i1:i2, julo:jhi]`) - E-W flux [mixing ratio]
- `qqv` (`inout`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from N-S advection [mixing ratio]
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Remarks
This routine is called from w/in a OpenMP parallel loop fro

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function fxppm!(
    ij::Int64,
    ilmt::Int64,
    crx::Matrix{Float64},
    dcx::Matrix{Float64},
    fx::Matrix{Float64},
    qqv::Matrix{Float64},
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64
)::Nothing
    crx_copy::Matrix{Float64} = copy(crx)
    dcx_copy::Matrix{Float64} = copy(dcx)

    il::Int64 = 0
    ilm1::Int64 = 0
    lenx::Int64 = 0
    rval::Float64 = 0.0

    a6::Vector{Float64} = zeros(ilo:ihi)
    al::Vector{Float64} = zeros(ilo:ihi)
    ar::Vector{Float64} = zeros(ilo:ihi)
    a61::Vector{Float64} = zeros((ihi - 1) - (ilo + 1) + 1)
    al1::Vector{Float64} = zeros((ihi - 1) - (ilo + 1) + 1)
    ar1::Vector{Float64} = zeros((ihi - 1) - (ilo + 1) + 1)
    dcxi1::Vector{Float64} = zeros((ihi - 1) - (ilo + 1) + 1)
    qqvi1::Vector{Float64} = zeros((ihi - 1) - (ilo + 1) + 1)

    r13::Float64 = 1.0 / 3.0
    r23::Float64 = 2.0 / 3.0

    for il ∈ ilo+1:ihi
        rval = 0.5 * (qqv[il-1, ij] + qqv[il, ij]) + (dcx_copy[il-1, ij] - dcx_copy[il, ij]) * r13

        al[il] = rval
        ar[il-1] = rval
    end

    for il ∈ ilo+1:ihi-1
        a6[il] = 3.0 * (qqv[il, ij] + qqv[il, ij] - (al[il] + ar[il]))
    end

    if ilmt <= 2
        a61 .= 0.0
        al1 .= 0.0
        ar1 .= 0.0

        dcxi1 .= 0.0
        qqvi1 .= 0.0

        lenx = 0

        for il ∈ ilo+1:ihi-1
            lenx += 1

            a61[lenx] = a6[il]
            al1[lenx] = al[il]
            ar1[lenx] = ar[il]

            dcxi1[lenx] = dcx_copy[il, ij]
            qqvi1[lenx] = qqv[il, ij]
        end

        lmtppm!(lenx, ilmt, a61, al1, ar1, dcxi1, qqvi1)

        lenx = 0

        for il ∈ ilo+1:ihi-1
            lenx += 1

            a6[il] = a61[lenx]
            al[il] = al1[lenx]
            ar[il] = ar1[lenx]

            dcx_copy[il, ij] = dcxi1[lenx]
            qqv[il, ij] = qqvi1[lenx]
        end

        # Populate ghost zones of qqv and dcx with new values (ccc, 11/20/08)
        for il ∈ -i2/3:0
            dcx_copy[il, ij] = dcx_copy[i2+il, ij]
            qqv[il, ij] = qqv[i2+il, ij]
        end

        for il ∈ i2+1:i2+i2/3
            dcx_copy[il, ij] = dcx_copy[il-i2, ij]
            qqv[il, ij] = qqv[il-i2, ij]
        end
    end

    for il ∈ i1+1:i2
        if crx_copy[il, ij] > 0.0
            ilm1 = il - 1

            fx[il, ij] = ar[ilm1] + 0.5 * crx_copy[il, ij] * (al[ilm1] - ar[ilm1] + a6[ilm1] * (1.0 - r23 * crx_copy[il, ij]))
        else
            fx[il, ij] = al[il] - 0.5 * crx_copy[il, ij] * (ar[il] - al[il] + a6[il] * (1.0 + r23 * crx_copy[il, ij]))
        end
    end

    # First box case (ccc, 11/20/08)
    if crx_copy[i1, ij] > 0.0
        ilm1 = i2

        fx[i1, ij] = ar[ilm1] + 0.5 * crx_copy[i1, ij] * (al[ilm1] - ar[ilm1] + a6[ilm1] * (1.0 - r23 * crx_copy[i1, ij]))
    else
        fx[i1, ij] = al[i1] - 0.5 * crx_copy[i1, ij] * (ar[i1] - al[i1] + a6[i1] * (1.0 + r23 * crx_copy[i1, ij]))
    end
end

"""
Enforces the full monotonic, semi-monotonic, or the positive-definite constraint to the sub-grid parabolic distribution of the Piecewise Parabolic Method (PPM).

## Arguments
- `lenx` (`in`) - Vector length
- `lmt` (`in`) - If 0 => full monotonicity; If 1 => semi-monotonic constraint (no undershoots); If 2 => positive-definite constraint
- `a6` (`inout`, `[lenx]`) - Curvature of the test parabola
- `al` (`inout`, `[lenx]`) - Left edge value of the test parabola
- `ar` (`inout`, `[lenx]`) - Right edge value of the test parabola
- `dc` (`inout`, `[lenx]`) - 0.5 * mismatch
- `qa` (`inout`, `[lenx]`) - Cell-averaged value

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function lmtppm!(
    lenx::Int64,
    lmt::Int64,
    a6::Vector{Float64},
    al::Vector{Float64},
    ar::Vector{Float64},
    dc::Vector{Float64},
    qa::Vector{Float64}
)::Nothing
    il::Int64 = 0
    a6da::Float64 = 0.0
    da1::Float64, da2::Float64 = 0.0, 0.0
    fmin::Float64, ftmp::Float64 = 0.0, 0.0
    r12::Float64 = 1.0 / 12.0

    if lmt == 0
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
- `jlmt` (`in`) - Controls various options in N-S advection
- `geofac_pc` (`in`) - special geometrical factor (geofac) for Polar cap
- `geofac` (`in`; `[ju1_gl:j2_gl]`) - geometrical factor for meridional advection; geofac uses correct spherical geometry, and replaces tpcore_fvdas_mod_acosp as the  meridional geometrical factor in tpcore
- `cry` (`in`; `[ilo:ihi, julo:jhi]`) - Courant number in N-S direction
- `dq1` (`inout`; `[ilo:ihi, julo:jhi]`) - Species density [hPa]
- `qqu` (`inout`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from E-W advection [mixing ratio]
- `qqv` (`inout`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from N-S advection [mixing ratio]
- `ymass` (`in`; `[ilo:ihi, julo:jhi]`) - Horizontal mass flux in N-S direction [hPa]
- `fy` (`out`; `[ilo:ihi, julo:jhi+1]`) - N-S flux [mixing ratio]
- `j1p` (`in`) - Global latitude indices at the edges of the S/N polar caps; J1P=JU1_GL+1 for a polar cap of 1 latitude band; J1P=JU1_GL+2 for a polar cap of 2 latitude bands
- `j2p` (`in`) - Global latitude indices at the edges of the S/N polar caps; J2P=J2_GL-1 for a polar cap of 1 latitude band; J2P=J2_GL-2 for a polar cap of 2 latitude bands
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `ilong` (`in`) - ???
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index
- `jord` (`in`) - N-S transport scheme (see module header for more info)

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history

"""
function ytp!(
    jlmt::Int64,
    geofac_pc::Float64,
    geofac::Vector{Float64},
    cry::Matrix{Float64},
    dq1::Matrix{Float64},
    qqu::Matrix{Float64},
    qqv::Matrix{Float64},
    ymass::Matrix{Float64},
    fy::Matrix{Float64},
    j1p::Int64,
    j2p::Int64,
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64
)::Nothing
    geofac_copy::Vector{Float64} = copy(geofac)
    cry_copy::Matrix{Float64} = copy(cry)
    ymass_copy::Matrix{Float64} = copy(ymass)

    il::Int64, ij::Int64 = 0, 0
    jv::Int64 = 0
    rj1p::Float64 = 0.0

    dcy::Matrix{Float64} = zeros(i1_gl:i2_gl, ju1_gl:j2_gl)

    dcy .= 0.0
    fy .= 0.0

    if jord == 1
        for ij ∈ j1p:j2p+1
            for il ∈ i1:i2
                jv = rj1p - cry_copy[il, ij]
                qqv[il, ij] = qqu[il, jv]
            end
        end
    else

        ymist!(4, dcy, qqu, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

        if jord <= 0 || jord >= 3
            fyppm!(jlmt, cry_copy, dcy, qqu, qqv, j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
        else
            for ij ∈ j1p:j2p+1
                for il = i1:i2
                    jv = rj1p - cry_copy[il, ij]
                    qqv[il, ij] = qqu[il, jv] + (sign(cry_copy[il, ij]) * 1.0 - cry_copy[il, ij]) * dcy[il, jv]
                end
            end
        end

    end

    for ij ∈ j1p:j2p+1
        qqv[i1:i2, ij] = qqv[i1:i2, ij] .* ymass_copy[i1:i2, ij]
    end

    #.sds.. save N-S species flux as diagnostic

    for ij ∈ i1:i2
        fy[ij, j1p:j2p+1] = qqv[ij, j1p:j2p+1] .* geofac_copy[j1p:j2p+1]
    end

    #... meridional flux update
    for ij ∈ j1p:j2p
        dq1[i1:i2, ij] .= dq1[i1:i2, ij] .+ (qqv[i1:i2, ij] .- qqv[i1:i2, ij+1]) * geofac_copy[ij]
    end

end

"""
Computes the linear tracer slope in the N-S direction.  It uses the Lin et. al. 1994 algorithm.

## Arguments
- `id` (`in`) - The "order" of the accuracy in the computed linear "slope" (or mismatch, Lin et al. 1994); it is either 2 or 4.
- `dcy` (`out`, `[ilo:ihi, julo:jhi]`) - Slope of concentration distribution in N-S direction [mixing ratio]
- `qqu` (`in`, `[ilo:ihi, julo:jhi]`) - Concentration contribution from E-W advection (mixing ratio)
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude index at the edge of the South polar cap; J1P=JU1_GL+1 for a polar cap of 1 latitude band; J1P=JU1_GL+2 for a polar cap of 2 latitude bands
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function ymist!(
    id::Int64,
    dcy::Matrix{Float64},
    qqu::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    qqu_copy::Matrix{Float64} = copy(qqu)

    il::Int64, ij::Int64 = 0, 0
    pmax::Float64, pmin::Float64 = 0.0, 0.0
    tmp::Float64 = 0.0

    qtmp::Matrix{Float64} = zeros(ilo:ihi, julo-2:jhi+2)

    r24::Float64 = 1.0 / 24.0

    # Populate qtmp
    qtmp .= 0.0
    for ij = ju1:j2
        qtmp[:, ij] .= qqu_copy[:, ij]
    end

    if id == 2
        for ij = ju1-1:j2-1
            for il = i1:i2
                tmp = 0.25 * (qtmp[il, ij+2] - qtmp[il, ij])
                pmax = max(qtmp[il, ij], qtmp[il, ij+1], qtmp[il, ij+2]) - qtmp[il, ij+1]
                pmin = qtmp[il, ij+1] - min(qtmp[il, ij], qtmp[il, ij+1], qtmp[il, ij+2])
                dcy[il, ij+1] = sign(tmp) * min(abs(tmp), pmin, pmax)
            end
        end
    else
        do_y_pole1_i2d2!(dcy, qtmp, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

        for ij = ju1-2:j2-2
            for il = i1:i2
                tmp = ((8.0 * (qtmp[il, ij+3] - qtmp[il, ij+1])) + qtmp[il, ij] - qtmp[il, ij+4]) * r24
                pmax = max(qtmp[il, ij+1], qtmp[il, ij+2], qtmp[il, ij+3]) - qtmp[il, ij+2]
                pmin = qtmp[il, ij+2] - min(qtmp[il, ij+1], qtmp[il, ij+2], qtmp[il, ij+3])
                dcy[il, ij+2] = sign(tmp) * min(abs(tmp), pmin, pmax)
            end
        end
    end

    do_y_pole2_i2d2!(dcy, qtmp, i1_gl, i2_gl, ju1_gl, j2_gl, j1p, ilo, ihi, julo, jhi, i1, i2, ju1, j2)
end

"""
Sets "dcy" at the Poles.

## Arguments
- `dcy` (`out`, `[ilo:ihi, julo:jhi]`) - Slope of concentration distribution in N-S direction [mixing ratio]
- `qqu` (`in`, `[ilo:ihi, julo-2:jhi+2]`) - Concentration contribution from E-W advection [mixing ratio]
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_ymist_pole1_i2d2!(
    dcy::Matrix{Float64},
    qqu::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    qqu_copy::Matrix{Float64} = copy(qqu)

    il::Int64 = 0
    pmax::Float64, pmin::Float64 = 0.0, 0.0
    tmp::Float64 = 0.0

    i2d2::Int64 = i2_gl / 2
    r24::Float64 = 1.0 / 24.0

    if ju1 == ju1_gl
        for il = i1:i2d2
            tmp = ((8.0 * (qqu_copy[il, ju1+2] - qqu_copy[il, ju1])) + qqu_copy[il+i2d2, ju1+1] - qqu_copy[il, ju1+3]) * r24
            pmax = max(qqu_copy[il, ju1], qqu_copy[il, ju1+1], qqu_copy[il, ju1+2]) - qqu_copy[il, ju1+1]
            pmin = qqu_copy[il, ju1+1] - min(qqu_copy[il, ju1], qqu_copy[il, ju1+1], qqu_copy[il, ju1+2])
            dcy[il, ju1+1] = sign(tmp) * min(abs(tmp), pmin, pmax)
        end
        for il = i1+i2d2:i2
            tmp = ((8.0 * (qqu_copy[il, ju1+2] - qqu_copy[il, ju1])) + qqu_copy[il-i2d2, ju1+1] - qqu_copy[il, ju1+3]) * r24
            pmax = max(qqu_copy[il, ju1], qqu_copy[il, ju1+1], qqu_copy[il, ju1+2]) - qqu_copy[il, ju1+1]
            pmin = qqu_copy[il, ju1+1] - min(qqu_copy[il, ju1], qqu_copy[il, ju1+1], qqu_copy[il, ju1+2])
            dcy[il, ju1+1] = sign(tmp) * min(abs(tmp), pmin, pmax)
        end
    end

    if j2 == j2_gl
        for il = i1:i2d2
            tmp = ((8.0 * (qqu_copy[il, j2] - qqu_copy[il, j2-2])) + qqu_copy[il, j2-3] - qqu_copy[il+i2d2, j2-1]) * r24
            pmax = max(qqu_copy[il, j2-2], qqu_copy[il, j2-1], qqu_copy[il, j2]) - qqu_copy[il, j2-1]
            pmin = qqu_copy[il, j2-1] - min(qqu_copy[il, j2-2], qqu_copy[il, j2-1], qqu_copy[il, j2])
            dcy[il, j2-1] = sign(tmp) * min(abs(tmp), pmin, pmax)
        end
        for il = i1+i2d2:i2
            tmp = ((8.0 * (qqu_copy[il, j2] - qqu_copy[il, j2-2])) + qqu_copy[il, j2-3] - qqu_copy[il-i2d2, j2-1]) * r24
            pmax = max(qqu_copy[il, j2-2], qqu_copy[il, j2-1], qqu_copy[il, j2]) - qqu_copy[il, j2-1]
            pmin = qqu_copy[il, j2-1] - min(qqu_copy[il, j2-2], qqu_copy[il, j2-1], qqu_copy[il, j2])
            dcy[il, j2-1] = sign(tmp) * min(abs(tmp), pmin, pmax)
        end
    end
end

"""
Sets "dcy" at the Poles.

## Arguments
- `dcy` (`out`; `[ilo:ihi, julo:jhi]`) - Slope of the concentration distribution in the N-S direction at the poles.
- `qqu` (`in`; `[ilo:ihi, julo-2:jhi+2]`) - Concentration contribution from E-W advection at the poles.
- `i1_gl` (`in`) - Global min longitude (I) index.
- `i2_gl` (`in`) - Global max longitude (I) index.
- `ju1_gl` (`in`) - Global latitude index at the edge of the South polar cap.
- `j2_gl` (`in`) - Global max latitude (J) index.
- `j1p` (`in`) - J1P=JU1_GL+1 for a polar cap of 1 latitude band.
- `i1` (`in`) - Local min longitude (I) index.
- `i2` (`in`) - Local max longitude (I) index.
- `ju1` (`in`) - Local latitude index at the edge of the South polar cap.
- `j2` (`in`) - Local max latitude (J) index.

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_ymist_pole2_i2d2!(
    dcy::Matrix{Float64},
    qqu::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    qqu_copy::Matrix{Float64} = copy(qqu)

    il::Int64 = 0
    pmax::Float64, pmin::Float64 = 0.0, 0.0
    tmp::Float64 = 0.0

    i2d2::Int64 = i2_gl / 2

    if ju1 == ju1_gl
        if j1p != ju1_gl + 1
            dcy[i1:i2, ju1] .= 0.0
        else
            # Determine slope in South Polar cap for scalars.
            for il = i1:i2d2
                tmp = 0.25 * (qqu_copy[il, ju1+1] - qqu_copy[il+i2d2, ju1+1])
                pmax = max(qqu_copy[il, ju1+1], qqu_copy[il, ju1], qqu_copy[il+i2d2, ju1+1]) - qqu_copy[il, ju1]
                pmin = qqu_copy[il, ju1] - min(qqu_copy[il, ju1+1], qqu_copy[il, ju1], qqu_copy[il+i2d2, ju1+1])
                dcy[il, ju1] = sign(tmp) * min(abs(tmp), pmax, pmin)
            end

            for il = i1+i2d2:i2
                dcy[il, ju1] = -dcy[il-i2d2, ju1]
            end
        end
    end

    if j2 == j2_gl
        if j1p != ju1_gl + 1
            dcy[i1:i2, j2] .= 0.0
        else
            # Determine slope in North Polar cap for scalars.
            for il = i1:i2d2
                tmp = 0.25 * (qqu_copy[il+i2d2, j2-1] - qqu_copy[il, j2-1])
                pmax = max(qqu_copy[il+i2d2, j2-1], qqu_copy[il, j2], qqu_copy[il, j2-1]) - qqu_copy[il, j2]
                pmin = qqu_copy[il, j2] - min(qqu_copy[il+i2d2, j2-1], qqu_copy[il, j2], qqu_copy[il, j2-1])
                dcy[il, j2] = sign(tmp) * min(abs(tmp), pmax, pmin)
            end

            for il = i1+i2d2:i2
                dcy[il, j2] = -dcy[il-i2d2, j2]
            end
        end
    end
end

"""
The 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the N-S direction.

## Arguments
- `jlmt` (`in`) - Number of latitude bands in the polar cap.
- `cry` (`in`; `[ilo:ihi, julo:jhi]`) - Flux in the N-S direction.
- `dcy` (`in`; `[ilo:ihi, julo:jhi]`) - Slope of the concentration distribution in the N-S direction.
- `qqu` (`in`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from E-W advection.
- `qqv` (`out`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from N-S advection.
- `j1p` (`in`) - J1P=JU1_GL+1 for a polar cap of 1 latitude band.
- `j2p` (`in`) - J2P=J2_GL-1 for a polar cap of 1 latitude band.
- `i1_gl` (`in`) - Global min longitude (I) index.
- `i2_gl` (`in`) - Global max longitude (I) index.
- `ju1_gl` (`in`) - Global latitude index at the edge of the South polar cap.
- `j2_gl` (`in`) - Global max latitude (J) index.
- `ilong` (`in`) - ???
- `ilo` (`in`) - Local min longitude (I) index.
- `ihi` (`in`) - Local max longitude (I) index.
- `julo` (`in`) - Local latitude index at the edge of the South polar cap.
- `jhi` (`in`) - Local max latitude (J) index.
- `i1` (`in`) - Local min longitude (I) index.
- `i2` (`in`) - Local max longitude (I) index.
- `ju1` (`in`) - Local latitude index at the edge of the South polar cap.
- `j2` (`in`) - Local max latitude (J) index.

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function fyppm!(
    jlmt::Int64,
    cry::Matrix{Float64},
    dcy::Matrix{Float64},
    qqu::Matrix{Float64},
    qqv::Matrix{Float64},
    j1p::Int64,
    j2p::Int64,
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilong::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    cry_copy::Matrix{Float64} = copy(cry)
    dcy_copy::Matrix{Float64} = copy(dcy)
    qqu_copy::Matrix{Float64} = copy(qqu)

    ijm1::Int64 = 0
    il::Int64, ij::Int64 = 0, 0
    lenx::Int64 = 0
    r13::Float64, r23::Float64 = 0.0, 0.0

    a61::Matrix{Float64} = zeros(ilong * ((jhi - 1) - (julo + 1) + 1))
    al1::Matrix{Float64} = zeros(ilong * ((jhi - 1) - (julo + 1) + 1))
    ar1::Matrix{Float64} = zeros(ilong * ((jhi - 1) - (julo + 1) + 1))
    dcy1::Matrix{Float64} = zeros(ilong * ((jhi - 1) - (julo + 1) + 1))
    qqu1::Matrix{Float64} = zeros(ilong * ((jhi - 1) - (julo + 1) + 1))
    a6::Matrix{Float64} = zeros(ihi - ilo + 1, jhi - julo + 1)
    al::Matrix{Float64} = zeros(ihi - ilo + 1, jhi - julo + 1)
    ar::Matrix{Float64} = zeros(ihi - ilo + 1, jhi - julo + 1)

    # NOTE: The code was writtein with I1:I2 as the first dimension of AL, AR, A6, AL1, A61, AR1. However, the limits should really should be ILO:IHI. In practice, however, for a global grid (and OpenMP parallelization) ILO=I1 and IHI=I2. Nevertheless, we will change the limits to ILO:IHI to be consistent and to avoid future problems. (bmy, 12/5/08)

    a6 .= 0.0
    al .= 0.0
    ar .= 0.0

    r13 = 1.0 / 3.0
    r23 = 2.0 / 3.0

    for ij ∈ julo+1:jhi
        for il ∈ ilo:ihi
            al[il, ij] = 0.5 * (qqu_copy[il, ij-1] + qqu_copy[il, ij]) + (dcy_copy[il, ij-1] - dcy_copy[il, ij]) * r13
            ar[il, ij-1] = al[il, ij]
        end
    end

    do_fyppm_pole_i2d2!(al, ar, i1_gl, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2)

    for ij ∈ julo+1:jhi-1
        for il ∈ ilo:ihi
            a6[il, ij] = 3.0 * (qqu_copy[il, ij] + qqu_copy[il, ij] - (al[il, ij] + ar[il, ij]))
        end
    end

    if jlmt <= 2
        lenx = 0
        for ij ∈ julo+1:jhi-1
            for il ∈ ilo:ihi
                lenx += 1
                a61[lenx] = a6[il, ij]
                al1[lenx] = al[il, ij]
                ar1[lenx] = ar[il, ij]
                dcy1[lenx] = dcy_copy[il, ij]
                qqu1[lenx] = qqu_copy[il, ij]
            end
        end
        lmtppm(lenx, jlmt, a61, al1, ar1, dcy1, qqu1)
        lenx = 0
        for ij ∈ julo+1:jhi-1
            for il ∈ ilo:ihi
                lenx += 1
                a6[il, ij] = a61[lenx]
                al[il, ij] = al1[lenx]
                ar[il, ij] = ar1[lenx]
            end
        end
    end

    for ij ∈ j1p:j2p+1
        ijm1 = ij - 1

        for il ∈ ilo:ihi
            if cry_copy[il, ij] > 0.0
                qqv[il, ij] = ar[il, ijm1] + 0.5 * cry_copy[il, ij] * (al[il, ijm1] - ar[il, ijm1] + (a6[il, ijm1] * (1.0 - (r23 * cry_copy[il, ij]))))
            else
                qqv[il, ij] = al[il, ij] - 0.5 * cry_copy[il, ij] * (ar[il, ij] - al[il, ij] + (a6[il, ij] * (1.0 + (r23 * cry_copy[il, ij]))))
            end
        end
    end
end

"""
Sets "al" & "ar" at the Poles.

## Arguments
- `al` (`inout`; `[ilo:ihi, julo:jhi]`) - Left edge values of the test parabola
- `ar` (`inout`; `[ilo:ihi, julo:jhi]`) - Right edge values of the test parabola
- `i1_gl` (`in`) - Global min longitude index
- `i2_gl` (`in`) - Global max longitude index
- `ju1_gl` (`in`) - Global min latitude index
- `j2_gl` (`in`) - Global max latitude index
- `ilo` (`in`) - Local min longitude index
- `ihi` (`in`) - Local max longitude index
- `julo` (`in`) - Local min latitude index
- `jhi` (`in`) - Local max latitude index
- `i1` (`in`) - Local min longitude index
- `i2` (`in`) - Local max longitude index
- `ju1` (`in`) - Local min latitude index
- `j2` (`in`) - Local max latitude index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_fyppm_pole_i2d2!(
    al::Matrix{Float64},
    ar::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    i2d2::Int64 = i2_gl / 2

    for il::Int64 ∈ i1:i2d2
        al[il, ju1] = al[il+i2d2, ju1+1]
        al[il+i2d2, ju1] = al[il, ju1+1]
        ar[il, j2] = ar[il+i2d2, j2-1]
        ar[il+i2d2, j2] = ar[il, j2-1]
    end
end

"""
Sets "dq1" at the Poles.

## Arguments
- `geofac_pc` (`in`) - Special geometrical factor (geofac) for Polar cap
- `dq1` (`inout`; `[ilo:ihi, julo:jhi]`) - Species density [hPa]
- `qqv` (`inout`; `[ilo:ihi, julo:jhi]`) - Concentration contribution from N-S advection [mixing ratio]
- `fy` (`inout`; `[ilo:ihi, julo:jhi+1]`) - N-S mass flux [mixing ratio]
- `i1_gl` (`in`) - Global min longitude (I) index
- `i2_gl` (`in`) - Global max longitude (I) index
- `ju1_gl` (`in`) - Global min latitude (J) index
- `j2_gl` (`in`) - Global max latitude (J) index
- `j1p` (`in`) - Global latitude index at the edge of the S polar cap
- `j2p` (`in`) - Global latitude index at the edge of the N polar cap
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_ytp_pole_sum!(
    geofac_pc::Float64,
    dq1::Matrix{Float64},
    qqv::Matrix{Float64},
    fy::Matrix{Float64},
    i1_gl::Int64,
    i2_gl::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    j1p::Int64,
    j2p::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64
)::Nothing
    il::Int64, ik::Int64 = 0, 0
    ri2::Float64 = 0.0
    dq_np::Float64 = 0.0
    dq_sp::Float64 = 0.0
    dqik::Vector{Float64} = zeros(2)
    sumnp::Float64 = 0.0
    sumsp::Float64 = 0.0

    ri2 = i2_gl
    dqik .= 0.0

    # Integrate N-S flux around polar cap lat circle for each level

    sumsp = 0.0
    sumnp = 0.0

    for il ∈ i1:i2
        sumsp += qqv[il, j1p]
        sumnp += qqv[il, j2p+1]
    end

    # wrap in E-W direction

    if i1 == i1_gl
        dqik[1] = dq1[i1, ju1]
        dqik[2] = dq1[i1, j2]
    end

    # normalize and set inside polar cap

    dq_sp = dqik[1] - (sumsp / ri2 * geofac_pc)
    dq_np = dqik[2] + (sumnp / ri2 * geofac_pc)

    for il = i1:i2
        dq1[il, ju1] = dq_sp
        dq1[il, j2] = dq_np
        # save polar flux
        fy[il, ju1] = -(sumsp / ri2 * geofac_pc)
        fy[il, j2+1] = (sumnp / ri2 * geofac_pc)
    end

    if j1p != ju1_gl + 1
        for il = i1:i2
            dq1[il, ju1+1] = dq_sp
            dq1[il, j2-1] = dq_np
            # save polar flux
            fy[il, ju1+1] = -(sumsp / ri2 * geofac_pc)
            fy[il, j2] = (sumnp / ri2 * geofac_pc)
        end
    end
end

"""
The 1D "outer" flux form operator based on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for computing the fluxes in the vertical direction.

Fzppm was modified by S.-J. Lin, 12/14/98, to allow the use of the KORD=7 (klmt=4) option.  KORD=7 enforces the 2nd monotonicity constraint of Huynh (1996).  Note that in Huynh's original scheme, two constraints are necessary for the preservation of monotonicity.  To use Huynh's algorithm, it was modified as follows.  The original PPM is still used to obtain the first guesses for the cell edges, and as such Huynh's 1st constraint is no longer needed.  Huynh's median function is also replaced by a simpler yet functionally equivalent in-line algorithm.

## Arguments
- `klmt` (`in`) - Controls various options in vertical advection
- `delp1` (`in`; `[ilo:ihi, julo:jhi, k1:k2]`) - Pressure thickness, the pseudo-density in a hydrostatic system at t1 [hPa]
- `wz` (`in`; `[i1:i2, ju1:j2, k1:k2]`) - Large scale mass flux (per time step tdt) in the vertical direction as diagnosed from the hydrostatic relationship [hPa]
- `dq1` (`inout`; `[ilo:ihi, julo:jhi, k1:k2]`) - Species density [hPa]
- `qq1` (`in`; `[:, :, :]`) - Species concentration [mixing ratio]
- `fz` (`out`; `[ilo:ihi, julo:jhi, k1:k2]`) - Vertical flux [mixing ratio]
- `j1p` (`in`) - Global latitude index at the edge of the South polar cap; J1P=JU1_GL+1 for a polar cap of 1 latitude band; J1P=JU1_GL+2 for a polar cap of 2 latitude bands
- `ju1_gl` (`in`) - Global latitude index at the edge of the North polar cap; JU1_GL=JULO for a polar cap of 1 latitude band; JU1_GL=JULO-1 for a polar cap of 2 latitude bands
- `j2_gl` (`in`) - Global max latitude (J) index
- `ilo` (`in`) - Global min longitude (I) index
- `ihi` (`in`) - Global max longitude (I) index
- `julo` (`in`) - Global min latitude (J) index
- `jhi` (`in`) - Global max latitude (J) index
- `ilong` (`in`) - Number of longitudes
- `ivert` (`in`) - Number of vertical levels
- `i1` (`in`) - Local min longitude (I) index
- `i2` (`in`) - Local max longitude (I) index
- `ju1` (`in`) - Local min latitude (J) index
- `j2` (`in`) - Local max latitude (J) index
- `k1` (`in`) - Local min vertical level (K) index
- `k2` (`in`) - Local max vertical level (K) index

## Author
Original code from Shian-Jiann Lin, DAO \\
John Tannahill, LLNL (jrt@llnl.gov)

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere. \\
See https://github.com/geoschem/geos-chem for complete history
"""
function fzppm!(
    klmt::Int64,
    delp1::Array{Float64,3},
    wz::Array{Float64,3},
    dq1::Array{Float64,3},
    qq1::Array{Float64,3},
    fz::Array{Float64,3},
    j1p::Int64,
    ju1_gl::Int64,
    j2_gl::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64,
    ilong::Int64,
    ivert::Int64,
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    k1::Int64,
    k2::Int64
)::Nothing
    delp1_copy::Array{Float64,3} = copy(delp1)
    wz_copy::Array{Float64,3} = copy(wz)
    qq1_copy::Array{Float64,3} = copy(qq1)

    il::Int64, ij::Int64, ik::Int64 = 0, 0, 0
    k1p1::Int64, k1p2::Int64 = 0, 0
    k2m1::Int64, k2m2::Int64 = 0, 0
    lenx::Int64 = 0
    a1::Float64, a2::Float64 = 0.0, 0.0
    aa::Float64, bb::Float64 = 0.0, 0.0
    c0::Float64, c1::Float64, c2::Float64 = 0.0, 0.0, 0.0
    cm::Float64, cp::Float64 = 0.0, 0.0
    fac1::Float64, fac2::Float64 = 0.0, 0.0
    lac::Float64 = 0.0
    qmax::Float64, qmin::Float64 = 0.0, 0.0
    qmp::Float64 = 0.0
    r13::Float64, r23::Float64 = 0.0, 0.0
    tmp::Float64 = 0.0

    a61::Vector{Float64} = zeros(ilong * (ivert - 4))
    al1::Vector{Float64} = zeros(ilong * (ivert - 4))
    ar1::Vector{Float64} = zeros(ilong * (ivert - 4))
    dca1::Vector{Float64} = zeros(ilong * (ivert - 4))
    qq1a1::Vector{Float64} = zeros(ilong * (ivert - 4))
    a6::Matrix{Float64} = zeros(i1:i2, k1:k2)
    al::Matrix{Float64} = zeros(i1:i2, k1:k2)
    ar::Matrix{Float64} = zeros(i1:i2, k1:k2)
    dca::Matrix{Float64} = zeros(i1:i2, k1:k2)
    dlp1a::Matrix{Float64} = zeros(i1:i2, k1:k2)
    qq1a::Matrix{Float64} = zeros(i1:i2, k1:k2)
    wza::Matrix{Float64} = zeros(i1:i2, k1:k2)
    dc::Array{Float64,3} = zeros(i1:i2, ju1:j2, k1:k2)

    dpi::Array{Float64,3} = zeros(i1:i2, ju1:j2, k1:k2)

    a6[:, :] = 0.0
    al[:, :] = 0.0
    ar[:, :] = 0.0
    dc[:, :, :] = 0.0
    # .sds... diagnostic vertical flux for species - set top to 0.0
    fz[:, :, :] = 0.0

    k1p1 = k1 + 1
    k1p2 = k1 + 2
    k2m1 = k2 - 1
    k2m2 = k2 - 2
    r13 = 1.0 / 3.0
    r23 = 2.0 / 3.0

    dpi[:, :, k1:k2m1] .=
        qq1_copy[i1:i2, ju1:j2, k1p1:k2] .-
        qq1_copy[i1:i2, ju1:j2, k1:k2m1]

    for ik = k1p1:k2m1
        for ij = ju1:j2
            for il = i1:i2
                c0 = delp1_copy[il, ij, ik] / (delp1_copy[il, ij, ik-1] + delp1_copy[il, ij, ik] + delp1_copy[il, ij, ik+1])
                c1 = (delp1_copy[il, ij, ik-1] + (0.5 * delp1_copy[il, ij, ik])) / (delp1_copy[il, ij, ik+1] + delp1_copy[il, ij, ik])
                c2 = (delp1_copy[il, ij, ik+1] + (0.5 * delp1_copy[il, ij, ik])) / (delp1_copy[il, ij, ik-1] + delp1_copy[il, ij, ik])
                tmp = c0 * ((c1 * dpi[il, ij, ik]) + (c2 * dpi[il, ij, ik-1]))
                qmax = max(qq1_copy[il, ij, ik-1], qq1_copy[il, ij, ik], qq1_copy[il, ij, ik+1]) - qq1_copy[il, ij, ik]
                qmin = qq1_copy[il, ij, ik] - min(qq1_copy[il, ij, ik-1], qq1_copy[il, ij, ik], qq1_copy[il, ij, ik+1])
                dc[il, ij, ik] = sign(tmp) * min(abs(tmp), qmax, qmin)
            end
        end
    end

    # Loop over latitudes (to save memory).

    for ij ∈ ju1:j2
        if (ij == ju1_gl + 1 || ij == j2_gl - 1) && (j1p != ju1_gl + 1)
            continue
        end

        for ik ∈ k1:k2
            for il ∈ i1:i2
                dca[il, ik] = dc[il, ij, ik]  # the monotone slope
                wza[il, ik] = wz_copy[il, ij, ik]
                dlp1a[il, ik] = delp1_copy[il, ij, ik]
                qq1a[il, ik] = qq1_copy[il, ij, ik]
            end
        end

        # Compute first guesses at cell interfaces. First guesses are required to be continuous. Three-cell parabolic subgrid distribution at model top; two-cell parabolic with zero gradient subgrid distribution at the surface.

        # First guess top edge value.

        for il ∈ i1:i2
            # Three-cell PPM; compute a, b, & c of q = aP^2 + bP + c using cell averages and dlp1a.
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
    end

    # Bottom.

    for il ∈ i1:i2
        # 2-cell PPM with zero gradient right at the surface.
        fac1 = dpi[il, ij, k2m1] * (dlp1a[il, k2] * dlp1a[il, k2]) / ((dlp1a[il, k2] + dlp1a[il, k2m1]) * (2.0 * dlp1a[il, k2] + dlp1a[il, k2m1]))
        ar[il, k2] = qq1a[il, k2] + fac1
        al[il, k2] = qq1a[il, k2] - (fac1 + fac1)

        if qq1a[il, k2] * ar[il, k2] <= 0.0
            ar[il, k2] = 0.0
        end

        dca[il, k2] = ar[il, k2] - qq1a[il, k2]
    end

    # 4th order interpolation in the interior.

    for ik ∈ k1p2:k2m1
        for il ∈ i1:i2
            c1 = (dpi[il, ij, ik-1] * dlp1a[il, ik-1]) / (dlp1a[il, ik-1] + dlp1a[il, ik])
            c2 = 2.0 / (dlp1a[il, ik-2] + dlp1a[il, ik-1] + dlp1a[il, ik] + dlp1a[il, ik+1])
            a1 = (dlp1a[il, ik-2] + dlp1a[il, ik-1]) / (2.0 * dlp1a[il, ik-1] + dlp1a[il, ik])
            a2 = (dlp1a[il, ik] + dlp1a[il, ik+1]) / (2.0 * dlp1a[il, ik] + dlp1a[il, ik-1])
            al[il, ik] = qq1a[il, ik-1] + c1 + c2 * (dlp1a[il, ik] * (c1 * (a1 - a2) + a2 * dca[il, ik-1]) - dlp1a[il, ik-1] * a1 * dca[il, ik])
        end
    end

    for ik ∈ k1:k2
        for il ∈ i1:i2
            ar[il, ik] = al[il, ik+1]
        end
    end

    # Top & Bottom 2 layers always monotonic.

    lenx = i2 - i1 + 1

    for ik ∈ k1:k1p1
        for il ∈ i1:i2
            a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (al[il, ik] + ar[il, ik]))
        end

        lmtppm!(lenx, 0, a6[i1, ik], al[i1, ik], ar[i1, ik], dca[i1, ik], qq1a[i1, ik])
    end

    for ik ∈ k2m1:k2
        for il ∈ i1:i2
            a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (al[il, ik] + ar[il, ik]))
        end

        lmtppm!(lenx, 0, a6[i1, ik], al[i1, ik], ar[i1, ik], dca[i1, ik], qq1a[i1, ik])
    end

    # Interior depending on klmt.
    if klmt == 4
        # KORD=7, Huynh's 2nd constraint.
        for ik ∈ k1p1:k2m1
            for il ∈ i1:i2
                dca[il, ik] = dpi[il, ij, ik] - dpi[il, ij, ik-1]
            end
        end

        for ik ∈ k1p2:k2m2
            for il ∈ i1:i2
                # Right edges.

                qmp = qq1a[il, ik] + (2.0 * dpi[il, ij, ik-1])
                lac = qq1a[il, ik] + (1.5 * dca[il, ik-1]) + (0.5 * dpi[il, ij, ik-1])
                qmin = min(qq1s[il, ik], qmp, lac)
                qmax = max(qq1s[il, ik], qmp, lac)
                ar[il, ik] = min(qmax, max(qmin, qq1a[il, ik]))

                # Left edges.

                qmp = qq1a[il, ik] - (2.0 * dpi[il, ij, ik-1])
                lac = qq1a[il, ik] - (1.5 * dca[il, ik-1]) - (0.5 * dpi[il, ij, ik-1])
                qmin = min(qq1s[il, ik], qmp, lac)
                qmax = max(qq1s[il, ik], qmp, lac)
                al[il, ik] = min(qmax, max(qmin, qq1a[il, ik]))

                # Recompute a6.

                a6[il, ik] = 3.0 * (qq1a[il, ik] + qq1a[il, ik] - (ar[il, ik] + al[il, ik]))
            end
        end
    elseif klmt <= 2
        lenx = 0

        for ik ∈ k1p2:k2m2
            for il ∈ i1:i2
                lenx += 1

                al1[lenx] = al[il, ik]
                ar1[lenx] = ar[il, ik]
                dca1[lenx] = dca[il, ik]
                qq1a1[lenx] = qq1a[il, ik]

                a61[lenx] = 3.0 * (qq1a1[lenx] + qq1a1[lenx] - (al1[lenx] + ar1[lenx]))
            end
        end

        lmtppm!(lenx, klmt, a61, al1, ar1, dca1, qq1a1)

        lenx = 0

        for ik ∈ k1p2:k2m2
            for il ∈ i1:i2
                lenx += 1

                a6[il, ik] = a61[lenx]
                al[il, ik] = al1[lenx]
                ar[il, ik] = ar1[lenx]
                dca[il, ik] = dca1[lenx]
                qq1a[il, ik] = qq1a1[lenx]
            end
        end

        for ik ∈ k1:k2m1
            for il ∈ i1:i2
                if wza[il, ik] > 0.0
                    cm = wza[il, ik] / dlp1a[il, ik]

                    dca[il, ik+1] = ar[il, ik] + 0.5 * cm * (al[il, ik] - ar[il, ik] + a6[il, ik] * (1.0 - r23 * cm))
                else
                    cp = wza[il, ik] / dlp1a[il, ik+1]

                    dca[il, ik+1] = al[il, ik+1] + 0.5 * cp * (al[il, ik+1] - ar[il, ik+1] - a6[il, ik+1] * (1.0 + r23 * cp))
                end
            end
        end

        for ik ∈ k1:k2m1
            for il ∈ i1:i2
                dca[il, ik+1] *= wza[il, ik]
                # .sds.. save vertical flux for species as diagnostic
                fz[il, ij, ik+1] = dca[il, ik+1]
            end
        end

        for il ∈ i1:i2
            dq1[il, ij, k1] -= dca[il, k1p1]
            dq1[il, ij, k2] += dca[il, k2]
        end

        for ik ∈ k1p1:k2m1
            for il ∈ i1:i2
                dq1[il, ij, ik] += dca[il, ik] - dca[il, ik+1]
            end
        end

    end
end

"""
Averages pressure at the Poles when the Polar cap is enlarged. It makes the last two latitudes equal.

## Arguments
- `area_1d` (`in`; `[ju1:j2]`) - Surface area of grid box
- `press` (`inout`; `[ilo:ihi, julo:jhi]`) - Surface pressure [hPa]
- `i1` (`in`) - Local min longitude (I)
- `i2` (`in`) - Local max longitude (I)
- `ju1` (`in`) - Local min latitude (J)
- `j2` (`in`) - Local max latitude (J)
- `ilo` (`in`) - Local min longitude (I) index
- `ihi` (`in`) - Local max longitude (I) index
- `julo` (`in`) - Local min latitude (J) index
- `jhi` (`in`) - Local max latitude (J) index

## Author
Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003) \\
Implemented into GEOS-Chem by Claire Carouge (ccarouge@seas.harvard.edu)

## Remarks
Subroutine from pjc_pfix. Call this one once everything is working fine.

## Revision History
05 Dec 2008 - C. Carouge - Replaced TPCORE routines by S-J Lin and Kevin Yeh with the TPCORE routines from GMI model. This eliminates the polar overshoot in the stratosphere.
See https://github.com/geoschem/geos-chem for complete history
"""
function average_press_poles!(
    area_1d::Vector{Float64},
    press::Matrix{Float64},
    i1::Int64,
    i2::Int64,
    ju1::Int64,
    j2::Int64,
    ilo::Int64,
    ihi::Int64,
    julo::Int64,
    jhi::Int64
)::Nothing
    area_1d_copy::Vector{Float64} = copy(area_1d)

    i::Int64, j::Int64 = 0, 0
    meanp::Float64 = 0.0
    rel_area::Vector{Float64} = zeros(ju1:j2)
    sum_area::Float64 = 0.0

    # Compute the sum of surface area
    sum_area = sum(area_1d_copy) * i2

    # Calculate rel_area for each lat. (ccc, 11/20/08)
    for j ∈ ju1:j2
        rel_area[j] = area_1d_copy[j] / sum_area
    end

    # South Pole

    # Surface area of the S. Polar cap
    sum_area = sum(rel_area[ju1:ju1+1]) * i2

    # Zero
    meanp = 0.0

    # Sum pressure * surface area over the S. Polar cap
    for j ∈ ju1:ju1+1
        for i ∈ i1:i2
            meanp += rel_area[j] * press[i, j]
        end
    end

    # Normalize pressure in all grid boxes w/in the S. Polar cap
    press[:, ju1:ju1+1] .= meanp / sum_area

    # North Pole

    # Surface area of the N. Polar cap
    sum_area = sum(rel_area[j2-1:j2]) * i2

    # Zero
    meanp = 0.0

    # Sum pressure * surface area over the N. Polar cap
    for j ∈ j2-1:j2
        for i ∈ i1:i2
            meanp += rel_area[j] * press[i, j]
        end
    end

    # Normalize pressure in all grid boxes w/in the N. Polar cap
    press[:, j2-1:j2] .= meanp / sum_area
end

end

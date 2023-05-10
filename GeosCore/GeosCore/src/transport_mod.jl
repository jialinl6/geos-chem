"""
Used to call the proper version of the TPCORE advection scheme for different met field datasets and their nested or global grids.

## Revision History

10 Mar 2003 - Y. Wang, R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
module transport_mod

using GeosCore.tpcore_fvdas_mod
# import ..tpcore_fvdas_mod: init_tpcore!, tpcore_fvdas!, cleanup_tpcore!

# using .GeosCore.tpcore_fvas_mod

import Headers
# import tpcore_fvas_mod

export cleanup_transport!, do_transport!, init_transport!, init_window_transport!

jfirst::Int64 = 0
# For fvDAS TPCORE
jlast::Int64 = 0
# For fvDAS TPCORE
ng::Int64 = 0
# For fvDAS TPCORE
mg::Int64 = 0
# For fvDAS TPCORE
n_adj::Int64 = 0
# Vertical coordinate array for TPCORE
ap::Vector{Float64} = []
# Vertical coordinate array for TPCORE
bp::Vector{Float64} = []
# Grid box surface areas [m2]
a_m2::Vector{Float64} = []

do_transport_first::Bool = true

mutable struct OptInput
end
mutable struct ChmState
end
mutable struct DgnState
end
mutable struct GrdState
end
mutable struct MetState
end

"""
Driver routine for the proper TPCORE program for GEOS-3, GEOS-4/GEOS-5, or window simulations.

## Arguments
- `input_opt` (`in`) - Input Options object
- `state_chm` (`inout`) - Chemistry State object
- `state_diag` (`inout`) - Diagnostics State object
- `state_grid` (`in`) - Grid State object
- `state_met` (`inout`) - Meteorology State object
- `rc` (`out`) - Success or failure?

## Revision History
10 Mar 2003 - R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_transport!(
    input_opt::OptInput,
    state_chm::ChmState,
    state_diag::DgnState,
    state_grid::GrdState,
    state_met::MetState,
    rc::Int64,
)::Nothing
    buff_size::Int64
    i0_w1::Int64
    j0_w1::Int64
    im_w1::Int64
    jm_w1::Int64
    ts_dyn::Int64
    dt_dyn::Float64
    errmsg::String
    thisloc::String
    rc = GC_SUCCESS

    errmsg = ""
    thisloc = " -> at do_transport! (in module GeosCore/transport_mod.F90)"

    # Transport budget diagnostics - Part 1 of 2
    if state_diag.archive_budgettransport
        # Get initial column masses (full, trop, PBL)
        compute_budget_diagnostics!(
            input_opt,
            state_chm,
            state_grid,
            state_met,
            state_diag.archive_budgettransportfull,
            nothing,
            state_diag.map_budgettransportfull,
            state_diag.archive_budgettransporttrop,
            nothing,
            state_diag.map_budgettransporttrop,
            state_diag.archive_budgettransportpbl,
            nothing,
            state_diag.map_budgettransportpbl,
            state_diag.budgetcolumnmass,
            true,
            rc,
        )

        if rc != GC_SUCCESS
            errmsg = "Transport budget diagnostics error 1"
            gc_error!(errmsg, rc, thisloc)
            return
        end
    end

    if first
        if state_grid.nestedgrid
            # All nested grid simulations
            init_window_transport!(input_opt, state_grid, rc)

            # Trap potential errors
            if rc != GC_SUCCESS
                errmsg = "Error encountered in init_window_transport!"
                gc_error!(errmsg, rc, thisloc)
                return
            end
        else
            # All global simulations
            init_transport!(input_opt, state_grid, rc)

            # Trap potential errors
            if rc != GC_SUCCESS
                errmsg = "Error encountered in init_transport!"
                gc_error!(errmsg, rc, thisloc)
                return
            end
        end
    end

    if state_grid.nestedgrid
        # Set buffer size
        buff_size = 2
        # (lzh, 09/01/2014)
        im_w1 = (state_grid.nx - state_grid.westbuffer - state_grid.eastbuffer) + 2 * buff_size
        jm_w1 = (state_grid.ny - state_grid.southbuffer - state_grid.northbuffer) + 2 * buff_size
        i0_w1 = state_grid.westbuffer - buff_size
        j0_w1 = state_grid.southbuffer - buff_size

        # Nested-grid simulation with GEOS-FP/MERRA2 met
        do_window_transport!(
            i0_w1,
            im_w1,
            j0_w1,
            jm_w1,
            input_opt,
            state_chm,
            state_diag,
            state_grid,
            state_met,
            rc,
        )

        # Trap potential errors
        if rc != GC_SUCCESS
            errmsg = "Error encountered in do_window_transport!"
            gc_error!(errmsg, rc, thisloc)
            return
        end
    else
        # Choose the proper version of TPCORE for global simulations

        do_global_adv!(
            input_opt,
            state_chm,
            state_diag,
            state_grid,
            state_met,
            rc,
        )

        # Trap potential errors
        if rc != GC_SUCCESS
            errmsg = "Error encountered in do_global_adv!"
            gc_error!(errmsg, rc, thisloc)
            return
        end
    end

    # Transport budget diagnostics - Part 2 of 2
    if state_diag.archive_budgettransport
        # Dynamic timestep [s]
        ts_dyn = get_ts_dyn()
        dt_dyn = Int64(ts_dyn)

        # Get final column masses (full, trop, PBL)
        compute_budget_diagnostics!(
            input_opt,
            state_chm,
            state_grid,
            state_met,
            state_diag.archive_budgettransportfull,
            state_diag.budgettransportfull,
            state_diag.map_budgettransportfull,
            state_diag.archive_budgettransporttrop,
            state_diag.budgettransporttrop,
            state_diag.map_budgettransporttrop,
            state_diag.archive_budgettransportpbl,
            state_diag.budgettransportpbl,
            state_diag.map_budgettransportpbl,
            state_diag.budgetcolumnmass,
            dt_dyn,
            rc,
        )

        if rc != GC_SUCCESS
            errmsg = "Transport budget diagnostics error 2"
            gc_error!(errmsg, rc, thisloc)
            return
        end
    end
end

"""
Driver routine for TPCORE with the GMAO GEOS-FP or MERRA-2 met fields.

## Arguments
- `input_opt` (`in`) - Input Options object
- `state_chm` (`inout`) - Chemistry State object
- `state_diag` (`inout`) - Diagnostics State object
- `state_grid` (`in`) - Grid State object
- `state_met` (`inout`) - Meteorology State object
- `rc` (`out`) - Success or failure?

## Remarks
As of July 2016, we assume that all of the advected species are listed first in the species database.  This is the easiest way to pass a slab to the TPCORE routine.  This may change in the future. (bmy, 7/13/16).

Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP) are incremented upside-down (level 1 = top of the atmosphere). The levels order is reversed only when written out to diagnostic output.

## Revision History
30 Oct 2007 - R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_global_adv!(
    input_opt::OptInput,
    state_chm::ChmState,
    state_diag::DgnState,
    state_grid::GrdState,
    state_met::MetState,
    rc::Int64
)::Nothing
    lfill::Bool
    prtdebug::Bool
    iord::Int64
    jord::Int64
    kord::Int64
    i::Int64
    j::Int64
    l::Int64
    l2::Int64
    n::Int64
    n_dyn::Int64
    na::Int64
    nadvect::Int64
    d_dyn::Float64

    p_tp1::Matrix{Float64} = zeros(state_grid.nx, state_grid.ny)
    p_tp2::Matrix{Float64} = zeros(state_grid.nx, state_grid.ny)
    p_temp::Matrix{Float64} = zeros(state_grid.nx, state_grid.ny)
    xmass::Array{Float64,3} = zeros(state_grid.nx, state_grid.ny, state_grid.nz)
    ymass::Array{Float64,3} = zeros(state_grid.nx, state_grid.ny, state_grid.nz)

    p_uwnd::Array{Float64,3}
    p_vwnd::Array{Float64,3}
    p_xmass::Array{Float64,3}
    p_ymass::Array{Float64,3}

    # Assume success
    rc = GC_SUCCESS

    lfill = input_opt.lfill
    prtdebug = input_opt.lprt && input_opt.amiroot
    iord = input_opt.tpcore_iord
    jord = input_opt.tpcore_jord
    kord = input_opt.tpcore_kord
    nadvect = state_chm.nadvect

    # Dynamic timestep [s]
    n_dyn = get_ts_dyn()
    d_dyn = n_dyn |> trunc |> Int64

    # Prepare variables for calls to PJC pressure-fixer and TPCORE
    # For hybrid grids, the pressure at the bottom edge of grid box (I,J,L) is given by:
    # P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    # where Psurface is the true surface pressure (i.e. not PS-PTOP).
    # and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)

    # DEBUG: Print a few global species sums
    # if prtdebug
    #     print_global_species_kg!(20, 20, 1, "O3", input_opt, state_chm, state_grid, state_met, "do_global_adv: pre-advection", rc)
    # end

    for j in 1:state_grid.ny
        for i in 1:state_grid.nx
            # Set true dry sfc pressure at midpoint of dynamic timestep [hPa]
            p_tp1[i, j] = get_pedge_dry(i, j, 1)
            # Set true dry sfc pressure at end of dynamic timestep [hPa]
            p_tp2[i, j] = state_met.psc2_dry[i, j]
        end
    end

    # Call the PJC/LLNL pressure fixer to get the adjusted air masses XMASS and YMASS. XMASS and YMASS need to be passed to TPCORE_FVDAS in order to ensure mass conservation.

    # NOTE: P_TP1 and P_TP2 are the true surface pressures!
    do_pjc_pfix!(
        state_grid,
        d_dyn,
        p_tp1,
        p_tp2,
        state_met.u,
        state_met.v,
        xmass,
        ymass,
        rc,
    )

    # Call TPCORE_FVDAS to perform the advection

    # Flip aray indices in the vertical using pointer storage
    p_uwnd = state_met.u[:, :, state_grid.nz:-1:1]
    p_vwnd = state_met.v[:, :, state_grid.nz:-1:1]
    p_xmass = xmass[:, :, state_grid.nz:-1:1]
    p_ymass = ymass[:, :, state_grid.nz:-1:1]

    # NOTE: For now, so as to avoid having to rewrite the internals of the TPCORE routines, just point to 1:nAdvect entries of State_Chm%Species.  This is OK for now, as of July 2016, all of the advected species are listed first.  This may change in the future, but we'll worry about that later.  The units of p_SPC will be converted to [kg/kg moist air] below. (bmy, 7/13/16)

    tpcore_fvdas!(
        d_dyn, re, state_grid.nx, state_grid.ny,
        state_grid.nz, jfirst, jlast, ng,
        mg, nadvect, ap, bp,
        p_uwnd, p_vwnd, p_tp1, p_tp2,
        p_temp, iord, jord, kord,
        n_adj, p_xmass, p_ymass, lfill,
        a_m2, state_chm, state_diag
    )

    # Reset surface pressure and ensure mass conservation

    # Update dry and wet floating pressures to the most recently interpolated values (State_Met%PSC2_DRY and State_Met%PSC2) (ewl, 7/6/16)
    set_floating_pressures!(state_grid, state_met, rc)

    # Update State_Met air quantities with new pressures and update tracer mixing ratio because after advection the mixing ratio values reflect the new air pressure (ewl, 7/6/16)
    airqnt!(input_opt, state_chm, state_grid, state_met, rc, false)

    # DEBUG: Print a few global species sums
    # if prtdebug
    #     print_global_species_kg!(20, 20, 1, "O3", input_opt, state_chm, state_grid, state_met, "do_global_adv: post-airqnt", rc)
    # end
end

"""
The driver program for the proper TPCORE program for the GEOS-FP/MERRA2 nested-grid simulations.

## Arguments
- `i0` (`in`) - The starting i-index for the window
- `im` (`in`) - The ending i-index for the window
- `j0` (`in`) - The starting j-index for the window
- `jm` (`in`) - The ending j-index for the window
- `input_opt` (`in`) - The input options object
- `state_chm` (`inout`) - The chemistry state object
- `state_diag` (`inout`) - The diagnostics state object
- `state_grid` (`in`) - The grid state object
- `state_met` (`inout`) - The meteorology state object
- `rc` (`out`) - The return code

## Remarks
As of July 2016, we assume that all of the advected species are listed first in the species database.  This is the easiest way to pass a slab to the TPCORE routine.  This may change in the future. (bmy, 7/13/16)

Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP) are incremented upside-down (level 1 = top of the atmosphere).  The levels order is reversed only when written out to diagnostic output.

## Revision History
10 Mar 2003 - R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
function do_window_transport!(
    i0::Int64,
    im::Int64,
    j0::Int64,
    jm::Int64,
    input_opt::OptInput,
    state_chm::ChmState,
    state_diag::DgnState,
    state_grid::GrdState,
    state_met::MetState,
    rc::Int64
)::Nothing
    na::Int64 = 0
    i::Int64, j::Int64, l::Int64, l2::Int64, n::Int64 = 0, 0, 0, 0, 0

    p_tp1::Matrix{Float64} = zeros(state_grid.nx, state_grid.ny)
    p_tp2::Matrix{Float64} = zeros(state_grid.nx, state_grid.ny)
    p_temp::Matrix{Float64} = zeros(state_grid.nx, state_grid.ny)
    xmass::Array{Float64,3} = zeros(state_grid.nx, state_grid.ny, state_grid.nz)
    ymass::Array{Float64,3} = zeros(state_grid.nx, state_grid.ny, state_grid.nz)
    q_spc::Array{Float64,4} = zeros(im, jm, state_grid.nz, state_chm.nadvect)

    p_a_m2::Vector{Float64} = zeros()
    p_p_tp1::Matrix{Float64} = zeros()
    p_p_tp2::Matrix{Float64} = zeros()
    p_p_temp::Matrix{Float64} = zeros()
    p_uwnd::Array{Float64,3} = zeros()
    p_vwnd::Array{Float64,3} = zeros()
    p_xmass::Array{Float64,3} = zeros()
    p_ymass::Array{Float64,3} = zeros()
    p_spc::Array{Float64,4} = zeros()

    lfill::Bool = input_opt.lfill
    iord::Int64 = input_opt.tpcore_iord
    jord::Int64 = input_opt.tpcore_jord
    kord::Int64 = input_opt.tpcore_kord
    nadvect::Int64 = state_chm.nadvect
    prtdebug::Bool = input_opt.lprt && input_opt.am_i_root

    # Dynamic timestep [s]
    # TODO:
    # get_ts_dyn is a custom function that returns the diagnostic timestep in seconds.
    n_dyn::Int64 = get_ts_dyn()
    d_dyn::Float64 = n_dyn |> Float64

    # Set start and end indices for the window
    ia::Int64 = i0 + 1
    ib::Int64 = i0 + im
    ja::Int64 = j0 + 1
    jb::Int64 = j0 + jm

    # Set local array for species concentrations in window
    for n ∈ 1:state_chm.nadvect
        q_spc[:, :, :, n] .= state_chm.species[n].conc[ia:ib, ja:jb, :]
    end

    # Prepare variables for calls to PJC pressure-fixer and TPCORE
    # For hybrid grids, the pressure at the bottom edge of grid box (I,J,L) is given by:
    # P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    # where Psurface is the true surface pressure (i.e. not PS-PTOP). and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
    # if prtdebug
    #     print_global_species_kg(20, 20, 1, "SPC_O3", input_opt, state_chm, state_grid, state_met, "do_window_transport: pre-advection", rc)
    # end

    for j ∈ 1:state_grid.ny
        for i ∈ 1:state_grid.nx
            # TODO: I'm not sure where get_pedge_dry is defined
            # Set true dry sfc pressure at midpoint of dynamic timestep [hPa]
            p_tp1[i, j] = get_pedge_dry(i, j, 1)
            # Set true dry sfc pressure at end of dynamic timestep [hPa]
            p_tp2[i, j] = state_met.psc2_dry[i, j]
        end
    end

    # Call the PJC/LLNL pressure fixer to get the adjusted air masses XMASS and YMASS.  XMASS and YMASS need to be passed to TPCORE_FVDAS in order to ensure mass conservation.

    # dan
    xmass .= 0.0
    ymass .= 0.0
    # NOTE: P_TP1 and P_TP2 are the true surface pressures!
    do_pjc_pfix_window(state_grid, d_dyn, p_tp1, p_tp2, state_met.u, state_met.v, xmass, ymass)

    if prtdebug
        debug_msg!("### FVDAS_WINDOW: b TPCORE_WINDOW")
    end

    # Flip array indices in the vertical using pointer storage
    # Exclude the buffer zone (lzh, 4/1/2015)
    p_a_m2 = a_m2[ja:jb]
    p_p_tp1 = p_tp1[ia:ib, ja:jb]
    p_p_tp2 = p_tp2[ia:ib, ja:jb]
    p_p_temp = p_temp[ia:ib, ja:jb]
    p_uwnd = state_met.u[ia:ib, ja:jb, state_grid.nz:-1:1]
    p_vwnd = state_met.v[ia:ib, ja:jb, state_grid.nz:-1:1]
    p_xmass = xmass[ia:ib, ja:jb, state_grid.nz:-1:1]
    p_ymass = ymass[ia:ib, ja:jb, state_grid.nz:-1:1]
    p_spc = q_spc[:, :, state_grid.nz:-1:1, :]

    # Do the advection
    tpcore_window!(
        d_dyn,
        re,
        im,
        jm,
        state_grid.nz,
        jfirst,
        jlast,
        ng,
        mg,
        nadvect,
        ap,
        bp,
        p_uwnd,
        p_vwnd,
        p_p_tp1,
        p_p_tp2,
        p_p_temp,
        p_spc,
        iord,
        jord,
        kord,
        n_adj,
        p_xmass,
        p_ymass,
        p_a_m2,
        state_chm,
        state_diag
    )

    # Update species concentrations from local array
    for n in 1:state_chm.nadvect
        state_chm.species[n].conc[ia:ib, ja:jb, :] .= q_spc[:, :, :, n]
    end

    # Reset surface pressure to the true surface pressure

    # Update dry and wet floating pressures to the most recently interpolated values (State_Met%PSC2_DRY and State_Met%PSC2) (ewl, 7/6/16)
    set_floating_pressures!(state_grid, state_met, rc)

    # Update State_Met air quantities with new pressures. Do not update tracer mixing ratio because after advection the mixing ratio values reflect the new air pressure (ewl, 3/31/15)
    airqnt!(input_opt, state_chm, state_grid, state_met, rc, false)

    # Debug
    # if prtdebug
    #     print_global_species_kg(20, 20, 1, "SPC_O3", input_opt, state_chm, state_grid, state_met, "do_window_transport: post-airqnt", rc)
    #     debug_msg!("### FVDAS_WINDOW: a TPCORE_WINDOW")
    # end
end

"""
Initializes all module variables and arrays.

## Revision History
10 Mar 2003 - R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
function init_transport!()::Nothing
    ltran::Bool = false
    j::Int, k::Int, l::Int, n_dyn::Int = 0, 0, 0, 0
    ymid_r::Vector{Float64} = zeros(state_grid.ny)
    real_n_dyn::Float64 = 0.0

    rc = GC_SUCCESS
    errmsg::String = ""
    thisloc::String = " -> at init_transport (in module GeosCore/transport_mod.F90)"

    # Allocate arrays for TPCORE vertical coordinates
    # For fvDAS TPCORE with for GEOS-FP or MERRA-2 met fields:
    #    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
    # Also here Ap, Bp will be flipped since both TPCORE versions index levels from the atm. top downwards (bdf, bmy, 10/30/07)
    bp = zeros(state_grid.nz + 1)

    # Flip Ap and Bp for TPCORE
    for l ∈ 1:state_grid.nz+1
        # As L runs from the surface up, K runs from the top down
        k = state_grid.nz + 1 - l + 1
        # Ap(L) is in [hPa]
        ap[l] = get_ap(k)
        bp[l] = get_bp(k)
    end

    # Allocate arrays for surface area and layer thickness
    a_m2 = zeros(state_grid.ny)

    # Surface area [m2]
    for j ∈ 1:state_grid.ny
        a_m2[j] = state_grid.area_m2[1, j]
    end

    # Additional setup for the fvDAS version of TPCORE

    n_dyn = get_ts_dyn()
    n_adj = 0
    ng = 0
    mg = 0

    # YMID_R is latitude of grid box center [radian]
    for j ∈ 1:state_grid.ny
        ymid_r[j] = state_grid.ymid_r[1, j]
    end

    real_n_dyn = n_dyn

    init_tpcore!(
        state_grid.nx,
        state_grid.ny,
        state_grid.nz,
        jfirst,
        jlast,
        ng,
        mg,
        real_n_dyn,
        re,
        ymid_r,
        rc
    )

    # Trap potential errors
    if rc != GC_SUCCESS
        err_msg = "Error encountered in \"Init_Tpcore\"!"
        gc_error(err_msg, rc, this_loc)
        return
    end
end

"""
Initializes all module variables and arrays for the GEOS-FP/MERRA2 nested grid simulation.

## Arguments
- `input_opt` (`in`) - Input options object
- `state_grid` (`in`) - Grid state object
- `rc` (`out`) - Success or failure?

## Revision History
06 Jun 2008 - D. Chen & R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
function init_window_transport!(
    input_opt::OptInput,
    state_grid::GrdState,
    rc::Int64,
)::Nothing
    ltran::Bool
    buff_size::Int64
    j::Int64
    k::Int64
    l::Int64
    n_dyn::Int64
    im_w1::Int64
    jm_w1::Int64
    i0_w1::Int64
    j0_w1::Int64
    ymid_r_w::Vector{Float64} = zeros(state_grid.ny + 2)

    # Assume success
    rc = GC_SUCCESS

    # Copy values from Input_Opt
    ltran = input_opt.ltran

    # Allocate arrays for TPCORE vertical coordinates (GEOS-FP/MERRA2 nested grid simulation only!!!)
    # For fvDAS TPCORE with for GEOS-FP/MERRA2 met fields:
    #    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
    # Also here Ap, Bp will be flipped since both TPCORE versions index levels from the atm. top downwards (bdf, bmy, 10/30/07)
    global ap = zeros(state_grid.nz + 1)
    global bp = zeros(state_grid.nz + 1)

    for l in 1:state_grid.nz+1
        # As L runs from the surface up, K runs from the top down
        k = state_grid.nz + 1 - l + 1
        # NOTE: `get_ap` and `get_bp` are not defined in this file
        ap[l] = get_ap(k)
        bp[l] = get_bp(k)
    end

    # Allocate arrays for surface area and layer thickness
    global a_m2 = zeros(state_grid.ny)

    # Surface area [m2]
    for j in 1:state_grid.ny
        a_m2[j] = state_grid.area_m2[1, j]
    end

    # Additional setup for the fvDAS version of TPCORE
    n_dyn = get_ts_dyn()
    n_adj = 0
    ng = 0
    mg = 0

    # (lzh, 4/1/2015)
    buff_size = 2
    im_w1 = (state_grid.nx - state_grid.westbuffer - state_grid.eastbuffer) + 2 * buff_size
    jm_w1 = (state_grid.ny - state_grid.southbuffer - state_grid.northbuffer) + 2 * buff_size
    i0_w1 = state_grid.westbuffer - buff_size
    j0_w1 = state_grid.southbuffer - buff_size

    # YMID_R is latitude of grid box center [radians]
    for j in 1:state_grid.ny
        ymid_r_w[j] = state_grid.ymid_r[1, j]
    end

    # Compute YMID_R_W at southern edge of nested region
    j = 0
    ymid_r_w[j] = state_grid.ymid_r[1, j+1] - (state_grid.dy * PI_180)

    # Compute YMID_R_W at northern edge of nested region
    j = state_grid.ny + 1
    ymid_r_w[j] = state_grid.ymid_r[1, j-1] + (state_grid.dy * PI_180)

    # Call INIT routine from "tpcore_window_mod.F90"
    init_window!(
        state_grid,
        im_w1,
        jm_w1,
        state_grid.nz,
        jfirst,
        jlast,
        ng,
        mg,
        n_dyn,
        re,
        ymid_r_w[j0_w1:j0_w1+jm_w1+1],
    )
end

"""
Deallocate all module arrays.

## Revision History
10 Mar 2003 - R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
function cleanup_transport()

    # IF ( ALLOCATED( Ap     ) ) DEALLOCATE( Ap     )
    # IF ( ALLOCATED( A_M2   ) ) DEALLOCATE( A_M2   )
    # IF ( ALLOCATED( Bp     ) ) DEALLOCATE( Bp     )
end

end

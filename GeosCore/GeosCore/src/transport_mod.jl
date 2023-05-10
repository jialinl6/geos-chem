# !------------------------------------------------------------------------------
# !                  GEOS-Chem Global Chemical Transport Model                  !
# !------------------------------------------------------------------------------
# !BOP
# !
# ! !MODULE: transport_mod.F90
# !
# ! !DESCRIPTION: Module TRANSPORT\_MOD is used to call the proper version of
# !  the TPCORE advection scheme for different met field datasets and their
# !  nested or global grids.
# !\\
# !\\
# ! !INTERFACE:
# !
# MODULE TRANSPORT_MOD
module transport_mod
# !
# ! !USES:
# !
#   USE PRECISION_MOD      ! For GEOS-Chem Precision (fp)
#   USE PRESSURE_MOD

#   IMPLICIT NONE
#   PRIVATE
# !
# ! !PUBLIC MEMBER FUNCTIONS:
# !
#   PUBLIC  :: CLEANUP_TRANSPORT
#   PUBLIC  :: DO_TRANSPORT
#   PUBLIC  :: INIT_TRANSPORT
#   PUBLIC  :: INIT_WINDOW_TRANSPORT
# !
# ! !PRIVATE MEMBER FUNCTIONS:
# !
#   PRIVATE :: DO_GLOBAL_ADV
#   PRIVATE :: DO_WINDOW_TRANSPORT
# !
# ! !REVISION HISTORY:
# !  10 Mar 2003 - Y. Wang, R. Yantosca - Initial version
# !  See https://github.com/geoschem/geos-chem for complete history
# !EOP
# !------------------------------------------------------------------------------
# !BOC
#   !=================================================================
#   ! MODULE VARIABLES:
#   !
#   ! (1 ) Ap     (REAL(fp) ) : Vertical coordinate array for TPCORE
#   ! (2 ) A_M2   (REAL(fp) ) : Grid box surface areas [m2]
#   ! (3 ) Bp     (REAL(fp) ) : Vertical coordinate array for TPCORE
#   ! (7 ) JLAST  (INTEGER)   : For fvDAS TPCORE
#   ! (8 ) MG     (INTEGER)   : For fvDAS TPCORE
#   ! (9 ) NG     (INTEGER)   : For fvDAS TPCORE
#   ! (10) N_ADJ  (INTEGER)   : For fvDAS TPCORE
#   !=================================================================
#   INTEGER                       :: JFIRST
#   INTEGER                       :: JLAST, NG,   MG,   N_ADJ
#   REAL(fp), ALLOCATABLE         :: Ap(:)
#   REAL(fp), ALLOCATABLE         :: Bp(:)
#   REAL(fp), ALLOCATABLE, TARGET :: A_M2(:)

# CONTAINS
# !EOC
# !------------------------------------------------------------------------------
# !                  GEOS-Chem Global Chemical Transport Model                  !
# !------------------------------------------------------------------------------
# !BOP
# !
# ! !IROUTINE: do_transport
# !
# ! !DESCRIPTION: Subroutine DO\_TRANSPORT is the driver routine for the proper
# !  TPCORE program for GEOS-3, GEOS-4/GEOS-5, or window simulations.
# !\\
# !\\
# ! !INTERFACE:
# !
#   SUBROUTINE DO_TRANSPORT( Input_Opt,  State_Chm, State_Diag, &
#                            State_Grid, State_Met, RC )
# !
# ! !USES:
# !
#     USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
#     USE ErrCode_Mod
#     USE Input_Opt_Mod,   ONLY : OptInput
#     USE State_Chm_Mod,   ONLY : ChmState
#     USE State_Diag_Mod,  ONLY : DgnState
#     USE State_Grid_Mod,  ONLY : GrdState
#     USE State_Met_Mod,   ONLY : MetState
#     USE TIME_MOD,        ONLY : GET_TS_DYN
# !
# ! !INPUT PARAMETERS:
# !
#     TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
#     TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
# !
# ! !INPUT/OUTPUT PARAMETERS:
# !
#     TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
#     TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
#     TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
# !
# ! !OUTPUT PARAMETERS:
# !
#     INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
# !
# ! !REVISION HISTORY:
# !  10 Mar 2003 - R. Yantosca - Initial version
# !  See https://github.com/geoschem/geos-chem for complete history
# !EOP
# !------------------------------------------------------------------------------
# !BOC

do_transport_first::Bool = true

function do_transport!(
    input_opt::optinput,
    state_chm::chmstate,
    state_diag::dgnstate,
    state_grid::grdstate,
    state_met::metstate,
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
    input_opt::optinput,
    state_chm::chmstate,
    state_diag::dgnstate,
    state_grid::grdstate,
    state_met::metstate,
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

    #     ! Assume success
    #     RC          =  GC_SUCCESS

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

    #     !!### DEBUG: Print a few global species sums
    #     !IF ( prtDebug ) THEN
    #     !   CALL Print_Global_Species_Kg( 20, 20, 1, 'O3',       &
    #     !                                 Input_Opt, State_Chm,  &
    #     !                                 State_Grid, State_Met, &
    #     !                                 "do_global_adv: pre-advection", &
    #     !                                 RC )
    #     !ENDIF

    #     !$OMP PARALLEL DO       &
    #     !$OMP DEFAULT( SHARED ) &
    #     !$OMP PRIVATE( I, J )
    #     DO J = 1, State_Grid%NY
    #     DO I = 1, State_Grid%NX

    #        ! Set true dry sfc pressure at midpoint of dynamic timestep [hPa]
    #        P_TP1(I,J) = GET_PEDGE_DRY(I,J,1)

    #        ! Set true dry sfc pressure at end of dynamic timestep [hPa]
    #        P_TP2(I,J) = State_Met%PSC2_DRY(I,J)

    #     ENDDO
    #     ENDDO
    #     !$OMP END PARALLEL DO

    #     !=================================================================
    #     ! Call the PJC/LLNL pressure fixer to get the adjusted air
    #     ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to
    #     ! TPCORE_FVDAS in order to ensure mass conservation.
    #     !=================================================================

    #     ! NOTE: P_TP1 and P_TP2 are the true surface pressures!
    #     CALL DO_PJC_PFIX( State_Grid,  D_DYN, P_TP1, P_TP2, &
    #                       State_Met%U, State_Met%V, XMASS, YMASS )

    #     !=================================================================
    #     ! Call TPCORE_FVDAS to perform the advection
    #     !=================================================================

    #     ! Flip array indices in the vertical using pointer storage
    #     p_UWND    => State_Met%U      (:,:,State_Grid%NZ:1:-1)
    #     p_VWND    => State_Met%V      (:,:,State_Grid%NZ:1:-1)
    #     p_XMASS   => XMASS            (:,:,State_Grid%NZ:1:-1)
    #     p_YMASS   => YMASS            (:,:,State_Grid%NZ:1:-1)

    #     ! NOTE: For now, so as to avoid having to rewrite the internals
    #     ! of the TPCORE routines, just point to 1:nAdvect entries of
    #     ! State_Chm%Species.  This is OK for now, as of July 2016, all of
    #     ! the advected species are listed first.  This may change in the
    #     ! future, but we'll worry about that later.  The units of p_SPC
    #     ! will be converted to [kg/kg moist air] below. (bmy, 7/13/16)

    #     ! Do the advection
    #     CALL TPCORE_FVDAS( D_DYN,         Re,        State_Grid%NX, State_Grid%NY, &
    #                        State_Grid%NZ, JFIRST,    JLAST,         NG,            &
    #                        MG,            nAdvect,   Ap,            Bp,            &
    #                        p_UWND,        p_VWND,    P_TP1,         P_TP2,         &
    #                        P_TEMP,        IORD,      JORD,          KORD,          &
    #                        N_ADJ,         p_XMASS,   p_YMASS,       LFILL,         &
    #                        A_M2,          State_Chm, State_Diag                   )

    #     ! Free pointer memory
    #     p_UWND  => NULL()
    #     p_VWND  => NULL()
    #     p_XMASS => NULL()
    #     p_YMASS => NULL()

    #     !=================================================================
    #     ! Reset surface pressure and ensure mass conservation
    #     !=================================================================

    #     ! Update dry and wet floating pressures to the most recently
    #     ! interpolated values (State_Met%PSC2_DRY and State_Met%PSC2)
    #     ! (ewl, 7/6/16)
    #     CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC)

    #     ! Update State_Met air quantities with new pressures.
    #     ! Do not update tracer mixing ratio because after advection
    #     ! the mixing ratio values reflect the new air pressure (ewl, 3/31/15)
    #     CALL AIRQNT( Input_Opt, State_Chm, State_Grid, State_Met, &
    #                  RC, update_mixing_ratio=.FALSE. )

    #     !!### DEBUG: Print a few global species sums
    #     !IF ( prtDebug ) THEN
    #     !   CALL Print_Global_Species_Kg( 20, 20, 1, 'O3',       &
    #     !                                 Input_Opt, State_Chm,  &
    #     !                                 State_Grid, State_Met, &
    #     !                                 "do_global_adv: post-airqnt", &
    #     !                                 RC )
    #     !   CALL DEBUG_MSG( '### G4_G5_GLOB_ADV: a TPCORE' )
    #     !ENDIF

    #   END SUBROUTINE DO_GLOBAL_ADV
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
    input_opt::optinput,
    state_chm::chmstate,
    state_diag::dgnstate,
    state_grid::grdstate,
    state_met::metstate,
    rc::Int64
)::Nothing
    lfill::Bool
    prtdebug::Bool
    iord::Int64
    jord::Int64
    kord::Int64
    ia::Int64
    ib::Int64
    ja::Int64
    jb::Int64
    na::Int64
    nadvect::Int64
    n_dyn::Int64
    i::Int64
    j::Int64
    l::Int64
    l2::Int64
    n::Int64
    d_dyn::Float64

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

    lfill = input_opt.lfill
    iord = input_opt.tpcore_iord
    jord = input_opt.tpcore_jord
    kord = input_opt.tpcore_kord
    nadvect = state_chm.nadvect
    prtdebug = input_opt.lprt && input_opt.am_i_root

    #     ! Initialize pointers
    #     p_A_M2      => NULL()
    #     p_P_TP1     => NULL()
    #     p_P_TP2     => NULL()
    #     p_P_TEMP    => NULL()
    #     p_UWND      => NULL()
    #     p_VWND      => NULL()
    #     p_XMASS     => NULL()
    #     p_YMASS     => NULL()
    #     p_Spc       => NULL()

    # Dynamic timestep [s]
    # TODO:
    # get_ts_dyn is a custom function that returns the diagnostic timestep in seconds.
    n_dyn = get_ts_dyn()
    d_dyn = n_dyn

    # Set start and end indices for the window
    ia = i0 + 1
    ib = i0 + im
    ja = j0 + 1
    jb = j0 + jm

    #     ! Set local array for species concentrations in window
    #     DO N = 1, State_Chm%nAdvect
    #        Q_Spc(:,:,:,N)= State_Chm%Species(N)%Conc(IA:IB,JA:JB,:)
    #     ENDDO

    # Set local array for species concentrations in window
    for n ∈ 1:nadvect
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

    #     !=================================================================
    #     ! Call the PJC/LLNL pressure fixer to get the adjusted air
    #     ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to
    #     ! TPCORE_FVDAS in order to ensure mass conservation.
    #     !=================================================================
    #     XMASS = 0e+0_fp !(dan)
    #     YMASS = 0e+0_fp
    #     ! NOTE: P_TP1 and P_TP2 are the true surface pressures!
    #     CALL DO_PJC_PFIX_WINDOW( State_Grid,  D_DYN, &
    #                              P_TP1,       P_TP2, &
    #                              State_Met%U, State_Met%V, &
    #                              XMASS,       YMASS )

    #     IF ( prtDebug ) CALL DEBUG_MSG( '### FVDAS_WINDOW: a PJC_PFIX_WINDOW')

    #     ! Flip array indices in the vertical using pointer storage
    #     ! Exclude the buffer zone (lzh, 4/1/2015)
    #     p_A_M2   => A_M2       ( JA:JB                               )
    #     p_P_TP1  => P_TP1      ( IA:IB, JA:JB                        )
    #     p_P_TP2  => P_TP2      ( IA:IB, JA:JB                        )
    #     p_P_TEMP => P_TEMP     ( IA:IB, JA:JB                        )
    #     p_UWND   => State_Met%U( IA:IB, JA:JB, State_Grid%NZ:1:-1    )
    #     p_VWND   => State_Met%V( IA:IB, JA:JB, State_Grid%NZ:1:-1    )
    #     p_XMASS  => XMASS      ( IA:IB, JA:JB, State_Grid%NZ:1:-1    )
    #     p_YMASS  => YMASS      ( IA:IB, JA:JB, State_Grid%NZ:1:-1    )
    #     p_Spc    => Q_Spc      ( :,     :,     State_Grid%NZ:1:-1, : )

    #     ! Do the advection
    #     CALL TPCORE_WINDOW(D_DYN,   Re,        IM,      JM,      State_Grid%NZ, &
    #                        JFIRST,  JLAST,     NG,      MG,      nAdvect,       &
    #                        Ap,      Bp,        p_UWND,  p_VWND,  p_P_TP1,       &
    #                        p_P_TP2, p_P_TEMP,  p_Spc,   IORD,    JORD,          &
    #                        KORD,    N_ADJ,     p_XMASS, p_YMASS,                &
    #                        p_A_M2,  State_Chm, State_Diag )


    #     ! Update species concentrations from local array
    #     DO N = 1, State_Chm%nAdvect
    #        State_Chm%Species(N)%Conc(IA:IB,JA:JB,:) = Q_Spc(:,:,:,N)
    #     ENDDO

    #     ! Free pointer memory
    #     p_UWND   => NULL()
    #     p_VWND   => NULL()
    #     p_Spc    => NULL()
    #     p_XMASS  => NULL()
    #     p_YMASS  => NULL()
    #     p_P_TP1  => NULL()
    #     p_P_TP2  => NULL()
    #     p_P_TEMP => NULL()
    #     p_A_M2   => NULL()

    #     !=================================================================
    #     ! Reset surface pressure and ensure mass conservation
    #     !=================================================================

    #     ! Update dry and wet floating pressures to the most recently
    #     ! interpolated values (State_Met%PSC2_DRY and State_Met%PSC2)
    #     ! (ewl, 7/6/16)
    #     CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC)

    #     ! Update State_Met air quantities with new pressures.
    #     ! Do not update tracer mixing ratio because after advection
    #     ! the mixing ratio values reflect the new air pressure (ewl, 3/31/15)
    #     CALL AIRQNT( Input_Opt, State_Chm, State_Grid, State_Met, RC, &
    #                  Update_Mixing_Ratio=.FALSE. )

    #     !!### Debug
    #     !IF ( prtDebug ) THEN
    #     !   CALL Print_Global_Species_Kg( 20, 20, 1, 'SPC_O3',   &
    #     !                                 Input_Opt, State_Chm,  &
    #     !                                 State_Grid, State_Met, &
    #     !                                 "do_window_transport: post-airqnt", &
    #     !                                 RC )
    #     !   CALL DEBUG_MSG( '### NESTED_ADV: a TPCORE' )
    #     !ENDIF

    #   END SUBROUTINE DO_WINDOW_TRANSPORT
end

# !EOC
# !------------------------------------------------------------------------------
# !                  GEOS-Chem Global Chemical Transport Model                  !
# !------------------------------------------------------------------------------
# !BOP
# !
# ! !IROUTINE: init_transport
# !
# ! !DESCRIPTION: Subroutine INIT\_TRANSPORT initializes all module variables
# !  and arrays.
# !\\
# !\\
# ! !INTERFACE:
# !
#   SUBROUTINE INIT_TRANSPORT( Input_Opt, State_Grid, RC )
# !
# ! !USES:
# !
#     USE ErrCode_Mod
#     USE ERROR_MOD,        ONLY : ALLOC_ERR
#     USE Input_Opt_Mod,    ONLY : OptInput
#     USE State_Grid_Mod,   ONLY : GrdState
#     USE PhysConstants          ! Re
#     USE TIME_MOD,         ONLY : GET_TS_DYN
#     USE TPCORE_FVDAS_MOD, ONLY : INIT_TPCORE
# !
# ! !INPUT PARAMETERS:
# !
#     TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
#     TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
# !
# ! !OUTPUT PARAMETERS:
# !
#     INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
# !
# ! !REVISION HISTORY:
# !  10 Mar 2003 - R. Yantosca - Initial version
# !  See https://github.com/geoschem/geos-chem for complete history
# !EOP
# !------------------------------------------------------------------------------
# !BOC
function init_transport!()::Nothing
    # !
    # ! !LOCAL VARIABLES:
    # !
    #     ! Scalars
    #     LOGICAL            :: LTRAN
    #     INTEGER            :: J, K, L, N_DYN
    #     REAL(fp)           :: YMID_R(State_Grid%NY)
    #     REAL(fp)           :: REAL_N_DYN

    #     ! Strings
    #     CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    #     !=================================================================
    #     ! Initialize
    #     !=================================================================
    #     RC      = GC_SUCCESS
    #     ErrMsg  = ''
    #     ThisLoc = ' -> at Init_Transport (in module GeosCore/transport_mod.F90)'

    #     !=================================================================
    #     ! Allocate arrays for TPCORE vertical coordinates
    #     !
    #     ! For fvDAS TPCORE with for GEOS-FP or MERRA-2 met fields:
    #     !
    #     !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
    #     !
    #     ! Also here Ap, Bp will be flipped since both TPCORE versions
    #     ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
    #     !=================================================================
    #     ALLOCATE( Ap( State_Grid%NZ+1 ), STAT=RC )
    #     CALL GC_CheckVar( 'transport_mod.F:Ap', 0, RC )
    #     IF ( RC /= GC_SUCCESS ) RETURN

    #     ALLOCATE( Bp( State_Grid%NZ+1 ), STAT=RC )
    #     CALL GC_CheckVar( 'transport_mod.F:Bp', 0, RC )
    #     IF ( RC /= GC_SUCCESS ) RETURN

    #     ! Flip Ap and Bp for TPCORE
    #     DO L = 1, State_Grid%NZ+1

    #        ! As L runs from the surface up,
    #        ! K runs from the top down
    #        K = ( State_Grid%NZ + 1 ) - L + 1

    #        Ap(L) = GET_AP(K)          ! Ap(L) is in [hPa]
    #        Bp(L) = GET_BP(K)
    #     ENDDO

    #     !=================================================================
    #     ! Allocate arrays for surface area and layer thickness
    #     !=================================================================
    #     ALLOCATE( A_M2( State_Grid%NY ), STAT=RC )
    #     CALL GC_CheckVar( 'transport_mod.F:A_m2', 0, RC )
    #     IF ( RC /= GC_SUCCESS ) RETURN

    #     ! Surface area [m2]
    #     DO J = 1, State_Grid%NY
    #        A_M2(J) = State_Grid%Area_M2(1,J)
    #     ENDDO

    #     !=================================================================
    #     ! Additional setup for the fvDAS version of TPCORE
    #     !=================================================================

    #     ! Initialize
    #     N_DYN = GET_TS_DYN()
    #     N_ADJ = 0
    #     NG    = 0
    #     MG    = 0

    #     ! YMID_R is latitude of grid box center [radian]
    #     DO J = 1,State_Grid%NY
    #        YMID_R(J) = State_Grid%YMid_R(1,J)
    #     ENDDO

    #     REAL_N_DYN = N_DYN

    #     ! Call INIT routine from "tpcore_fvdas_mod.f"
    #     CALL INIT_TPCORE( State_Grid%NX, State_Grid%NY,  State_Grid%NZ, &
    #                       JFIRST, JLAST, NG, MG,         REAL_N_DYN,    &
    #                       Re,    YMID_R, RC                            )

    #     ! Trap potential errors
    #     IF ( RC /= GC_SUCCESS ) THEN
    #        ErrMsg = 'Error encountered in "Init_Tpcore"!'
    #        CALL GC_Error( ErrMsg, RC, ThisLoc )
    #        RETURN
    #     ENDIF

    #   END SUBROUTINE INIT_TRANSPORT
    # !EOC
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
    input_opt::optinput,
    state_grid::grdstate,
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

    #     !=================================================================
    #     ! Initialize
    #     !=================================================================

    # TODO: Check if this is correct
    # Assume success
    # RC = GC_SUCCESS

    # Copy values from Input_Opt
    ltran = input_opt.ltran

    #     !=================================================================
    #     ! Allocate arrays for TPCORE vertical coordinates
    #     ! GEOS-FP/MERRA2 nested grid simulation only!!!
    #     !
    #     ! For fvDAS TPCORE with for GEOS-FP/MERRA2 met fields:
    #     !
    #     !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
    #     !
    #     ! Also here Ap, Bp will be flipped since both TPCORE versions
    #     ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
    #     !=================================================================
    #     ALLOCATE( Ap( State_Grid%NZ+1 ), STAT=RC )

    # Allocate arrays for TPCORE vertical coordinates (GEOS-FP/MERRA2 nested grid simulation only!!!)
    # For fvDAS TPCORE with for GEOS-FP/MERRA2 met fields:
    #    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
    # Also here Ap, Bp will be flipped since both TPCORE versions index levels from the atm. top downwards (bdf, bmy, 10/30/07)
    global ap = zeros(state_grid.nz + 1)
    global bp = zeros(state_grid.nz + 1)

    for l in 1:state_grid.nz+1
        k = state_grid.nz + 1 - l + 1
        ap[l] = get_ap(k)
        bp[l] = get_bp(k)
    end

    #     !=================================================================
    #     ! Allocate arrays for surface area and layer thickness
    #     !=================================================================
    #     ALLOCATE( A_M2( State_Grid%NY ), STAT=RC )
    #     IF ( RC /= 0 ) CALL ALLOC_ERR( 'A_M2' )

    #     ! Surface area [m2]
    #     DO J = 1, State_Grid%NY
    #        A_M2(J) = State_Grid%Area_M2(1,J)
    #     ENDDO

    #     !=================================================================
    #     ! Additional setup for the fvDAS version of TPCORE
    #     !=================================================================

    #     ! Initialize
    #     N_DYN  = GET_TS_DYN()
    #     N_ADJ  = 0
    #     NG     = 0
    #     MG     = 0

    #     ! (lzh, 4/1/2015)
    #     BUFF_SIZE = 2
    #     IM_W1       =  ( State_Grid%NX - State_Grid%WestBuffer - &
    #                      State_Grid%EastBuffer  ) + 2 * BUFF_SIZE
    #     JM_W1       =  ( State_Grid%NY - State_Grid%SouthBuffer - &
    #                      State_Grid%NorthBuffer ) + 2 * BUFF_SIZE
    #     I0_W1     = State_Grid%WestBuffer  - BUFF_SIZE
    #     J0_W1     = State_Grid%SouthBuffer - BUFF_SIZE

    #     ! YMID_R is latitude of grid box center [radians]
    #     DO J = 1, State_Grid%NY
    #        YMID_R_W(J) = State_Grid%YMid_R(1,J)
    #     ENDDO

    #     ! Compute YMID_R_W at southern edge of nested region
    #     J = 0
    #     YMID_R_W(J) = State_Grid%YMid_R(1,J+1) - (State_Grid%DY * PI_180)

    #     ! Compute YMID_R_W at northern edge of nested region
    #     J = State_Grid%NY+1
    #     YMID_R_W(J) = State_Grid%YMid_R(1,J-1) + (State_Grid%DY * PI_180)

    #     ! Call INIT routine from "tpcore_window_mod.F90"
    #     CALL INIT_WINDOW(                                                        &
    #          State_Grid = State_Grid,                                            &
    #          IM         = IM_W1,                                                 &
    #          JM         = JM_W1,                                                 &
    #          KM         = State_Grid%NZ,                                         &
    #          JFIRST     = JFIRST,                                                &
    #          JLAST      = JLAST,                                                 &
    #          NG         = NG,                                                    &
    #          MG         = MG,                                                    &
    #          DT         = REAL( N_DYN, fp ),                                     &
    #          AE         = REAL( Re,    fp ),                                     &
    #          CLAT       = YMID_R_W( J0_W1:(J0_W1+JM_W1+1) )                     )

    #   END SUBROUTINE INIT_WINDOW_TRANSPORT
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

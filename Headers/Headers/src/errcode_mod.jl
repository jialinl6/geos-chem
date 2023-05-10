"""
Contains the error codes (i.e. that report success or failure) returned by GEOS-Chem routines.

## Remarks
The error codes are returned by routines at various levels of GEOS-Chem.

## Revision History
19 Oct 2012 - R. Yantosca - Initial version \\
See https://github.com/geoschem/geos-chem for complete history
"""
module errcode_mod
# !
# ! !USES:
# !
#   IMPLICIT NONE
#   PRIVATE
# !
# ! !PUBLIC MEMBER FUNCTIONS:
# !
#   PUBLIC :: GC_Error
#   PUBLIC :: GC_Warning
#   PUBLIC :: GC_CheckVar

const GC_SUCCESS = 0
const GC_FAILURE = -1


# !------------------------------------------------------------------------------
# !                  GEOS-Chem Global Chemical Transport Model                  !
# !------------------------------------------------------------------------------
# !BOP
# !
# ! !IROUTINE: GC_Error
# !
# ! !DESCRIPTION: Subroutine GC\_Error prints an error message and sets RC to
# !  GC\_FAILURE. Note that this routine does not stop a run, but it will cause
# !  a stop at a higher level if you add a catch for RC /= GC\_SUCCESS.
# !\\
# !\\
# ! !INTERFACE:
# !
#   SUBROUTINE GC_Error( ErrMsg, RC, ThisLoc, Instr )
"""
Prints an error message and sets RC to GC_FAILURE. Note that this routine does not stop a run, but it will cause a stop at a higher level if you add a catch for RC /= GC_SUCCESS.

## Arguments
- `errmsg` (`in`) - Error message to display
- `rc` (`inout`) - Error code
- `thisloc` (`in`, optional) - Location of error
- `instr` (`in`, optional) - Other instructions
"""
function gc_error!(
    errmsg::String,
    rc::Int64,
    thisloc::String="",
    instr::String="",
)::Nothing
    # !
    # ! !USES:
    # !
    #     USE Charpak_Mod,    ONLY : WordWrapPrint
    # #if defined( MODEL_CESM )
    #     USE CAM_ABORTUTILS, ONLY : ENDRUN
    # #endif
    # #if defined( ESMF_ )
    # #include "MAPL_Generic.h"
    #     USE ESMF
    #     USE MAPL_Mod
    # #endif
    # !
    # ! !INPUT PARAMETERS:
    # !
    #     CHARACTER(LEN=*), INTENT(IN   )            :: ErrMsg  ! Message to display
    #     CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: ThisLoc ! Location of error
    #     CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: Instr   ! Other instructions
    # !
    # ! !INPUT/OUTPUT PARAMETERS:
    # !
    #     INTEGER,          INTENT(INOUT)            :: RC      ! Error code
    # !
    # ! !REVISION HISTORY:
    # !  13 Aug 2015 - E. Lundgren - Initial version, based on C. Keller's HCO_ERROR
    # !  See https://github.com/geoschem/geos-chem for complete history
    # !EOP
    # !------------------------------------------------------------------------------
    # !BOC

    #     CHARACTER(LEN=1000) :: Message
    message::String
    # #if defined( ESMF_)
    #     INTEGER             :: localPET, STATUS
    #     CHARACTER(4)        :: localPETchar
    #     CHARACTER(LEN=1023) :: MSG, MSG1, MSG2
    #     TYPE(ESMF_VM)       :: VM
    # #endif

    @static if @isdefined(ESMF_)
        localpet::Int64
        status::Int64
        localpetchar::String
        msg::String
        msg1::String
        msg2::String
        vm::ESMF_VM

        # Get current thread number
        esmf_vmgetcurrent(vm, status)
        esmf_vmget(vm, localpet, __rc__)
        localpetchar = string(localpet)
        msg1 = "GEOS-Chem ERROR [" * strip(localpetchar) * "]: " * strip(errmsg)
        msg2 = ""
        if !empty(strip(thisloc))
            msg2 = "\n --> LOCATION: " * strip(thisloc)
        end
        msg = "\n" * strip(msg1) * strip(msg2)

        # Print error message
        print(msg)
    else
        println("="^79)
        message = "GEOS-Chem ERROR: " * strip(errmsg)
        # TODO: idk where this is defined
        WordWrapPrint!(message, 78)

        if !empty(strip(thisloc))
            message = "ERROR LOCATION: " * strip(thisloc)
            println(message)
        end

        if !empty(strip(instr))
            message = "INSTRUCTIONS: " * strip(instr)
            println(message)
        end

        println("="^79)
        println()
    end

    # Force the message to be flushed to the log file
    flush(stdout)

    @static if @isdefined(MODEL_CESM)
        endrun("GEOS-Chem failure!")
    end

    # Return with failure, but preserve existing error code
    if rc == GC_SUCCESS
        rc = GC_FAILURE
    end
end

#   END SUBROUTINE GC_Error
# !EOC
# !------------------------------------------------------------------------------
# !                  GEOS-Chem Global Chemical Transport Model                  !
# !------------------------------------------------------------------------------
# !BOP
# !
# ! !IROUTINE: GC_Warning
# !
# ! !DESCRIPTION: Subroutine GC\_Warning prints an warning (i.e. non-fatal
# !  error message) and sets RC to GC\_SUCCESS.
# !\\
# !\\
# ! !INTERFACE:
# !
#   SUBROUTINE GC_Warning( WarnMsg, RC, ThisLoc, Instr )
# !
# ! !USES:
# !
#     USE Charpak_Mod, ONLY : WordWrapPrint
# !!
# ! !INPUT PARAMETERS:
# !
#     CHARACTER(LEN=*), INTENT(IN   )            :: WarnMsg ! Message to display
#     CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: ThisLoc ! Location of warning
#     CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: Instr   ! Other instructions
# !
# ! !INPUT/OUTPUT PARAMETERS:
# !
#     INTEGER,          INTENT(INOUT)            :: RC
# !
# ! !REVISION HISTORY:
# !  13 Aug 2015 - E. Lundgren - Initial version, based on C. Keller's HCO_ERROR
# !  See https://github.com/geoschem/geos-chem for complete history
# !EOP
# !------------------------------------------------------------------------------
# !BOC

#     CHARACTER(LEN=1000) :: Message

#     !=======================================================================
#     ! GC_ERROR begins here
#     !=======================================================================

#     ! Separator
#     WRITE( 6, '(a)' ) REPEAT( '=', 79 )

#     ! Print error message to log
#     Message =  'GEOS-Chem WARNING: ' // TRIM( WarnMsg )
#     CALL WordWrapPrint( Message, 78 )

#     ! Print error location to log
#     IF ( PRESENT( ThisLoc ) ) THEN
#        Message = 'WARNING LOCATION: ' // TRIM( ThisLoc )
#        WRITE( 6, '(a)' ) TRIM( ThisLoc )
#     ENDIF

#     ! Print additional instructions to log
#     IF ( PRESENT( Instr ) ) THEN
#        WRITE( 6, '(a)' )
#        CALL WordWrapPrint( Instr, 78 )
#     ENDIF

#     ! Separators
#     WRITE( 6, '(a)' ) REPEAT( '=', 79 )
#     WRITE( 6, '(a)' ) ''

#     ! Force the message to be flushed to the log file
#     CALL Flush( 6 )

#     ! Return with success, since this is only a warning message
#     RC = GC_SUCCESS

#   END SUBROUTINE GC_Warning
# !EOC
# !------------------------------------------------------------------------------
# !                  GEOS-Chem Global Chemical Transport Model                  !
# !------------------------------------------------------------------------------
# !BOP
# !
# ! !IROUTINE: GC_CheckVar
# !
# ! !DESCRIPTION: Wrapper routine for GC\_Error.  Prints an error message
# !  if there is an allocation or registration error.  This is intended to
# !  be called from the state initialization method (e.g. Init\_State\_Met).
# !\\
# !\\
# ! !INTERFACE:
# !
#   SUBROUTINE GC_CheckVar( Variable, Operation, RC )
# !
# ! !INPUT PARAMETERS:
# !
#     CHARACTER(LEN=*), INTENT(IN)    :: Variable   ! Name of variable to check
#     INTEGER,          INTENT(IN)    :: Operation  ! 0=Allocate
#                                                   ! 1=Register
#                                                   ! 2=Deallocate
# !
# ! !OUTPUT PARAMETERS:
# !
#     INTEGER,          INTENT(INOUT) :: RC         ! Success or failure
# !
# ! !REMARKS:
# !  You also need to add an
# !    IF ( RC /= GC_SUCCESS ) RETURN
# !  from the calling routine for proper error handling.
# !
# ! !REVISION HISTORY:
# !  27 Jun 2017 - R. Yantosca - Initial version
# !  See https://github.com/geoschem/geos-chem for complete history
# !EOP
# !------------------------------------------------------------------------------
# !BOC
# !
# ! !LOCAL VARIABLES:
# !
#   ! Strings
#   CHARACTER(LEN=255) :: ErrMsg, ThisLoc

#   !=========================================================================
#   ! Initialize
#   !=========================================================================

#   ! Define error message
#   SELECT CASE( Operation )
#      CASE( 1 )
#         ErrMsg = 'Could not register '   // TRIM( Variable ) // '!'
#      CASE( 2 )
#         ErrMsg = 'Could not deallocate ' // TRIM( Variable ) // '!'
#      CASE DEFAULT
#         ErrMsg = 'Could not allocate '   // TRIM( Variable ) // '!'
#   END SELECT

#   ! Define location string
#   ThisLoc   = ' -> at GC_CheckVar (in Headers/errcode_mod.F90)'

#   !=========================================================================
#   ! Display error message if necessary
#   !=========================================================================
#   IF ( RC /= GC_SUCCESS ) THEN
#      CALL GC_Error( ErrMsg, RC, ThisLoc )
#   ENDIF

#   END SUBROUTINE GC_CheckVar
# !EOC
end
# END MODULE ErrCode_Mod


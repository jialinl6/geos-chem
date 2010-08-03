! $Id: logical_mod.f,v 1.2 2010/03/15 19:33:19 ccarouge Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: logical_mod.f
!
! !DESCRIPTION: Module LOGICAL_MOD contains all of the logical switches 
!  used by GEOS-CHEM.
!\\
!\\
! !INTERFACE:
!
      MODULE LOGICAL_MOD
!
! !USES:
!
      IMPLICIT NONE
!
! !REVISION HISTORY:
!  05 Nov 2004 - R. Yantosca - Added LNEI99 switch to toggle EPA/NEI emissions 
!  20 Dec 2004 - R. Yantosca - Added LAVHRRLAI switch for AVHRR LAI fields 
!  20 Oct 2005 - T-M Fu.     - Added LMEGAN switch for MEGAN biogenics
!  01 Nov 2005 - B. Field    - Added LEMEP switch
!  26 Feb 2006 - R. Yantosca - Added LDYNOCEAN switch for online ocean Hg model
!  05 Apr 2006 - R. Yantosca - Added LGFED2BB switch for GFED2 BIOMASS BURNING
!  05 May 2006 - L. Murray   - Added LCTH, LMFLUX, LPRECON for lightning
!  30 May 2006 - S. Wu       - Added LFUTURE  
!  26 Jun 2006 - R. Park     - Added LBRAVO  
!  06 Jul 2006 - Aaron van D.- Added LEDGAR, LEDGARNOx, LEDGARCO, LEDGARSHIP, 
!                              LEDGARSOx switches for EDGAR emissions
!  17 Aug 2006 - R. Yantosca - Added LSTREETS for David Streets' emissions
!  21 Aug 2006 - P. Le Sager - Added LVARTROP for variable tropopause
!  31 Jan 2007 - L. Murray   - Added LOTDREG, LOTDLOC for regional or local 
!                              OTD-LIS redistribution of lightning flashes
!  31 Jan 2007 - L. Murray   - Added LOTDSCALE
!  08 Mar 2008 - Aaron van D.- Added LCAC, LARCSHIP, LEMEPSHIP
!  24 Nov 2008 - Aaron van D.- Added LVISTAS
!  16 Oct 2009 - Y. Chen     - Added L8DAYBB, L3HRBB and LSYNOPBB for 
!                              8-day and 3-hr GFED BB emissions
!  26 Jan 2009 - P. Le Sager - Added LICARTT to account for Hudman 
!                              corrections to EPA/NEI99
!  12 Feb 2009 - D. Henze    - Added LSVCSPEC 
!  10 Mar 2009 - T-M Fu      - Added LMEGANMONO
!  10 Mar 2009 - T-M Fu      - Added LDICARB
!  29 May 2009 - J. Lin      - Add LNLPBL, LARPBLH and LDEPBCK (non-local PBL)
!  18 May 2009 - P. Le Sager - Added LCOOKE
!  28 May 2009 - P. Le Sager - Added LKPP 
!  16 Oct 2009 - C. Lee      - Added LICOADSSHIP
!  18 Aug 2009 - K. Wecht    - Added switches for CH4 emissions & budget  
!  16 Oct 2009 - R. Yantosca - Added LLINOZ switch for Linoz O3 strat chem
!  16 Oct 2009 - R. Yantosca - Added ProTeX header
!  30 Oct 2009 - Aaron van D - Added LNEI2005
!  19 Nov 2009 - M. Barkley  - Added LMODISLAI and LPECCA
!  18 Dec 2009 - Aaron van D - Added HDF5 logical switches
!  18 Dec 2009 - Aaron van D - Added logicals for NA, EU, CH, CU nested grids
!  18 Dec 2009 - Aaron van D - Added logical for 2 x 2.5 TPCORE BC's 
!  22 Jan 2010 - R. Yantosca - Added LTOMAS switch
!  29 Jan 2009 - F. Paulot   - Added LFERTILIZERNOX.
!  26 Feb 2010 - R. Yantosca - Remove obsolete LEMBED flag
!  18 May 2010 - R. Nassar   - Add logical flags for CO2 offline simulation
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! The logical switches below are for turning options ON or OFF
      !=================================================================

      !%%%% TOMAS aerosol microphysics %%%%
      LOGICAL :: LTOMAS

      !%%%% Aerosols %%%%
      LOGICAL :: LATEQ           ! ??     
      LOGICAL :: LCARB           ! Use carbon aerosol tracers?
      LOGICAL :: LCRYST          ! Use Crystalline aerosols? (not implemented)
      LOGICAL :: LCOOKE          ! Use the C
      LOGICAL :: LDEAD           ! Use the DEAD/Zender dust algorithm?
      LOGICAL :: LDUST           ! Use dust aerosol tracers?
      LOGICAL :: LSULF           ! Use sulfate aerosol tracers?
      LOGICAL :: LSOA            ! Use secondary organic aerosol tracers?
      LOGICAL :: LSSALT          ! Use sea-salt aerosol tracers?
      LOGICAL :: LDICARB         ! Use dicarbonyl chemistry mechanism?

      !%%%% Chemistry %%%%
      LOGICAL :: LCHEM           ! Use chemistry?
      LOGICAL :: LKPP            ! Use KPP solver instead of SMVGEAR?

      !%%%% Cloud convection %%%%
      LOGICAL :: LCONV           ! Use cloud convection?

      !%%%% Diagnostics %%%%
      LOGICAL :: LDBUG
      LOGICAL :: LDIAG           ! Use diagnostics?
      LOGICAL :: LPRT            ! Save debug output to log file?
      LOGICAL :: LSTDRUN         ! Save out Ox.mass files for benchmarks

      !%%%% Deposition %%%%
      LOGICAL :: LDRYD           ! Use dry deposition?
      LOGICAL :: LWETD           ! Use wet deposition?

      !%%%% Emissions %%%%
      LOGICAL :: LAIRNOX         ! Use aircraft NOx emissions?
      LOGICAL :: LANTHRO         ! Turns all anthropogenic emissions on/off
      LOGICAL :: LBBSEA          ! Use seasonal biomass burning emissions?
      LOGICAL :: LBIONOX         ! <-- deprecated: replace w/ LBIOMASS soon
      LOGICAL :: LBIOMASS        ! Use biomass emissions?
      LOGICAL :: LBIOFUEL        ! Use biofuel emissions?
      LOGICAL :: LBIOGENIC       ! Use biogenic emissions?
      LOGICAL :: LCAC            ! Use CAC Canadian regional emissions
      LOGICAL :: LBRAVO          ! Use BRAVO Mexican regional emissions?
      LOGICAL :: LEDGAR          ! Use EDGAR emissions?
      LOGICAL :: LEDGARNOx       ! Use EDGAR NOx emissions?
      LOGICAL :: LEDGARCO        ! Use EDGAR CO emissions?
      LOGICAL :: LEDGARSHIP      ! Use EDGAR ship emissions?
      LOGICAL :: LEDGARSOx       ! Use EDGAR SOx emissions?
      LOGICAL :: LEMEP           ! Use EMEP European regional emissions?
      LOGICAL :: LEMIS           ! Use emissions in GEOS-Chem (main switch)
      LOGICAL :: LFFNOX          ! Use anthropogenic NOx emissions
      LOGICAL :: LFOSSIL         ! <-- deprecated: replace w/ LANTHRO soon
      LOGICAL :: LSTREETS        ! Use David Streets SE Asian emissions?
      LOGICAL :: LICARTT         ! Use ICARTT fix to EPA emissions?
      LOGICAL :: LICOADSSHIP     ! Use ICOADS ship emissions inventory
      LOGICAL :: LLIGHTNOX       ! Use lightning NOx emissions?
      LOGICAL :: LOTDREG         ! Use OTD-LIS regional flash redistribution?
      LOGICAL :: LOTDLOC         ! Use OTD-LIS local flash redistribution?
      LOGICAL :: LOTDSCALE       ! Scale LNOx totals to OTD-LIS climatology
      LOGICAL :: LCTH            ! Use cloud-top height option in lightning?
      LOGICAL :: LMFLUX          ! Use mass-flux option in lightning?
      LOGICAL :: LPRECON         ! Use precip option in lightning?
      LOGICAL :: LMEGAN          ! Use MEGAN biogenic emissions?
      LOGICAL :: LMEGANMONO      ! Use MEGAN monoterpenes?
      LOGICAL :: LMONOT          ! Use old 
      LOGICAL :: LNEI99          ! Use EPA 1999 regional emissions?
      LOGICAL :: LNEI05          ! Use EPA 2005 regional emissions?
      LOGICAL :: LSHIPSO2        ! Use SO2 from ship emissions?
      LOGICAL :: LSOILNOX        ! Use soil NOx emissions
      !FP_ISOP (6/2009) LFERTILIZERNOX
      LOGICAL :: LFERTILIZERNOX  ! Use fertilizer NOx emissions
      LOGICAL :: LTOMSAI         ! Scale biomass burning to TOMS AI index?
      LOGICAL :: LWOODCO         ! <-- deprecated: replace w/ LBIOFUEL soon
      LOGICAL :: LAVHRRLAI       ! Use AVHRR leaf-area-indices?
      LOGICAL :: LGFED2BB        ! Use GFED2 biomass burning?
      LOGICAL :: LFUTURE         ! Use future-years emission scaling (GCAP)?
      LOGICAL :: LARCSHIP        ! Use ARCTAS ship emissions?
      LOGICAL :: LEMEPSHIP       ! Use EMEP ship emissions?
      LOGICAL :: LVISTAS         ! Use VISTAS NOx emissions?
      LOGICAL :: L8DAYBB         ! Use GFED2 8-day biomass burning?
      LOGICAL :: L3HRBB          ! Use GFED2 3-hr biomass burning?
      LOGICAL :: LSYNOPBB        ! Use GFED2 synoptic biomass burning
      LOGICAL :: LMODISLAI       ! MODIS LAI (mpb, 2009)
      LOGICAL :: LPECCA          ! PCEEA BVOC emission model (mpb,2009)

      !%%%% Transport and strat BC's %%%%
      LOGICAL :: LFILL           ! Fill negative values in TPCORE?
      LOGICAL :: LMFCT           ! Turns TPCORE MFCT option on/off
      LOGICAL :: LTRAN           ! Turns advection on/off
      LOGICAL :: LTPFV           ! Are we using the GEOS-4 
      LOGICAL :: LUPBD           ! Use stratospheric O3, NOY bdry conditions
      LOGICAL :: LLINOZ          ! Use LINOZ chemistry in the stratosphere?
      LOGICAL :: LWINDO          ! Use nested-grid simulation?
      LOGICAL :: LWINDO2x25      ! Use nested-grid BC's @ 2 x 2.5 resolution?
      LOGICAL :: LWINDO_NA       
      LOGICAL :: LWINDO_EU
      LOGICAL :: LWINDO_CH
      LOGICAL :: LWINDO_CU

      !%%%% Met fields %%%%
      LOGICAL :: LUNZIP          ! Unzip met fields on-the-fly?
      LOGICAL :: LWAIT           ! Wait for met fields to be unzipped?

      !%%%% PBL mixing %%%%
      LOGICAL :: LTURB           ! Use PBL mixing?
      LOGICAL :: LNLPBL          ! Use non-local PBL scheme instead of default?
      LOGICAL :: LARPBLH         ! --> pblh_ar in vdiff_mod.f
      LOGICAL :: LDEPBCK         ! --> drydep_back_cons in vdiff_mod.f

      !%%%% Restart files %%%%
      LOGICAL :: LSVGLB          ! Save tracer restart file?
      LOGICAL :: LSVCSPEC        ! Save CSPEC chemical species restart file?
  
      !%%%% Tagged simulations %%%%
      LOGICAL :: LSPLIT          ! Are we using tagged tracers?

      !%%%% Variable Tropopause %%%%
      LOGICAL :: LVARTROP        ! Use dynamic tropopause option?

      !%%%% Dynamic ocean Hg model %%%%
      LOGICAL :: LDYNOCEAN       ! Use dynamic ocean Hg model?

      !%%%% For the CH4 offline simulation only %%%%
      LOGICAL :: LGAO            ! Use gas & oil emissions?
      LOGICAL :: LCOL            ! Use coal emissions?
      LOGICAL :: LLIV            ! Use livestock emissions?
      LOGICAL :: LWAST           ! Use waste emissions?
      LOGICAL :: LRICE           ! Use rice emissions?
      LOGICAL :: LOTANT          ! Use other anthropogenic emissions?
      LOGICAL :: LWETL           ! Use wetland emissions?
      LOGICAL :: LSOABS          ! Use soil absorption?
      LOGICAL :: LOTNAT          ! Use other natural emissions?
      LOGICAL :: LBFCH4          ! Use CH4 biofuel emissions?
      LOGICAL :: LBMCH4          ! Use CH4 biomass emissions?
      LOGICAL :: LCH4BUD         ! Use computing CH4 budget

      !%%%% For the CO2 offline simulation only %%%%
      LOGICAL :: LGENFF      
      LOGICAL :: LANNFF      
      LOGICAL :: LMONFF      
      LOGICAL :: LSEASBB
      LOGICAL :: LBIONETORIG  
      LOGICAL :: LBIONETCLIM  
      LOGICAL :: LBIODAILY
      LOGICAL :: LBIODIURNAL    
      LOGICAL :: LOCEAN     
      LOGICAL :: LFFBKGRD
      LOGICAL :: LBIOSPHTAG     
      LOGICAL :: LFOSSILTAG     
      LOGICAL :: LOCN1997
      LOGICAL :: LOCN2009ANN
      LOGICAL :: LOCN2009MON
      LOGICAL :: LSHIPEDG
      LOGICAL :: LSHIPICO
      LOGICAL :: LSHIPSCALE
      LOGICAL :: LSHIPTAG
      LOGICAL :: LPLANE
      LOGICAL :: LPLANESCALE
      LOGICAL :: LPLANETAG
      LOGICAL :: LCHEMCO2

      !%%%% HDF5 output %%%%
      LOGICAL :: LND50_HDF       ! Save ND50  diagnostic in HDF5?
      LOGICAL :: LND51_HDF       ! Save ND51  diagnostic in HDF5?
      LOGICAL :: LND51b_HDF      ! Save ND51b diagnostic in HDF5?

      END MODULE LOGICAL_MOD
!EOC

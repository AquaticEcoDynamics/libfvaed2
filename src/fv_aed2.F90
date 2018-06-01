!###############################################################################
!#                                                                             #
!# fv_aed2.F90                                                                 #
!#                                                                             #
!# Interface for FV (Finite Volume) Model to AED2 modules (libaed2).           #
!#   Designed for TUFLOW-FV, released by BMT-WBM:                              #
!#   http://www.tuflow.com/Tuflow%20FV.aspx                                    #
!#                                                                             #
!# This is the main interface module that manages the connection with the      #
!# host hydrodynamic model; done through the 4 PUBLIC functions listed below.  #
!#                                                                             #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created Aug 2013                                                            #
!# Update Mar 2016: Add new env variables and feedback links                   #
!# Update Jun 2016: Add riparian exchange functionality                        #
!# Update Aug 2016: Add links for getting wave stress                          #
!# Update Oct 2016: Add variable mbility support                               #
!#                                                                             #
!###############################################################################

#include "aed2.h"

#define FV_AED_VERS "1.0.0dev"

#ifndef DEBUG
#define DEBUG      0
#endif
#define _NO_ODE_   1

!###############################################################################
MODULE fv_aed2
!-------------------------------------------------------------------------------
   USE aed2_common
   USE fv_zones

   IMPLICIT NONE

   PUBLIC init_aed2_models, &
          init_var_aed2_models, &
          set_env_aed2_models, &
          do_aed2_models, &
          clean_aed2_models


   !#--------------------------------------------------------------------------#
   !# Module Data

   AED_REAL :: Kw, Ksed
   INTEGER  :: solution_method

   !# Main arrays storing/pointing to the state and diagnostic variables
   AED_REAL,DIMENSION(:,:),POINTER :: cc, cc_diag

   !# Array pointing to the lagrangian particle masses and diagnostic properties
   AED_REAL,DIMENSION(:,:),POINTER :: pp, pp_diag

   !# Name of files being used to load initial values for benthic
   !  or benthic_diag vars, and the horizontal routing table for riparian flows
   CHARACTER(len=128) :: init_values_file = ''
   CHARACTER(len=128) :: route_table_file = ''

   !# Maps of surface, bottom and wet/dry (active) cells
   INTEGER,DIMENSION(:),POINTER :: surf_map, benth_map
   LOGICAL,DIMENSION(:),POINTER :: active

   !# Maps to nearest cell with water (for riparian exchange)
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: nearest_active
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: nearest_depth
   INTEGER, DIMENSION(:), ALLOCATABLE       :: route_table

   !# Arrays for work, vertical movement (ws), and cross-boundary fluxes
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: flux, ws
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: total
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: Fsed_setl
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: min_, max_

   !# Arrays for environmental variables (used if they are not supplied externally)
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: nir
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: par
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: uva
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: uvb
   AED_REAL,DIMENSION(:),POINTER :: lpar
   AED_REAL,TARGET :: col_taub  ! a temp var for bottom stress (computed from ustar_bed)

   !# External variables
   AED_REAL :: dt
   INTEGER, DIMENSION(:,:),POINTER :: mat
   AED_REAL,DIMENSION(:,:),POINTER :: rad
   AED_REAL,DIMENSION(:),  POINTER :: temp, salt, rho, nuh, h, z
   AED_REAL,DIMENSION(:),  POINTER :: extcoeff, tss, bio_drag
   AED_REAL,DIMENSION(:),  POINTER :: I_0, wnd, air_temp, rain
   AED_REAL,DIMENSION(:),  POINTER :: area, bathy, shadefrac, rainloss
   AED_REAL,DIMENSION(:),  POINTER :: ustar_bed
   AED_REAL,DIMENSION(:),  POINTER :: wv_uorb, wv_t

   !# maximum single precision real is 2**128 = 3.4e38
   AED_REAL :: glob_min = -1.0e38
   AED_REAL :: glob_max =  1.0e38

   !# Misc variables/options
   AED_REAL :: min_water_depth =  0.0401
   AED_REAL :: wave_factor =  1.0
   LOGICAL  :: old_zones = .TRUE.
   LOGICAL  :: do_zone_averaging = .FALSE.
   LOGICAL  :: request_nearest = .FALSE.
   LOGICAL  :: reinited = .FALSE.
   INTEGER  :: ThisStep = 0, n_equil_substep = 1

   !# Switches for configuring model operation and active links with the host model
   LOGICAL :: link_water_clarity, link_water_density, &
              link_bottom_drag, link_surface_drag, &
              link_ext_par, ext_tss_extinction, &
              link_solar_shade, link_rain_loss, &
              do_limiter, no_glob_lim, do_2d_atm_flux, do_particle_bgc, &
              link_wave_stress, display_minmax

   !# Integers storing number of variables being simulated
   INTEGER :: n_aed2_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet


CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE init_aed2_models(namlst,dname,nwq_var,nben_var,ndiag_var,names,bennames,diagnames)
!-------------------------------------------------------------------------------
! This routine is called by TuflowFV to define numbers and names of variables.
! TuflowFV will allocate the variables after return from this routine.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,          INTENT(in)  :: namlst
   INTEGER,          INTENT(out) :: nwq_var,nben_var,ndiag_var
   CHARACTER(len=*), INTENT(in)  :: dname
   CHARACTER(len=30),ALLOCATABLE,INTENT(out) :: names(:)
   CHARACTER(len=30),ALLOCATABLE,INTENT(out) :: bennames(:)
   CHARACTER(len=30),ALLOCATABLE,INTENT(out) :: diagnames(:)
!
!LOCALS
   TYPE(aed2_variable_t),POINTER :: tvar
   CHARACTER(len=128)            :: aed2_nml_file
   CHARACTER(len=128)            :: tname
   AED_REAL                      :: base_par_extinction, tss_par_extinction
   INTEGER                       :: status, n_sd, i, j

   CHARACTER(len=64) :: models(64)
   NAMELIST /aed2_models/ models

   NAMELIST /aed2_bio/ solution_method, link_bottom_drag,                      &
                       link_surface_drag, link_water_density,                  &
                       link_water_clarity, aed2_nml_file,                      &
                       link_ext_par, base_par_extinction,                      &
                       ext_tss_extinction, tss_par_extinction,                 &
                       do_particle_bgc, do_2d_atm_flux, do_zone_averaging,     &
                       link_solar_shade, link_rain_loss, init_values_file,     &
                       do_limiter, glob_min, glob_max, no_glob_lim,            &
                       route_table_file, n_equil_substep, min_water_depth,     &
                       link_wave_stress, wave_factor, display_minmax
!
!-------------------------------------------------------------------------------
!BEGIN
   print *, " "
   print *, "    using fv_aed2 version ", TRIM(FV_AED_VERS)

   ! Set default AED2 link options
   aed2_nml_file = 'aed2.nml'
   solution_method = 1
   link_bottom_drag = .false.
   link_surface_drag = .false.
   link_water_density = .false.
   link_water_clarity = .false.
   link_solar_shade = .true.
   link_rain_loss = .false.
   link_ext_par = .false.
   base_par_extinction = 0.1
   ext_tss_extinction = .false.
   tss_par_extinction = 0.02
   do_2d_atm_flux = .TRUE.
   do_limiter = .false.
   no_glob_lim = .false.
   do_particle_bgc = .false.
   min_water_depth = 0.0401
   link_wave_stress = .false.
   display_minmax = .false.

   ! Process input file (aed2.nml) to get run options
   print *, "    initialise aed2_core "
   IF ( aed2_init_core(dname) /= 0 ) STOP "Initialisation of aed2_core failed"
   CALL aed2_print_version
   tname = TRIM(dname)//'aed2.nml'
   print *,"    reading fv_aed2 config from ",TRIM(tname)
   OPEN(namlst,file=tname,action='read',status='old',iostat=status)
   IF ( status /= 0 ) CALL STOPIT("Cannot open file " // TRIM(tname))
   READ(namlst,nml=aed2_bio,iostat=status)
   IF ( status /= 0 ) STOP "Cannot read namelist entry aed2_bio"
   Kw = base_par_extinction
   Ksed = tss_par_extinction
   IF ( do_zone_averaging ) old_zones = .FALSE.
   print *,'    link options configured between TFV & AED2 - '
   print *,'        link_ext_par       :  ',link_ext_par
   print *,'        link_water_clarity :  ',link_water_clarity
   print *,'        link_surface_drag  :  ',link_surface_drag,' (not active)'
   print *,'        link_bottom_drag   :  ',link_bottom_drag
   print *,'        link_wave_stress   :  ',link_wave_stress
   print *,'        link_solar_shade   :  ',link_solar_shade
   print *,'        link_rain_loss     :  ',link_rain_loss
   print *,'        link_particle_bgc  :  ',do_particle_bgc,' (not active)'
   print *,'        link_water_density :  ',link_water_density,' (not active)'

   ! Process input file (aed2.nml) to get selected models
   print *,"    reading aed2_models config from ",TRIM(tname)
   models = ''
   READ(namlst, nml=aed2_models, iostat=status)
   IF ( status /= 0 ) STOP "Cannot read namelist entry aed2_models"

   ! Process each model define/setup block
   DO i=1,size(models)
      IF (models(i)=='') EXIT
      CALL aed2_define_model(models(i), namlst)
   ENDDO

   ! Set number of configured variables
   n_aed2_vars = aed2_core_status(nwq_var, nben_var, ndiag_var, n_sd)
   ndiag_var = ndiag_var + n_sd
   n_vars = nwq_var
   n_vars_ben = nben_var
   n_vars_diag = ndiag_var
   n_vars_diag_sheet = n_sd

#if DEBUG
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tvar) ) THEN
            print *,"AED2 var ", i, tvar%sheet, tvar%diag, tvar%extern, TRIM(tvar%name)
         ELSE
            print *,"AED2 var ", i, " is empty"
         ENDIF
      ENDDO
#endif
#if DEBUG
   print*,'    init_aed2_models : n_aed2_vars = ',n_aed2_vars,&
          ' nwq_var = ',nwq_var,' nben_var ',nben_var
#endif

   CALL check_data

   !# names = grab the names from info
   ALLOCATE(names(1:nwq_var),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): ERROR allocating (names)'
   ALLOCATE(bennames(1:nben_var),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): ERROR allocating (bennames)'
   IF ( .NOT. ALLOCATED(diagnames) ) ALLOCATE(diagnames(ndiag_var))
   IF (status /= 0) STOP 'allocate_memory(): ERROR allocating (diagnames)'

   ALLOCATE(min_(1:nwq_var+nben_var)) ; ALLOCATE(max_(1:nwq_var+nben_var))

   print *,"    configured variables - "
   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( .NOT. (tvar%sheet .OR. tvar%diag .OR. tvar%extern) ) THEN
            j = j + 1
            names(j) = TRIM(tvar%name)
            min_(j) = tvar%minimum
            max_(j) = tvar%maximum
            print *,"     S(",j,") AED2 pelagic(3D) variable: ", TRIM(names(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( tvar%sheet .AND. .NOT. (tvar%diag .OR. tvar%extern) ) THEN
            j = j + 1
            bennames(j) = TRIM(tvar%name)
            min_(nwq_var+j) = tvar%minimum
            max_(nwq_var+j) = tvar%maximum
            print *,"     B(",j,") AED2 benthic(2D) variable: ", TRIM(bennames(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( tvar%diag ) THEN
            j = j + 1
            diagnames(j) = TRIM(tvar%name)
            print *,"     D(",j,") AED2 diagnostic variable:  ", TRIM(diagnames(j))
         ENDIF
      ENDIF
   ENDDO

   CLOSE(namlst)
END SUBROUTINE init_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE init_var_aed2_models(nCells, cc_, cc_diag_, nwq, nwqben, sm, bm)
!-------------------------------------------------------------------------------
! Points the AED2 main variable arrays to those provided by the host model.
! At this point TuflowFV should have allocated the variable space.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                         :: nCells
   AED_REAL,POINTER,DIMENSION(:,:),INTENT(in) :: cc_, cc_diag_
   INTEGER,INTENT(inout)                      :: nwq, nwqben
   INTEGER,POINTER,DIMENSION(:),INTENT(in)    :: sm, bm
!
!LOCALS
   INTEGER :: rc, av, v, sv, d, sd
   TYPE(aed2_variable_t),POINTER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   nwq = n_vars
   nwqben = n_vars_ben

   print *,'    init_var_aed2_models : nwq = ',nwq,' nwqben = ',nwqben

   cc => cc_
   cc_diag => cc_diag_
   surf_map => sm
   benth_map => bm

   ! Allocate state variable array
   IF ( .NOT. ASSOCIATED(cc) ) STOP ' ERROR : no association for (cc)'
   cc = 0.

   IF (.not. ASSOCIATED(cc_diag) ) STOP ' ERROR : no association for (cc_diag)'
   cc_diag = 0.

   ! Allocate array with vertical movement rates (m/s, positive for upwards)
   ! These will be set the values provided by the modules
   ALLOCATE(ws(1:nCells,1:n_aed2_vars),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (ws)'
   ws = 0.

   !!# place holder for lagranigan particles
   !IF(do_particle_bgc) THEN
   !  pp => pp_
   !END IF

   ! Allocate array for photosynthetically active radiation (PAR).
   ! This will be calculated internally during each time step.
   ALLOCATE(par(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (par)'
   par = 0.
   ALLOCATE(nir(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (nir)'
   nir = 0.
   ALLOCATE(uva(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (uva)'
   uva = 0.
   ALLOCATE(uvb(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (uvb)'
   uvb = 0.

   !# Allocate array for sedimentation fluxes and initialize these to zero (no flux).
   ALLOCATE(Fsed_setl(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (Fsed_setl)'
   Fsed_setl = 0.

   !# Now set initial values
   v = 0 ; sv = 0;
   DO av=1,n_aed2_vars
      IF ( .NOT.  aed2_get_var(av, tv) ) STOP "ERROR getting variable info"
      IF ( .NOT. ( tv%extern .OR. tv%diag) ) THEN  !# neither global nor diagnostic variable
         IF ( tv%sheet ) THEN
            sv = sv + 1
            cc(n_vars+sv, :) = tv%initial
         ELSE
            v = v + 1
            cc(v,:) = tv%initial
         ENDIF
      ENDIF
   ENDDO

   IF ( init_values_file /= '' ) CALL set_initial_from_file
   IF ( route_table_file /= '' ) CALL load_route_table(ubound(bm, 1))

   ALLOCATE(flux(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (flux)'
#if !_NO_ODE_
   ALLOCATE(flux2(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (flux2)'
   ALLOCATE(flux3(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (flux3)'
   ALLOCATE(flux4(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (flux4)'
   ALLOCATE(cc1(n_vars, nCells),stat=rc)   ; IF (rc /= 0) STOP 'allocate_memory(): ERROR allocating (cc1)'
#endif
!
!-------------------------------------------------------------------------------
CONTAINS

   !############################################################################
   CHARACTER FUNCTION tolower(c)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      CHARACTER, INTENT(in) :: c
   !LOCALS
      INTEGER :: ic
   !BEGIN
   !----------------------------------------------------------------------------
      ic = ichar(c)
      if (ic >= 65 .and. ic < 90) ic = (ic+32)
      tolower = char(ic)
   END FUNCTION tolower
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !############################################################################
   FUNCTION same_str_icase(a, b) RESULT(res)
   !----------------------------------------------------------------------------
   !ARGUMENTS
      CHARACTER(len=*), INTENT(in) :: a,b
   !LOCALS
      INTEGER :: len, i
      LOGICAL :: res
   !
   !BEGIN
   !----------------------------------------------------------------------------
      res = .FALSE.
      len = LEN_TRIM(a)
      IF ( len /= LEN_TRIM(b) ) RETURN
      DO i=1, len
         if (tolower(a(i:i)) /= tolower(b(i:i)) ) RETURN
      ENDDO
      res = .TRUE.
   END FUNCTION same_str_icase
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !############################################################################
   SUBROUTINE set_initial_from_file
   !----------------------------------------------------------------------------
   USE aed2_csv_reader
   !
   !LOCALS
      INTEGER :: unit, nccols, ccol
      CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
      TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
      INTEGER :: idx_col = 0, numv = 0, numd = 0, t
      INTEGER,DIMENSION(:),ALLOCATABLE :: vars, vmap
      INTEGER,DIMENSION(:),ALLOCATABLE :: dvar, dmap
      LOGICAL,DIMENSION(:),ALLOCATABLE :: vsheet, dsheet
      LOGICAL :: meh
   !
   !BEGIN
   !----------------------------------------------------------------------------
      unit = aed_csv_read_header(init_values_file, csvnames, nccols)
      IF (unit <= 0) RETURN !# No file found
      print *,'    spatial AED2 var initialisation from file: '
      print *,'        ', TRIM(init_values_file)

      DO ccol=1,nccols
         IF ( csvnames(ccol) == "ID" ) THEN
            idx_col = ccol
            EXIT
         ENDIF
      ENDDO

      ALLOCATE(vars(nccols))   ; ALLOCATE(vmap(nccols))
      ALLOCATE(dvar(nccols))   ; ALLOCATE(dmap(nccols))
      ALLOCATE(vsheet(nccols)) ; ALLOCATE(dsheet(nccols))
      ALLOCATE(values(nccols))
      vmap = 0 ; dmap = 0

      IF ( idx_col > 0 ) THEN
         v = 0 ; sv = 0; d = 0; sd = 0
         DO av=1,n_aed2_vars
            IF ( .NOT. aed2_get_var(av, tv) ) STOP "ERROR getting variable info"
            IF ( .NOT. ( tv%extern ) ) THEN  !#  dont do environment vars
               IF (tv%diag) THEN
                  d = d + 1
               ELSE
                  IF ( tv%sheet ) THEN
                     sv = sv + 1
                  ELSE
                     v = v + 1
                  ENDIF
               ENDIF
               DO ccol=1,nccols
                  IF ( same_str_icase(tv%name, csvnames(ccol)) ) THEN
                     IF (tv%diag) THEN
                        numd = numd + 1
                        dmap(numd) = ccol
                      ! IF ( same_str_icase(tv%name, "LND_phreatic") ) THEN
                      ! phreat_id = av ; phreat_col = ccol ; phreat_var = d ; ENDIF
                        dvar(numd) = d
                        dsheet(numd) = tv%sheet
                     ELSE
                        numv = numv + 1
                        vmap(numv) = ccol
                        IF ( tv%sheet ) THEN
                           vars(numv) = n_vars + sv
                        ELSE
                           vars(numv) = v
                        ENDIF
                        vsheet(numv) = tv%sheet
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         DO WHILE ( aed_csv_read_row(unit, values) )
            t = extract_integer(values(idx_col))
            DO v=1,numv
               IF ( vmap(v) == 0 ) CYCLE
               If ( vsheet(v) ) THEN
                  cc(vars(v), bm(t)) = extract_double(values(vmap(v)))
               ELSE
                  cc(vars(v), sm(t):bm(t)) = extract_double(values(vmap(v)))
               ENDIF
            ENDDO
            DO v=1,numd
               IF ( dmap(v) == 0 ) CYCLE
               ! IF (dmap(v) == phreat_col ) &
               ! print*, " XXX setting phreat_col ", phreat_var
               If ( vsheet(v) ) THEN
                  cc_diag(dvar(v), bm(t)) = extract_double(values(dmap(v)))
               ELSE
                  cc_diag(dvar(v), sm(t):bm(t)) = extract_double(values(dmap(v)))
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      meh = aed_csv_close(unit)
      !# don't care if close fails

      IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
      IF (ALLOCATED(values))    DEALLOCATE(values)
      IF (ALLOCATED(vars))      DEALLOCATE(vars)
      IF (ALLOCATED(vmap))      DEALLOCATE(vmap)
      IF (ALLOCATED(dvar))      DEALLOCATE(dvar)
      IF (ALLOCATED(dmap))      DEALLOCATE(dmap)
   END SUBROUTINE set_initial_from_file
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !############################################################################
   SUBROUTINE load_route_table(nrows)
   !----------------------------------------------------------------------------
   USE aed2_csv_reader
   !ARGUMENTS
      INTEGER,INTENT(in) :: nrows
   !
   !LOCALS
      INTEGER :: unit, nccols, ccol, crow
      CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
      TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
      INTEGER :: idx_col = 0, t
      LOGICAL :: meh
   !
   !BEGIN
   !----------------------------------------------------------------------------
      unit = aed_csv_read_header(route_table_file, csvnames, nccols)
      IF (unit <= 0) RETURN !# No file found
      print *,'    riparian cell routing set from file: '
      print *,'        ', TRIM(route_table_file)
!# The format of the file should be me, "lowest ajoining" - ie always 2 colums
!# and always in the order - and we dont really care about the header, but being
!# csv it should have it so we read but ignore it.
!     DO ccol=1,nccols
!        IF ( csvnames(ccol) == "ID" ) THEN
!           idx_col = ccol
!           EXIT
!        ENDIF
!     ENDDO
!     IF (idx_col == 0) THEN
!        print*,"Could not find column 'ID'"
!        RETURN
!     ENDIF
      idx_col = 1

      ALLOCATE(values(nccols))
      ALLOCATE(route_table(nrows))
      ALLOCATE(nearest_active(nrows))
      ALLOCATE(nearest_depth(nrows))
      route_table = 0

      crow = 0
      DO WHILE ( aed_csv_read_row(unit, values) )
         crow = crow + 1
         IF ( crow > nrows ) THEN
            print*, "        NOTE: routing table has more rows than expected - extras ignored"
         ENDIF
         t = extract_integer(values(idx_col))
         route_table(crow) = extract_integer(values(2))

         !MH PUT A CHECK HERE TO MAKE SURE NO CIRCULAR REFERENCE
      ENDDO

      IF ( crow < nrows ) &
      print*, "        NOTE: routing table has less rows than expected? ",crow,"/",nrows

      meh = aed_csv_close(unit)
      !# don't care if close fails

      IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
      IF (ALLOCATED(values))    DEALLOCATE(values)
   END SUBROUTINE load_route_table
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END SUBROUTINE init_var_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE set_env_aed2_models(dt_,              &
                            ! 3D env variables
                               temp_,            &
                               salt_,            &
                               rho_,             &
                               h_,               &
                               tss_,             &
                               rad_,             &
                            ! 3D feedback arrays
                               extcoeff_,        &
                            ! 2D env variables
                               area_,            &
                               I_0_,             &
                               wnd_,             &
                               rain_,            &
                               air_temp_,        &
                               ustar_bed_,       &
                               ustar_surf_,      &
                               wv_uorb_,         &
                               wv_t_,            &
                               z_,               &
                               bathy_,           &
                               mat_id_,          &
                               active_,          &
                            ! 2D feedback arrays
                               biodrag_,         &
                               solarshade_,      &
                               rainloss_         &
                               )
!-------------------------------------------------------------------------------
! Provide environmental information from TuflowFV and set feedback arrays
!-------------------------------------------------------------------------------
!ARGUMENTS
   DOUBLETYPE, INTENT(in) :: dt_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: temp_, salt_, rho_, h_,    &
                                                    area_, tss_, extcoeff_, z_
   AED_REAL, INTENT(in), DIMENSION(:,:), POINTER :: rad_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: I_0_, wnd_, ustar_bed_, ustar_surf_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: wv_uorb_, wv_t_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: rain_, bathy_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: air_temp_
   INTEGER,  INTENT(in), DIMENSION(:,:), POINTER :: mat_id_
   LOGICAL,  INTENT(in), DIMENSION(:),   POINTER :: active_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: biodrag_, solarshade_, rainloss_
!
!LOCALS
   INTEGER :: i, j
   INTEGER :: nTypes, cType
   INTEGER, DIMENSION(:),ALLOCATABLE :: mat_t
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,'    set_env_aed2_models : linking to host environment vars '

   !# Provide pointers to arrays with environmental variables to AED2.
   dt = dt_

   !# 2D (sheet) variables being pointed to
   area => area_
   I_0 => I_0_
   wnd => wnd_
   ustar_bed => ustar_bed_
   mat => mat_id_
   bathy => bathy_
   rain  => rain_
   shadefrac => solarshade_
   rainloss => rainloss_
   bio_drag => biodrag_
   air_temp => air_temp_
   IF(link_wave_stress)THEN
     wv_uorb => wv_uorb_
     wv_t => wv_t_
   END IF

   !# 3D variables being pointed to
   h => h_               !# layer heights [1d array] needed for advection, diffusion
   z => z_               !# depth [1d array], used to calculate local pressure
   extcoeff => extcoeff_ !# biogeochemical light attenuation coefficients [1d array],
                         !# output of biogeochemistry, input for physics
   salt => salt_
   temp => temp_

   rho => rho_
   tss => tss_
   active => active_

   IF (link_ext_par) lpar => rad_(1,:)

!  ALLOCATE(pactive(size(active)))
!  pactive = active
#if DEBUG
!if ( .not. associated(area) ) print*, " No association for area"
!if ( .not. associated(I_0) ) print*, " No association for I_0"
!if ( .not. associated(wnd) ) print*, " No association for wnd"
!if ( .not. associated(ustar_bed) ) print*, " No association for ustar_bed"
!if ( .not. associated(mat) ) print*, " No association for mat"
!if ( .not. associated(bathy) ) print*, " No association for bathy"
!if ( .not. associated(rain) ) print*, " No association for rain"
!if ( .not. associated(shadefrac) ) print*, " No association for shadefrac"
!if ( .not. associated(rainloss) ) print*, " No association for rainloss"
!if ( .not. associated(bio_drag) ) print*, " No association for bio_drag"
!if ( .not. associated(air_temp) ) print*, " No association for air_temp"
!if ( .not. associated(h) ) print*, " No association for h"
!if ( .not. associated(z) ) print*, " No association for z"
!if ( .not. associated(extcoeff) ) print*, " No association for extcoeff"
!if ( .not. associated(salt) ) print*, " No association for salt"
!if ( .not. associated(temp) ) print*, " No association for temp"
!if ( .not. associated(rho) ) print*, " No association for rho"
!if ( .not. associated(tss) ) print*, " No association for tss"
!if ( .not. associated(active) ) print*, " No association for active"
#endif
! CALL CheckPhreatic

   IF (old_zones) THEN
      !# We allocate the full size array because that's how the indices are presented
      ALLOCATE(zone(ubound(mat_id_,2))) ! JC
      zone = 1.
      DO i=1, ubound(mat_id_,2) ! JC
         !# use the bottom index to fill the array
         zone(i) = mat_id_(1,i)! JC
      ENDDO
   ELSE
      CALL init_zones(ubound(mat_id_, 2), mat_id_, n_aed2_vars, n_vars, n_vars_ben)
   ENDIF

END SUBROUTINE set_env_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! !###############################################################################
! SUBROUTINE CheckPhreatic
! !-------------------------------------------------------------------------------
! !ARGUMENTS
! !
! !LOCALS
!    AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben),       &
!                                             flux_rip(n_vars+n_vars_ben)
!    TYPE (aed2_column_t) :: column(n_aed2_vars)
!    INTEGER :: col, bot
! !
! !BEGIN
!    DO col=1, size(active)
!       bot = benth_map(col)
!       !# set column data structure from global arrays
!       CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben, flux_rip)
!
!       print*," ^^^", column(phreat_id)%cell_sheet, cc_diag(phreat_var, bot)
!    ENDDO
! END SUBROUTINE CheckPhreatic


!###############################################################################
SUBROUTINE check_data
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: av, i, top, bot
   INTEGER :: v, d, sv, sd, ev, err_count
   TYPE(aed2_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   err_count = 0

   DO av=1,n_aed2_vars
      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "ERROR getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; tvar%found = .true.
            CASE ( 'salinity' )    ; tvar%found = .true.
            CASE ( 'density' )     ; tvar%found = .true.
            CASE ( 'layer_ht' )    ; tvar%found = .true.
            CASE ( 'layer_area' )  ; tvar%found = .true.
            CASE ( 'rain' )        ; tvar%found = .true.
            CASE ( 'rainloss' )    ; tvar%found = .true.
            CASE ( 'material' )    ; tvar%found = .true.
            CASE ( 'bathy' )       ; tvar%found = .true.
            CASE ( 'extc_coef' )   ; tvar%found = .true.
            CASE ( 'tss' )         ; tvar%found = .true.
            CASE ( 'par' )         ; tvar%found = .true.
            CASE ( 'nir' )         ; tvar%found = .true.
            CASE ( 'uva' )         ; tvar%found = .true.
            CASE ( 'uvb' )         ; tvar%found = .true.
            CASE ( 'sed_zone' )    ; tvar%found = .true.
            CASE ( 'wind_speed' )  ; tvar%found = .true.
            CASE ( 'par_sf' )      ; tvar%found = .true.
            CASE ( 'taub' )        ; tvar%found = .true.
            CASE ( 'air_temp' )    ; tvar%found = .true.
            CASE ( 'nearest_active' ) ; tvar%found = .true. ; request_nearest = .true.
            CASE ( 'nearest_depth' )  ; tvar%found = .true. ; request_nearest = .true.
         !  CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
         ELSE
            d = d + 1
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
         ELSE
            v = v + 1
         ENDIF
      ENDIF
      IF ( .NOT. tvar%found ) THEN
         print *, "ERROR: Undefined variable ", trim(tvar%name)
         err_count = err_count + 1
      ENDIF
   ENDDO

   if ( n_vars < v ) print *,"More vars than expected"
   if ( n_vars_ben < sv ) print *,"More sheet vars than expected"
   if ( n_vars_diag < sd + d ) print *,"More diag vars than expected"
   if ( n_vars_diag_sheet < sd ) print *,"More sheet diag vars than expected"

   IF ( err_count > 0 ) CALL STOPIT("*** ERRORs in configuration")
END SUBROUTINE check_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_column(column, col, cc, cc_diag, flux_pel, flux_atm, flux_ben, flux_rip)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in) :: col
   AED_REAL, TARGET, INTENT(in) :: cc(:,:)       !# (n_vars, n_layers)
   AED_REAL, TARGET, INTENT(in) :: cc_diag(:,:)  !# (n_vars, n_layers)
   AED_REAL, TARGET, INTENT(inout) :: flux_pel(:,:) !# (n_vars, n_layers)
   AED_REAL, TARGET, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_ben(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_rip(:)   !# (n_vars)
!
!LOCALS
   INTEGER :: av, i, top, bot
   INTEGER :: v, d, sv, sd, ev
   TYPE(aed2_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   top = surf_map(col)
   bot = benth_map(col)

   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed2_vars

      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "ERROR getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => temp(top:bot)
            CASE ( 'salinity' )    ; column(av)%cell => salt(top:bot)
            CASE ( 'density' )     ; column(av)%cell => rho(top:bot)
            CASE ( 'layer_ht' )    ; column(av)%cell => h(top:bot)
            CASE ( 'layer_area' )  ; column(av)%cell_sheet => area(col)
            CASE ( 'rain' )        ; column(av)%cell_sheet => rain(col)
            CASE ( 'rainloss' )    ; column(av)%cell_sheet => rainloss(col)
            CASE ( 'material' )    ; column(av)%cell_sheet => zone(col)
            CASE ( 'bathy' )       ; column(av)%cell_sheet => bathy(col)
            CASE ( 'extc_coef' )   ; column(av)%cell => extcoeff(top:bot)
            CASE ( 'tss' )         ; column(av)%cell => tss(top:bot)
            CASE ( 'par' )         ; IF (link_ext_par) THEN
                                        column(av)%cell => lpar(top:bot)
                                     ELSE
                                        column(av)%cell => par(top:bot)
                                     ENDIF
            CASE ( 'nir' )         ; column(av)%cell => nir(top:bot)
            CASE ( 'uva' )         ; column(av)%cell => uva(top:bot)
            CASE ( 'uvb' )         ; column(av)%cell => uvb(top:bot)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => zone(col)
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd(col)
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0(col)
            CASE ( 'taub' )        ; column(av)%cell_sheet => col_taub
            CASE ( 'air_temp' )    ; column(av)%cell_sheet => air_temp(col)

            CASE ( 'nearest_active' ) ; column(av)%cell_sheet => nearest_active(col);
            CASE ( 'nearest_depth' )  ; column(av)%cell_sheet => nearest_depth(col);
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         d = d + 1
         IF ( tvar%sheet ) THEN
            column(av)%cell_sheet => cc_diag(d, bot)
         ELSE
            column(av)%cell => cc_diag(d,top:bot)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => cc(n_vars+sv, bot)
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => cc(n_vars+sv, top)
            ENDIF
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
            column(av)%flux_rip => flux_rip(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => cc(v,top:bot)
            column(av)%flux_pel => flux_pel(v,top:bot)
            column(av)%flux_ben => flux_ben(v)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_rip => flux_rip(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE calculate_fluxes(column, count, flux_pel, flux_atm, flux_ben, flux_rip, h)
!-------------------------------------------------------------------------------
! Checks the current values of all state variables and repairs these
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in) :: count
   AED_REAL, INTENT(inout) :: flux_pel(:,:) !# (n_vars, n_layers)
   AED_REAL, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, INTENT(inout) :: flux_ben(:)   !# (n_vars)
   AED_REAL, INTENT(inout) :: flux_rip(:)   !# (n_vars)
   AED_REAL, INTENT(inout) :: h(:)          !# (n_layers)
!
!LOCALS
   INTEGER :: i
!-------------------------------------------------------------------------------
!BEGIN
   flux_pel = zero_
   flux_atm = zero_
   flux_ben = zero_
   flux_rip = zero_

   !# Calculate temporal derivatives due to air-water exchange.
   CALL aed2_calculate_surface(column, 1)

   !# Distribute the fluxes into pelagic surface layer
   IF ( do_2d_atm_flux .OR. count > 1 ) &
      flux_pel(:,1) = flux_pel(:,1) + flux_atm(:)/h(1)

   IF ( .NOT. do_zone_averaging ) THEN
      !# Calculate temporal derivatives due to benthic exchange processes.
      CALL aed2_calculate_benthic(column, count)

      !# Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
      flux_pel(:,count) = flux_pel(:,count)/h(count)
   ENDIF

   !# Add pelagic sink and source terms for all depth levels.
   DO i=1,count
      CALL aed2_calculate(column, i)
   ENDDO
END SUBROUTINE calculate_fluxes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_states(column, top, bot)
!-------------------------------------------------------------------------------
!USES
   USE IEEE_ARITHMETIC
!
!ARGUMENTS
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: top, bot
!
!LOCALS
   TYPE(aed2_variable_t),POINTER :: tv
   INTEGER i,v,d,lev
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lev=top, bot
      !CALL aed2_equilibrate(column, lev)
      v = 0; d = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               IF ( do_limiter ) THEN
                  IF ( .NOT. isnan(min_(v)) ) THEN
                     IF ( cc(v, lev) < min_(v) ) cc(v, lev) = min_(v)
                  ELSE IF (.NOT. no_glob_lim) THEN
                     IF ( cc(v, lev) < glob_min ) THEN
                        cc(v, lev) = MISVAL
                        print*, "Variable ", v, " below glob_min"
                     ENDIF
                  ENDIF
                  IF ( .NOT. isnan(max_(v)) ) THEN
                     IF ( cc(v, lev) > max_(v) ) cc(v, lev) = max_(v)
                  ELSE IF (.NOT. no_glob_lim) THEN
                     IF ( cc(v, lev) > glob_max ) THEN
                        cc(v, lev) = MISVAL
                        print*, "Variable ", v, " above glob_max"
                     ENDIF
                  ENDIF
               ENDIF
            ELSE IF ( tv%diag .AND. .NOT. no_glob_lim ) THEN
               d = d + 1
               IF ( cc_diag(d, lev) < glob_min .OR. cc_diag(d, lev) > glob_max ) THEN
                  print *, "Diagnostic ", d, " exceeded global bounds", cc_diag(d, lev)
                  cc_diag(d, lev) = MISVAL
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO
END SUBROUTINE check_states
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE fill_nearest(nCols)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in) :: nCols
!
!LOCALS
   INTEGER  :: k, col, prev, next
!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( ALLOCATED(route_table) ) THEN
      DO col=1, nCols
         IF (active(col) .AND. h(benth_map(col))>=min_water_depth) THEN
            nearest_active(col) = col
            nearest_depth(col) = h(benth_map(col)) + bathy(col)
         ELSE
            k = route_table(col)
            DO WHILE ( .NOT. active(k) .OR. h(benth_map(k))<min_water_depth)
               IF ( k == route_table(k) ) EXIT
               k = route_table(k)
            ENDDO
            nearest_active(col) = k
            nearest_depth(col) = h(benth_map(k)) + bathy(k)
            ! this needs fixing to sum over top:bot, as h is layer thicknesses, not references to datum
         ENDIF
      ENDDO
   ELSE
      nearest_active = 0.
   ENDIF
END SUBROUTINE fill_nearest
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE do_aed2_models(nCells, nCols)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in) :: nCells, nCols
!
!LOCALS
   TYPE(aed2_variable_t),POINTER :: tv
   AED_REAL :: tr

   AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben),       &
               flux_rip(n_vars+n_vars_ben)
   TYPE (aed2_column_t) :: column(n_aed2_vars)

   INTEGER  :: i, j, col, lev, top, bot, v, na, d
   AED_REAL :: rain_loss
   LOGICAL  :: aed_active_col
   AED_REAL,DIMENSION(:),POINTER :: tpar
!
!-------------------------------------------------------------------------------
!BEGIN
! print *," START do_aed2_models"

   !#--------------------------------------------------------------------
   !# START-UP JOBS
   rainloss = zero_

   IF ( request_nearest ) CALL fill_nearest(nCols)

   IF ( .not. reinited ) CALL re_initialize()

   ThisStep = ThisStep + 1

   IF ( do_zone_averaging ) THEN
      IF (link_ext_par) THEN
         tpar => lpar
      ELSE
         tpar => par
      ENDIF
      CALL calc_zone_areas(nCols, active, temp, salt, h, area, wnd, rho,       &
                 extcoeff, I_0, tpar, tss, rain, rainloss, air_temp, bathy, col_taub)
   ENDIF

!!$OMP DO
   !#--------------------------------------------------------------------
   !# LOOP THROUGH COLUMNS DOING JOBS PRIOR TO THE KINETICS BEING SOLVED
   DO col=1, nCols
      ! move to next column if dry
      IF (.NOT. active(col)) CYCLE

      ! identify cell indicies within the domain
      top = surf_map(col)
      bot = benth_map(col)

      ! set column data
      CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben, flux_rip)

      ! compute vertical settling/mobility
      v = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               ! only for state_vars that are not sheet
               IF ( .NOT. isnan(tv%mobility) ) THEN
                  ! default to ws that was set during initialisation
                  ws(top:bot,i) = tv%mobility
               ELSE
                  ! zero nan values
                  ws(top:bot,i) = zero_
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      DO lev = top, bot
        ! update ws for modules that use the mobility method
        CALL aed2_mobility(column, lev-top+1, ws(lev,:))
      ENDDO
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) .AND. SUM(ABS(ws(top:bot,i)))>zero_ ) THEN
               CALL Settling(bot-top+1, dt, h(top:bot), ws(top:bot,i), Fsed_setl(col), column(i)%cell)
            ENDIF
         ENDIF
      ENDDO
      CALL check_states(column, top, bot)

      !# populate local light/extc arrays one column at a time
      IF (.NOT. link_ext_par) &  !#MH check link_ext_par logic
         CALL Light(column, bot-top+1, I_0(col), extcoeff(top:bot), par(top:bot), h(top:bot))

      !# non-PAR light fudge
      nir(top:bot) = (par(top:bot)/0.45) * 0.510
      uva(top:bot) = (par(top:bot)/0.45) * 0.035
      uvb(top:bot) = (par(top:bot)/0.45) * 0.005
   ENDDO
!!$OMP END DO

   IF ( do_zone_averaging ) THEN
      CALL copy_to_zone(nCols, cc, area, active, benth_map)
      CALL compute_zone_benthic_fluxes(n_aed2_vars, dt)
   ENDIF

!!$OMP DO
   !#--------------------------------------------------------------------
   !# THIS IS THE MAIN WQ SOLUTION LOOP
   DO col=1, nCols

      !# find top and bottom cell indicies based on maps provided by the host
      top = surf_map(col)
      bot = benth_map(col)

      !# compute bottom shear stress for this column based on ustar from host
      col_taub = rho(bot)*(ustar_bed(col)*ustar_bed(col))
      IF ( link_wave_stress .AND. ASSOCIATED(wv_uorb) .AND. ASSOCIATED(wv_t) ) THEN
         CALL Stress(h(bot),rho(bot),col_taub,ustar_bed(col),wv_uorb(col),wv_t(col))
      ELSE
         CALL Stress(h(bot),rho(bot),col_taub,ustar_bed(col))
      ENDIF

      !# set column data structure from global arrays
      CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben, flux_rip)

      !# See if there has been a change in active(wet/dry) state
!     IF (active(col) .NE. pactive(col)) THEN
!        IF (active(col)) THEN
!        ELSE
!        ENDIF
!     ENDIF

      !# do riparian interfaces for this column and update fluxes
      flux_rip = zero_
      shadefrac(col) = one_  ! zero_
      rainloss(col) = zero_
      aed_active_col = active(col)
      IF( h(benth_map(col))<min_water_depth ) aed_active_col = .false.  ! MH TUFLOWFV 4cm dry cells
      IF ( .NOT.  Riparian(column, aed_active_col, shadefrac(col), rainloss(col)) ) THEN
         IF ( request_nearest ) THEN
            na = nearest_active(col)
            ! Check for cells that are routed to dry pools
            IF ( h(benth_map(na)) >= min_water_depth ) THEN
               cc(:,benth_map(na))=cc(:,benth_map(na))+dt*flux_rip(:) &
                              * MIN((area(col)/area(na)),1e2)/h(benth_map(na))
            ENDIF
         ENDIF
         CYCLE
      ENDIF

      !# do non-kinetic updates to BGC variables (eq equilibration)
      IF ( ThisStep >= n_equil_substep ) CALL Update(column, bot-top+1)

      !# find the particles in this column and update particle bgc
      IF (do_particle_bgc) CALL Particles(column, bot-top+1, h(top:bot))

#if _NO_ODE_
      !# for this column, do the main kinetic/bgc flux calculation
      !# (this includes water column, surface and benthic interfaces)
      CALL calculate_fluxes(column, bot-top+1, flux(:,top:bot), flux_atm, flux_ben, flux_rip, h(top:bot))

      !# now go forth and solve the ODE (Euler! - link to fv_ode would be nice)
      DO lev = top, bot
         DO i = 1, n_vars
            cc(i,lev)=cc(i,lev)+dt*flux(i,lev)
#if DEBUG>1
            !# check for NaNs
            IF ( isnan(cc(i,lev)) ) THEN
               print*,'Nan at i = ', i, ' lev = ', lev
               print*,'h(lev) = ', h(lev), ' flux(i,lev) = ', flux(i,lev)
               print*,'Top of column @ ', top, ' bottom of column @ ', bot
               call STOPIT('NaN value')
            ENDIF
#endif
         ENDDO ! vars
      ENDDO  ! levels

      !# benthic state variables
      DO i = n_vars+1, n_vars+n_vars_ben
        cc(i,bot)=cc(i,bot)+dt*flux_ben(i)
      ENDDO

      !# add riparian flux
      IF ( do_zone_averaging ) THEN ! Untested
         DO i = n_vars+1, n_vars+n_vars_ben
            cc(i,bot)=cc(i,bot)+dt*flux(i,bot)
#if DEBUG>1
            !# check for NaNs
            IF ( isnan(cc(i,bot)) ) THEN
               print*,'Nan at i = ', i, ' bot = ', bot
               print*,'h(bot) = ', h(bot), ' flux(i,bot) = ', flux(i,bot)
               call STOPIT('NaN value')
            ENDIF
#endif
         ENDDO ! ben vars
      ENDIF
#endif


      !# now the bgc updates are complete, update links to host model
      CALL BioDrag(column, bot-top+1, bio_drag(col))
      CALL BioExtinction(column, bot-top+1, extcoeff(top:bot))
      !CALL BioDensity()

      CALL check_states(column, top, bot)
!     IF (active(col) .NE. pactive(col)) THEN
!        IF (active(col)) THEN
!        ELSE
!        ENDIF
!        pactive(col) = active(col)
!     ENDIF
   ENDDO ! cols
!!$OMP END DO

   IF ( do_zone_averaging ) &
      CALL copy_from_zone(nCols, cc, area, active, benth_map)

   IF ( ThisStep >= n_equil_substep ) ThisStep = 0

   IF ( display_minmax ) THEN
      v = 0; d = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               WRITE(*,'(1X,"VarLims: ",I0,1X,"<=> ",f15.8,f15.8," : ",A)')v,MINVAL(cc(v,:)),MAXVAL(cc(v,:)),TRIM(tv%name)
               !print *,'VarLims',v,TRIM(tv%name),MINVAL(cc(v,:)),MAXVAL(cc(v,:))
            ELSE IF ( tv%diag .AND. .NOT. no_glob_lim ) THEN
               d = d + 1
               WRITE(*,'(1X,"DiagLim: ",I0,1X,"<=> ",f15.8,f15.8," : ",A)')d,MINVAL(cc_diag(d,:)),MAXVAL(cc_diag(d,:)),TRIM(tv%name)
               !print *,'DiagLim',d,TRIM(tv%name),MINVAL(cc_diag(d,:)),MAXVAL(cc_diag(d,:))
            ENDIF
         ENDIF
      ENDDO
    ENDIF

CONTAINS

   !###############################################################################
   SUBROUTINE re_initialize()
   !-------------------------------------------------------------------------------
   !ARGUMENTS
   !
   !LOCALS
      INTEGER  :: i, col, lev, top, bot, count, nCols
   !
   !-------------------------------------------------------------------------------
   !BEGIN
      nCols = ubound(active, 1)
      DO col=1, nCols
         top = surf_map(col)
         bot = benth_map(col)
         count = top-bot+1
         CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben, flux_rip)
         DO lev=1, count
            CALL aed2_initialize(column, lev)
         ENDDO
      ENDDO
      reinited = .TRUE.
   END SUBROUTINE re_initialize

END SUBROUTINE do_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE clean_aed2_models
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Deallocate internal arrays
   IF (allocated(ws))             deallocate(ws)
   IF (allocated(total))          deallocate(total)
   IF (allocated(nir))            deallocate(nir)
   IF (allocated(par))            deallocate(par)
   IF (allocated(uva))            deallocate(uva)
   IF (allocated(uvb))            deallocate(uvb)
!  IF (allocated(pactive))        deallocate(pactive)
END SUBROUTINE clean_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Light(column, count, Io, extc, par_, h_)
!-------------------------------------------------------------------------------
!
! Calculate photosynthetically active radiation over entire column
! based on surface radiation, and background and biotic extinction.
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(in)    :: Io
   AED_REAL, INTENT(inout) :: extc(:)
   AED_REAL, INTENT(inout) :: par_(:)
   AED_REAL, INTENT(inout) :: h_(:)
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: zz, localext, localshade
!
!-------------------------------------------------------------------------------
!BEGIN
   zz = zero_
   localext = zero_

   CALL BioExtinction(column,count,extc)

   localext = extc(1)
   zz = 0.001 !0.5*h_(1)    !MH: assume top of layer
   par_(1) = 0.45 * Io * EXP( -(localext) * zz )

   IF (count <= 1) RETURN

   DO i = 2, count
      localext = extc(i)

      !zz = zz + 0.5*h_(i)
      zz = h_(i)
      par_(i) = par_(i-1) * EXP( -(localext) * zz )
   ENDDO
END SUBROUTINE Light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Settling(N,dt,h,wvel,Fsed,Y)
!-------------------------------------------------------------------------------
!
! Update settling of AED2 state variables in a given column
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)     :: N       !# number of vertical layers
   AED_REAL,INTENT(in)    :: dt      !# time step (s)
   AED_REAL,INTENT(in)    :: h(:)    !# layer thickness (m)
   AED_REAL,INTENT(in)    :: wvel(:) !# vertical advection speed
   AED_REAL,INTENT(inout) :: Fsed    !# value of sediment input due to settling
   AED_REAL,INTENT(inout) :: Y(:)
!
!CONSTANTS
   INTEGER,PARAMETER :: itmax=100
!
!LOCALS
   INTEGER  :: i,k,it
   AED_REAL :: step_dt
   AED_REAL :: Yc
   AED_REAL :: c,cmax
   AED_REAL :: cu(0:N)
!
!-------------------------------------------------------------------------------
!BEGIN
   Fsed = 0. !# initialize sediment settling fluxes with zero
   cu   = 0. !# initialize interface fluxes with zero
   cmax = 0. !# initialize maximum Courant number

   !# compute maximum Courant number
   !      calculated as number of layers that the particles will travel based
   !      on settling or buoyancy velocity.
   !      This number is then used to split the vertical movement
   !      calculations to limit movement across a single layer
   DO k=1,N-1
      !# sinking particles
      c=abs(wvel(k+1))*dt/(0.5*(h(k+1)+h(k)))
      IF (c > cmax) cmax=c
      !# rising particles
      c=abs(wvel(k))*dt/(0.5*(h(k+1)+h(k)))
      IF (c > cmax) cmax=c
   ENDDO

   it=min(itmax,int(cmax)+1)
   step_dt = dt / float(it);

   !# splitting loop
   DO i = 1,it
      !# vertical loop
      DO k=1,N-1
         !# compute the slope ratio
         IF (wvel(k) > 0.) THEN !# Particle is rising
            Yc=Y(k)       !# central value
         ELSE !# negative speed Particle is sinking
            Yc=Y(k+1)     !# central value
         ENDIF

         !# compute the limited flux
         cu(k)=wvel(k) * Yc
      ENDDO

      !# do the upper boundary conditions
      cu(N) = zero_       !# limit flux into the domain from atmosphere

      !# do the lower boundary conditions
      IF (wvel(1) > 0.) THEN !# Particle is rising
         cu(0) = 0.  !flux from benthos is zero
      ELSE  !# Particle is settling
         cu(0) = wvel(1)*Y(1)
         Fsed = cu(0) * step_dt !# flux settled into the sediments per sub time step
      ENDIF
      !# do the vertical advection step including positive migration
      !# and settling of suspended matter.
      DO k=1,N
          Y(k)=Y(k) - step_dt * ((cu(k) - cu(k-1)) / h(k))
      ENDDO
   ENDDO !# end of the iteration loop
   Fsed = Fsed / dt !# Average flux rate for full time step used in AED2
END SUBROUTINE Settling
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
LOGICAL FUNCTION Riparian(column, actv, shade_frac, rain_loss)
!-------------------------------------------------------------------------------
!
! Do riparian functionality, including operations in dry and fringing cells &
! populate feedback arrays to the host model associated with riparian effects
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   LOGICAL,  INTENT(in)    :: actv
   AED_REAL, INTENT(inout) :: shade_frac, rain_loss
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: localshade
   AED_REAL :: localrainl
!
!-------------------------------------------------------------------------------
!BEGIN
   !# compute the methods relevant to either DRY or WET cells
   IF (.NOT. actv ) THEN
      CALL aed2_calculate_dry(column, 1);
      CALL aed2_calculate_riparian(column, 1, zero_);
   ELSE
      CALL aed2_calculate_riparian(column, 1, one_);
   ENDIF

   !# update feedback arrays to host model, to reduce rain (or if -ve then add flow)
   CALL aed2_rain_loss(column, 1, localrainl);
   IF (link_rain_loss) rain_loss = localrainl


   !# update feedback arrays to shade the water (ie reduce incoming light, Io)
   CALL aed2_light_shading(column, 1, localshade)
   IF (link_solar_shade) shade_frac = localshade

   Riparian = actv

END FUNCTION Riparian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Update(column,count)
!-------------------------------------------------------------------------------
!
! Do non-kinetic (eg equilibrium) updates to state variables in AED2 modules
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in) :: count
!
!LOCAL VARIABLES:
   INTEGER :: lev
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lev=1,count
      CALL aed2_equilibrate(column, lev)
   ENDDO
END SUBROUTINE Update
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Particles(column, count, h_)
!-------------------------------------------------------------------------------
!
! Calculate biogeochemical transformations on particles !TO BE COMPELTED!
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(inout) :: h_(:)
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: zz
!
!-------------------------------------------------------------------------------
!BEGIN
   zz = zero_

 ! CALL aed2_particle_bgc(column,1,ppid ...)
   IF (count <= 1) RETURN

   DO i = 2, count
  !   CALL aed2_particle_bgc(column,count,ppid ...)
   ENDDO
END SUBROUTINE Particles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE BioExtinction(column,count,extc)
!-------------------------------------------------------------------------------
!
! Calculate the specific light attenuation additions due to AED2 modules
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(inout) :: extc(:)
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: localext
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_

   CALL aed2_light_extinction(column, 1, localext)
   IF (link_water_clarity) THEN
     extc(1) = localext
   ELSE
     extc(1) = localext + Kw
   END IF

   IF (count <= 1) RETURN

   DO i = 2, count
      CALL aed2_light_extinction(column, i, localext)
     IF (link_water_clarity) THEN
       extc(i) = localext
     ELSE
       extc(i) = localext + Kw
     END IF
   ENDDO
END SUBROUTINE BioExtinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE BioDrag(column,count,bdrag)
!-------------------------------------------------------------------------------
!
! Calculate the drag addition to be returned to the host model due to vegetation
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(inout) :: bdrag
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: localdrag
!
!-------------------------------------------------------------------------------
!BEGIN
   bdrag = zero_
   localdrag = zero_

   CALL aed2_bio_drag(column, count, localdrag)

   IF (link_bottom_drag) bdrag = localdrag
END SUBROUTINE BioDrag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE BioDensity(column,count,bio_density)
!-------------------------------------------------------------------------------
!
! Calculate the density addition to be returned to the host model due to WQ
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(inout) :: bio_density(:)
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: zz, localdensity
!
!-------------------------------------------------------------------------------
!BEGIN
   RETURN
END SUBROUTINE BioDensity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE Stress(h,rho,taub,ustar,uorb,wvperiod)
!-------------------------------------------------------------------------------
!
! Calculate the density addition to be returned to the host model due to WQ
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL, INTENT(in) :: h,rho,ustar
   AED_REAL, INTENT(in),OPTIONAL :: uorb,wvperiod
   AED_REAL, INTENT(inout) :: taub
!
!LOCAL VARIABLES:
   AED_REAL,PARAMETER :: pi = 4.0*ATAN(1.0)
   AED_REAL,PARAMETER :: nuw = 1.05e-6
   AED_REAL,PARAMETER :: ksw = 0.001
   AED_REAL,PARAMETER :: kappa = 0.41
   AED_REAL :: Aw,Rew,fwr,fws,fw,tauw
!
!-------------------------------------------------------------------------------
!BEGIN


   ! Current shear stress
   taub = rho*(ustar**2)

   IF( .NOT.PRESENT(uorb) .OR. .NOT.PRESENT(wvperiod) ) RETURN

   ! Shear stress due to wave-induced orbital velocity
   IF (h<0.05 .OR. uorb<0.001 .OR. wvperiod<0.01) THEN
   ELSE
      Aw = uorb*wvperiod/(2.*pi)
      Rew = uorb*Aw/nuw
      ! Smooth friction factor
      fws = 0.035*Rew**(-0.16)
      ! Turbulent friction factor
      fwr = EXP(5.21*(ksw/Aw)**0.194-5.98)
      fw = MAX(fws,fwr)
      ! Calculate wave stress (wave_factor allows for user scaling)
      tauw = (0.5*rho*fw*uorb**2)*wave_factor

      ! Total current + wave stress
      !-- a) mean bed shear
      taub = taub*( 1.+1.2*(tauw/(taub+tauw+1e-10))**3.2 )
      !-- b) RMS bed shear
      taub = SQRT( taub**2 + 0.5*tauw**2 )
   END IF
END SUBROUTINE Stress
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE fv_aed2

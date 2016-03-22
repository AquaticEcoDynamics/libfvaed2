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
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created Aug 2013                                                            #
!# Update Mar 2016: Add new env variables and feedback links                   #
!#                                                                             #
!###############################################################################

#include "aed2.h"

#define FV_AED_VERS "0.9.21"

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

   !# Name of file being used to load initial values for benthic or benthic_diag vars
   CHARACTER(len=128) :: init_values_file = ''

   !# Maps of surface, bottom and wet/dry (active) cells
   LOGICAL,DIMENSION(:),POINTER :: active
   INTEGER,DIMENSION(:),POINTER :: surf_map, benth_map

   !# Arrays for work, vertical movement, and cross-boundary fluxes
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: flux
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: ws, total
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: min_, max_

   !# Arrays for environmental variables (used if they are not supplied externally)
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: nir
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: par
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: uva
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: uvb
   AED_REAL,DIMENSION(:),POINTER :: lpar

   AED_REAL,TARGET :: col_taub  ! a temp var for the taub for column (computed from ustar_bed)

   !# External variables
   AED_REAL :: dt
   AED_REAL,DIMENSION(:),  POINTER :: temp, salt, rho, nuh, h, z
   AED_REAL,DIMENSION(:),  POINTER :: extcoeff, tss, bio_drag
   AED_REAL,DIMENSION(:),  POINTER :: I_0, wnd, air_temp, rain
   AED_REAL,DIMENSION(:),  POINTER :: area, bathy, shadefrac, rainloss
   AED_REAL,DIMENSION(:),  POINTER :: Fsed_setl, ustar_bed
   AED_REAL,DIMENSION(:,:),POINTER :: rad
   INTEGER, DIMENSION(:,:),POINTER :: mat

   !# Switches for configuring model operation and active links with the host model
   LOGICAL :: link_water_clarity, link_water_density,                          &
              link_bottom_drag, link_surface_drag,                             &
              link_ext_par, ext_tss_extinction,                                &
              link_solar_shade, link_rain_loss,                                &
              do_limiter, do_2d_atm_flux, do_particle_bgc

   LOGICAL :: old_zones = .TRUE.
   LOGICAL :: do_zone_averaging = .FALSE.

   !# Integers storing number of variables being simulated
   INTEGER :: n_aed2_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet

CONTAINS
!===============================================================================


!###############################################################################
!SUBROUTINE STOPIT(message)
!!-------------------------------------------------------------------------------
!!ARGUMENTS
!   CHARACTER(*) :: message
!!-------------------------------------------------------------------------------
!   PRINT *,message
!   STOP
!END SUBROUTINE STOPIT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE init_aed2_models(namlst,dname,nwq_var,nben_var,ndiag_var,names,bennames,diagnames)
!-------------------------------------------------------------------------------
! This routine is called by tuflowfv to define numbers of variables. TuflowFV
! will allocate the variables after return from this routine.
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
   CHARACTER(len=128) :: aed2_nml_file
   CHARACTER(len=128) :: tname
   AED_REAL           :: base_par_extinction, tss_par_extinction
   INTEGER            :: status
   INTEGER            :: n_sd, i, j
   TYPE(aed2_variable_t),POINTER :: tvar

   CHARACTER(len=64) :: models(64)
   NAMELIST /aed2_models/ models

   NAMELIST /aed2_bio/ solution_method, link_bottom_drag,                 &
                      link_surface_drag, link_water_density,              &
                      link_water_clarity, aed2_nml_file,                  &
                      link_ext_par, base_par_extinction,                  &
                      ext_tss_extinction, tss_par_extinction,             &
                      do_limiter, do_2d_atm_flux, do_zone_averaging,      &
                      link_solar_shade, link_rain_loss, init_values_file, &
                      do_particle_bgc
!
!-------------------------------------------------------------------------------
!BEGIN
   print *, "*** Using fv_aed2 version ", TRIM(FV_AED_VERS)

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
   do_particle_bgc = .false.

   IF ( aed2_init_core(dname) /= 0 ) STOP "Initialisation of aed2_core failed"
   tname = TRIM(dname)//'aed2.nml'
   print *,"Reading aed2_models config from ",TRIM(tname)
   OPEN(namlst,file=tname,action='read',status='old',iostat=status)
   IF ( status /= 0 ) CALL STOPIT("Cannot open file " // TRIM(tname))
   READ(namlst,nml=aed2_bio,iostat=status)
   IF ( status /= 0 ) STOP "Cannot read namelist entry aed2_bio"
   Kw = base_par_extinction
   Ksed = tss_par_extinction
   

   models = ''
   READ(namlst, nml=aed2_models, iostat=status)
   IF ( status /= 0 ) STOP "Cannot read namelist entry aed2_models"

   IF ( do_zone_averaging ) old_zones = .FALSE.

   DO i=1,size(models)
      IF (models(i)=='') EXIT
      CALL aed2_define_model(models(i), namlst)
   ENDDO

   n_aed2_vars = aed2_core_status(nwq_var, nben_var, ndiag_var, n_sd)

#if DEBUG
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         print *,"AED2 var ", i, tvar%sheet, tvar%diag, tvar%extern, TRIM(tvar%name)
      ELSE
         print *,"AED2 var ", i, " is empty"
      ENDIF
   ENDDO
#endif

   ndiag_var = ndiag_var + n_sd
   n_vars = nwq_var
   n_vars_ben = nben_var
   n_vars_diag = ndiag_var
   n_vars_diag_sheet = n_sd

#if DEBUG
   print*,'AED2 init_aed2_models : n_aed2_vars = ',n_aed2_vars,' nwq_var = ',nwq_var,' nben_var ',nben_var
#endif

   CALL check_data

   !# names = grab the names from info
   ALLOCATE(names(1:nwq_var),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (names)'
   ALLOCATE(bennames(1:nben_var),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (bennames)'
   IF ( .NOT. ALLOCATED(diagnames) ) ALLOCATE(diagnames(ndiag_var))
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (diagnames)'

   ALLOCATE(min_(1:nwq_var+nben_var)) ; ALLOCATE(max_(1:nwq_var+nben_var))

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( .NOT. (tvar%sheet .OR. tvar%diag .OR. tvar%extern) ) THEN
            j = j + 1
            names(j) = TRIM(tvar%name)
            min_(j) = tvar%minimum
            max_(j) = tvar%maximum
            print *,"AED2 var name(",j,") : ", TRIM(names(j))
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
            print *,"AED2 var_ben name(",j,") : ", TRIM(bennames(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( tvar%diag ) THEN
            j = j + 1
            diagnames(j) = TRIM(tvar%name)
            print *,"AED2 diag name(",j,") : ", trim(diagnames(j))
         ENDIF
      ENDIF
   ENDDO

   CLOSE(namlst)
END SUBROUTINE init_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE init_var_aed2_models(nCells, cc_, cc_diag_, nwq, nwqben, sm, bm)
!-------------------------------------------------------------------------------
! At this point TuflowFV should have allocated the variable space.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: nCells
   AED_REAL,POINTER,DIMENSION(:,:),INTENT(in) :: cc_, cc_diag_
   INTEGER,INTENT(inout)                   :: nwq, nwqben
   INTEGER,POINTER,DIMENSION(:),INTENT(in) :: sm, bm
!
!LOCALS
   INTEGER :: rc, av, v, sv, d, sd
   TYPE(aed2_variable_t),POINTER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   nwq = n_vars
   nwqben = n_vars_ben

   print *,'init_var_aed2_models : nwq = ',nwq,' nwqben = ',nwqben

   cc => cc_
   cc_diag => cc_diag_
   surf_map => sm
   benth_map => bm
   
   ! Allocate state variable array
   IF ( .NOT. ASSOCIATED(cc) ) STOP ' Error : no association for (cc)'
   cc = 0.

   IF (.not. ASSOCIATED(cc_diag) ) STOP ' Error : no association for (cc_diag)'
   cc_diag = 0.

   ! Allocate array with vertical movement rates (m/s, positive for upwards),
   ! and set these to the values provided by the model.
   ALLOCATE(ws(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   ws = 0.

   !!# place holder for lagranigan particles
   !IF(do_particle_bgc) THEN
   !  pp => pp_
   !END IF

   ! Allocate array for photosynthetically active radiation (PAR).
   ! This will be calculated internally during each time step.
   ALLOCATE(par(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (par)'
   par = 0.
   ALLOCATE(nir(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (nir)'
   nir = 0.
   ALLOCATE(uva(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (uva)'
   uva = 0.
   ALLOCATE(uvb(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (uvb)'
   uvb = 0.

   !# Allocate array for sedimentation fluxes and initialize these to zero (no flux).
   ALLOCATE(Fsed_setl(1:nCells),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (Fsed_setl)'
   Fsed_setl = 0.

   !# Now set initial values
   v = 0 ; sv = 0;
   DO av=1,n_aed2_vars
      IF ( .NOT.  aed2_get_var(av, tv) ) STOP "Error getting variable info"
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

!   CALL fv_initialize()

   ALLOCATE(flux(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux)'
#if !_NO_ODE_
   ALLOCATE(flux2(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux2)'
   ALLOCATE(flux3(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux3)'
   ALLOCATE(flux4(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux4)'
   ALLOCATE(cc1(n_vars, nCells),stat=rc)   ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc1)'
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
   USE fv_aed2_csv_reader
   !
   !LOCALS
      INTEGER :: unit, nccols, ccol
      CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
      AED_REAL,DIMENSION(:),ALLOCATABLE :: values
      INTEGER :: idx_col = 0, numv = 0, numd = 0, t
      INTEGER,DIMENSION(:),ALLOCATABLE :: vars, vmap
      INTEGER,DIMENSION(:),ALLOCATABLE :: dvar, dmap
      LOGICAL,DIMENSION(:),ALLOCATABLE :: vsheet, dsheet
      LOGICAL :: meh
   !
   !BEGIN
   !----------------------------------------------------------------------------
      unit = aed_csv_read_header(init_values_file, csvnames, nccols)
      print *,'benthic vars initialised from file : ', csvnames
      IF (unit <= 0) RETURN !# No file found
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

      IF ( idx_col > 0 ) THEN
         v = 0 ; sv = 0; d = 0; sd = 0
         DO av=1,n_aed2_vars
            IF ( .NOT. aed2_get_var(av, tv) ) STOP "Error getting variable info"
            IF ( .NOT. ( tv%extern ) ) THEN  !#  dont do environment vars
               DO ccol=1,nccols
                  IF ( same_str_icase(tv%name, csvnames(ccol)) ) THEN
                     IF (tv%diag) THEN
                        numd = numd + 1
                        dmap(numd) = ccol
                        dvar(numd) = numd
                        dsheet(numd) = tv%sheet
                     ELSE
                        numv = numv + 1
                        vmap(numv) = ccol
                        IF ( tv%sheet ) THEN
                           sv = sv + 1
                           vars(numv) = n_vars + sv
                        ELSE
                           v = v + 1
                           vars(numv) = v
                        ENDIF
                        vsheet(numd) = tv%sheet
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         DO WHILE ( aed_csv_read_row(unit, values) )
            t = INT(values(idx_col))
 !             print *,'values',values
            DO v=1,numv
               If ( vsheet(v) ) THEN
                  cc(vars(v), bm(t)) = values(vmap(v))
               ELSE
                  cc(vars(v), sm(t):bm(t)) = values(vmap(v))
               ENDIF
            ENDDO
            DO v=1,numd
               If ( vsheet(v) ) THEN
 !             print *,'dmpa',dmap(v)
                  cc_diag(dvar(v), bm(t)) = values(dmap(v))
               ELSE
                  cc_diag(dvar(v), sm(t):bm(t)) = values(dmap(v))
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



END SUBROUTINE init_var_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE set_env_aed2_models(dt_,              &   !
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
! Provide information about tuflowfv data not declared by aed2_bio
!-------------------------------------------------------------------------------
!ARGUMENTS
   DOUBLETYPE, INTENT(in) :: dt_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: temp_, salt_, rho_, h_,    &
                                                    area_, tss_, extcoeff_, z_
   AED_REAL, INTENT(in), DIMENSION(:,:), POINTER :: rad_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: I_0_, wnd_, ustar_bed_
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: rain_, bathy_ ! JC    bathy,
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: biodrag_, ustar_surf_, solarshade_, rainloss_ ! JC    bathy,
   AED_REAL, INTENT(in), DIMENSION(:),   POINTER :: air_temp_
   INTEGER,  INTENT(in), DIMENSION(:,:), POINTER :: mat_id_
   LOGICAL,  INTENT(in), DIMENSION(:),   POINTER :: active_
!
!LOCALS
   INTEGER :: i, j
   INTEGER :: nTypes, cType
   INTEGER, DIMENSION(:),ALLOCATABLE :: mat_t
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Provide pointers to arrays with environmental variables to AED2.
   print *,'set_env_aed2_models  '

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

   !# 3D variables being pointed to
   h    => h_            !# layer heights [1d array] needed for advection, diffusion
   z    => z_            !# depth [1d array], used to calculate local pressure
   extcoeff => extcoeff_ !# biogeochemical light attenuation coefficients [1d array],
                         !# output of biogeochemistry, input for physics
   salt => salt_
   temp => temp_

   rho => rho_
   tss => tss_
   active => active_

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

   dt = dt_

   IF (link_ext_par) THEN
      lpar => rad_(1,:)
   ENDIF

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
   
   print *,'aed2 re_initialize starting ... '
   CALL re_initialize()
   print *,'aed2 re_initialize completed'

   CONTAINS 
   
   !###############################################################################
   SUBROUTINE re_initialize()
   !-------------------------------------------------------------------------------
   !ARGUMENTS
   !
   !LOCALS
      TYPE(aed2_variable_t),POINTER :: tv
      AED_REAL :: tr

      AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben)
      TYPE (aed2_column_t) :: column(n_aed2_vars)

      INTEGER  :: i, col, lev, top, bot, count, nCols
   !
   !-------------------------------------------------------------------------------
   !BEGIN
      nCols = ubound(active, 1)
      DO col=1, nCols
         top = surf_map(col)
         bot = benth_map(col)
         count = top-bot+1
         CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben)
         DO lev=1, count
            CALL aed2_initialize(column, lev)
         ENDDO
      ENDDO
   END SUBROUTINE re_initialize


END SUBROUTINE set_env_aed2_models
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "Error getting variable info"

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
         print *, "ERROR: Undefined variable ", tvar%name
         err_count = err_count + 1
      ENDIF
   ENDDO

   if ( n_vars < v ) print *,"More vars than expected"
   if ( n_vars_ben < sv ) print *,"More sheet vars than expected"
   if ( n_vars_diag < sd + d ) print *,"More diag vars than expected"
   if ( n_vars_diag_sheet < sd ) print *,"More sheet diag vars than expected"

   IF ( err_count > 0 ) CALL STOPIT("*** Errors in configuration")
END SUBROUTINE check_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_column(column, col, cc, cc_diag, flux_pel, flux_atm, flux_ben)
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

      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => temp(top:bot)
            CASE ( 'salinity' )    ; column(av)%cell => salt(top:bot)
            CASE ( 'density' )     ; column(av)%cell => rho(top:bot)
            CASE ( 'layer_ht' )    ; column(av)%cell => h(top:bot)
            CASE ( 'layer_area' )  ; column(av)%cell_sheet => area(col)
            CASE ( 'rain' )        ; column(av)%cell_sheet => rain(col)! Does it mean this var is already pointed to the correct variable of TFFV? JC
            CASE ( 'rainloss' )    ; column(av)%cell_sheet => rainloss(col)! JC
            CASE ( 'material' )    ; column(av)%cell_sheet => zone(col)
            CASE ( 'bathy' )       ; column(av)%cell_sheet => bathy(col)! Does it mean this var is already pointed to the correct variable of TFFV? JC
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
         ELSE
            v = v + 1
            column(av)%cell => cc(v,top:bot)
            column(av)%flux_pel => flux_pel(v,top:bot)
            column(av)%flux_ben => flux_ben(v)
            column(av)%flux_atm => flux_atm(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE calculate_fluxes(column, count, flux_pel, flux_atm, flux_ben, h)
!-------------------------------------------------------------------------------
! Checks the current values of all state variables and repairs these
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in) :: count
   AED_REAL, INTENT(inout) :: flux_pel(:,:) !# (n_vars, n_layers)
   AED_REAL, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, INTENT(inout) :: flux_ben(:)   !# (n_vars)
   AED_REAL, INTENT(inout) :: h(:)          !# (n_layers)
!
!LOCALS
   INTEGER :: i
!-------------------------------------------------------------------------------
!BEGIN
   flux_pel = zero_
   flux_atm = zero_
   flux_ben = zero_

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
   INTEGER i,v,lev
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lev=top, bot
      !CALL aed2_equilibrate(column, lev)
      v = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               IF ( do_limiter ) THEN
                  IF ( .NOT. isnan(min_(v)) ) THEN
                     IF ( cc(v, lev) < min_(v) ) cc(v, lev) = min_(v)
                  ENDIF
                  IF ( .NOT. isnan(max_(v)) ) THEN
                     IF ( cc(v, lev) > max_(v) ) cc(v, lev) = max_(v)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO
END SUBROUTINE check_states
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

   AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben)
   TYPE (aed2_column_t) :: column(n_aed2_vars)

   INTEGER  :: i, j, col, lev, top, bot, v
   AED_REAL :: rain_loss
!
!-------------------------------------------------------------------------------
!BEGIN
 print *," START do_aed2_models"

   IF ( do_zone_averaging ) THEN
      IF (link_ext_par) THEN
         CALL calc_zone_areas(nCols, temp, salt, h, area, wnd, rho, extcoeff, I_0, par, tss, active, rain)
      ELSE
         CALL calc_zone_areas(nCols, temp, salt, h, area, wnd, rho, extcoeff, I_0, lpar, tss, active, rain)
      ENDIF
   ENDIF

!!$OMP DO
   !#--------------------------------------------------------------------
   !# LOOP THROUGH COLUMNS DOING JOBS PRIOR TO THE KINETICS BEING SOLVED
   DO col=1, nCols
      IF (.NOT. active(col)) CYCLE

      top = surf_map(col)
      bot = benth_map(col)

      CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben)

      !# Firstly run through ALL vars and select state vars
      v = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               !# only for state_vars that are not sheet
               IF ( .NOT. isnan(tv%mobility) ) THEN
                  ws(top:bot) = tv%mobility
                  CALL Settling(bot-top+1, dt, h(top:bot), ws(top:bot), Fsed_setl(col), column(i)%cell)
               ENDIF
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
!print*,"stage 2"

   IF ( do_zone_averaging ) THEN
      CALL copy_to_zone(nCols, cc, area, active, benth_map)
      CALL compute_zone_benthic_fluxes(n_aed2_vars, dt)
   ENDIF

!!$OMP DO
!print*,"stage 3"
   !#--------------------------------------------------------------------
   !# THIS IS THE MAIN WQ SOLUTION LOOP
   DO col=1, nCols

      !# find top and bottom cell indicies based on maps provided by the host
      top = surf_map(col)
      bot = benth_map(col)

      !# compute bottom shear stress for this column based on ustar from host
      col_taub = rho(bot)*(ustar_bed(col)*ustar_bed(col))

      !# set column data structure from global arrays
      CALL define_column(column, col, cc, cc_diag, flux, flux_atm, flux_ben)

!      IF (.NOT. active(col)) THEN
!         CALL aed2_calculate_riparian(column, bot-top+1, zero_);
!         CALL aed2_calculate_dry(column, bot-top+1);
!         CALL aed2_rain_loss(column, bot-top+1, rain_loss);! JC
!         CYCLE
!      ELSE
!         CALL aed2_calculate_riparian(column, bot-top+1, 1.0);
!      ENDIF


#if _NO_ODE_
      !# for this column, do the main kinetic/bgc flux calculation
      !# (this includes water column, surface and benthic interfaces)
      CALL calculate_fluxes(column, bot-top+1, flux(:,top:bot), flux_atm, flux_ben, h(top:bot))

      !# find the particles in this column and update particle bgc
      if(do_particle_bgc) CALL Particles(column, bot-top+1, h(top:bot))

      !# do riparian interfaces for this column and update fluxes
      IF ( .NOT.  Riparian(column, active(col), shadefrac(col), rainloss(col)) ) &
         CYCLE

!print*,"stage 3b"
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

      !# do non-kinetic updates to BGC variables (eq equilibration)
      CALL Update(column, bot-top+1)

      !# now the bgc updates are complete, update links to host model
      CALL BioDrag(column, bot-top+1, bio_drag(col))
      CALL BioExtinction(column, bot-top+1, extcoeff(top:bot))
      !CALL BioDensity()

      CALL check_states(column, top, bot)
   ENDDO ! cols
!!$OMP END DO

   IF ( do_zone_averaging ) &
      CALL copy_to_zone(nCols, cc, area, active, benth_map)

 print*,"DONE do_aed2_models"
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
   par_(1) = 0.45 * Io * EXP( -(Kw+localext) * zz )

   IF (count <= 1) RETURN

   DO i = 2, count
      localext = extc(i)

      !zz = zz + 0.5*h_(i)
      zz = h_(i)
      par_(i) = par_(i-1) * EXP( -(Kw+localext) * zz )
   ENDDO
END SUBROUTINE Light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Settling(N,dt,h,ww,Fsed,Y)
!-------------------------------------------------------------------------------
!
! Update settling of AED2 state variables in a given column
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)     :: N      !# number of vertical layers
   AED_REAL,INTENT(in)    :: dt     !# time step (s)
   AED_REAL,INTENT(in)    :: h(:)   !# layer thickness (m)
   AED_REAL,INTENT(in)    :: ww(:)  !# vertical advection speed
   AED_REAL,INTENT(inout) :: Fsed   !# value of sediment flux due to settling
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
   !# calculated as number of layers that the particles will travel based on settling or
   !# buoyancy velocity
   !# this number is then used to split the vertical movement calculations to limit
   !# movement across a single layer
   DO k=1,N-1
      !# sinking particles
      c=abs(ww(k+1))*dt/(0.5*(h(k+1)+h(k)))
      IF (c > cmax) cmax=c
      !# rising particles
      c=abs(ww(k))*dt/(0.5*(h(k+1)+h(k)))
      IF (c > cmax) cmax=c
   ENDDO

   it=min(itmax,int(cmax)+1)
   step_dt = dt / float(it);

   !# splitting loop
   DO i=1, it
      !# vertical loop
      DO k=1,N-1
         !# compute the slope ration
         IF (ww(k) > 0.) THEN !# Particle is rising
            Yc=Y(k  )     !# central value
         ELSE !# negative speed Particle is sinking
            Yc=Y(k+1)     !# central value
         ENDIF

         !# compute the limited flux
         cu(k)=ww(k) * Yc
      ENDDO

      !# do the upper boundary conditions
      cu(N) = zero_       !# flux into the domain from atmosphere

      !# do the lower boundary conditions
      IF (ww(1) > 0.) THEN !# Particle is rising
         cu(0) = 0.  !flux from benthos is zero
      ELSE  !# Particle is settling
         cu(0) = ww(1)*Y(1)
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
! Do riparian functionality, including operations in dry and fringing cells
! Populate feedback arrays to the hoist model associated with riparian effects
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
   IF (.NOT. actv) THEN
      CALL aed2_calculate_riparian(column, 1, zero_);
      CALL aed2_calculate_dry(column, 1);
      Riparian = .FALSE.
   ELSE
      CALL aed2_calculate_riparian(column, 1, 1.0);
   ENDIF

   !# update feedback arrays to host model, to reduce rain (or if -ve then add flow)
   CALL aed2_rain_loss(column, 1, localrainl);
   IF (link_rain_loss) rain_loss = localrainl
   

   !# update feedback arrays to shade the water (ie reduce incoming light, Io) 
   CALL aed2_light_shading(column, 1, localshade)
   IF (link_solar_shade) shade_frac = localshade

   Riparian = .TRUE.
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
   INTEGER,  INTENT(in)    :: count
!
!LOCAL VARIABLES:
   INTEGER :: i, lev
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
   IF (link_water_clarity) extc(1) = localext
   IF (count <= 1) RETURN

   DO i = 2, count
      CALL aed2_light_extinction(column, i, localext)
      IF (link_water_clarity) extc(i) = localext
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


!===============================================================================
END MODULE fv_aed2

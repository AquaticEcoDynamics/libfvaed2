
! Copyright 2014 by BMT WBM Pty Ltd under the GNU Public License - www.gnu.org

MODULE tuflowfv_external_wq

INCLUDE 'COMPILER_DIRECTIVES.FI'

! MODULE USE STATEMENTS
!DEC$ IF DEFINED(AED2)
  USE fv_aed2
!DEC$ END IF

  IMPLICIT NONE
! MODULE ACCESS
  PRIVATE
  PUBLIC :: fvwq, wqrk, wqdk
  PUBLIC :: wq
  PUBLIC :: tuflowfv_init_extern_wq, tuflowfv_construct_extern_wq, tuflowfv_do_extern_wq, tuflowfv_destruct_extern_wq

! MODULE PARAMETERS
  INTEGER,PARAMETER :: WQFileNum = 20
  INTEGER,PARAMETER :: wqrk = 8
  INTEGER,PARAMETER :: wqdk = 8

!DEC$ IF (PLATFORM==1) ! Windows
  CHARACTER(LEN=1),PARAMETER :: slash = '\'
!DEC$ ELSE IF (PLATFORM==2) ! Linux
  CHARACTER(LEN=1),PARAMETER :: slash = '/'
!DEC$ END IF

! MODULE TYPE DEFINITIONS

! WQ TYPE
TYPE :: fvwq
    LOGICAL :: init = .FALSE.                                   ! WQ INITIALISED
    LOGICAL :: disable = .FALSE.                                ! WQ CALCULATION FLAG
    INTEGER :: typ                                              ! WQ MODEL ID
    CHARACTER(LEN=30) :: model                                  ! WQ MODEL DESCRIPTION
    LOGICAL :: updated                                          ! UPDATED STATUS
    REAL(wqdk) :: dt_update                                     ! UPDATE TIMESTEP
    REAL(wqdk) :: t_update                                      ! NEXT UPDATE TIME
    INTEGER :: NC2                                              ! NUMBER OF 2D CELLS
    INTEGER :: NC3                                              ! NUMBER OF 3D CELLS
    INTEGER :: Nwq                                              ! NUMBER OF WQ CONSTITUENTS
    INTEGER :: Nben                                             ! NUMBER OF BENTHIC WQ CONSTITUENTS
    INTEGER :: Ndiag                                            ! NUMBER OF WQ DIAGNOSTIC VARIABLES
    CHARACTER(LEN=30),ALLOCATABLE,DIMENSION(:) :: names         ! WQ PELAGIC CONSTITUENT NAMES
    CHARACTER(LEN=30),ALLOCATABLE,DIMENSION(:) :: ben_names     ! WQ BENTHIC CONSTITUENT NAMES
    CHARACTER(LEN=30),ALLOCATABLE,DIMENSION(:) :: diag_names    ! WQ DIAGNOSTIC VARIABLE NAMES
    INTEGER,POINTER,DIMENSION(:) :: surf_map                ! SURFACE CELL MAP (NC2)
    INTEGER,POINTER,DIMENSION(:) :: benth_map               ! BOTTOM/BENTHIC LAYER MAP (NC2)
    INTEGER,POINTER,DIMENSION(:) :: NL                      ! NUMBER OF LAYERS (NC2)
    INTEGER,POINTER,DIMENSION(:,:) :: mat_id                ! MATERIAL ID (NMG,NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: thick                ! CELL THICKNESS (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: depth                ! LOCAL MID-CELL DEPTH (NC3)
    REAL(wqrk),POINTER,DIMENSION(:,:) :: dcdt               ! TEMPORAL DERIVATIVE OF WQ CONSTITUENTS (NWQ,NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: sal                  ! SALINITY POINTER (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: temp                 ! TEMPERATURE POINTER (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: tss                  ! TOTAL SUSPENDED SOLIDS POINTER (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: vvel                 ! VERTICAL VELOCITIES (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: cvel                 ! CELL VELOCITIES (NC3)
    REAL(wqrk),POINTER,DIMENSION(:,:) :: par                ! NET SHORTWAVE RADIATION (NC3)
    REAL(wqrk),POINTER,DIMENSION(:,:) :: cc                 ! WQ CONSTITUENT CONCENTRATIONS (NWQ,NC3)
    REAL(wqrk),POINTER,DIMENSION(:,:) :: diag               ! DIAGNOSTIC WQ VARIABLES (NDIAG,NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: density              ! ABSOLUTE DENSITY (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: I_0                  ! NET SURFACE IRRADIANCE (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: wind                 ! 10M WINDSPEED (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: precip               ! RAIN (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: ustar_bed            ! BED FRICTION VELOCITY (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: ustar_surf           ! SURFACE FRICTION VELOCITY (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: air_temp             ! AIR TEMPERATURE (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: wv_uorb              ! Wave Stress ?
    REAL(wqrk),POINTER,DIMENSION(:) :: wv_t                 !
    ! Arrays that control feedbacks between the models
    REAL(wqrk),POINTER,DIMENSION(:) :: bioshade             ! BIOGEOCHEMICAL LIGHT EXTINCTION COEFFICIENT RETURNED FROM WQ (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: biodrag              ! ADDITIONAL DRAG ON FLOW FROM BIOLOGY, RETURNED FROM WQ (NC3)
    REAL(wqrk),POINTER,DIMENSION(:) :: solarshade           ! REDUCTION OF SOLAR RADIATION DUE TO SHADING RETURNED FROM WQ (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: rainloss             ! LOSS OF RAINFALL INTO EXPOSED SEDIMENT RETURNED FROM WQ (NC2)
    ! Variables required for AED2 dry cell models
    LOGICAL,POINTER,DIMENSION(:) :: active                  ! COLUMN ACTIVE STATUS (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: area                 ! CELL AREA (NC2)
    REAL(wqrk),POINTER,DIMENSION(:) :: bathy                ! HEIGHT OF COLUMN BOTTOM (NC2)
    ! Variables related to particle tracking
    INTEGER :: NG                                           ! Number of Particle Groups
    TYPE(partgroup),DIMENSION(:),ALLOCATABLE :: parts       ! Particle groups
END TYPE

! MODULE OBJECTS
  TYPE(fvwq) :: wq
!DEC$ ATTRIBUTES DLLEXPORT :: wq

CONTAINS

!DEC$ IF DEFINED(AED2)

!******************************************************************************
! AED2 LIBRARY INTERFACE ******************************************************
!******************************************************************************
SUBROUTINE tuflowfv_init_extern_wq(nlog, wqDir, wqMod, Nwqvars, Nwqben, Nwqdiags, varnames, bennames, diagnames)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_init_extern_wq

  IMPLICIT NONE

! SUBROUTINE ARGUMENTS
  INTEGER,INTENT(IN) :: nlog
  CHARACTER(LEN=*),INTENT(IN) :: wqDir
  CHARACTER(LEN=*),INTENT(OUT) :: wqMod
  INTEGER,INTENT(OUT) :: Nwqvars, Nwqben, Nwqdiags
  CHARACTER(LEN=*),ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: varnames, bennames, diagnames

! LOCAL VARIABLES
  CHARACTER(LEN=30) :: sub = 'tuflowfv_construct_extern_wq'
  CHARACTER(LEN=200) :: wqFilename
  INTEGER :: i
  LOGICAL :: openstat

  INQUIRE(UNIT=nlog,OPENED=openstat)

! BEGIN
  ! GET AED2 CONFIGURATION
  WRITE(*,'(a\)') 'Configuring "AED2" external module... '
  IF (openstat) WRITE(nlog,'(a\)') 'Configuring "AED2" external module... '

  i = LEN_TRIM(wqDir)
  IF (i>0) THEN
    IF (wqDir(i:i)==slash) THEN
!     WQFileName = TRIM(wqDir)//'aed.nml'
      WQFileName = TRIM(wqDir)
    ELSE
!     WQFileName = TRIM(wqDir)//slash//'aed.nml'
      WQFileName = TRIM(wqDir)//slash
    END IF
  ELSE
!   WQFileName = 'aed.nml'
    WQFileName = ''
  END IF

  wqmod = 'AED2'
  CALL init_aed2_models(WQFileNum,WQFileName,Nwqvars,Nwqben,Nwqdiags,varnames,bennames,diagnames)

  WRITE(*,'(a)') 'Successful.'
  IF (openstat) WRITE(nlog,'(a)') 'Successful.'

END SUBROUTINE tuflowfv_init_extern_wq
!******************************************************************************
SUBROUTINE tuflowfv_construct_extern_wq(nlog)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_construct_extern_wq

  IMPLICIT NONE

! SUBROUTINE ARGUMENTS
  INTEGER,INTENT(IN) :: nlog

! LOCAL VARIABLES
  CHARACTER(LEN=30) :: sub = 'tuflowfv_construct_extern_wq'
  integer i,j
  LOGICAL :: openstat

  INQUIRE(UNIT=nlog,OPENED=openstat)

! BEGIN
  ! INITIALISE AED2 LINKAGE
  WRITE(*,'(a)') 'Initialising "AED2" external module:'
  IF (openstat) WRITE(nlog,'(a)') 'Initialising "AED2" external module:'

  if ( .not. associated(wq%cc) .or. .not. associated(wq%diag) ) then
    write(*, '(a)') 'Water quality variables not associated'
    IF (openstat) write(nlog, '(a)') 'Water quality variables not associated'
    return
  endif

  ! ALLOCATE LOCAL MEMORY BLOCK FOR AED2 WQ VARIABLES
  CALL init_var_aed2_models(wq%nc3,wq%cc,wq%diag,wq%nwq,wq%nben, &
                                    wq%surf_map,wq%benth_map)

  CALL set_env_aed2_models( wq%dt_update,       &
                            ! 3D env variables
                            wq%temp,            &
                            wq%sal,             &
                            wq%density,         &
                            wq%thick,           &
                            wq%tss,             &
                            wq%par,             &
                            wq%vvel,            &
                            wq%cvel,            &
                            ! 3D feedback arrays
                            wq%bioshade,        &
                            ! 2D env variables
                            wq%area,            &
                            wq%I_0,             &
                            wq%wind,            &
                            wq%precip,          &
                            wq%air_temp,        &
                            wq%ustar_bed,       &
                            wq%ustar_surf,      &
                            wq%wv_uorb,         &
                            wq%wv_t,            &
                            wq%depth,           &
                            wq%bathy,           &
                            wq%mat_id,          &
                            wq%active,          &
                            ! 2D feedback arrays
                            wq%biodrag,         &
                            wq%solarshade,      &
                            wq%rainloss)


  WRITE(*,'(a)') 'Successful.'
  IF (openstat) WRITE(nlog,'(a)') 'Successful.'
  wq%init = .TRUE.

END SUBROUTINE tuflowfv_construct_extern_wq
!******************************************************************************
SUBROUTINE tuflowfv_do_extern_wq(nlog)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_do_extern_wq
  IMPLICIT NONE

! SUBROUTINE ARGUMENTS
  INTEGER,INTENT(IN) :: nlog

! LOCAL VARIABLES
  CHARACTER(LEN=30) :: sub = 'tuflowfv_do_extern_wq'
  LOGICAL :: openstat
!
!BEGIN

  INQUIRE(UNIT=nlog,OPENED=openstat)

  IF (.NOT. wq%init) THEN
     WRITE(*,'(a)') 'ERROR "AED2" external module not initialised.'
     IF (openstat) WRITE(nlog,'(a)') 'ERROR "AED2" external module not initialised.'
     WRITE(*,'(a)') 'Stopping in '//TRIM(sub)//'.'
     IF (openstat) WRITE(nlog,'(a)') 'Stopping in '//TRIM(sub)//'.'
     STOP
  ENDIF

  IF (do_particle_bgc) CALL set_env_particles(wq%NG,wq%parts)

  CALL do_aed2_models(wq%nc3,wq%nc2)
END SUBROUTINE tuflowfv_do_extern_wq

!******************************************************************************

SUBROUTINE tuflowfv_destruct_extern_wq(nlog)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_destruct_extern_wq
  IMPLICIT NONE

! SUBROUTINE ARGUMENTS
  INTEGER,INTENT(IN) :: nlog

! LOCAL VARIABLES
  CHARACTER(LEN=30) :: sub = 'tuflowfv_destruct_extern_wq'
  LOGICAL :: openstat

  INQUIRE(UNIT=nlog,OPENED=openstat)

! BEGIN
  IF (.NOT. wq%init) RETURN
  WRITE(*,'(a)') 'Cleaning "AED2" external module:'
  IF (openstat) WRITE(nlog,'(a)') 'Cleaning "AED2" external module:'

  CALL clean_aed2_models()

  WRITE(*,'(a)') 'Successful.'
  IF (openstat) WRITE(nlog,'(a)') 'Successful.'
  wq%init = .FALSE.

END SUBROUTINE tuflowfv_destruct_extern_wq
!******************************************************************************
! END AED2 LIBRARY INTERFACE **************************************************
!******************************************************************************


!DEC$ ELSE


!******************************************************************************
! BLANK LIBRARY INTERFACE *****************************************************
!******************************************************************************
SUBROUTINE tuflowfv_init_extern_wq(nlog, wqDir, wqMod, Nwqvars, Nwqben, Nwqdiags, varnames, bennames, diagnames)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_init_extern_wq
IMPLICIT NONE
! SUBROUTINE ARGUMENTS
INTEGER,INTENT(IN) :: nlog
CHARACTER(LEN=*),INTENT(IN) :: wqDir
CHARACTER(LEN=*),INTENT(OUT) :: wqMod
INTEGER,INTENT(OUT) :: Nwqvars, Nwqben, Nwqdiags
CHARACTER(LEN=*),ALLOCATABLE,DIMENSION(:),INTENT(OUT) :: varnames, bennames, diagnames
END SUBROUTINE tuflowfv_init_extern_wq
!******************************************************************************
SUBROUTINE tuflowfv_construct_extern_wq(nlog)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_construct_extern_wq
IMPLICIT NONE
! SUBROUTINE ARGUMENTS
INTEGER,INTENT(IN) :: nlog
! LOCAL VARIABLES
LOGICAL :: openstat
INQUIRE(UNIT=nlog,OPENED=openstat)
WRITE(*,*) 'External wq library is empty.'
WRITE(*,*) 'Replace external wq library ''tuflowfv_external_wq''.'
WRITE(*,*) 'Exiting.'
IF (openstat) WRITE(nlog,*) 'External wq library is empty.'
IF (openstat) WRITE(nlog,*) 'Replace external wq library ''tuflowfv_external_wq''.'
IF (openstat) WRITE(nlog,*) 'Exiting.'
STOP
END SUBROUTINE tuflowfv_construct_extern_wq
!******************************************************************************
SUBROUTINE tuflowfv_do_extern_wq(nlog)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_do_extern_wq
IMPLICIT NONE
! SUBROUTINE ARGUMENTS
INTEGER,INTENT(IN) :: nlog
END SUBROUTINE tuflowfv_do_extern_wq
!******************************************************************************
SUBROUTINE tuflowfv_destruct_extern_wq(nlog)
!DEC$ ATTRIBUTES DLLEXPORT :: tuflowfv_destruct_extern_wq
IMPLICIT NONE
! SUBROUTINE ARGUMENTS
INTEGER,INTENT(IN) :: nlog
END SUBROUTINE tuflowfv_destruct_extern_wq
!******************************************************************************
! END BLANK LIBRARY INTERFACE *************************************************
!******************************************************************************
!DEC$ END IF

END MODULE tuflowfv_external_wq

!###############################################################################
!#                                                                             #
!# fv_ode.F90                                                                  #
!#                                                                             #
!# Interface for FV (Finite Volume) Model to AED2 modules.                     #
!#                                                                             #
!# This is a support module to managing solution of the BGC ODE set            #
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
!# Created Apr 2015                                                            #
!#                                                                             #
!###############################################################################

#include "aed2.h"

#ifndef DEBUG
#define DEBUG      0
#endif

!###############################################################################
MODULE fv_ode
!-------------------------------------------------------------------------------
   USE aed2_common

   IMPLICIT NONE

   PUBLIC init_ode, do_ode

   !#--------------------------------------------------------------------------#
   !# Module Data

   INTEGER :: solution_method

   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: flux, flux2, flux3, flux4
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: cc1

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE init_ode(ode_type, n_vars, nCells)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: ode_type, n_vars, nCells
!
!LOCALS
   INTEGER :: rc
!
!-------------------------------------------------------------------------------
!BEGIN
   solution_method = ode_type
   SELECT CASE (solution_method)
      CASE (1)   !# euler forward
      CASE (2)   !# runge_kutta_2
         ALLOCATE(flux2(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux2)'
         ALLOCATE(cc1(n_vars, nCells),stat=rc)   ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc1)'
      CASE (3)   !# runge_kutta_4
         ALLOCATE(flux2(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux2)'
         ALLOCATE(flux3(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux3)'
         ALLOCATE(flux4(n_vars, nCells),stat=rc) ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (flux4)'
         ALLOCATE(cc1(n_vars, nCells),stat=rc)   ; IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc1)'
      CASE DEFAULT
         STOP "no valid solution_method specified!"
      END SELECT
END SUBROUTINE init_ode
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE do_ode(col, top, bot)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: col, top, bot
!
!LOCALS
   INTEGER  :: lev, i
   AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben)
   TYPE (aed2_column_t) :: column(n_aed2_vars)
   TYPE (aed2_column_t) :: column2(n_aed2_vars)
   TYPE (aed2_column_t) :: column3(n_aed2_vars)
   TYPE (aed2_column_t) :: column4(n_aed2_vars)
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Time-integrate one biological time step
   SELECT CASE (solution_method)
      CASE (1)   !# This is what the euler forward would do ....
         CALL calculate_fluxes(column, bot-top+1, flux(:,top:bot), flux_atm, flux_ben, h(top:bot))
         DO lev = top, bot
            DO i = 1, n_vars
               cc(i,lev)=cc(i,lev)+dt*flux(i,lev)
#if DEBUG>1
               IF ( isnan(cc(i,lev)) ) THEN
                  print*,'Nan at i = ', i, ' lev = ', lev
                  print*,'h(lev) = ', h(lev), ' flux(i,lev) = ', flux(i,lev)
                  print*,'Top of column @ ', top, ' bottom of column @ ', bot
                  call STOPIT('NaN value')
               ENDIF
#endif
            ENDDO
         ENDDO
      CASE (2)   !# This is what the runge_kutta_2 would do ....
         CALL calculate_fluxes(column, bot-top+1, flux(:,top:bot), flux_atm, flux_ben, h(top:bot))

         DO lev = top, bot
            DO i = 1, n_vars
               cc1(i,lev)=cc(i,lev)+dt*flux(i,lev)
            ENDDO
         ENDDO

         CALL define_column(column2, col, cc1, cc_diag, flux2, flux_atm, flux_ben)
         CALL calculate_fluxes(column2, bot-top+1, flux2(:,top:bot), flux_atm, flux_ben, h(top:bot))

         DO lev = top, bot
            DO i = 1, n_vars
               cc(i,lev)=cc(i,lev)+dt*0.5*(flux(i,lev)+flux1(i,lev))
            ENDDO
         ENDDO
      CASE (3)   !# This is what the runge_kutta_4 would do ....
         CALL calculate_fluxes(column, bot-top+1, flux(:,top:bot), flux_atm, flux_ben, h(top:bot))

         DO lev = top, bot
            DO i = 1, n_vars
               cc1(i,lev)=cc(i,lev)+dt*flux(i,lev)
            ENDDO
         ENDDO

         CALL define_column(column2, col, cc1, cc_diag, flux2, flux_atm, flux_ben)
         CALL calculate_fluxes(column2, bot-top+1, flux2(:,top:bot), flux_atm, flux_ben, h(top:bot))

         DO lev = top, bot
            DO i = 1, n_vars
               cc1(i,lev)=cc(i,lev)+dt*flux1(i,lev)
            ENDDO
         ENDDO

         CALL define_column(column3, col, cc1, cc_diag, flux3, flux_atm, flux_ben)
         CALL calculate_fluxes(column3, bot-top+1, flux3(:,top:bot), flux_atm, flux_ben, h(top:bot))

         DO lev = top, bot
            DO i = 1, n_vars
               cc1(i,lev)=cc(i,lev)+dt*flux2(i,lev)
            ENDDO
         ENDDO

         CALL define_column(column4, col, cc1, cc_diag, flux4, flux_atm, flux_ben)
         CALL calculate_fluxes(column4, bot-top+1, flux4(:,top:bot), flux_atm, flux_ben, h(top:bot))

         DO lev = top, bot
            DO i = 1, n_vars
               cc(i,lev)=cc(i,lev)+dt*1./3.*(0.5*flux(i,lev)+flux2(i,lev)+flux3(i,lev)+0.5*flux4(i,lev))
            ENDDO
         ENDDO

      CASE DEFAULT
         STOP "no valid solution_method specified!"
   END SELECT
END SUBROUTINE do_ode
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE fv_ode

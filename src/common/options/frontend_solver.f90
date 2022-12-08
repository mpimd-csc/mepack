!
! THIS PROGRAM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
! IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
! THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
! (AT YOUR OPTION) ANY LATER VERSION.
!
! THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
! BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
! MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
! GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!
! YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
! ALONG WITH THIS PROGRAM; IF NOT, SEE <HTTP://WWW.GNU.ORG/LICENSES/>.
!
! Copyright (C) Martin Koehler, 2017-2022
!

!> @addtogroup frontend_solver
!> @{

!> @brief Select the solver for the frontend methods.
!>
!> The MEPACK_OPTIONS_FRONTEND_SOLVER module contains the selection of the solvers used inside the frontend routines.
MODULE MEPACK_OPTIONS_FRONTEND_SOLVER
    USE ISO_C_BINDING
    IMPLICIT NONE

    !> Select the default frontend solver. (LEVEL-3)
    INTEGER, PARAMETER :: FRONTEND_SOLVER_DEFAULT = 1

    !> Select the level-3 frontend solver
    INTEGER, PARAMETER :: FRONTEND_SOLVER_LEVEL3 = 1

    !> Select the level-2 frontend solver.
    INTEGER, PARAMETER :: FRONTEND_SOLVER_LEVEL2 = 2

    !> Select the DAG frontend solver.
    INTEGER, PARAMETER :: FRONTEND_SOLVER_DAG = 3

    !> Select the 2STAG frontend solver.
    INTEGER, PARAMETER :: FRONTEND_SOLVER_2STAGE = 4

    !> Select the recursive blocking solver.
    INTEGER, PARAMETER :: FRONTEND_SOLVER_RECURSIVE = 5

    !> Select the Gardiner-Laub solver, only TGSYLV
    INTEGER, PARAMETER :: FRONTEND_SOLVER_GARDINER_LAUB = 6

    !> Select the solver from LAPACK, only available in DLA_GESYLV and SLA_GESYLV
    INTEGER, PARAMETER :: FRONTEND_SOLVER_LAPACK = 7

    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TRSYLV = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TRSYLV2 = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TGSYLV = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TGCSYLV = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TGCSYLV_DUAL = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TRLYAP = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TRSTEIN = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TGLYAP = 0
    INTEGER, SAVE, PRIVATE :: FRONTEND_SOLVER_TGSTEIN = 0

    CONTAINS

        !> @brief Get the frontend solver for TRSYLV.
        !> @return Return the number for the inner solver for TRSYLV
        !>
        !> The TRSYLV_FRONTEND_SOLVER function returns the solver used in the
        !> TRSYLV frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !> @see mepack_options_frontend_solver::frontend_solver_lapack
        !>
        FUNCTION TRSYLV_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TRSYLV .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TRSYLV
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TRSYLV_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TRSYLV.
        !>
        !> The TRSYLV_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TRSYLV frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !> @see mepack_options_frontend_solver::frontend_solver_lapack
        !>
        SUBROUTINE TRSYLV_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TRSYLV = FSOLVE
            ELSE
                FRONTEND_SOLVER_TRSYLV = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TRSYLV_FRONTEND_SOLVER_SET

        !> @brief Get the frontend solver for TRSYLV2.
        !> @return Return the number for the inner solver for TRSYLV2
        !>
        !> The TRSYLV2_FRONTEND_SOLVER function returns the solver used in the
        !> TRSYLV2 frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        FUNCTION TRSYLV2_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TRSYLV2 .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TRSYLV2
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TRSYLV2_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TRSYLV2.
        !>
        !> The TRSYLV2_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TRSYLV2 frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        SUBROUTINE TRSYLV2_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TRSYLV2 = FSOLVE
            ELSE
                FRONTEND_SOLVER_TRSYLV2 = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TRSYLV2_FRONTEND_SOLVER_SET

        !> @brief Get the frontend solver for TGSYLV.
        !> @return Return the number for the inner solver for TGSYLV
        !>
        !> The TGSYLV_FRONTEND_SOLVER function returns the solver used in the
        !> TGSYLV frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !> @see mepack_options_frontend_solver::frontend_solver_gardiner_laub
        !>
        FUNCTION TGSYLV_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TGSYLV .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TGSYLV
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TGSYLV_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TGSYLV.
        !>
        !> The TGSYLV_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TGSYLV frontend.
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !> @see mepack_options_frontend_solver::frontend_solver_gardiner_laub
        !>
        !>
        SUBROUTINE TGSYLV_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 6 ) THEN
                FRONTEND_SOLVER_TGSYLV = FSOLVE
            ELSE
                FRONTEND_SOLVER_TGSYLV = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TGSYLV_FRONTEND_SOLVER_SET


        !> @brief Get the frontend solver for TGCSYLV.
        !> @return Return the number for the inner solver for TGCSYLV
        !>
        !> The TGCSYLV_FRONTEND_SOLVER function returns the solver used in the
        !> TGCSYLV frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        FUNCTION TGCSYLV_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TGCSYLV .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TGCSYLV
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TGCSYLV_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TGCSYLV.
        !>
        !> The TGCSYLV_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TGCSYLV frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        SUBROUTINE TGCSYLV_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TGCSYLV = FSOLVE
            ELSE
                FRONTEND_SOLVER_TGCSYLV = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TGCSYLV_FRONTEND_SOLVER_SET

        !> @brief Get the frontend solver for TGCSYLV_DUAL.
        !> @return Return the number for the inner solver for TGCSYLV_DUAL
        !>
        !> The TGCSYLV_DUAL_FRONTEND_SOLVER function returns the solver used in the
        !> TGCSYLV frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        FUNCTION TGCSYLV_DUAL_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TGCSYLV_DUAL .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TGCSYLV_DUAL
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TGCSYLV_DUAL_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TGCSYLV_DUAL.
        !>
        !> The TGCSYLV_DUAL_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TGCSYLV frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        SUBROUTINE TGCSYLV_DUAL_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TGCSYLV_DUAL = FSOLVE
            ELSE
                FRONTEND_SOLVER_TGCSYLV_DUAL = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TGCSYLV_DUAL_FRONTEND_SOLVER_SET


        !> @brief Get the frontend solver for TRLYAP.
        !> @return Return the number for the inner solver for TRLYAP
        !>
        !> The TRLYAP_FRONTEND_SOLVER function returns the solver used in the
        !> TRLYAP frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        FUNCTION TRLYAP_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TRLYAP .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TRLYAP
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TRLYAP_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TRLYAP.
        !>
        !> The TRLYAP_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TRLYAP frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        SUBROUTINE TRLYAP_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TRLYAP = FSOLVE
            ELSE
                FRONTEND_SOLVER_TRLYAP = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TRLYAP_FRONTEND_SOLVER_SET

        !> @brief Get the frontend solver for TGLYAP.
        !> @return Return the number for the inner solver for TGLYAP
        !>
        !> The TGLYAP_FRONTEND_SOLVER function returns the solver used in the
        !> TGLYAP frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        FUNCTION TGLYAP_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TGLYAP .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TGLYAP
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TGLYAP_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TGLYAP.
        !>
        !> The TGLYAP_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TGLYAP frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        SUBROUTINE TGLYAP_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TGLYAP = FSOLVE
            ELSE
                FRONTEND_SOLVER_TGLYAP = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TGLYAP_FRONTEND_SOLVER_SET

        !> @brief Get the frontend solver for TRSTEIN.
        !> @return Return the number for the inner solver for TRSTEIN
        !>
        !> The TRSTEIN_FRONTEND_SOLVER function returns the solver used in the
        !> TRSTEIN frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        FUNCTION TRSTEIN_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TRSTEIN .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TRSTEIN
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TRSTEIN_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TRSTEIN.
        !>
        !> The TRSTEIN_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TRSTEIN frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        SUBROUTINE TRSTEIN_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TRSTEIN = FSOLVE
            ELSE
                FRONTEND_SOLVER_TRSTEIN = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TRSTEIN_FRONTEND_SOLVER_SET

        !> @brief Get the frontend solver for TGSTEIN.
        !> @return Return the number for the inner solver for TGSTEIN
        !>
        !> The TGSTEIN_FRONTEND_SOLVER function returns the solver used in the
        !> TGSTEIN frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        FUNCTION TGSTEIN_FRONTEND_SOLVER() RESULT(FSOLVE)
            IMPLICIT NONE
            INTEGER :: FSOLVE

            IF ( FRONTEND_SOLVER_TGSTEIN .GT. 0 ) THEN
                FSOLVE = FRONTEND_SOLVER_TGSTEIN
            ELSE
                FSOLVE = FRONTEND_SOLVER_DEFAULT
            END IF
        END FUNCTION TGSTEIN_FRONTEND_SOLVER

        !> @brief Set the frontend solver for TGSTEIN.
        !>
        !> The TGSTEIN_FRONTEND_SOLVER_SET subroutine set the solver used in the
        !> TGSTEIN frontend.
        !>
        !> @see mepack_options_frontend_solver::frontend_solver_default
        !> @see mepack_options_frontend_solver::frontend_solver_level3
        !> @see mepack_options_frontend_solver::frontend_solver_level2
        !> @see mepack_options_frontend_solver::frontend_solver_dag
        !> @see mepack_options_frontend_solver::frontend_solver_recursive
        !> @see mepack_options_frontend_solver::frontend_solver_2stage
        !>
        !>
        SUBROUTINE TGSTEIN_FRONTEND_SOLVER_SET(FSOLVE)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: FSOLVE

            IF ( FSOLVE .GE. 1 .AND. FSOLVE .LE. 5 ) THEN
                FRONTEND_SOLVER_TGSTEIN = FSOLVE
            ELSE
                FRONTEND_SOLVER_TGSTEIN = FRONTEND_SOLVER_DEFAULT
            END IF
        END SUBROUTINE TGSTEIN_FRONTEND_SOLVER_SET

END MODULE MEPACK_OPTIONS_FRONTEND_SOLVER
!> @}


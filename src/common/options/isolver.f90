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
! Copyright (C) Martin Koehler, 2017-2023
!

!> @addtogroup isolver
!> @{

!> @brief Select the internal solver for the level-3 algorithms.
!>
!> The MEPACK_OPTIONS_ISOLVER module contains the selection of the inner solvers
!> for the level-3 solvers.
MODULE MEPACK_OPTIONS_ISOLVER
    USE ISO_C_BINDING
    IMPLICIT NONE

    !> Select the default inner solver.
    INTEGER, PARAMETER :: ISOLVER_DEFAULT = 1
    !> Select the vectorized with local copies inner solver.
    INTEGER, PARAMETER :: ISOLVER_VECTOR_LOCAL_COPY = 1
    !> Select the local copies inner solver.
    INTEGER, PARAMETER :: ISOLVER_LOCAL_COPY = 2
    !> Select the reordered level 2 inner solver.
    INTEGER, PARAMETER :: ISOLVER_REORDER = 3
    !> Select the level 2 inner solver.
    INTEGER, PARAMETER :: ISOLVER_LEVEL2 = 4
    !> Select the naive implementation of the level 2 inner solver.
    INTEGER, PARAMETER :: ISOLVER_LEVEL2_NAIVE = 5
    !> Select the recursive blocking as inner solver.
    INTEGER, PARAMETER :: ISOLVER_RECURSIVE = 6


    INTEGER, SAVE, PRIVATE :: ISOLVER_TRSYLV = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TRSYLV2 = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TGSYLV = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TGCSYLV = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TGCSYLV_DUAL = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TRLYAP = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TRSTEIN = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TGLYAP = 0
    INTEGER, SAVE, PRIVATE :: ISOLVER_TGSTEIN = 0

    CONTAINS

        !> @brief Get the inner solver for TRSYLV.
        !> @return Return the number for the inner solver for TRSYLV
        !>
        !> The TRSYLV_ISOLVER function returns the number of the inner solver
        !> for the TRSYLV level-3 and DAG solvers.
        !>
        FUNCTION TRSYLV_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TRSYLV .GT. 0 ) THEN
                ISOLVE = ISOLVER_TRSYLV
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TRSYLV_ISOLVER

        !> @brief Set the inner solver for TRSYLV
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TRSYLV_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TRSYLV level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TRSYLV_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TRSYLV = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TRSYLV = ISOLVE
            END IF
        END SUBROUTINE

        !> @brief Get the inner solver for TRSYLV2.
        !> @return Return the number for the inner solver for TRSYLV2
        !>
        !> The TRSYLV2_ISOLVER function returns the number of the inner solver
        !> for the TRSYLV2 level-3 and DAG solvers.
        !>
        FUNCTION TRSYLV2_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TRSYLV2 .GT. 0 ) THEN
                ISOLVE = ISOLVER_TRSYLV2
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TRSYLV2_ISOLVER

        !> @brief Set the inner solver for TRSYLV2
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TRSYLV2_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TRSYLV2 level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TRSYLV2_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TRSYLV2 = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TRSYLV2 = ISOLVE
            END IF
        END SUBROUTINE

        !> @brief Get the inner solver for TGSYLV.
        !> @return Return the number for the inner solver for TGSYLV
        !>
        !> The TGSYLV_ISOLVER function returns the number of the inner solver
        !> for the TGSYLV level-3 and DAG solvers.
        !>
        FUNCTION TGSYLV_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TGSYLV .GT. 0 ) THEN
                ISOLVE = ISOLVER_TGSYLV
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TGSYLV_ISOLVER

        !> @brief Set the inner solver for TGSYLV
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TGSYLV_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TGSYLV level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TGSYLV_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TGSYLV = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TGSYLV = ISOLVE
            END IF
        END SUBROUTINE

        !> @brief Get the inner solver for TGCSYLV.
        !> @return Return the number for the inner solver for TGCSYLV
        !>
        !> The TGCSYLV_ISOLVER function returns the number of the inner solver
        !> for the TGCSYLV level-3 and DAG solvers.
        !>
        FUNCTION TGCSYLV_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TGCSYLV .GT. 0 ) THEN
                ISOLVE = ISOLVER_TGCSYLV
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TGCSYLV_ISOLVER

        !> @brief Set the inner solver for TGCSYLV
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TGCSYLV_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TGCSYLV level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TGCSYLV_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TGCSYLV = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TGCSYLV = ISOLVE
            END IF
        END SUBROUTINE

        !> @brief Get the inner solver for TGCSYLV_DUAL.
        !> @return Return the number for the inner solver for TGCSYLV_DUAL
        !>
        !> The TGCSYLV_DUAL_ISOLVER function returns the number of the inner solver
        !> for the TGCSYLV_DUAL level-3 and DAG solvers.
        !>
        FUNCTION TGCSYLV_DUAL_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TGCSYLV_DUAL .GT. 0 ) THEN
                ISOLVE = ISOLVER_TGCSYLV_DUAL
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TGCSYLV_DUAL_ISOLVER

        !> @brief Set the inner solver for TGCSYLV_DUAL
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TGCSYLV_DUAL_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TGCSYLV level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TGCSYLV_DUAL_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TGCSYLV_DUAL = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TGCSYLV_DUAL = ISOLVE
            END IF
        END SUBROUTINE


        !> @brief Get the inner solver for TRLYAP.
        !> @return Return the number for the inner solver for TRLYAP
        !>
        !> The TRLYAP_ISOLVER function returns the number of the inner solver
        !> for the TRLYAP level-3 and DAG solvers.
        !>
        FUNCTION TRLYAP_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TRLYAP .GT. 0 ) THEN
                ISOLVE = ISOLVER_TRLYAP
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TRLYAP_ISOLVER

        !> @brief Set the inner solver for TRLYAP
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TRLYAP_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TRLYAP level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TRLYAP_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TRLYAP = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TRLYAP = ISOLVE
            END IF
        END SUBROUTINE

        !> @brief Get the inner solver for TGLYAP.
        !> @return Return the number for the inner solver for TGLYAP
        !>
        !> The TGLYAP_ISOLVER function returns the number of the inner solver
        !> for the TGLYAP level-3 and DAG solvers.
        !>
        FUNCTION TGLYAP_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TGLYAP .GT. 0 ) THEN
                ISOLVE = ISOLVER_TGLYAP
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TGLYAP_ISOLVER

        !> @brief Set the inner solver for TGLYAP
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TGLYAP_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TGLYAP level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TGLYAP_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TGLYAP = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TGLYAP = ISOLVE
            END IF
        END SUBROUTINE

        !> @brief Get the inner solver for TGSTEIN.
        !> @return Return the number for the inner solver for TGSTEIN
        !>
        !> The TGSTEIN_ISOLVER function returns the number of the inner solver
        !> for the TGSTEIN level-3 and DAG solvers.
        !>
        FUNCTION TGSTEIN_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TGSTEIN .GT. 0 ) THEN
                ISOLVE = ISOLVER_TGSTEIN
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TGSTEIN_ISOLVER

        !> @brief Set the inner solver for TGSTEIN
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TGSTEIN_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TGSTEIN level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TGSTEIN_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TGSTEIN = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TGSTEIN = ISOLVE
            END IF
        END SUBROUTINE

        !> @brief Get the inner solver for TRSTEIN.
        !> @return Return the number for the inner solver for TRSTEIN
        !>
        !> The TRSTEIN_ISOLVER function returns the number of the inner solver
        !> for the TRSTEIN level-3 and DAG solvers.
        !>
        FUNCTION TRSTEIN_ISOLVER() RESULT(ISOLVE)
            IMPLICIT NONE
            INTEGER :: ISOLVE

            IF ( ISOLVER_TRSTEIN .GT. 0 ) THEN
                ISOLVE = ISOLVER_TRSTEIN
            ELSE
                ISOLVE = 1
            END IF
        END FUNCTION TRSTEIN_ISOLVER

        !> @brief Set the inner solver for TRSTEIN
        !> @param[in]   ISOLVE      Number of the inner solver
        !>
        !> The TRSTEIN_ISOLVER_SET subroutine sets the inner solver to be used
        !> by the TRSTEIN level-3 and DAG solvers.
        !> If ISOLVE is out of range the default inner solver is set.
        !>
        SUBROUTINE TRSTEIN_ISOLVER_SET(ISOLVE)
            INTEGER, INTENT(IN) :: ISOLVE
            IF ( ISOLVE .LT. 1 .OR. ISOLVE .GT. 6) THEN
                ISOLVER_TRSTEIN = ISOLVER_DEFAULT
            ELSE
                ISOLVER_TRSTEIN = ISOLVE
            END IF
        END SUBROUTINE



END MODULE MEPACK_OPTIONS_ISOLVER
!> @}

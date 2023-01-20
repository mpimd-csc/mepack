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

!> @addtogroup blocksize
!> @{

!> @brief Blocksize for Double Precision Calculations
!>
!> The MEPACK_BLOCKSIZE_DOUBLE module contains the blocksize parameters specific parameters for
!> double prevent calculations.
MODULE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE ISO_C_BINDING
    USE CSC_LUA
    IMPLICIT NONE

    INTEGER, SAVE, PRIVATE :: TRSYLV_BLOCKSIZE_MB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TRSYLV_BLOCKSIZE_NB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TRSYLV2_BLOCKSIZE_MB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TRSYLV2_BLOCKSIZE_NB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGSYLV_BLOCKSIZE_MB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGSYLV_BLOCKSIZE_NB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGCSYLV_BLOCKSIZE_MB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGCSYLV_BLOCKSIZE_NB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGCSYLV_DUAL_BLOCKSIZE_MB_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGCSYLV_DUAL_BLOCKSIZE_NB_FIXED = 0

    INTEGER, SAVE, PRIVATE :: TRLYAP_BLOCKSIZE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGLYAP_BLOCKSIZE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TRSTEIN_BLOCKSIZE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGSTEIN_BLOCKSIZE_FIXED = 0

    INTEGER, SAVE, PRIVATE :: TRLYAP_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGLYAP_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TRSTEIN_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGSTEIN_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TRSYLV_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TRSYLV2_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGSYLV_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGCSYLV_BLOCKSIZE_2STAGE_FIXED = 0
    INTEGER, SAVE, PRIVATE :: TGCSYLV_DUAL_BLOCKSIZE_2STAGE_FIXED = 0

    !> \brief Interface to access the LUA interpreter
    !> \internal
    INTERFACE GET_LUA
        !> \brief Interface for mepack_lua_state
        !> \internal
        FUNCTION GET_LUAX() BIND(C, name = "mepack_lua_state")
            IMPORT
            TYPE(C_PTR) :: GET_LUAX
        END FUNCTION
    END INTERFACE

    !> \brief Interface to check if the LUA interpreter is loaded
    !> \internal
    INTERFACE LUA_LOADED
        !> \brief Interface for mepack_lua_loaded
        !> \internal
        FUNCTION LUA_LOADEX() BIND(C, name = "mepack_lua_loaded")
            IMPORT
            INTEGER(KIND = C_INT) :: LUA_LOADEX
        END FUNCTION
    END INTERFACE


    PRIVATE :: GET_LUA, LUA_LOADED

    CONTAINS

        !> @brief Return the row blocksize for TRSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TRSYLV with respect to the row partitioning.
        !>
        !> The TRSYLV_BLOCKSIZE_MB function returns the blocksize for the TRSYLV level 3
        !> solver with respect to the number of rows.
        FUNCTION TRSYLV_BLOCKSIZE_MB(M,N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO

            IF (TRSYLV_BLOCKSIZE_MB_FIXED .GT. 0 ) THEN
                MB = TRSYLV_BLOCKSIZE_MB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trsylv_double_mb", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trsylv_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF
            END IF
        END FUNCTION TRSYLV_BLOCKSIZE_MB

        !> @brief Return the column blocksize for TRSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TRSYLV with respect to the column partitioning.
        !>
        !> The TRSYLV_BLOCKSIZE_NB function returns the blocksize for the TRSYLV level 3
        !> solver with respect to the number of columns.
        FUNCTION TRSYLV_BLOCKSIZE_NB(M,N) RESULT(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: NB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRSYLV_BLOCKSIZE_NB_FIXED .GT. 0 ) THEN
                NB = TRSYLV_BLOCKSIZE_NB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trsylv_double_nb", M, N, NB, INFO)
                    IF (NB.LT.1 .OR. INFO .NE. 0) THEN
                        NB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trsylv_double_nb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    NB = 32
                END IF
            END IF
        END FUNCTION TRSYLV_BLOCKSIZE_NB


        !> @brief Set the row blocksize for TRSYLV
        !> @param[in]  MB   New row blocksize
        !>
        !> Set the row blocksize for the TRSYLV level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSYLV_BLOCKSIZE_MB_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRSYLV_BLOCKSIZE_MB_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Set the column blocksize for TRSYLV
        !> @param[in]  NB   New column blocksize
        !>
        !> Set the column blocksize for the TRSYLV level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSYLV_BLOCKSIZE_NB_SET(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NB

            IF ( NB .GE. 0 ) THEN
                TRSYLV_BLOCKSIZE_NB_FIXED = NB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the row blocksize for TRSYLV2
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TRSYLV2 with respect to the row partitioning.
        !>
        !> The TRSYLV2_BLOCKSIZE_MB function returns the blocksize for the TRSYLV2 level 3
        !> solver with respect to the number of rows.
        FUNCTION TRSYLV2_BLOCKSIZE_MB(M,N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRSYLV2_BLOCKSIZE_MB_FIXED .GT. 0 ) THEN
                MB = TRSYLV2_BLOCKSIZE_MB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trsylv2_double_mb", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trsylv2_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF
            END IF
        END FUNCTION TRSYLV2_BLOCKSIZE_MB

        !> @brief Return the column blocksize for TRSYLV2
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TRSYLV2 with respect to the column partitioning.
        !>
        !> The TRSYLV2_BLOCKSIZE_NB function returns the blocksize for the TRSYLV2 level 3
        !> solver with respect to the number of columns.
        FUNCTION TRSYLV2_BLOCKSIZE_NB(M,N) RESULT(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: NB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRSYLV2_BLOCKSIZE_NB_FIXED .GT. 0 ) THEN
                NB = TRSYLV2_BLOCKSIZE_NB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trsylv2_double_nb", M, N, NB, INFO)
                    IF (NB.LT.1 .OR. INFO .NE. 0) THEN
                        NB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trsylv2_double_nb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    NB = 32
                END IF

            END IF
        END FUNCTION TRSYLV2_BLOCKSIZE_NB


        !> @brief Set the row blocksize for TRSYLV2
        !> @param[in]  MB   New row blocksize
        !>
        !> Set the row blocksize for the TRSYLV2 level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSYLV2_BLOCKSIZE_MB_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRSYLV2_BLOCKSIZE_MB_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Set the column blocksize for TRSYLV2
        !> @param[in]  NB   New column blocksize
        !>
        !> Set the column blocksize for the TRSYLV2 level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSYLV2_BLOCKSIZE_NB_SET(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NB

            IF ( NB .GE. 0 ) THEN
                TRSYLV2_BLOCKSIZE_NB_FIXED = NB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Return the row blocksize for TGSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGSYLV with respect to the row partitioning.
        !>
        !> The TGSYLV_BLOCKSIZE_MB function returns the blocksize for the TGSYLV level 3
        !> solver with respect to the number of rows.
        FUNCTION TGSYLV_BLOCKSIZE_MB(M,N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGSYLV_BLOCKSIZE_MB_FIXED .GT. 0 ) THEN
                MB = TGSYLV_BLOCKSIZE_MB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgsylv_double_mb", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgsylv_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF

            END IF
        END FUNCTION TGSYLV_BLOCKSIZE_MB

        !> @brief Return the column blocksize for TGSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGSYLV with respect to the column partitioning.
        !>
        !> The TGSYLV_BLOCKSIZE_NB function returns the blocksize for the TGSYLV level 3
        !> solver with respect to the number of columns.
        FUNCTION TGSYLV_BLOCKSIZE_NB(M,N) RESULT(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: NB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGSYLV_BLOCKSIZE_NB_FIXED .GT. 0 ) THEN
                NB = TGSYLV_BLOCKSIZE_NB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgsylv_double_nb", M, N, NB, INFO)
                    IF (NB.LT.1 .OR. INFO .NE. 0) THEN
                        NB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgsylv_double_nb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    NB = 32
                END IF
            END IF
        END FUNCTION TGSYLV_BLOCKSIZE_NB

        !> @brief Set the row blocksize for TGSYLV
        !> @param[in]  MB   New row blocksize
        !>
        !> Set the row blocksize for the TGSYLV level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGSYLV_BLOCKSIZE_MB_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGSYLV_BLOCKSIZE_MB_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Set the column blocksize for TGSYLV
        !> @param[in]  NB   New column blocksize
        !>
        !> Set the column blocksize for the TGSYLV level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGSYLV_BLOCKSIZE_NB_SET(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NB

            IF ( NB .GE. 0 ) THEN
                TGSYLV_BLOCKSIZE_NB_FIXED = NB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Return the row blocksize for TGCSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGCSYLV with respect to the row partitioning.
        !>
        !> The TGCSYLV_BLOCKSIZE_MB function returns the blocksize for the TGCSYLV level 3
        !> solver with respect to the number of rows.
        FUNCTION TGCSYLV_BLOCKSIZE_MB(M,N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGCSYLV_BLOCKSIZE_MB_FIXED .GT. 0 ) THEN
                MB = TGCSYLV_BLOCKSIZE_MB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgcsylv_double_mb", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgcsylv_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF

            END IF
        END FUNCTION TGCSYLV_BLOCKSIZE_MB

        !> @brief Return the column blocksize for TGCSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGCSYLV with respect to the column partitioning.
        !>
        !> The TGCSYLV_BLOCKSIZE_NB function returns the blocksize for the TGCSYLV level 3
        !> solver with respect to the number of columns.
        FUNCTION TGCSYLV_BLOCKSIZE_NB(M,N) RESULT(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: NB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGCSYLV_BLOCKSIZE_NB_FIXED .GT. 0 ) THEN
                NB = TGCSYLV_BLOCKSIZE_NB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgcsylv_double_nb", M, N, NB, INFO)
                    IF (NB.LT.1 .OR. INFO .NE. 0) THEN
                        NB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgcsylv_double_nb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    NB = 32
                END IF
            END IF
        END FUNCTION TGCSYLV_BLOCKSIZE_NB

        !> @brief Set the row blocksize for TGCSYLV
        !> @param[in]  MB   New row blocksize
        !>
        !> Set the row blocksize for the TGCSYLV level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGCSYLV_BLOCKSIZE_MB_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGCSYLV_BLOCKSIZE_MB_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Set the column blocksize for TGCSYLV
        !> @param[in]  NB   New column blocksize
        !>
        !> Set the column blocksize for the TGCSYLV level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGCSYLV_BLOCKSIZE_NB_SET(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NB

            IF ( NB .GE. 0 ) THEN
                TGCSYLV_BLOCKSIZE_NB_FIXED = NB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Return the row blocksize for TGCSYLV_DUAL
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGCSYLV_DUAL with respect to the row partitioning.
        !>
        !> The TGCSYLV_DUAL_BLOCKSIZE_MB function returns the blocksize for the TGCSYLV_DUAL level 3
        !> solver with respect to the number of rows.
        FUNCTION TGCSYLV_DUAL_BLOCKSIZE_MB(M,N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGCSYLV_DUAL_BLOCKSIZE_MB_FIXED .GT. 0 ) THEN
                MB = TGCSYLV_DUAL_BLOCKSIZE_MB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgcsylv_dual_double_mb", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgcsylv_dual_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF

            END IF
        END FUNCTION TGCSYLV_DUAL_BLOCKSIZE_MB

        !> @brief Return the column blocksize for TGCSYLV_DUAL
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGCSYLV_DUAL with respect to the column partitioning.
        !>
        !> The TGCSYLV_DUAL_BLOCKSIZE_NB function returns the blocksize for the TGCSYLV_DUAL level 3
        !> solver with respect to the number of columns.
        FUNCTION TGCSYLV_DUAL_BLOCKSIZE_NB(M,N) RESULT(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: NB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGCSYLV_DUAL_BLOCKSIZE_NB_FIXED .GT. 0 ) THEN
                NB = TGCSYLV_DUAL_BLOCKSIZE_NB_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgcsylv_dual_double_nb", M, N, NB, INFO)
                    IF (NB.LT.1 .OR. INFO .NE. 0) THEN
                        NB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgcsylv_dual_double_nb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    NB = 32
                END IF
            END IF
        END FUNCTION TGCSYLV_DUAL_BLOCKSIZE_NB

        !> @brief Set the row blocksize for TGCSYLV_DUAL
        !> @param[in]  MB   New row blocksize
        !>
        !> Set the row blocksize for the TGCSYLV_DUAL level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGCSYLV_DUAL_BLOCKSIZE_MB_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGCSYLV_DUAL_BLOCKSIZE_MB_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE

        !> @brief Set the column blocksize for TGCSYLV_DUAL
        !> @param[in]  NB   New column blocksize
        !>
        !> Set the column blocksize for the TGCSYLV_DUAL level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGCSYLV_DUAL_BLOCKSIZE_NB_SET(NB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NB

            IF ( NB .GE. 0 ) THEN
                TGCSYLV_DUAL_BLOCKSIZE_NB_FIXED = NB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the blocksize for TRLYAP
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TRLYAP.
        !>
        !> The TRLYAP_BLOCKSIZE function returns the blocksize for the TRLYAP level 3
        !> solver.
        FUNCTION TRLYAP_BLOCKSIZE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRLYAP_BLOCKSIZE_FIXED .GT. 0 ) THEN
                MB = TRLYAP_BLOCKSIZE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trlyap_double_mb", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trlyap_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF

            END IF
        END FUNCTION TRLYAP_BLOCKSIZE

        !> @brief Set the blocksize for TRLYAP
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TRLYAP level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRLYAP_BLOCKSIZE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRLYAP_BLOCKSIZE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the blocksize for TRSTEIN
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TRSTEIN.
        !>
        !> The TRSTEIN_BLOCKSIZE function returns the blocksize for the TRSTEIN level 3
        !> solver.
        FUNCTION TRSTEIN_BLOCKSIZE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRSTEIN_BLOCKSIZE_FIXED .GT. 0 ) THEN
                MB = TRSTEIN_BLOCKSIZE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trstein_double_mb", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trstein_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF


            END IF
        END FUNCTION TRSTEIN_BLOCKSIZE

        !> @brief Set the blocksize for TRSTEIN
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TRSTEIN level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSTEIN_BLOCKSIZE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRSTEIN_BLOCKSIZE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the blocksize for TGLYAP
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TGLYAP.
        !>
        !> The TGLYAP_BLOCKSIZE function returns the blocksize for the TGLYAP level 3
        !> solver.
        FUNCTION TGLYAP_BLOCKSIZE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGLYAP_BLOCKSIZE_FIXED .GT. 0 ) THEN
                MB = TGLYAP_BLOCKSIZE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tglyap_double_mb", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tglyap_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF


            END IF
        END FUNCTION TGLYAP_BLOCKSIZE

        !> @brief Set the blocksize for TGLYAP
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TGLYAP level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGLYAP_BLOCKSIZE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGLYAP_BLOCKSIZE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the blocksize for TGSTEIN
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TGSTEIN.
        !>
        !> The TGSTEIN_BLOCKSIZE function returns the blocksize for the TGSTEIN level 3
        !> solver.
        FUNCTION TGSTEIN_BLOCKSIZE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGSTEIN_BLOCKSIZE_FIXED .GT. 0 ) THEN
                MB = TGSTEIN_BLOCKSIZE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgstein_double_mb", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 32
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgstein_double_mb) returned a blocksize lower than 1 or failed. Use 32."
                    END IF
                ELSE
                    MB = 32
                END IF


            END IF
        END FUNCTION TGSTEIN_BLOCKSIZE

        !> @brief Set the blocksize for TGSTEIN
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TGSTEIN level-3 and DAG algorithms. If it
        !> is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGSTEIN_BLOCKSIZE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGSTEIN_BLOCKSIZE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the big blocksize for 2-stage TRLYAP
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TRLYAP.
        !>
        !> The TRLYAP_BLOCKSIZE function returns the blocksize for the two stage
        !> TRLYAP level 3 solver.
        FUNCTION TRLYAP_BLOCKSIZE_2STAGE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRLYAP_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TRLYAP_BLOCKSIZE_2STAGE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trlyap_double_2stage", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trlyap_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF


            END IF
        END FUNCTION TRLYAP_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TRLYAP two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TRLYAP level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRLYAP_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRLYAP_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the big blocksize for 2-stage TRSTEIN
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TRSTEIN.
        !>
        !> The TRSTEIN_BLOCKSIZE function returns the blocksize for the two stage
        !> TRSTEIN level 3 solver.
        FUNCTION TRSTEIN_BLOCKSIZE_2STAGE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRSTEIN_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TRSTEIN_BLOCKSIZE_2STAGE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trstein_double_2stage", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trstein_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF
            END IF
        END FUNCTION TRSTEIN_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TRSTEIN two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TRSTEIN level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSTEIN_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRSTEIN_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the big blocksize for 2-stage TGLYAP
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TGLYAP.
        !>
        !> The TGLYAP_BLOCKSIZE function returns the blocksize for the two stage
        !> TGLYAP level 3 solver.
        FUNCTION TGLYAP_BLOCKSIZE_2STAGE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGLYAP_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TGLYAP_BLOCKSIZE_2STAGE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tglyap_double_2stage", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tglyap_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF
            END IF
        END FUNCTION TGLYAP_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TGLYAP two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TGLYAP level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGLYAP_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGLYAP_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE



        !> @brief Return the big blocksize for 2-stage TGSTEIN
        !> @param[in]   M       Number of rows
        !> @return The blocksize for TGSTEIN.
        !>
        !> The TGSTEIN_BLOCKSIZE function returns the blocksize for the two stage
        !> TGSTEIN level 3 solver.
        FUNCTION TGSTEIN_BLOCKSIZE_2STAGE(M) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGSTEIN_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TGSTEIN_BLOCKSIZE_2STAGE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgstein_double_2stage", M, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgstein_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF


            END IF
        END FUNCTION TGSTEIN_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TGSTEIN two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TGSTEIN level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGSTEIN_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGSTEIN_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the big blocksize for 2-stage TRSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TRSYLV.
        !>
        !> The TRSYLV_BLOCKSIZE function returns the blocksize for the two stage
        !> TRSYLV level 3 solver.
        FUNCTION TRSYLV_BLOCKSIZE_2STAGE(M, N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRSYLV_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TRSYLV_BLOCKSIZE_2STAGE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trsylv_double_2stage", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trsylv_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF

            END IF
        END FUNCTION TRSYLV_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TRSYLV two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TRSYLV level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSYLV_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRSYLV_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the big blocksize for 2-stage TRSYLV2
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TRSYLV2.
        !>
        !> The TRSYLV2_BLOCKSIZE function returns the blocksize for the two stage
        !> TRSYLV2 level 3 solver.
        FUNCTION TRSYLV2_BLOCKSIZE_2STAGE(M, N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TRSYLV2_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TRSYLV2_BLOCKSIZE_2STAGE_FIXED
            ELSE
                IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "trsylv2_double_2stage", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(trsylv2_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF
            END IF
        END FUNCTION TRSYLV2_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TRSYLV2 two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TRSYLV2 level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TRSYLV2_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TRSYLV2_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the big blocksize for 2-stage TGSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of rows
        !> @return The blocksize for TGSYLV.
        !>
        !> The TGSYLV_BLOCKSIZE function returns the blocksize for the two stage
        !> TGSYLV level 3 solver.
        FUNCTION TGSYLV_BLOCKSIZE_2STAGE(M, N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGSYLV_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TGSYLV_BLOCKSIZE_2STAGE_FIXED
            ELSE
                 IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgsylv_double_2stage", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgsylv_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF


            END IF
        END FUNCTION TGSYLV_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TGSYLV two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TGSYLV level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGSYLV_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGSYLV_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE


        !> @brief Return the big blocksize for 2-stage TGCSYLV
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGCSYLV.
        !>
        !> The TGCSYLV_BLOCKSIZE function returns the blocksize for the two stage
        !> TGCSYLV level 3 solver.
        FUNCTION TGCSYLV_BLOCKSIZE_2STAGE(M,N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGCSYLV_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TGCSYLV_BLOCKSIZE_2STAGE_FIXED
            ELSE
                 IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgcsylv_double_2stage", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgcsylv_double_2stage) returned a blocksize lower than 1 or failed. Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF

            END IF
        END FUNCTION TGCSYLV_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TGCSYLV two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TGCSYLV level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGCSYLV_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGCSYLV_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE



        !> @brief Return the big blocksize for 2-stage TGCSYLV_DUAL
        !> @param[in]   M       Number of rows
        !> @param[in]   N       Number of columns
        !> @return The blocksize for TGCSYLV_DUAL.
        !>
        !> The TGCSYLV_DUAL_BLOCKSIZE function returns the blocksize for the two stage
        !> TGCSYLV_DUAL level 3 solver.
        FUNCTION TGCSYLV_DUAL_BLOCKSIZE_2STAGE(M,N) RESULT(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: M, N
            INTEGER :: MB
            TYPE(C_PTR) :: LUA
            INTEGER :: INFO


            IF (TGCSYLV_DUAL_BLOCKSIZE_2STAGE_FIXED .GT. 0 ) THEN
                MB = TGCSYLV_DUAL_BLOCKSIZE_2STAGE_FIXED
            ELSE
                 IF ( LUA_LOADED() .EQ. 1 ) THEN
                    LUA = GET_LUA()
                    CALL CSC_LUA_CALL_RETURN_INT(LUA, "tgcsylv_dual_double_2stage", M, N, MB, INFO)
                    IF (MB.LT.1 .OR. INFO .NE. 0) THEN
                        MB = 2048
                        WRITE(*,*) "MEPACK_LUA_CONFIG(tgcsylv_dual_double_2stage) returned a blocksize lower than 1 or failed."
                        WRITE(*,*) "--> Use 2048."
                    END IF
                ELSE
                    MB = 2048
                END IF

            END IF
        END FUNCTION TGCSYLV_DUAL_BLOCKSIZE_2STAGE

        !> @brief Set the blocksize for TGCSYLV_DUAL two-stage algorithms
        !> @param[in]  MB   New blocksize
        !>
        !> Set the blocksize for the TGCSYLV_DUAL level-3 two stage algorithm,
        !> If it is set to 0 the automatic selection will be enabled.
        !>
        SUBROUTINE TGCSYLV_DUAL_BLOCKSIZE_2STAGE_SET(MB)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: MB

            IF ( MB .GE. 0 ) THEN
                TGCSYLV_DUAL_BLOCKSIZE_2STAGE_FIXED = MB
            END IF
            RETURN
        END SUBROUTINE



END MODULE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
!> @}

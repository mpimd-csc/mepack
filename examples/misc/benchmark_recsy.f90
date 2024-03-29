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
! COPYRIGHT (C) Martin Koehler, 2017-2022
!

PROGRAM MAIN
    USE CSC_COMMON
    USE CSC_TABLE
    USE CSC_INIFILE
    USE OMP_LIB
    IMPLICIT NONE


    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A, B, C, D, X, Y, YBAK
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK

    DOUBLE PRECISION, DIMENSION(20) :: MACHINE

    INTEGER :: M, MSTART, MSTOP, MSTEP, RUNS, R, MAT, NMAT
    INTEGER, DIMENSION(4) :: ISEED, ISEED2
    INTEGER :: SORT
    CHARACTER(1) :: TRANSA, TRANSB
    DOUBLE PRECISION, PARAMETER :: DONE = 1.0D0
    DOUBLE PRECISION :: SCAL
    INTEGER :: INFO, LWORK
    INTEGER, PARAMETER :: NB = 64, MB = 64
    DOUBLE PRECISION :: TS, TE

    CHARACTER(LEN=1024) :: OUTPUT

    EXTERNAL BENCHMARK_INIT_F
    EXTERNAL BENCHMARK_RANDOM_GEVP_DOUBLE_F
    EXTERNAL BENCHMARK_RHS_GSYLV_DOUBLE_F
    EXTERNAL BENCHMARK_CHECK_X_DOUBLE_F
    DOUBLE PRECISION BENCHMARK_CHECK_X_DOUBLE_F
    EXTERNAL DGEMM
    EXTERNAL DLA_TGSYLV_L3
    EXTERNAL RECGSYL, RECSY_MACHINE
    EXTERNAL RECGSYL_P
    EXTERNAL MEPACK_INITIALIZE
    INTEGER :: THREADS
    ! Table
    TYPE(C_PTR) :: TAB
    INTEGER :: COLM
    INTEGER :: COLRECSY, COLRECSYP, COLSPD
    INTEGER :: COLRESRECSY, COLRESRECSYP


    ! TIMES
    DOUBLE PRECISION :: TRECSY, TRECSYP, T1, T2
    DOUBLE PRECISION :: RRECSY, RRECSYP

    CALL MEPACK_INITIALIZE()
    TAB = CSC_TABLE_NEW(1)
    COLM          = CSC_TABLE_ADD_COLUMN(TAB, "M", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT)
    COLRECSY      = CSC_TABLE_ADD_COLUMN(TAB, "Time RECSY", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLRESRECSY   = CSC_TABLE_ADD_COLUMN(TAB, "Res. RECSY", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLRECSYP     = CSC_TABLE_ADD_COLUMN(TAB, "Time RECSY Par.", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLRESRECSYP  = CSC_TABLE_ADD_COLUMN(TAB, "Res. RECSY Par.", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLSPD        = CSC_TABLE_ADD_COLUMN(TAB, "Speed Up.", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    CALL CSC_TABLE_COMMENT_ALLINFO(TAB)

    CALL RECSY_MACHINE(MACHINE)
    THREADS = OMP_GET_MAX_THREADS()


    MSTART = 512
    MSTOP  = 6144
    MSTEP  = 512
    RUNS   = 10
    NMAT = 5
    MAT = 0

    ISEED = 1
    SORT = 0

    TRANSA = 'N'
    TRANSB = 'N'

    WRITE(*,'("# ", A )') "Benchmark the RECSY and parallel RECSY codes for the generalized Sylvester Equation"

    CALL BENCHMARK_INIT_F()

    DO M = MSTART, MSTOP, MSTEP
        ISEED = 1
        TRECSY  = 0.0D0
        TRECSYP = 0.0D0
        RRECSY  = 0.0D0
        RRECSYP = 0.0D0

        CALL CSC_TABLE_NEW_ROW(TAB)

        DO MAT = 1, NMAT
            ALLOCATE(A(M,M), B(M,M), C(M,M), D(M,M), X(M,M), Y(M,M), YBAK(M, M) )
            X = DONE
            LWORK =  MAX(2*MAX(NB,MB),MAX(MB*NB, M*NB), M*M+4)
            ALLOCATE(WORK(LWORK))

            MACHINE(4) = MB
            MACHINE(5) = NB

            ISEED2 = ISEED
            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED, A, C, SORT)
            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED, B, D, SORT)

            CALL BENCHMARK_RHS_GSYLV_DOUBLE_F(TRANSA, TRANSB, DONE, M, M, A, M, B, M, C, M, D, M, X, M, Y, M)
            YBAK = Y

            ! L2 Aligned Copy
            TE = 0.0D0
            DO R = 1, RUNS
                Y = YBAK
                TS = CSC_WTIME()
                CALL RECGSYL(0, SCAL, M , M, A, M, B, M, C, M, D, M, Y, M, INFO, MACHINE, WORK, LWORK)
                TE = TE + (CSC_WTIME() - TS)
            END DO
            TRECSY = TRECSY + TE / DBLE (RUNS)
            RRECSY = RRECSY + BENCHMARK_CHECK_X_DOUBLE_F(M, M, Y, M, X, M)


            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED2, A, C, MB)
            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED2, B, D, NB)
            ! CALL DGEMM("N", "N", M, M, M, DONE, A, M, C, M, DONE, Y, M)

            CALL BENCHMARK_RHS_GSYLV_DOUBLE_F(TRANSA, TRANSB, DONE, M, M, A, M, B, M, C, M, D, M, X, M, Y, M)
            YBAK = Y

            ! L2 Aligned Copy
            MACHINE(6) = 3D6
            TE = 0.0D0
            DO R = 1, RUNS
                Y = YBAK
                TS = CSC_WTIME()
                INFO = 0
                CALL RECGSYL_P(THREADS, 0, SCAL, M , M, A, M, B, M, C, M, D, M, Y, M, INFO, MACHINE, WORK, LWORK)
                TE = TE + (CSC_WTIME() - TS)
            END DO
            TRECSYP = TRECSYP + TE / DBLE (RUNS)
            RRECSYP = RRECSYP + BENCHMARK_CHECK_X_DOUBLE_F(M, M, Y, M, X, M)



            DEALLOCATE(A, B, C, D, X, Y, YBAK)
            DEALLOCATE(WORK)

        END DO

        T1 = TRECSY/DBLE(NMAT)
        T2 = TRECSYP/DBLE(NMAT)
        CALL CSC_TABLE_SET_ENTRY_INTEGER(TAB, COLM, M)
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLRECSY, TRECSY/DBLE(NMAT))
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLRESRECSY, RRECSY/DBLE(NMAT))
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLRECSYP, TRECSYP/DBLE(NMAT))
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLRESRECSYP, RRECSYP/DBLE(NMAT))
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLSPD, T1 / T2)
    END DO

    CALL CSC_TABLE_COMMENT_CMD(TAB)
    CALL CSC_TABLE_COMMENT_DATE(TAB)
    CALL CSC_TABLE_COMMENT_TEXT(TAB, "MSTART:  " // INT2STR(MSTART))
    CALL CSC_TABLE_COMMENT_TEXT(TAB, "MSTOP:   " // INT2STR(MSTOP))
    CALL CSC_TABLE_COMMENT_TEXT(TAB, "MSTEP:   " // INT2STR(MSTEP))
    CALL CSC_TABLE_COMMENT_TEXT(TAB, "RUNS:    " // INT2STR(RUNS))
    CALL CSC_TABLE_COMMENT_TEXT(TAB, "NMAT:    " // INT2STR(NMAT))


    CALL CSC_TABLE_PRINT(TAB, "    ")
    IF ( COMMAND_ARGUMENT_COUNT() .EQ. 1) THEN
        CALL GET_COMMAND_ARGUMENT(1, OUTPUT)
        CALL CSC_TABLE_SAVE_ASCII(TRIM(OUTPUT), TAB, "   ")
    END IF

    CALL CSC_TABLE_DESTROY(TAB)


END PROGRAM MAIN




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
    IMPLICIT NONE


    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A, B, C, D, X, Y, YBAK
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK


    INTEGER :: M, MSTART, MSTOP, MSTEP, RUNS, R, MAT, NMAT
    INTEGER, DIMENSION(4) :: ISEED, ISEED2
    INTEGER :: SORT
    CHARACTER(1) :: TRANSA, TRANSB
    DOUBLE PRECISION, PARAMETER :: DONE = 1.0D0
    DOUBLE PRECISION :: SCAL
    INTEGER :: INFO, LWORK
    INTEGER, PARAMETER :: NB = 64, MB = 64
    DOUBLE PRECISION :: TS, TE, T1, T2

    CHARACTER(LEN=1024) :: OUTPUT

    EXTERNAL BENCHMARK_INIT_F
    EXTERNAL BENCHMARK_RANDOM_GEVP_DOUBLE_F
    EXTERNAL BENCHMARK_RHS_GSYLV_DOUBLE_F
    EXTERNAL BENCHMARK_CHECK_X_DOUBLE_F
    DOUBLE PRECISION BENCHMARK_CHECK_X_DOUBLE_F
    EXTERNAL DGEMM
    EXTERNAL DLA_TGSYLV_L3
    EXTERNAL MEPACK_INITIALIZE
    ! Table
    TYPE(C_PTR) :: TAB
    INTEGER :: COLM
    INTEGER :: COLL2ACOPY, COLRESL2ACOPY
    INTEGER :: COLL2ACOPYA, COLRESL2ACOPYA
    INTEGER :: COLSPD

    ! TIMES
    DOUBLE PRECISION :: TL2ACOPY, TL2ACOPYA
    DOUBLE PRECISION :: RL2ACOPY, RL2ACOPYA


    TAB = CSC_TABLE_NEW(1)
    COLM           = CSC_TABLE_ADD_COLUMN(TAB, "M", CSC_TABLE_INTEGER, CSC_TABLE_RIGHT)
    COLL2ACOPY      = CSC_TABLE_ADD_COLUMN(TAB, "Time No Reorder", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLRESL2ACOPY   = CSC_TABLE_ADD_COLUMN(TAB, "Res. No Reorder", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLL2ACOPYA      = CSC_TABLE_ADD_COLUMN(TAB, "Time Reorder", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLRESL2ACOPYA   = CSC_TABLE_ADD_COLUMN(TAB, "Res. Reorder", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    COLSPD           = CSC_TABLE_ADD_COLUMN(TAB, "Speed Up", CSC_TABLE_FLOAT, CSC_TABLE_RIGHT)
    CALL CSC_TABLE_COMMENT_ALLINFO(TAB)

    CALL MEPACK_INITIALIZE()

    MSTART = 10
    MSTOP  = 2100
    MSTEP  = 10
    RUNS   = 10
    NMAT = 5
    MAT = 0

    ISEED = 1
    SORT = 0

    TRANSA = 'N'
    TRANSB = 'N'

    WRITE(*,'("# ", A )') "Benchmark the influence of reordered eigenvalues on the execution time of the"
    WRITE(*,'("# ", A )') "aligned local copy solvers for the generalized Sylvester equation."

    CALL BENCHMARK_INIT_F()

    DO M = MSTART, MSTOP, MSTEP
        ISEED = 1
        TL2ACOPY = 0.0D0
        TL2ACOPYA = 0.0D0
        RL2ACOPY = 0.0D0
        RL2ACOPYA = 0.0D0

        CALL CSC_TABLE_NEW_ROW(TAB)

        DO MAT = 1, NMAT
            ALLOCATE(A(M,M), B(M,M), C(M,M), D(M,M), X(M,M), Y(M,M), YBAK(M, M) )
            X = DONE
            LWORK =  MAX(2*MAX(NB,MB),MAX(MB*NB, M*NB))
            ALLOCATE(WORK(LWORK))

            ISEED2 = ISEED
            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED, A, C, SORT)
            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED, D, B, SORT)
            ! CALL DGEMM("N", "N", M, M, M, DONE, A, M, C, M, DONE, Y, M)

            CALL BENCHMARK_RHS_GSYLV_DOUBLE_F(TRANSA, TRANSB, DONE, M, M, A, M, B, M, C, M, D, M, X, M, Y, M)
            YBAK = Y

            ! Without reordering
            TE = 0.0D0
            DO R = 1, RUNS
                Y = YBAK
                TS = CSC_WTIME()
                CALL DLA_TGSYLV_L3(TRANSA, TRANSB, DONE, M, M, A, M, B, M, C, M, D, M, Y, M, SCAL, WORK, INFO)
                TE = TE + (CSC_WTIME() - TS)
            END DO
            TL2ACOPY = TL2ACOPY + TE / DBLE (RUNS)
            RL2ACOPY = RL2ACOPY + BENCHMARK_CHECK_X_DOUBLE_F(M, M, Y, M, X, M)


            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED2, A, C, MB)
            CALL BENCHMARK_RANDOM_GEVP_DOUBLE_F(M, ISEED2, D, B, NB)
            ! CALL DGEMM("N", "N", M, M, M, DONE, A, M, C, M, DONE, Y, M)

            CALL BENCHMARK_RHS_GSYLV_DOUBLE_F(TRANSA, TRANSB, DONE, M, M, A, M, B, M, C, M, D, M, X, M, Y, M)
            YBAK = Y

            ! With reordering
            TE = 0.0D0
            DO R = 1, RUNS
                Y = YBAK
                TS = CSC_WTIME()
                CALL DLA_TGSYLV_L3(TRANSA, TRANSB, DONE, M, M, A, M, B, M, C, M, D, M, Y, M, SCAL, WORK, INFO)
                TE = TE + (CSC_WTIME() - TS)
            END DO
            TL2ACOPYA = TL2ACOPYA + TE / DBLE (RUNS)
            RL2ACOPYA = RL2ACOPYA + BENCHMARK_CHECK_X_DOUBLE_F(M, M, Y, M, X, M)



            DEALLOCATE(A, B, C, D, X, Y, YBAK)
            DEALLOCATE(WORK)

        END DO

        T1 =  TL2ACOPY/DBLE(NMAT)
        T2 = TL2ACOPYA/DBLE(NMAT)

        CALL CSC_TABLE_SET_ENTRY_INTEGER(TAB, COLM, M)

        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLL2ACOPY,T1)
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLRESL2ACOPY, RL2ACOPY/DBLE(NMAT))
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLL2ACOPYA, T2)
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLRESL2ACOPYA, RL2ACOPYA/DBLE(NMAT))
        CALL CSC_TABLE_SET_ENTRY_FLOAT(TAB, COLSPD, T1/T2)
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




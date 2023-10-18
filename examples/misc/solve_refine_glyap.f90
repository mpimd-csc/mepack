!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) Martin Koehler, 2017-2023

!
! Solve a generalized Lyapunov Equation with iterative refinement.
!
PROGRAM RUN_REFINE_GLYAP
    USE CSC_HDF5
    USE CSC_COMMON
    USE CSC_SYSINFO
    USE MEPACK_OPTIONS_FRONTEND_SOLVER

    IMPLICIT NONE


    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A, C, AS, CS, Y, X, Q, Z, W
    DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
    INTEGER(CSC_HDF5_T) :: H5FILE1
    DOUBLE PRECISION :: TS, TE, TALL
    DOUBLE PRECISION, DIMENSION(100) :: CONVLOG
    INTEGER :: NARGS, M1, M2, M, INFO, MAXLEN
    DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
    DOUBLE PRECISION, DIMENSION(2) :: TIMES1
    DOUBLE PRECISION :: NRMX, RES, RES2, NRMRHS, TMP
    INTEGER :: LDWORK, K
    CHARACTER(1) TRANSA, TRANSB
    EXTERNAL DLANGE
    DOUBLE PRECISION DLANGE
    DOUBLE PRECISION :: TAU
    INTEGER :: MAXIT, SOLVER
    EXTERNAL DGEMM, DLA_GGLYAP_REFINE, DLA_RESIDUAL_GLYAP
    DOUBLE PRECISION :: DLA_RESIDUAL_GLYAP
    EXTERNAL MEPACK_INITIALIZE

    CHARACTER(LEN=1024) :: FILENAME1, FN

    CALL MEPACK_INITIALIZE()

    TRANSA ='T'
    TRANSB ='N'
    MAXIT = 30
    TAU = 0.001

    NARGS = COMMAND_ARGUMENT_COUNT()
    IF ( NARGS .NE. 2 ) THEN
        WRITE(*,'(A)') "Usage: solve_refine_glyap pair1.h5 SOLVER"
        WRITE(*,'(A)') ""
        WRITE(*,'(A)') "Thereby SOLVER is one of the following values :"
        WRITE(*,'(A)') " 0 - OpenMP4 - DAG solver"
        WRITE(*,'(A)') " 1 - Standard Level 3 Solver"
        WRITE(*,'(A)') " 2 - Level-3 2-Stage solver"
        WRITE(*,'(A)') " 3 - Recursive Blocking"
        STOP
    END IF

    CALL GET_COMMAND_ARGUMENT(1, value=FN, length=MAXLEN)
    FILENAME1 = FN(1:LEN_TRIM(FN))
    CALL GET_COMMAND_ARGUMENT(2, value=FN, length=MAXLEN)
    READ(FN, *) SOLVER

    CALL CSC_HDF5_OPEN(H5FILE1, TRIM(FILENAME1), "r")

    ! Read Matrix META
    CALL CSC_HDF5_MATRIX_SIZE_F(H5FILE1, "ASCHUR", M1, M2)

    WRITE(*,'(A)') "# Solve a generalized Lyapunov equation with Iterative Refinement/MEPACK"
    WRITE(*,'(A, I6)') "# Dimension of the Problem - M = ", M1
    M = M1
    WRITE(*,'(A, A)') "# Input filename = ", TRIM(FILENAME1)

    ALLOCATE(A(M,M), AS(M,M), C(M,M), CS(M,M))
    ALLOCATE(Q(M,M), Z(M,M))
    ALLOCATE(X(M, M), Y(M,M), W(M,M))

    SELECT CASE (SOLVER)
        CASE(0)
            WRITE(*,'(A)') "# Solver = OpenMP4 - DAG "
            CALL TGLYAP_FRONTEND_SOLVER_SET(FRONTEND_SOLVER_DAG)
        CASE(1)
            WRITE(*,'(A)') "# Solver = Level-3"
            CALL TGLYAP_FRONTEND_SOLVER_SET(FRONTEND_SOLVER_LEVEL3)
        CASE(2)
            WRITE(*,'(A)') "# Solver = Level-4 2-Stage"
            CALL TGLYAP_FRONTEND_SOLVER_SET(FRONTEND_SOLVER_2STAGE)
        CASE(3)
            WRITE(*,'(A)') "# Solver = Recursive Blocking"
            CALL TGLYAP_FRONTEND_SOLVER_SET(FRONTEND_SOLVER_RECURSIVE)
        CASE DEFAULT
            WRITE(*,'(A)') "# Solver is unknown, please select a proper one"
            STOP
    END SELECT

    ! Workspace Query
    LDWORK = -1
    CALL DLA_GGLYAP_REFINE(TRANSA, "N", M, A, M, C, M, Y, M, W, M,&
        & AS, M, CS, M, Q, M , Z, M, MAXIT, TAU, CONVLOG, TMP, LDWORK, INFO)

    IF ( INFO .EQ. 0 .AND. LDWORK .GE. 0 ) THEN
        WRITE(*,'(A, I10)') '# WORKSPACE = ', LDWORK
    ELSE
        WRITE(*,'(A)') '# Workspace query failed. '
        STOP
    END IF

    ALLOCATE(WORK(LDWORK))

    CALL CSC_HDF5_READ_MATRIX(H5FILE1, "A", M, M, A, M,INFO)
    CALL CSC_HDF5_READ_MATRIX(H5FILE1, "B", M, M, C, M,INFO)
    CALL CSC_HDF5_READ_MATRIX(H5FILE1, "ASCHUR", M, M, AS, M, INFO)
    CALL CSC_HDF5_READ_MATRIX(H5FILE1, "BSCHUR", M, M, CS, M, INFO)
    CALL CSC_HDF5_READ_MATRIX(H5FILE1, "Q", M, M, Q, M, INFO)
    CALL CSC_HDF5_READ_MATRIX(H5FILE1, "Z", M, M, Z, M, INFO)
    CALL CSC_HDF5_READ_VECTOR(H5FILE1, "TIMES", 2, TIMES1, INFO)
    CALL CSC_HDF5_CLOSE(H5FILE1)

    DO K = 1, M
        IF ( CS(K,K) .EQ. 0.0D0) THEN
            WRITE(*,*) "B is singular."
            STOP
        END IF
    END DO

    X = 1.0D0
    Y = 0.0D0

    NRMX = DLANGE("F", M, M, X, M, WORK)

    CALL DGEMM(TRANSA, "N", M, M, M, ONE, A, M, X, M, ZERO, WORK, M)
    CALL DGEMM("N", TRANSB, M, M, M, ONE, WORK, M, C, M, ONE, Y, M)
    CALL DGEMM(TRANSA, "N", M, M, M, ONE, C, M, X, M, ZERO, WORK, M)
    CALL DGEMM("N", TRANSB, M, M, M, ONE, WORK, M, A, M, ONE, Y, M)
    NRMRHS = DLANGE("F", M, M, Y, M, WORK)

    WORK = 0.0D0
    W = Y
    X = Y

    TS =CSC_WTIME()
    CALL DLA_GGLYAP_REFINE(TRANSA, "N", M, A, M, C, M, Y, M, W, M,&
        & AS, M, CS, M, Q, M , Z, M, MAXIT, TAU, CONVLOG, WORK(1), LDWORK, INFO)
    TE = CSC_WTIME()

    TALL = TE - TS + TIMES1(1)

    RES2 = DLA_RESIDUAL_GLYAP(TRANSA, M, A, M, C, M, Y, M, X, M, ONE)

    Y = Y - 1.0D0
    RES = DLANGE("F", M, M, Y, M, WORK) / NRMX

    WRITE(*,'(A, F15.6)') "# Refinement LYAP TIME = ", TE-TS
    WRITE(*,'(A, I3)') "# Number of Iterations = ", MAXIT
    WRITE(*,'(A, F15.6)') "# QZ TIME     = ", TIMES1(1)
    WRITE(*,'(A)') "# M, IT,  Time, Forward, Residual"
    WRITE(*,'(I5," ",I3, " ", D20.13," ",D20.13," ",D20.13)') M, MAXIT, TALL, RES/NRMX, RES2

    WRITE(*,'(A)') "#"
    WRITE(*,'(A)') "# Refinement Log"
    WRITE(*,'(A)') "# IT, Relative Residual"
    DO K = 1, MAXIT
        WRITE(*,'(I3, D22.15)') K, CONVLOG(K)
    END DO

    DEALLOCATE(A, AS, C, CS, X, Y, Q, Z, WORK, W)

END PROGRAM


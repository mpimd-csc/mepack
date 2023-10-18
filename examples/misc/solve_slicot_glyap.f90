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
! Solve a generalized Lyapunov Equation with SLICOT.
!


PROGRAM RUN_SLICOT_GLYAP
    USE CSC_HDF5
    USE CSC_COMMON
    USE CSC_SYSINFO
    IMPLICIT NONE


    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A, C, AS, CS, Y, X, Q, Z
    DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
    INTEGER(CSC_HDF5_T) :: H5FILE1
    DOUBLE PRECISION :: TS, TE, TALL

    INTEGER :: NARGS, M1, M2, M, INFO, MAXLEN
    DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
    DOUBLE PRECISION, DIMENSION(2) :: TIMES1
    DOUBLE PRECISION :: SCALE, NRMX, RES, NRMRHS, RES2
    INTEGER :: LDWORK, K
    CHARACTER(1) TRANSA, TRANSB
    EXTERNAL DLANGE
    DOUBLE PRECISION DLANGE
    EXTERNAL SG03AD
    DOUBLE PRECISION FERR,SEP
    EXTERNAL DGEMM, DLA_GGLYAP_L3, DLA_RESIDUAL_GLYAP
    DOUBLE PRECISION :: DLA_RESIDUAL_GLYAP

    CHARACTER(LEN=1024) :: FILENAME1, FN

    EXTERNAL MEPACK_INITIALIZE

    CALL MEPACK_INITIALIZE()

    NARGS = COMMAND_ARGUMENT_COUNT()
    IF ( NARGS .NE. 1 ) THEN
        WRITE(*,*) "Usage: solve_slicot_glyap input.mat"
        STOP
    END IF


    CALL GET_COMMAND_ARGUMENT(1, value=FN, length=MAXLEN)
    FILENAME1 = FN(1:LEN_TRIM(FN))

    CALL CSC_HDF5_OPEN(H5FILE1, TRIM(FILENAME1), "r")

    ! Read Matrix META
    CALL CSC_HDF5_MATRIX_SIZE_F(H5FILE1, "ASCHUR", M1, M2)

    WRITE(*,'(A)') "# Solve a generalized Lyapunov equation with SLICOT"
    WRITE(*,'(A, I6)') "# Dimension of the Problem - M = ", M1
    M = M1

    ALLOCATE(A(M,M), AS(M,M), C(M,M), CS(M,M))
    ALLOCATE(Q(M,M), Z(M,M))
    ALLOCATE(X(M, M), Y(M,M))

    LDWORK = M*M+128
    ALLOCATE(WORK(LDWORK))
    WRITE(*,'(A, A)') "# Input Filename = ", TRIM(FILENAME1)

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

    TRANSA ='N'
    TRANSB ='T'
    NRMX = DLANGE("F", M, M, X, M, WORK)

    CALL DGEMM(TRANSA, "N", M, M, M, ONE, A, M, X, M, ZERO, WORK, M)
    CALL DGEMM("N", TRANSB, M, M, M, ONE, WORK, M, C, M, ONE, Y, M)
    CALL DGEMM(TRANSA, "N", M, M, M, ONE, C, M, X, M, ZERO, WORK, M)
    CALL DGEMM("N", TRANSB, M, M, M, ONE, WORK, M, A, M, ONE, Y, M)

    X = Y

    NRMRHS = DLANGE("F", M, M, Y, M, WORK)

    WORK = 0.0D0

    TS = CSC_WTIME()
    CALL SG03AD("C", "X", "F", TRANSB, "U", M, AS, M, CS, M, Q, M, Z, M, Y, M, SCALE, SEP, FERR, &
        WORK(1), WORK(2), WORK(3), K, WORK, LDWORK, INFO)
    TE = CSC_WTIME()



    RES2 = DLA_RESIDUAL_GLYAP(TRANSA, M, A, M, C, M, Y, M, X, M, SCALE)
    Y = Y - 1.0D0

    RES = DLANGE("F", M, M, Y, M, WORK)
    TALL = TE - TS + TIMES1(1)

    WRITE(*,'( A, F15.6)') "# SLICOT TIME = ", TE-TS
    WRITE(*,'( A, F15.6)') "# QZ TIME     = ", TIMES1(1)

    WRITE(*,'(A)') "# M, Time, Forward, Residual"
    WRITE(*,'(I5," ",D17.7," ",D20.13," ",D20.13)') M, TALL, RES/NRMX, RES2

    DEALLOCATE(A, AS, C, CS, X, Y, Q, Z, WORK)

END PROGRAM

!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) Martin Koehler, 2017-2023
!


!> \brief Kronecker-product based Kernel for the generalized Sylvester equation
!
!  Definition:
!  ===========
!      DLA_TGSYLV_KRON_KERNEL(TRANSA, TRANSB, SGN, M, N, FA, LDA, FB, LDB, FC, LDC, FD, LDD, X, LDX, SCAL, INFO)
!                CHARACTER(1) TRANSA, TRANSB
!                INTEGER  M,N,LDA,LDB,LDC,LDD,LDX
!                INTEGER  INFO
!                DOUBLE PRECISION FA(LDA,*)
!                DOUBLE PRECISION FB(LDB,*)
!                DOUBLE PRECISION FC(LDC,*)
!                DOUBLE PRECISION FD(LDD,*)
!                DOUBLE PRECISION X(LDX,*)
!                DOUBLE PRECISION SCAL, SGN
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> To solve for X either the reduced generalized Sylvester equation
!>
!>      opA(FA) * X * opB(FB) + opA(FC) * X * opB(FD) = SCAL * Y        (1)
!>
!> or
!>
!>      opA(FA) * X * opB(FB) - opA(FC) * X * opB(FD) = SCAL * Y        (2)
!>
!> where A and C are full matrices and B and D are at most upper
!> Hessenberg. Especially this includes the case where  (A,C) and
!> (D,B) are in generalized Schur form. The matrices
!> A and C or respectively B and D are transposed implicitly inside
!> the algorithm. SCALE is an output scale factor, set to avoid
!> overflow in X.
!>
!> The routine is limited by 0 <= M <= 4 and 0 <= N <= 4. For solving
!> large generalized Sylvester equations please use a different solver.
!>
!> The solution X of (1) or (2) is computed via the Kronecker product representation
!> of (1) or (2). Therefore, this routine works only on small matrices up to dimension
!> 4 x 4.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A and C:
!>          == 'N':  opA(A) = A, opA(C) = C (No transpose for A and C)
!>          == 'T':  opA(A) = A**T, opA(C) = C **T (Transpose A and C)
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to B and D:
!>          == 'N':  opB(B) = B, opB(D) = D (No transpose for B and D)
!>          == 'T':  opB(B) = B**T, opB(D) = D **T (Transpose B and D)
!> \endverbatim
!>
!> \param[in] SGN
!> \verbatim
!>          SGN is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between the two parts of the Sylvester equation.
!>          = 1 :  Solve Equation (1)
!>          == -1:  Solve Equation (2)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  4 >= M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices B and D.  4 >= N >= 0.
!> \endverbatim
!>
!> \param[in] FA
!> \verbatim
!>          FA is DOUBLE PRECISION array, dimension (LDA,M)
!>          The matrix FA must be (quasi-) upper triangular.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array FA.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] FB
!> \verbatim
!>          FB is DOUBLE PRECISION array, dimension (LDB,N)
!>          The matrix FB must be (quasi-) upper triangular. If the matrix D is already
!>          quasi-upper triangular the matrix B must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array FB.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in] FC
!> \verbatim
!>          FC is DOUBLE PRECISION array, dimension (LDC,M)
!>          The matrix FC must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array FC.  LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in] FD
!> \verbatim
!>          FD is DOUBLE PRECISION array, dimension (LDD,N)
!>          The matrix FD must be (quasi-) upper triangular. If the matrix B is already
!>          quasi-upper triangular the matrix FD must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array FD.  LDD >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1) or (2)
!>          as selected by TRANSA, TRANSB, and SGN.
!>          Right hand side Y and the solution X are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] SCAL
!> \verbatim
!>          SCAL is DOUBLE PRECISION
!>          SCAL is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCAL <= 1 holds true.
!> \endverbatim
!>
!> \param[in,out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          On input:
!>          = 1 :  Skip the input data checks
!>          <> 1:  Check input data like normal LAPACK like routines.
!>          On output:
!>          == 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  The equation is not solved correctly. One of the arising inner
!>                system got singular.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup dbltghelp
!
SUBROUTINE DLA_TGSYLV_KRON_KERNEL(TRANSA, TRANSB, SGN, M, N, FA, LDA, FB, LDB, FC, LDC, FD, LDD, X, LDX, SCAL, INFO)

    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Input Parameters
    CHARACTER(1) TRANSA, TRANSB
    INTEGER  M,N,LDA,LDB,LDC,LDD,LDX
    INTEGER  INFO
    DOUBLE PRECISION FA(LDA,*)
    DOUBLE PRECISION FB(LDB,*)
    DOUBLE PRECISION FC(LDC,*)
    DOUBLE PRECISION FD(LDD,*)
    DOUBLE PRECISION X(LDX,*)
    DOUBLE PRECISION SCAL, SGN

    ! External Functions
    EXTERNAL DGETC2X, DGESC2
    EXTERNAL XERBLA
    EXTERNAL LSAME
    LOGICAL LSAME

    INTRINSIC MAX

    ! Internal variables
    INTEGER LDZ, DIMZ, I, J, IINFO
    INTEGER IPIV(16), JPIV(16)
    DOUBLE PRECISION Z(16,16)
    DOUBLE PRECISION RHS(16)
    DOUBLE PRECISION LSCAL
    DOUBLE PRECISION ZERO, ONE
    PARAMETER( ZERO = 0.0d0, ONE = 1.0d0)
    DOUBLE PRECISION A(4,4), B(4,4), C(4,4), D(4,4)

    IF ( INFO .NE. 1 ) THEN
        INFO = 0
        IF ( .NOT. LSAME(TRANSA,'N') .AND. .NOT. LSAME(TRANSA,'T')) THEN
            INFO = -1
        ELSE IF ( .NOT. LSAME(TRANSB,'N') .AND. .NOT. LSAME(TRANSB,'T')) THEN
            INFO = -2
        ELSE IF ( SGN .NE. ONE .AND. SGN .NE. -ONE ) THEN
            INFO = -3
        ELSE IF ( M < 0 .OR. M > 4) THEN
            INFO = -4
        ELSE IF ( N < 0 .OR. N > 4) THEN
            INFO = -5
        ELSE IF ( LDA < MAX(1, M) ) THEN
            INFO = -7
        ELSE IF ( LDB < MAX(1, N) ) THEN
            INFO = -9
        ELSE IF ( LDC < MAX(1, M) ) THEN
            INFO = -11
        ELSE IF ( LDD < MAX(1, N) ) THEN
            INFO = -13
        ELSE IF ( LDX < MAX(1, M) ) THEN
            INFO = -14
        END IF

        IF ( INFO .NE. 0 ) THEN
            CALL XERBLA('DLA_TGSYLV_KRON_KERNEL',-INFO)
            RETURN
        END IF
    ELSE
        INFO = 0
    ENDIF

    ! Quick Return
    IF ( M == 0 .OR. N == 0 ) THEN
        RETURN
    ENDIF

    ! Prepare
    SCAL = ONE;
    LSCAL = ONE;

    LDZ = 16
    Z(1:16,1:16) = ZERO

    ! Create a Local Copy of A,B,C,D to ensure that they are in memory for the Kronecker products.
    IF ( LSAME(TRANSA, 'N' )) THEN
        A(1:M,1:M) = FA(1:M, 1:M)
        C(1:M,1:M) = SGN*FC(1:M, 1:M)
    ELSE
        A(1:M,1:M) = TRANSPOSE(FA(1:M, 1:M))
        C(1:M,1:M) = SGN*TRANSPOSE(FC(1:M, 1:M))
    ENDIF
    B(1:N,1:N) = FB(1:N, 1:N)
    D(1:N,1:N) = FD(1:N, 1:N)


    ! Main
    DIMZ = M * N
    IF ( LSAME (TRANSB, 'N') .EQV. .TRUE. ) THEN
        IF ( M == 1 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  2 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  2 ) = Z(  1,  2 ) + B(  2,  1) * A(  1,  1 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  2 ) = Z(  1,  2 ) + D(  2,  1) * C(  1,  1 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  2 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  3,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  3,  2 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  3,  3 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  2 ) = Z(  1,  2 ) + B(  2,  1) * A(  1,  1 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  2,  3 ) = Z(  2,  3 ) + B(  3,  2) * A(  1,  1 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  2 ) = Z(  1,  2 ) + D(  2,  1) * C(  1,  1 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  2,  3 ) = Z(  2,  3 ) + D(  3,  2) * C(  1,  1 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  2 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  3,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  3,  2 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  3,  3 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  4,  1 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z(  4,  2 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z(  4,  3 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z(  4,  4 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  2 ) = Z(  1,  2 ) + B(  2,  1) * A(  1,  1 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  2,  3 ) = Z(  2,  3 ) + B(  3,  2) * A(  1,  1 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  3,  4 ) = Z(  3,  4 ) + B(  4,  3) * A(  1,  1 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  2 ) = Z(  1,  2 ) + D(  2,  1) * C(  1,  1 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  2,  3 ) = Z(  2,  3 ) + D(  3,  2) * C(  1,  1 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  3,  4 ) = Z(  3,  4 ) + D(  4,  3) * C(  1,  1 )
                END IF
            END IF
        ELSE IF ( M == 2 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  4,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  4,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  3 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  4,  3 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  4,  4 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  3 ) = Z(  1,  3 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  3 ) = Z(  2,  3 ) + B(  2,  1) * A(  2,  1 )
                    Z(  1,  4 ) = Z(  1,  4 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  4 ) = Z(  2,  4 ) + B(  2,  1) * A(  2,  2 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  3 ) = Z(  1,  3 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  3 ) = Z(  2,  3 ) + D(  2,  1) * C(  2,  1 )
                    Z(  1,  4 ) = Z(  1,  4 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  4 ) = Z(  2,  4 ) + D(  2,  1) * C(  2,  2 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  4,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  4,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  3 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  4,  3 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  4,  4 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  5,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  6,  1 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  5,  2 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  6,  2 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  5,  3 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  6,  3 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  5,  4 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  6,  4 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  5,  5 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  5,  6 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  3 ) = Z(  1,  3 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  3 ) = Z(  2,  3 ) + B(  2,  1) * A(  2,  1 )
                    Z(  1,  4 ) = Z(  1,  4 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  4 ) = Z(  2,  4 ) + B(  2,  1) * A(  2,  2 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  3,  2) * A(  1,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + B(  3,  2) * A(  2,  1 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  3,  2) * A(  1,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + B(  3,  2) * A(  2,  2 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  3 ) = Z(  1,  3 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  3 ) = Z(  2,  3 ) + D(  2,  1) * C(  2,  1 )
                    Z(  1,  4 ) = Z(  1,  4 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  4 ) = Z(  2,  4 ) + D(  2,  1) * C(  2,  2 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  3,  2) * C(  1,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + D(  3,  2) * C(  2,  1 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  3,  2) * C(  1,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + D(  3,  2) * C(  2,  2 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  4,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  4,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  3 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  4,  3 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  4,  4 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  5,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  6,  1 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  5,  2 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  6,  2 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  5,  3 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  6,  3 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  5,  4 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  6,  4 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  5,  5 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  5,  6 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z(  7,  1 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z(  8,  1 ) = B(  1,  4) * A(  2,  1 ) + D(  1,  4 ) * C(  2,  1 )
                Z(  7,  2 ) = B(  1,  4) * A(  1,  2 ) + D(  1,  4 ) * C(  1,  2 )
                Z(  8,  2 ) = B(  1,  4) * A(  2,  2 ) + D(  1,  4 ) * C(  2,  2 )
                Z(  7,  3 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z(  8,  3 ) = B(  2,  4) * A(  2,  1 ) + D(  2,  4 ) * C(  2,  1 )
                Z(  7,  4 ) = B(  2,  4) * A(  1,  2 ) + D(  2,  4 ) * C(  1,  2 )
                Z(  8,  4 ) = B(  2,  4) * A(  2,  2 ) + D(  2,  4 ) * C(  2,  2 )
                Z(  7,  5 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z(  8,  5 ) = B(  3,  4) * A(  2,  1 ) + D(  3,  4 ) * C(  2,  1 )
                Z(  7,  6 ) = B(  3,  4) * A(  1,  2 ) + D(  3,  4 ) * C(  1,  2 )
                Z(  8,  6 ) = B(  3,  4) * A(  2,  2 ) + D(  3,  4 ) * C(  2,  2 )
                Z(  7,  7 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )
                Z(  8,  7 ) = B(  4,  4) * A(  2,  1 ) + D(  4,  4 ) * C(  2,  1 )
                Z(  7,  8 ) = B(  4,  4) * A(  1,  2 ) + D(  4,  4 ) * C(  1,  2 )
                Z(  8,  8 ) = B(  4,  4) * A(  2,  2 ) + D(  4,  4 ) * C(  2,  2 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  3 ) = Z(  1,  3 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  3 ) = Z(  2,  3 ) + B(  2,  1) * A(  2,  1 )
                    Z(  1,  4 ) = Z(  1,  4 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  4 ) = Z(  2,  4 ) + B(  2,  1) * A(  2,  2 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  3,  2) * A(  1,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + B(  3,  2) * A(  2,  1 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  3,  2) * A(  1,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + B(  3,  2) * A(  2,  2 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  5,  7 ) = Z(  5,  7 ) + B(  4,  3) * A(  1,  1 )
                    Z(  6,  7 ) = Z(  6,  7 ) + B(  4,  3) * A(  2,  1 )
                    Z(  5,  8 ) = Z(  5,  8 ) + B(  4,  3) * A(  1,  2 )
                    Z(  6,  8 ) = Z(  6,  8 ) + B(  4,  3) * A(  2,  2 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  3 ) = Z(  1,  3 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  3 ) = Z(  2,  3 ) + D(  2,  1) * C(  2,  1 )
                    Z(  1,  4 ) = Z(  1,  4 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  4 ) = Z(  2,  4 ) + D(  2,  1) * C(  2,  2 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  3,  2) * C(  1,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + D(  3,  2) * C(  2,  1 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  3,  2) * C(  1,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + D(  3,  2) * C(  2,  2 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  5,  7 ) = Z(  5,  7 ) + D(  4,  3) * C(  1,  1 )
                    Z(  6,  7 ) = Z(  6,  7 ) + D(  4,  3) * C(  2,  1 )
                    Z(  5,  8 ) = Z(  5,  8 ) + D(  4,  3) * C(  1,  2 )
                    Z(  6,  8 ) = Z(  6,  8 ) + D(  4,  3) * C(  2,  2 )
                END IF
            END IF
        ELSE IF ( M == 3 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  5,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  6,  1 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  4,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  5,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  6,  2 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  4,  3 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  5,  3 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  6,  3 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  4 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  5,  4 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  6,  4 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  5,  5 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  6,  5 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  5,  6 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  6,  6 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  4 ) = Z(  1,  4 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  4 ) = Z(  2,  4 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  4 ) = Z(  3,  4 ) + B(  2,  1) * A(  3,  1 )
                    Z(  1,  5 ) = Z(  1,  5 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  5 ) = Z(  2,  5 ) + B(  2,  1) * A(  2,  2 )
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  2,  1) * A(  3,  2 )
                    Z(  1,  6 ) = Z(  1,  6 ) + B(  2,  1) * A(  1,  3 )
                    Z(  2,  6 ) = Z(  2,  6 ) + B(  2,  1) * A(  2,  3 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  2,  1) * A(  3,  3 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  4 ) = Z(  1,  4 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  4 ) = Z(  2,  4 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  4 ) = Z(  3,  4 ) + D(  2,  1) * C(  3,  1 )
                    Z(  1,  5 ) = Z(  1,  5 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  5 ) = Z(  2,  5 ) + D(  2,  1) * C(  2,  2 )
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  2,  1) * C(  3,  2 )
                    Z(  1,  6 ) = Z(  1,  6 ) + D(  2,  1) * C(  1,  3 )
                    Z(  2,  6 ) = Z(  2,  6 ) + D(  2,  1) * C(  2,  3 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  2,  1) * C(  3,  3 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  5,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  6,  1 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  4,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  5,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  6,  2 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  4,  3 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  5,  3 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  6,  3 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  4 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  5,  4 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  6,  4 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  5,  5 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  6,  5 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  5,  6 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  6,  6 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  7,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  8,  1 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  9,  1 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z(  7,  2 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  8,  2 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  9,  2 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z(  7,  3 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z(  8,  3 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z(  9,  3 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z(  7,  4 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  8,  4 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  9,  4 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z(  7,  5 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  8,  5 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  9,  5 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z(  7,  6 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z(  8,  6 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z(  9,  6 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z(  7,  7 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  8,  7 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  9,  7 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z(  7,  8 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  8,  8 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z(  9,  8 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z(  7,  9 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z(  8,  9 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z(  9,  9 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  4 ) = Z(  1,  4 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  4 ) = Z(  2,  4 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  4 ) = Z(  3,  4 ) + B(  2,  1) * A(  3,  1 )
                    Z(  1,  5 ) = Z(  1,  5 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  5 ) = Z(  2,  5 ) + B(  2,  1) * A(  2,  2 )
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  2,  1) * A(  3,  2 )
                    Z(  1,  6 ) = Z(  1,  6 ) + B(  2,  1) * A(  1,  3 )
                    Z(  2,  6 ) = Z(  2,  6 ) + B(  2,  1) * A(  2,  3 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  2,  1) * A(  3,  3 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  4,  7 ) = Z(  4,  7 ) + B(  3,  2) * A(  1,  1 )
                    Z(  5,  7 ) = Z(  5,  7 ) + B(  3,  2) * A(  2,  1 )
                    Z(  6,  7 ) = Z(  6,  7 ) + B(  3,  2) * A(  3,  1 )
                    Z(  4,  8 ) = Z(  4,  8 ) + B(  3,  2) * A(  1,  2 )
                    Z(  5,  8 ) = Z(  5,  8 ) + B(  3,  2) * A(  2,  2 )
                    Z(  6,  8 ) = Z(  6,  8 ) + B(  3,  2) * A(  3,  2 )
                    Z(  4,  9 ) = Z(  4,  9 ) + B(  3,  2) * A(  1,  3 )
                    Z(  5,  9 ) = Z(  5,  9 ) + B(  3,  2) * A(  2,  3 )
                    Z(  6,  9 ) = Z(  6,  9 ) + B(  3,  2) * A(  3,  3 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  4 ) = Z(  1,  4 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  4 ) = Z(  2,  4 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  4 ) = Z(  3,  4 ) + D(  2,  1) * C(  3,  1 )
                    Z(  1,  5 ) = Z(  1,  5 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  5 ) = Z(  2,  5 ) + D(  2,  1) * C(  2,  2 )
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  2,  1) * C(  3,  2 )
                    Z(  1,  6 ) = Z(  1,  6 ) + D(  2,  1) * C(  1,  3 )
                    Z(  2,  6 ) = Z(  2,  6 ) + D(  2,  1) * C(  2,  3 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  2,  1) * C(  3,  3 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  4,  7 ) = Z(  4,  7 ) + D(  3,  2) * C(  1,  1 )
                    Z(  5,  7 ) = Z(  5,  7 ) + D(  3,  2) * C(  2,  1 )
                    Z(  6,  7 ) = Z(  6,  7 ) + D(  3,  2) * C(  3,  1 )
                    Z(  4,  8 ) = Z(  4,  8 ) + D(  3,  2) * C(  1,  2 )
                    Z(  5,  8 ) = Z(  5,  8 ) + D(  3,  2) * C(  2,  2 )
                    Z(  6,  8 ) = Z(  6,  8 ) + D(  3,  2) * C(  3,  2 )
                    Z(  4,  9 ) = Z(  4,  9 ) + D(  3,  2) * C(  1,  3 )
                    Z(  5,  9 ) = Z(  5,  9 ) + D(  3,  2) * C(  2,  3 )
                    Z(  6,  9 ) = Z(  6,  9 ) + D(  3,  2) * C(  3,  3 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  5,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  6,  1 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  4,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  5,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  6,  2 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  4,  3 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  5,  3 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  6,  3 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  4 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  5,  4 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  6,  4 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  5,  5 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  6,  5 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  5,  6 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  6,  6 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  7,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  8,  1 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  9,  1 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z(  7,  2 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  8,  2 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  9,  2 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z(  7,  3 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z(  8,  3 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z(  9,  3 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z(  7,  4 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  8,  4 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  9,  4 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z(  7,  5 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  8,  5 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  9,  5 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z(  7,  6 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z(  8,  6 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z(  9,  6 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z(  7,  7 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  8,  7 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  9,  7 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z(  7,  8 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  8,  8 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z(  9,  8 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z(  7,  9 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z(  8,  9 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z(  9,  9 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )
                Z( 10,  1 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z( 11,  1 ) = B(  1,  4) * A(  2,  1 ) + D(  1,  4 ) * C(  2,  1 )
                Z( 12,  1 ) = B(  1,  4) * A(  3,  1 ) + D(  1,  4 ) * C(  3,  1 )
                Z( 10,  2 ) = B(  1,  4) * A(  1,  2 ) + D(  1,  4 ) * C(  1,  2 )
                Z( 11,  2 ) = B(  1,  4) * A(  2,  2 ) + D(  1,  4 ) * C(  2,  2 )
                Z( 12,  2 ) = B(  1,  4) * A(  3,  2 ) + D(  1,  4 ) * C(  3,  2 )
                Z( 10,  3 ) = B(  1,  4) * A(  1,  3 ) + D(  1,  4 ) * C(  1,  3 )
                Z( 11,  3 ) = B(  1,  4) * A(  2,  3 ) + D(  1,  4 ) * C(  2,  3 )
                Z( 12,  3 ) = B(  1,  4) * A(  3,  3 ) + D(  1,  4 ) * C(  3,  3 )
                Z( 10,  4 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z( 11,  4 ) = B(  2,  4) * A(  2,  1 ) + D(  2,  4 ) * C(  2,  1 )
                Z( 12,  4 ) = B(  2,  4) * A(  3,  1 ) + D(  2,  4 ) * C(  3,  1 )
                Z( 10,  5 ) = B(  2,  4) * A(  1,  2 ) + D(  2,  4 ) * C(  1,  2 )
                Z( 11,  5 ) = B(  2,  4) * A(  2,  2 ) + D(  2,  4 ) * C(  2,  2 )
                Z( 12,  5 ) = B(  2,  4) * A(  3,  2 ) + D(  2,  4 ) * C(  3,  2 )
                Z( 10,  6 ) = B(  2,  4) * A(  1,  3 ) + D(  2,  4 ) * C(  1,  3 )
                Z( 11,  6 ) = B(  2,  4) * A(  2,  3 ) + D(  2,  4 ) * C(  2,  3 )
                Z( 12,  6 ) = B(  2,  4) * A(  3,  3 ) + D(  2,  4 ) * C(  3,  3 )
                Z( 10,  7 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z( 11,  7 ) = B(  3,  4) * A(  2,  1 ) + D(  3,  4 ) * C(  2,  1 )
                Z( 12,  7 ) = B(  3,  4) * A(  3,  1 ) + D(  3,  4 ) * C(  3,  1 )
                Z( 10,  8 ) = B(  3,  4) * A(  1,  2 ) + D(  3,  4 ) * C(  1,  2 )
                Z( 11,  8 ) = B(  3,  4) * A(  2,  2 ) + D(  3,  4 ) * C(  2,  2 )
                Z( 12,  8 ) = B(  3,  4) * A(  3,  2 ) + D(  3,  4 ) * C(  3,  2 )
                Z( 10,  9 ) = B(  3,  4) * A(  1,  3 ) + D(  3,  4 ) * C(  1,  3 )
                Z( 11,  9 ) = B(  3,  4) * A(  2,  3 ) + D(  3,  4 ) * C(  2,  3 )
                Z( 12,  9 ) = B(  3,  4) * A(  3,  3 ) + D(  3,  4 ) * C(  3,  3 )
                Z( 10, 10 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )
                Z( 11, 10 ) = B(  4,  4) * A(  2,  1 ) + D(  4,  4 ) * C(  2,  1 )
                Z( 12, 10 ) = B(  4,  4) * A(  3,  1 ) + D(  4,  4 ) * C(  3,  1 )
                Z( 10, 11 ) = B(  4,  4) * A(  1,  2 ) + D(  4,  4 ) * C(  1,  2 )
                Z( 11, 11 ) = B(  4,  4) * A(  2,  2 ) + D(  4,  4 ) * C(  2,  2 )
                Z( 12, 11 ) = B(  4,  4) * A(  3,  2 ) + D(  4,  4 ) * C(  3,  2 )
                Z( 10, 12 ) = B(  4,  4) * A(  1,  3 ) + D(  4,  4 ) * C(  1,  3 )
                Z( 11, 12 ) = B(  4,  4) * A(  2,  3 ) + D(  4,  4 ) * C(  2,  3 )
                Z( 12, 12 ) = B(  4,  4) * A(  3,  3 ) + D(  4,  4 ) * C(  3,  3 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  4 ) = Z(  1,  4 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  4 ) = Z(  2,  4 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  4 ) = Z(  3,  4 ) + B(  2,  1) * A(  3,  1 )
                    Z(  1,  5 ) = Z(  1,  5 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  5 ) = Z(  2,  5 ) + B(  2,  1) * A(  2,  2 )
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  2,  1) * A(  3,  2 )
                    Z(  1,  6 ) = Z(  1,  6 ) + B(  2,  1) * A(  1,  3 )
                    Z(  2,  6 ) = Z(  2,  6 ) + B(  2,  1) * A(  2,  3 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  2,  1) * A(  3,  3 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  4,  7 ) = Z(  4,  7 ) + B(  3,  2) * A(  1,  1 )
                    Z(  5,  7 ) = Z(  5,  7 ) + B(  3,  2) * A(  2,  1 )
                    Z(  6,  7 ) = Z(  6,  7 ) + B(  3,  2) * A(  3,  1 )
                    Z(  4,  8 ) = Z(  4,  8 ) + B(  3,  2) * A(  1,  2 )
                    Z(  5,  8 ) = Z(  5,  8 ) + B(  3,  2) * A(  2,  2 )
                    Z(  6,  8 ) = Z(  6,  8 ) + B(  3,  2) * A(  3,  2 )
                    Z(  4,  9 ) = Z(  4,  9 ) + B(  3,  2) * A(  1,  3 )
                    Z(  5,  9 ) = Z(  5,  9 ) + B(  3,  2) * A(  2,  3 )
                    Z(  6,  9 ) = Z(  6,  9 ) + B(  3,  2) * A(  3,  3 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  7, 10 ) = Z(  7, 10 ) + B(  4,  3) * A(  1,  1 )
                    Z(  8, 10 ) = Z(  8, 10 ) + B(  4,  3) * A(  2,  1 )
                    Z(  9, 10 ) = Z(  9, 10 ) + B(  4,  3) * A(  3,  1 )
                    Z(  7, 11 ) = Z(  7, 11 ) + B(  4,  3) * A(  1,  2 )
                    Z(  8, 11 ) = Z(  8, 11 ) + B(  4,  3) * A(  2,  2 )
                    Z(  9, 11 ) = Z(  9, 11 ) + B(  4,  3) * A(  3,  2 )
                    Z(  7, 12 ) = Z(  7, 12 ) + B(  4,  3) * A(  1,  3 )
                    Z(  8, 12 ) = Z(  8, 12 ) + B(  4,  3) * A(  2,  3 )
                    Z(  9, 12 ) = Z(  9, 12 ) + B(  4,  3) * A(  3,  3 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  4 ) = Z(  1,  4 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  4 ) = Z(  2,  4 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  4 ) = Z(  3,  4 ) + D(  2,  1) * C(  3,  1 )
                    Z(  1,  5 ) = Z(  1,  5 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  5 ) = Z(  2,  5 ) + D(  2,  1) * C(  2,  2 )
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  2,  1) * C(  3,  2 )
                    Z(  1,  6 ) = Z(  1,  6 ) + D(  2,  1) * C(  1,  3 )
                    Z(  2,  6 ) = Z(  2,  6 ) + D(  2,  1) * C(  2,  3 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  2,  1) * C(  3,  3 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  4,  7 ) = Z(  4,  7 ) + D(  3,  2) * C(  1,  1 )
                    Z(  5,  7 ) = Z(  5,  7 ) + D(  3,  2) * C(  2,  1 )
                    Z(  6,  7 ) = Z(  6,  7 ) + D(  3,  2) * C(  3,  1 )
                    Z(  4,  8 ) = Z(  4,  8 ) + D(  3,  2) * C(  1,  2 )
                    Z(  5,  8 ) = Z(  5,  8 ) + D(  3,  2) * C(  2,  2 )
                    Z(  6,  8 ) = Z(  6,  8 ) + D(  3,  2) * C(  3,  2 )
                    Z(  4,  9 ) = Z(  4,  9 ) + D(  3,  2) * C(  1,  3 )
                    Z(  5,  9 ) = Z(  5,  9 ) + D(  3,  2) * C(  2,  3 )
                    Z(  6,  9 ) = Z(  6,  9 ) + D(  3,  2) * C(  3,  3 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  7, 10 ) = Z(  7, 10 ) + D(  4,  3) * C(  1,  1 )
                    Z(  8, 10 ) = Z(  8, 10 ) + D(  4,  3) * C(  2,  1 )
                    Z(  9, 10 ) = Z(  9, 10 ) + D(  4,  3) * C(  3,  1 )
                    Z(  7, 11 ) = Z(  7, 11 ) + D(  4,  3) * C(  1,  2 )
                    Z(  8, 11 ) = Z(  8, 11 ) + D(  4,  3) * C(  2,  2 )
                    Z(  9, 11 ) = Z(  9, 11 ) + D(  4,  3) * C(  3,  2 )
                    Z(  7, 12 ) = Z(  7, 12 ) + D(  4,  3) * C(  1,  3 )
                    Z(  8, 12 ) = Z(  8, 12 ) + D(  4,  3) * C(  2,  3 )
                    Z(  9, 12 ) = Z(  9, 12 ) + D(  4,  3) * C(  3,  3 )
                END IF
            END IF
        ELSE IF ( M == 4 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )
                Z(  5,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  6,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  7,  1 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  8,  1 ) = B(  1,  2) * A(  4,  1 ) + D(  1,  2 ) * C(  4,  1 )
                Z(  5,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  6,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  7,  2 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  8,  2 ) = B(  1,  2) * A(  4,  2 ) + D(  1,  2 ) * C(  4,  2 )
                Z(  5,  3 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  6,  3 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  7,  3 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  8,  3 ) = B(  1,  2) * A(  4,  3 ) + D(  1,  2 ) * C(  4,  3 )
                Z(  5,  4 ) = B(  1,  2) * A(  1,  4 ) + D(  1,  2 ) * C(  1,  4 )
                Z(  6,  4 ) = B(  1,  2) * A(  2,  4 ) + D(  1,  2 ) * C(  2,  4 )
                Z(  7,  4 ) = B(  1,  2) * A(  3,  4 ) + D(  1,  2 ) * C(  3,  4 )
                Z(  8,  4 ) = B(  1,  2) * A(  4,  4 ) + D(  1,  2 ) * C(  4,  4 )
                Z(  5,  5 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  7,  5 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  8,  5 ) = B(  2,  2) * A(  4,  1 ) + D(  2,  2 ) * C(  4,  1 )
                Z(  5,  6 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  7,  6 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  8,  6 ) = B(  2,  2) * A(  4,  2 ) + D(  2,  2 ) * C(  4,  2 )
                Z(  5,  7 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  6,  7 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  7,  7 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  8,  7 ) = B(  2,  2) * A(  4,  3 ) + D(  2,  2 ) * C(  4,  3 )
                Z(  5,  8 ) = B(  2,  2) * A(  1,  4 ) + D(  2,  2 ) * C(  1,  4 )
                Z(  6,  8 ) = B(  2,  2) * A(  2,  4 ) + D(  2,  2 ) * C(  2,  4 )
                Z(  7,  8 ) = B(  2,  2) * A(  3,  4 ) + D(  2,  2 ) * C(  3,  4 )
                Z(  8,  8 ) = B(  2,  2) * A(  4,  4 ) + D(  2,  2 ) * C(  4,  4 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  5 ) = Z(  1,  5 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  5 ) = Z(  2,  5 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  2,  1) * A(  3,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + B(  2,  1) * A(  4,  1 )
                    Z(  1,  6 ) = Z(  1,  6 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  6 ) = Z(  2,  6 ) + B(  2,  1) * A(  2,  2 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  2,  1) * A(  3,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + B(  2,  1) * A(  4,  2 )
                    Z(  1,  7 ) = Z(  1,  7 ) + B(  2,  1) * A(  1,  3 )
                    Z(  2,  7 ) = Z(  2,  7 ) + B(  2,  1) * A(  2,  3 )
                    Z(  3,  7 ) = Z(  3,  7 ) + B(  2,  1) * A(  3,  3 )
                    Z(  4,  7 ) = Z(  4,  7 ) + B(  2,  1) * A(  4,  3 )
                    Z(  1,  8 ) = Z(  1,  8 ) + B(  2,  1) * A(  1,  4 )
                    Z(  2,  8 ) = Z(  2,  8 ) + B(  2,  1) * A(  2,  4 )
                    Z(  3,  8 ) = Z(  3,  8 ) + B(  2,  1) * A(  3,  4 )
                    Z(  4,  8 ) = Z(  4,  8 ) + B(  2,  1) * A(  4,  4 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  5 ) = Z(  1,  5 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  5 ) = Z(  2,  5 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  2,  1) * C(  3,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + D(  2,  1) * C(  4,  1 )
                    Z(  1,  6 ) = Z(  1,  6 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  6 ) = Z(  2,  6 ) + D(  2,  1) * C(  2,  2 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  2,  1) * C(  3,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + D(  2,  1) * C(  4,  2 )
                    Z(  1,  7 ) = Z(  1,  7 ) + D(  2,  1) * C(  1,  3 )
                    Z(  2,  7 ) = Z(  2,  7 ) + D(  2,  1) * C(  2,  3 )
                    Z(  3,  7 ) = Z(  3,  7 ) + D(  2,  1) * C(  3,  3 )
                    Z(  4,  7 ) = Z(  4,  7 ) + D(  2,  1) * C(  4,  3 )
                    Z(  1,  8 ) = Z(  1,  8 ) + D(  2,  1) * C(  1,  4 )
                    Z(  2,  8 ) = Z(  2,  8 ) + D(  2,  1) * C(  2,  4 )
                    Z(  3,  8 ) = Z(  3,  8 ) + D(  2,  1) * C(  3,  4 )
                    Z(  4,  8 ) = Z(  4,  8 ) + D(  2,  1) * C(  4,  4 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )
                Z(  5,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  6,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  7,  1 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  8,  1 ) = B(  1,  2) * A(  4,  1 ) + D(  1,  2 ) * C(  4,  1 )
                Z(  5,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  6,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  7,  2 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  8,  2 ) = B(  1,  2) * A(  4,  2 ) + D(  1,  2 ) * C(  4,  2 )
                Z(  5,  3 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  6,  3 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  7,  3 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  8,  3 ) = B(  1,  2) * A(  4,  3 ) + D(  1,  2 ) * C(  4,  3 )
                Z(  5,  4 ) = B(  1,  2) * A(  1,  4 ) + D(  1,  2 ) * C(  1,  4 )
                Z(  6,  4 ) = B(  1,  2) * A(  2,  4 ) + D(  1,  2 ) * C(  2,  4 )
                Z(  7,  4 ) = B(  1,  2) * A(  3,  4 ) + D(  1,  2 ) * C(  3,  4 )
                Z(  8,  4 ) = B(  1,  2) * A(  4,  4 ) + D(  1,  2 ) * C(  4,  4 )
                Z(  5,  5 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  7,  5 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  8,  5 ) = B(  2,  2) * A(  4,  1 ) + D(  2,  2 ) * C(  4,  1 )
                Z(  5,  6 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  7,  6 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  8,  6 ) = B(  2,  2) * A(  4,  2 ) + D(  2,  2 ) * C(  4,  2 )
                Z(  5,  7 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  6,  7 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  7,  7 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  8,  7 ) = B(  2,  2) * A(  4,  3 ) + D(  2,  2 ) * C(  4,  3 )
                Z(  5,  8 ) = B(  2,  2) * A(  1,  4 ) + D(  2,  2 ) * C(  1,  4 )
                Z(  6,  8 ) = B(  2,  2) * A(  2,  4 ) + D(  2,  2 ) * C(  2,  4 )
                Z(  7,  8 ) = B(  2,  2) * A(  3,  4 ) + D(  2,  2 ) * C(  3,  4 )
                Z(  8,  8 ) = B(  2,  2) * A(  4,  4 ) + D(  2,  2 ) * C(  4,  4 )
                Z(  9,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z( 10,  1 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z( 11,  1 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z( 12,  1 ) = B(  1,  3) * A(  4,  1 ) + D(  1,  3 ) * C(  4,  1 )
                Z(  9,  2 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z( 10,  2 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z( 11,  2 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z( 12,  2 ) = B(  1,  3) * A(  4,  2 ) + D(  1,  3 ) * C(  4,  2 )
                Z(  9,  3 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z( 10,  3 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z( 11,  3 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z( 12,  3 ) = B(  1,  3) * A(  4,  3 ) + D(  1,  3 ) * C(  4,  3 )
                Z(  9,  4 ) = B(  1,  3) * A(  1,  4 ) + D(  1,  3 ) * C(  1,  4 )
                Z( 10,  4 ) = B(  1,  3) * A(  2,  4 ) + D(  1,  3 ) * C(  2,  4 )
                Z( 11,  4 ) = B(  1,  3) * A(  3,  4 ) + D(  1,  3 ) * C(  3,  4 )
                Z( 12,  4 ) = B(  1,  3) * A(  4,  4 ) + D(  1,  3 ) * C(  4,  4 )
                Z(  9,  5 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z( 10,  5 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z( 11,  5 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z( 12,  5 ) = B(  2,  3) * A(  4,  1 ) + D(  2,  3 ) * C(  4,  1 )
                Z(  9,  6 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z( 10,  6 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z( 11,  6 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z( 12,  6 ) = B(  2,  3) * A(  4,  2 ) + D(  2,  3 ) * C(  4,  2 )
                Z(  9,  7 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z( 10,  7 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z( 11,  7 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z( 12,  7 ) = B(  2,  3) * A(  4,  3 ) + D(  2,  3 ) * C(  4,  3 )
                Z(  9,  8 ) = B(  2,  3) * A(  1,  4 ) + D(  2,  3 ) * C(  1,  4 )
                Z( 10,  8 ) = B(  2,  3) * A(  2,  4 ) + D(  2,  3 ) * C(  2,  4 )
                Z( 11,  8 ) = B(  2,  3) * A(  3,  4 ) + D(  2,  3 ) * C(  3,  4 )
                Z( 12,  8 ) = B(  2,  3) * A(  4,  4 ) + D(  2,  3 ) * C(  4,  4 )
                Z(  9,  9 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z( 10,  9 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z( 11,  9 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z( 12,  9 ) = B(  3,  3) * A(  4,  1 ) + D(  3,  3 ) * C(  4,  1 )
                Z(  9, 10 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z( 10, 10 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z( 11, 10 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z( 12, 10 ) = B(  3,  3) * A(  4,  2 ) + D(  3,  3 ) * C(  4,  2 )
                Z(  9, 11 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z( 10, 11 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z( 11, 11 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )
                Z( 12, 11 ) = B(  3,  3) * A(  4,  3 ) + D(  3,  3 ) * C(  4,  3 )
                Z(  9, 12 ) = B(  3,  3) * A(  1,  4 ) + D(  3,  3 ) * C(  1,  4 )
                Z( 10, 12 ) = B(  3,  3) * A(  2,  4 ) + D(  3,  3 ) * C(  2,  4 )
                Z( 11, 12 ) = B(  3,  3) * A(  3,  4 ) + D(  3,  3 ) * C(  3,  4 )
                Z( 12, 12 ) = B(  3,  3) * A(  4,  4 ) + D(  3,  3 ) * C(  4,  4 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  5 ) = Z(  1,  5 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  5 ) = Z(  2,  5 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  2,  1) * A(  3,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + B(  2,  1) * A(  4,  1 )
                    Z(  1,  6 ) = Z(  1,  6 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  6 ) = Z(  2,  6 ) + B(  2,  1) * A(  2,  2 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  2,  1) * A(  3,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + B(  2,  1) * A(  4,  2 )
                    Z(  1,  7 ) = Z(  1,  7 ) + B(  2,  1) * A(  1,  3 )
                    Z(  2,  7 ) = Z(  2,  7 ) + B(  2,  1) * A(  2,  3 )
                    Z(  3,  7 ) = Z(  3,  7 ) + B(  2,  1) * A(  3,  3 )
                    Z(  4,  7 ) = Z(  4,  7 ) + B(  2,  1) * A(  4,  3 )
                    Z(  1,  8 ) = Z(  1,  8 ) + B(  2,  1) * A(  1,  4 )
                    Z(  2,  8 ) = Z(  2,  8 ) + B(  2,  1) * A(  2,  4 )
                    Z(  3,  8 ) = Z(  3,  8 ) + B(  2,  1) * A(  3,  4 )
                    Z(  4,  8 ) = Z(  4,  8 ) + B(  2,  1) * A(  4,  4 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  9 ) = Z(  5,  9 ) + B(  3,  2) * A(  1,  1 )
                    Z(  6,  9 ) = Z(  6,  9 ) + B(  3,  2) * A(  2,  1 )
                    Z(  7,  9 ) = Z(  7,  9 ) + B(  3,  2) * A(  3,  1 )
                    Z(  8,  9 ) = Z(  8,  9 ) + B(  3,  2) * A(  4,  1 )
                    Z(  5, 10 ) = Z(  5, 10 ) + B(  3,  2) * A(  1,  2 )
                    Z(  6, 10 ) = Z(  6, 10 ) + B(  3,  2) * A(  2,  2 )
                    Z(  7, 10 ) = Z(  7, 10 ) + B(  3,  2) * A(  3,  2 )
                    Z(  8, 10 ) = Z(  8, 10 ) + B(  3,  2) * A(  4,  2 )
                    Z(  5, 11 ) = Z(  5, 11 ) + B(  3,  2) * A(  1,  3 )
                    Z(  6, 11 ) = Z(  6, 11 ) + B(  3,  2) * A(  2,  3 )
                    Z(  7, 11 ) = Z(  7, 11 ) + B(  3,  2) * A(  3,  3 )
                    Z(  8, 11 ) = Z(  8, 11 ) + B(  3,  2) * A(  4,  3 )
                    Z(  5, 12 ) = Z(  5, 12 ) + B(  3,  2) * A(  1,  4 )
                    Z(  6, 12 ) = Z(  6, 12 ) + B(  3,  2) * A(  2,  4 )
                    Z(  7, 12 ) = Z(  7, 12 ) + B(  3,  2) * A(  3,  4 )
                    Z(  8, 12 ) = Z(  8, 12 ) + B(  3,  2) * A(  4,  4 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  5 ) = Z(  1,  5 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  5 ) = Z(  2,  5 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  2,  1) * C(  3,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + D(  2,  1) * C(  4,  1 )
                    Z(  1,  6 ) = Z(  1,  6 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  6 ) = Z(  2,  6 ) + D(  2,  1) * C(  2,  2 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  2,  1) * C(  3,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + D(  2,  1) * C(  4,  2 )
                    Z(  1,  7 ) = Z(  1,  7 ) + D(  2,  1) * C(  1,  3 )
                    Z(  2,  7 ) = Z(  2,  7 ) + D(  2,  1) * C(  2,  3 )
                    Z(  3,  7 ) = Z(  3,  7 ) + D(  2,  1) * C(  3,  3 )
                    Z(  4,  7 ) = Z(  4,  7 ) + D(  2,  1) * C(  4,  3 )
                    Z(  1,  8 ) = Z(  1,  8 ) + D(  2,  1) * C(  1,  4 )
                    Z(  2,  8 ) = Z(  2,  8 ) + D(  2,  1) * C(  2,  4 )
                    Z(  3,  8 ) = Z(  3,  8 ) + D(  2,  1) * C(  3,  4 )
                    Z(  4,  8 ) = Z(  4,  8 ) + D(  2,  1) * C(  4,  4 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  9 ) = Z(  5,  9 ) + D(  3,  2) * C(  1,  1 )
                    Z(  6,  9 ) = Z(  6,  9 ) + D(  3,  2) * C(  2,  1 )
                    Z(  7,  9 ) = Z(  7,  9 ) + D(  3,  2) * C(  3,  1 )
                    Z(  8,  9 ) = Z(  8,  9 ) + D(  3,  2) * C(  4,  1 )
                    Z(  5, 10 ) = Z(  5, 10 ) + D(  3,  2) * C(  1,  2 )
                    Z(  6, 10 ) = Z(  6, 10 ) + D(  3,  2) * C(  2,  2 )
                    Z(  7, 10 ) = Z(  7, 10 ) + D(  3,  2) * C(  3,  2 )
                    Z(  8, 10 ) = Z(  8, 10 ) + D(  3,  2) * C(  4,  2 )
                    Z(  5, 11 ) = Z(  5, 11 ) + D(  3,  2) * C(  1,  3 )
                    Z(  6, 11 ) = Z(  6, 11 ) + D(  3,  2) * C(  2,  3 )
                    Z(  7, 11 ) = Z(  7, 11 ) + D(  3,  2) * C(  3,  3 )
                    Z(  8, 11 ) = Z(  8, 11 ) + D(  3,  2) * C(  4,  3 )
                    Z(  5, 12 ) = Z(  5, 12 ) + D(  3,  2) * C(  1,  4 )
                    Z(  6, 12 ) = Z(  6, 12 ) + D(  3,  2) * C(  2,  4 )
                    Z(  7, 12 ) = Z(  7, 12 ) + D(  3,  2) * C(  3,  4 )
                    Z(  8, 12 ) = Z(  8, 12 ) + D(  3,  2) * C(  4,  4 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )
                Z(  5,  1 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  6,  1 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  7,  1 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  8,  1 ) = B(  1,  2) * A(  4,  1 ) + D(  1,  2 ) * C(  4,  1 )
                Z(  5,  2 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  6,  2 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  7,  2 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  8,  2 ) = B(  1,  2) * A(  4,  2 ) + D(  1,  2 ) * C(  4,  2 )
                Z(  5,  3 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  6,  3 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  7,  3 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  8,  3 ) = B(  1,  2) * A(  4,  3 ) + D(  1,  2 ) * C(  4,  3 )
                Z(  5,  4 ) = B(  1,  2) * A(  1,  4 ) + D(  1,  2 ) * C(  1,  4 )
                Z(  6,  4 ) = B(  1,  2) * A(  2,  4 ) + D(  1,  2 ) * C(  2,  4 )
                Z(  7,  4 ) = B(  1,  2) * A(  3,  4 ) + D(  1,  2 ) * C(  3,  4 )
                Z(  8,  4 ) = B(  1,  2) * A(  4,  4 ) + D(  1,  2 ) * C(  4,  4 )
                Z(  5,  5 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  7,  5 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  8,  5 ) = B(  2,  2) * A(  4,  1 ) + D(  2,  2 ) * C(  4,  1 )
                Z(  5,  6 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  7,  6 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  8,  6 ) = B(  2,  2) * A(  4,  2 ) + D(  2,  2 ) * C(  4,  2 )
                Z(  5,  7 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  6,  7 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  7,  7 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  8,  7 ) = B(  2,  2) * A(  4,  3 ) + D(  2,  2 ) * C(  4,  3 )
                Z(  5,  8 ) = B(  2,  2) * A(  1,  4 ) + D(  2,  2 ) * C(  1,  4 )
                Z(  6,  8 ) = B(  2,  2) * A(  2,  4 ) + D(  2,  2 ) * C(  2,  4 )
                Z(  7,  8 ) = B(  2,  2) * A(  3,  4 ) + D(  2,  2 ) * C(  3,  4 )
                Z(  8,  8 ) = B(  2,  2) * A(  4,  4 ) + D(  2,  2 ) * C(  4,  4 )
                Z(  9,  1 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z( 10,  1 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z( 11,  1 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z( 12,  1 ) = B(  1,  3) * A(  4,  1 ) + D(  1,  3 ) * C(  4,  1 )
                Z(  9,  2 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z( 10,  2 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z( 11,  2 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z( 12,  2 ) = B(  1,  3) * A(  4,  2 ) + D(  1,  3 ) * C(  4,  2 )
                Z(  9,  3 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z( 10,  3 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z( 11,  3 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z( 12,  3 ) = B(  1,  3) * A(  4,  3 ) + D(  1,  3 ) * C(  4,  3 )
                Z(  9,  4 ) = B(  1,  3) * A(  1,  4 ) + D(  1,  3 ) * C(  1,  4 )
                Z( 10,  4 ) = B(  1,  3) * A(  2,  4 ) + D(  1,  3 ) * C(  2,  4 )
                Z( 11,  4 ) = B(  1,  3) * A(  3,  4 ) + D(  1,  3 ) * C(  3,  4 )
                Z( 12,  4 ) = B(  1,  3) * A(  4,  4 ) + D(  1,  3 ) * C(  4,  4 )
                Z(  9,  5 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z( 10,  5 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z( 11,  5 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z( 12,  5 ) = B(  2,  3) * A(  4,  1 ) + D(  2,  3 ) * C(  4,  1 )
                Z(  9,  6 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z( 10,  6 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z( 11,  6 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z( 12,  6 ) = B(  2,  3) * A(  4,  2 ) + D(  2,  3 ) * C(  4,  2 )
                Z(  9,  7 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z( 10,  7 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z( 11,  7 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z( 12,  7 ) = B(  2,  3) * A(  4,  3 ) + D(  2,  3 ) * C(  4,  3 )
                Z(  9,  8 ) = B(  2,  3) * A(  1,  4 ) + D(  2,  3 ) * C(  1,  4 )
                Z( 10,  8 ) = B(  2,  3) * A(  2,  4 ) + D(  2,  3 ) * C(  2,  4 )
                Z( 11,  8 ) = B(  2,  3) * A(  3,  4 ) + D(  2,  3 ) * C(  3,  4 )
                Z( 12,  8 ) = B(  2,  3) * A(  4,  4 ) + D(  2,  3 ) * C(  4,  4 )
                Z(  9,  9 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z( 10,  9 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z( 11,  9 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z( 12,  9 ) = B(  3,  3) * A(  4,  1 ) + D(  3,  3 ) * C(  4,  1 )
                Z(  9, 10 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z( 10, 10 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z( 11, 10 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z( 12, 10 ) = B(  3,  3) * A(  4,  2 ) + D(  3,  3 ) * C(  4,  2 )
                Z(  9, 11 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z( 10, 11 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z( 11, 11 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )
                Z( 12, 11 ) = B(  3,  3) * A(  4,  3 ) + D(  3,  3 ) * C(  4,  3 )
                Z(  9, 12 ) = B(  3,  3) * A(  1,  4 ) + D(  3,  3 ) * C(  1,  4 )
                Z( 10, 12 ) = B(  3,  3) * A(  2,  4 ) + D(  3,  3 ) * C(  2,  4 )
                Z( 11, 12 ) = B(  3,  3) * A(  3,  4 ) + D(  3,  3 ) * C(  3,  4 )
                Z( 12, 12 ) = B(  3,  3) * A(  4,  4 ) + D(  3,  3 ) * C(  4,  4 )
                Z( 13,  1 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z( 14,  1 ) = B(  1,  4) * A(  2,  1 ) + D(  1,  4 ) * C(  2,  1 )
                Z( 15,  1 ) = B(  1,  4) * A(  3,  1 ) + D(  1,  4 ) * C(  3,  1 )
                Z( 16,  1 ) = B(  1,  4) * A(  4,  1 ) + D(  1,  4 ) * C(  4,  1 )
                Z( 13,  2 ) = B(  1,  4) * A(  1,  2 ) + D(  1,  4 ) * C(  1,  2 )
                Z( 14,  2 ) = B(  1,  4) * A(  2,  2 ) + D(  1,  4 ) * C(  2,  2 )
                Z( 15,  2 ) = B(  1,  4) * A(  3,  2 ) + D(  1,  4 ) * C(  3,  2 )
                Z( 16,  2 ) = B(  1,  4) * A(  4,  2 ) + D(  1,  4 ) * C(  4,  2 )
                Z( 13,  3 ) = B(  1,  4) * A(  1,  3 ) + D(  1,  4 ) * C(  1,  3 )
                Z( 14,  3 ) = B(  1,  4) * A(  2,  3 ) + D(  1,  4 ) * C(  2,  3 )
                Z( 15,  3 ) = B(  1,  4) * A(  3,  3 ) + D(  1,  4 ) * C(  3,  3 )
                Z( 16,  3 ) = B(  1,  4) * A(  4,  3 ) + D(  1,  4 ) * C(  4,  3 )
                Z( 13,  4 ) = B(  1,  4) * A(  1,  4 ) + D(  1,  4 ) * C(  1,  4 )
                Z( 14,  4 ) = B(  1,  4) * A(  2,  4 ) + D(  1,  4 ) * C(  2,  4 )
                Z( 15,  4 ) = B(  1,  4) * A(  3,  4 ) + D(  1,  4 ) * C(  3,  4 )
                Z( 16,  4 ) = B(  1,  4) * A(  4,  4 ) + D(  1,  4 ) * C(  4,  4 )
                Z( 13,  5 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z( 14,  5 ) = B(  2,  4) * A(  2,  1 ) + D(  2,  4 ) * C(  2,  1 )
                Z( 15,  5 ) = B(  2,  4) * A(  3,  1 ) + D(  2,  4 ) * C(  3,  1 )
                Z( 16,  5 ) = B(  2,  4) * A(  4,  1 ) + D(  2,  4 ) * C(  4,  1 )
                Z( 13,  6 ) = B(  2,  4) * A(  1,  2 ) + D(  2,  4 ) * C(  1,  2 )
                Z( 14,  6 ) = B(  2,  4) * A(  2,  2 ) + D(  2,  4 ) * C(  2,  2 )
                Z( 15,  6 ) = B(  2,  4) * A(  3,  2 ) + D(  2,  4 ) * C(  3,  2 )
                Z( 16,  6 ) = B(  2,  4) * A(  4,  2 ) + D(  2,  4 ) * C(  4,  2 )
                Z( 13,  7 ) = B(  2,  4) * A(  1,  3 ) + D(  2,  4 ) * C(  1,  3 )
                Z( 14,  7 ) = B(  2,  4) * A(  2,  3 ) + D(  2,  4 ) * C(  2,  3 )
                Z( 15,  7 ) = B(  2,  4) * A(  3,  3 ) + D(  2,  4 ) * C(  3,  3 )
                Z( 16,  7 ) = B(  2,  4) * A(  4,  3 ) + D(  2,  4 ) * C(  4,  3 )
                Z( 13,  8 ) = B(  2,  4) * A(  1,  4 ) + D(  2,  4 ) * C(  1,  4 )
                Z( 14,  8 ) = B(  2,  4) * A(  2,  4 ) + D(  2,  4 ) * C(  2,  4 )
                Z( 15,  8 ) = B(  2,  4) * A(  3,  4 ) + D(  2,  4 ) * C(  3,  4 )
                Z( 16,  8 ) = B(  2,  4) * A(  4,  4 ) + D(  2,  4 ) * C(  4,  4 )
                Z( 13,  9 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z( 14,  9 ) = B(  3,  4) * A(  2,  1 ) + D(  3,  4 ) * C(  2,  1 )
                Z( 15,  9 ) = B(  3,  4) * A(  3,  1 ) + D(  3,  4 ) * C(  3,  1 )
                Z( 16,  9 ) = B(  3,  4) * A(  4,  1 ) + D(  3,  4 ) * C(  4,  1 )
                Z( 13, 10 ) = B(  3,  4) * A(  1,  2 ) + D(  3,  4 ) * C(  1,  2 )
                Z( 14, 10 ) = B(  3,  4) * A(  2,  2 ) + D(  3,  4 ) * C(  2,  2 )
                Z( 15, 10 ) = B(  3,  4) * A(  3,  2 ) + D(  3,  4 ) * C(  3,  2 )
                Z( 16, 10 ) = B(  3,  4) * A(  4,  2 ) + D(  3,  4 ) * C(  4,  2 )
                Z( 13, 11 ) = B(  3,  4) * A(  1,  3 ) + D(  3,  4 ) * C(  1,  3 )
                Z( 14, 11 ) = B(  3,  4) * A(  2,  3 ) + D(  3,  4 ) * C(  2,  3 )
                Z( 15, 11 ) = B(  3,  4) * A(  3,  3 ) + D(  3,  4 ) * C(  3,  3 )
                Z( 16, 11 ) = B(  3,  4) * A(  4,  3 ) + D(  3,  4 ) * C(  4,  3 )
                Z( 13, 12 ) = B(  3,  4) * A(  1,  4 ) + D(  3,  4 ) * C(  1,  4 )
                Z( 14, 12 ) = B(  3,  4) * A(  2,  4 ) + D(  3,  4 ) * C(  2,  4 )
                Z( 15, 12 ) = B(  3,  4) * A(  3,  4 ) + D(  3,  4 ) * C(  3,  4 )
                Z( 16, 12 ) = B(  3,  4) * A(  4,  4 ) + D(  3,  4 ) * C(  4,  4 )
                Z( 13, 13 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )
                Z( 14, 13 ) = B(  4,  4) * A(  2,  1 ) + D(  4,  4 ) * C(  2,  1 )
                Z( 15, 13 ) = B(  4,  4) * A(  3,  1 ) + D(  4,  4 ) * C(  3,  1 )
                Z( 16, 13 ) = B(  4,  4) * A(  4,  1 ) + D(  4,  4 ) * C(  4,  1 )
                Z( 13, 14 ) = B(  4,  4) * A(  1,  2 ) + D(  4,  4 ) * C(  1,  2 )
                Z( 14, 14 ) = B(  4,  4) * A(  2,  2 ) + D(  4,  4 ) * C(  2,  2 )
                Z( 15, 14 ) = B(  4,  4) * A(  3,  2 ) + D(  4,  4 ) * C(  3,  2 )
                Z( 16, 14 ) = B(  4,  4) * A(  4,  2 ) + D(  4,  4 ) * C(  4,  2 )
                Z( 13, 15 ) = B(  4,  4) * A(  1,  3 ) + D(  4,  4 ) * C(  1,  3 )
                Z( 14, 15 ) = B(  4,  4) * A(  2,  3 ) + D(  4,  4 ) * C(  2,  3 )
                Z( 15, 15 ) = B(  4,  4) * A(  3,  3 ) + D(  4,  4 ) * C(  3,  3 )
                Z( 16, 15 ) = B(  4,  4) * A(  4,  3 ) + D(  4,  4 ) * C(  4,  3 )
                Z( 13, 16 ) = B(  4,  4) * A(  1,  4 ) + D(  4,  4 ) * C(  1,  4 )
                Z( 14, 16 ) = B(  4,  4) * A(  2,  4 ) + D(  4,  4 ) * C(  2,  4 )
                Z( 15, 16 ) = B(  4,  4) * A(  3,  4 ) + D(  4,  4 ) * C(  3,  4 )
                Z( 16, 16 ) = B(  4,  4) * A(  4,  4 ) + D(  4,  4 ) * C(  4,  4 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  5 ) = Z(  1,  5 ) + B(  2,  1) * A(  1,  1 )
                    Z(  2,  5 ) = Z(  2,  5 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  5 ) = Z(  3,  5 ) + B(  2,  1) * A(  3,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + B(  2,  1) * A(  4,  1 )
                    Z(  1,  6 ) = Z(  1,  6 ) + B(  2,  1) * A(  1,  2 )
                    Z(  2,  6 ) = Z(  2,  6 ) + B(  2,  1) * A(  2,  2 )
                    Z(  3,  6 ) = Z(  3,  6 ) + B(  2,  1) * A(  3,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + B(  2,  1) * A(  4,  2 )
                    Z(  1,  7 ) = Z(  1,  7 ) + B(  2,  1) * A(  1,  3 )
                    Z(  2,  7 ) = Z(  2,  7 ) + B(  2,  1) * A(  2,  3 )
                    Z(  3,  7 ) = Z(  3,  7 ) + B(  2,  1) * A(  3,  3 )
                    Z(  4,  7 ) = Z(  4,  7 ) + B(  2,  1) * A(  4,  3 )
                    Z(  1,  8 ) = Z(  1,  8 ) + B(  2,  1) * A(  1,  4 )
                    Z(  2,  8 ) = Z(  2,  8 ) + B(  2,  1) * A(  2,  4 )
                    Z(  3,  8 ) = Z(  3,  8 ) + B(  2,  1) * A(  3,  4 )
                    Z(  4,  8 ) = Z(  4,  8 ) + B(  2,  1) * A(  4,  4 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  9 ) = Z(  5,  9 ) + B(  3,  2) * A(  1,  1 )
                    Z(  6,  9 ) = Z(  6,  9 ) + B(  3,  2) * A(  2,  1 )
                    Z(  7,  9 ) = Z(  7,  9 ) + B(  3,  2) * A(  3,  1 )
                    Z(  8,  9 ) = Z(  8,  9 ) + B(  3,  2) * A(  4,  1 )
                    Z(  5, 10 ) = Z(  5, 10 ) + B(  3,  2) * A(  1,  2 )
                    Z(  6, 10 ) = Z(  6, 10 ) + B(  3,  2) * A(  2,  2 )
                    Z(  7, 10 ) = Z(  7, 10 ) + B(  3,  2) * A(  3,  2 )
                    Z(  8, 10 ) = Z(  8, 10 ) + B(  3,  2) * A(  4,  2 )
                    Z(  5, 11 ) = Z(  5, 11 ) + B(  3,  2) * A(  1,  3 )
                    Z(  6, 11 ) = Z(  6, 11 ) + B(  3,  2) * A(  2,  3 )
                    Z(  7, 11 ) = Z(  7, 11 ) + B(  3,  2) * A(  3,  3 )
                    Z(  8, 11 ) = Z(  8, 11 ) + B(  3,  2) * A(  4,  3 )
                    Z(  5, 12 ) = Z(  5, 12 ) + B(  3,  2) * A(  1,  4 )
                    Z(  6, 12 ) = Z(  6, 12 ) + B(  3,  2) * A(  2,  4 )
                    Z(  7, 12 ) = Z(  7, 12 ) + B(  3,  2) * A(  3,  4 )
                    Z(  8, 12 ) = Z(  8, 12 ) + B(  3,  2) * A(  4,  4 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  9, 13 ) = Z(  9, 13 ) + B(  4,  3) * A(  1,  1 )
                    Z( 10, 13 ) = Z( 10, 13 ) + B(  4,  3) * A(  2,  1 )
                    Z( 11, 13 ) = Z( 11, 13 ) + B(  4,  3) * A(  3,  1 )
                    Z( 12, 13 ) = Z( 12, 13 ) + B(  4,  3) * A(  4,  1 )
                    Z(  9, 14 ) = Z(  9, 14 ) + B(  4,  3) * A(  1,  2 )
                    Z( 10, 14 ) = Z( 10, 14 ) + B(  4,  3) * A(  2,  2 )
                    Z( 11, 14 ) = Z( 11, 14 ) + B(  4,  3) * A(  3,  2 )
                    Z( 12, 14 ) = Z( 12, 14 ) + B(  4,  3) * A(  4,  2 )
                    Z(  9, 15 ) = Z(  9, 15 ) + B(  4,  3) * A(  1,  3 )
                    Z( 10, 15 ) = Z( 10, 15 ) + B(  4,  3) * A(  2,  3 )
                    Z( 11, 15 ) = Z( 11, 15 ) + B(  4,  3) * A(  3,  3 )
                    Z( 12, 15 ) = Z( 12, 15 ) + B(  4,  3) * A(  4,  3 )
                    Z(  9, 16 ) = Z(  9, 16 ) + B(  4,  3) * A(  1,  4 )
                    Z( 10, 16 ) = Z( 10, 16 ) + B(  4,  3) * A(  2,  4 )
                    Z( 11, 16 ) = Z( 11, 16 ) + B(  4,  3) * A(  3,  4 )
                    Z( 12, 16 ) = Z( 12, 16 ) + B(  4,  3) * A(  4,  4 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  1,  5 ) = Z(  1,  5 ) + D(  2,  1) * C(  1,  1 )
                    Z(  2,  5 ) = Z(  2,  5 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  5 ) = Z(  3,  5 ) + D(  2,  1) * C(  3,  1 )
                    Z(  4,  5 ) = Z(  4,  5 ) + D(  2,  1) * C(  4,  1 )
                    Z(  1,  6 ) = Z(  1,  6 ) + D(  2,  1) * C(  1,  2 )
                    Z(  2,  6 ) = Z(  2,  6 ) + D(  2,  1) * C(  2,  2 )
                    Z(  3,  6 ) = Z(  3,  6 ) + D(  2,  1) * C(  3,  2 )
                    Z(  4,  6 ) = Z(  4,  6 ) + D(  2,  1) * C(  4,  2 )
                    Z(  1,  7 ) = Z(  1,  7 ) + D(  2,  1) * C(  1,  3 )
                    Z(  2,  7 ) = Z(  2,  7 ) + D(  2,  1) * C(  2,  3 )
                    Z(  3,  7 ) = Z(  3,  7 ) + D(  2,  1) * C(  3,  3 )
                    Z(  4,  7 ) = Z(  4,  7 ) + D(  2,  1) * C(  4,  3 )
                    Z(  1,  8 ) = Z(  1,  8 ) + D(  2,  1) * C(  1,  4 )
                    Z(  2,  8 ) = Z(  2,  8 ) + D(  2,  1) * C(  2,  4 )
                    Z(  3,  8 ) = Z(  3,  8 ) + D(  2,  1) * C(  3,  4 )
                    Z(  4,  8 ) = Z(  4,  8 ) + D(  2,  1) * C(  4,  4 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  9 ) = Z(  5,  9 ) + D(  3,  2) * C(  1,  1 )
                    Z(  6,  9 ) = Z(  6,  9 ) + D(  3,  2) * C(  2,  1 )
                    Z(  7,  9 ) = Z(  7,  9 ) + D(  3,  2) * C(  3,  1 )
                    Z(  8,  9 ) = Z(  8,  9 ) + D(  3,  2) * C(  4,  1 )
                    Z(  5, 10 ) = Z(  5, 10 ) + D(  3,  2) * C(  1,  2 )
                    Z(  6, 10 ) = Z(  6, 10 ) + D(  3,  2) * C(  2,  2 )
                    Z(  7, 10 ) = Z(  7, 10 ) + D(  3,  2) * C(  3,  2 )
                    Z(  8, 10 ) = Z(  8, 10 ) + D(  3,  2) * C(  4,  2 )
                    Z(  5, 11 ) = Z(  5, 11 ) + D(  3,  2) * C(  1,  3 )
                    Z(  6, 11 ) = Z(  6, 11 ) + D(  3,  2) * C(  2,  3 )
                    Z(  7, 11 ) = Z(  7, 11 ) + D(  3,  2) * C(  3,  3 )
                    Z(  8, 11 ) = Z(  8, 11 ) + D(  3,  2) * C(  4,  3 )
                    Z(  5, 12 ) = Z(  5, 12 ) + D(  3,  2) * C(  1,  4 )
                    Z(  6, 12 ) = Z(  6, 12 ) + D(  3,  2) * C(  2,  4 )
                    Z(  7, 12 ) = Z(  7, 12 ) + D(  3,  2) * C(  3,  4 )
                    Z(  8, 12 ) = Z(  8, 12 ) + D(  3,  2) * C(  4,  4 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  9, 13 ) = Z(  9, 13 ) + D(  4,  3) * C(  1,  1 )
                    Z( 10, 13 ) = Z( 10, 13 ) + D(  4,  3) * C(  2,  1 )
                    Z( 11, 13 ) = Z( 11, 13 ) + D(  4,  3) * C(  3,  1 )
                    Z( 12, 13 ) = Z( 12, 13 ) + D(  4,  3) * C(  4,  1 )
                    Z(  9, 14 ) = Z(  9, 14 ) + D(  4,  3) * C(  1,  2 )
                    Z( 10, 14 ) = Z( 10, 14 ) + D(  4,  3) * C(  2,  2 )
                    Z( 11, 14 ) = Z( 11, 14 ) + D(  4,  3) * C(  3,  2 )
                    Z( 12, 14 ) = Z( 12, 14 ) + D(  4,  3) * C(  4,  2 )
                    Z(  9, 15 ) = Z(  9, 15 ) + D(  4,  3) * C(  1,  3 )
                    Z( 10, 15 ) = Z( 10, 15 ) + D(  4,  3) * C(  2,  3 )
                    Z( 11, 15 ) = Z( 11, 15 ) + D(  4,  3) * C(  3,  3 )
                    Z( 12, 15 ) = Z( 12, 15 ) + D(  4,  3) * C(  4,  3 )
                    Z(  9, 16 ) = Z(  9, 16 ) + D(  4,  3) * C(  1,  4 )
                    Z( 10, 16 ) = Z( 10, 16 ) + D(  4,  3) * C(  2,  4 )
                    Z( 11, 16 ) = Z( 11, 16 ) + D(  4,  3) * C(  3,  4 )
                    Z( 12, 16 ) = Z( 12, 16 ) + D(  4,  3) * C(  4,  4 )
                END IF
            END IF
        END IF
    ELSE
        IF ( M == 1 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  1,  2 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  2 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  2,  1 ) = Z(  2,  1 ) + B(  2,  1) * A(  1,  1 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  2,  1 ) = Z(  2,  1 ) + D(  2,  1) * C(  1,  1 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  1,  2 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  2 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  1,  3 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  3 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  3,  3 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  2,  1 ) = Z(  2,  1 ) + B(  2,  1) * A(  1,  1 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  2 ) = Z(  3,  2 ) + B(  3,  2) * A(  1,  1 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  2,  1 ) = Z(  2,  1 ) + D(  2,  1) * C(  1,  1 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  2 ) = Z(  3,  2 ) + D(  3,  2) * C(  1,  1 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  1,  2 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  2 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  1,  3 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  3 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  3,  3 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  1,  4 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z(  2,  4 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z(  3,  4 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z(  4,  4 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  2,  1 ) = Z(  2,  1 ) + B(  2,  1) * A(  1,  1 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  2 ) = Z(  3,  2 ) + B(  3,  2) * A(  1,  1 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  4,  3 ) = Z(  4,  3 ) + B(  4,  3) * A(  1,  1 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  2,  1 ) = Z(  2,  1 ) + D(  2,  1) * C(  1,  1 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  3,  2 ) = Z(  3,  2 ) + D(  3,  2) * C(  1,  1 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  4,  3 ) = Z(  4,  3 ) + D(  4,  3) * C(  1,  1 )
                END IF
            END IF
        ELSE IF ( M == 2 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  1,  3 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  3 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  1,  4 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  4 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  3 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  4,  3 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  4,  4 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  3,  1 ) = Z(  3,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  4,  1 ) = Z(  4,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  2 ) = Z(  3,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  4,  2 ) = Z(  4,  2 ) + B(  2,  1) * A(  2,  2 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  3,  1 ) = Z(  3,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  4,  1 ) = Z(  4,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  2 ) = Z(  3,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  4,  2 ) = Z(  4,  2 ) + D(  2,  1) * C(  2,  2 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  1,  3 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  3 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  1,  4 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  4 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  3 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  4,  3 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  4,  4 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  1,  5 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  5 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  1,  6 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  2,  6 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  3,  5 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  4,  5 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  3,  6 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  4,  6 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  5,  5 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  5,  6 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  3,  1 ) = Z(  3,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  4,  1 ) = Z(  4,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  2 ) = Z(  3,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  4,  2 ) = Z(  4,  2 ) + B(  2,  1) * A(  2,  2 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  3,  2) * A(  1,  1 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  3,  2) * A(  2,  1 )
                    Z(  5,  4 ) = Z(  5,  4 ) + B(  3,  2) * A(  1,  2 )
                    Z(  6,  4 ) = Z(  6,  4 ) + B(  3,  2) * A(  2,  2 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  3,  1 ) = Z(  3,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  4,  1 ) = Z(  4,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  2 ) = Z(  3,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  4,  2 ) = Z(  4,  2 ) + D(  2,  1) * C(  2,  2 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  3,  2) * C(  1,  1 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  3,  2) * C(  2,  1 )
                    Z(  5,  4 ) = Z(  5,  4 ) + D(  3,  2) * C(  1,  2 )
                    Z(  6,  4 ) = Z(  6,  4 ) + D(  3,  2) * C(  2,  2 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  1,  3 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  3 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  1,  4 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  4 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  3 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  4,  3 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  4,  4 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  1,  5 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  5 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  1,  6 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  2,  6 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  3,  5 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  4,  5 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  3,  6 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  4,  6 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  5,  5 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  5,  6 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z(  1,  7 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z(  2,  7 ) = B(  1,  4) * A(  2,  1 ) + D(  1,  4 ) * C(  2,  1 )
                Z(  1,  8 ) = B(  1,  4) * A(  1,  2 ) + D(  1,  4 ) * C(  1,  2 )
                Z(  2,  8 ) = B(  1,  4) * A(  2,  2 ) + D(  1,  4 ) * C(  2,  2 )
                Z(  3,  7 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z(  4,  7 ) = B(  2,  4) * A(  2,  1 ) + D(  2,  4 ) * C(  2,  1 )
                Z(  3,  8 ) = B(  2,  4) * A(  1,  2 ) + D(  2,  4 ) * C(  1,  2 )
                Z(  4,  8 ) = B(  2,  4) * A(  2,  2 ) + D(  2,  4 ) * C(  2,  2 )
                Z(  5,  7 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z(  6,  7 ) = B(  3,  4) * A(  2,  1 ) + D(  3,  4 ) * C(  2,  1 )
                Z(  5,  8 ) = B(  3,  4) * A(  1,  2 ) + D(  3,  4 ) * C(  1,  2 )
                Z(  6,  8 ) = B(  3,  4) * A(  2,  2 ) + D(  3,  4 ) * C(  2,  2 )
                Z(  7,  7 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )
                Z(  8,  7 ) = B(  4,  4) * A(  2,  1 ) + D(  4,  4 ) * C(  2,  1 )
                Z(  7,  8 ) = B(  4,  4) * A(  1,  2 ) + D(  4,  4 ) * C(  1,  2 )
                Z(  8,  8 ) = B(  4,  4) * A(  2,  2 ) + D(  4,  4 ) * C(  2,  2 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  3,  1 ) = Z(  3,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  4,  1 ) = Z(  4,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  3,  2 ) = Z(  3,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  4,  2 ) = Z(  4,  2 ) + B(  2,  1) * A(  2,  2 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  3,  2) * A(  1,  1 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  3,  2) * A(  2,  1 )
                    Z(  5,  4 ) = Z(  5,  4 ) + B(  3,  2) * A(  1,  2 )
                    Z(  6,  4 ) = Z(  6,  4 ) + B(  3,  2) * A(  2,  2 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  7,  5 ) = Z(  7,  5 ) + B(  4,  3) * A(  1,  1 )
                    Z(  8,  5 ) = Z(  8,  5 ) + B(  4,  3) * A(  2,  1 )
                    Z(  7,  6 ) = Z(  7,  6 ) + B(  4,  3) * A(  1,  2 )
                    Z(  8,  6 ) = Z(  8,  6 ) + B(  4,  3) * A(  2,  2 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  3,  1 ) = Z(  3,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  4,  1 ) = Z(  4,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  3,  2 ) = Z(  3,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  4,  2 ) = Z(  4,  2 ) + D(  2,  1) * C(  2,  2 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  3,  2) * C(  1,  1 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  3,  2) * C(  2,  1 )
                    Z(  5,  4 ) = Z(  5,  4 ) + D(  3,  2) * C(  1,  2 )
                    Z(  6,  4 ) = Z(  6,  4 ) + D(  3,  2) * C(  2,  2 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z(  7,  5 ) = Z(  7,  5 ) + D(  4,  3) * C(  1,  1 )
                    Z(  8,  5 ) = Z(  8,  5 ) + D(  4,  3) * C(  2,  1 )
                    Z(  7,  6 ) = Z(  7,  6 ) + D(  4,  3) * C(  1,  2 )
                    Z(  8,  6 ) = Z(  8,  6 ) + D(  4,  3) * C(  2,  2 )
                END IF
            END IF
        ELSE IF ( M == 3 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  1,  4 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  4 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  1,  5 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  5 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  5 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  1,  6 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  2,  6 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  3,  6 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  4 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  5,  4 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  6,  4 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  5,  5 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  6,  5 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  5,  6 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  6,  6 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  4,  1 ) = Z(  4,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  5,  1 ) = Z(  5,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + B(  2,  1) * A(  3,  1 )
                    Z(  4,  2 ) = Z(  4,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  5,  2 ) = Z(  5,  2 ) + B(  2,  1) * A(  2,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + B(  2,  1) * A(  3,  2 )
                    Z(  4,  3 ) = Z(  4,  3 ) + B(  2,  1) * A(  1,  3 )
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  2,  1) * A(  2,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  2,  1) * A(  3,  3 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  4,  1 ) = Z(  4,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  5,  1 ) = Z(  5,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + D(  2,  1) * C(  3,  1 )
                    Z(  4,  2 ) = Z(  4,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  5,  2 ) = Z(  5,  2 ) + D(  2,  1) * C(  2,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + D(  2,  1) * C(  3,  2 )
                    Z(  4,  3 ) = Z(  4,  3 ) + D(  2,  1) * C(  1,  3 )
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  2,  1) * C(  2,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  2,  1) * C(  3,  3 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  1,  4 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  4 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  1,  5 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  5 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  5 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  1,  6 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  2,  6 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  3,  6 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  4 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  5,  4 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  6,  4 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  5,  5 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  6,  5 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  5,  6 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  6,  6 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  1,  7 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  7 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  3,  7 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z(  1,  8 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  2,  8 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  3,  8 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z(  1,  9 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z(  2,  9 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z(  3,  9 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z(  4,  7 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  5,  7 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  6,  7 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z(  4,  8 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  5,  8 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  6,  8 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z(  4,  9 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z(  5,  9 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z(  6,  9 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z(  7,  7 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  8,  7 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  9,  7 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z(  7,  8 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  8,  8 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z(  9,  8 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z(  7,  9 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z(  8,  9 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z(  9,  9 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  4,  1 ) = Z(  4,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  5,  1 ) = Z(  5,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + B(  2,  1) * A(  3,  1 )
                    Z(  4,  2 ) = Z(  4,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  5,  2 ) = Z(  5,  2 ) + B(  2,  1) * A(  2,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + B(  2,  1) * A(  3,  2 )
                    Z(  4,  3 ) = Z(  4,  3 ) + B(  2,  1) * A(  1,  3 )
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  2,  1) * A(  2,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  2,  1) * A(  3,  3 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  7,  4 ) = Z(  7,  4 ) + B(  3,  2) * A(  1,  1 )
                    Z(  8,  4 ) = Z(  8,  4 ) + B(  3,  2) * A(  2,  1 )
                    Z(  9,  4 ) = Z(  9,  4 ) + B(  3,  2) * A(  3,  1 )
                    Z(  7,  5 ) = Z(  7,  5 ) + B(  3,  2) * A(  1,  2 )
                    Z(  8,  5 ) = Z(  8,  5 ) + B(  3,  2) * A(  2,  2 )
                    Z(  9,  5 ) = Z(  9,  5 ) + B(  3,  2) * A(  3,  2 )
                    Z(  7,  6 ) = Z(  7,  6 ) + B(  3,  2) * A(  1,  3 )
                    Z(  8,  6 ) = Z(  8,  6 ) + B(  3,  2) * A(  2,  3 )
                    Z(  9,  6 ) = Z(  9,  6 ) + B(  3,  2) * A(  3,  3 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  4,  1 ) = Z(  4,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  5,  1 ) = Z(  5,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + D(  2,  1) * C(  3,  1 )
                    Z(  4,  2 ) = Z(  4,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  5,  2 ) = Z(  5,  2 ) + D(  2,  1) * C(  2,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + D(  2,  1) * C(  3,  2 )
                    Z(  4,  3 ) = Z(  4,  3 ) + D(  2,  1) * C(  1,  3 )
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  2,  1) * C(  2,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  2,  1) * C(  3,  3 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  7,  4 ) = Z(  7,  4 ) + D(  3,  2) * C(  1,  1 )
                    Z(  8,  4 ) = Z(  8,  4 ) + D(  3,  2) * C(  2,  1 )
                    Z(  9,  4 ) = Z(  9,  4 ) + D(  3,  2) * C(  3,  1 )
                    Z(  7,  5 ) = Z(  7,  5 ) + D(  3,  2) * C(  1,  2 )
                    Z(  8,  5 ) = Z(  8,  5 ) + D(  3,  2) * C(  2,  2 )
                    Z(  9,  5 ) = Z(  9,  5 ) + D(  3,  2) * C(  3,  2 )
                    Z(  7,  6 ) = Z(  7,  6 ) + D(  3,  2) * C(  1,  3 )
                    Z(  8,  6 ) = Z(  8,  6 ) + D(  3,  2) * C(  2,  3 )
                    Z(  9,  6 ) = Z(  9,  6 ) + D(  3,  2) * C(  3,  3 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  1,  4 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  4 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  4 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  1,  5 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  5 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  5 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  1,  6 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  2,  6 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  3,  6 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  4 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  5,  4 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  6,  4 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  5,  5 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  6,  5 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  5,  6 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  6,  6 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  1,  7 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  7 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  3,  7 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z(  1,  8 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  2,  8 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  3,  8 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z(  1,  9 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z(  2,  9 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z(  3,  9 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z(  4,  7 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  5,  7 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  6,  7 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z(  4,  8 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  5,  8 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  6,  8 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z(  4,  9 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z(  5,  9 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z(  6,  9 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z(  7,  7 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z(  8,  7 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z(  9,  7 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z(  7,  8 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z(  8,  8 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z(  9,  8 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z(  7,  9 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z(  8,  9 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z(  9,  9 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )
                Z(  1, 10 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z(  2, 10 ) = B(  1,  4) * A(  2,  1 ) + D(  1,  4 ) * C(  2,  1 )
                Z(  3, 10 ) = B(  1,  4) * A(  3,  1 ) + D(  1,  4 ) * C(  3,  1 )
                Z(  1, 11 ) = B(  1,  4) * A(  1,  2 ) + D(  1,  4 ) * C(  1,  2 )
                Z(  2, 11 ) = B(  1,  4) * A(  2,  2 ) + D(  1,  4 ) * C(  2,  2 )
                Z(  3, 11 ) = B(  1,  4) * A(  3,  2 ) + D(  1,  4 ) * C(  3,  2 )
                Z(  1, 12 ) = B(  1,  4) * A(  1,  3 ) + D(  1,  4 ) * C(  1,  3 )
                Z(  2, 12 ) = B(  1,  4) * A(  2,  3 ) + D(  1,  4 ) * C(  2,  3 )
                Z(  3, 12 ) = B(  1,  4) * A(  3,  3 ) + D(  1,  4 ) * C(  3,  3 )
                Z(  4, 10 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z(  5, 10 ) = B(  2,  4) * A(  2,  1 ) + D(  2,  4 ) * C(  2,  1 )
                Z(  6, 10 ) = B(  2,  4) * A(  3,  1 ) + D(  2,  4 ) * C(  3,  1 )
                Z(  4, 11 ) = B(  2,  4) * A(  1,  2 ) + D(  2,  4 ) * C(  1,  2 )
                Z(  5, 11 ) = B(  2,  4) * A(  2,  2 ) + D(  2,  4 ) * C(  2,  2 )
                Z(  6, 11 ) = B(  2,  4) * A(  3,  2 ) + D(  2,  4 ) * C(  3,  2 )
                Z(  4, 12 ) = B(  2,  4) * A(  1,  3 ) + D(  2,  4 ) * C(  1,  3 )
                Z(  5, 12 ) = B(  2,  4) * A(  2,  3 ) + D(  2,  4 ) * C(  2,  3 )
                Z(  6, 12 ) = B(  2,  4) * A(  3,  3 ) + D(  2,  4 ) * C(  3,  3 )
                Z(  7, 10 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z(  8, 10 ) = B(  3,  4) * A(  2,  1 ) + D(  3,  4 ) * C(  2,  1 )
                Z(  9, 10 ) = B(  3,  4) * A(  3,  1 ) + D(  3,  4 ) * C(  3,  1 )
                Z(  7, 11 ) = B(  3,  4) * A(  1,  2 ) + D(  3,  4 ) * C(  1,  2 )
                Z(  8, 11 ) = B(  3,  4) * A(  2,  2 ) + D(  3,  4 ) * C(  2,  2 )
                Z(  9, 11 ) = B(  3,  4) * A(  3,  2 ) + D(  3,  4 ) * C(  3,  2 )
                Z(  7, 12 ) = B(  3,  4) * A(  1,  3 ) + D(  3,  4 ) * C(  1,  3 )
                Z(  8, 12 ) = B(  3,  4) * A(  2,  3 ) + D(  3,  4 ) * C(  2,  3 )
                Z(  9, 12 ) = B(  3,  4) * A(  3,  3 ) + D(  3,  4 ) * C(  3,  3 )
                Z( 10, 10 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )
                Z( 11, 10 ) = B(  4,  4) * A(  2,  1 ) + D(  4,  4 ) * C(  2,  1 )
                Z( 12, 10 ) = B(  4,  4) * A(  3,  1 ) + D(  4,  4 ) * C(  3,  1 )
                Z( 10, 11 ) = B(  4,  4) * A(  1,  2 ) + D(  4,  4 ) * C(  1,  2 )
                Z( 11, 11 ) = B(  4,  4) * A(  2,  2 ) + D(  4,  4 ) * C(  2,  2 )
                Z( 12, 11 ) = B(  4,  4) * A(  3,  2 ) + D(  4,  4 ) * C(  3,  2 )
                Z( 10, 12 ) = B(  4,  4) * A(  1,  3 ) + D(  4,  4 ) * C(  1,  3 )
                Z( 11, 12 ) = B(  4,  4) * A(  2,  3 ) + D(  4,  4 ) * C(  2,  3 )
                Z( 12, 12 ) = B(  4,  4) * A(  3,  3 ) + D(  4,  4 ) * C(  3,  3 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  4,  1 ) = Z(  4,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  5,  1 ) = Z(  5,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + B(  2,  1) * A(  3,  1 )
                    Z(  4,  2 ) = Z(  4,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  5,  2 ) = Z(  5,  2 ) + B(  2,  1) * A(  2,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + B(  2,  1) * A(  3,  2 )
                    Z(  4,  3 ) = Z(  4,  3 ) + B(  2,  1) * A(  1,  3 )
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  2,  1) * A(  2,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  2,  1) * A(  3,  3 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  7,  4 ) = Z(  7,  4 ) + B(  3,  2) * A(  1,  1 )
                    Z(  8,  4 ) = Z(  8,  4 ) + B(  3,  2) * A(  2,  1 )
                    Z(  9,  4 ) = Z(  9,  4 ) + B(  3,  2) * A(  3,  1 )
                    Z(  7,  5 ) = Z(  7,  5 ) + B(  3,  2) * A(  1,  2 )
                    Z(  8,  5 ) = Z(  8,  5 ) + B(  3,  2) * A(  2,  2 )
                    Z(  9,  5 ) = Z(  9,  5 ) + B(  3,  2) * A(  3,  2 )
                    Z(  7,  6 ) = Z(  7,  6 ) + B(  3,  2) * A(  1,  3 )
                    Z(  8,  6 ) = Z(  8,  6 ) + B(  3,  2) * A(  2,  3 )
                    Z(  9,  6 ) = Z(  9,  6 ) + B(  3,  2) * A(  3,  3 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z( 10,  7 ) = Z( 10,  7 ) + B(  4,  3) * A(  1,  1 )
                    Z( 11,  7 ) = Z( 11,  7 ) + B(  4,  3) * A(  2,  1 )
                    Z( 12,  7 ) = Z( 12,  7 ) + B(  4,  3) * A(  3,  1 )
                    Z( 10,  8 ) = Z( 10,  8 ) + B(  4,  3) * A(  1,  2 )
                    Z( 11,  8 ) = Z( 11,  8 ) + B(  4,  3) * A(  2,  2 )
                    Z( 12,  8 ) = Z( 12,  8 ) + B(  4,  3) * A(  3,  2 )
                    Z( 10,  9 ) = Z( 10,  9 ) + B(  4,  3) * A(  1,  3 )
                    Z( 11,  9 ) = Z( 11,  9 ) + B(  4,  3) * A(  2,  3 )
                    Z( 12,  9 ) = Z( 12,  9 ) + B(  4,  3) * A(  3,  3 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  4,  1 ) = Z(  4,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  5,  1 ) = Z(  5,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + D(  2,  1) * C(  3,  1 )
                    Z(  4,  2 ) = Z(  4,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  5,  2 ) = Z(  5,  2 ) + D(  2,  1) * C(  2,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + D(  2,  1) * C(  3,  2 )
                    Z(  4,  3 ) = Z(  4,  3 ) + D(  2,  1) * C(  1,  3 )
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  2,  1) * C(  2,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  2,  1) * C(  3,  3 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  7,  4 ) = Z(  7,  4 ) + D(  3,  2) * C(  1,  1 )
                    Z(  8,  4 ) = Z(  8,  4 ) + D(  3,  2) * C(  2,  1 )
                    Z(  9,  4 ) = Z(  9,  4 ) + D(  3,  2) * C(  3,  1 )
                    Z(  7,  5 ) = Z(  7,  5 ) + D(  3,  2) * C(  1,  2 )
                    Z(  8,  5 ) = Z(  8,  5 ) + D(  3,  2) * C(  2,  2 )
                    Z(  9,  5 ) = Z(  9,  5 ) + D(  3,  2) * C(  3,  2 )
                    Z(  7,  6 ) = Z(  7,  6 ) + D(  3,  2) * C(  1,  3 )
                    Z(  8,  6 ) = Z(  8,  6 ) + D(  3,  2) * C(  2,  3 )
                    Z(  9,  6 ) = Z(  9,  6 ) + D(  3,  2) * C(  3,  3 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z( 10,  7 ) = Z( 10,  7 ) + D(  4,  3) * C(  1,  1 )
                    Z( 11,  7 ) = Z( 11,  7 ) + D(  4,  3) * C(  2,  1 )
                    Z( 12,  7 ) = Z( 12,  7 ) + D(  4,  3) * C(  3,  1 )
                    Z( 10,  8 ) = Z( 10,  8 ) + D(  4,  3) * C(  1,  2 )
                    Z( 11,  8 ) = Z( 11,  8 ) + D(  4,  3) * C(  2,  2 )
                    Z( 12,  8 ) = Z( 12,  8 ) + D(  4,  3) * C(  3,  2 )
                    Z( 10,  9 ) = Z( 10,  9 ) + D(  4,  3) * C(  1,  3 )
                    Z( 11,  9 ) = Z( 11,  9 ) + D(  4,  3) * C(  2,  3 )
                    Z( 12,  9 ) = Z( 12,  9 ) + D(  4,  3) * C(  3,  3 )
                END IF
            END IF
        ELSE IF ( M == 4 ) THEN
            IF ( N == 1 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )


            ELSE IF ( N == 2 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )
                Z(  1,  5 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  5 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  5 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  1,  2) * A(  4,  1 ) + D(  1,  2 ) * C(  4,  1 )
                Z(  1,  6 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  6 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  6 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  1,  2) * A(  4,  2 ) + D(  1,  2 ) * C(  4,  2 )
                Z(  1,  7 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  2,  7 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  3,  7 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  7 ) = B(  1,  2) * A(  4,  3 ) + D(  1,  2 ) * C(  4,  3 )
                Z(  1,  8 ) = B(  1,  2) * A(  1,  4 ) + D(  1,  2 ) * C(  1,  4 )
                Z(  2,  8 ) = B(  1,  2) * A(  2,  4 ) + D(  1,  2 ) * C(  2,  4 )
                Z(  3,  8 ) = B(  1,  2) * A(  3,  4 ) + D(  1,  2 ) * C(  3,  4 )
                Z(  4,  8 ) = B(  1,  2) * A(  4,  4 ) + D(  1,  2 ) * C(  4,  4 )
                Z(  5,  5 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  7,  5 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  8,  5 ) = B(  2,  2) * A(  4,  1 ) + D(  2,  2 ) * C(  4,  1 )
                Z(  5,  6 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  7,  6 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  8,  6 ) = B(  2,  2) * A(  4,  2 ) + D(  2,  2 ) * C(  4,  2 )
                Z(  5,  7 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  6,  7 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  7,  7 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  8,  7 ) = B(  2,  2) * A(  4,  3 ) + D(  2,  2 ) * C(  4,  3 )
                Z(  5,  8 ) = B(  2,  2) * A(  1,  4 ) + D(  2,  2 ) * C(  1,  4 )
                Z(  6,  8 ) = B(  2,  2) * A(  2,  4 ) + D(  2,  2 ) * C(  2,  4 )
                Z(  7,  8 ) = B(  2,  2) * A(  3,  4 ) + D(  2,  2 ) * C(  3,  4 )
                Z(  8,  8 ) = B(  2,  2) * A(  4,  4 ) + D(  2,  2 ) * C(  4,  4 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  5,  1 ) = Z(  5,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  7,  1 ) = Z(  7,  1 ) + B(  2,  1) * A(  3,  1 )
                    Z(  8,  1 ) = Z(  8,  1 ) + B(  2,  1) * A(  4,  1 )
                    Z(  5,  2 ) = Z(  5,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + B(  2,  1) * A(  2,  2 )
                    Z(  7,  2 ) = Z(  7,  2 ) + B(  2,  1) * A(  3,  2 )
                    Z(  8,  2 ) = Z(  8,  2 ) + B(  2,  1) * A(  4,  2 )
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  2,  1) * A(  1,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  2,  1) * A(  2,  3 )
                    Z(  7,  3 ) = Z(  7,  3 ) + B(  2,  1) * A(  3,  3 )
                    Z(  8,  3 ) = Z(  8,  3 ) + B(  2,  1) * A(  4,  3 )
                    Z(  5,  4 ) = Z(  5,  4 ) + B(  2,  1) * A(  1,  4 )
                    Z(  6,  4 ) = Z(  6,  4 ) + B(  2,  1) * A(  2,  4 )
                    Z(  7,  4 ) = Z(  7,  4 ) + B(  2,  1) * A(  3,  4 )
                    Z(  8,  4 ) = Z(  8,  4 ) + B(  2,  1) * A(  4,  4 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  5,  1 ) = Z(  5,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  7,  1 ) = Z(  7,  1 ) + D(  2,  1) * C(  3,  1 )
                    Z(  8,  1 ) = Z(  8,  1 ) + D(  2,  1) * C(  4,  1 )
                    Z(  5,  2 ) = Z(  5,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + D(  2,  1) * C(  2,  2 )
                    Z(  7,  2 ) = Z(  7,  2 ) + D(  2,  1) * C(  3,  2 )
                    Z(  8,  2 ) = Z(  8,  2 ) + D(  2,  1) * C(  4,  2 )
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  2,  1) * C(  1,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  2,  1) * C(  2,  3 )
                    Z(  7,  3 ) = Z(  7,  3 ) + D(  2,  1) * C(  3,  3 )
                    Z(  8,  3 ) = Z(  8,  3 ) + D(  2,  1) * C(  4,  3 )
                    Z(  5,  4 ) = Z(  5,  4 ) + D(  2,  1) * C(  1,  4 )
                    Z(  6,  4 ) = Z(  6,  4 ) + D(  2,  1) * C(  2,  4 )
                    Z(  7,  4 ) = Z(  7,  4 ) + D(  2,  1) * C(  3,  4 )
                    Z(  8,  4 ) = Z(  8,  4 ) + D(  2,  1) * C(  4,  4 )
                END IF
            ELSE IF ( N == 3 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )
                Z(  1,  5 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  5 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  5 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  1,  2) * A(  4,  1 ) + D(  1,  2 ) * C(  4,  1 )
                Z(  1,  6 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  6 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  6 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  1,  2) * A(  4,  2 ) + D(  1,  2 ) * C(  4,  2 )
                Z(  1,  7 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  2,  7 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  3,  7 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  7 ) = B(  1,  2) * A(  4,  3 ) + D(  1,  2 ) * C(  4,  3 )
                Z(  1,  8 ) = B(  1,  2) * A(  1,  4 ) + D(  1,  2 ) * C(  1,  4 )
                Z(  2,  8 ) = B(  1,  2) * A(  2,  4 ) + D(  1,  2 ) * C(  2,  4 )
                Z(  3,  8 ) = B(  1,  2) * A(  3,  4 ) + D(  1,  2 ) * C(  3,  4 )
                Z(  4,  8 ) = B(  1,  2) * A(  4,  4 ) + D(  1,  2 ) * C(  4,  4 )
                Z(  5,  5 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  7,  5 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  8,  5 ) = B(  2,  2) * A(  4,  1 ) + D(  2,  2 ) * C(  4,  1 )
                Z(  5,  6 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  7,  6 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  8,  6 ) = B(  2,  2) * A(  4,  2 ) + D(  2,  2 ) * C(  4,  2 )
                Z(  5,  7 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  6,  7 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  7,  7 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  8,  7 ) = B(  2,  2) * A(  4,  3 ) + D(  2,  2 ) * C(  4,  3 )
                Z(  5,  8 ) = B(  2,  2) * A(  1,  4 ) + D(  2,  2 ) * C(  1,  4 )
                Z(  6,  8 ) = B(  2,  2) * A(  2,  4 ) + D(  2,  2 ) * C(  2,  4 )
                Z(  7,  8 ) = B(  2,  2) * A(  3,  4 ) + D(  2,  2 ) * C(  3,  4 )
                Z(  8,  8 ) = B(  2,  2) * A(  4,  4 ) + D(  2,  2 ) * C(  4,  4 )
                Z(  1,  9 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  9 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  3,  9 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z(  4,  9 ) = B(  1,  3) * A(  4,  1 ) + D(  1,  3 ) * C(  4,  1 )
                Z(  1, 10 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  2, 10 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  3, 10 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z(  4, 10 ) = B(  1,  3) * A(  4,  2 ) + D(  1,  3 ) * C(  4,  2 )
                Z(  1, 11 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z(  2, 11 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z(  3, 11 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z(  4, 11 ) = B(  1,  3) * A(  4,  3 ) + D(  1,  3 ) * C(  4,  3 )
                Z(  1, 12 ) = B(  1,  3) * A(  1,  4 ) + D(  1,  3 ) * C(  1,  4 )
                Z(  2, 12 ) = B(  1,  3) * A(  2,  4 ) + D(  1,  3 ) * C(  2,  4 )
                Z(  3, 12 ) = B(  1,  3) * A(  3,  4 ) + D(  1,  3 ) * C(  3,  4 )
                Z(  4, 12 ) = B(  1,  3) * A(  4,  4 ) + D(  1,  3 ) * C(  4,  4 )
                Z(  5,  9 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  6,  9 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  7,  9 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z(  8,  9 ) = B(  2,  3) * A(  4,  1 ) + D(  2,  3 ) * C(  4,  1 )
                Z(  5, 10 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  6, 10 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  7, 10 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z(  8, 10 ) = B(  2,  3) * A(  4,  2 ) + D(  2,  3 ) * C(  4,  2 )
                Z(  5, 11 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z(  6, 11 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z(  7, 11 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z(  8, 11 ) = B(  2,  3) * A(  4,  3 ) + D(  2,  3 ) * C(  4,  3 )
                Z(  5, 12 ) = B(  2,  3) * A(  1,  4 ) + D(  2,  3 ) * C(  1,  4 )
                Z(  6, 12 ) = B(  2,  3) * A(  2,  4 ) + D(  2,  3 ) * C(  2,  4 )
                Z(  7, 12 ) = B(  2,  3) * A(  3,  4 ) + D(  2,  3 ) * C(  3,  4 )
                Z(  8, 12 ) = B(  2,  3) * A(  4,  4 ) + D(  2,  3 ) * C(  4,  4 )
                Z(  9,  9 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z( 10,  9 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z( 11,  9 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z( 12,  9 ) = B(  3,  3) * A(  4,  1 ) + D(  3,  3 ) * C(  4,  1 )
                Z(  9, 10 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z( 10, 10 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z( 11, 10 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z( 12, 10 ) = B(  3,  3) * A(  4,  2 ) + D(  3,  3 ) * C(  4,  2 )
                Z(  9, 11 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z( 10, 11 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z( 11, 11 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )
                Z( 12, 11 ) = B(  3,  3) * A(  4,  3 ) + D(  3,  3 ) * C(  4,  3 )
                Z(  9, 12 ) = B(  3,  3) * A(  1,  4 ) + D(  3,  3 ) * C(  1,  4 )
                Z( 10, 12 ) = B(  3,  3) * A(  2,  4 ) + D(  3,  3 ) * C(  2,  4 )
                Z( 11, 12 ) = B(  3,  3) * A(  3,  4 ) + D(  3,  3 ) * C(  3,  4 )
                Z( 12, 12 ) = B(  3,  3) * A(  4,  4 ) + D(  3,  3 ) * C(  4,  4 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  5,  1 ) = Z(  5,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  7,  1 ) = Z(  7,  1 ) + B(  2,  1) * A(  3,  1 )
                    Z(  8,  1 ) = Z(  8,  1 ) + B(  2,  1) * A(  4,  1 )
                    Z(  5,  2 ) = Z(  5,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + B(  2,  1) * A(  2,  2 )
                    Z(  7,  2 ) = Z(  7,  2 ) + B(  2,  1) * A(  3,  2 )
                    Z(  8,  2 ) = Z(  8,  2 ) + B(  2,  1) * A(  4,  2 )
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  2,  1) * A(  1,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  2,  1) * A(  2,  3 )
                    Z(  7,  3 ) = Z(  7,  3 ) + B(  2,  1) * A(  3,  3 )
                    Z(  8,  3 ) = Z(  8,  3 ) + B(  2,  1) * A(  4,  3 )
                    Z(  5,  4 ) = Z(  5,  4 ) + B(  2,  1) * A(  1,  4 )
                    Z(  6,  4 ) = Z(  6,  4 ) + B(  2,  1) * A(  2,  4 )
                    Z(  7,  4 ) = Z(  7,  4 ) + B(  2,  1) * A(  3,  4 )
                    Z(  8,  4 ) = Z(  8,  4 ) + B(  2,  1) * A(  4,  4 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  9,  5 ) = Z(  9,  5 ) + B(  3,  2) * A(  1,  1 )
                    Z( 10,  5 ) = Z( 10,  5 ) + B(  3,  2) * A(  2,  1 )
                    Z( 11,  5 ) = Z( 11,  5 ) + B(  3,  2) * A(  3,  1 )
                    Z( 12,  5 ) = Z( 12,  5 ) + B(  3,  2) * A(  4,  1 )
                    Z(  9,  6 ) = Z(  9,  6 ) + B(  3,  2) * A(  1,  2 )
                    Z( 10,  6 ) = Z( 10,  6 ) + B(  3,  2) * A(  2,  2 )
                    Z( 11,  6 ) = Z( 11,  6 ) + B(  3,  2) * A(  3,  2 )
                    Z( 12,  6 ) = Z( 12,  6 ) + B(  3,  2) * A(  4,  2 )
                    Z(  9,  7 ) = Z(  9,  7 ) + B(  3,  2) * A(  1,  3 )
                    Z( 10,  7 ) = Z( 10,  7 ) + B(  3,  2) * A(  2,  3 )
                    Z( 11,  7 ) = Z( 11,  7 ) + B(  3,  2) * A(  3,  3 )
                    Z( 12,  7 ) = Z( 12,  7 ) + B(  3,  2) * A(  4,  3 )
                    Z(  9,  8 ) = Z(  9,  8 ) + B(  3,  2) * A(  1,  4 )
                    Z( 10,  8 ) = Z( 10,  8 ) + B(  3,  2) * A(  2,  4 )
                    Z( 11,  8 ) = Z( 11,  8 ) + B(  3,  2) * A(  3,  4 )
                    Z( 12,  8 ) = Z( 12,  8 ) + B(  3,  2) * A(  4,  4 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  5,  1 ) = Z(  5,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  7,  1 ) = Z(  7,  1 ) + D(  2,  1) * C(  3,  1 )
                    Z(  8,  1 ) = Z(  8,  1 ) + D(  2,  1) * C(  4,  1 )
                    Z(  5,  2 ) = Z(  5,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + D(  2,  1) * C(  2,  2 )
                    Z(  7,  2 ) = Z(  7,  2 ) + D(  2,  1) * C(  3,  2 )
                    Z(  8,  2 ) = Z(  8,  2 ) + D(  2,  1) * C(  4,  2 )
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  2,  1) * C(  1,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  2,  1) * C(  2,  3 )
                    Z(  7,  3 ) = Z(  7,  3 ) + D(  2,  1) * C(  3,  3 )
                    Z(  8,  3 ) = Z(  8,  3 ) + D(  2,  1) * C(  4,  3 )
                    Z(  5,  4 ) = Z(  5,  4 ) + D(  2,  1) * C(  1,  4 )
                    Z(  6,  4 ) = Z(  6,  4 ) + D(  2,  1) * C(  2,  4 )
                    Z(  7,  4 ) = Z(  7,  4 ) + D(  2,  1) * C(  3,  4 )
                    Z(  8,  4 ) = Z(  8,  4 ) + D(  2,  1) * C(  4,  4 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  9,  5 ) = Z(  9,  5 ) + D(  3,  2) * C(  1,  1 )
                    Z( 10,  5 ) = Z( 10,  5 ) + D(  3,  2) * C(  2,  1 )
                    Z( 11,  5 ) = Z( 11,  5 ) + D(  3,  2) * C(  3,  1 )
                    Z( 12,  5 ) = Z( 12,  5 ) + D(  3,  2) * C(  4,  1 )
                    Z(  9,  6 ) = Z(  9,  6 ) + D(  3,  2) * C(  1,  2 )
                    Z( 10,  6 ) = Z( 10,  6 ) + D(  3,  2) * C(  2,  2 )
                    Z( 11,  6 ) = Z( 11,  6 ) + D(  3,  2) * C(  3,  2 )
                    Z( 12,  6 ) = Z( 12,  6 ) + D(  3,  2) * C(  4,  2 )
                    Z(  9,  7 ) = Z(  9,  7 ) + D(  3,  2) * C(  1,  3 )
                    Z( 10,  7 ) = Z( 10,  7 ) + D(  3,  2) * C(  2,  3 )
                    Z( 11,  7 ) = Z( 11,  7 ) + D(  3,  2) * C(  3,  3 )
                    Z( 12,  7 ) = Z( 12,  7 ) + D(  3,  2) * C(  4,  3 )
                    Z(  9,  8 ) = Z(  9,  8 ) + D(  3,  2) * C(  1,  4 )
                    Z( 10,  8 ) = Z( 10,  8 ) + D(  3,  2) * C(  2,  4 )
                    Z( 11,  8 ) = Z( 11,  8 ) + D(  3,  2) * C(  3,  4 )
                    Z( 12,  8 ) = Z( 12,  8 ) + D(  3,  2) * C(  4,  4 )
                END IF
            ELSE IF ( N == 4 ) THEN
                Z(  1,  1 ) = B(  1,  1) * A(  1,  1 ) + D(  1,  1 ) * C(  1,  1 )
                Z(  2,  1 ) = B(  1,  1) * A(  2,  1 ) + D(  1,  1 ) * C(  2,  1 )
                Z(  3,  1 ) = B(  1,  1) * A(  3,  1 ) + D(  1,  1 ) * C(  3,  1 )
                Z(  4,  1 ) = B(  1,  1) * A(  4,  1 ) + D(  1,  1 ) * C(  4,  1 )
                Z(  1,  2 ) = B(  1,  1) * A(  1,  2 ) + D(  1,  1 ) * C(  1,  2 )
                Z(  2,  2 ) = B(  1,  1) * A(  2,  2 ) + D(  1,  1 ) * C(  2,  2 )
                Z(  3,  2 ) = B(  1,  1) * A(  3,  2 ) + D(  1,  1 ) * C(  3,  2 )
                Z(  4,  2 ) = B(  1,  1) * A(  4,  2 ) + D(  1,  1 ) * C(  4,  2 )
                Z(  1,  3 ) = B(  1,  1) * A(  1,  3 ) + D(  1,  1 ) * C(  1,  3 )
                Z(  2,  3 ) = B(  1,  1) * A(  2,  3 ) + D(  1,  1 ) * C(  2,  3 )
                Z(  3,  3 ) = B(  1,  1) * A(  3,  3 ) + D(  1,  1 ) * C(  3,  3 )
                Z(  4,  3 ) = B(  1,  1) * A(  4,  3 ) + D(  1,  1 ) * C(  4,  3 )
                Z(  1,  4 ) = B(  1,  1) * A(  1,  4 ) + D(  1,  1 ) * C(  1,  4 )
                Z(  2,  4 ) = B(  1,  1) * A(  2,  4 ) + D(  1,  1 ) * C(  2,  4 )
                Z(  3,  4 ) = B(  1,  1) * A(  3,  4 ) + D(  1,  1 ) * C(  3,  4 )
                Z(  4,  4 ) = B(  1,  1) * A(  4,  4 ) + D(  1,  1 ) * C(  4,  4 )
                Z(  1,  5 ) = B(  1,  2) * A(  1,  1 ) + D(  1,  2 ) * C(  1,  1 )
                Z(  2,  5 ) = B(  1,  2) * A(  2,  1 ) + D(  1,  2 ) * C(  2,  1 )
                Z(  3,  5 ) = B(  1,  2) * A(  3,  1 ) + D(  1,  2 ) * C(  3,  1 )
                Z(  4,  5 ) = B(  1,  2) * A(  4,  1 ) + D(  1,  2 ) * C(  4,  1 )
                Z(  1,  6 ) = B(  1,  2) * A(  1,  2 ) + D(  1,  2 ) * C(  1,  2 )
                Z(  2,  6 ) = B(  1,  2) * A(  2,  2 ) + D(  1,  2 ) * C(  2,  2 )
                Z(  3,  6 ) = B(  1,  2) * A(  3,  2 ) + D(  1,  2 ) * C(  3,  2 )
                Z(  4,  6 ) = B(  1,  2) * A(  4,  2 ) + D(  1,  2 ) * C(  4,  2 )
                Z(  1,  7 ) = B(  1,  2) * A(  1,  3 ) + D(  1,  2 ) * C(  1,  3 )
                Z(  2,  7 ) = B(  1,  2) * A(  2,  3 ) + D(  1,  2 ) * C(  2,  3 )
                Z(  3,  7 ) = B(  1,  2) * A(  3,  3 ) + D(  1,  2 ) * C(  3,  3 )
                Z(  4,  7 ) = B(  1,  2) * A(  4,  3 ) + D(  1,  2 ) * C(  4,  3 )
                Z(  1,  8 ) = B(  1,  2) * A(  1,  4 ) + D(  1,  2 ) * C(  1,  4 )
                Z(  2,  8 ) = B(  1,  2) * A(  2,  4 ) + D(  1,  2 ) * C(  2,  4 )
                Z(  3,  8 ) = B(  1,  2) * A(  3,  4 ) + D(  1,  2 ) * C(  3,  4 )
                Z(  4,  8 ) = B(  1,  2) * A(  4,  4 ) + D(  1,  2 ) * C(  4,  4 )
                Z(  5,  5 ) = B(  2,  2) * A(  1,  1 ) + D(  2,  2 ) * C(  1,  1 )
                Z(  6,  5 ) = B(  2,  2) * A(  2,  1 ) + D(  2,  2 ) * C(  2,  1 )
                Z(  7,  5 ) = B(  2,  2) * A(  3,  1 ) + D(  2,  2 ) * C(  3,  1 )
                Z(  8,  5 ) = B(  2,  2) * A(  4,  1 ) + D(  2,  2 ) * C(  4,  1 )
                Z(  5,  6 ) = B(  2,  2) * A(  1,  2 ) + D(  2,  2 ) * C(  1,  2 )
                Z(  6,  6 ) = B(  2,  2) * A(  2,  2 ) + D(  2,  2 ) * C(  2,  2 )
                Z(  7,  6 ) = B(  2,  2) * A(  3,  2 ) + D(  2,  2 ) * C(  3,  2 )
                Z(  8,  6 ) = B(  2,  2) * A(  4,  2 ) + D(  2,  2 ) * C(  4,  2 )
                Z(  5,  7 ) = B(  2,  2) * A(  1,  3 ) + D(  2,  2 ) * C(  1,  3 )
                Z(  6,  7 ) = B(  2,  2) * A(  2,  3 ) + D(  2,  2 ) * C(  2,  3 )
                Z(  7,  7 ) = B(  2,  2) * A(  3,  3 ) + D(  2,  2 ) * C(  3,  3 )
                Z(  8,  7 ) = B(  2,  2) * A(  4,  3 ) + D(  2,  2 ) * C(  4,  3 )
                Z(  5,  8 ) = B(  2,  2) * A(  1,  4 ) + D(  2,  2 ) * C(  1,  4 )
                Z(  6,  8 ) = B(  2,  2) * A(  2,  4 ) + D(  2,  2 ) * C(  2,  4 )
                Z(  7,  8 ) = B(  2,  2) * A(  3,  4 ) + D(  2,  2 ) * C(  3,  4 )
                Z(  8,  8 ) = B(  2,  2) * A(  4,  4 ) + D(  2,  2 ) * C(  4,  4 )
                Z(  1,  9 ) = B(  1,  3) * A(  1,  1 ) + D(  1,  3 ) * C(  1,  1 )
                Z(  2,  9 ) = B(  1,  3) * A(  2,  1 ) + D(  1,  3 ) * C(  2,  1 )
                Z(  3,  9 ) = B(  1,  3) * A(  3,  1 ) + D(  1,  3 ) * C(  3,  1 )
                Z(  4,  9 ) = B(  1,  3) * A(  4,  1 ) + D(  1,  3 ) * C(  4,  1 )
                Z(  1, 10 ) = B(  1,  3) * A(  1,  2 ) + D(  1,  3 ) * C(  1,  2 )
                Z(  2, 10 ) = B(  1,  3) * A(  2,  2 ) + D(  1,  3 ) * C(  2,  2 )
                Z(  3, 10 ) = B(  1,  3) * A(  3,  2 ) + D(  1,  3 ) * C(  3,  2 )
                Z(  4, 10 ) = B(  1,  3) * A(  4,  2 ) + D(  1,  3 ) * C(  4,  2 )
                Z(  1, 11 ) = B(  1,  3) * A(  1,  3 ) + D(  1,  3 ) * C(  1,  3 )
                Z(  2, 11 ) = B(  1,  3) * A(  2,  3 ) + D(  1,  3 ) * C(  2,  3 )
                Z(  3, 11 ) = B(  1,  3) * A(  3,  3 ) + D(  1,  3 ) * C(  3,  3 )
                Z(  4, 11 ) = B(  1,  3) * A(  4,  3 ) + D(  1,  3 ) * C(  4,  3 )
                Z(  1, 12 ) = B(  1,  3) * A(  1,  4 ) + D(  1,  3 ) * C(  1,  4 )
                Z(  2, 12 ) = B(  1,  3) * A(  2,  4 ) + D(  1,  3 ) * C(  2,  4 )
                Z(  3, 12 ) = B(  1,  3) * A(  3,  4 ) + D(  1,  3 ) * C(  3,  4 )
                Z(  4, 12 ) = B(  1,  3) * A(  4,  4 ) + D(  1,  3 ) * C(  4,  4 )
                Z(  5,  9 ) = B(  2,  3) * A(  1,  1 ) + D(  2,  3 ) * C(  1,  1 )
                Z(  6,  9 ) = B(  2,  3) * A(  2,  1 ) + D(  2,  3 ) * C(  2,  1 )
                Z(  7,  9 ) = B(  2,  3) * A(  3,  1 ) + D(  2,  3 ) * C(  3,  1 )
                Z(  8,  9 ) = B(  2,  3) * A(  4,  1 ) + D(  2,  3 ) * C(  4,  1 )
                Z(  5, 10 ) = B(  2,  3) * A(  1,  2 ) + D(  2,  3 ) * C(  1,  2 )
                Z(  6, 10 ) = B(  2,  3) * A(  2,  2 ) + D(  2,  3 ) * C(  2,  2 )
                Z(  7, 10 ) = B(  2,  3) * A(  3,  2 ) + D(  2,  3 ) * C(  3,  2 )
                Z(  8, 10 ) = B(  2,  3) * A(  4,  2 ) + D(  2,  3 ) * C(  4,  2 )
                Z(  5, 11 ) = B(  2,  3) * A(  1,  3 ) + D(  2,  3 ) * C(  1,  3 )
                Z(  6, 11 ) = B(  2,  3) * A(  2,  3 ) + D(  2,  3 ) * C(  2,  3 )
                Z(  7, 11 ) = B(  2,  3) * A(  3,  3 ) + D(  2,  3 ) * C(  3,  3 )
                Z(  8, 11 ) = B(  2,  3) * A(  4,  3 ) + D(  2,  3 ) * C(  4,  3 )
                Z(  5, 12 ) = B(  2,  3) * A(  1,  4 ) + D(  2,  3 ) * C(  1,  4 )
                Z(  6, 12 ) = B(  2,  3) * A(  2,  4 ) + D(  2,  3 ) * C(  2,  4 )
                Z(  7, 12 ) = B(  2,  3) * A(  3,  4 ) + D(  2,  3 ) * C(  3,  4 )
                Z(  8, 12 ) = B(  2,  3) * A(  4,  4 ) + D(  2,  3 ) * C(  4,  4 )
                Z(  9,  9 ) = B(  3,  3) * A(  1,  1 ) + D(  3,  3 ) * C(  1,  1 )
                Z( 10,  9 ) = B(  3,  3) * A(  2,  1 ) + D(  3,  3 ) * C(  2,  1 )
                Z( 11,  9 ) = B(  3,  3) * A(  3,  1 ) + D(  3,  3 ) * C(  3,  1 )
                Z( 12,  9 ) = B(  3,  3) * A(  4,  1 ) + D(  3,  3 ) * C(  4,  1 )
                Z(  9, 10 ) = B(  3,  3) * A(  1,  2 ) + D(  3,  3 ) * C(  1,  2 )
                Z( 10, 10 ) = B(  3,  3) * A(  2,  2 ) + D(  3,  3 ) * C(  2,  2 )
                Z( 11, 10 ) = B(  3,  3) * A(  3,  2 ) + D(  3,  3 ) * C(  3,  2 )
                Z( 12, 10 ) = B(  3,  3) * A(  4,  2 ) + D(  3,  3 ) * C(  4,  2 )
                Z(  9, 11 ) = B(  3,  3) * A(  1,  3 ) + D(  3,  3 ) * C(  1,  3 )
                Z( 10, 11 ) = B(  3,  3) * A(  2,  3 ) + D(  3,  3 ) * C(  2,  3 )
                Z( 11, 11 ) = B(  3,  3) * A(  3,  3 ) + D(  3,  3 ) * C(  3,  3 )
                Z( 12, 11 ) = B(  3,  3) * A(  4,  3 ) + D(  3,  3 ) * C(  4,  3 )
                Z(  9, 12 ) = B(  3,  3) * A(  1,  4 ) + D(  3,  3 ) * C(  1,  4 )
                Z( 10, 12 ) = B(  3,  3) * A(  2,  4 ) + D(  3,  3 ) * C(  2,  4 )
                Z( 11, 12 ) = B(  3,  3) * A(  3,  4 ) + D(  3,  3 ) * C(  3,  4 )
                Z( 12, 12 ) = B(  3,  3) * A(  4,  4 ) + D(  3,  3 ) * C(  4,  4 )
                Z(  1, 13 ) = B(  1,  4) * A(  1,  1 ) + D(  1,  4 ) * C(  1,  1 )
                Z(  2, 13 ) = B(  1,  4) * A(  2,  1 ) + D(  1,  4 ) * C(  2,  1 )
                Z(  3, 13 ) = B(  1,  4) * A(  3,  1 ) + D(  1,  4 ) * C(  3,  1 )
                Z(  4, 13 ) = B(  1,  4) * A(  4,  1 ) + D(  1,  4 ) * C(  4,  1 )
                Z(  1, 14 ) = B(  1,  4) * A(  1,  2 ) + D(  1,  4 ) * C(  1,  2 )
                Z(  2, 14 ) = B(  1,  4) * A(  2,  2 ) + D(  1,  4 ) * C(  2,  2 )
                Z(  3, 14 ) = B(  1,  4) * A(  3,  2 ) + D(  1,  4 ) * C(  3,  2 )
                Z(  4, 14 ) = B(  1,  4) * A(  4,  2 ) + D(  1,  4 ) * C(  4,  2 )
                Z(  1, 15 ) = B(  1,  4) * A(  1,  3 ) + D(  1,  4 ) * C(  1,  3 )
                Z(  2, 15 ) = B(  1,  4) * A(  2,  3 ) + D(  1,  4 ) * C(  2,  3 )
                Z(  3, 15 ) = B(  1,  4) * A(  3,  3 ) + D(  1,  4 ) * C(  3,  3 )
                Z(  4, 15 ) = B(  1,  4) * A(  4,  3 ) + D(  1,  4 ) * C(  4,  3 )
                Z(  1, 16 ) = B(  1,  4) * A(  1,  4 ) + D(  1,  4 ) * C(  1,  4 )
                Z(  2, 16 ) = B(  1,  4) * A(  2,  4 ) + D(  1,  4 ) * C(  2,  4 )
                Z(  3, 16 ) = B(  1,  4) * A(  3,  4 ) + D(  1,  4 ) * C(  3,  4 )
                Z(  4, 16 ) = B(  1,  4) * A(  4,  4 ) + D(  1,  4 ) * C(  4,  4 )
                Z(  5, 13 ) = B(  2,  4) * A(  1,  1 ) + D(  2,  4 ) * C(  1,  1 )
                Z(  6, 13 ) = B(  2,  4) * A(  2,  1 ) + D(  2,  4 ) * C(  2,  1 )
                Z(  7, 13 ) = B(  2,  4) * A(  3,  1 ) + D(  2,  4 ) * C(  3,  1 )
                Z(  8, 13 ) = B(  2,  4) * A(  4,  1 ) + D(  2,  4 ) * C(  4,  1 )
                Z(  5, 14 ) = B(  2,  4) * A(  1,  2 ) + D(  2,  4 ) * C(  1,  2 )
                Z(  6, 14 ) = B(  2,  4) * A(  2,  2 ) + D(  2,  4 ) * C(  2,  2 )
                Z(  7, 14 ) = B(  2,  4) * A(  3,  2 ) + D(  2,  4 ) * C(  3,  2 )
                Z(  8, 14 ) = B(  2,  4) * A(  4,  2 ) + D(  2,  4 ) * C(  4,  2 )
                Z(  5, 15 ) = B(  2,  4) * A(  1,  3 ) + D(  2,  4 ) * C(  1,  3 )
                Z(  6, 15 ) = B(  2,  4) * A(  2,  3 ) + D(  2,  4 ) * C(  2,  3 )
                Z(  7, 15 ) = B(  2,  4) * A(  3,  3 ) + D(  2,  4 ) * C(  3,  3 )
                Z(  8, 15 ) = B(  2,  4) * A(  4,  3 ) + D(  2,  4 ) * C(  4,  3 )
                Z(  5, 16 ) = B(  2,  4) * A(  1,  4 ) + D(  2,  4 ) * C(  1,  4 )
                Z(  6, 16 ) = B(  2,  4) * A(  2,  4 ) + D(  2,  4 ) * C(  2,  4 )
                Z(  7, 16 ) = B(  2,  4) * A(  3,  4 ) + D(  2,  4 ) * C(  3,  4 )
                Z(  8, 16 ) = B(  2,  4) * A(  4,  4 ) + D(  2,  4 ) * C(  4,  4 )
                Z(  9, 13 ) = B(  3,  4) * A(  1,  1 ) + D(  3,  4 ) * C(  1,  1 )
                Z( 10, 13 ) = B(  3,  4) * A(  2,  1 ) + D(  3,  4 ) * C(  2,  1 )
                Z( 11, 13 ) = B(  3,  4) * A(  3,  1 ) + D(  3,  4 ) * C(  3,  1 )
                Z( 12, 13 ) = B(  3,  4) * A(  4,  1 ) + D(  3,  4 ) * C(  4,  1 )
                Z(  9, 14 ) = B(  3,  4) * A(  1,  2 ) + D(  3,  4 ) * C(  1,  2 )
                Z( 10, 14 ) = B(  3,  4) * A(  2,  2 ) + D(  3,  4 ) * C(  2,  2 )
                Z( 11, 14 ) = B(  3,  4) * A(  3,  2 ) + D(  3,  4 ) * C(  3,  2 )
                Z( 12, 14 ) = B(  3,  4) * A(  4,  2 ) + D(  3,  4 ) * C(  4,  2 )
                Z(  9, 15 ) = B(  3,  4) * A(  1,  3 ) + D(  3,  4 ) * C(  1,  3 )
                Z( 10, 15 ) = B(  3,  4) * A(  2,  3 ) + D(  3,  4 ) * C(  2,  3 )
                Z( 11, 15 ) = B(  3,  4) * A(  3,  3 ) + D(  3,  4 ) * C(  3,  3 )
                Z( 12, 15 ) = B(  3,  4) * A(  4,  3 ) + D(  3,  4 ) * C(  4,  3 )
                Z(  9, 16 ) = B(  3,  4) * A(  1,  4 ) + D(  3,  4 ) * C(  1,  4 )
                Z( 10, 16 ) = B(  3,  4) * A(  2,  4 ) + D(  3,  4 ) * C(  2,  4 )
                Z( 11, 16 ) = B(  3,  4) * A(  3,  4 ) + D(  3,  4 ) * C(  3,  4 )
                Z( 12, 16 ) = B(  3,  4) * A(  4,  4 ) + D(  3,  4 ) * C(  4,  4 )
                Z( 13, 13 ) = B(  4,  4) * A(  1,  1 ) + D(  4,  4 ) * C(  1,  1 )
                Z( 14, 13 ) = B(  4,  4) * A(  2,  1 ) + D(  4,  4 ) * C(  2,  1 )
                Z( 15, 13 ) = B(  4,  4) * A(  3,  1 ) + D(  4,  4 ) * C(  3,  1 )
                Z( 16, 13 ) = B(  4,  4) * A(  4,  1 ) + D(  4,  4 ) * C(  4,  1 )
                Z( 13, 14 ) = B(  4,  4) * A(  1,  2 ) + D(  4,  4 ) * C(  1,  2 )
                Z( 14, 14 ) = B(  4,  4) * A(  2,  2 ) + D(  4,  4 ) * C(  2,  2 )
                Z( 15, 14 ) = B(  4,  4) * A(  3,  2 ) + D(  4,  4 ) * C(  3,  2 )
                Z( 16, 14 ) = B(  4,  4) * A(  4,  2 ) + D(  4,  4 ) * C(  4,  2 )
                Z( 13, 15 ) = B(  4,  4) * A(  1,  3 ) + D(  4,  4 ) * C(  1,  3 )
                Z( 14, 15 ) = B(  4,  4) * A(  2,  3 ) + D(  4,  4 ) * C(  2,  3 )
                Z( 15, 15 ) = B(  4,  4) * A(  3,  3 ) + D(  4,  4 ) * C(  3,  3 )
                Z( 16, 15 ) = B(  4,  4) * A(  4,  3 ) + D(  4,  4 ) * C(  4,  3 )
                Z( 13, 16 ) = B(  4,  4) * A(  1,  4 ) + D(  4,  4 ) * C(  1,  4 )
                Z( 14, 16 ) = B(  4,  4) * A(  2,  4 ) + D(  4,  4 ) * C(  2,  4 )
                Z( 15, 16 ) = B(  4,  4) * A(  3,  4 ) + D(  4,  4 ) * C(  3,  4 )
                Z( 16, 16 ) = B(  4,  4) * A(  4,  4 ) + D(  4,  4 ) * C(  4,  4 )

                IF ( B (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  5,  1 ) = Z(  5,  1 ) + B(  2,  1) * A(  1,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + B(  2,  1) * A(  2,  1 )
                    Z(  7,  1 ) = Z(  7,  1 ) + B(  2,  1) * A(  3,  1 )
                    Z(  8,  1 ) = Z(  8,  1 ) + B(  2,  1) * A(  4,  1 )
                    Z(  5,  2 ) = Z(  5,  2 ) + B(  2,  1) * A(  1,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + B(  2,  1) * A(  2,  2 )
                    Z(  7,  2 ) = Z(  7,  2 ) + B(  2,  1) * A(  3,  2 )
                    Z(  8,  2 ) = Z(  8,  2 ) + B(  2,  1) * A(  4,  2 )
                    Z(  5,  3 ) = Z(  5,  3 ) + B(  2,  1) * A(  1,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + B(  2,  1) * A(  2,  3 )
                    Z(  7,  3 ) = Z(  7,  3 ) + B(  2,  1) * A(  3,  3 )
                    Z(  8,  3 ) = Z(  8,  3 ) + B(  2,  1) * A(  4,  3 )
                    Z(  5,  4 ) = Z(  5,  4 ) + B(  2,  1) * A(  1,  4 )
                    Z(  6,  4 ) = Z(  6,  4 ) + B(  2,  1) * A(  2,  4 )
                    Z(  7,  4 ) = Z(  7,  4 ) + B(  2,  1) * A(  3,  4 )
                    Z(  8,  4 ) = Z(  8,  4 ) + B(  2,  1) * A(  4,  4 )
                END IF
                IF ( B (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  9,  5 ) = Z(  9,  5 ) + B(  3,  2) * A(  1,  1 )
                    Z( 10,  5 ) = Z( 10,  5 ) + B(  3,  2) * A(  2,  1 )
                    Z( 11,  5 ) = Z( 11,  5 ) + B(  3,  2) * A(  3,  1 )
                    Z( 12,  5 ) = Z( 12,  5 ) + B(  3,  2) * A(  4,  1 )
                    Z(  9,  6 ) = Z(  9,  6 ) + B(  3,  2) * A(  1,  2 )
                    Z( 10,  6 ) = Z( 10,  6 ) + B(  3,  2) * A(  2,  2 )
                    Z( 11,  6 ) = Z( 11,  6 ) + B(  3,  2) * A(  3,  2 )
                    Z( 12,  6 ) = Z( 12,  6 ) + B(  3,  2) * A(  4,  2 )
                    Z(  9,  7 ) = Z(  9,  7 ) + B(  3,  2) * A(  1,  3 )
                    Z( 10,  7 ) = Z( 10,  7 ) + B(  3,  2) * A(  2,  3 )
                    Z( 11,  7 ) = Z( 11,  7 ) + B(  3,  2) * A(  3,  3 )
                    Z( 12,  7 ) = Z( 12,  7 ) + B(  3,  2) * A(  4,  3 )
                    Z(  9,  8 ) = Z(  9,  8 ) + B(  3,  2) * A(  1,  4 )
                    Z( 10,  8 ) = Z( 10,  8 ) + B(  3,  2) * A(  2,  4 )
                    Z( 11,  8 ) = Z( 11,  8 ) + B(  3,  2) * A(  3,  4 )
                    Z( 12,  8 ) = Z( 12,  8 ) + B(  3,  2) * A(  4,  4 )
                END IF
                IF ( B (  4,  3) .NE. 0.0d0 ) THEN
                    Z( 13,  9 ) = Z( 13,  9 ) + B(  4,  3) * A(  1,  1 )
                    Z( 14,  9 ) = Z( 14,  9 ) + B(  4,  3) * A(  2,  1 )
                    Z( 15,  9 ) = Z( 15,  9 ) + B(  4,  3) * A(  3,  1 )
                    Z( 16,  9 ) = Z( 16,  9 ) + B(  4,  3) * A(  4,  1 )
                    Z( 13, 10 ) = Z( 13, 10 ) + B(  4,  3) * A(  1,  2 )
                    Z( 14, 10 ) = Z( 14, 10 ) + B(  4,  3) * A(  2,  2 )
                    Z( 15, 10 ) = Z( 15, 10 ) + B(  4,  3) * A(  3,  2 )
                    Z( 16, 10 ) = Z( 16, 10 ) + B(  4,  3) * A(  4,  2 )
                    Z( 13, 11 ) = Z( 13, 11 ) + B(  4,  3) * A(  1,  3 )
                    Z( 14, 11 ) = Z( 14, 11 ) + B(  4,  3) * A(  2,  3 )
                    Z( 15, 11 ) = Z( 15, 11 ) + B(  4,  3) * A(  3,  3 )
                    Z( 16, 11 ) = Z( 16, 11 ) + B(  4,  3) * A(  4,  3 )
                    Z( 13, 12 ) = Z( 13, 12 ) + B(  4,  3) * A(  1,  4 )
                    Z( 14, 12 ) = Z( 14, 12 ) + B(  4,  3) * A(  2,  4 )
                    Z( 15, 12 ) = Z( 15, 12 ) + B(  4,  3) * A(  3,  4 )
                    Z( 16, 12 ) = Z( 16, 12 ) + B(  4,  3) * A(  4,  4 )
                END IF

                IF ( D (  2,  1) .NE. 0.0d0 ) THEN
                    Z(  5,  1 ) = Z(  5,  1 ) + D(  2,  1) * C(  1,  1 )
                    Z(  6,  1 ) = Z(  6,  1 ) + D(  2,  1) * C(  2,  1 )
                    Z(  7,  1 ) = Z(  7,  1 ) + D(  2,  1) * C(  3,  1 )
                    Z(  8,  1 ) = Z(  8,  1 ) + D(  2,  1) * C(  4,  1 )
                    Z(  5,  2 ) = Z(  5,  2 ) + D(  2,  1) * C(  1,  2 )
                    Z(  6,  2 ) = Z(  6,  2 ) + D(  2,  1) * C(  2,  2 )
                    Z(  7,  2 ) = Z(  7,  2 ) + D(  2,  1) * C(  3,  2 )
                    Z(  8,  2 ) = Z(  8,  2 ) + D(  2,  1) * C(  4,  2 )
                    Z(  5,  3 ) = Z(  5,  3 ) + D(  2,  1) * C(  1,  3 )
                    Z(  6,  3 ) = Z(  6,  3 ) + D(  2,  1) * C(  2,  3 )
                    Z(  7,  3 ) = Z(  7,  3 ) + D(  2,  1) * C(  3,  3 )
                    Z(  8,  3 ) = Z(  8,  3 ) + D(  2,  1) * C(  4,  3 )
                    Z(  5,  4 ) = Z(  5,  4 ) + D(  2,  1) * C(  1,  4 )
                    Z(  6,  4 ) = Z(  6,  4 ) + D(  2,  1) * C(  2,  4 )
                    Z(  7,  4 ) = Z(  7,  4 ) + D(  2,  1) * C(  3,  4 )
                    Z(  8,  4 ) = Z(  8,  4 ) + D(  2,  1) * C(  4,  4 )
                END IF
                IF ( D (  3,  2) .NE. 0.0d0 ) THEN
                    Z(  9,  5 ) = Z(  9,  5 ) + D(  3,  2) * C(  1,  1 )
                    Z( 10,  5 ) = Z( 10,  5 ) + D(  3,  2) * C(  2,  1 )
                    Z( 11,  5 ) = Z( 11,  5 ) + D(  3,  2) * C(  3,  1 )
                    Z( 12,  5 ) = Z( 12,  5 ) + D(  3,  2) * C(  4,  1 )
                    Z(  9,  6 ) = Z(  9,  6 ) + D(  3,  2) * C(  1,  2 )
                    Z( 10,  6 ) = Z( 10,  6 ) + D(  3,  2) * C(  2,  2 )
                    Z( 11,  6 ) = Z( 11,  6 ) + D(  3,  2) * C(  3,  2 )
                    Z( 12,  6 ) = Z( 12,  6 ) + D(  3,  2) * C(  4,  2 )
                    Z(  9,  7 ) = Z(  9,  7 ) + D(  3,  2) * C(  1,  3 )
                    Z( 10,  7 ) = Z( 10,  7 ) + D(  3,  2) * C(  2,  3 )
                    Z( 11,  7 ) = Z( 11,  7 ) + D(  3,  2) * C(  3,  3 )
                    Z( 12,  7 ) = Z( 12,  7 ) + D(  3,  2) * C(  4,  3 )
                    Z(  9,  8 ) = Z(  9,  8 ) + D(  3,  2) * C(  1,  4 )
                    Z( 10,  8 ) = Z( 10,  8 ) + D(  3,  2) * C(  2,  4 )
                    Z( 11,  8 ) = Z( 11,  8 ) + D(  3,  2) * C(  3,  4 )
                    Z( 12,  8 ) = Z( 12,  8 ) + D(  3,  2) * C(  4,  4 )
                END IF
                IF ( D (  4,  3) .NE. 0.0d0 ) THEN
                    Z( 13,  9 ) = Z( 13,  9 ) + D(  4,  3) * C(  1,  1 )
                    Z( 14,  9 ) = Z( 14,  9 ) + D(  4,  3) * C(  2,  1 )
                    Z( 15,  9 ) = Z( 15,  9 ) + D(  4,  3) * C(  3,  1 )
                    Z( 16,  9 ) = Z( 16,  9 ) + D(  4,  3) * C(  4,  1 )
                    Z( 13, 10 ) = Z( 13, 10 ) + D(  4,  3) * C(  1,  2 )
                    Z( 14, 10 ) = Z( 14, 10 ) + D(  4,  3) * C(  2,  2 )
                    Z( 15, 10 ) = Z( 15, 10 ) + D(  4,  3) * C(  3,  2 )
                    Z( 16, 10 ) = Z( 16, 10 ) + D(  4,  3) * C(  4,  2 )
                    Z( 13, 11 ) = Z( 13, 11 ) + D(  4,  3) * C(  1,  3 )
                    Z( 14, 11 ) = Z( 14, 11 ) + D(  4,  3) * C(  2,  3 )
                    Z( 15, 11 ) = Z( 15, 11 ) + D(  4,  3) * C(  3,  3 )
                    Z( 16, 11 ) = Z( 16, 11 ) + D(  4,  3) * C(  4,  3 )
                    Z( 13, 12 ) = Z( 13, 12 ) + D(  4,  3) * C(  1,  4 )
                    Z( 14, 12 ) = Z( 14, 12 ) + D(  4,  3) * C(  2,  4 )
                    Z( 15, 12 ) = Z( 15, 12 ) + D(  4,  3) * C(  3,  4 )
                    Z( 16, 12 ) = Z( 16, 12 ) + D(  4,  3) * C(  4,  4 )
                END IF
            END IF
        END IF
    ENDIF
    ! Get RHS
    DO J = 1, N
        DO I = 1, M
            RHS( I + (J-1) * M) = X(I, J)
        END DO
    END DO

    ! Solve
    CALL DGETC2X(DIMZ,Z,LDZ,IPIV,JPIV,IINFO)
    IF (IINFO .NE. 0 ) THEN
        INFO = 4
        RETURN
    END IF
    CALL DGESC2(DIMZ,Z,LDZ,RHS,IPIV,JPIV,LSCAL)
    IF (SCAL .NE. ONE) THEN
        INFO = 1
    ENDIF
    SCAL = LSCAL

    ! Restore RHS
    DO J = 1, N
        DO I = 1, M
            X(I,J) = RHS( I + (J-1) * M)
        END DO
    END DO

    RETURN
END SUBROUTINE


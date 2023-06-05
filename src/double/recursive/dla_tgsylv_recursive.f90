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


!> \brief Level-3 Recursive Blocked Solver for the generalized Sylvester equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TGSYLV_RECURSIVE ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           DOUBLE PRECISION SGN, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TGSYLV_RECURSIVE solves a generalized Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side E and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by DGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is a rewrite of RECSY with some optimizations.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A and C:
!>          == 'N':  op1(A) = A, op1(C) = C (No transpose for A and C)
!>          == 'T':  op1(A) = A**T, op1(C) = C **T (Transpose A and C)
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to B and D:
!>          == 'N':  op2(B) = B, op2(D) = D (No transpose for B and D)
!>          == 'T':  op2(B) = B**T, op2(D) = D **T (Transpose B and D)
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
!>          The order of the matrices A and C.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices B and D.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          The matrix A must be (quasi-) upper triangular.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The matrix B must be (quasi-) upper triangular. If the matrix D is already
!>          quasi-upper triangular the matrix B must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,M)
!>          The matrix C must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (LDD,N)
!>          The matrix D must be (quasi-) upper triangular. If the matrix B is already
!>          quasi-upper triangular the matrix D must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D.  LDD >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution X of Equation (1) or (2)
!>          as selected by TRANSA, TRANSB, and SGN.
!>          Right hand side Y and the solution X are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,M).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension M*N
!>          Workspace for the algorithm
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          == 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  The equation is not solved correctly. One of the arising inner
!>                system got singular.
!> \endverbatim
!>
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date June 2023
!> \ingroup dbltgsylv
!

RECURSIVE SUBROUTINE DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, &
        & SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Input
    CHARACTER TRANSA(1), TRANSB(1)
    INTEGER M, N, LDA, LDB, LDC, LDD, LDX , INFO
    DOUBLE PRECISION SGN
    DOUBLE PRECISION A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), X(LDX, *)
    DOUBLE PRECISION SCALE
    DOUBLE PRECISION WORK(*)

    ! External Functions
    EXTERNAL DLA_GEAXB
    EXTERNAL DLA_TGSYLV_KERNEL_44NN
    EXTERNAL DLA_TGSYLV_KERNEL_44NT
    EXTERNAL DLA_TGSYLV_KERNEL_44TN
    EXTERNAL DLA_TGSYLV_KERNEL_44TT
    EXTERNAL DLA_TGSYLV_L2
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGSYLV_L2_REORDER
    EXTERNAL DLA_TGSYLV_L2_UNOPT
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    ! Local Variables
    INTEGER ININFO
    CHARACTER TYPEB, TYPED
    LOGICAL BTRANSA, BTRANSB
    INTEGER M1, M2, N1, N2, IINFO, I
    DOUBLE PRECISION ZERO, ONE
    PARAMETER(ZERO = 0.0d0, ONE = 1.0d0)
    DOUBLE PRECISION SCAL
    INTEGER BOUND
    ! INTEGER MB, NB, MAXNB, ISOLVER
    PARAMETER(BOUND = 4)

    ININFO = INFO
    INFO = 0
    BTRANSA = LSAME(TRANSA, 'T')
    BTRANSB = LSAME(TRANSB, 'T')

    ! CHECK Input
    IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'N')) THEN
        INFO = -1
    ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'N')) THEN
        INFO = -2
    ELSE IF ( SGN .NE. -ONE .AND. SGN .NE. ONE ) THEN
        INFO = -3
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -4
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -5
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -7
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -9
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -11
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -13
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -15
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_TGSYLV_RECURSIVE', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        INFO = MAX(1, M*N)
        RETURN
    END IF


    TYPEB = 'T'
    TYPED = 'T'
    DO I = 1, N-1
        IF ( B(I+1,I) .NE. ZERO) TYPEB = 'Q'
        IF ( D(I+1,I) .NE. ZERO) TYPED = 'Q'
    END DO

    ! Quick Return
    IF (M .LE. 0 .OR. N.LE.0 ) THEN
        RETURN
    ENDIF

    SCALE = ONE

    IF (M .LE. BOUND .AND. N .LE. BOUND) THEN
        ! IF ( BOUND .LE. 4 ) THEN
            IINFO = 0
            IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB) THEN
                !(N,N)
                CALL DLA_TGSYLV_KERNEL_44NN(SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE, IINFO)
            ELSE IF ( .NOT. BTRANSA .AND. BTRANSB) THEN
                !(N,T)
                CALL DLA_TGSYLV_KERNEL_44NT(SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE,IINFO)
            ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
                !(T,N)
                CALL DLA_TGSYLV_KERNEL_44TN(SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE,IINFO)

            ELSE IF ( BTRANSA .AND. BTRANSB) THEN
                !(T,T)
                CALL DLA_TGSYLV_KERNEL_44TT(SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, SCALE,IINFO)

            ENDIF
            INFO = IINFO
    !     ELSE
    !         SCAL = ONE;
    !         SELECT CASE(ISOLVER)
    !         CASE(1)
    !             IF ( MAXNB .LE. 32 ) THEN
    !                 CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                     & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !             ELSE IF ( MAXNB .LE. 64) THEN
    !                 CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                     & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !             ELSE IF ( MAXNB .LE. 96) THEN
    !                 CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                     & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !             ELSE IF ( MAXNB .LE. 128) THEN
    !                 CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                     & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !             ELSE
    !                 CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                     & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !             END IF
    !         CASE(2)
    !             IF ( M .GT. 128 .OR. N .GT. 128) THEN
    !                 CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                     & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !             ELSE
    !                 CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                     & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !             END IF
    !         CASE(3)
    !             CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                 & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !         CASE(4)
    !             CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                 & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !         CASE(5)
    !             CALL DLA_TGSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                 & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !         CASE(6)
    !             ! Do not start the recursive algorithm again.
    !             CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
    !                 & C, LDC, D, LDD,  X, LDX, SCAL, WORK, IINFO)
    !         END SELECT
    !         IF ( IINFO.NE. 0 ) THEN
    !             INFO = IINFO
    !         END IF
    !         IF ( SCAL .NE. ONE) THEN
    !             SCALE = SCAL
    !         END IF
    !
    !     END IF
        RETURN
    END IF


    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB ) THEN
        ! (N,N)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('N', 'N', SGN, M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D,LDD, &
                & X(M1+1,1), LDX, SCALE, WORK, INFO )

            CALL DLA_GEAXB('N', 'N', 'F', TYPEB, M1, M2, N, N, -ONE, A(1,M1+1), LDA, X(M1+1,1), LDX, B, LDB, SCALE, X(1,1), &
                & LDX, WORK)
            CALL DLA_GEAXB('N', 'N', 'F', TYPED, M1, M2, N, N, -SGN, C(1,M1+1), LDC, X(M1+1,1), LDX, D, LDD, ONE, X(1,1), &
                & LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('N', 'N', SGN, M1, N, A(1,1), LDA, B, LDB, C(1,1), LDC, D, LDD, X(1,1), LDX, &
                & SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(M1+1:M, 1:N) = SCAL*X(M1+1:M, 1:N)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('N', 'N', SGN, M, N1, A, LDA, B(1,1), LDB, C, LDC, D(1,1), LDD, X(1,1), LDX, &
                & SCALE, WORK, INFO)
            CALL DLA_GEAXB('N', 'N', 'Q', 'F', M, M, N1, N2, -ONE, A, LDA, X(1,1), LDX, B(1,N1+1), LDB, SCALE, X(1,N1+1), &
                & LDX, WORK)
            CALL DLA_GEAXB('N', 'N', 'T', 'F', M, M, N1, N2, -SGN, C, LDC, X(1,1), LDX, D(1,N1+1), LDD, ONE, X(1,N1+1), &
                & LDX, WORK)
            ! X2
            CALL DLA_TGSYLV_RECURSIVE('N', 'N', SGN, M, N2, A, LDA, B(N1+1,N1+1), LDB, C, LDC, D(N1+1,N1+1), LDD, &
                & X(1,N1+1), LDX, SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, 1:N1) = SCAL*X(1:M, 1:N1)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('N','N', SGN, M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, C(M1+1,M1+1), LDC, D(1,1), LDD, &
                & X(M1+1,1), LDX, SCALE, WORK, INFO)

            CALL DLA_GEAXB('N', 'N', 'F', TYPEB, M1, M2, N1, N1, -ONE, A(1,M1+1), LDA, X(M1+1,1), LDX, B(1,1), LDB, &
                & SCALE, X(1,1), LDX, WORK)
            CALL DLA_GEAXB('N', 'N', 'F', TYPED, M1, M2, N1, N1, -SGN, C(1,M1+1), LDC, X(M1+1,1), LDX, D(1,1), LDD, &
                & ONE, X(1,1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('N','N', SGN, M1, N1, A(1,1), LDA, B(1,1), LDB, C(1,1), LDC, D(1, 1), LDD, &
                & X(1,1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(M1+1:M, 1:N1) = SCAL * X(M1+1:M, 1:N1)
            ENDIF

            CALL DLA_GEAXB('N', 'N', 'Q', 'F', M2, M2, N1, N2, -ONE, A(M1+1, M1+1), LDA, X(M1+1,1), LDX, B(1,N1+1), LDB, &
                & SCALE, X(M1+1,N1+1), LDX, WORK)
            CALL DLA_GEAXB('N', 'N', 'T', 'F', M2, M2, N1, N2, -SGN, C(M1+1, M1+1), LDC, X(M1+1,1), LDX, D(1,N1+1), LDD, &
                & ONE, X(M1+1,N1+1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('N','N', SGN, M2, N2, A(M1+1, M1+1), LDA, B(N1+1,N1+1), LDB, C(M1+1, M1+1), LDC,&
                & D(N1+1,N1+1), LDD, X(M1+1, N1+1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M, 1:N1) = SCAL * X(1:M, 1:N1)
            ENDIF

            CALL DLA_GEAXB('N', 'N', 'Q', 'F', M1, M, N1, N2, -ONE, A(1,1), LDA, X(1,1),LDX, B(1,N1+1),LDB, &
                & SCALE,  X(1, N1+1), LDX, WORK)

            CALL DLA_GEAXB('N', 'N', 'F', TYPEB, M1, M2, N2,N2, -ONE, A(1,M1+1), LDA, X(M1+1,N1+1), LDX, B(N1+1, N1+1), LDB, &
                & ONE, X(1,N1+1), LDX, WORK)

            CALL DLA_GEAXB('N', 'N', 'T', 'F', M1, M, N1, N2, -SGN, C(1,1), LDC, X(1,1),LDX, D(1,N1+1),LDD, &
                & ONE, X(1, N1+1), LDX, WORK)
            CALL DLA_GEAXB('N', 'N', 'F', TYPED, M1, M2, N2,N2, -SGN, C(1,M1+1), LDC, X(M1+1,N1+1), LDX, D(N1+1, N1+1), LDD, &
                & ONE, X(1,N1+1), LDX, WORK)

            CALL DLA_TGSYLV_RECURSIVE('N', 'N',SGN,  M1, N2, A(1,1), LDA, B(N1+1, N1+1), LDB, C(1,1), LDC, D(N1+1,N1+1),&
                & LDD, X(1,N1+1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M, 1:N1) = SCAL * X(1:M, 1:N1)
                X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF
        END IF


    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N, T)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('N', 'T', SGN, M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D,LDD, &
                X(M1+1,1), LDX, SCALE, WORK, INFO )

            CALL DLA_GEAXB('N', 'T', 'F', TYPEB, M1, M2, N, N, -ONE, A(1,M1+1), LDA, X(M1+1,1), LDX, B, LDB, SCALE, X(1,1), &
                & LDX, WORK)
            CALL DLA_GEAXB('N', 'T', 'F', TYPED, M1, M2, N, N, -SGN, C(1,M1+1), LDC, X(M1+1,1), LDX, D, LDD, ONE, X(1,1), &
                & LDX, WORK)
            ! X22
            CALL DLA_TGSYLV_RECURSIVE('N', 'T', SGN, M1, N, A(1,1), LDA, B, LDB, C(1,1), LDC, D, LDD, X(1,1), LDX, &
                & SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(M1+1:M, 1:N) = SCAL*X(M1+1:M, 1:N)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('N', 'T', SGN, M, N2, A, LDA, B(N1+1, N1+1), LDB, C, LDC, D(N1+1,N1+1), LDD, &
                X(1,N1+1), LDX, SCALE, WORK, INFO)
            CALL DLA_GEAXB('N', 'T', 'Q', 'F', M, M, N2, N1, -ONE, A, LDA, X(1,N1+1), LDX, B(1,N1+1), LDB, SCALE, X(1,1), &
                & LDX, WORK)
            CALL DLA_GEAXB('N', 'T', 'T', 'F', M, M, N2, N1, -SGN, C, LDC, X(1,N1+1), LDX, D(1,N1+1), LDD, ONE, X(1,1), &
                & LDX, WORK)

            CALL DLA_TGSYLV_RECURSIVE('N', 'T',  SGN, M, N1, A, LDA, B(1,1), LDB, C, LDC, D(1,1), LDD, X(1,1), LDX, &
                & SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, N1+1:N) = SCAL*X(1:M, N1+1:N)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('N','T', SGN, M2, N2, A(M1+1,M1+1), LDA, B(N1+1,N1+1), LDB, C(M1+1,M1+1), LDC, &
                & D(N1+1,N1+1), LDD, X(M1+1,N1+1), LDX, SCALE, WORK, INFO)

            CALL DLA_GEAXB('N', 'T', 'Q', 'F', M2, M2, N2, N1, -ONE, A(M1+1,M1+1), LDA, X(M1+1,N1+1), LDX, B(1,N1+1), LDB, &
                 & SCALE, X(M1+1,1), LDX, WORK)
            CALL DLA_GEAXB('N', 'T', 'T', 'F', M2, M2, N2, N1, -SGN, C(M1+1,M1+1), LDC, X(M1+1,N1+1), LDX, D(1,N1+1), LDD, &
                 & ONE, X(M1+1,1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('N','T', SGN,  M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, C(M1+1,M1+1), LDC, D(1, 1), LDD, &
                 & X(M1+1,1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                 SCALE = SCAL * SCALE
                 X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF

            CALL DLA_GEAXB('N', 'T', 'F', TYPEB, M1, M2, N2, N2, -ONE, A(1, M1+1), LDA, X(M1+1,N1+1), LDX, B(N1+1,N1+1), LDB, &
                & SCALE, X(1,N1+1), LDX, WORK)
            CALL DLA_GEAXB('N', 'T', 'F', TYPED, M1, M2, N2, N2, -SGN, C(1, M1+1), LDC, X(M1+1,N1+1), LDX, D(N1+1,N1+1), LDD, &
                & ONE, X(1,N1+1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('N','T', SGN,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, C(1,1), LDC, D(N1+1,N1+1), LDD, &
                & X(1, N1+1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(M1+1, 1:N) = SCAL * X(M1+1, 1:N)
            ENDIF

            CALL DLA_GEAXB('N', 'T', 'F', TYPEB, M1, M2, N1, N1, -ONE, A(1,M1+1), LDA, X(M1+1,1),LDX, B(1,1),LDB, &
                & SCALE,  X(1, 1), LDX, WORK)
            CALL DLA_GEAXB('N', 'T', 'Q', 'F',   M1, M,  N2, N1, -ONE, A(1,1), LDA, X(1,N1+1), LDX, B(1, N1+1), LDB, &
                & ONE, X(1,1), LDX, WORK)
            CALL DLA_GEAXB('N', 'T', 'F', TYPED, M1, M2, N1, N1, -SGN, C(1,M1+1), LDC, X(M1+1,1),LDX, D(1,1),LDD, &
                & ONE,  X(1, 1), LDX, WORK)
            CALL DLA_GEAXB('N', 'T', 'T', 'F',   M1, M,  N2, N1, -SGN, C(1,1), LDC, X(1,N1+1), LDX, D(1, N1+1), LDD, &
                & ONE, X(1,1), LDX, WORK)


            CALL DLA_TGSYLV_RECURSIVE('N', 'T',SGN,  M1, N1, A(1,1), LDA, B(1,1), LDB, C(1,1), LDC, D(1,1),&
                & LDD, X(1,1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M, 1:N1) = SCAL * X(1:M, 1:N1)
                X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF
        END IF


    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
        ! (T, N)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('T', 'N', SGN,  M1, N, A(1,1), LDA, B, LDB, C(1,1), LDC, D,LDD, X(1,1), LDX, &
                & SCALE, WORK, INFO )

            CALL DLA_GEAXB('T', 'N', 'F', TYPEB, M2, M1, N, N, -ONE, A(1,M1+1), LDA, X(1,1), LDX, B, LDB, SCALE, X(M1+1,1), &
                & LDX, WORK)
            CALL DLA_GEAXB('T', 'N', 'F', TYPED, M2, M1, N, N, -SGN, C(1,M1+1), LDC, X(1,1), LDX, D, LDD, ONE, X(M1+1,1), &
                & LDX, WORK)

            CALL DLA_TGSYLV_RECURSIVE('T', 'N', SGN,  M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D, LDD, &
                & X(M1+1,1), LDX, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M1, 1:N) = SCAL*X(1:M1, 1:N)
                SCALE = SCAL * SCALE
            ENDIF

        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF ( B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('T', 'N', SGN,  M, N1, A, LDA, B(1,1), LDB, C, LDC, D(1,1), LDD, X(1,1), LDX, &
                & SCALE, WORK, INFO)
            CALL DLA_GEAXB('T', 'N', 'Q', 'F', M, M, N1, N2, -ONE, A, LDA, X(1,1), LDX, B(1,N1+1), LDB, SCALE, X(1,N1+1), &
                & LDX, WORK)
            CALL DLA_GEAXB('T', 'N', 'T', 'F', M, M, N1, N2, -SGN, C, LDC, X(1,1), LDX, D(1,N1+1), LDD, ONE, X(1,N1+1), &
                & LDX, WORK)

            CALL DLA_TGSYLV_RECURSIVE('T', 'N', SGN,  M, N2, A, LDA, B(N1+1,N1+1), LDB, C, LDC, D(N1+1,N1+1), LDD, &
                & X(1,N1+1), LDX, SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, 1:N1) = SCAL*X(1:M, 1:N1)
                SCALE = SCAL * SCALE
            ENDIF
        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF ( B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF


            CALL DLA_TGSYLV_RECURSIVE('T','N', SGN,  M1, N1, A(1,1), LDA, B(1,1), LDB, C(1,1), LDC, D(1,1), LDD, X(1,1), LDX, &
                & SCALE, WORK, INFO)

            CALL DLA_GEAXB('T', 'N', 'Q', 'F', M1, M1, N1, N2, -ONE, A(1,1), LDA, X(1,1), LDX, B(1,N1+1), LDB, &
                & SCALE, X(1,N1+1), LDX, WORK)
            CALL DLA_GEAXB('T', 'N', 'T', 'F', M1, M1, N1, N2, -SGN, C(1,1), LDC, X(1,1), LDX, D(1,N1+1), LDD, &
                & ONE, X(1,N1+1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('T','N', SGN,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, C(1,1), LDC, D(N1+1, N1+1), LDD, &
                & X(1,N1+1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N1) = SCAL * X(1:M1, 1:N1)
            ENDIF

            ! X21
            CALL DLA_GEAXB('T', 'N', 'F', TYPEB, M2, M1, N1, N1, -ONE, A(1, M1+1), LDA, X(1,1), LDX, B(1,1), LDB, &
                & SCALE, X(M1+1,1), LDX, WORK)
            CALL DLA_GEAXB('T', 'N', 'F', TYPED, M2, M1, N1, N1, -SGN, C(1, M1+1), LDC, X(1,1), LDX, D(1,1), LDD, &
                & ONE, X(M1+1,1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('T','N', SGN,  M2, N1, A(M1+1, M1+1), LDA, B(1,1), LDB, C(M1+1, M1+1), LDC, D(1,1), LDD, &
                & X(M1+1, 1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
            ENDIF

            !X22
            CALL DLA_GEAXB('T', 'N', 'Q', 'F', M2, M, N1, N2, -ONE, A(1,M1+1), LDA, X(1,1),LDX, B(1,N1+1),LDB, &
                & SCALE,  X(M1+1, N1+1), LDX, WORK)

            CALL DLA_GEAXB('T', 'N', 'F', TYPEB, M2, M1, N2,N2, -ONE, A(1,M1+1), LDA, X(1,N1+1), LDX, B(N1+1, N1+1), LDB, &
                & ONE, X(M1+1,N1+1), LDX, WORK)

            CALL DLA_GEAXB('T', 'N', 'T', 'F', M2, M, N1, N2, -SGN, C(1,M1+1), LDC, X(1,1),LDX, D(1,N1+1),LDD, &
                & ONE, X(M1+1, N1+1), LDX, WORK)
            CALL DLA_GEAXB('T', 'N', 'F', TYPED, M2, M1, N2,N2, -SGN, C(1,M1+1), LDC, X(1,N1+1), LDX, D(N1+1, N1+1), LDD, &
                & ONE, X(M1+1,N1+1), LDX, WORK)

            CALL DLA_TGSYLV_RECURSIVE('T', 'N', SGN,  M2, N2, A(M1+1,M1+1), LDA, B(N1+1, N1+1), LDB, C(M1+1,M1+1), LDC, &
                & D(N1+1,N1+1), LDD, X(M1+1,N1+1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
                X(M1+1:M, 1:N1) = SCAL * X(M1+1:M, 1:N1)
            ENDIF
        END IF
    ELSE
        ! (T,T)
        IF ( N .LE. M/2 .OR. N .LE. BOUND ) THEN
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF

            !X11
            CALL DLA_TGSYLV_RECURSIVE('T', 'T', SGN,  M1, N, A(1,1), LDA, B, LDB, C(1,1), LDC, D,LDD, X(1,1), LDX, &
                & SCALE, WORK, INFO )

            CALL DLA_GEAXB('T', 'T', 'F', TYPEB, M2, M1, N, N, -ONE, A(1,M1+1), LDA, X(1,1), LDX, B, LDB, SCALE, X(M1+1,1), &
                & LDX, WORK)
            CALL DLA_GEAXB('T', 'T', 'F', TYPED, M2, M1, N, N, -SGN, C(1,M1+1), LDC, X(1,1), LDX, D, LDD, ONE, X(M1+1,1), &
                & LDX, WORK)
            ! X22
            CALL DLA_TGSYLV_RECURSIVE('T', 'T', SGN,  M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D, LDD, &
                & X(M1+1,1), LDX, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M1, 1:N) = SCAL*X(1:M1, 1:N)
                SCALE = SCAL * SCALE
            ENDIF

        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL DLA_TGSYLV_RECURSIVE('T', 'T', SGN,  M, N2, A, LDA, B(N1+1, N1+1), LDB, C, LDC, D(N1+1,N1+1), LDD, &
                & X(1,N1+1), LDX, SCALE, WORK, INFO)
            CALL DLA_GEAXB('T', 'T', 'Q', 'F', M, M, N2, N1, -ONE, A, LDA, X(1,N1+1), LDX, B(1,N1+1), LDB, SCALE, X(1,1), &
                & LDX, WORK)
            CALL DLA_GEAXB('T', 'T', 'T', 'F', M, M, N2, N1, -SGN, C, LDC, X(1,N1+1), LDX, D(1,N1+1), LDD, ONE, X(1,1), &
                & LDX, WORK)
            !X1
            CALL DLA_TGSYLV_RECURSIVE('T', 'T', SGN,  M, N1, A, LDA, B(1,1), LDB, C, LDC, D(1,1), LDD, X(1,1), LDX, &
                & SCAL,  WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                X(1:M, N1+1:N) = SCAL*X(1:M, N1+1:N)
                SCALE = SCAL * SCALE
            ENDIF

        ELSE
            M1 = M/2
            M2 = M - M1
            IF ( A(M1+1, M1) .NE. ZERO) THEN
                M1 = M1 +1
                M2 = M2 -1
            ENDIF
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            !X12
            CALL DLA_TGSYLV_RECURSIVE('T','T', SGN,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, C(1,1), LDC, D(N1+1,N1+1), LDD, &
                & X(1,N1+1), LDX, SCALE, WORK, INFO)

            !X11
            CALL DLA_GEAXB('T', 'T', 'Q', 'F', M1, M1, N2, N1, -ONE, A(1,1), LDA, X(1,N1+1), LDX, B(1,N1+1), LDB, &
                & SCALE, X(1,1), LDX, WORK)
            CALL DLA_GEAXB('T', 'T', 'T', 'F', M1, M1, N2, N1, -SGN, C(1,1), LDC, X(1,N1+1), LDX, D(1,N1+1), LDD, &
                & ONE, X(1,1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('T','T', SGN,  M1, N1, A(1,1), LDA, B(1,1), LDB, C(1,1), LDC, D(1, 1), LDD, &
                & X(1,1), LDX, SCAL, WORK, INFO)

            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, N1+1) = SCAL * X(1:M1, N1+1)
            ENDIF

            !X22
            CALL DLA_GEAXB('T', 'T', 'F', TYPEB, M2, M1, N2, N2, -ONE, A(1, M1+1), LDA, X(1,N1+1), LDX, B(N1+1,N1+1), LDB, &
                & SCALE, X(M1+1,N1+1), LDX, WORK)
            CALL DLA_GEAXB('T', 'T', 'F', TYPED, M2, M1, N2, N2, -SGN, C(1, M1+1), LDC, X(1,N1+1), LDX, D(N1+1,N1+1), LDD, &
                & ONE, X(M1+1,N1+1), LDX, WORK)
            CALL DLA_TGSYLV_RECURSIVE('T','T', SGN,  M2, N2, A(M1+1, M1+1), LDA, B(N1+1,N1+1), LDB, C(M1+1, M1+1), LDC,  &
                & D(N1+1,N1+1), LDD, X(M1+1, N1+1), LDX, SCAL, WORK, INFO)
            IF (SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
            ENDIF

            ! X21

            CALL DLA_GEAXB('T', 'T', 'F', TYPEB, M2, M1, N1, N1, -ONE, A(1,M1+1), LDA, X(1,1),LDX, B(1,1),LDB, &
                & SCALE,  X(M1+1, 1), LDX, WORK)
            CALL DLA_GEAXB('T', 'T', 'Q', 'F', M2, M, N2, N1, -ONE, A(1,M1+1), LDA, X(1,N1+1),LDX, B(1,N1+1),LDB, &
                & ONE,  X(M1+1, 1), LDX, WORK)

            CALL DLA_GEAXB('T', 'T', 'F', TYPED, M2, M1, N1, N1, -SGN, C(1,M1+1), LDC, X(1,1),LDX, D(1,1),LDD, &
                & ONE,  X(M1+1, 1), LDX, WORK)
            CALL DLA_GEAXB('T', 'T', 'T', 'F', M2, M, N2, N1, -SGN, C(1,M1+1), LDC, X(1,N1+1),LDX, D(1,N1+1),LDD, &
                & ONE,  X(M1+1, 1), LDX, WORK)


            CALL DLA_TGSYLV_RECURSIVE('T', 'T',SGN,  M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, C(M1+1,M1+1), LDC, D(1,1),&
                & LDD, X(M1+1,1), LDX, SCAL, WORK, INFO)
            IF ( SCAL .NE. ONE) THEN
                SCALE = SCAL * SCALE
                X(1:M1, 1:N) = SCAL * X(1:M1, 1:N)
                X(M1+1:M, N1+1:N) = SCAL * X(M1+1:M, N1+1:N)
            ENDIF

        ENDIF

    ENDIF
END SUBROUTINE


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

!> \brief Recursive Algorithm for the generalized coupled Sylvester equation
!
!  Definition:
!  ===========
!       SUBROUTINE SLA_TGCSYLV_RECURSIVE ( TRANSA, TRANSB, SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
!                                       & E, LDE, F, LDF,  SCALE, WORK, INFO)
!           CHARACTER TRANSA(1), TRANSB(1)
!           REAL SGN1, SGN2, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
!           REAL A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TGCSYLV_RECURSIVE solves a generalized coupled  Sylvester equation of the following form
!>
!>    op1(A) * R  + SGN1 * L  * op2(B) = SCALE * E                              (1)
!>    op1(C) * R  + SGN2 * L  * op2(D) = SCALE * F
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by SGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using BLAS level 3 operations.
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
!> \param[in] SGN1
!> \verbatim
!>          SGN1 is REAL, allowed values: +/-1
!>          Specifies the sign between in the first equation.
!> \endverbatim
!>
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is REAL, allowed values: +/-1
!>          Specifies the sign between in the second equation.
!> \endverbatim
!>
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
!>          A is REAL array, dimension (LDA,M)
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
!>          B is REAL array, dimension (LDB,N)
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
!>          C is REAL array, dimension (LDC,M)
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
!>          D is REAL array, dimension (LDD,N)
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
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (LDE,N)
!>          On input, the matrix E contains the right hand side E.
!>          On output, the matrix E contains the solution R of Equation (1)
!>          as selected by TRANSA, TRANSB, and SGN1/SGN2.
!>          Right hand side E and the solution R are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of the array E.  LDE >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is REAL array, dimension (LDF,N)
!>          On input, the matrix F contains the right hand side F.
!>          On output, the matrix F contains the solution L of Equation (1)
!>          as selected by TRANSA, TRANSB, and SGN1/SGN2.
!>          Right hand side F and the solution L are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F.  LDE >= max(1,M).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is REAL array, dimension 1
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
!> \ingroup sgltgsylv
!
RECURSIVE SUBROUTINE SLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE,  F, LDF, &
        & SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Input
    CHARACTER TRANSA(1), TRANSB(1)
    INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF , INFO
    REAL SGN1, SGN2
    REAL A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)
    REAL SCALE
    REAL WORK(*)

    ! External Functions
    EXTERNAL SGEMM
    EXTERNAL SLA_TGCSYLV_KERNEL_44NN
    EXTERNAL SLA_TGCSYLV_KERNEL_44NT
    EXTERNAL SLA_TGCSYLV_KERNEL_44TN
    EXTERNAL SLA_TGCSYLV_KERNEL_44TT
    EXTERNAL SLA_TGCSYLV_L2
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TGCSYLV_L2_REORDER
    EXTERNAL SLA_TGCSYLV_L2_UNOPT
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME


    ! Local Variables
    INTEGER ININFO
    CHARACTER TYPEB, TYPED
    LOGICAL BTRANSA, BTRANSB
    INTEGER M1, M2, N1, N2, IINFO, I
    REAL ZERO, ONE
    PARAMETER(ZERO = 0.0d0, ONE = 1.0d0)
    REAL SCAL
    INTEGER BOUND
    ! INTEGER  MB, NB, MAXNB, ISOLVER
    PARAMETER(BOUND = 4)

    ! CHECK Input
    ININFO = INFO
    INFO = 0
    BTRANSA = LSAME(TRANSA, 'T')
    BTRANSB = LSAME(TRANSB, 'T')

    IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'N')) THEN
        INFO = -1
    ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'N')) THEN
        INFO = -2
    ELSE IF ( SGN1 .NE. -ONE .AND. SGN1 .NE. ONE ) THEN
        INFO = -3
    ELSE IF ( SGN2 .NE. -ONE .AND. SGN2 .NE. ONE ) THEN
        INFO = -4
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -5
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -6
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -8
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -10
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -12
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -14
    ELSE IF ( LDE .LT. MAX(1, M)) THEN
        INFO = -16
    ELSE IF ( LDF .LT. MAX(1, M)) THEN
        INFO = -18
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TGCSYLV_RECURSIVE', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        INFO = 1
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
                CALL SLA_TGCSYLV_KERNEL_44NN(SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, IINFO)
            ELSE IF ( .NOT. BTRANSA .AND. BTRANSB) THEN
                !(N,T)
                CALL SLA_TGCSYLV_KERNEL_44NT(SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, IINFO)
            ELSE IF ( BTRANSA .AND. .NOT. BTRANSB) THEN
                !(T,N)
                CALL SLA_TGCSYLV_KERNEL_44TN(SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, IINFO)

            ELSE IF ( BTRANSA .AND. BTRANSB) THEN
                !(T,T)
                CALL SLA_TGCSYLV_KERNEL_44TT(SGN1, SGN2,  M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, IINFO)
            ENDIF
            INFO = IINFO
        ! ELSE
        !     SCAL = ONE;
        !     SELECT CASE(ISOLVER)
        !     CASE(1)
        !         IF ( MAXNB .LE. 32 ) THEN
        !             CALL SLA_TGCSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !                 & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !         ELSE IF ( MAXNB .LE. 64) THEN
        !             CALL SLA_TGCSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !                 & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !         ELSE IF ( MAXNB .LE. 96) THEN
        !             CALL SLA_TGCSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !                 & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !         ELSE IF ( MAXNB .LE. 128) THEN
        !             CALL SLA_TGCSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !                 & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !         ELSE
        !             CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !                 & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !         END IF
        !     CASE(2)
        !         IF ( M .GT. 128 .OR. N .GT. 128) THEN
        !             CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !                 & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !         ELSE
        !             CALL SLA_TGCSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !                 & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !         END IF
        !     CASE(3)
        !         CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !             & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !     CASE(4)
        !         CALL SLA_TGCSYLV_L2(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !             & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !     CASE(5)
        !         CALL SLA_TGCSYLV_L2_UNOPT(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !             & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !     CASE(6)
        !         ! Do not start the recursive algorithm again.
        !         CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
        !             & C, LDC, D, LDD,  E, LDE, F, LDF, SCAL, WORK, IINFO)
        !     END SELECT
        !     IF ( IINFO.NE. 0 ) THEN
        !         INFO = IINFO
        !     END IF
        !     IF ( SCAL .NE. ONE) THEN
        !         SCALE = SCAL
        !     END IF
        !
        ! END IF
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


            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D,LDD, &
                & E(M1+1,1), LDE, F(M1+1,1), LDF,  SCALE, WORK, INFO )

            CALL SGEMM("N", "N", M1, N, M2, -ONE, A(1,M1+1), LDA, E(M1+1,1), LDE, SCALE, E(1,1), LDE)
            CALL SGEMM("N", "N", M1, N, M2, -ONE, C(1,M1+1), LDC, E(M1+1,1), LDE, SCALE, F(1,1), LDF)

            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M1, N, A(1,1), LDA, B, LDB, C(1,1), LDC, D,LDD, &
                & E(1,1), LDE, F(1,1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(M1+1:M, 1:N) = SCAL*E(M1+1:M, 1:N)
                F(M1+1:M, 1:N) = SCAL*F(M1+1:M, 1:N)
                SCALE = SCAL * SCALE
            ENDIF


        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M, N1, A, LDA, B(1,1), LDB, C, LDC, D(1,1),LDD, &
                & E(1,1), LDE, F(1,1), LDF,  SCALE, WORK, INFO)

            CALL SGEMM("N", "N", M, N2, N1, -SGN1, F(1,1), LDF, B(1,N1+1), LDB, SCALE, E(1,N1+1), LDE)
            CALL SGEMM("N", "N", M, N2, N1, -SGN2, F(1,1), LDF, D(1,N1+1), LDD, SCALE, F(1,N1+1), LDF)

            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M, N2, A, LDA, B(N1+1,N1+1), LDB, C, LDC, D(N1+1,N1+1),LDD, &
                & E(1,N1+1), LDE, F(1,N1+1), LDF,  SCAL, WORK, INFO  )

            IF (SCAL.NE.ONE) THEN
                E(1:M, 1:N1) = SCAL*E(1:M, 1:N1)
                F(1:M, 1:N1) = SCAL*F(1:M, 1:N1)
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

            ! E21
            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M2, N1, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D,LDD, &
                & E(M1+1,1), LDE, F(M1+1,1), LDF,  SCALE, WORK, INFO  )

            ! E11
            CALL SGEMM("N", "N", M1, N1, M2, -ONE, A(1,M1+1), LDA, E(M1+1,1), LDE, SCALE, E(1,1), LDE)
            CALL SGEMM("N", "N", M1, N1, M2, -ONE, C(1,M1+1), LDC, E(M1+1,1), LDE, SCALE, F(1,1), LDF)

            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M1, N1, A, LDA, B, LDB, C, LDC, D,LDD, &
                & E, LDE, F, LDF,  SCAL, WORK, INFO )
            IF (SCAL.NE.ONE) THEN
                E(M1+1:M, 1:N1) = SCAL*E(M1+1:M, 1:N1)
                F(M1+1:M, 1:N1) = SCAL*F(M1+1:M, 1:N1)
                SCALE = SCAL * SCALE
            ENDIF

            ! E22
            CALL SGEMM("N","N", M2, N2, N1, -SGN1, F(M1+1,1), LDF, B(1,N1+1), LDB, SCALE, E(M1+1,N1+1), LDE)
            CALL SGEMM("N","N", M2, N2, N1, -SGN2, F(M1+1,1), LDF, D(1,N1+1), LDD, SCALE, F(M1+1,N1+1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M2, N2, A(M1+1,M1+1), LDA, B(N1+1,N1+1), LDB, C(M1+1,M1+1), LDC, &
                & D(N1+1,N1+1),LDD, E(M1+1,N1+1), LDE, F(M1+1,N1+1), LDF,  SCAL, WORK, INFO )
            IF (SCAL.NE.ONE) THEN
                E(1:M, 1:N1) = SCAL*E(1:M, 1:N1)
                F(1:M, 1:N1) = SCAL*F(1:M, 1:N1)
                SCALE = SCAL * SCALE
            ENDIF


            ! E12
            CALL SGEMM("N", "N", M1, N2, M2, -ONE, A(1,M1+1), LDA, E(M1+1,N1+1), LDE, SCALE, E(1,N1+1), LDE)
            CALL SGEMM("N", "N", M1, N2, M2, -ONE, C(1,M1+1), LDC, E(M1+1,N1+1), LDE, SCALE, F(1,N1+1), LDF)
            CALL SGEMM("N", "N", M1, N2, N1, -SGN1, F(1,1), LDF, B(1,N1+1), LDB, ONE, E(1,N1+1), LDE)
            CALL SGEMM("N", "N", M1, N2, N1, -SGN2, F(1,1), LDF, D(1,N1+1), LDD, ONE, F(1,N1+1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('N', 'N', SGN1, SGN2, M1, N2, A, LDA, B(N1+1,N1+1), LDB, C, LDC, &
                & D(N1+1,N1+1),LDD, E(1,N1+1), LDE, F(1,N1+1), LDF,  SCAL, WORK, INFO )
            IF (SCAL.NE.ONE) THEN
                E(1:M, 1:N1) = SCAL*E(1:M, 1:N1)
                F(1:M, 1:N1) = SCAL*F(1:M, 1:N1)
                E(M1+1:M, N1+1:N) = SCAL*E(M1+1:M,N1+1:N)
                F(M1+1:M, N1+1:N) = SCAL*F(M1+1:M,N1+1:N)
                SCALE = SCAL * SCALE
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

            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D,LDD, &
                & E(M1+1,1), LDE, F(M1+1,1), LDF,  SCALE, WORK, INFO )

            CALL SGEMM("N", "N", M1, N, M2, -ONE, A(1,M1+1), LDA, E(M1+1,1), LDE, SCALE, E(1,1), LDE)
            CALL SGEMM("N", "N", M1, N, M2, -ONE, C(1,M1+1), LDC, E(M1+1,1), LDE, SCALE, F(1,1), LDF)

            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M1, N, A(1,1), LDA, B, LDB, C(1,1), LDC, D,LDD, &
                & E(1,1), LDE, F(1,1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(M1+1:M, 1:N) = SCAL*E(M1+1:M, 1:N)
                F(M1+1:M, 1:N) = SCAL*F(M1+1:M, 1:N)
                SCALE = SCAL * SCALE
            ENDIF


        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M, N2, A, LDA, B(N1+1,N1+1), LDB, C, LDC, D(N1+1,N1+1),LDD, &
                & E(1,N1+1), LDE, F(1,N1+1), LDF,  SCALE, WORK, INFO )

            CALL SGEMM("N", "T", M, N1, N2, -SGN1, F(1,N1+1), LDF, B(1,N1+1), LDB, SCALE, E, LDE)
            CALL SGEMM("N", "T", M, N1, N2, -SGN2, F(1,N1+1), LDF, D(1,N1+1), LDD, SCALE, F, LDF)
            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M, N1, A, LDA, B, LDB, C, LDC, D,LDD, &
                & E(1,1), LDE, F(1,1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(1:M, N1+1:N) = SCAL*E(1:M, N1+1:N)
                F(1:M, N1+1:N) = SCAL*F(1:M, N1+1:N)
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

            ! E22
            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M2, N2, A(M1+1,M1+1), LDA, B(N1+1,N1+1), LDB, C(M1+1,M1+1), LDC, &
                & D(N1+1,N1+1),LDD, E(M1+1,N1+1), LDE, F(M1+1,N1+1), LDF,  SCALE, WORK, INFO)

            ! E12
            CALL SGEMM("N", "N", M1, N2, M2, -ONE, A(1,M1+1), LDA, E(M1+1,N1+1), LDE, SCALE, E(1,N1+1), LDE)
            CALL SGEMM("N", "N", M1, N2, M2, -ONE, C(1,M1+1), LDC, E(M1+1,N1+1), LDE, SCALE, F(1,N1+1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, C(1,1), LDC, &
                & D(N1+1,N1+1),LDD, E(1,N1+1), LDE, F(1,N1+1), LDF,  SCAL, WORK, INFO)

            IF (SCAL.NE.ONE) THEN
                E(M1+1:M, N1+1:N) = SCAL*E(M1+1:M, N1+1:N)
                F(M1+1:M, N1+1:N) = SCAL*F(M1+1:M, N1+1:N)
                SCALE = SCAL * SCALE
            ENDIF

            ! E21
            CALL SGEMM("N", "T", M2, N1, N2, -SGN1, F(M1+1,N1+1), LDF, B(1,N1+1), LDB, SCALE, E(M1+1,1), LDE)
            CALL SGEMM("N", "T", M2, N1, N2, -SGN2, F(M1+1,N1+1), LDF, D(1,N1+1), LDD, SCALE, F(M1+1,1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, C(M1+1,M1+1), LDC, &
                & D(1,1),LDD, E(M1+1,1), LDE, F(M1+1,1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(1:M, N1+1:N) = SCAL*E(1:M, N1+1:N)
                F(1:M, N1+1:N) = SCAL*F(1:M, N1+1:N)
                SCALE = SCAL * SCALE
            ENDIF


            CALL SGEMM("N", "N", M1, N1, M2, -ONE, A(1,M1+1), LDA, E(M1+1,1), LDE, SCALE, E(1,1), LDE)
            CALL SGEMM("N", "N", M1, N1, M2, -ONE, C(1,M1+1), LDC, E(M1+1,1), LDE, SCALE, F(1,1), LDF)
            CALL SGEMM("N", "T", M1, N1, N2, -SGN1, F(1,N1+1), LDF, B(1,N1+1), LDB, ONE, E(1,1), LDE)
            CALL SGEMM("N", "T", M1, N1, N2, -SGN2, F(1,N1+1), LDF, D(1,N1+1), LDD, ONE, F(1,1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('N', 'T', SGN1, SGN2, M1, N1, A(1,1), LDA, B(1,1), LDB, C(1,1), LDC, &
                & D(1,1),LDD, E(1,1), LDE, F(1,1), LDF,  SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                E(1:M, N1+1:N) = SCAL*E(1:M, N1+1:N)
                F(1:M, N1+1:N) = SCAL*F(1:M, N1+1:N)
                E(M1+1:M, 1:N1) = SCAL*E(M1+1:M, 1:N1)
                F(M1+1:M, 1:N1) = SCAL*F(M1+1:M, 1:N1)
                SCALE = SCAL * SCALE
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

            CALL SLA_TGCSYLV_RECURSIVE('T', 'N', SGN1, SGN2, M1, N, A, LDA, B, LDB, C, LDC, D,LDD, &
                & E, LDE, F, LDF,  SCALE, WORK, INFO )

            CALL SGEMM("T", "N", M2, N, M1, -ONE, A(1,M1+1), LDA, E(1,1), LDE, SCALE, E(M1+1,1), LDE)
            CALL SGEMM("T", "N", M2, N, M1, -ONE, C(1,M1+1), LDC, E(1,1), LDE, SCALE, F(M1+1,1), LDF)

            CALL SLA_TGCSYLV_RECURSIVE('T', 'N', SGN1, SGN2, M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D,LDD, &
                & E(M1+1,1), LDE, F(M1+1,1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(1:M1, 1:N) = SCAL*E(1:M1, 1:N)
                F(1:M1, 1:N) = SCAL*F(1:M1, 1:N)

                SCALE = SCAL * SCALE
            ENDIF
        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF ( B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL SLA_TGCSYLV_RECURSIVE('T', 'N', SGN1, SGN2, M, N1, A, LDA, B(1,1), LDB, C, LDC, D(1,1),LDD, &
                & E(1,1), LDE, F(1,1), LDF,  SCALE, WORK, INFO )

            CALL SGEMM("N", "N", M, N2, N1, -SGN1, F(1,1), LDF, B(1,N1+1), LDB, SCALE, E(1,N1+1), LDE)
            CALL SGEMM("N", "N", M, N2, N1, -SGN2, F(1,1), LDF, D(1,N1+1), LDD, SCALE, F(1,N1+1), LDF)

            CALL SLA_TGCSYLV_RECURSIVE('T', 'N', SGN1, SGN2, M, N2, A, LDA, B(N1+1,N1+1), LDB, C, LDC, D(N1+1,N1+1),LDD, &
                & E(1,N1+1), LDE, F(1,N1+1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(1:M, 1:N1) = SCAL*E(1:M, 1:N1)
                F(1:M, 1:N1) = SCAL*F(1:M, 1:N1)
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

            ! E11
            CALL SLA_TGCSYLV_RECURSIVE('T','N', SGN1, SGN2,  M1, N1, A(1,1), LDA, B(1,1), LDB, C(1,1), LDC, D(1,1), LDD, &
                & E(1,1), LDE, F(1,1), LDF, SCALE, WORK, INFO)

            ! E21
            CALL SGEMM("T", "N", M2, N1, M1, -ONE, A(1,M1+1), LDA, E(1,1), LDE, SCALE, E(M1+1,1), LDE)
            CALL SGEMM("T", "N", M2, N1, M1, -ONE, C(1,M1+1), LDC, E(1,1), LDE, SCALE, F(M1+1,1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('T','N', SGN1, SGN2,  M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, C(M1+1,M1+1), LDC, &
                & D(1,1), LDD, E(M1+1,1), LDE, F(M1+1,1), LDF, SCAL, WORK, INFO)

            IF (SCAL.NE.ONE) THEN
                E(1:M1, 1:N1) = SCAL * E(1:M1, 1:N1)
                F(1:M1, 1:N1) = SCAL * F(1:M1, 1:N1)
                SCALE = SCALE * SCAL
            END IF

            ! E12
            CALL SGEMM("N", "N", M1, N2, N1, -SGN1, F(1,1), LDF, B(1,N1+1), LDB, SCALE, E(1,N1+1), LDE)
            CALL SGEMM("N", "N", M1, N2, N1, -SGN2, F(1,1), LDF, D(1,N1+1), LDD, SCALE, F(1,N1+1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('T','N', SGN1, SGN2,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, C(1,1), LDC, &
                & D(N1+1,N1+1), LDD,  E(1,N1+1), LDE, F(1,N1+1), LDF, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                E(1:M, 1:N1) = SCAL * E(1:M, 1:N1)
                F(1:M, 1:N1) = SCAL * F(1:M, 1:N1)
                SCALE = SCALE * SCAL
            END IF

            ! E22
            CALL SGEMM("T", "N", M2, N2, M1, -ONE, A(1,M1+1), LDA, E(1,N1+1), LDE, SCALE, E(M1+1,N1+1), LDE)
            CALL SGEMM("T", "N", M2, N2, M1, -ONE, C(1,M1+1), LDC, E(1,N1+1), LDE, SCALE, F(M1+1,N1+1), LDF)
            CALL SGEMM("N", "N", M2, N2, N1, -SGN1, F(M1+1,1), LDF, B(1,N1+1), LDB, ONE , E(M1+1,N1+1), LDE)
            CALL SGEMM("N", "N", M2, N2, N1, -SGN2, F(M1+1,1), LDF, D(1,N1+1), LDD, ONE , F(M1+1,N1+1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('T','N', SGN1, SGN2,  M2, N2, A(M1+1,M1+1), LDA, B(N1+1,N1+1), LDB, C(M1+1,M1+1), LDC, &
                & D(N1+1,N1+1), LDD,  E(M1+1,N1+1), LDE, F(M1+1,N1+1), LDF, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                E(1:M, 1:N1) = SCAL * E(1:M, 1:N1)
                F(1:M, 1:N1) = SCAL * F(1:M, 1:N1)
                E(1:M1, N1+1:N) = SCAL * E(1:M1, N1+1:N)
                F(1:M1, N1+1:N) = SCAL * F(1:M1, N1+1:N)
                SCALE = SCALE * SCAL
            END IF
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

            CALL SLA_TGCSYLV_RECURSIVE('T', 'T', SGN1, SGN2, M1, N, A, LDA, B, LDB, C, LDC, D,LDD, &
                & E, LDE, F, LDF,  SCALE, WORK, INFO )

            CALL SGEMM("T", "N", M2, N, M1, -ONE, A(1,M1+1), LDA, E(1,1), LDE, SCALE, E(M1+1,1), LDE)
            CALL SGEMM("T", "N", M2, N, M1, -ONE, C(1,M1+1), LDC, E(1,1), LDE, SCALE, F(M1+1,1), LDF)

            CALL SLA_TGCSYLV_RECURSIVE('T', 'T', SGN1, SGN2, M2, N, A(M1+1,M1+1), LDA, B, LDB, C(M1+1,M1+1), LDC, D,LDD, &
                & E(M1+1,1), LDE, F(M1+1,1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(1:M1, 1:N) = SCAL*E(1:M1, 1:N)
                F(1:M1, 1:N) = SCAL*F(1:M1, 1:N)

                SCALE = SCAL * SCALE
            ENDIF


        ELSE IF ( M .LE. N/2  .OR. M .LE. BOUND ) THEN
            N1 = N/2
            N2 = N - N1
            IF (B(N1+1,N1) .NE. ZERO .OR. D(N1+1,N1) .NE. ZERO) THEN
                N1 = N1 +1
                N2 = N2 -1
            ENDIF

            CALL SLA_TGCSYLV_RECURSIVE('T', 'T', SGN1, SGN2, M, N2, A, LDA, B(N1+1,N1+1), LDB, C, LDC, D(N1+1,N1+1),LDD, &
                & E(1,N1+1), LDE, F(1,N1+1), LDF,  SCALE, WORK, INFO )

            CALL SGEMM("N", "T", M, N1, N2, -SGN1, F(1,N1+1), LDF, B(1,N1+1), LDB, SCALE, E, LDE)
            CALL SGEMM("N", "T", M, N1, N2, -SGN2, F(1,N1+1), LDF, D(1,N1+1), LDD, SCALE, F, LDF)
            CALL SLA_TGCSYLV_RECURSIVE('T', 'T', SGN1, SGN2, M, N1, A, LDA, B, LDB, C, LDC, D,LDD, &
                & E(1,1), LDE, F(1,1), LDF,  SCAL, WORK, INFO )

            IF (SCAL.NE.ONE) THEN
                E(1:M, N1+1:N) = SCAL*E(1:M, N1+1:N)
                F(1:M, N1+1:N) = SCAL*F(1:M, N1+1:N)
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

            ! E12
            CALL SLA_TGCSYLV_RECURSIVE('T','T', SGN1, SGN2,  M1, N2, A(1,1), LDA, B(N1+1,N1+1), LDB, C(1,1), LDC, &
                & D(N1+1,N1+1), LDD, E(1,N1+1), LDE, F(1,N1+1), LDF, SCALE, WORK, INFO)

            ! E22
            CALL SGEMM("T", "N", M2, N2, M1, -ONE, A(1,M1+1), LDA, E(1,N1+1), LDE, SCALE, E(M1+1,N1+1), LDE)
            CALL SGEMM("T", "N", M2, N2, M1, -ONE, C(1,M1+1), LDC, E(1,N1+1), LDE, SCALE, F(M1+1,N1+1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('T','T', SGN1, SGN2,  M2, N2, A(M1+1,M1+1), LDA, B(N1+1,N1+1), LDB, C(M1+1,M1+1), LDC, &
                & D(N1+1,N1+1), LDD, E(M1+1,N1+1), LDE, F(M1+1,N1+1), LDF, SCAL, WORK, INFO)

            IF (SCAL.NE.ONE) THEN
                E(1:M1, N1+1:N) = SCAL*E(1:M1, N1+1:N)
                F(1:M1, N1+1:N) = SCAL*F(1:M1, N1+1:N)
                SCALE = SCALE * SCAL
            END IF

            ! E11
            CALL SGEMM("N", "T", M1, N1, N2, -SGN1, F(1,N1+1), LDF, B(1,N1+1), LDB, SCALE, E(1,1), LDE)
            CALL SGEMM("N", "T", M1, N1, N2, -SGN2, F(1,N1+1), LDF, D(1,N1+1), LDD, SCALE, F(1,1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('T','T', SGN1, SGN2,  M1, N1, A(1,1), LDA, B(1,1), LDB, C(1,1), LDC, &
                & D(1,1), LDD, E(1,1), LDE, F(1,1), LDF, SCAL, WORK, INFO)
            IF (SCAL.NE.ONE) THEN
                E(1:M, N1+1:N) = SCAL*E(1:M, N1+1:N)
                F(1:M, N1+1:N) = SCAL*F(1:M, N1+1:N)
                SCALE = SCAL * SCALE
            ENDIF

            !E21

            CALL SGEMM("T", "N", M2, N1, M1, -ONE, A(1,M1+1), LDA, E(1,1), LDE, SCALE, E(M1+1,1), LDE)
            CALL SGEMM("T", "N", M2, N1, M1, -ONE, C(1,M1+1), LDC, E(1,1), LDE, SCALE, F(M1+1,1), LDF)
            CALL SGEMM("N", "T", M2, N1, N2, -SGN1, F(M1+1,N1+1), LDF, B(1,N1+1), LDB, ONE, E(M1+1,1), LDE)
            CALL SGEMM("N", "T", M2, N1, N2, -SGN2, F(M1+1,N1+1), LDF, D(1,N1+1), LDD, ONE, F(M1+1,1), LDF)
            CALL SLA_TGCSYLV_RECURSIVE('T','T', SGN1, SGN2,  M2, N1, A(M1+1,M1+1), LDA, B(1,1), LDB, C(M1+1,M1+1), LDC, &
                & D(1,1), LDD, E(M1+1,1), LDE, F(M1+1,1), LDF, SCAL, WORK, INFO)

            IF (SCAL.NE.ONE) THEN
                E(1:M, N1+1:N) = SCAL*E(1:M, N1+1:N)
                F(1:M, N1+1:N) = SCAL*F(1:M, N1+1:N)
                E(1:M1, 1:N1) = SCAL*E(1:M1, 1:N1)
                F(1:M1, 1:N1) = SCAL*F(1:M1, 1:N1)
                SCALE = SCAL * SCALE
            ENDIF

        ENDIF

    ENDIF
END SUBROUTINE


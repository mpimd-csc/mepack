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

!> \brief Recursive Blocking Algorithm for the generalized Lyapunov equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLA_TGLYAP_RECURSIVE ( TRANS, M,  A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
!           CHARACTER TRANS(1)
!           REAL SCALE
!           INTEGER M, N, LDA, LDB, LDX, INFO
!           REAL A(LDA, *), B(LDB, *) , X(LDX, *)
!           REAL WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_TGLYAP_RECURSIVE solves a generalized Lyapunov  equation of the following forms
!>
!>    A * X * B^T + B * X * A^T = SCALE * Y                                              (1)
!>
!> or
!>
!>    A^T * X * B + B^T * X * A =  SCALE * Y                                             (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, B is a M-by-M upper triangular,
!> and X and Y are symmetric  M-by-M matrices.
!> Typically the matrix pencil (A,B) is created by SGGES from LAPACK.
!> \endverbatim
!> \attention The algorithm is implemented using recursive blocking.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A and B:
!>          == 'N':  Equation (1) is solved.
!>          == 'T':  Equation (2) is solved.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and B.  M >= 0.
!> \endverbatim
!>
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
!>          B is REAL array, dimension (LDB,M)
!>          The matrix B must be upper triangular.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1) or (2)
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
!>          SCALE is REAL
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is REAL array, dimension M*M
!>          Workspace for the algorithm.
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
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date October 2023
!> \ingroup sgltglyap
!
RECURSIVE SUBROUTINE SLA_TGLYAP_RECURSIVE ( TRANS, M, A, LDA, B, LDB, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANS(1)
    REAL SCALE
    INTEGER M, LDA, LDX, INFO, LDB
    REAL A(LDA, *), X(LDX, *), B(LDB,*)
    REAL WORK(*)


    ! Local Variables
    INTEGER ININFO
    REAL SCAL
    REAL MAT(4,4)
    REAL RHS(4)
    INTEGER IPIV(4), JPIV(4)

    LOGICAL BTRANS
    REAL ZERO , ONE, SGN, HALF, TWO
    PARAMETER (ZERO = 0.0, ONE = 1.0, SGN = 1.0, HALF = 0.5, TWO = 2.0 )
    INTEGER IZERO, IONE
    PARAMETER(IZERO = 0, IONE = 1 )

    INTEGER M1, M2, IINFO, I, INFO1


    INTEGER BOUND
    PARAMETER(BOUND = 2)

    ! EXTERNAL FUNCTIONS
    EXTERNAL SGEMM
    EXTERNAL SGESC2
    EXTERNAL SGETC2X
    EXTERNAL SLA_STRMM_FIX
    EXTERNAL SLA_ITRANSPOSE
    EXTERNAL SLA_SMALL_SOLVE4
    EXTERNAL SLA_TGSYLV_RECURSIVE
    EXTERNAL SLA_STRMM_FIX_RIGHT
    EXTERNAL SSYR2K
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME

    ! Check Input
    ININFO = INFO
    INFO = 0
    INFO1 = 1
    BTRANS = LSAME(TRANS, 'T')

    IF ( .NOT. BTRANS .AND. .NOT. LSAME(TRANS, 'N')) THEN
        INFO = -1
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -2
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -4
    ELSE IF ( LDB .LT. MAX(1, M)) THEN
        INFO = -6
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -8
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_TGLYAP_RECURSIVE', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1 ) THEN
        INFO = MAX(1, M*M)
        RETURN
    END IF


    ! Quick Return
    SCALE = ONE
    SCAL = ONE
    INFO = 0
    IINFO = 0
    IF ( M .EQ. 0 ) THEN
        RETURN
    END IF


    SCALE = ONE
    ! Solution Kernels for Small Systems.
    IF (M .LE. BOUND ) THEN
        IINFO = 0
        IF ( .NOT. BTRANS ) THEN
            SELECT CASE(M)
            CASE(1)
                X(1,1) = (X(1,1) / (2*A(1,1)*B(1,1)))
            CASE(2)
                IF ( A(2,1) .EQ. ZERO ) THEN
                    X(2,2) = (X(2,2) / (TWO*A(2,2)*B(2,2)))
                    X(1,2) = (X(1,2) - A(1,2)*X(2,2)*B(2,2) - B(1,2)*X(2,2)*A(2,2))/(A(1,1)*B(2,2)+B(1,1)*A(2,2))
                    X(1,1) = (X(1,1) - TWO*(A(1,1)*X(1,2)+A(1,2)*X(2,2))*B(1,2) &
                        &    - TWO*B(1,1)*X(1,2)*A(1,2)) / (TWO*A(1,1)*B(1,1))
                    X(2,1) = X(1,2)
                ELSE
                    I = 0
                    DO
                        MAT(1,1) = B(1,1)*A(1,1) + A(1,1)*B(1,1)
                        MAT(1,2) = B(1,1)*A(1,2) + A(1,1)*B(1,2)
                        MAT(1,3) = B(1,2)*A(1,1) + A(1,2)*B(1,1)
                        MAT(1,4) = B(1,2)*A(1,2) + A(1,2)*B(1,2)

                        MAT(2,1) = B(1,1)*A(2,1)
                        MAT(2,2) = B(1,1)*A(2,2) + A(1,1)*B(2,2)
                        MAT(2,3) = B(1,2)*A(2,1)
                        MAT(2,4) = B(1,2)*A(2,2) + A(1,2)*B(2,2)

                        MAT(3,1) = A(2,1)*B(1,1)
                        MAT(3,2) = A(2,1)*B(1,2)
                        MAT(3,3) = B(2,2)*A(1,1) + A(2,2)*B(1,1)
                        MAT(3,4) = B(2,2)*A(1,2) + A(2,2)*B(1,2)

                        MAT(4,1) = ZERO
                        MAT(4,2) = A(2,1)*B(2,2)
                        MAT(4,3) = B(2,2)*A(2,1)
                        MAT(4,4) = B(2,2)*A(2,2) + A(2,2)*B(2,2)

                        RHS(1) = X(1,1)
                        RHS(2) = X(1,2)
                        RHS(3) = X(1,2)
                        RHS(4) = X(2,2)

                        IF ( I .EQ. 0 ) THEN
                            CALL SLA_SMALL_SOLVE4( 4, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                            IF ( INFO1 .EQ. 0 ) EXIT
                            I = 1
                        ELSE
                            CALL SGETC2X(4, MAT, 4, IPIV, JPIV, INFO)
                            IF ( INFO .NE. 0 ) EXIT
                            CALL SGESC2(4, MAT, 4, RHS, IPIV, JPIV, SCAL)
                            IF ( SCAL .NE. ONE ) THEN
                                X(1:M, 1:M) = SCAL*X(1:M, 1:M)
                                SCALE = SCALE * SCAL
                                INFO = 1
                            END IF
                            EXIT
                        END IF

                    END DO
                    X(1,1) = RHS(1)
                    X(1,2) = RHS(2)
                    X(2,1) = RHS(2)
                    X(2,2) = RHS(4)
                END IF
            END SELECT
        ELSE IF ( BTRANS ) THEN
            SELECT CASE(M)
            CASE(1)
                X(1,1) = X(1,1) / (TWO*A(1,1)*B(1,1))
            CASE(2)
                IF ( A(2,1) .EQ. ZERO ) THEN
                    X(1,1) = X(1,1) / (TWO*A(1,1)*B(1,1))
                    X(2,1) = (X(2,1) - A(1,2)*X(1,1)*B(1,1) - B(1,2) *X(1,1) * A(1,1)) / (A(1,1)*B(2,2)+B(1,1)*A(2,2))
                    X(1,2) = X(2,1)
                    X(2,2) = (X(2,2) - TWO * B(1,2) * (X(1,1)*A(1,2) + X(2,1) * A(2,2)) &
                         &           - TWO * A(1,2) * X(2,1) *B(2,2)) / (TWO*A(2,2)*B(2,2))
                ELSE
                    I = 0
                    DO
                        MAT(1,1) = B(1,1)*A(1,1) + A(1,1)*B(1,1)
                        MAT(1,2) = B(1,1)*A(2,1)
                        MAT(1,3) = A(2,1)*B(1,1)
                        MAT(1,4) = ZERO

                        MAT(2,1) = B(1,1)*A(1,2) + A(1,1)*B(1,2)
                        MAT(2,2) = B(1,1)*A(2,2) + A(1,1)*B(2,2)
                        MAT(2,3) = A(2,1)*B(1,2)
                        MAT(2,4) = A(2,1)*B(2,2)

                        MAT(3,1) = B(1,2)*A(1,1) + A(1,2)*B(1,1)
                        MAT(3,2) = B(1,2)*A(2,1)
                        MAT(3,3) = B(2,2)*A(1,1) + A(2,2)*B(1,1)
                        MAT(3,4) = B(2,2)*A(2,1)

                        MAT(4,1) = B(1,2)*A(1,2) + A(1,2)*B(1,2)
                        MAT(4,2) = B(1,2)*A(2,2) + A(1,2)*B(2,2)
                        MAT(4,3) = B(2,2)*A(1,2) + A(2,2)*B(1,2)
                        MAT(4,4) = B(2,2)*A(2,2) + A(2,2)*B(2,2)

                        RHS(1) = X(1,1)
                        RHS(2) = X(2,1)
                        RHS(3) = X(2,1)
                        RHS(4) = X(2,2)

                        IF ( I .EQ. 0 ) THEN
                            CALL SLA_SMALL_SOLVE4( 4, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                            IF ( INFO1 .EQ. 0 ) EXIT
                            I = 1
                        ELSE
                            CALL SGETC2X(4, MAT, 4, IPIV, JPIV, INFO)
                            IF ( INFO .NE. 0 ) EXIT
                            CALL SGESC2(4, MAT, 4, RHS, IPIV, JPIV, SCAL)
                            IF ( SCAL .NE. ONE ) THEN
                                X(1:M, 1:M) = SCAL*X(1:M, 1:M)
                                SCALE = SCALE * SCAL
                                INFO = 1
                            END IF
                            EXIT
                        END IF

                    END DO
                    X(1,1) = RHS(1)
                    X(2,1) = RHS(2)
                    X(1,2) = RHS(2)
                    X(2,2) = RHS(4)


                END IF

            END SELECT
        ENDIF
        RETURN
    ENDIF

    ! Recursive Blocked Algorithm.
    M1 = M/2
    M2 = M - M1
    IF ( A(M1+1, M1) .NE. ZERO) THEN
        M1 = M1 +1
        M2 = M2 -1
    ENDIF

    IF ( .NOT. BTRANS ) THEN
         ! (N, T)
         CALL SLA_TGLYAP_RECURSIVE('N', M2, A(M1+1,M1+1), LDA, B(M1+1,M1+1), LDB,  X(M1+1,M1+1), LDX, SCALE, WORK, INFO )

         ! CALL SGEMM("N", "T", M2, M2, M2,  ONE, X(M1+1,M1+1), LDX, B(M1+1,M1+1), LDB, ZERO, WORK, M2)
         CALL SLA_STRMM_FIX_RIGHT("N", "T", M2, M2,  ONE, X(M1+1,M1+1), LDX, B(M1+1,M1+1), LDB, ZERO, WORK, M2)
         CALL SGEMM("N","N",  M1, M2, M2, -ONE, A(1,M1+1), LDA, WORK, M2, SCALE, X(1,M1+1), LDX)

         ! CALL SGEMM("N", "T", M2, M2, M2,  ONE, X(M1+1,M1+1), LDX, A(M1+1,M1+1), LDA, ZERO, WORK, M2)
         CALL SLA_STRMM_FIX_RIGHT("N", "T", M2, M2,  ONE, X(M1+1,M1+1), LDX, A(M1+1,M1+1), LDA, ZERO, WORK, M2)
         CALL SGEMM("N","N",  M1, M2, M2, -ONE, B(1,M1+1), LDB, WORK, M2, ONE, X(1,M1+1), LDX)

         CALL SLA_TGSYLV_RECURSIVE('N', 'T', ONE, M1, M2, A(1,1), LDA, B(M1+1,M1+1),LDB, &
             & B(1,1), LDB, A(M1+1,M1+1), LDA, X(1,M1+1), LDX, SCAL, WORK, INFO)
         IF (SCAL .NE. ONE) THEN
             SCALE = SCALE * SCAL
             X(M1+1:M, M1+1:M) = SCAL * X(M1+1:M, M1+1:M)
         END IF

         CALL SLA_ITRANSPOSE(1,M1,M1+1,M, X, LDX)

         ! CALL SGEMM("N", "N", M1, M2, M1, ONE, A(1,1), LDA, X(1,M1+1), LDX, ZERO, WORK, M1)
         CALL SLA_STRMM_FIX("N", M1, M2, ONE, A(1,1), LDA, X(1,M1+1), LDX, ZERO, WORK, M1)
         CALL SGEMM("N", "N", M1, M2, M2, ONE, A(1,M1+1), LDA, X(M1+1,M1+1), LDX, ONE, WORK, M1)


         CALL SSYR2K("U", "N", M1, M2, -ONE, WORK, M1, B(1,M1+1), LDB, ONE, X(1,1), LDX)

         ! CALL SGEMM("N", "N", M1, M2, M1, ONE, B(1,1), LDB, X(1,M1+1), LDX, ZERO, WORK, M1)
         CALL SLA_STRMM_FIX("N", M1, M2, ONE, B(1,1), LDB, X(1,M1+1), LDX, ZERO, WORK, M1)
         CALL SSYR2K("U", "N", M1, M2, -ONE, WORK, M1, A(1,M1+1), LDA, ONE, X(1,1), LDX)

         CALL SLA_TGLYAP_RECURSIVE('N', M1, A(1,1), LDA, B(1,1), LDB, X(1,1), LDX, SCAL, WORK, INFO )
         IF (SCAL .NE. ONE) THEN
             SCALE = SCALE * SCAL
             X(1:M, M1+1:M) = SCAL * X(1:M, M1+1:M)
             X(M1+1:M, 1:M1) = SCAL * X(M1+1:M, 1:M1)
         END IF

     ELSE
         CALL SLA_TGLYAP_RECURSIVE('T', M1, A(1,1), LDA, B(1,1), LDB, X(1,1), LDX, SCALE, WORK, INFO )

         ! CALL SGEMM("N", "N", M1, M1, M1, ONE, X(1,1), LDX, B(1,1), LDB, ZERO, WORK, M1)
         CALL SLA_STRMM_FIX_RIGHT("N", "N", M1,  M1, ONE, X(1,1), LDX, B(1,1), LDB, ZERO, WORK, M1)
         CALL SGEMM("T", "N", M2, M1, M1, -ONE, A(1,M1+1), LDA, WORK, M1, SCALE, X(M1+1,1), LDX)

         ! CALL SGEMM("N", "N", M1, M1, M1, ONE, X(1,1), LDX, A(1,1), LDA, ZERO, WORK, M1)
         CALL SLA_STRMM_FIX_RIGHT("N", "N",  M1,  M1, ONE, X(1,1), LDX, A(1,1), LDA, ZERO, WORK, M1)
         CALL SGEMM("T", "N", M2, M1, M1, -ONE, B(1,M1+1), LDB, WORK, M1, ONE, X(M1+1,1), LDX)


         CALL SLA_TGSYLV_RECURSIVE('T', 'N', ONE, M2, M1, A(M1+1,M1+1), LDA, B(1,1),LDB, &
             & B(M1+1,M1+1), LDB, A(1,1), LDA, X(M1+1,1), LDX, SCAL, WORK, INFO)
         IF (SCAL .NE. ONE) THEN
             SCALE = SCALE * SCAL
             X(1:M1, 1:M1) = SCAL * X(1:M1, 1:M1)
         END IF

         CALL SLA_ITRANSPOSE(M1+1,M, 1, M1,  X, LDX)

         CALL SGEMM("T", "N", M1, M2, M1, ONE, X(1,1), LDX, A(1,M1+1), LDA, ZERO, WORK, M1)
         CALL SGEMM("T", "N", M1, M2, M2, ONE, X(M1+1,1), LDX, A(M1+1,M1+1), LDA, ONE, WORK, M1)
         ! CALL SLA_STRMM_FIX_RIGHT("T", "N",  M1, M2, ONE, X(M1+1,1), LDX, A(M1+1,M1+1), LDA, ONE, WORK, M1)

         CALL SSYR2K("L", "T", M2, M1, -ONE, WORK, M1, B(1,M1+1), LDB, SCALE, X(M1+1,M1+1), LDX)
         CALL SGEMM("T", "N", M1, M2, M2, ONE, X(M1+1, 1), LDX, B(M1+1,M1+1), LDB, ZERO, WORK, M1)
         ! CALL SLA_STRMM_FIX_RIGHT("T", "N",  M1, M2, ONE, X(M1+1,1), LDX, B(M1+1,M1+1), LDB, ONE, WORK, M1)

         CALL SSYR2K("L", "T", M2, M1, -ONE, WORK, M1,  A(1,M1+1), LDA, ONE, X(M1+1,M1+1), LDX)

         CALL SLA_TGLYAP_RECURSIVE('T', M2, A(M1+1,M1+1), LDA, B(M1+1,M1+1), LDB, X(M1+1,M1+1), LDX, SCAL, WORK, INFO )

         IF (SCAL .NE. ONE) THEN
             SCALE = SCALE * SCAL
             X(1:M, 1:M1) = SCAL * X(1:M, 1:M1)
             X(1:M1,M1+1:M) = SCAL * X(1:M1,M1+1:M)
         END IF

     END IF
     RETURN

     END SUBROUTINE
     SUBROUTINE SLA_STRMM_FIX_RIGHT(TRANSA, TRANSB, M, N , ALPHA, X, LDX, B, LDB, BETA, C, LDC)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
         IMPLICIT NONE
         CHARACTER(1) TRANSA,TRANSB
         INTEGER M, N, LDX, LDB, LDC
         REAL X(LDX, *) , B(LDB, *), C(LDC, *), ALPHA, BETA

         EXTERNAL SGEMM
         EXTERNAL LSAME
         LOGICAL LSAME
         INTEGER NB, K, KH, KB
         REAL ZERO
         PARAMETER(NB = 64, ZERO=0.0)

         IF ( LSAME(TRANSA, "N")) THEN

             IF ( LSAME(TRANSB, "T")) THEN
                 !$omp parallel private(K,KH,KB) if (N/NB .GT. 6)
                 !$omp do schedule(dynamic)
                 DO K = 1, N, NB
                     KH = MIN(N, K + NB  - 1)
                     KB = KH - K + 1

                     IF ( K .GT. 1 ) THEN
                         IF ( B(K,K-1) .NE. ZERO )  THEN
                             CALL SGEMM("N", "T", M, KB, N-K+2 ,  ALPHA, X(1,K-1), LDX, B(K,K-1), LDB, BETA, C(1,K), LDC)
                         ELSE
                             CALL SGEMM("N", "T", M, KB, N-K+1 ,  ALPHA, X(1,K), LDX, B(K,K), LDB, BETA, C(1,K), LDC)
                         END IF

                     ELSE
                         CALL SGEMM("N", "T", M, KB, N-K+1 ,  ALPHA, X(1,K), LDX, B(K,K), LDB, BETA, C(1,K), LDC)

                     END IF
                 END DO
                 !$omp end do
                 !$omp end parallel
             ELSE
                 !$omp parallel private(K,KH,KB) if (N/NB .GT. 6)
                 !$omp do schedule(dynamic)
                 DO K = 1, N, NB
                     KH = MIN(N, K + NB  - 1)
                     KB = KH - K + 1

                     IF ( K .LT. N ) THEN
                         IF ( B(KH+1,KH) .NE. ZERO ) THEN
                             CALL SGEMM("N", "N", M, KB, KH+1, ALPHA , X, LDX, B(1,K), LDB, BETA, C(1,K), LDC)
                         ELSE
                             CALL SGEMM("N", "N", M, KB, KH, ALPHA, X, LDX, B(1,K), LDB, BETA, C(1,K), LDC)

                         END IF
                     ELSE
                         CALL SGEMM("N", "N", M, KB, KH, ALPHA, X, LDX, B(1,K), LDB, BETA, C(1,K), LDC)
                     END IF
                 END DO
                 !$omp end do
                 !$omp end parallel

             END IF
         END IF

 END SUBROUTINE

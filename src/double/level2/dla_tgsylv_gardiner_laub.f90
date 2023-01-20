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

!> \brief Level-2 Gardiner-Laub algorithm for the generalized Sylvester equation
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TGSYLV_GARDINER_LAUB ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, &
!                  & LDC, D, LDD, X, LDX, SCALE, WORK, INFO)
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
!> DLA_TGSYLV_GARDINER_LAUB solves a generalized Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) + op1(C) * X * op2(D) = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X * op2(B) - op1(C) * X * op2(D) = SCALE * Y                              (2)
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by DGGES from LAPACK.
!> \endverbatim
!> \remark The algorithm implements the ideas by Gardiner et. al. on top of BLAS level 2 operations.
!> \attention Scaling is not supported by this algorithm so SCALE == 1.
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
!>          WORK is DOUBLE PRECISION array, dimension LWORK
!>          Workspace for the algorithm.
!>          The workspace needs to queried before the running the computation.
!>          The query is performed by calling the subroutine with INFO == -1 on input.
!>          The required workspace is then returned in INFO.
!> \endverbatim
!>
!> \param[inout] INFO
!> \verbatim
!>          INFO is INTEGER
!>
!>          On input:
!>            == -1 : Perform a workspace query
!>            <> -1: normal operation
!>
!>          On exit, workspace query:
!>            < 0 :  if INFO = -i, the i-th argument had an illegal value
!>            >= 0:  The value of INFO is the required number of elements in the workspace.
!>
!>          On exit, normal operation:
!>            == 0:  successful exit
!>            < 0:  if INFO = -i, the i-th argument had an illegal value
!>            > 0:  The equation is not solved correctly. One of the arising inner
!>                  system got singular.
!> \endverbatim
!>
!>
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Januar 2023
!> \ingroup dbltgsylv
!
SUBROUTINE DLA_TGSYLV_GARDINER_LAUB ( TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, &
        & D, LDD, X, LDX, SCALE, WORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    CHARACTER TRANSA(1), TRANSB(1)
    DOUBLE PRECISION SGN, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDX, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), X(LDX, *)
    DOUBLE PRECISION WORK(2*M,*)


    ! Local Variables
    INTEGER ININFO
    LOGICAL BTRANSA, BTRANSB
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
    INTEGER IONE, ITWO
    PARAMETER(IONE = 1, ITWO = 2 )

    INTEGER I,IINFO
    INTEGER K, L
    DOUBLE PRECISION AA,BB,CC,SS
    DOUBLE PRECISION MAT(4,4),RHS(4)
    LOGICAL LTMP


    ! EXTERNAL FUNCTIONS
    EXTERNAL DGEMV
    EXTERNAL DTRSV
    EXTERNAL DROT
    EXTERNAL DROTG
    EXTERNAL DGETC2X
    EXTERNAL DLA_SMALL_SOLVE4
    EXTERNAL DTRMV
    EXTERNAL DCOPY
    EXTERNAL DSCAL
    EXTERNAL DAXPY
    EXTERNAL XERROR_HANDLER
    EXTERNAL LSAME
    LOGICAL LSAME



    ! Check Input
    ININFO = INFO
    INFO = 0
    BTRANSA = LSAME(TRANSA, 'T')
    BTRANSB = LSAME(TRANSB, 'T')

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
        CALL XERROR_HANDLER('DLA_TGSYLV_GARDINER_LAUB', -INFO)
        RETURN
    END IF

    IF ( ININFO .EQ. -1) THEN
        INFO = 4*M*M+2*M
        RETURN
    END IF

    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    IF ( .NOT. BTRANSA .AND. .NOT. BTRANSB) THEN
        ! ( N, N)
        K = 1
        DO WHILE ( K <= N)
            IF ( K .LT. N ) THEN
                LTMP = (D(K+1,K) .EQ. ZERO .AND. B(K+1, K) .EQ. ZERO)
            ELSE
                LTMP = .FALSE.
            END IF
            IF ( LTMP  .OR. (K .EQ. N)) THEN
                ! X(:,k)=(B(k,k)*A+D(k,k)*C)'\X(:,k);
                ! WORK(1:2*M,1:2*M) = ZERO
                DO I = 1, M
                    CALL DCOPY(M, A(1,I), IONE, WORK(1,I), IONE)
                    CALL DSCAL(M, B(K,K), WORK(1,I), IONE)
                    CALL DAXPY(M, SGN*D(K,K), C(1,I), IONE, WORK(1,I), IONE)
                END DO

                L = 1
                DO WHILE ( L < M )
                    IF ( WORK(L+1,L) .NE. ZERO ) THEN
                        AA = WORK(L,L)
                        BB = WORK(L+1,L)
                        CALL DROTG(AA,BB,CC,SS)
                        CALL DROT(M,WORK(L,L),2*M,WORK(L+1,L),2*M,CC,SS)
                        CALL DROT(IONE,X(L,K),IONE,X(L+1,K),IONE,CC,SS)
                        L = L +2
                    ELSE
                        L = L +1
                    END IF
                END DO
                CALL DTRSV("U","N","N",M,WORK,2*M,X(1,K),IONE)


                CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("N",M,M,SGN,C,LDA, X(1,K),IONE, ZERO, WORK(1,2),IONE)

                DO L = K+1,N
                    CALL DAXPY(M,-B(K,L),WORK(1,1),IONE,X(1,L),IONE)
                    CALL DAXPY(M,-D(K,L),WORK(1,2),IONE,X(1,L),IONE)
                END DO
                K = K + 1
            ELSE
                ! WORK(1:2*M, 1:2*M+1) = ZERO
                ! TMP = [ (B(k,k)*A+D(k,k)*C) ,(B(k,k+1)*A+D(k,k+1)*C) ;
                !         D(k+1,k)*C, (B(k+1,k+1)*A+D(k+1,k+1)*C)]'\reshape(X(:,k:k+1),2*m,1);
                DO L = 1, M
                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L-1),ITWO)
                    CALL DSCAL(M,B(K,K),WORK(1,2*L-1),ITWO)
                    CALL DAXPY(M,SGN*D(K,K),C(1,L),IONE,WORK(1,2*L-1),ITWO)

                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L),ITWO)
                    CALL DSCAL(M,B(K+1,K),WORK(1,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K+1,K),C(1,L),IONE,WORK(1,2*L),ITWO)

                    IF ( B(K+1,K) .NE. ZERO ) THEN
                        CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,B(K,K+1),WORK(2,2*L-1),ITWO)
                    ELSE
                        CALL DCOPY(M,C(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,SGN*D(K,K+1),WORK(2,2*L-1),ITWO)
                    END IF

                    CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L),ITWO)
                    CALL DSCAL(M,B(K+1,K+1),WORK(2,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K+1,K+1),C(1,L),IONE,WORK(2,2*L),ITWO)
                END DO


                CALL DCOPY(M, X(1,K), IONE, WORK(1,2*M+1),ITWO)
                CALL DCOPY(M, X(1,K+1),IONE, WORK(2,2*M+1),ITWO)

                L = 2*M
                DO WHILE ( L >= 1  )
                    LTMP = .FALSE.
                    IF ( L .GE. 4 ) THEN
                        IF ( WORK(L,L-2) .NE. ZERO) THEN
                            LTMP = .TRUE.
                        END IF
                    END IF

                    IF ( LTMP  ) THEN
                        MAT(1:4,1:4) = WORK(L-3:L,L-3:L)
                        RHS(1:4) = WORK(L-3:L,2*M+1)
                        CALL DLA_SMALL_SOLVE4( 4, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L-3:L,2*M+1) = RHS(1:4)
                        IF ( IINFO .NE. 0 )  INFO = 100

                        CALL DGEMV("N",L-4, 4, -ONE, WORK(1,L-3),2*M,WORK(L-3,2*M+1), &
                            & IONE,ONE,WORK(1,2*M+1),IONE)
                        L = L - 4
                    ELSE
                        MAT(1:2,1:2) = WORK(L-1:L,L-1:L)
                        RHS(1:2) = WORK(L-1:L,2*M+1)

                        CALL DLA_SMALL_SOLVE4( 2, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L-1:L,2*M+1) = RHS(1:2)
                        IF ( IINFO .NE. 0 )  INFO = 100


                        CALL DGEMV("N",L-2,2, -ONE, WORK(1,L-1),2*M,WORK(L-1,2*M+1), &
                            & IONE,ONE,WORK(1,2*M+1),IONE)
                        L = L - 2
                    END IF
                END DO


                CALL DCOPY(M, WORK(1,2*M+1),ITWO, X(1,K), IONE)
                CALL DCOPY(M, WORK(2,2*M+1),ITWO, X(1,K+1),IONE)

                ! Precompute the Update January 2021
                CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("N",M,M,SGN,C,LDA, X(1,K),IONE, ZERO, WORK(M+1,1),IONE)
                CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K+1),IONE, ZERO, WORK(1,2),IONE)
                CALL DGEMV("N",M,M,SGN,C,LDA, X(1,K+1),IONE, ZERO, WORK(M+1,2),IONE)

                DO L = K+2,N
                    CALL DAXPY(M,-B(K,L),WORK(1,1),IONE, X(1,L),IONE)
                    CALL DAXPY(M,-D(K,L),WORK(M+1,1),IONE, X(1,L),IONE)

                    CALL DAXPY(M,-B(K+1,L),WORK(1,2),IONE, X(1,L),IONE)
                    CALL DAXPY(M,-D(K+1,L),WORK(M+1,2),IONE, X(1,L),IONE)
                END DO

                K = K + 2
            END IF
        END DO

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN
        ! (T,N)
        K = 1
        DO WHILE ( K <= N)
            IF ( K .LT. N ) THEN
                LTMP = D(K+1,K) .EQ. ZERO .AND. B(K+1, K) .EQ. ZERO
            ELSE
                LTMP = .FALSE.
            END IF
            IF ( LTMP  .OR. (K .EQ. N)) THEN
                ! X(:,k)=(B(k,k)*A+D(k,k)*C)'\X(:,k);
                ! WORK(1:2*M,1:2*M) = ZERO
                DO I = 1, M
                    CALL DCOPY(M, A(1,I), IONE, WORK(1,I), IONE)
                    CALL DSCAL(M, B(K,K), WORK(1,I), IONE)
                    CALL DAXPY(M, SGN*D(K,K), C(1,I), IONE, WORK(1,I), IONE)
                END DO

                L = 1
                DO WHILE ( L < M )
                    IF ( WORK(L+1,L) .NE. ZERO ) THEN
                        AA = WORK(L,L)
                        BB = WORK(L+1,L)
                        CALL DROTG(AA,BB,WORK(2*L-1,2*M+1),WORK(2*L,2*M+1))
                        CALL DROT(M,WORK(L,L),2*M,WORK(L+1,L),2*M,WORK(2*L-1,2*M+1),WORK(2*L,2*M+1))
                        WORK(L+1,L) = ONE
                        L = L +2
                    ELSE
                        L = L +1
                    END IF
                END DO
                CALL DTRSV("U","T","N",M,WORK,2*M,X(1,K),IONE)

                L = 1
                DO WHILE ( L < M )
                    IF ( WORK(L+1,L) .EQ. ONE ) THEN
                        CALL DROT(IONE,X(L,K),IONE,X(L+1,K),IONE,WORK(2*L-1,2*M+1),-WORK(2*L,2*M+1))
                        L = L +2
                    ELSE
                        L = L +1
                    END IF
                END DO

                !Update January 2021
                !Precompute A^T X_k and C^T X_k
                CALL DGEMV("T",M,M,ONE, A, LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("T",M,M,SGN, C, LDA, X(1,K),IONE, ZERO, WORK(1,2),IONE)
                DO L = K+1,N
                    CALL DAXPY(M,-B(K,L),WORK(1,1),IONE,X(1,L),IONE)
                    CALL DAXPY(M,-D(K,L),WORK(1,2),IONE,X(1,L),IONE)
                END DO
                K = K + 1
            ELSE
                ! WORK(1:2*M, 1:2*M+1) = ZERO
                ! TMP = [ (B(k,k)*A+D(k,k)*C) ,(B(k,k+1)*A+D(k,k+1)*C) ;
                !         D(k+1,k)*C, (B(k+1,k+1)*A+D(k+1,k+1)*C)]'\reshape(X(:,k:k+1),2*m,1);
                DO L = 1, M
                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L-1),ITWO)
                    CALL DSCAL(M,B(K,K),WORK(1,2*L-1),ITWO)
                    CALL DAXPY(M,SGN*D(K,K),C(1,L),IONE,WORK(1,2*L-1),ITWO)

                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L),ITWO)
                    CALL DSCAL(M,B(K,K+1),WORK(1,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K,K+1),C(1,L),IONE,WORK(1,2*L),ITWO)

                    IF ( B(K+1,K) .NE. ZERO ) THEN
                        CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,B(K+1,K),WORK(2,2*L-1),ITWO)
                    ELSE
                        CALL DCOPY(M,C(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,SGN*D(K+1,K),WORK(2,2*L-1),ITWO)
                    END IF

                    CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L),ITWO)
                    CALL DSCAL(M,B(K+1,K+1),WORK(2,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K+1,K+1),C(1,L),IONE,WORK(2,2*L),ITWO)
                END DO

                CALL DCOPY(M, X(1,K), IONE, WORK(1,2*M+1),ITWO)
                CALL DCOPY(M, X(1,K+1),IONE, WORK(2,2*M+1),ITWO)

                L = 1
                DO WHILE ( L <= 2*M  )
                    IF ( (L+3 <= 2*M) .AND. ( WORK(L+2,L) .NE. ZERO) ) THEN
                        MAT(1:4,1:4) = TRANSPOSE(WORK(L:L+3,L:L+3))
                        RHS(1:4) = WORK(L:L+3,2*M+1)
                        CALL DLA_SMALL_SOLVE4( 4, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L:L+3,2*M+1) = RHS(1:4)
                        IF ( IINFO .NE. 0 )  INFO = 100

                        IF ( L+4 .LE. 2*M ) THEN
                            CALL DGEMV("T",4, 2*M-L-3, -ONE, WORK(L,L+4),2*M,WORK(L,2*M+1), &
                                & IONE,ONE,WORK(L+4,2*M+1),IONE)
                        END IF
                        L = L + 4
                    ELSE
                        MAT(1:2,1:2) = TRANSPOSE(WORK(L:L+1,L:L+1))
                        RHS(1:2) = WORK(L:L+1,2*M+1)
                        CALL DLA_SMALL_SOLVE4( 2, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L:L+1,2*M+1) = RHS(1:2)

                        IF ( IINFO .NE. 0 )  INFO = 100
                        CALL DGEMV("T",2, 2*M-L-1, -ONE, WORK(L,L+2),2*M,WORK(L,2*M+1), &
                            & IONE,ONE,WORK(L+2,2*M+1),IONE)
                        L = L + 2
                    END IF
                END DO

                CALL DCOPY(M, WORK(1,2*M+1),ITWO, X(1,K), IONE)
                CALL DCOPY(M, WORK(2,2*M+1),ITWO, X(1,K+1),IONE)

                ! Precompute the Update January 2021
                CALL DGEMV("T",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("T",M,M,SGN,C,LDA, X(1,K),IONE, ZERO, WORK(M+1,1),IONE)
                CALL DGEMV("T",M,M,ONE,A,LDA, X(1,K+1),IONE, ZERO, WORK(1,2),IONE)
                CALL DGEMV("T",M,M,SGN,C,LDA, X(1,K+1),IONE, ZERO, WORK(M+1,2),IONE)

                DO L = K+2,N
                    CALL DAXPY(M,-B(K,L),WORK(1,1),IONE, X(1,L),IONE)
                    CALL DAXPY(M,-D(K,L),WORK(M+1,1),IONE, X(1,L),IONE)

                    CALL DAXPY(M,-B(K+1,L),WORK(1,2),IONE, X(1,L),IONE)
                    CALL DAXPY(M,-D(K+1,L),WORK(M+1,2),IONE, X(1,L),IONE)
                END DO
                K = K + 2
            END IF
        END DO
    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN
        ! (N,T)
        K = N
        DO WHILE ( K >= 1 )
            IF ( K .GT. 1 ) THEN
                LTMP = B(K,K-1) .EQ. ZERO .AND. D(K,K-1) .EQ. ZERO
            ELSE
                LTMP = .FALSE.
            END IF
            IF ( K .EQ. 1 .OR. LTMP ) THEN
                ! WORK(1:M,1:M) = ZERO
                DO I = 1, M
                    CALL DCOPY(M,A(1,I), IONE, WORK(1,I), IONE)
                    CALL DSCAL(M,B(K,K), WORK(1,I),IONE)
                    CALL DAXPY(M,SGN*D(K,K),C(1,I),IONE,WORK(1,I),IONE)
                END DO

                ! New one using givens rotations
                L = 1
                DO WHILE ( L < M )
                    IF ( WORK(L+1,L) .NE. ZERO ) THEN
                        AA = WORK(L,L)
                        BB = WORK(L+1,L)
                        CALL DROTG(AA,BB,CC,SS)
                        CALL DROT(M,WORK(L,L),2*M,WORK(L+1,L),2*M,CC,SS)
                        CALL DROT(IONE,X(L,K),IONE,X(L+1,K),IONE,CC,SS)
                        L = L +2
                    ELSE
                        L = L +1
                    END IF
                END DO
                CALL DTRSV("U","N","N",M,WORK,2*M,X(1,K),IONE)

                CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("N",M,M,SGN,C,LDA, X(1,K),IONE, ZERO, WORK(1,2),IONE)

                !Update January 2021
                DO L = 1,K-1
                    CALL DAXPY(M,-B(L,K),WORK(1,1),IONE,X(1,L),IONE)
                    CALL DAXPY(M,-D(L,K),WORK(1,2),IONE,X(1,L),IONE)
                END DO
                K = K - 1
            ELSE
                WORK(1:2*M, 1:2*M+1) = ZERO

                ! New Idea
                DO L = 1, M
                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L-1),ITWO)
                    CALL DSCAL(M,B(K-1,K-1),WORK(1,2*L-1),ITWO)
                    CALL DAXPY(M,SGN*D(K-1,K-1),C(1,L),IONE,WORK(1,2*L-1),ITWO)

                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L),ITWO)
                    CALL DSCAL(M,B(K-1,K),WORK(1,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K-1,K),C(1,L),IONE,WORK(1,2*L),ITWO)

                    IF( B(K,K-1) .NE. ZERO) THEN
                        CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,B(K,K-1),WORK(2,2*L-1),ITWO)
                    ELSE
                        CALL DCOPY(M,C(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,SGN*D(K,K-1),WORK(2,2*L-1),ITWO)
                    END IF

                    CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L),ITWO)
                    CALL DSCAL(M,B(K,K),WORK(2,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K,K),C(1,L),IONE,WORK(2,2*L),ITWO)

                END DO

                CALL DCOPY(M, X(1,K-1), IONE, WORK(1,2*M+1),ITWO)
                CALL DCOPY(M, X(1,K),IONE, WORK(2,2*M+1),ITWO)

                L = 2*M
                DO WHILE ( L >= 1  )
                    LTMP = .FALSE.
                    IF ( (L >= 4)) THEN
                        IF (  WORK(L,L-2) .NE. ZERO)  THEN
                            LTMP = .TRUE.
                        END IF
                    END IF

                    IF ( LTMP ) THEN
                        MAT(1:4,1:4) = WORK(L-3:L,L-3:L)
                        RHS(1:4) = WORK(L-3:L,2*M+1)
                        CALL DLA_SMALL_SOLVE4( 4, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L-3:L,2*M+1) = RHS(1:4)
                        IF ( IINFO .NE. 0 ) INFO = 100

                        CALL DGEMV("N",L-4, 4, -ONE, WORK(1,L-3),2*M,WORK(L-3,2*M+1), &
                            & IONE,ONE,WORK(1,2*M+1),IONE)
                        L = L - 4
                    ELSE
                        MAT(1:2,1:2) = WORK(L-1:L,L-1:L)
                        RHS(1:2) = WORK(L-1:L,2*M+1)

                        CALL DLA_SMALL_SOLVE4( 2, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L-1:L,2*M+1) = RHS(1:2)
                        IF ( IINFO .NE. 0 ) INFO = 100

                        CALL DGEMV("N",L-2,2, -ONE, WORK(1,L-1),2*M,WORK(L-1,2*M+1), &
                            & IONE,ONE,WORK(1,2*M+1),IONE)
                        L = L - 2
                    END IF
                END DO

                CALL DCOPY(M, WORK(1,2*M+1),ITWO, X(1,K-1), IONE)
                CALL DCOPY(M, WORK(2,2*M+1),ITWO, X(1,K),IONE)


                CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("N",M,M,SGN,C,LDA, X(1,K),IONE, ZERO, WORK(M+1,1),IONE)
                CALL DGEMV("N",M,M,ONE,A,LDA, X(1,K-1),IONE, ZERO, WORK(1,2),IONE)
                CALL DGEMV("N",M,M,SGN,C,LDA, X(1,K-1),IONE, ZERO, WORK(M+1,2),IONE)

                DO L = 1,K-2
                    CALL DAXPY(M,-B(L,K),  WORK(1,1),IONE,X(1,L),IONE)
                    CALL DAXPY(M,-D(L,K),  WORK(M+1,1),IONE, X(1,L),IONE)

                    CALL DAXPY(M,-B(L,K-1),WORK(1,2),IONE, X(1,L),IONE)
                    CALL DAXPY(M,-D(L,K-1),WORK(M+1,2),IONE, X(1,L),IONE)
                END DO

                K = K - 2
            END IF
        END DO

    ELSE IF ( BTRANSA .AND. BTRANSB ) THEN
        ! (T,T)
        K = N
        DO WHILE ( K >= 1 )
            IF ( K .GT. 1 ) THEN
                LTMP = B(K,K-1) .EQ. ZERO .AND. D(K,K-1) .EQ. ZERO
            ELSE
                LTMP = .FALSE.
            END IF

            IF ( K .EQ. 1 .OR. LTMP) THEN
                ! WORK(1:M,1:M) = ZERO
                DO I = 1, M
                    CALL DCOPY(M, A(1,I), IONE, WORK(1,I), IONE)
                    CALL DSCAL(M, B(K,K), WORK(1,I), IONE)
                    CALL DAXPY(M, SGN*D(K,K), C(1,I), IONE, WORK(1,I),IONE)
                END DO

                L = 1
                DO WHILE ( L < M )
                    IF ( WORK(L+1,L) .NE. ZERO ) THEN
                        AA = WORK(L,L)
                        BB = WORK(L+1,L)
                        CALL DROTG(AA,BB,WORK(2*L-1,2*M+1),WORK(2*L,2*M+1))
                        CALL DROT(M,WORK(L,L),2*M,WORK(L+1,L),2*M,WORK(2*L-1,2*M+1),WORK(2*L,2*M+1))
                        WORK(L+1,L) = ONE
                        L = L +2
                    ELSE
                        L = L +1
                    END IF
                END DO
                CALL DTRSV("U","T","N",M,WORK,2*M,X(1,K),IONE)

                L = 1
                DO WHILE ( L < M )
                    IF ( WORK(L+1,L) .EQ. ONE ) THEN
                        CALL DROT(IONE,X(L,K),IONE,X(L+1,K),IONE,WORK(2*L-1,2*M+1),-WORK(2*L,2*M+1))
                        L = L +2
                    ELSE
                        L = L +1
                    END IF
                END DO


                !Update January 2021
                CALL DGEMV("T",M,M,ONE, A, LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("T",M,M,SGN, C, LDA, X(1,K),IONE, ZERO, WORK(1,2),IONE)

                DO L = 1,K-1
                    CALL DAXPY(M,-B(L,K),WORK(1,1),IONE,X(1,L),IONE)
                    CALL DAXPY(M,-D(L,K),WORK(1,2),IONE,X(1,L),IONE)
                END DO
                K = K - 1
            ELSE
                ! WORK(1:2*M, 1:2*M+1) = ZERO

                ! New Idea
                DO L = 1, M
                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L-1),ITWO)
                    CALL DSCAL(M,B(K-1,K-1),WORK(1,2*L-1),ITWO)
                    CALL DAXPY(M,SGN*D(K-1,K-1),C(1,L),IONE,WORK(1,2*L-1),ITWO)

                    CALL DCOPY(M,A(1,L),IONE,WORK(1,2*L),ITWO)
                    CALL DSCAL(M,B(K,K-1),WORK(1,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K,K-1),C(1,L),IONE,WORK(1,2*L),ITWO)

                    IF ( B(K,K-1) .NE. ZERO ) THEN
                        CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,B(K-1,K),WORK(2,2*L-1),ITWO)
                    ELSE
                        CALL DCOPY(M,C(1,L),IONE,WORK(2,2*L-1),ITWO)
                        CALL DSCAL(M,SGN*D(K-1,K),WORK(2,2*L-1),ITWO)
                    END IF

                    CALL DCOPY(M,A(1,L),IONE,WORK(2,2*L),ITWO)
                    CALL DSCAL(M,B(K,K),WORK(2,2*L),ITWO)
                    CALL DAXPY(M,SGN*D(K,K),C(1,L),IONE,WORK(2,2*L),ITWO)

                END DO


                CALL DCOPY(M, X(1,K-1), IONE, WORK(1,2*M+1),ITWO)
                CALL DCOPY(M, X(1,K),IONE, WORK(2,2*M+1),ITWO)

                L = 1
                DO WHILE ( L <= 2*M  )
                    IF ( (L+3 <= 2*M) .AND. ( WORK(L+2,L) .NE. ZERO) ) THEN
                        MAT(1:4,1:4) = TRANSPOSE(WORK(L:L+3,L:L+3))
                        RHS(1:4) = WORK(L:L+3,2*M+1)
                        CALL DLA_SMALL_SOLVE4( 4, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L:L+3,2*M+1) = RHS(1:4)

                        IF ( IINFO .NE. 0 ) INFO = 100

                        IF ( L + 4 .LE. 2*M ) THEN
                            CALL DGEMV("T",4, 2*M-L-3, -ONE, WORK(L,L+4),2*M,WORK(L,2*M+1), &
                                & IONE,ONE,WORK(L+4,2*M+1),IONE)
                        END IF
                        L = L + 4
                    ELSE

                        MAT(1:2,1:2) = TRANSPOSE(WORK(L:L+1,L:L+1))
                        RHS(1:2) = WORK(L:L+1,2*M+1)
                        CALL DLA_SMALL_SOLVE4( 2, MAT, RHS, IINFO, MACHINE_PRECISION, SMALL_NUMBER)
                        WORK(L:L+1,2*M+1) = RHS(1:2)

                        IF ( IINFO .NE. 0 ) INFO = 100

                        CALL DGEMV("T",2, 2*M-L-1, -ONE, WORK(L,L+2),2*M,WORK(L,2*M+1), &
                            & IONE,ONE,WORK(L+2,2*M+1),IONE)
                        L = L + 2
                    END IF
                END DO


                CALL DCOPY(M, WORK(1,2*M+1),ITWO, X(1,K-1), IONE)
                CALL DCOPY(M, WORK(2,2*M+1),ITWO, X(1,K),IONE)

                ! Precompute the Update January 2021
                CALL DGEMV("T",M,M,ONE,A,LDA, X(1,K),IONE, ZERO, WORK(1,1),IONE)
                CALL DGEMV("T",M,M,SGN,C,LDA, X(1,K),IONE, ZERO, WORK(M+1,1),IONE)
                CALL DGEMV("T",M,M,ONE,A,LDA, X(1,K-1),IONE, ZERO, WORK(1,2),IONE)
                CALL DGEMV("T",M,M,SGN,C,LDA, X(1,K-1),IONE, ZERO, WORK(M+1,2),IONE)

                DO L = 1,K-2
                    CALL DAXPY(M,-B(L,K),  WORK(1,1),IONE,X(1,L),IONE)
                    CALL DAXPY(M,-D(L,K),  WORK(M+1,1),IONE, X(1,L),IONE)

                    CALL DAXPY(M,-B(L,K-1),WORK(1,2),IONE, X(1,L),IONE)
                    CALL DAXPY(M,-D(L,K-1),WORK(M+1,2),IONE, X(1,L),IONE)
                END DO
                K = K - 2
            END IF
        END DO
    END IF


END SUBROUTINE

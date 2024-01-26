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

!> \brief Solver for a 4x4 generalized coupled Sylvester equation (TRANSA = N, TRANSB = N)
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLA_TGCSYLV_KERNEL_44NN (SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, INFO)
!           DOUBLE PRECISION SGN1, SGN2, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TGCSYLV_KERNEL_44NN solves a generalized Sylvester equation of the following form
!>
!>    A * R  + SGN1 * L * B = SCALE * E                              (1)
!>    C * R  + SGN2 * L * D = SCALE * F                              (1)
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand sides (E,F)  and the solutions
!> (R,L) M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by DGGES from LAPACK. The algorithm is implemented using BLAS level 2
!> operations. Thereby the order of M and N is at most 4. Furthermore, for fast execution
!> the function does not check the input arguments.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SGN1
!> \verbatim
!>          SGN1 is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign in the first equation of the Sylvester equation.
!> \endverbatim
!>
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign in the second equation of the Sylvester equation.
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
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (LDE,N)
!>          On input, the matrix E contains the right hand side E.
!>          On output, the matrix E contains the solution R of Equation (1).
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
!>          F is DOUBLE PRECISION array, dimension (LDF,N)
!>          On input, the matrix F contains the right hand side F.
!>          On output, the matrix F contains the solution L of Equation (1).
!>          Right hand side F and the solution L are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F.  LDF >= max(1,M).
!> \endverbatim
!
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          SCALE is a scaling factor to prevent the overflow in the result.
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
SUBROUTINE DLA_TGCSYLV_KERNEL_44NN (SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    DOUBLE PRECISION SGN1, SGN2, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)


    ! Local Variables
    INTEGER K, L, KB, LB, KH, LH, KE, IT
    INTEGER DIMMAT, INFO1, IPIV(8), JPIV(8)
    DOUBLE PRECISION MAT(8,8), RHS(8)
    DOUBLE PRECISION SCAL

    DOUBLE PRECISION ZERO , ONE
    PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

    ! EXTERNAL FUNCTIONS
    EXTERNAL DGESC2
    EXTERNAL DGETC2X
    EXTERNAL DLA_SMALL_SOLVE8


    ! Quick Return
    SCALE = ONE
    IF ( M .EQ. 0 .OR. N .EQ. 0 ) THEN
        RETURN
    END IF

    ! Prepare
    INFO = 0

    ! ( N, N )
    K = M
    DO WHILE (K .GT. 0)
        KB = 1
        KH = K
        IF ( K .GT. 1 ) THEN
            IF ( A(K,K-1) .NE. ZERO ) THEN
                KB = 2
                K = K - 1
            END IF
        END IF
        KE = K - 1

        L = 1
        DO WHILE ( L .LE. N )
            LB = 1
            IF ( L .LT. N ) THEN
                IF ( ( B(L+1,L) .NE. ZERO .OR. D(L+1,L) .NE. ZERO ) ) THEN
                    LB = 2
                END IF
            END IF
            LH = L + LB - 1

            IT = 0
            DO
                ! Solve inner system
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT(1,1) = A(K,K)
                    MAT(1,2) = SGN1 * B(L,L)
                    MAT(2,1) = C(K,K)
                    MAT(2,2) = SGN2 * D(L,L)
                    RHS(1) = E(K,L)
                    RHS(2) = F(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO
                    MAT(1,1) = A(K,K)
                    MAT(1,3) = SGN1 * B(L,L)
                    MAT(1,4) = SGN1 * B(LH, L)

                    MAT(2,2) = A(K,K)
                    MAT(2,3) = SGN1 * B(L,LH)
                    MAT(2,4) = SGN1 * B(LH, LH)

                    MAT(3,1) = C(K,K)
                    MAT(3,3) = SGN2 * D(L,L)
                    MAT(3,4) = SGN2 * D(LH, L)

                    MAT(4,2) = C(K,K)
                    MAT(4,3) = SGN2 * D(L,LH)
                    MAT(4,4) = SGN2 * D(LH, LH)

                    RHS(1) = E(K, L)
                    RHS(2) = E(K, LH)
                    RHS(3) = F(K, L)
                    RHS(4) = F(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = A(K,K)
                    MAT(1,2) = A(K,KH)
                    MAT(1,3) = SGN1 * B(L,L)

                    MAT(2,1) = A(KH,K)
                    MAT(2,2) = A(KH, KH)
                    MAT(2,4) = SGN1 * B(L,L)

                    MAT(3,1) = C(K,K)
                    MAT(3,2) = C(K,KH)
                    MAT(3,3) = SGN2 * D(L,L)

                    MAT(4,2) = C(KH, KH)
                    MAT(4,4) = SGN2 * D(L,L)

                    RHS(1) = E(K, L)
                    RHS(2) = E(KH, L)
                    RHS(3) = F(K, L)
                    RHS(4) = F(KH, L)
                ELSE
                    DIMMAT = 8

                    MAT(1:8,1:8) = ZERO
                    MAT(1,1) = A(K,K)
                    MAT(1,2) = A(K,KH)
                    MAT(1,5) = SGN1 * B(L,L)
                    MAT(1,7) = SGN1 * B(LH, L)

                    MAT(2,1) = A(KH,K)
                    MAT(2,2) = A(KH,KH)
                    MAT(2,6) = SGN1 * B(L,L)
                    MAT(2,8) = SGN1 * B(LH, L)

                    MAT(3,3) = A(K,K)
                    MAT(3,4) = A(K,KH)
                    MAT(3,5) = SGN1 * B(L,LH)
                    MAT(3,7) = SGN1 * B(LH,LH)

                    MAT(4,3) = A(KH,K)
                    MAT(4,4) = A(KH,KH)
                    MAT(4,6) = SGN1 * B(L,LH)
                    MAT(4,8) = SGN1 * B(LH,LH)

                    MAT(5,1) = C(K,K)
                    MAT(5,2) = C(K,KH)
                    MAT(5,5) = SGN2 * D(L,L)
                    MAT(5,7) = SGN2 * D(LH, L)

                    MAT(6,2) = C(KH,KH)
                    MAT(6,6) = SGN2 * D(L,L)
                    MAT(6,8) = SGN2 * D(LH, L)

                    MAT(7,3) = C(K,K)
                    MAT(7,4) = C(K,KH)
                    MAT(7,5) = SGN2 * D(L,LH)
                    MAT(7,7) = SGN2 * D(LH,LH)

                    MAT(8,4) = C(KH,KH)
                    MAT(8,6) = SGN2 * D(L,LH)
                    MAT(8,8) = SGN2 * D(LH,LH)


                    RHS(1) = E(K,L)
                    RHS(2) = E(KH,L)
                    RHS(3) = E(K,LH)
                    RHS(4) = E(KH,LH)

                    RHS(5) = F(K,L)
                    RHS(6) = F(KH,L)
                    RHS(7) = F(K,LH)
                    RHS(8) = F(KH,LH)
                END IF

                IF ( IT .EQ. 0 ) THEN
                    CALL DLA_SMALL_SOLVE8( DIMMAT, MAT, RHS, INFO1, MACHINE_PRECISION, SMALL_NUMBER)
                    IF ( INFO1 .EQ. 0 ) EXIT
                    IT = 1
                ELSE
                    CALL DGETC2X(DIMMAT, MAT, 8, IPIV, JPIV, INFO1)
                    IF ( INFO1 .NE. 0 ) INFO = INFO1
                    CALL DGESC2(DIMMAT, MAT, 8, RHS, IPIV, JPIV, SCAL)
                    IF ( SCAL .NE. ONE ) THEN
                        E(1:M, 1:N) = SCAL*E(1:M, 1:N)
                        F(1:M, 1:N) = SCAL*F(1:M, 1:N)

                        SCALE = SCALE * SCAL
                        INFO = 1
                    END IF
                    EXIT
                END IF

            END DO

            IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                E(K,L) = RHS(1)
                F(K,L) = RHS(2)

            ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN

                E(K, L)  = RHS(1)
                E(K, LH)  = RHS(2)
                F(K, L)  = RHS(3)
                F(K, LH)  = RHS(4)

            ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                E(K, L)  = RHS(1)
                E(KH, L)  = RHS(2)
                F(K, L)  = RHS(3)
                F(KH, L)  = RHS(4)
            ELSE
                E(K,L) = RHS(1)
                E(KH,L)  = RHS(2)
                E(K,LH)  = RHS(3)
                E(KH,LH)  = RHS(4)

                F(K,L) = RHS(5)
                F(KH,L)  = RHS(6)
                F(K,LH)  = RHS(7)
                F(KH,LH)  = RHS(8)
            END IF


            ! Update January 2021
            ! in current row
            IF ( LH .LT. N ) THEN
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    E(K,LH+1:N) = E(K,LH+1:N) - SGN1 *  F(K,L) * B(L,LH+1:N)
                    F(K,LH+1:N) = F(K,LH+1:N) - SGN2 *  F(K,L) * D(L,LH+1:N)
                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    E(K,LH+1:N) = E(K,LH+1:N) - SGN1* F(K,L) * B(L,LH+1:N) - SGN1 * F(K,LH) * B(LH,LH+1:N)
                    F(K,LH+1:N) = F(K,LH+1:N) - SGN2* F(K,L) * D(L,LH+1:N) - SGN2 * F(K,LH) * D(LH,LH+1:N)
                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    E(K,LH+1:N)  = E(K,LH+1:N)  - SGN1 * F(K,L) * B(L,LH+1:N)
                    E(KH,LH+1:N) = E(KH,LH+1:N) - SGN1 * F(KH,L) * B(L,LH+1:N)
                    F(K,LH+1:N)  = F(K,LH+1:N)  - SGN2 * F(K,L) *  D(L,LH+1:N)
                    F(KH,LH+1:N) = F(KH,LH+1:N) - SGN2 * F(KH,L) * D(L,LH+1:N)
                ELSE
                    E(K,LH+1:N)   = E(K,LH+1:N)   - SGN1 * F(K,L)  * B(L,LH+1:N) - SGN1 * F(K,LH)  * B(LH,LH+1:N)
                    E(KH,LH+1:N)  = E(KH,LH+1:N)  - SGN1 * F(KH,L) * B(L,LH+1:N) - SGN1 * F(KH,LH) * B(LH,LH+1:N)

                    F(K,LH+1:N)   = F(K,LH+1:N)   - SGN2 * F(K,L)  * D(L,LH+1:N) - SGN2 * F(K,LH)  * D(LH,LH+1:N)
                    F(KH,LH+1:N)  = F(KH,LH+1:N)  - SGN2 * F(KH,L) * D(L,LH+1:N) - SGN2 * F(KH,LH) * D(LH,LH+1:N)
                END IF

            END IF

            L = LH + 1
        END DO

        ! Update January 2021
        IF ( K .GT. 1 ) THEN
            IF ( KB .EQ. 1 ) THEN
                SELECT CASE(N)
                CASE(1)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1)


                CASE(2)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1)
                    E(1:KE,2) = E(1:KE, 2) -  A(1:KE,K) * E(K,2)
                    F(1:KE,2) = F(1:KE, 2) -  C(1:KE,K) * E(K,2)

                CASE(3)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1)
                    E(1:KE,2) = E(1:KE, 2) -  A(1:KE,K) * E(K,2)
                    F(1:KE,2) = F(1:KE, 2) -  C(1:KE,K) * E(K,2)
                    E(1:KE,3) = E(1:KE, 3) -  A(1:KE,K) * E(K,3)
                    F(1:KE,3) = F(1:KE, 3) -  C(1:KE,K) * E(K,3)

                CASE(4)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1)
                    E(1:KE,2) = E(1:KE, 2) -  A(1:KE,K) * E(K,2)
                    F(1:KE,2) = F(1:KE, 2) -  C(1:KE,K) * E(K,2)
                    E(1:KE,3) = E(1:KE, 3) -  A(1:KE,K) * E(K,3)
                    F(1:KE,3) = F(1:KE, 3) -  C(1:KE,K) * E(K,3)
                    E(1:KE,4) = E(1:KE, 4) -  A(1:KE,K) * E(K,4)
                    F(1:KE,4) = F(1:KE, 4) -  C(1:KE,K) * E(K,4)


                END SELECT
            ELSE
                SELECT CASE(N)
                CASE(1)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1) - A(1:KE,KH) * E(KH,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1) - C(1:KE,KH) * E(KH,1)
                CASE(2)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1) - A(1:KE,KH) * E(KH,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1) - C(1:KE,KH) * E(KH,1)
                    E(1:KE,2) = E(1:KE, 2) -  A(1:KE,K) * E(K,2) - A(1:KE,KH) * E(KH,2)
                    F(1:KE,2) = F(1:KE, 2) -  C(1:KE,K) * E(K,2) - C(1:KE,KH) * E(KH,2)

                CASE(3)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1) - A(1:KE,KH) * E(KH,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1) - C(1:KE,KH) * E(KH,1)
                    E(1:KE,2) = E(1:KE, 2) -  A(1:KE,K) * E(K,2) - A(1:KE,KH) * E(KH,2)
                    F(1:KE,2) = F(1:KE, 2) -  C(1:KE,K) * E(K,2) - C(1:KE,KH) * E(KH,2)
                    E(1:KE,3) = E(1:KE, 3) -  A(1:KE,K) * E(K,3) - A(1:KE,KH) * E(KH,3)
                    F(1:KE,3) = F(1:KE, 3) -  C(1:KE,K) * E(K,3) - C(1:KE,KH) * E(KH,3)

                CASE(4)
                    E(1:KE,1) = E(1:KE, 1) -  A(1:KE,K) * E(K,1) - A(1:KE,KH) * E(KH,1)
                    F(1:KE,1) = F(1:KE, 1) -  C(1:KE,K) * E(K,1) - C(1:KE,KH) * E(KH,1)
                    E(1:KE,2) = E(1:KE, 2) -  A(1:KE,K) * E(K,2) - A(1:KE,KH) * E(KH,2)
                    F(1:KE,2) = F(1:KE, 2) -  C(1:KE,K) * E(K,2) - C(1:KE,KH) * E(KH,2)
                    E(1:KE,3) = E(1:KE, 3) -  A(1:KE,K) * E(K,3) - A(1:KE,KH) * E(KH,3)
                    F(1:KE,3) = F(1:KE, 3) -  C(1:KE,K) * E(K,3) - C(1:KE,KH) * E(KH,3)
                    E(1:KE,4) = E(1:KE, 4) -  A(1:KE,K) * E(K,4) - A(1:KE,KH) * E(KH,4)
                    F(1:KE,4) = F(1:KE, 4) -  C(1:KE,K) * E(K,4) - C(1:KE,KH) * E(KH,4)
                END SELECT
            END IF
        END IF

        K = K - 1
    END DO

END SUBROUTINE




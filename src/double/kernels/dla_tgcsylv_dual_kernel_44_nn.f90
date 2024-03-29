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

!> \brief Solver for a 4x4 dual generalized coupled Sylvester equation (TRANSA = N, TRANSB = N)
!
!  Definition:
!  ===========
!
!    SUBROUTINE DLA_TGCSYLV_DUAL_KERNEL_44NN (SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, INFO)
!           DOUBLE PRECISION SGN1, SGN2, SCALE
!           INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
!           DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_TGCSYLV_DUAL_KERNEL_44NN solves a generalized coupled  Sylvester equation of the following form
!>
!>    op1(A)**T * R  + op1(C)**T * L               = SCALE * E                             (1)
!>    SGN1 * R * op2(B)**T + SGN2 * L * op2(D)** T = SCALE * F
!>
!> where A is a M-by-M quasi upper triangular matrix, C is a M-by-M upper triangular
!> matrix and B and D are N-by-N upper triangular matrices. Furthermore either B or
!> D can be quasi upper triangular as well. The right hand side Y and the solution X
!> M-by-N matrices. Typically the matrix pencils (A,C) and (B,D) ( or (D,B) ) are
!> created by DGGES from LAPACK.
!> The equation (1) is the dual to the generalized coupled Sylvester equation
!>
!>    op1(A) * R + SGN1 * L * op2(B)  = SCALE * E                                          (2)
!>    op1(C) * R + SGN2 * L * op2(D)  = SCALE * F
!>
!> The equation (1) is the dual one to equation (2) with respect to the underlying linear system.
!> Let Z be the matrix formed by rewriting (2) into its Kronecker form. This yields
!>
!>         | kron(I, op1(A))   SGN1*kron(op2(B)**T, I) | | Vec R |   | Vec E |
!>   Z X = |                                           |*|       | = |       |
!>         | kron(I, op1(C))   SGN2*kron(op2(D)**T, I) | | Vec L |   | Vec F |
!>
!> Regarding Z**T one obtains
!>
!>            | kron(I, op1(A)**T )    kron(I, op1(C)**T)   | | Vec R |   | Vec E |
!>   Z**T X = |                                             |*|       | = |       |
!>            | SGN1*kron(op2(B), I)   SGN2*kron(op2(D), I) | | Vec L |   | Vec F |
!>
!> which belongs to the Sylvester equation (1). For this reason the parameters TRANSA and TRANSB
!> are expressed in terms of the Sylvester equation (2).
!>
!> \endverbatim
!>
!> \attention This routine only allows to solve equations that are at most of size 4 x 4. The operator arguments TRANSA and TRANSB
!>            are fixed to op1(H) = H and op2(H) = H. The primary use of this function is to end the recursion in
!>            \ref dla_tgcsylv_dual_recursive.
!
!  Arguments:
!  ==========
!
!>
!> \param[in] SGN1
!> \verbatim
!>          SGN1 is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between in the first equation.
!> \endverbatim
!>
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between in the second equation.
!> \endverbatim
!>
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  128 >= M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices B and D.  128 >= N >= 0.
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
!>          F is DOUBLE PRECISION array, dimension (LDF,N)
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
!>          SCALE is DOUBLE PRECISION
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[inout] INFO
!> \verbatim
!>          INFO is INTEGER
!>
!>          On exit:
!>            == 0:  successful exit
!>            > 0:  The equation is not solved correctly. One of the arising inner
!>                  system got singular.
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup dbltgsylv
!
SUBROUTINE DLA_TGCSYLV_DUAL_KERNEL_44NN (SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, SCALE, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE
    ! Arguments
    DOUBLE PRECISION SGN1, SGN2, SCALE
    INTEGER M, N, LDA, LDB, LDC, LDD, LDE, LDF, INFO
    DOUBLE PRECISION A(LDA, *), B(LDB, *) , C(LDC, *), D(LDD, *), E(LDE, *), F(LDF,*)


    ! Local Variables
    INTEGER K, L, KB, LB, KH, LH, LE, IT
    INTEGER DIMMAT, INFO1, IPIV(8), JPIV(8)
    DOUBLE PRECISION MAT(8,8), RHS(8)
    DOUBLE PRECISION SCAL
    DOUBLE PRECISION AL(4,4), CL(4,4)

    INTEGER IONE
    DOUBLE PRECISION ZERO , ONE
    PARAMETER (IONE = 1, ZERO = 0.0D0, ONE = 1.0D0)

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

    DO L = 1, M
        DO K = MAX(1,L-1), M
            AL(K,L) = A(L,K)
            CL(K,L) = C(L,K)
        END DO
    END DO

    L = N
    DO WHILE ( L .GT. 0 )
        LB = 1
        LH = L
        IF ( L .GT. 1 ) THEN
            IF (B(L,L-1) .NE. ZERO) THEN
                LB = 2
                L = L - 1
            ELSE IF ( D(L,L-1) .NE. ZERO)  THEN
                LB = 2
                L = L - 1
            END IF
        END IF
        LE = L - 1


        K = 1
        DO WHILE ( K .LE. M )
            KB = 1
            IF ( K .LT. M ) THEN
                IF ( ( AL(K,K+1) .NE. ZERO )) THEN
                    KB = 2
                END IF
            END IF
            KH = K + KB - 1

            IT = 0
            DO
                ! Solve inner system
                IF ( KB .EQ. 1 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 2
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = CL(K,K)
                    MAT(2,1) = SGN1 * B(L,L)
                    MAT(2,2) = SGN2 * D(L,L)
                    RHS(1) = E(K,L)
                    RHS(2) = F(K,L)

                ELSE IF ( KB .EQ. 1 .AND. LB .EQ. 2 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,3) = CL(K,K)

                    MAT(2,2) = AL(K,K)
                    MAT(2,4) = CL(K,K)

                    MAT(3,1) = SGN1 * B(L,L)
                    MAT(3,2) = SGN1 * B(L, LH)
                    MAT(3,3) = SGN2 * D(L,L)
                    MAT(3,4) = SGN2 * D(L, LH)

                    MAT(4,1) = SGN1 * B(LH,L)
                    MAT(4,2) = SGN1 * B(LH, LH)
                    MAT(4,3) = SGN2 * D(LH,L)
                    MAT(4,4) = SGN2 * D(LH, LH)

                    RHS(1) = E(K, L)
                    RHS(2) = E(K, LH)
                    RHS(3) = F(K, L)
                    RHS(4) = F(K, LH)

                ELSE IF ( KB .EQ. 2 .AND. LB .EQ. 1 ) THEN
                    DIMMAT = 4
                    MAT(1:4,1:4) = ZERO

                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,3) = CL(K,K)
                    MAT(1,4) = CL(K,KH)

                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,3) = CL(KH,K)
                    MAT(2,4) = CL(KH,KH)

                    MAT(3,1) = SGN1 * B(L,L)
                    MAT(3,3) =  SGN2 * D(L,L)
                    MAT(4,2) = SGN1 * B(L,L)
                    MAT(4,4) =  SGN2 * D(L,L)

                    RHS(1) = E(K, L)
                    RHS(2) = E(KH, L)
                    RHS(3) = F(K, L)
                    RHS(4) = F(KH, L)
                ELSE
                    DIMMAT = 8

                    MAT(1:8,1:8) = ZERO
                    MAT(1,1) = AL(K,K)
                    MAT(1,2) = AL(K,KH)
                    MAT(1,5) = CL(K,K)
                    MAT(1,6) = CL(K,KH)
                    MAT(2,1) = AL(KH,K)
                    MAT(2,2) = AL(KH,KH)
                    MAT(2,5) = CL(KH,K)
                    MAT(2,6) = CL(KH,KH)

                    MAT(3,3) = AL(K,K)
                    MAT(3,4) = AL(K,KH)
                    MAT(3,7) = CL(K,K)
                    MAT(3,8) = CL(K,KH)
                    MAT(4,3) = AL(KH,K)
                    MAT(4,4) = AL(KH,KH)
                    MAT(4,7) = CL(KH,K)
                    MAT(4,8) = CL(KH,KH)

                    MAT(5,1) = SGN1 * B(L,L)
                    MAT(5,3) = SGN1 * B(L,LH)
                    MAT(5,5) =  SGN2 * D(L,L)
                    MAT(5,7) =  SGN2 * D(L,LH)

                    MAT(6,2) = SGN1 * B(L,L)
                    MAT(6,4) = SGN1 * B(L,LH)
                    MAT(6,6) =  SGN2 * D(L,L)
                    MAT(6,8) =  SGN2 * D(L,LH)

                    MAT(7,1) = SGN1 * B(LH,L)
                    MAT(7,3) = SGN1 * B(LH,LH)
                    MAT(7,5) =  SGN2 * D(LH,L)
                    MAT(7,7) =  SGN2 * D(LH,LH)

                    MAT(8,2) = SGN1 * B(LH,L)
                    MAT(8,4) = SGN1 * B(LH,LH)
                    MAT(8,6) =  SGN2 * D(LH,L)
                    MAT(8,8) =  SGN2 * D(LH,LH)

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

            IF ( KH .LT. M ) THEN
                IF (KB .EQ. IONE .AND. LB .EQ. IONE) THEN
                    E(KH+1:M, L) = E(KH+1:M, L) - AL(KH+1:M,K) * E(K,L)
                    E(KH+1:M, L) = E(KH+1:M, L) - CL(KH+1:M,K) * F(K,L)

                ELSE IF ( KB .EQ. IONE .AND. LB .NE. IONE) THEN
                    E(KH+1:M, L)  = E(KH+1:M, L) -  AL(KH+1:M,K) * E(K,L)
                    E(KH+1:M, LH) = E(KH+1:M, LH) - AL(KH+1:M,K) * E(K,LH)

                    E(KH+1:M, L)  = E(KH+1:M, L)  - CL(KH+1:M,K) * F(K,L)
                    E(KH+1:M, LH) = E(KH+1:M, LH) - CL(KH+1:M,K) * F(K,LH)

                ELSE IF ( KB .NE. IONE .AND. LB .EQ. IONE) THEN
                    E(KH+1:M,L) = E(KH+1:M,L) - AL(KH+1:M,K)*E(K,L) - AL(KH+1:M,KH)*E(KH,L)
                    E(KH+1:M,L) = E(KH+1:M,L) - CL(KH+1:M,K)*F(K,L) - CL(KH+1:M,KH)*F(KH,L)
                ELSE
                    E(KH+1:M,L)  = E(KH+1:M,L)  - AL(KH+1:M,K) * E(K,L)  - AL(KH+1:M,KH) * E(KH,L)
                    E(KH+1:M,LH) = E(KH+1:M,LH) - AL(KH+1:M,K) * E(K,LH) - AL(KH+1:M,KH) * E(KH,LH)
                    E(KH+1:M,L)  = E(KH+1:M,L)  - CL(KH+1:M,K) * F(K,L)  - CL(KH+1:M,KH) * F(KH,L)
                    E(KH+1:M,LH) = E(KH+1:M,LH) - CL(KH+1:M,K) * F(K,LH) - CL(KH+1:M,KH) * F(KH,LH)
                END IF
            END IF

            K = KH + 1
        END DO

        ! Update January 2021
        IF ( L .GT. 1 ) THEN
            IF ( LB .EQ. IONE ) THEN
                DO IT = 1,LE
                    F(1:M, IT) = F(1:M,IT) - E(1:M,L) * (SGN1 * B(IT,L))
                    F(1:M, IT) = F(1:M,IT) - F(1:M,L) * (SGN2 * D(IT,L))
                END DO
            ELSE
                DO IT = 1, LE
                    F(1:M, IT) = F(1:M,IT) - E(1:M,L) * (SGN1 * B(IT,L)) - E(1:M,LH) * (SGN1 * B(IT,LH))
                    F(1:M, IT) = F(1:M,IT) - F(1:M,L) * (SGN2 * D(IT,L)) - F(1:M,LH) * (SGN2 * D(IT,LH))
                END DO
            ENDIF
        END IF

        L = L - 1
    END DO


END SUBROUTINE




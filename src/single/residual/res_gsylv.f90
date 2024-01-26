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

!> \brief Compute the relative residual for the generalized Sylvester equation.
!
!  Definition:
!  ===========
!       FUNCTION SLA_RESIDUAL_GSYLV(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, Y, LDY, SCALE)
!           CHARACTER TRANSA(1), TRANSB(1)
!           REAL SCALE, SGN
!           INTEGER M, N, LDA, LDB. LDC, LDD, LDX, LDY
!           REAL A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), X(LDX, *), Y(LDY, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Compute the relative residual of the generalized Sylvester equation
!>
!>  RelRes = || SCALE*Y - opA(A)*X*opB(B) - SGN * opA(C) * X opB(D) ) || / ( (||A||*||B||+||C||*||D||)*||X|| + SCALE * ||Y|| )  (1)
!>
!> The norms are evaluated in terms of the Frobenius norm.
!> \endverbatim
!> \remark Auxiliary memory is managed by the function itself.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A :
!>          == 'N':  opA(A) = A (No transpose for A)
!>          == 'T':  opA(A) = A**T (Transpose A)
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to B :
!>          == 'N':  opB(B) = B (No transpose for B)
!>          == 'T':  opB(B) = B**T (Transpose B)
!> \endverbatim
!
!> \param[in] SGN
!> \verbatim
!>          SGN is REAL, allowed values: +/-1
!>          Specifies the sign between the two parts of the Sylvester equation.
!> \endverbatim
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,M)
!>          The coefficient matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          The coefficient matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array A.  LDB >= max(1,N).
!> \endverbatim
!
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (LDC,M)
!>          The coefficient matrix C.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (LDD,N)
!>          The coefficient matrix D.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array A.  LDD >= max(1,N).
!> \endverbatim

!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (LDX,N)
!>          Solution of the Sylvester Equation.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,M).
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is REAL array, dimension (LDY,N)
!>          The right hand side of the Sylvester Equation.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of the array Y.  LDY >= max(1,M).
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is REAL
!>          SCALE is a scaling factor to prevent the overflow in the result.
!> \endverbatim
!>
!> \return
!> \verbatim
!>          The function returns the relative residual as defined in (1). If an error
!>          occurs a negative integer I is returned, signalizing that the parameter -I
!>          has a wrong/illegal value. If I = -1000, the auxiliary memory could not be allocated.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup residual
!
FUNCTION SLA_RESIDUAL_GSYLV(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, Y, LDY, SCALE)
    IMPLICIT NONE
    CHARACTER :: TRANSA, TRANSB
    INTEGER :: M, N, LDA, LDB, LDC, LDD, LDX, LDY
    REAL :: A(LDA, *), B(LDB, *), C(LDC,*), D(LDD, *), X(LDX, *), Y(LDY, *)
    REAL :: SCALE, SGN, SLA_RESIDUAL_GSYLV

    REAL :: SLANGE
    EXTERNAL :: SGEMM, SLANGE, LSAME, SLACPY, XERROR_HANDLER
    LOGICAL :: LSAME

    INTEGER :: INFO, ALLOC_STATUS
    REAL :: NRMA, NRMX, NRMY, NRMR, NRMB, NRMC, NRMD
    REAL, ALLOCATABLE :: WORK(:,:), WORK2(:,:)
    REAL, PARAMETER :: DONE = 1.0
    REAL, PARAMETER :: DZERO = 0.0
    REAL, PARAMETER :: MONE = -1.0
    REAL :: DUMMY(1)



    INFO = 0
    IF ( .NOT. LSAME(TRANSA, "N") .AND. .NOT. LSAME(TRANSA, "T" )) THEN
        INFO = -1
    ELSE IF ( .NOT. LSAME(TRANSB, "N") .AND. .NOT. LSAME(TRANSB, "T" )) THEN
        INFO = -2
    ELSE IF ( SGN.NE.DONE .AND. SGN.NE.MONE) THEN
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
    ELSE IF ( LDY .LT. MAX(1, M)) THEN
        INFO = -17
    END IF

    IF (INFO .NE. 0 ) THEN
        SLA_RESIDUAL_GSYLV = REAL(INFO)
        CALL XERROR_HANDLER("SLA_RESIDUAL_GSYLV", -INFO)
        RETURN
    END IF

    IF ( M .EQ. 0 .OR. N .EQ. 0) THEN
        SLA_RESIDUAL_GSYLV = 0.0
        RETURN
    END IF

    ALLOCATE(WORK(M,N), STAT = ALLOC_STATUS)
    IF ( ALLOC_STATUS .NE. 0) THEN
        SLA_RESIDUAL_GSYLV = REAL(-1000)
        RETURN
    END IF

    ALLOCATE(WORK2(M,N), STAT = ALLOC_STATUS)
    IF ( ALLOC_STATUS .NE. 0) THEN
        SLA_RESIDUAL_GSYLV = REAL(-1000)
        DEALLOCATE(WORK)
        RETURN
    END IF


    NRMA = SLANGE("F", M, M, A, LDA, WORK)
    NRMB = SLANGE("F", N, N, B, LDB, WORK)
    NRMC = SLANGE("F", M, M, C, LDC, WORK)
    NRMD = SLANGE("F", N, N, D, LDD, WORK)
    NRMX = SLANGE("F", M, N, X, LDX, WORK)
    NRMY = SLANGE("F", M, N, Y, LDY, WORK)

    CALL SLACPY("A", M, N, Y, LDY, WORK(1,1), M)
    CALL SGEMM(TRANSA, "N", M, N, M, DONE, A, LDA, X, LDX, DZERO, WORK2, M)
    CALL SGEMM("N", TRANSB, M, N, N, MONE , WORK2, M, B, LDB, SCALE, WORK, M)
    CALL SGEMM(TRANSA, "N", M, N, M, DONE, C, LDC, X, LDX, DZERO, WORK2, M)
    CALL SGEMM("N", TRANSB, M, N, N, SGN*MONE , WORK2, M, D, LDD, DONE, WORK, M)

    NRMR = SLANGE("F", M, N, WORK(1,1), M, DUMMY)

    SLA_RESIDUAL_GSYLV = NRMR / ( (NRMA  * NRMB + NRMC*NRMD ) * NRMX + SCALE * NRMY)

    DEALLOCATE(WORK)
    DEALLOCATE(WORK2)

    RETURN

END FUNCTION SLA_RESIDUAL_GSYLV

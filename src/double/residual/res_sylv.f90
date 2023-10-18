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

!> \brief Compute the relative residual for the Sylvester equation.
!
!  Definition:
!  ===========
!       FUNCTION DLA_RESIDUAL_SYLV(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, Y, LDY, SCALE)
!           CHARACTER TRANSA(1), TRANSB(1)
!           DOUBLE PRECISION SCALE, SGN
!           INTEGER M, N, LDA, LDB. LDX, LDY
!           DOUBLE PRECISION A(LDA, *), B(LDB, *), X(LDX, *), Y(LDY, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Compute the relative residual of the Sylvester equation
!>
!>     RelRes = || SCALE*Y - opA(A)*X - SGN * X * opB(B) || / ( (||A||+||B||)*||X|| + SCALE * ||Y|| )     (1)
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
!>          SGN is DOUBLE PRECISION, allowed values: +/-1
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
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
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
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The coefficient matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array A.  LDB >= max(1,N).
!> \endverbatim
!
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
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
!>          Y is DOUBLE PRECISION array, dimension (LDY,N)
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
!>          SCALE is DOUBLE PRECISION
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
!> \date October 2023
!> \ingroup residual
!
FUNCTION DLA_RESIDUAL_SYLV(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, Y, LDY, SCALE)
    IMPLICIT NONE
    CHARACTER :: TRANSA, TRANSB
    INTEGER :: M, N, LDA, LDB, LDX, LDY
    DOUBLE PRECISION :: A(LDA, *), B(LDB, *), X(LDX, *), Y(LDY, *)
    DOUBLE PRECISION :: SCALE, SGN, DLA_RESIDUAL_SYLV

    DOUBLE PRECISION :: DLANGE
    EXTERNAL :: DGEMM, DLANGE, LSAME, DLACPY, XERROR_HANDLER
    LOGICAL :: LSAME

    INTEGER :: INFO, ALLOC_STATUS
    DOUBLE PRECISION :: NRMA, NRMX, NRMY, NRMR, NRMB
    DOUBLE PRECISION, ALLOCATABLE :: WORK(:,:)
    DOUBLE PRECISION, PARAMETER :: DONE = 1.0D0
    DOUBLE PRECISION, PARAMETER :: MONE = -1.0D0
    DOUBLE PRECISION :: DUMMY(1)



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
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -11
    ELSE IF ( LDY .LT. MAX(1, M)) THEN
        INFO = -13
    END IF

    IF (INFO .NE. 0 ) THEN
        DLA_RESIDUAL_SYLV = DBLE(INFO)
        CALL XERROR_HANDLER("DLA_RESIDUAL_SYLV", -INFO)
        RETURN
    END IF

    IF ( M .EQ. 0 .OR. N .EQ. 0) THEN
        DLA_RESIDUAL_SYLV = 0.0D0
        RETURN
    END IF

    ALLOCATE(WORK(M,N), STAT = ALLOC_STATUS)
    IF ( ALLOC_STATUS .NE. 0) THEN
        DLA_RESIDUAL_SYLV = DBLE(-1000)
        RETURN
    END IF


    NRMA = DLANGE("F", M, M, A, LDA, WORK)
    NRMB = DLANGE("F", N, N, B, LDB, WORK)
    NRMX = DLANGE("F", M, N, X, LDX, WORK)
    NRMY = DLANGE("F", M, N, Y, LDY, WORK)

    CALL DLACPY("A", M, N, Y, LDY, WORK(1,1), M)
    CALL DGEMM(TRANSA, "N", M, N, M, MONE, A, LDA, X, LDX, SCALE, WORK, M)
    CALL DGEMM("N", TRANSB, M, N, N, MONE * SGN , X, LDX, B, LDB, DONE, WORK, M)

    NRMR = DLANGE("F", M, N, WORK(1,1), M, DUMMY)

    DLA_RESIDUAL_SYLV = NRMR / ( (NRMA  + NRMB ) * NRMX + SCALE * NRMY)

    DEALLOCATE(WORK)

    RETURN

END FUNCTION DLA_RESIDUAL_SYLV

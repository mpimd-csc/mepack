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

!> \brief Compute the relative residual for the Lyapunov equation.
!
!  Definition:
!  ===========
!       FUNCTION DLA_RESIDUAL_LYAP(TRANS, M, A, LDA, X, LDX, Y, LDY, SCALE)
!           CHARACTER TRANS(!)
!           DOUBLE PRECISION SCALE
!           INTEGER M, LDA, LDX, LDY
!           DOUBLE PRECISION A(LDA, *), X(LDX, *), Y(LDY, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Compute the relative residual of the Lyapunov equation
!>
!>     RelRes = || SCALE*Y - op(A)*X - X * op(A)**T || / ( 2*||A||*||X|| + SCALE * ||Y|| )     (1)
!>
!> The norms are evaluated in terms of the Frobenius norm.
!> \endverbatim
!> \remark Auxiliary memory is managed by the function itself.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A :
!>          == 'N':  op(A) = A (No transpose for A)
!>          == 'T':  op(A) = A**T (Transpose A)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrix A.  M >= 0.
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
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,M)
!>          Solution of the Lyapunov Equation.
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
!>          Y is DOUBLE PRECISION array, dimension (LDY,M)
!>          The right hand side of the Lyapunov Equation.
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
FUNCTION DLA_RESIDUAL_LYAP(TRANS, M, A, LDA, X, LDX, Y, LDY, SCALE)
    IMPLICIT NONE
    CHARACTER :: TRANS
    INTEGER :: M, LDA, LDX, LDY
    DOUBLE PRECISION :: A(LDA, *), X(LDX, *), Y(LDY, *)
    DOUBLE PRECISION :: SCALE, DLA_RESIDUAL_LYAP

    DOUBLE PRECISION :: DLANGE
    EXTERNAL :: DGEMM, DLANGE, LSAME, DLACPY, XERROR_HANDLER
    LOGICAL :: LSAME

    INTEGER :: INFO, ALLOC_STATUS
    CHARACTER :: TRANS2
    DOUBLE PRECISION :: NRMA, NRMX, NRMY, NRMR
    DOUBLE PRECISION, ALLOCATABLE :: WORK(:,:)
    DOUBLE PRECISION, PARAMETER :: DONE = 1.0D0
    DOUBLE PRECISION, PARAMETER :: MONE = -1.0D0
    DOUBLE PRECISION :: DUMMY(1)



    INFO = 0
    IF ( .NOT. LSAME(TRANS, "N") .AND. .NOT. LSAME(TRANS, "T" )) THEN
        INFO = -1
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -2
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -4
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -6
    ELSE IF ( LDY .LT. MAX(1, M)) THEN
        INFO = -8
    END IF
    IF (INFO .NE. 0 ) THEN
        DLA_RESIDUAL_LYAP = DBLE(INFO)
        CALL XERROR_HANDLER("DLA_RESIDUAL_LYAP", -INFO)
        RETURN
    END IF

    IF ( M .EQ. 0 ) THEN
        DLA_RESIDUAL_LYAP = 0.0D0
        RETURN
    END IF

    ALLOCATE(WORK(M,M), STAT = ALLOC_STATUS)
    IF ( ALLOC_STATUS .NE. 0) THEN
        DLA_RESIDUAL_LYAP = DBLE(-1000)
        RETURN
    END IF


    NRMA = DLANGE("F", M, M, A, LDA, WORK)
    NRMX = DLANGE("F", M, M, X, LDX, WORK)
    NRMY = DLANGE("F", M, M, Y, LDY, WORK)

    IF (LSAME (TRANS, "N")) THEN
        TRANS2 = 'T'
    ELSE
        TRANS2 = 'N'
    END IF

    CALL DLACPY("A", M, M, Y, LDY, WORK(1,1), M)
    CALL DGEMM(TRANS, "N", M, M, M, MONE, A, LDA, X, LDX, SCALE, WORK, M)
    CALL DGEMM("N", TRANS2, M, M, M, MONE, X, LDX, A, LDA, DONE, WORK, M)

    NRMR = DLANGE("F", M, M, WORK(1,1), M, DUMMY)

    DLA_RESIDUAL_LYAP = NRMR / ( 2.0D0 * NRMA * NRMX + SCALE * NRMY)

    DEALLOCATE(WORK)

    RETURN

END FUNCTION DLA_RESIDUAL_LYAP

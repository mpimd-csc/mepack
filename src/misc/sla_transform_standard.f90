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

!> \brief Transform the Right-Hand Side or Solution for standard projections.
!
!  Definition:
!  ===========
!
!    SUBROUTINE SLA_TRANSFORM_STANDARD(TRANS, M, X, LDX, Q, LDQ, WORK)
!          CHARACTER(1) TRANS,
!          INTEGER M, LDQ, LDX
!          REAL Q(LDQ, *), X(LDX, *), WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SLA_TRANSFORM_STANDARD computes either
!>
!>    X <-  Q*X*Q**T                              (1)
!>
!> or
!>
!>    X <- Q**T*X*Q                               (2)
!>
!> where Q and X are M-by-M matrices.
!> \endverbatim
!> \remark The function does not perform any check on the input arguments.
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER(1)
!>          Specifies which transformation is applied:
!>          == 'N':  Transformation (1)
!>          == 'T':  Transformation (2)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices Q and X.  M >= 0.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,M)
!>          The matrix X is a general M-by-M matrix overwritten by (1) or (2)
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,M).
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,M)
!>          The matrix Q is M-by-M the transformation matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,M).
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is REAL array, dimension (M**2)
!>          Workspace for the computation.
!> \endverbatim
!>
!>
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup sglaux
!
SUBROUTINE SLA_TRANSFORM_STANDARD(TRANS, M, X, LDX, Q, LDQ, WORK)
    IMPLICIT NONE

    ! Arguments
    CHARACTER TRANS(*)
    INTEGER M, LDQ, LDX
    REAL Q(LDQ, *), X(LDX, *), WORK(*)

    ! Local Data
    REAL ONE, ZERO
    PARAMETER (ONE = 1.0, ZERO = 0.0)

    ! External Functions
    EXTERNAL LSAME, SGEMM
    LOGICAL LSAME


    IF (LSAME(TRANS, 'T')) THEN
        ! Transform Y <- Q**T * Y * Q
        CALL SGEMM("T", "N", M, M, M, ONE, Q, LDQ, X, LDX, ZERO, WORK, M)
        CALL SGEMM("N", "N", M, M, M, ONE, WORK, M, Q, LDQ, ZERO, X, LDX)
    ELSE
        ! Transform X <- Q * X * Q**T
        CALL SGEMM("N", "N", M, M, M, ONE, Q, LDQ, X, LDX, ZERO, WORK, M)
        CALL SGEMM("N", "T", M, M, M, ONE, WORK, M, Q, LDQ, ZERO, X, LDX)
    END IF
    RETURN
END SUBROUTINE

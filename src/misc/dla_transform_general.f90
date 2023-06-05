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

!> \brief Transform the Right-Hand Side or Solution for generalized projections.
!
!  Definition:
!  ===========
!
!    SUBROUTINE DLA_TRANSFORM_GENERAL(TRANS, M, N,  X, LDX, QA, LDQA, QB, LDQB,  WORK)
!          CHARACTER(1) TRANS,
!          INTEGER M, LDQ, LDX
!          DOUBLE PRECISION Q(LDQ, *), X(LDX, *), WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DLA_TRANSFORM_STANDARD computes either
!>
!>    X <-  QA*X*QB**T                              (1)
!>
!> or
!>
!>    X <- QA**T*X*QB                               (2)
!>
!> where QA is a M-by-M matrices, QB is a N-by-N matrix, and X is a M-by-N matrix.
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
!>          The order of the matrices QA and the number of rows of X.  M >= 0.
!> \endverbatim
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices QB and the number of columns of X.  N >= 0.
!> \endverbatim
!
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>          The matrix X is a general M-by-N matrix overwritten by (1) or (2)
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,M).
!> \endverbatim
!>
!> \param[in] QA
!> \verbatim
!>          QA is DOUBLE PRECISION array, dimension (LDQA,M)
!>          The matrix QA is M-by-M the transformation matrix.
!> \endverbatim
!>
!> \param[in] LDQA
!> \verbatim
!>          LDQA is INTEGER
!>          The leading dimension of the array QA.  LDQA >= max(1,M).
!> \endverbatim
!>
!> \param[in] QB
!> \verbatim
!>          QB is DOUBLE PRECISION array, dimension (LDQB,N)
!>          The matrix QB is N-by-N the transformation matrix.
!> \endverbatim
!>
!> \param[in] LDQB
!> \verbatim
!>          LDQB is INTEGER
!>          The leading dimension of the array QB.  LDQB >= max(1,N).
!> \endverbatim
!>
!> \param[in] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (M*N)
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
!> \date June 2023
!> \ingroup dblaux
!
SUBROUTINE DLA_TRANSFORM_GENERAL(TRANS, M, N, X, LDX, QA, LDQA, QB, LDQB, WORK )
    IMPLICIT NONE

    ! Arguments
    CHARACTER TRANS(*)
    INTEGER M, LDQA, LDQB, LDX, N
    DOUBLE PRECISION QA(LDQA, *), QB(LDQB, *), X(LDX, *), WORK(*)

    ! Local Data
    DOUBLE PRECISION ONE, ZERO
    PARAMETER (ONE = 1.0D0, ZERO = 0.0D0)

    ! External Functions
    EXTERNAL LSAME, DGEMM
    LOGICAL LSAME


    IF (LSAME(TRANS, 'T')) THEN
        ! Transform Y <- Q**T * Y * Q
        CALL DGEMM("T", "N", M, N, M, ONE, QA, LDQA, X, LDX, ZERO, WORK, M)
        CALL DGEMM("N", "N", M, N, N, ONE, WORK, M, QB, LDQB, ZERO, X, LDX)
    ELSE
        ! Transform X <- Q * X * Q**T
        CALL DGEMM("N", "N", M, N, M, ONE, QA, LDQA, X, LDX, ZERO, WORK, M)
        CALL DGEMM("N", "T", M, N, N, ONE, WORK, M, QB, LDQB, ZERO, X, LDX)
    END IF
    RETURN
END SUBROUTINE

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

!> \brief Block Reordering of generalized Eigenvalues.
!
!  Definition:
!  ===========
!       SUBROUTINE DLA_SORT_EV( N, A, LDA, Q, LDQ, NB, WORK, INFO )
!             IMPLICIT NONE
!             INTEGER INFO, LDA, LDQ, N, NB
!             DOUBLE PRECISION   A( LDA, * ), Q( LDQ, * ), WORK( * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>     Reorder the eigenvalues of a matrix A such that
!>     the complex eigenvalues pairs will not cause a change of the block
!>     size in blocked linear matrix equation solver codes.
!>
!>     The function uses DTREXC from LAPACK for moving the eigenvalues
!>     on the diagonal.
!>
!> \endverbatim
!>
!
!  Arguments:
!  ==========
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and Q, N >= 0.
!> \endverbatim
!
!> \param[in,out] A
!> \verbatim
!>             A  is a DOUBLE PRECISION array, dimension (LDA,N)
!>             On entry, the leading N-by-N upper
!>             Hessenberg part of this array must contain the
!>             generalized Schur factor A_s of the matrix A.
!>             On exit, the leading N-by-N part of this array contains
!>             the generalized Schur factor A_s of the matrix A with
!>             no complex eigenvalues pairs an A(k*NB,k*NB), where k is
!>             a non negative integer.
!> \endverbatim
!
!> \param[in] LDA
!> \verbatim
!>             LDA is INTEGER
!>             The leading dimension of the array A.  LDA >= MAX(1,N).
!> \endverbatim
!
!> \param[in,out] Q
!> \verbatim
!>             Q is a DOUBLE PRECISION array, dimension (LDQ,N)
!>             On entry, the leading N-by-N part of
!>             this array must contain the orthogonal matrix Q from
!>             the generalized Schur factorization.
!>             On exit, the leading N-by-N part of this array contains
!>             the update January 2021
!>             factorization with integrated eigenvalue reordering.
!> \endverbatim
!
!> \param[in] LDQ
!> \verbatim
!>             LDQ  is  INTEGER
!>             The leading dimension of the array Q.  LDQ >= MAX(1,N).
!> \endverbatim
!
!> \param[in] NB
!> \verbatim
!>            NB is INTEGER
!>             Block size of the solver planed to use. Typical values
!>             are 16, 32, or 64. The block size must be an even positive
!>             integer otherwise an error is returned.
!> \endverbatim
!
!  Workspace
!  =========
!
!> \param[in]  WORK
!> \verbatim
!>             WORK is INTEGER array, dimension (N)
!>             Workspace
!> \endverbatim
!
!   Error indicator
!   ===============
!
!> \param[out] INFO
!> \verbatim
!>             INFO is INTEGER
!>             == 0:  successful exit;
!>             < 0:  if INFO = -i, the i-th argument had an illegal
!>                   value;
!>             = 1:  The internal DTGEX2 caused an error. ]
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Januar 2023
!> \ingroup dblaux
!

SUBROUTINE DLA_SORT_EV( N, A, LDA, Q, LDQ, NB, WORK, INFO )
    IMPLICIT NONE

    ! .. Arguments ..
    INTEGER INFO, LDA, LDQ, N, NB
    DOUBLE PRECISION   A( LDA, * ), Q( LDQ, * ), WORK( * )

    ! .. Local Variables ..
    INTEGER K, POS, SJ, J, INFO2

    DOUBLE PRECISION ZERO
    PARAMETER (ZERO = 0.0D0)

    ! .. External Functions
    INTRINSIC MAX, MIN, MOD
    EXTERNAL XERBLA, DTREXC


    ! .. Check Parameters
    INFO = 0
    INFO2 = 0
    IF ( N .LT. 0) THEN
        INFO = -1
    ELSE IF ( LDA .LT. MAX(1, N))  THEN
        INFO = -3
    ELSE IF ( LDQ .LT. MAX(1, N))  THEN
        INFO = -5
    ELSE IF ( MOD(NB,2) .EQ. 1 ) THEN
        INFO = -6
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -6
    END IF

    IF ( INFO .NE. 0 ) THEN
        CALL XERBLA('DLA_SORT_EV', -INFO)
    END IF


    ! .. Sort the Eigenvalues on the Block Boundary..
    DO K = 1, N, NB
        POS = K+NB-1
        IF ( POS .GE. N ) THEN
            EXIT
        END IF

        IF ( A(POS+1, POS) .NE. ZERO) THEN
            SJ = K
            ! .. Search for the last previous real eigenvalue ..
            DO J = POS-2, K, -2
                IF( J .GT. 0 .AND. A(J+1,J) .EQ. ZERO) THEN
                    SJ = J + 1
                    EXIT
                END IF
            END DO

            ! .. Move the real eigenvalue to the blocksize boundary ..
            CALL DTREXC("Vector", N, A, LDA, Q, LDQ, SJ, POS, WORK, INFO2)
            IF ( INFO2 .NE. 0 ) THEN
                INFO = 1
                RETURN
            END IF
        END IF
    END DO


END SUBROUTINE

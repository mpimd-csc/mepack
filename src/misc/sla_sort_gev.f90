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
! Copyright (C) Martin Koehler, 2017-2022
!

!> \brief Block Reordering of generalized Eigenvalues.
!
!  Definition:
!  ===========
!       SUBROUTINE SLA_SORT_GEV( N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, NB, WORK, LWORK, INFO )
!             IMPLICIT NONE
!             INTEGER INFO, LDA, LDB, LDQ, LDZ, LWORK, N, NB
!             REAL   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>     Reorder the eigenvalues of a matrix pencil (A,B) such that
!>     the complex eigenvalues pairs will not cause a change of the block
!>     size in blocked linear matrix equation solver codes.
!>
!>     The function uses STGEX2 from LAPACK for moving the eigenvalues
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
!>          The order of the matrices A, B, Q, and Z, N >= 0.
!> \endverbatim
!
!> \param[in,out] A
!> \verbatim
!>             A  is a REAL array, dimension (LDA,N)
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
!> \param[in,out] B
!> \verbatim
!>             B is a REAL array, dimension (LDB,N)
!>             On entry, the leading N-by-N upper
!>             triangular part of this array must contain the
!>             generalized Schur factor B_s of the matrix B,
!>             On exit, the leading N-by-N part of this array contains
!>             the generalized Schur factor B_s of the matrix B
!>             with reordered eigenvalues.
!> \endverbatim
!
!> \param[in] LDB
!> \verbatim
!>             LDB is INTEGER
!>             The leading dimension of the array B.  LDB >= MAX(1,N).
!> \endverbatim
!
!> \param[in,out] Q
!> \verbatim
!>             Q is a REAL array, dimension (LDQ,N)
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
!> \param[in,out] Z
!> \verbatim
!>              Z is a REAL array, dimension (LDZ,N)
!>              On entry, the leading N-by-N part of
!>              this array must contain the orthogonal matrix Z from
!>              the generalized Schur factorization.
!>              On exit, the leading N-by-N part of this array contains
!>              the update January 2021
!>              factorization with integrated eigenvalue reordering.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>             LDZ is INTEGER
!>             The leading dimension of the array Z.  LDZ >= MAX(1,N).
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
!>             WORK is INTEGER array, dimension (LWORK)
!>             Workspace
!> \endverbatim
!
!> \param[in] LWORK
!> \verbatim
!>            LWORK  is INTEGER
!>            Size of the workspace.
!>            LWORK >= MAX(1, 4*N, 256)
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
!>             = 1:  The internal STGEX2 caused an error. ]
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup sglaux
!

SUBROUTINE SLA_SORT_GEV( N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, NB, WORK, LWORK, INFO )
    IMPLICIT NONE

    ! .. Arguments ..
    INTEGER INFO, LDA, LDB, LDQ, LDZ, LWORK, N, NB
    REAL   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )

    ! .. Local Variables ..
    INTEGER K, POS, SJ, J, INFO2

    REAL ZERO
    PARAMETER (ZERO = 0.0d0)
    LOGICAL WANT
    PARAMETER (WANT = .TRUE.)

    ! .. External Functions
    INTRINSIC MAX, MIN, MOD
    EXTERNAL XERBLA, STGEX2


    ! .. Check Parameters
    INFO = 0
    IF ( N .LT. 0) THEN
        INFO = -1
    ELSE IF ( LDA .LT. MAX(1, N))  THEN
        INFO = -3
    ELSE IF ( LDB .LT. MAX(1, N))  THEN
        INFO = -5
    ELSE IF ( LDQ .LT. MAX(1, N))  THEN
        INFO = -7
    ELSE IF ( LDZ .LT. MAX(1, N))  THEN
        INFO = -9
    ELSE IF ( MOD(NB,2) .EQ. 1 ) THEN
        INFO = -10
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -10
    ELSE IF ( LWORK .LT. MAX(1, N*4, 256)) THEN
        INFO = -12
    END IF

    IF ( INFO .NE. 0 ) THEN
        CALL XERBLA('SLA_SORT_GEV', -INFO)
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
            DO WHILE ( SJ .LT. POS )
                CALL STGEX2(WANT, WANT, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, SJ, 1, 2, WORK, LWORK, INFO2)
                SJ = SJ + 2
                IF ( INFO2 .NE. 0 ) THEN
                    INFO = 1
                    RETURN
                END IF
            END DO
        END IF
    END DO


END SUBROUTINE

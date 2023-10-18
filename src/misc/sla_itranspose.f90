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
!> \brief Solver for 4-by-4 Linear Systems
!
!  Definition:
!  ===========
!
!   SUBROUTINE SLA_ITRANSPOSE(K,KH,L,LH, X,LDX)
!            INTEGER K,KH,L,LH,LDX
!             REAL X(LDX,*)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_ITRANSPOSE performs the implicit transpose
!>
!>    X(L:LH,K:KH ) <- X(K:KH,L:LH)**T                                                   (1)
!>
!> \endverbatim
!>
!
!  Arguments:
!  ==========
!
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          First index of the transpose operations as in (1)
!> \endverbatim
!>
!> \param[in] KH
!> \verbatim
!>          KH is INTEGER
!>          Second index of the transpose operations as in (1)
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          Third index of the transpose operations as in (1)
!> \endverbatim
!>
!> \param[in] LH
!> \verbatim
!>          LH is INTEGER
!>          Fourth index of the transpose operations as in (1)
!> \endverbatim
!>
!> \param[inout] X
!> \verbatim
!>          X is REAL array, dimension (LDX,*)
!>          X is the matrix to perform the transpose on.
!> \endverbatim
!
!> \param[in] LDX
!> \verbatim
!>             LDX  is  INTEGER
!>             The leading dimension of the array X.  LDX >= MAX(1,KH,LH).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date October 2023
!> \ingroup sglaux
!
SUBROUTINE SLA_ITRANSPOSE(K,KH,L,LH, X,LDX)
    IMPLICIT NONE
    INTEGER K,KH,L,LH,LDX
    REAL X(LDX,*)

    INTEGER I, J

    !$omp parallel private(I) default(shared) if (KH-K+1 >= 192)
    !$omp do
    DO I = K,KH
        DO J = L, LH
            X(J, I) = X(I, J)
        END DO
        ! The vectorized statement creates a temporary array.
        ! X(L:LH,I) = X(I,L:LH)
    END DO
    !$omp end do
    !$omp end parallel

END SUBROUTINE


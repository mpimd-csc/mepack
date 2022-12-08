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
!> \brief Solver for 8-by-8 Linear Systems
!
!  Definition:
!  ===========
!       SUBROUTINE SLA_SMALL_SOLVE8( N, A, RHS, INFO, EPS, SMLNUM )
!             IMPLICIT NONE
!             INTEGER            INFO, N
!             REAL   A( 8, 8 )
!             REAL   RHS(8)
!             REAL   EPS
!             REAL   SMLNUM
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLA_SMALL_SOLVE8 solves a linear system
!>
!>    A * x = b                                                                          (1)
!>
!> where A is a most a 8 by 8 matrix with leading dimension 8.
!> \endverbatim
!>
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A  N = {1,2,4,8} .
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (8,8)
!>          The matrix A must be stored in a 4-by-4 matrix even if it is smaller.
!>          On output the matrix  A is destroyed.
!> \endverbatim
!>
!> \param[in,out] RHS
!> \verbatim
!>          RHS is REAL array, dimension (8)
!>          The matrix RHS contains the right hand side on input and the solution X on output if INFO == 0
!> \endverbatim
!>
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          == 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  The equation is not solved correctly. One of the arising inner
!>                system got singular.
!> \endverbatim
!>
!> \param[in] EPS
!> \verbatim
!>          EPS is REAL
!>          EPS is the machine epsilon determined by DLAMCH
!> \endverbatim
!>
!> \param[in] SMLNUM
!> \verbatim
!>          SMLNUM is REAL
!>          SMLNUM is the smallest number in REAL which can be represented without underflow.
!> \endverbatim
!>
!>
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup sglaux
!


SUBROUTINE SLA_SMALL_SOLVE8( N, A, RHS, INFO, EPS, SMLNUM )
    IMPLICIT NONE
    INTEGER            INFO, N
    REAL   A( 8, 8 )
    REAL   RHS(8)

    !     .. Parameters ..
    REAL   ZERO, ONE, TWO
    PARAMETER          ( ZERO = 0.0, ONE = 1.0, TWO = 2.0 )
    !     .. Local Scalars ..
    INTEGER            I, IP, IPV
    REAL   EPS, SMIN, SMLNUM, XMAX, TEMP, TEMPX(8)
    !dir$ attributes align: 64:: TEMPX


    !     .. Intrinsic Functions ..
    INTRINSIC          ABS, MAX, MAXVAL

    INFO = 0

    TEMPX(1:8) = ZERO
    IF (N.EQ.1) THEN
        XMAX = ABS(RHS(1))
        IF( TWO*SMLNUM*XMAX.GT.ABS( A( 1, 1 ) ) ) INFO = 1
        RHS(1) = RHS(1)/A(1,1)
        RETURN
    ELSE IF (N.EQ.2) THEN
        IPV = 1
        TEMPX(1) = ABS(A(1,1))
        TEMPX(2) = ABS(A(2,1))
        IF ( TEMPX(1) .LT. TEMPX(2) ) IPV = 2
        SMIN = MAX( EPS*TEMPX(IPV), SMLNUM )
        IF( IPV.NE.1 ) THEN
            TEMPX(1:2)  = A(IPV,1:2)
            A(IPV,1:2) = A(1,1:2)
            A(1,1:2) = TEMPX(1:2)
            TEMP = RHS(1)
            RHS(1) = RHS(2)
            RHS(2) = TEMP

        ENDIF

        ! Check for singularity
        IF( ABS( A( 1, 1 ) ).LT.SMIN ) THEN
            INFO = 1
            A( 1, 1 ) = SMIN
        END IF
        A( 2, 1 ) = A( 2, 1 ) / A( 1, 1 )
        A( 2, 2 ) = A( 2, 2) - A(2,1) * A(1,2)

        ! Solve
        RHS( 2 ) = RHS( 2 ) - A( 2, 1 )*RHS( 1 )

        XMAX = MAXVAL(ABS(RHS(1:2)),1)


        IF( TWO*SMLNUM*XMAX .GT. ABS( A( 2, 2 ) ) ) THEN
            INFO = 1
        END IF

        RHS( 2 ) = RHS( 2 )/A(2,2) ! *TEMP
        RHS( 1 ) = (RHS( 1 ) - RHS( 2 )*A( 1, 2 ))/A(1,1)

        RETURN
    ELSE IF ( N .EQ. 4 ) THEN
        IPV = 1
        SMIN = ZERO
        DO I = 1, 3
            !        Find max element in matrix A
            XMAX = ZERO
            DO  IP = I, 4
                TEMP = ABS(A(IP,I))
                IF( TEMP .GE. XMAX ) THEN
                    XMAX = TEMP
                    IPV = IP
                END IF
            END DO
            IF( I.EQ.1 )  SMIN = MAX( EPS*XMAX, SMLNUM )
            !        Swap rows
            IF( IPV.NE.I ) THEN
                TEMPX(1:4) = A(IPV,1:4)
                A(IPV,1:4) = A(I,1:4)
                A(I,1:4) = TEMPX(1:4)
                TEMP = RHS(IPV)
                RHS(IPV) = RHS(I)
                RHS(I) = TEMP
            ENDIF

            !        Check for singularity
            IF( ABS( A( I, I ) ).LT.SMIN ) THEN
                INFO = I
                A( I, I ) = SMIN
            END IF
            A(I+1:4,I) = A(I+1:4,I) / A(I,I)
            DO IP = I+1, 4
                A(I+1:N,IP) =  A(I+1:N,IP) - A(I+1:N,I)  * A(I,IP)
            END DO
        END DO

        ! Solve
        DO I = 2,4
            RHS(I:4) = RHS(I:4) - A(I:4,I-1) * RHS(I-1)
        END DO

        XMAX = MAXVAL(ABS(RHS(1:4)),1)

        IF( TWO*SMLNUM*XMAX .GT. ABS( A( 4, 4 ) ) ) THEN
            INFO = 1
        END IF

        DO I = 4,2, -1
            RHS(I) = RHS(I)/A(I,I)
            RHS(1:I-1) = RHS(1:I-1) - A(1:I-1,I) * RHS(I)
        END DO
        RHS(1) = RHS(1) / A(1,1)

        RETURN
    ELSE IF ( N.EQ.8 ) THEN
        IPV = 1
        SMIN = ZERO
        DO I = 1, 7
            !        Find max element in matrix A
            XMAX = ZERO
            DO  IP = I, 8
                TEMP = ABS(A(IP,I))
                IF( TEMP .GE. XMAX ) THEN
                    XMAX = TEMP
                    IPV = IP
                END IF
            END DO
            IF( I.EQ.1 )  SMIN = MAX( EPS*XMAX, SMLNUM )
            !        Swap rows
            IF( IPV.NE.I ) THEN
                TEMPX(1:8) = A(IPV,1:8)
                A(IPV,1:8) = A(I,1:8)
                A(I,1:8) = TEMPX(1:8)
                TEMP = RHS(IPV)
                RHS(IPV) = RHS(I)
                RHS(I) = TEMP
            ENDIF

            !        Check for singularity
            IF( ABS( A( I, I ) ).LT.SMIN ) THEN
                INFO = I
                A( I, I ) = SMIN
            END IF
            A(I+1:8,I) = A(I+1:8,I) / A(I,I)
            DO IP = I+1, 8
                IPV = I+1
                A(IPV:8,IP) =  A(IPV:8,IP) - A(IPV:8,I)  * A(I,IP)
            END DO
        END DO

        ! Solve
        DO I = 2,8
            RHS(I:8) = RHS(I:8) - A(I:8,I-1) * RHS(I-1)
        END DO

        XMAX = MAXVAL(ABS(RHS(1:8)),1)

        IF( TWO*SMLNUM*XMAX .GT. ABS( A( 4, 4 ) ) ) THEN
            INFO = 1
        END IF

        DO I = 8,2, -1
            RHS(I) = RHS(I)/A(I,I)
            RHS(1:I-1) = RHS(1:I-1) - A(1:I-1,I) * RHS(I)
        END DO
        RHS(1) = RHS(1) / A(1,1)

        RETURN
    ELSE
        INFO = 1
    END IF
    RETURN
END SUBROUTINE

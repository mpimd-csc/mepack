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
!       SUBROUTINE DLA_SMALL_SOLVE4( N, A, RHS, INFO, EPS, SMLNUM )
!             IMPLICIT NONE
!             INTEGER            INFO, N
!             DOUBLE PRECISION   A( 4, 4 )
!             DOUBLE PRECISION   RHS(4)
!             DOUBLE PRECISION   EPS
!             DOUBLE PRECISION   SMLNUM
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLA_SMALL_SOLVE4 solves a linear system
!>
!>    A * x = b                                                                          (1)
!>
!> where A is a most a 4 by 4 matrix.
!> \endverbatim
!>
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A  N = (1,2,4).
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (4,4)
!>          The matrix A must be stored in a 4-by-4 matrix even if it is smaller.
!>          On output the matrix  A is destroyed.
!> \endverbatim
!>
!> \param[in,out] RHS
!> \verbatim
!>          RHS is DOUBLE PRECISION array, dimension (4)
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
!>          EPS is DOUBLE PRECISION
!>          EPS is the machine epsilon determined by DLAMCH
!> \endverbatim
!>
!> \param[in] SMLNUM
!> \verbatim
!>          SMLNUM is DOUBLE PRECISION
!>          SMLNUM is the smallest number in DOUBLE PRECISION which can be represented without underflow.
!> \endverbatim
!>
!>
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date October 2023
!> \ingroup dblaux
!


SUBROUTINE DLA_SMALL_SOLVE4( N, A, RHS, INFO, EPS, SMLNUM )
    IMPLICIT NONE
    INTEGER            INFO, N
    DOUBLE PRECISION   A( 4, 4 )
    DOUBLE PRECISION   RHS(4)

    !     .. Parameters ..
    DOUBLE PRECISION   ZERO, ONE, TWO
    PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
    !     .. Local Scalars ..
    INTEGER            I, IP, IPV
    DOUBLE PRECISION   EPS, SMIN, SMLNUM, XMAX, TEMP, TEMPX(4), RHSX(4)
    !dir$ attributes align: 64:: TEMPX
    !dir$ attributes align: 64:: RHSX

    !     .. Intrinsic Functions ..
    INTRINSIC          ABS, MAX, MAXVAL

    INFO = 0

    TEMPX(1:4) = ZERO

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

        RHS( 2 ) = RHS( 2 )/A(2,2)
        RHS( 1 ) = (RHS( 1 ) - RHS( 2 )*A( 1, 2 ))/A(1,1)

        RETURN
    ELSE IF ( N .EQ. 3 ) THEN
        IPV = 1
        SMIN = ZERO
        DO I = 1, 2

            XMAX = ZERO
            DO  IP = I, 3
                TEMP = ABS(A(IP,I))
                IF( TEMP .GE. XMAX ) THEN
                    XMAX = TEMP
                    IPV = IP
                END IF
            END DO
            IF( I.EQ.1 )  SMIN = MAX( EPS*XMAX, SMLNUM )
            !        Swap rows
            IF( IPV.NE.I ) THEN
                TEMPX(1:3) = A(IPV,1:3)
                A(IPV,1:3) = A(I,1:3)
                A(I,1:3) = TEMPX(1:3)
                TEMP = RHS(IPV)
                RHS(IPV) = RHS(I)
                RHS(I) = TEMP
            ENDIF

            !        Check for singularity
            IF( ABS( A( I, I ) ).LT.SMIN ) THEN
                INFO = I
                A( I, I ) = SMIN
            END IF
            A(I+1:3,I) = A(I+1:3,I) / A(I,I)
            DO IP = I+1, 3
                A(I+1:N,IP) =  A(I+1:N,IP) - A(I+1:N,I)  * A(I,IP)
            END DO
        END DO

        DO I = 2,3
            RHS(I:3) = RHS(I:3) - A(I:3,I-1) * RHS(I-1)
        END DO

        XMAX = MAXVAL(ABS(RHS(1:3)),1)

        IF( TWO*SMLNUM*XMAX .GT. ABS( A( 3, 3 ) ) ) THEN
            INFO = 1
        END IF

        DO I = 3,2, -1
            RHS(I) = RHS(I)/A(I,I)
            RHS(1:I-1) = RHS(1:I-1) - A(1:I-1,I) * RHS(I)
        END DO
        RHS(1) = RHS(1) / A(1,1)

        RETURN

    ELSE IF ( N .EQ. 4 ) THEN
        IPV = 1
        SMIN = ZERO
        RHSX(1:4) = RHS(1:4)
        ! I = 1
        IPV = MAXLOC(ABS(A(1:4,1)),1)
        XMAX = ABS(A(IPV,1))
        SMIN = MAX( EPS*XMAX, SMLNUM )
        !        Swap rows
        IF( IPV.NE. 1 ) THEN
            TEMPX = A(IPV,1:4)
            A(IPV,1:4) = A(1,1:4)
            A(1,1:4) = TEMPX(1:4)
            TEMP = RHSX(IPV)
            RHSX(IPV) = RHSX(1)
            RHSX(1) = TEMP
        ENDIF

        !        Check for singularity
        IF( ABS( A( 1, 1 ) ).LT.SMIN ) THEN
            INFO = 1
            A( 1, 1 ) = SMIN
        END IF
        A(2:4,1) = A(2:4,1) / A(1,1)
        A(2:4,2) =  A(2:4,2) - A(2:4,1)  * A(1,2)
        A(2:4,3) =  A(2:4,3) - A(2:4,1)  * A(1,3)
        A(2:4,4) =  A(2:4,4) - A(2:4,1)  * A(1,4)

        ! I = 2
        IPV = MAXLOC(ABS(A(2:4,2)),1) + 1
        XMAX = ABS(A(IPV,2))
        !        Swap rows
        IF( IPV.NE.2 ) THEN
            TEMPX = A(IPV,1:4)
            A(IPV,1:4) = A(2,1:4)
            A(2,1:4) = TEMPX(1:4)
            TEMP = RHSX(IPV)
            RHSX(IPV) = RHSX(2)
            RHSX(2) = TEMP
        ENDIF

        !        Check for singularity
        IF( ABS( A( 2, 2 ) ).LT.SMIN ) THEN
            INFO = 2
            A( 2, 2 ) = SMIN
        END IF
        A(3:4,2) =  A(3:4,2) / A(2,2)
        A(3:4,3) =  A(3:4,3) - A(3:4,2)  * A(2,3)
        A(3:4,4) =  A(3:4,4) - A(3:4,2)  * A(2,4)

        ! I = 3
        IPV = MAXLOC(ABS(A(3:4,3)),1) + 2
        XMAX = ABS(A(IPV,3))
        !        Swap rows
        IF( IPV.NE.3 ) THEN
            TEMPX = A(IPV,1:4)
            A(IPV,1:4) = A(3,1:4)
            A(3,1:4) = TEMPX(1:4)
            TEMP = RHSX(IPV)
            RHSX(IPV) = RHSX(3)
            RHSX(3) = TEMP
        ENDIF

        !        Check for singularity
        IF( ABS( A( 3, 3 ) ).LT.SMIN ) THEN
            INFO = 3
            A( 3, 3 ) = SMIN
        END IF
        A(4,3) =  A(4,3) / A(3,3)
        A(4,4) =  A(4,4) - A(4,3)  * A(3,4)

        ! Forward / Backward
        RHSX(2:4) = RHSX(2:4) - A(2:4,1) * RHSX(1)
        RHSX(3:4) = RHSX(3:4) - A(3:4,2) * RHSX(2)
        RHSX(4) = RHSX(4) - A(4,3) * RHSX(3)

        XMAX = MAXVAL(ABS(RHSX(1:4)),1)

        IF( TWO*SMLNUM*XMAX .GT. ABS( A( 4, 4 ) ) ) THEN
            INFO = 1
        END IF

        RHSX(4)   = RHSX(4)/A(4,4)
        RHSX(1:3) = RHSX(1:3) - A(1:3,4) * RHSX(4)
        RHSX(3)   = RHSX(3)/A(3,3)
        RHSX(1:2) = RHSX(1:2) - A(1:2,3) * RHSX(3)
        RHSX(2)   = RHSX(2)/A(2,2)
        RHSX(1)   = RHSX(1) - A(1,2) * RHSX(2)
        RHSX(1)   = RHSX(1) / A(1,1)

        RHS(1:4) = RHSX(1:4)
        RETURN
    ELSE
        INFO = -1
    END IF
    RETURN
END SUBROUTINE

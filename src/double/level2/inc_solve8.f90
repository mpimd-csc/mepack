!
! Inline code to solve a 8x8 Linear system

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


INFO1 = 0
DO J = 1, 8
    JP = J-1+MAXLOC(ABS(MAT(J:8,J)),1)
    IF( MAT( JP, J ).NE. 0.0D0 ) THEN
        IF( JP.NE.J ) THEN
            TEMPX(1:12) = MAT(J,1:12)
            MAT(J,1:12) = MAT(JP,1:12)
            MAT(JP,1:12) = TEMPX(1:12)
        END IF

        IF( J.LT.8 ) THEN
!            !$omp simd safelen(8) simdlen(8)
            DO I = 1, 8-J
                MAT( J+I, J ) = MAT( J+I, J ) / MAT( J, J )
            END DO
!            !$omp end simd
        END IF
    ELSE
        INFO1= 1
    END IF

    IF( J.LT. 8 ) THEN
        DO I = J+1, 9
!            !$omp simd safelen(8) simdlen(8)
            DO IX = J+1, 8
                MAT(IX, I) = MAT(IX, I) - MAT(IX,J) * MAT(J,I)
            END DO
!            !$omp end simd
        END DO
    END IF
END DO

! J = 8
MAT(8,9) = MAT(8,9)/MAT(8,8)
MAT(1:7,9) = MAT(1:7,9) - MAT(8,9)*MAT(1:7,8)
! J = 7
MAT(7,9) = MAT(7,9)/MAT(7,7)
MAT(1:6,9) = MAT(1:6,9) - MAT(7,9)*MAT(1:6,7)
! J = 6
MAT(6,9) = MAT(6,9)/MAT(6,6)
MAT(1:5,9) = MAT(1:5,9) - MAT(6,9)*MAT(1:5,6)
! J = 5
MAT(5,9) = MAT(5,9)/MAT(5,5)
MAT(1:4,9) = MAT(1:4,9) - MAT(5,9)*MAT(1:4,5)
! J = 4
MAT(4,9) = MAT(4,9)/MAT(4,4)
MAT(1:3,9) = MAT(1:3,9) - MAT(4,9)*MAT(1:3,4)
! J = 3
MAT(3,9) = MAT(3,9)/MAT(3,3)
MAT(1:2,9) = MAT(1:2,9) - MAT(3,9)*MAT(1:2,3)
! J = 2
MAT(2,9) = MAT(2,9)/MAT(2,2)
MAT(1,9) = MAT(1,9) - MAT(2,9)*MAT(1,2)
! J = 1
MAT(1,9) = MAT(1,9)/MAT(1,1)


! DO J = 8,1,-1
!     IF (MAT(J,9).NE.0.0D0) THEN
!         MAT(J,9) = MAT(J,9)/MAT(J,J)
!         TEMP   = MAT(J,9)
!         ! MAT(1:J-1,9) = MAT(1:J-1,9) - TEMP*MAT(1:J-1,J)
!         !$omp simd safelen(8) simdlen(8)
!         DO I = 1, J-1
!             MAT(I, 9 )  = MAT(I, 9) - TEMP * MAT(I, J)
!         END DO
!         !$omp end simd
!     END IF
! END DO



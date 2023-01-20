!
! THIS PROGRAM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
! IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
! THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
! (AT YOUR OPTION) ANY LATER VERSION.
!
! THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
! BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
! MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
! GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!
! YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
! ALONG WITH THIS PROGRAM; IF NOT, SEE <HTTP://WWW.GNU.ORG/LICENSES/>.
!
! Copyright (C) Martin Koehler, 2017-2023
!

!> @addtogroup machine
!> @{

!> @brief Machine specific parameters for Single Precision Calculations
!>
!> The MEPACK_OPTIONS_MACHINE_SINGLE module contains the machine specific parameters for
!> single prevent calculations like the machine precision.
MODULE MEPACK_OPTIONS_MACHINE_SINGLE
    IMPLICIT NONE

    PRIVATE
    INTERFACE
        FUNCTION SLAMCH(WHAT)
            CHARACTER, DIMENSION(1), INTENT(IN) :: WHAT
            REAL :: SLAMCH
        END FUNCTION
        SUBROUTINE SLABAD(X, Y)
            REAL, INTENT(INOUT) :: X, Y
        END SUBROUTINE
    END INTERFACE

    REAL, PARAMETER :: ONE = 1.0


    !> @brief The machine precision in single precision.
    REAL, PUBLIC, SAVE :: MACHINE_PRECISION = 0.0
    !> @brief The smallest save single precision number.
    REAL, PUBLIC, SAVE :: SMALL_NUMBER = 0.0
    !> @brief The biggest save single precision number.
    REAL, PUBLIC, SAVE :: BIG_NUMBER = 0.0

    PUBLIC :: MEPACK_OPTIONS_MACHINE_SINGLE_INIT

    CONTAINS
        !> @brief Initialize the Double Precision Machine Specific Values
        !>
        !> The MEPACK_OPTIONS_MACHINE_DOUBLE_INIT subroutine initializes the machine specific
        !> parameters, like the machine precision, in single precision.
        !> The subroutine must be called once before the first MEPACK routine
        !> handling single precision data is called.
        SUBROUTINE MEPACK_OPTIONS_MACHINE_SINGLE_INIT()
            IMPLICIT NONE

            MACHINE_PRECISION = SLAMCH( 'E' )
            SMALL_NUMBER      = SLAMCH( 'S' ) / MACHINE_PRECISION
            BIG_NUMBER        = ONE / SMALL_NUMBER
            CALL SLABAD( SMALL_NUMBER, BIG_NUMBER )
        END SUBROUTINE
END MODULE MEPACK_OPTIONS_MACHINE_SINGLE

!> @}

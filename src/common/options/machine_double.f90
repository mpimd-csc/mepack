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
! Copyright (C) Martin Koehler, 2017-2022
!

!> @addtogroup machine
!> @{

!> @brief Machine specific parameters for Double Precision Calculations
!>
!> The MEPACK_OPTIONS_MACHINE_DOUBLE module contains the machine specific parameters for
!> double prevent calculations like the machine precision.
MODULE MEPACK_OPTIONS_MACHINE_DOUBLE
    IMPLICIT NONE

    PRIVATE
    INTERFACE
        FUNCTION DLAMCH(WHAT)
            CHARACTER, DIMENSION(1), INTENT(IN) :: WHAT
            DOUBLE PRECISION :: DLAMCH
        END FUNCTION
        SUBROUTINE DLABAD(X, Y)
            DOUBLE PRECISION, INTENT(INOUT) :: X, Y
        END SUBROUTINE
    END INTERFACE

    DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0


    !> @brief The machine precision in double precision.
    DOUBLE PRECISION, PUBLIC :: MACHINE_PRECISION = 0.0D0
    !> @brief The smallest save double precision number.
    DOUBLE PRECISION, PUBLIC :: SMALL_NUMBER = 0.0D0
    !> @brief The biggest save double precision number.
    DOUBLE PRECISION, PUBLIC :: BIG_NUMBER = 0.0D0

    PUBLIC :: MEPACK_OPTIONS_MACHINE_DOUBLE_INIT

    SAVE

    CONTAINS
        !> @brief Initialize the Double Precision Machine Specific Values
        !>
        !> The MEPACK_OPTIONS_MACHINE_DOUBLE_INIT subroutine initializes the machine specific
        !> parameters, like the machine precision, in double precision.
        !> The subroutine must be called once before the first MEPACK routine
        !> handling double precision data is called.
        SUBROUTINE MEPACK_OPTIONS_MACHINE_DOUBLE_INIT()
            IMPLICIT NONE

            MACHINE_PRECISION = DLAMCH( 'E' )
            SMALL_NUMBER      = DLAMCH( 'S' ) / MACHINE_PRECISION
            BIG_NUMBER        = ONE / SMALL_NUMBER
            CALL DLABAD( SMALL_NUMBER, BIG_NUMBER )
        END SUBROUTINE
END MODULE MEPACK_OPTIONS_MACHINE_DOUBLE
!> @}

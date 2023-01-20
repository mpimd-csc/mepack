!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) Martin Koehler, 2017-2023

!> @addtogroup auxoptions
!> @{

!> @brief Control the verbosity of MEPACK
!>
!> MEPACK can output some details about the used environment, for example which tuning file is loaded.
!> Therefore, the verbosity can be turned on or off using the environment variable \b MEPACK_VERBOSE or
!> the subroutines from this module. The verbosity does not effect any computations performed by MEPACK.
!>
MODULE MEPACK_OPTIONS_VERBOSE

    !> @internal
    INTEGER, PRIVATE ::  VERBOSITY_LEVEL = 0

CONTAINS

    !> @brief Initializes the module.
    !>
    !> The VERBOSE_INIT subroutine is automatically called when MEPACK is initialized. It reads
    !> the \b MEPACK_VERBOSE environment variable and set the internal verbosity state accordingly.
    !>
    !> @see mepack_options_verbose::verbose_set
    SUBROUTINE VERBOSE_INIT()
        IMPLICIT NONE

        CHARACTER(LEN=100) :: VALUE
        INTEGER :: LEN, STA, STATE

        CALL GET_ENVIRONMENT_VARIABLE(NAME="MEPACK_VERBOSE", VALUE=VALUE, LENGTH=LEN, STATUS=STA)

        IF ( STA .EQ. 0 ) THEN
            READ(VALUE, *) STATE
        ELSE
            STATE = 0
        END IF

        VERBOSITY_LEVEL = STATE
    END SUBROUTINE

    !> @brief Set the verbosity level of MEPACK
    !> @param[in]   LEVEL   New verbosity level.
    !>
    !> The VERBOSE_SET subroutine sets the verbosity level of MEPACK. At the moment this is either zero if
    !> MEPACK should be quiet or anything non-zero, which means MEPACK is showing up with some additional messages.
    SUBROUTINE VERBOSE_SET(LEVEL)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: LEVEL

        VERBOSITY_LEVEL = LEVEL
    END SUBROUTINE

    !> @brief Returns the current verbosity level of MEPACK.
    !> @return Integer value defining the current verbosity level.
    !>
    !> The VERBOSE_GET function returns the current verbosity level. Zero means
    !> nothing should be written to the screen. Everything non-zero is interpreted as print some information on the screen.
    FUNCTION VERBOSE_GET() RESULT(LEVEL)
        IMPLICIT NONE
        INTEGER :: LEVEL

        LEVEL = VERBOSITY_LEVEL

    END FUNCTION

END MODULE
!> @}

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

!> @addtogroup auxoptions
!> @{

!> @brief Enable/Disable to use of OpenMP accelerated routines.
!>
!> The MEPACK_OPTIONS_OPENMP module provide a global option to enable
!> or disable the OpenMP accelerated solvers. This is required if
!> MEPACK is executed from MATLAB.
MODULE MEPACK_OPTIONS_OPENMP
    IMPLICIT NONE

    INTEGER, SAVE, PRIVATE :: USE_OPENMP = 1

    CONTAINS

        !> @brief Return if OpenMP is enabled.
        !> @return A non zero value if OpenMP is enabled.
        !>
        !> The OPENMP_ENABLED function returns a non zero value
        !> if the OpenMP accelerated solvers should be used. This only
        !> influences the automatic solver selection.
        FUNCTION OPENMP_ENABLED() RESULT(OPENMP)
            IMPLICIT NONE
            INTEGER :: OPENMP

            OPENMP = USE_OPENMP
        END FUNCTION OPENMP_ENABLED

        !> @brief Enable or disable the use of OpenMP accelerated routines.
        !> @param[in]  OPENMP     Integer identifying if OpenMP is enabled or not.
        !>
        !> The OPENMP_ENABLE function enables or disables the use of OpenMP in the automatic
        !> solver selection. If the OpenMP solver gets selected and OpenMP is not allowed,
        !> the standard level-3 solver is used. This can be beneficial in MATLAB if the OpenMP
        !> routines crash.
        !>
        SUBROUTINE OPENMP_ENABLE(OPENMP)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: OPENMP
            USE_OPENMP = OPENMP
            RETURN
        END SUBROUTINE OPENMP_ENABLE


END MODULE MEPACK_OPTIONS_OPENMP
!> @}

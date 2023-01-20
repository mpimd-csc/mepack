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

!> @addtogroup sortev
!> @{

!> @brief Options to sort eigenvalues in the frontend.
!>
!> The MEPACK_OPTIONS_SORTEV module provide a global option to enable
!> or disable the sort of the eigenvalues before a solver
!> is called from the frontend routines.
MODULE MEPACK_OPTIONS_SORTEV
    IMPLICIT NONE

    INTEGER, SAVE, PRIVATE :: SORTEV_ENABLED_VAR = 0

    CONTAINS

        !> @brief Return if the sorting of the eigenvalues is enabled.
        !> @return A non zero value if the eigenvalues should be sorted.
        !>
        !> The SORTEV_ENABLED function returns a non zero value if the eigenvalues
        !> of the coefficient matrices should be sorted before the equations are solved.
        FUNCTION SORTEV_ENABLED() RESULT(SORTEV)
            IMPLICIT NONE
            INTEGER :: SORTEV

            SORTEV = SORTEV_ENABLED_VAR
        END FUNCTION SORTEV_ENABLED

        !> @brief Enable or disable the sorting of the eigenvalues.
        !> @param[in]  SORT         Integer identifying if the eigenvalues should be sorted or not. #
        !>
        !> The SORTEV_ENABLE subroutine enables or disable the sorting of the eigenvalues
        !> before a triangular matrix equation solver is called in a frontend routine.
        !>
        SUBROUTINE SORTEV_ENABLE(SORT)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: SORT
            SORTEV_ENABLED_VAR = SORT
            RETURN
        END SUBROUTINE SORTEV_ENABLE


END MODULE MEPACK_OPTIONS_SORTEV
!> @}

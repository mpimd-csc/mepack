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

!> @brief Get the MEPACK version number
!> @param[out] MAJOR        MEPACK major version number
!> @param[out] MINOR        MEPACK minor version number
!> @param[out] PATCHLEVEL   MEPACK patch level
!>
!> The MEPACK_GET_VERSION subroutine returns the MEPACK version number
!> split into its components.
!>
!> @see mepack_version
!> @ingroup aux
SUBROUTINE MEPACK_GET_VERSION(MAJOR, MINOR, PATCHLEVEL)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: MAJOR, MINOR, PATCHLEVEL

    MAJOR = MEPACK_MAJOR
    MINOR = MEPACK_MINOR
    PATCHLEVEL = MEPACK_PATCH_LEVEL

END SUBROUTINE MEPACK_GET_VERSION

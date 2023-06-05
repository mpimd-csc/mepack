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

!>  \brief Error handler callback subroutine 
!> 
!   Definition:
!   ===========
!
!       SUBROUTINE XERROR_HANDLER ( NAME, INFO ) 
!            CHARACTER NAME(*) 
!            INTEGER INFO 
!
!  
!>  \par Purpose:
!   =============
!>
!>  \verbatim 
!>  The XERROR_HANDLER subroutine is the callback function used by the sanity checks
!>  of the computational subroutines. It either displays an error message of 
!>  if an alternative error handler is set before via XERROR_SET_HANDLER_C or 
!>  XERROR_SET_HANDLER_F it delegates the error to this function. 
!>  \endverbatim 
!> 
!>  \see xerror_set_handler_c
!>  \see xerror_set_handler_f
!  
!>  \par Arguments:
!   ===============
!>  
!>  \param[in] NAME
!>  \verbatim
!>           NAME is CHARACTER(*) 
!>           Specifies the name of the function where the error appeared. 
!>  \endverbatim
!>  
!>  \param[in] INFO 
!>  \verbatim
!>           INFO is INTEGER 
!>           Specifies the error which appeared in the subroutine NAME. 
!>  \endverbatim
!>  
!    Authors:
!    ========
!>  
!>  \author Martin Koehler, MPI Magdeburg 
!>  
!>  \date June 2023
!>  \ingroup auxerror
SUBROUTINE XERROR_HANDLER(NAME, INFO)
    RETURN 
END SUBROUTINE

!>  \brief Set a C routine as error handling callback.
!> 
!   Definition:
!   ===========
!
!       SUBROUTINE XERROR_SET_HANDLER_C ( ROUTINE ) 
!            EXTERNAL ROUTINE 
!
!  
!>  \par Purpose:
!   =============
!>
!>  \verbatim 
!>  The XERROR_SET_HANDLER_C subroutine sets the error callback in XERROR_HANDLER to 
!>  a C function with the following signature: 
!> 
!>      void error_handler(const char *name, int info); 
!> 
!>  where name is a 0 terminated C compatible string. In order to avoid memory leaks
!>  the function must return to its call or exit the whole program.  
!>  \endverbatim 
!> 
!>  \see xerror_handler 
!>  \see xerror_set_handler_f
!  
!>  \par Arguments:
!   ===============
!>  
!>  \param[in] ROUTINE
!>  \verbatim
!>           ROUTINE is an EXTERNAL SUBROUTINE 
!>           Specifies a C subroutine which is executed as alternative error handler. 
!>  \endverbatim
!>  
!>   Authors:
!    ========
!>  
!>  \author Martin Koehler, MPI Magdeburg 
!>  
!>  \date June 2023
!>  \ingroup auxerror
SUBROUTINE XERROR_SET_HANDLER_C(ROUTINE)
    RETURN 
END SUBROUTINE

!>  \brief Set a Fortran routine as error handling callback.
!> 
!   Definition:
!   ===========
!
!       SUBROUTINE XERROR_SET_HANDLER_F ( ROUTINE ) 
!            EXTERNAL ROUTINE 
!
!  
!>  \par Purpose:
!   =============
!>
!>  \verbatim 
!>  The XERROR_SET_HANDLER_F subroutine sets the error callback in XERROR_HANDLER to 
!>  a Fortran function with the following signature: 
!> 
!>      SUBROUTINE ERROR_HANDLER(NAME, INFO) 
!>          CHARACTER NAME(*) 
!>          INTEGER INFO 
!> 
!>  where NAME is a Fortran-like string. In order to avoid memory leaks
!>  the function must return to its call or exit the whole program.  
!>  \endverbatim 
!> 
!>  \see xerror_handler
!>  \see xerror_set_handler_c
!  
!>  \par Arguments:
!   ===============
!>  
!>  \param[in] ROUTINE
!>  \verbatim
!>           ROUTINE is an EXTERNAL SUBROUTINE 
!>           Specifies a Fortran subroutine which is executed as alternative error handler. 
!>  \endverbatim
!>  
!>   Authors:
!    ========
!>  
!>  \author Martin Koehler, MPI Magdeburg 
!>  
!>  \date June 2023
!>  \ingroup auxerror
SUBROUTINE XERROR_SET_HANDLER_F(ROUTINE)
    RETURN 
END SUBROUTINE


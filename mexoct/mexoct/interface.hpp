/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) Martin Koehler, 2019
 */

/**
 * @file interface.hpp
 */

#pragma once
#ifndef MEXOCT_EXTERN_C
#  ifdef __cplusplus
#    define MEXOCT_EXTERN_C extern "C"
#  else
#    define MEXOCT_EXTERN_C extern
#  endif
#endif
#if defined(_WIN32) || defined(_MSC_VER)
#define MEXOCT_DLL_EXPORT_SYM __declspec(dllexport)
#define MEXOCT_DLL_IMPORT_SYM __declspec(dllimport)
#elif __GNUC__ >= 4
#define MEXOCT_DLL_EXPORT_SYM __attribute__((visibility("default")))
#define MEXOCT_DLL_IMPORT_SYM __attribute__((visibility("default")))
#else
#define MEXOCT_DLL_EXPORT_SYM
#define MEXOCT_DLL_IMPORT_SYM
#endif

#if __cplusplus <= 199711L
  #error This library needs at least a C++11 compliant compiler
#endif

#include <cstdint>
#ifdef MEXOCT_OCTAVE
    /* Octave Interface  */
    #include <octave/oct.h>

    #if ( OCTAVE_MAJOR_VERSION == 4 && OCTAVE_MINOR_VERSION == 0 ) || OCTAVE_MAJOR_VERSION < 4
        #error At least OCTAVE 4.2.x is required.
    #endif
    #define mexoct_int octave_idx_type
    #ifdef F77_INT
    #define mexoct_fortran_int F77_INT
    #else
    #define mexoct_fortran_int int32_t
    #endif
#else
    /* Matlab Interface */
    #include "mex.h"
    #define mexoct_int mwIndex
    #define mexoct_fortran_int mwSignedIndex
#endif

    #include "mexoct/utils.hpp"
    #include "mexoct/errors.hpp"
    #include "mexoct/paramtypes.hpp"
    #include "mexoct/parameter.hpp"

#ifdef MEXOCT_OCTAVE
    #include "mexoct/octave/octave_types.hpp"
#else
    #include "mexoct/matlab/matlab_types.hpp"
#endif

    #include "mexoct/scalar_type.hpp"
    #include "mexoct/array_type.hpp"
    #include "mexoct/traits.hpp"
    #include "mexoct/string.hpp"
    #include "mexoct/print.hpp"
    #include "mexoct/memory.hpp"
    #include "mexoct/function_signature.hpp"
    #include "mexoct/argparser.hpp"

/**
 * \defgroup macros Macro Defintions
 *
 * The macro defintions are required to build the interface according to Octave's or
 * MATLAB's interface definition. By preprocessesor defines from outside, the macros
 * distinguish between GNU Octave and MATLAB.
 *
 * @{
 */

/**
 * \def MEXOCT_INIT
 * \brief Initialize the mexoct interface.
 *
 * The #MEXOCT_INIT() macro initializes mexoct and needs to be called at first in the \ref MEXOCT_ENTRY function.
 * Depending on the interface build, if MATLAB or OCTAVE, it sets up the necessary variables.
 *
 */

/**
 * \def MEXOCT_ENTRY
 * \brief Defines the entry point of a Mex or an Oct file.
 * \param name      Name of the MATLAB/OCTAVE function
 * \param desc      The help text, including the description, of the function.
 *
 *
 * The #MEXOCT_ENTRY(name,desc) macro defines the entry point of a mex-file or an oct-file. It is expanded to
 * its required form using a set of macros. The name argument, which defines the name of the function, should be
 * equal to the name of the file (without its extension). The desc argumment is used to document the file. In case
 * of Octave this is passed to Octave's internal documentation feature. If the interface is compiled for MATLAB, the
 * description text is included as a static string in the mex file and it extraced at compile time to generate
 * an m-file including the corresponding documentation block. A linebreak in the ouput needs to be written as \c \n
 * and a line break in the source file needs a trailing \ in each line.
 *
 * \b Example:
\code{.cpp}
MEXOCT_ENTRY(lu_nopiv, "\
lu_nopiv  LU decomposition without pivoting\n\
\n\
    [L,U] = LU(A) stores an upper triangular matrix in U and a\n\
    \"psychologically lower triangular matrix\" (i.e. a product of lower\n\
    triangular and permutation matrices) in L, so that A = L*U. A can be\n\
    rectangular.\n")
{
    MEXOCT_INIT();

    ....

    MEXOCT_RETURN;
}
\endcode
 */

/**
 * \def MEXOCT_PARSE
 * \brief Execute the parser and execute the matched function.
 * \param parser  Parser object to be used.
 *
 * The #MEXOCT_PARSE(parser) macro parses the functions arguments given from MATLAB or OCTAVE and
 * feeds them to the addressed parser. The parser object then executes the desired piece of code.
 * If an error derived from \ref MexOctException is thrown, it is catched and \ref mexoct_error is
 * called on it. The parser is created with the help of the \ref mexoct::makeArgParser function.
 *
 * @see MEXOCT_PARSE_NOEXCEPT
 */


/**
 * \def MEXOCT_PARSE_NOEXCEPT
 * \brief Execute the parser and execute the matched function.
 * \param parser  Parser object to be used.
 *
 * The #MEXOCT_PARSE_NOEXCEPT(parser) macro parses the functions arguments given from MATLAB or OCTAVE and
 * feeds them to the addressed parser. The parser object then executes the desired piece of code.
 *
 * @attention In contrast to the \ref MEXOCT_PARSE macro this one does not enclose the call to the
 * parser with an exception handling block.
 *
 * @see MEXOCT_PARSE
 */

/**
 * \def MEXOCT_RETURN
 * \brief Returns from a mex or oct function
 *
 * The #MEXOCT_RETURN macro defines the end of a mex-function or an oct-function. It is necessary to
 * call this macro at the end of each mex-function or oct-function. If called before, the execution
 * is stopped and the control flow returns to MATLAB or Octave.
 */
#ifdef MEXOCT_OCTAVE
    #define MEXOCT_INIT()           octave_value_list * __mexoct_retval;
    #define MEXOCT_ENTRY(name,desc) char const * octHelpTxt = desc; \
                                    MEXOCT_EXTERN_C char const * octHelp(void) { return octHelpTxt; } \
                                    DEFUN_DLD ( name , __mexoct_args, __mexoct_nargout, desc )
    #define MEXOCT_PARSE(parser)    try { \
                                        (parser).parse(__mexoct_nargout, &__mexoct_retval, __mexoct_args.length(),  __mexoct_args); \
                                    } catch ( MexOctException & e ) { \
                                        mexoct_error(e); \
                                    }
    #define MEXOCT_PARSE_NOEXCEPT(parser) do { (parser).parse(__mexoct_nargout, &__mexoct_retval, __mexoct_args.length(),  __mexoct_args);} while(0)

    #define MEXOCT_RETURN return *__mexoct_retval;

#else
    #define MEXOCT_INIT()
    #define MEXOCT_ENTRY(name, desc) char const * mexHelp = desc; \
                                     MEXOCT_EXTERN_C MEXOCT_DLL_EXPORT_SYM char const * mexoct_gethelptext(); \
                                     char const * mexoct_gethelptext() { return mexHelp;  }\
                                     MEXOCT_EXTERN_C MEXOCT_DLL_EXPORT_SYM void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
    #define MEXOCT_PARSE(parser)     try {\
                                        (parser).parse(nlhs, plhs, nrhs, prhs); \
                                     } catch ( MexOctException &e ) {\
                                         mexoct_error( e ); \
                                     }
    #define MEXOCT_PARSE_NOEXCEPT(parser)  do { (parser).parse(nlhs, plhs, nrhs, prhs); } while(0);

    #define MEXOCT_RETURN do { return; } while(0)

#endif

/**
 * @}
 */

/*  #undef MEXOCT_EXTERN_C
#undef MEXOCT_DLL_IMPORT_SYM
#undef MEXOCT_DLL_EXPORT_SYM */

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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
 * Copyright (C) Martin Koehler, 2017-2023
 */


#include "FCMangle.h"
#include "mepack.h"
#include "mepack_internal.h"

/**
 * @addtogroup auxoptions
 * @{
 */

/**
 * @brief C Interface for MEPACK_OPTIONS_VERBOSE::VERBOSE_INIT
 *
 * The mepack_verbose_init function wraps the Fortran subroutine MEPACK_OPTIONS_VERBOSE::VERBOSE_INIT such that
 * it can be called from C easily.
 *
 * @see mepack_options_verbose::verbose_init
 */
void mepack_verbose_init() {

   FC_MODULE_(mepack_options_verbose,verbose_init,MEPACK_OPTIONS_VERBOSE,VERBOSE_INIT)();
}

/**
 * @brief C Interface for MEPACK_OPTIONS_VERBOSE::VERBOSE_SET
 * @param[in] level     New verbosity level
 *
 * The mepack_verbose_set function is a wrapper around MEPACK_OPTIONS_VERBOSE::VERBOSE_SET. It set
 * the verbosity level of MEPACK.
 *
 * @see mepack_options_verbose::verbose_set
 */
void mepack_verbose_set(int level) {
    Int _l = level;
    FC_MODULE_(mepack_options_verbose,verbose_set,MEPACK_OPTIONS_VERBOSE,VERBOSE_SET)(&_l);
}

/**
 * @brief C Interface for MEPACK_OPTIONS_VERBOSE::VERBOSE_GET
 * @return The current verbosity level of MEPACK
 *
 * The mepack_verbose_get function is a wrapper around MEPACK_OPTIONS_VERBOSE::VERBOSE_GET.
 * It returns the current verbosity level of MEPACK.
 *
 * @see mepack_options_verbose::verbose_get
 */
int mepack_verbose_get() {
    return FC_MODULE_(mepack_options_verbose,verbose_get,MEPACK_OPTIONS_VERBOSE,VERBOSE_GET)();
}
/**
 * @}
 */


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
 * @brief Return if the OpenMP accelerated subroutines are enabled.
 * @return A non-zero if the OpenMP accelerated subroutines are enabled.
 *
 * The mepack_openmp_enabled function returns a true value if the
 * OpenMP accelerated subroutines are used.
 * This is beneficial in MATLAB if the OpenMP accelerated solvers cause MATLAB to crash.
 */
int mepack_openmp_enabled(void)
{
    int ret;
    ret = (int)  FC_MODULE_(mepack_options_openmp,openmp_enabled,MEPACK_OPTIONS_OPENMP,OPENMP_ENABLED)();
    return ret;
}

/**
 * @brief Enable or disable the OpenMP accelerated solvers.
 * @param[in] openmp    Boolean integer value to enable or disable the OpenMP accelerated solvers.
 *
 * The mepack_openmp_enable function enables or disables the OpenMP accelerated subroutines.
 *
 */
void mepack_openmp_enable(int openmp)
{
    Int _openmp = openmp;
    FC_MODULE_(mepack_options_openmp,openmp_enable,MEPACK_OPTIONS_OPENMP,OPENMP_ENABLE)(&_openmp);
    return;
}


/**
 * @}
 */


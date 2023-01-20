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
 * @addtogroup sortev
 * @{
 */

/**
 * @brief Return if the eigenvalue sorting is enabled or not.
 * @return A non-zero if the eigenvalue sorting is enabled.
 *
 * The mepack_sortev_enabled function returns a true value if the eigenvalue
 * sorting is enabled. That means, frontend routines sort the eigenvalues
 * before calling a triangular solver in order to match the block structure.
 */
int mepack_sortev_enabled()
{
    int ret;
    ret = (int)  FC_MODULE_(mepack_options_sortev,sortev_enabled,MEPACK_OPTIONS_SORTEV,SORTEV_ENABLED)();
    return ret;
}

/**
 * @brief Enable or disable the eigenvalues sorting in the frontend routines.
 * @param[in] sort      Boolean integer value to enable or disable the eigenvalue sorting.
 *
 * The mepack_sortev_enable function enables or disables the sorting of the eigenvalues
 * before a triangular solver is called by one of the frontend routines.
 *
 */
void mepack_sortev_enable(int sort)
{
    Int _sort = sort;
    FC_MODULE_(mepack_options_sortev,sortev_enable,MEPACK_OPTIONS_SORTEV,SORTEV_ENABLE)(&_sort);
    return;
}


/**
 * @}
 */

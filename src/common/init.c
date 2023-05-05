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
 * @brief Initialize MEPACK from C
 * @ingroup auxinit
 *
 * The mepack_init function is a wrapper around the \ref mepack_initialize routine
 * to be called from C easily. If accidentally both, the C and the Fortran initialization routine,
 * are called only the one which is called first is executed. If the function is called twice, nothing
 * happens.
 *
 * @see mepack_initialize
 */
void mepack_init(void)
{

    FC_GLOBAL_(mepack_initialize,MEPACK_INITIALIZE) ();
    return;
}


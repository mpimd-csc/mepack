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
 * Copyright (C) Martin Koehler, 2017-2022
 */


#include "FCMangle.h"
#include "mepack.h"
#include "mepack_internal.h"


/**
 * @brief Get the version information of MEPACK
 * @param[out]  MAJOR       Major version number of MEPACK
 * @param[out]  MINOR       Minor version number of MEPACK
 * @param[out]  PATCHLEVEL  Patchlevel of MEPACK
 * @ingroup aux
 *
 * The mepack_version function returns the version of MEPACK. The version number is returned in
 * its three components MAJOR, MINOR and PATCHLEVEL.
 */
void mepack_version(int *MAJOR, int *MINOR, int *PATCHLEVEL)
{
    *MAJOR = MEPACK_MAJOR;
    *MINOR = MEPACK_MINOR;
    *PATCHLEVEL = MEPACK_PATCH_LEVEL;
    return;
}


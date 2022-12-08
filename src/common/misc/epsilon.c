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

#ifdef INTEGER8
#include <stdint.h>
#define Int int64_t
#endif
#ifndef Int
#define Int int
#endif
#include <stdlib.h>

#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif

double FC_GLOBAL(dlamch,DLAMCH)(char *str, fortran_charlen_t len);
float  FC_GLOBAL(slamch,SLAMCH)(char *str, fortran_charlen_t len);



double mepack_double_epsilon()
{
    return FC_GLOBAL(dlamch,DLAMCH)("E", 1);
}


float mepack_single_epsilon()
{
    return FC_GLOBAL(slamch,SLAMCH)("E",1);
}




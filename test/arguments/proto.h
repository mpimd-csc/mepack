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

#ifdef __cplusplus
extern "C" {
#endif
#include <ctype.h>
#include <math.h>
#include "FCMangle.h"
#ifdef INTEGER8
#include <stdint.h>
#define Int int64_t
#endif
#ifndef Int
#define Int int
#endif

#include "mepack.h"
#include "mepack_internal.h"

#define MAX(A,B) (((A)>(B))?(A):(B))

    /*-----------------------------------------------------------------------------
     *  Error Handler
     *-----------------------------------------------------------------------------*/
    void FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C) ( void (*eh)(const char*,Int));
    void FC_GLOBAL_(xerror_set_handler_f,XERROR_SET_HANDLER_F) ( void (*eh)(const char*,Int*, Int));
    /*-----------------------------------------------------------------------------
     *  MEPACK Double
     *-----------------------------------------------------------------------------*/


    extern void FC_GLOBAL_(dla_sort_gev,DLA_SORT_GEV)( Int * M,  double * A, Int * lDA, double *C, Int *ldc, double * Q, Int *ldq, double *Z, Int *ldz, Int *NB, double * work, Int *ldwork, Int *info );
    extern void FC_GLOBAL_(dla_sort_ev,DLA_SORT_EV)( Int * M,  double * A, Int * lDA, double * Q, Int *ldq, Int *NB, double * work, Int *ldwork, Int *info );
    extern void FC_GLOBAL_(sla_sort_ev,SLA_SORT_EV)( Int * M,  float * A, Int * lDA, float * Q, Int *ldq, Int *NB, float * work, Int *ldwork, Int *info );
    extern void FC_GLOBAL_(sla_sort_gev,SLA_SORT_GEV)( Int * M,  float * A, Int * lDA, float *C, Int *ldc, float * Q, Int *ldq, float *Z, Int *ldz, Int *NB, float * work, Int *ldwork, Int *info );

    #ifdef __cplusplus
};
#endif

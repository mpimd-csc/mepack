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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @addtogroup auxmem
 * @{
 */

/**
 * \brief Query the size of the workspace for triangular solvers.
 *
 * \par Purpose:
 *
 * The mepack_memory function queries the size of the required routine. It wraps around the Fortran
 * style workspace query and adds a typecast to avoid integer overflows. This can happen if MEPACK
 * is compiled with the support for 64-bit integers in Fortran. If solver only relies on one dimension
 * argument, e.g. a Lyapunov or a Stein equation is planed to be solved, only the \b M parameter is used.
 * The value of \b N is neglected in this case. The frontend routines, i.e. the ones dealing with the
 * transformation matrices as well, have to use \ref mepack_memory_frontend to query their memory
 * requirements.
 *
 * \see mepack_memory_frontend
 *
 * \remark
 *    This function should be used from C. If the workspace query is performed from a Fortran code,
 *    the workspace is queried by passing `INFO == -1` to the corresponding solver.
 *
 * \param[in] func
 *    The name of the solver for which the size of the workspace should be queried. The argument is not
 *    case sensitive.
 *
 * \param[in] M
 *    The number of rows of the solution, the right hand side, and the order of the matrices from
 *    the left side. If the equation only depends on one dimension argument, this is given here.
 *
 * \param[in] N
 *    The number of columns of the solution, the right hand side, and the order of the matrices from
 *    the right side. If the equation only depends on one dimension argument, this parameter is not used.
 *
 * \return If the return value is greater or equal to zero, it is the workspace required by the desired
 *    routine counted in floating point number of the desired precision. That means the return value
 *    either needs to be multiplied by `sizeof(float)` or `sizeof(double)`. A negative return value
 *    indicates an error. This error code is explained in the description of the corresponding routine
 *    for which the memory query was performed.
 */
ssize_t mepack_memory(const char *func, int M, int N)
{
    Int info = -1;

    /*
     * TRLYAP (DOUBLE)
     */
    if ((strcasecmp("DLA_TRLYAP_L3", func) == 0)) {
        FC_GLOBAL_(dla_trlyap_l3,DLA_TRLYAP_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TRLYAP_L3_2S", func) == 0)) {
        FC_GLOBAL_(dla_trlyap_l3_2s,DLA_TRLYAP_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TRLYAP_DAG", func) == 0)){
        FC_GLOBAL_(dla_trlyap_dag,DLA_TRLYAP_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TRLYAP_L2_OPT", func) == 0)) {
        FC_GLOBAL_(dla_trlyap_l2_opt,DLA_TRLYAP_L2_OPT)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TRLYAP_L2", func) == 0) ) {
        FC_GLOBAL_(dla_trlyap_l2,DLA_TRLYAP_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TRLYAP_RECURSIVE", func) == 0))  {
        FC_GLOBAL_(dla_trlyap_recursive,DLA_TRLYAP_RECURSIVE)("N", &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    /*
     * TRLYAP (SINGLE)
     */
    if ((strcasecmp("SLA_TRLYAP_L3", func) == 0)) {
        FC_GLOBAL_(sla_trlyap_l3,SLA_TRLYAP_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TRLYAP_L3_2S", func) == 0)) {
        FC_GLOBAL_(sla_trlyap_l3_2s,SLA_TRLYAP_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TRLYAP_DAG", func) == 0)) {
        FC_GLOBAL_(sla_trlyap_dag,SLA_TRLYAP_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TRLYAP_L2_OPT", func) == 0))  {
        FC_GLOBAL_(sla_trlyap_l2_opt,SLA_TRLYAP_L2_OPT)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TRLYAP_L2", func) == 0) ) {
        FC_GLOBAL_(sla_trlyap_l2,SLA_TRLYAP_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TRLYAP_RECURSIVE", func) == 0)){
        FC_GLOBAL_(sla_trlyap_recursive,SLA_TRLYAP_RECURSIVE)("N", &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }



    /*
     * TRSTEIN (DOUBLE)
     */
    if ((strcasecmp("DLA_TRSTEIN_L3", func) == 0)  ) {
        FC_GLOBAL_(dla_trstein_l3,DLA_TRSTEIN_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TRSTEIN_L3_2S", func) == 0) ) {
        FC_GLOBAL_(dla_trstein_l3_2s,DLA_TRSTEIN_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TRSTEIN_DAG", func) == 0)  ) {
        FC_GLOBAL_(dla_trstein_dag,DLA_TRSTEIN_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TRSTEIN_L2", func) == 0)  ) {
        FC_GLOBAL_(dla_trstein_l2,DLA_TRSTEIN_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TRSTEIN_RECURSIVE", func) == 0)  ) {
        FC_GLOBAL_(dla_trstein_recursive,DLA_TRSTEIN_RECURSIVE)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    /*
     * TRSTEIN (SINGLE)
     */
    if ((strcasecmp("SLA_TRSTEIN_L3", func) == 0)  ) {
        FC_GLOBAL_(sla_trstein_l3,SLA_TRSTEIN_L3)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TRSTEIN_L3_2S", func) == 0) ) {
        FC_GLOBAL_(sla_trstein_l3_2s,SLA_TRSTEIN_L3_2S)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TRSTEIN_DAG", func) == 0)  ) {
        FC_GLOBAL_(sla_trstein_dag,SLA_TRSTEIN_DAG)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TRSTEIN_L2", func) == 0) ) {
        FC_GLOBAL_(sla_trstein_l2,SLA_TRSTEIN_L2)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TRSTEIN_RECURSIVE", func) == 0)  ) {
        FC_GLOBAL_(sla_trstein_recursive,SLA_TRSTEIN_RECURSIVE)("N",  &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*X*/, &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    /*
     * TRSYLV (DOUBLE)
     */
    if ( strcasecmp("DLA_TRSYLV_L3", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l3,DLA_TRSYLV_L3)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L3_UNOPT", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l3_unopt,DLA_TRSYLV_L3_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_DAG", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_dag,DLA_TRSYLV_DAG)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L3_2S", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l3_2s,DLA_TRSYLV_L3_2S)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2,DLA_TRSYLV_L2)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_UNOPT", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2_unopt,DLA_TRSYLV_L2_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_REORDER", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2_reorder,DLA_TRSYLV_L2_REORDER)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2_local_copy,DLA_TRSYLV_L2_LOCAL_COPY)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_32", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2_local_copy_32,DLA_TRSYLV_L2_LOCAL_COPY_32)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_64", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2_local_copy_64,DLA_TRSYLV_L2_LOCAL_COPY_64)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_96", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2_local_copy_96,DLA_TRSYLV_L2_LOCAL_COPY_96)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_L2_LOCAL_COPY_128", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_l2_local_copy_128,DLA_TRSYLV_L2_LOCAL_COPY_128)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV_RECURSIVE", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv_recursive,DLA_TRSYLV_RECURSIVE)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TRSYLV (Single)
     */
    if ( strcasecmp("SLA_TRSYLV_L3", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l3,SLA_TRSYLV_L3)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L3_UNOPT", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l3_unopt,SLA_TRSYLV_L3_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_DAG", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_dag,SLA_TRSYLV_DAG)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L3_2S", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l3_2s,SLA_TRSYLV_L3_2S)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2,SLA_TRSYLV_L2)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_UNOPT", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2_unopt,SLA_TRSYLV_L2_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_REORDER", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2_reorder,SLA_TRSYLV_L2_REORDER)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2_local_copy,SLA_TRSYLV_L2_LOCAL_COPY)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_32", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2_local_copy_32,SLA_TRSYLV_L2_LOCAL_COPY_32)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_64", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2_local_copy_64,SLA_TRSYLV_L2_LOCAL_COPY_64)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_96", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2_local_copy_96,SLA_TRSYLV_L2_LOCAL_COPY_96)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_L2_LOCAL_COPY_128", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_l2_local_copy_128,SLA_TRSYLV_L2_LOCAL_COPY_128)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV_RECURSIVE", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv_recursive,SLA_TRSYLV_RECURSIVE)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TRSYLV2 (DOUBLE)
     */
    if ( strcasecmp("DLA_TRSYLV2_L3", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l3,DLA_TRSYLV2_L3)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L3_UNOPT", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l3_unopt,DLA_TRSYLV2_L3_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_DAG", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_dag,DLA_TRSYLV2_DAG)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L3_2S", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l3_2s,DLA_TRSYLV2_L3_2S)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2,DLA_TRSYLV2_L2)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_UNOPT", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2_unopt,DLA_TRSYLV2_L2_UNOPT)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_REORDER", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2_reorder,DLA_TRSYLV2_L2_REORDER)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2_local_copy,DLA_TRSYLV2_L2_LOCAL_COPY)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_32", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2_local_copy_32,DLA_TRSYLV2_L2_LOCAL_COPY_32)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_64", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2_local_copy_64,DLA_TRSYLV2_L2_LOCAL_COPY_64)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_96", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2_local_copy_96,DLA_TRSYLV2_L2_LOCAL_COPY_96)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_L2_LOCAL_COPY_128", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_l2_local_copy_128,DLA_TRSYLV2_L2_LOCAL_COPY_128)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TRSYLV2_RECURSIVE", func) == 0 ) {
         FC_GLOBAL_(dla_trsylv2_recursive,DLA_TRSYLV2_RECURSIVE)("N", "N", &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TRSYLV2 (Single)
     */
    if ( strcasecmp("SLA_TRSYLV2_L3", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l3,SLA_TRSYLV2_L3)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L3_UNOPT", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l3_unopt,SLA_TRSYLV2_L3_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_DAG", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_dag,SLA_TRSYLV2_DAG)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L3_2S", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l3_2s,SLA_TRSYLV2_L3_2S)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2,SLA_TRSYLV2_L2)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_UNOPT", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2_unopt,SLA_TRSYLV2_L2_UNOPT)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_REORDER", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2_reorder,SLA_TRSYLV2_L2_REORDER)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2_local_copy,SLA_TRSYLV2_L2_LOCAL_COPY)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_32", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2_local_copy_32,SLA_TRSYLV2_L2_LOCAL_COPY_32)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_64", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2_local_copy_64,SLA_TRSYLV2_L2_LOCAL_COPY_64)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_96", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2_local_copy_96,SLA_TRSYLV2_L2_LOCAL_COPY_96)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_L2_LOCAL_COPY_128", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_l2_local_copy_128,SLA_TRSYLV2_L2_LOCAL_COPY_128)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TRSYLV2_RECURSIVE", func) == 0 ) {
         FC_GLOBAL_(sla_trsylv2_recursive,SLA_TRSYLV2_RECURSIVE)("N", "N", &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }


    /*
     * TGLYAP (DOUBLE)
     */
    if ( strcasecmp("DLA_TGLYAP_L3", func) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l3,DLA_TGLYAP_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGLYAP_L3_2S", func) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l3_2s,DLA_TGLYAP_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGLYAP_DAG", func) == 0 ) {
        FC_GLOBAL_(dla_tglyap_dag,DLA_TGLYAP_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGLYAP_L2", func) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l2,DLA_TGLYAP_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGLYAP_L2_UNOPT", func) == 0 ) {
        FC_GLOBAL_(dla_tglyap_l2_unopt,DLA_TGLYAP_L2_UNOPT)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGLYAP_RECURSIVE", func) == 0 ) {
        FC_GLOBAL_(dla_tglyap_recursive,DLA_TGLYAP_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }




    /*
     * TGLYAP (SINGLE)
     */
    if ( strcasecmp("SLA_TGLYAP_L3", func) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l3,SLA_TGLYAP_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGLYAP_L3_2S", func) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l3_2s,SLA_TGLYAP_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGLYAP_DAG", func) == 0 ) {
        FC_GLOBAL_(sla_tglyap_dag,SLA_TGLYAP_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGLYAP_L2", func) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l2,SLA_TGLYAP_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGLYAP_L2_UNOPT", func) == 0 ) {
        FC_GLOBAL_(sla_tglyap_l2_unopt,SLA_TGLYAP_L2_UNOPT)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGLYAP_RECURSIVE", func) == 0 ) {
        FC_GLOBAL_(sla_tglyap_recursive,SLA_TGLYAP_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }

    /*
     * TGSTEIN (DOUBLE)
     */
    if ( strcasecmp("DLA_TGSTEIN_L3", func) == 0 ) {
        FC_GLOBAL_(dla_tgstein_l3,DLA_TGSTEIN_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGSTEIN_L3_2S", func) == 0 ) {
        FC_GLOBAL_(dla_tgstein_l3_2s,DLA_TGSTEIN_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGSTEIN_DAG", func) == 0 ) {
        FC_GLOBAL_(dla_tgstein_dag,DLA_TGSTEIN_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGSTEIN_L2", func) == 0 ) {
        FC_GLOBAL_(dla_tgstein_l2,DLA_TGSTEIN_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGSTEIN_RECURSIVE", func) == 0 ) {
        FC_GLOBAL_(dla_tgstein_recursive,DLA_TGSTEIN_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }



    /*
     * TGSTEIN(SINGLE)
     */
    if ( strcasecmp("SLA_TGSTEIN_L3", func) == 0 ) {
        FC_GLOBAL_(sla_tgstein_l3,SLA_TGSTEIN_L3)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGSTEIN_L3_2S", func) == 0 ) {
        FC_GLOBAL_(sla_tgstein_l3_2s,SLA_TGSTEIN_L3_2S)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGSTEIN_DAG", func) == 0 ) {
        FC_GLOBAL_(sla_tgstein_dag,SLA_TGSTEIN_DAG)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGSTEIN_L2", func) == 0 ) {
        FC_GLOBAL_(sla_tgstein_l2,SLA_TGSTEIN_L2)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGSTEIN_RECURSIVE", func) == 0 ) {
        FC_GLOBAL_(sla_tgstein_recursive,SLA_TGSTEIN_RECURSIVE)("N", &(Int) {M}, NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){M}, NULL /*X*/ , &(Int){M},
            NULL /*SCALE*/, NULL /*WORK*/, &info, 1);
        return (ssize_t) info;
    }





    /*
     * TGSYLV (DOUBLE)
     */
    if ((strcasecmp("DLA_TGSYLV_L3", func ) == 0 ) ) {
        FC_GLOBAL_(dla_tgsylv_l3,DLA_TGSYLV_L3)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L3_ROWWISE", func ) == 0 )) {
        FC_GLOBAL_(dla_tgsylv_l3_rowwise,DLA_TGSYLV_L3_ROWWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L3_COLWISE", func ) == 0 )) {
        FC_GLOBAL_(dla_tgsylv_l3_colwise,DLA_TGSYLV_L3_COLWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L3_2S", func ) == 0 )) {
        FC_GLOBAL_(dla_tgsylv_l3_2s,DLA_TGSYLV_L3_2S)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_DAG", func ) == 0 ) ) {
        FC_GLOBAL_(dla_tgsylv_dag,DLA_TGSYLV_DAG)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_32", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_32,DLA_TGSYLV_L2_LOCAL_COPY_32)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_64", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_64,DLA_TGSYLV_L2_LOCAL_COPY_64)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_96", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_96,DLA_TGSYLV_L2_LOCAL_COPY_96)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY_128", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy_128,DLA_TGSYLV_L2_LOCAL_COPY_128)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_LOCAL_COPY", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_local_copy,DLA_TGSYLV_L2_LOCAL_COPY)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_REORDER", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_reorder,DLA_TGSYLV_L2_REORDER)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TGSYLV_L2", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2,DLA_TGSYLV_L2)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_ROWWISE", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_rowwise,DLA_TGSYLV_L2_ROWWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("DLA_TGSYLV_L2_COLWISE", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_colwise,DLA_TGSYLV_L2_COLWISE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TGSYLV_L2_UNOPT", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_l2_unopt,DLA_TGSYLV_L2_UNOPT)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TGSYLV_RECURSIVE", func ) == 0)) {
        FC_GLOBAL_(dla_tgsylv_recursive,DLA_TGSYLV_RECURSIVE)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("DLA_TGSYLV_GARDINER_LAUB", func ) == 0 ) || (strcasecmp("SLA_TGSYLV_GARDINER_LAUB", func) == 0)) {
        FC_GLOBAL_(dla_tgsylv_gardiner_laub,DLA_TGSYLV_GARDINER_LAUB)("N", "N", &(double){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TGSYLV (SINGLE)
     */
    if ((strcasecmp("SLA_TGSYLV_L3", func ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3", func) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3,SLA_TGSYLV_L3)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L3_ROWWISE", func ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3_ROWISE", func) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3_rowwise,SLA_TGSYLV_L3_ROWWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L3_COLWISE", func ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3_COLWISE", func) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3_colwise,SLA_TGSYLV_L3_COLWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L3_2S", func ) == 0 ) || (strcasecmp("SLA_TGSYLV_L3_2S", func) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l3_2s,SLA_TGSYLV_L3_2S)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_DAG", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_dag,SLA_TGSYLV_DAG)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_32", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_32,SLA_TGSYLV_L2_LOCAL_COPY_32)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_64", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_64,SLA_TGSYLV_L2_LOCAL_COPY_64)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_96", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_96,SLA_TGSYLV_L2_LOCAL_COPY_96)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY_128", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy_128,SLA_TGSYLV_L2_LOCAL_COPY_128)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_LOCAL_COPY", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_local_copy,SLA_TGSYLV_L2_LOCAL_COPY)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_REORDER", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_reorder,SLA_TGSYLV_L2_REORDER)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TGSYLV_L2", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2,SLA_TGSYLV_L2)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_ROWWISE", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_rowwise,SLA_TGSYLV_L2_ROWWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ((strcasecmp("SLA_TGSYLV_L2_COLWISE", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_colwise,SLA_TGSYLV_L2_COLWISE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TGSYLV_L2_UNOPT", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_l2_unopt,SLA_TGSYLV_L2_UNOPT)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TGSYLV_RECURSIVE", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_recursive,SLA_TGSYLV_RECURSIVE)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    if ((strcasecmp("SLA_TGSYLV_GARDINER_LAUB", func ) == 0)) {
        FC_GLOBAL_(sla_tgsylv_gardiner_laub,SLA_TGSYLV_GARDINER_LAUB)("N", "N", &(float){-1},  &(Int){M}, &(Int){N},
                NULL  /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},
                NULL /*SCALE*/, NULL /* WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TGCSYLV(DOUBLE)
     */
    if ( strcasecmp("DLA_TGCSYLV_L2",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2,DLA_TGCSYLV_L2)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_UNOPT",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_unopt,DLA_TGCSYLV_L2_UNOPT)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_REORDER",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_reorder,DLA_TGCSYLV_L2_REORDER)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy,DLA_TGCSYLV_L2_LOCAL_COPY)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_32",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_32,DLA_TGCSYLV_L2_LOCAL_COPY_32)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_64",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_64,DLA_TGCSYLV_L2_LOCAL_COPY_64)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_96",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_96,DLA_TGCSYLV_L2_LOCAL_COPY_96)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L2_LOCAL_COPY_128",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l2_local_copy_128,DLA_TGCSYLV_L2_LOCAL_COPY_128)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L3",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l3,DLA_TGCSYLV_L3)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_L3_UNOPT",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l3_unopt,DLA_TGCSYLV_L3_UNOPT)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
     if ( strcasecmp("DLA_TGCSYLV_L3_2S",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_l3_2s,DLA_TGCSYLV_L3_2S)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DAG",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dag,DLA_TGCSYLV_DAG)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_RECURSIVE",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_recursive,DLA_TGCSYLV_RECURSIVE)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TGCSYLV_DUAL_DUAL(DOUBLE)
     */
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2,DLA_TGCSYLV_DUAL_L2)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_32,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_64,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_96,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l2_local_copy_128,DLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_L3",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l3,DLA_TGCSYLV_DUAL_L3)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
     if ( strcasecmp("DLA_TGCSYLV_DUAL_L3_2S",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_l3_2s,DLA_TGCSYLV_DUAL_L3_2S)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_DAG",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_dag,DLA_TGCSYLV_DUAL_DAG)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("DLA_TGCSYLV_DUAL_RECURSIVE",func) == 0) {
        FC_GLOBAL_(dla_tgcsylv_dual_recursive,DLA_TGCSYLV_DUAL_RECURSIVE)("N", "N", &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TGCSYLV(SINGLE)
     */
    if ( strcasecmp("SLA_TGCSYLV_L2",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2,SLA_TGCSYLV_L2)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_UNOPT",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_unopt,SLA_TGCSYLV_L2_UNOPT)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_REORDER",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_reorder,SLA_TGCSYLV_L2_REORDER)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy,SLA_TGCSYLV_L2_LOCAL_COPY)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_32",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_32,SLA_TGCSYLV_L2_LOCAL_COPY_32)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_64",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_64,SLA_TGCSYLV_L2_LOCAL_COPY_64)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_96",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_96,SLA_TGCSYLV_L2_LOCAL_COPY_96)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L2_LOCAL_COPY_128",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l2_local_copy_128,SLA_TGCSYLV_L2_LOCAL_COPY_128)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L3",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l3,SLA_TGCSYLV_L3)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_L3_UNOPT",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l3_unopt,SLA_TGCSYLV_L3_UNOPT)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
     if ( strcasecmp("SLA_TGCSYLV_L3_2S",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_l3_2s,SLA_TGCSYLV_L3_2S)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DAG",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dag,SLA_TGCSYLV_DAG)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_RECURSIVE",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_recursive,SLA_TGCSYLV_RECURSIVE)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }

    /*
     * TGCSYLV_DUAL(SINGLE)
     */
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2,SLA_TGCSYLV_DUAL_L2)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_32,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_32)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_64,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_64)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_96,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_96)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l2_local_copy_128,SLA_TGCSYLV_DUAL_L2_LOCAL_COPY_128)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L3",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l3,SLA_TGCSYLV_DUAL_L3)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_L3_2S",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_l3_2s,SLA_TGCSYLV_DUAL_L3_2S)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_DAG",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_dag,SLA_TGCSYLV_DUAL_DAG)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }
    if ( strcasecmp("SLA_TGCSYLV_DUAL_RECURSIVE",func) == 0) {
        FC_GLOBAL_(sla_tgcsylv_dual_recursive,SLA_TGCSYLV_DUAL_RECURSIVE)("N", "N", &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int) {N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int) {N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int) {M}, NULL /*SCALE*/, NULL /*WORK*/, &info, 1, 1);
        return (ssize_t) info;
    }



    return -1;
}


/**
 * @}
 */

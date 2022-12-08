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
 * \brief Query the size of the workspace for frontend routines.
 *
 * \par Purpose:
 *
 * The mepack_memory_frontend function queries the size of the required routine. It wraps around the Fortran
 * style workspace query and adds a typecast to avoid integer overflows. This can happen if MEPACK
 * is compiled with the support for 64-bit integers in Fortran. If solver only relies on one dimension
 * argument, e.g. a Lyapunov or a Stein equation is planed to be solved, only the \b M and the \b factA
 * parameters are used. The values of \b N and \b factB are neglected in this case. The triangular solvers,
 * i.e. the ones which only aware of the triangular coefficient matrices, have to use \ref mepack_memory
 * to query their memory requirements.
 *
 * \see mepack_memory
 *
 * \remark
 *    This function should be used from C. If the workspace query is performed from a Fortran code,
 *    the workspace is queried by passing `INFO == -1` to the corresponding solver.
 *
 * \param[in] func
 *    The name of the solver for which the size of the workspace should be queried. The argument is not
 *    case sensitive.
 *
 * \param[in] factA
 *    Indicates whether the first matrix pair is already factorized or not. If "N" is passed the
 *    matrix pair is not factorized. If "F" is passed it is. For Lyapunov and Stein equations only
 *    this parameter is used.
 *
 * \param[in] factB
 *    Indicates whether the second matrix pair is already factorized or not. If "N" is passed the
 *    matrix pair is not factorized. If "F" is passed it is. For Lyapunov and Stein equations this
 *    parameter is not used.
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
ssize_t mepack_memory_frontend(const char *func, const char * factA, const char *factB, int M, int N)
{
    Int info = 0;
    Int ldwork = -1;

    /*
     *  GELYAP
     */
    if ( strcasecmp("DLA_GELYAP", func) == 0 ) {
        FC_GLOBAL_(dla_gelyap,DLA_GELYAP)(factA, "N", &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*Q*/, &(Int){M},
                NULL /*X*/, &(Int) {M}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GELYAP", func) == 0 ) {
        FC_GLOBAL_(sla_gelyap,SLA_GELYAP)(factA, "N", &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*Q*/, &(Int){M},
                NULL /*X*/, &(Int) {M}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * GESTEIN
     */
    if ( strcasecmp("DLA_GESTEIN", func) == 0 ) {
        FC_GLOBAL_(dla_gestein,DLA_GESTEIN)(factA, "N", &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*Q*/, &(Int){M},
                NULL /*X*/, &(Int) {M}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GESTEIN", func) == 0 ) {
        FC_GLOBAL_(sla_gestein,SLA_GESTEIN)(factA, "N", &(Int){M}, NULL /*A*/, &(Int){M}, NULL /*Q*/, &(Int){M},
                NULL /*X*/, &(Int) {M}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }


    /*
     * GESYLV
     */
    if ( strcasecmp("DLA_GESYLV", func) == 0 ) {
        FC_GLOBAL_(dla_gesylv,DLA_GESYLV)(factA, factB, "N", "N", &(double){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M},  NULL /* B */, &(Int){N}, NULL /*QA*/, &(Int){M}, NULL /*QB*/, &(Int){N}, NULL /* X */, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GESYLV", func) == 0 ) {
        FC_GLOBAL_(sla_gesylv,SLA_GESYLV)(factA, factB, "N", "N", &(float){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M},  NULL /* B */, &(Int){N}, NULL /*QA*/, &(Int){M}, NULL /*QB*/, &(Int){N}, NULL /* X */, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * GESYLV2
     */
    if ( strcasecmp("DLA_GESYLV2", func) == 0 ) {
        FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)(factA, factB, "N", "N", &(double){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M},  NULL /* B */, &(Int){N}, NULL /*QA*/, &(Int){M}, NULL /*QB*/, &(Int){N}, NULL /* X */, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GESYLV2", func) == 0 ) {
        FC_GLOBAL_(sla_gesylv2,SLA_GESYLV2)(factA, factB, "N", "N", &(float){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M},  NULL /* B */, &(Int){N}, NULL /*QA*/, &(Int){M}, NULL /*QB*/, &(Int){N}, NULL /* X */, &(Int){M},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }


    /*
     * CGSYLV
     */
    if ( strcasecmp("DLA_GGCSYLV", func) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)(factA, factB, "N", "N",
                &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N},  NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*QA*/, &(Int){M}, NULL /*ZA*/, &(Int){M},  NULL /*QB*/, &(Int){N}, NULL /*ZB*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GGCSYLV", func) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv,SLA_GGCSYLV)(factA, factB, "N", "N",
                &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N},  NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*QA*/, &(Int){M}, NULL /*ZA*/, &(Int){M},  NULL /*QB*/, &(Int){N}, NULL /*ZB*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * CGSYLV DUAL
     */
    if ( strcasecmp("DLA_GGCSYLV_DUAL", func) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)(factA, factB, "N", "N",
                &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N},  NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*QA*/, &(Int){M}, NULL /*ZA*/, &(Int){M},  NULL /*QB*/, &(Int){N}, NULL /*ZB*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GGCSYLV_DUAL", func) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv_dual,SLA_GGCSYLV_DUAL)(factA, factB, "N", "N",
                &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N},  NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*QA*/, &(Int){M}, NULL /*ZA*/, &(Int){M},  NULL /*QB*/, &(Int){N}, NULL /*ZB*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }


    /*
     * GGSYLV
     */
    if ( strcasecmp("DLA_GGSYLV", func) == 0 ) {
       FC_GLOBAL_(dla_ggsylv,DLA_GGSYLV)(factA, factB, "N", "N",
                &(double) { 1 },  &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*QA*/, &(Int){M}, NULL /*ZA*/, &(Int){M}, NULL /*QB*/, &(Int){N}, NULL /*ZB*/, &(Int){N},
                NULL /*X*/, &(Int){M}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GGSYLV", func) == 0 ) {
       FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)(factA, factB, "N", "N",
                &(float) { 1 },  &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*QA*/, &(Int){M}, NULL /*ZA*/, &(Int){M}, NULL /*QB*/, &(Int){N}, NULL /*ZB*/, &(Int){N},
                NULL /*X*/, &(Int){M}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * GGLYAP
     */
    if ( strcasecmp("DLA_GGLYAP", func) == 0 ) {
        FC_GLOBAL_(dla_gglyap,DLA_GGLYAP)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M}, NULL /*C*/, &(Int){M}, NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int) {M},
                NULL /*X*/,  &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GGLYAP", func) == 0 ) {
        FC_GLOBAL_(sla_gglyap,SLA_GGLYAP)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M}, NULL /*C*/, &(Int){M}, NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int) {M},
                NULL /*X*/,  &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * GGSTEIN
     */
    if ( strcasecmp("DLA_GGSTEIN", func) == 0 ) {
        FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M}, NULL /*C*/, &(Int){M}, NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int) {M},
                NULL /*X*/,  &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }
    if ( strcasecmp("SLA_GGSTEIN", func) == 0 ) {
        FC_GLOBAL_(sla_ggstein,SLA_GGSTEIN)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M}, NULL /*C*/, &(Int){M}, NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int) {M},
                NULL /*X*/,  &(Int){M}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }



    /*
     * Refine GELYAP
     */
    if ( strcasecmp("DLA_GELYAP_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_gelyap_refine, DLA_GELYAP_REFINE) ("N", "N", &(Int){M},
                NULL /*A*/, &(Int) {M}, NULL /*X*/, &(Int){M},
                NULL /*Y*/, &(Int) {M},
                NULL /*AS*/, &(Int) {M}, NULL /*Q*/, &(Int){M},
                &(Int) {2}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GELYAP_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_gelyap_refine, SLA_GELYAP_REFINE) ("N", "N", &(Int){M},
                NULL /*A*/, &(Int) {M}, NULL /*X*/, &(Int){M},
                NULL /*Y*/, &(Int) {M},
                NULL /*AS*/, &(Int) {M}, NULL /*Q*/, &(Int){M},
                &(Int) {2}, &(float){1.0},NULL,  NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * Refine GESTEIN
     */
    if ( strcasecmp("DLA_GESTEIN_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_gestein_refine, DLA_GESTEIN_REFINE) ("N", "N", &(Int){M},
                NULL /*A*/, &(Int) {M}, NULL /*X*/, &(Int){M}, NULL, &(Int){M},
                NULL /*AS*/, &(Int) {M}, NULL /*Q*/, &(Int){M},
                &(Int) {2}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GESTEIN_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_gestein_refine, sla_gestein_refine) ("N","N", &(Int){M},
                NULL /*A*/, &(Int) {M}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int) {M}, NULL /*Q*/, &(Int){M},
                &(Int) {2}, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }


    /*
     * Refine GESYLV
     */
    if ( strcasecmp("DLA_GESYLV_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_gesylv_refine, DLA_GESYLV_REFINE) ("N", "N","N", &(double){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int){M}, NULL /*BS*/, &(Int){N}, NULL /*Q*/, &(Int){M}, NULL /*U*/, &(Int){N},
                &(Int){10}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GESYLV_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_gesylv_refine, SLA_GESYLV_REFINE) ("N", "N","N", &(float){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int){M}, NULL /*BS*/, &(Int){N}, NULL /*Q*/, &(Int){M}, NULL /*U*/, &(Int){N},
                &(Int){10}, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * Refine GESYLV2
     */
    if ( strcasecmp("DLA_GESYLV2_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_gesylv2_refine, DLA_GESYLV2_REFINE) ("N", "N","N", &(double){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int){M}, NULL /*BS*/, &(Int){N}, NULL /*Q*/, &(Int){M}, NULL /*U*/, &(Int){N},
                &(Int){10}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GESYLV2_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_gesylv2_refine, SLA_GESYLV2_REFINE) ("N", "N","N", &(float){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int){M}, NULL /*BS*/, &(Int){N}, NULL /*Q*/, &(Int){M}, NULL /*U*/, &(Int){N},
                &(Int){10}, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * GGLYAP
     */
    if ( strcasecmp("DLA_GGLYAP_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_gglyap_refine, DLA_GGLYAP_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M}, NULL /*B*/, &(Int){M}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int) {M}, NULL /*BS*/, &(Int) {M},NULL /*Q*/, &(Int) {M},NULL /*Z*/, &(Int) {M},
                &(Int) {10}, &(double) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GGLYAP_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_gglyap_refine, SLA_GGLYAP_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M}, NULL /*B*/, &(Int){M}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int) {M}, NULL /*BS*/, &(Int) {M},NULL /*Q*/, &(Int) {M},NULL /*Z*/, &(Int) {M},
                &(Int) {10}, &(float) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * GGSTEIN
     */
    if ( strcasecmp("DLA_GGSTEIN_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_ggstein_refine, DLA_GGSTEIN_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M}, NULL /*B*/, &(Int){M}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int) {M}, NULL /*BS*/, &(Int) {M},NULL /*Q*/, &(Int) {M},NULL /*Z*/, &(Int) {M},
                &(Int) {10}, &(double) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GGSTEIN_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_ggstein_refine, SLA_GGSTEIN_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M}, NULL /*B*/, &(Int){M}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*AS*/, &(Int) {M}, NULL /*BS*/, &(Int) {M},NULL /*Q*/, &(Int) {M},NULL /*Z*/, &(Int) {M},
                &(Int) {10}, &(float) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * GGSYLV
     */
    if ( strcasecmp("DLA_GGSYLV_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_ggsylv_refine, DLA_GGSYLV_REFINE) ("N", "N", "N",&(double) {1.0}, &(Int){M}, &(Int) {N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*As*/, &(Int){M}, NULL /*Bs*/, &(Int){N}, NULL /*Cs*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int){M}, NULL /*U*/, &(Int){N}, NULL /*V*/, &(Int){N},
                &(Int){ 10}, &(double) {1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GGSYLV_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_ggsylv_refine, SLA_GGSYLV_REFINE) ("N", "N", "N",&(float) {1.0}, &(Int){M}, &(Int) {N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N}, NULL /*X*/, &(Int){M},NULL, &(Int){M},
                NULL /*As*/, &(Int){M}, NULL /*Bs*/, &(Int){N}, NULL /*Cs*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int){M}, NULL /*U*/, &(Int){N}, NULL /*V*/, &(Int){N},
                &(Int){ 10}, &(float) {1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * CGSYLV
     */
    if ( strcasecmp("DLA_GGCSYLV_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv_refine, DLA_GGCSYLV_REFINE) ("N", "N", "N",&(double) { 1 }, &(double) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M},NULL, &(Int){M}, NULL, &(Int){M},
                NULL /*As*/, &(Int){M}, NULL /*Bs*/, &(Int){N}, NULL /*Cs*/, &(Int){M}, NULL /*Ds*/, &(Int){N},
                NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int){M}, NULL /*U*/, &(Int){N}, NULL /*V*/, &(Int){N},
                &(Int){ 10 }, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GGCSYLV_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv_refine, SLA_GGCSYLV_REFINE) ("N", "N", "N",&(float) { 1 }, &(float) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M},NULL, &(Int){M}, NULL, &(Int){M},
                NULL /*As*/, &(Int){M}, NULL /*Bs*/, &(Int){N}, NULL /*Cs*/, &(Int){M}, NULL /*Ds*/, &(Int){N},
                NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int){M}, NULL /*U*/, &(Int){N}, NULL /*V*/, &(Int){N},
                &(Int){ 10 }, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    /*
     * CGSYLV DUAL
     */
    if ( strcasecmp("DLA_GGCSYLV_DUAL_REFINE", func) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv_dual_refine, DLA_GGCSYLV_DUAL_REFINE) ("N", "N", "N",&(double) { 1 }, &(double) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M},NULL, &(Int){M}, NULL, &(Int){M},
                NULL /*As*/, &(Int){M}, NULL /*Bs*/, &(Int){N}, NULL /*Cs*/, &(Int){M}, NULL /*Ds*/, &(Int){N},
                NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int){M}, NULL /*U*/, &(Int){N}, NULL /*V*/, &(Int){N},
                &(Int){ 10 }, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }

    if ( strcasecmp("SLA_GGCSYLV_DUAL_REFINE", func) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv_dual_refine, SLA_GGCSYLV_DUAL_REFINE) ("N", "N", "N",&(float) { 1 }, &(float) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M}, NULL /*B*/, &(Int){N}, NULL /*C*/, &(Int){M}, NULL /*D*/, &(Int){N},
                NULL /*E*/, &(Int){M}, NULL /*F*/, &(Int){M},NULL, &(Int){M}, NULL, &(Int){M},
                NULL /*As*/, &(Int){M}, NULL /*Bs*/, &(Int){N}, NULL /*Cs*/, &(Int){M}, NULL /*Ds*/, &(Int){N},
                NULL /*Q*/, &(Int){M}, NULL /*Z*/, &(Int){M}, NULL /*U*/, &(Int){N}, NULL /*V*/, &(Int){N},
                &(Int){ 10 }, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) return (ssize_t) -1;
        return (ssize_t) ldwork;
    }




    return -1;

}


/**
 * @}
 */

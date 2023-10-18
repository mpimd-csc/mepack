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
 * \return memory requirements of the routine counted in floating point number of the desired precision.
 *    That means the return value either needs to be multiplied by `sizeof(float)` or `sizeof(double)`.
 *    A negative return value  indicates an error. This error code is explained in the description of
 *    the corresponding routine for which the memory query was performed.
 */
ssize_t mepack_memory_frontend(const char *func, const char * factA, const char *factB, int M, int N)
{
    Int info = 0;
    Int ldwork = -1;
    ssize_t return_value = -1;

    char *function_name_internal = NULL;

    if (strncasecmp(func, "DLA", 3) == 0 || strncasecmp(func, "SLA", 3) == 0)
    {
        function_name_internal = strdup(func);
    }
    else
    {
        if ( strcasecmp("mepack_double_gelyap", func ) == 0 ) {
            function_name_internal = strdup("DLA_GELYAP");
        } else if ( strcasecmp("mepack_single_gelyap", func ) == 0 ) {
            function_name_internal = strdup("SLA_GELYAP");
        } else if ( strcasecmp("mepack_double_gestein", func ) == 0 ) {
            function_name_internal = strdup("DLA_GESTEIN");
        } else if ( strcasecmp("mepack_single_gestein", func ) == 0 ) {
            function_name_internal = strdup("SLA_GESTEIN");
        } else if ( strcasecmp("mepack_double_gesylv", func ) == 0 ) {
            function_name_internal = strdup("DLA_GESYLV");
        } else if ( strcasecmp("mepack_single_gesylv", func ) == 0 ) {
            function_name_internal = strdup("SLA_GESYLV");
        } else if ( strcasecmp("mepack_double_gesylv2", func ) == 0 ) {
            function_name_internal = strdup("DLA_GESYLV2");
        } else if ( strcasecmp("mepack_single_gesylv2", func ) == 0 ) {
            function_name_internal = strdup("SLA_GESYLV2");
        } else if ( strcasecmp("mepack_double_gglyap", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGLYAP");
        } else if ( strcasecmp("mepack_single_gglyap", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGLYAP");
        } else if ( strcasecmp("mepack_double_ggstein", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGSTEIN");
        } else if ( strcasecmp("mepack_single_ggstein", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGSTEIN");
        } else if ( strcasecmp("mepack_double_ggsylv", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGSYLV");
        } else if ( strcasecmp("mepack_single_ggsylv", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGSYLV");
        } else if ( strcasecmp("mepack_double_ggcsylv", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGCSYLV");
        } else if ( strcasecmp("mepack_single_ggcsylv", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGCSYLV");
        } else if ( strcasecmp("mepack_double_ggcsylv_dual", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGCSYLV_DUAL");
        } else if ( strcasecmp("mepack_single_ggcsylv_dual", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGCSYLV_DUAL");
        } else if ( strcasecmp("mepack_double_gelyap_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GELYAP_REFINE");
        } else if ( strcasecmp("mepack_single_gelyap_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GELYAP_REFINE");
        } else if ( strcasecmp("mepack_double_gestein_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GESTEIN_REFINE");
        } else if ( strcasecmp("mepack_single_gestein_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GESTEIN_REFINE");
        } else if ( strcasecmp("mepack_double_gglyap_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGLYAP_REFINE");
        } else if ( strcasecmp("mepack_single_gglyap_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGLYAP_REFINE");
        } else if ( strcasecmp("mepack_double_ggstein_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGSTEIN_REFINE");
        } else if ( strcasecmp("mepack_single_ggstein_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGSTEIN_REFINE");
        } else if ( strcasecmp("mepack_double_gesylv_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GESYLV_REFINE");
        } else if ( strcasecmp("mepack_single_gesylv_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GESYLV_REFINE");
        } else if ( strcasecmp("mepack_double_gesylv2_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GESYLV2_REFINE");
        } else if ( strcasecmp("mepack_single_gesylv2_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GESYLV2_REFINE");
        } else if ( strcasecmp("mepack_double_ggsylv_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGSYLV_REFINE");
        } else if ( strcasecmp("mepack_single_ggsylv_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGSYLV_REFINE");
        } else if ( strcasecmp("mepack_double_ggcsylv_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGCSYLV_REFINE");
        } else if ( strcasecmp("mepack_single_ggcsylv_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGCSYLV_REFINE");
        } else if ( strcasecmp("mepack_double_ggcsylv_dual_refine", func ) == 0 ) {
            function_name_internal = strdup("DLA_GGCSYLV_DUAL_REFINE");
        } else if ( strcasecmp("mepack_single_ggcsylv_dual_refine", func ) == 0 ) {
            function_name_internal = strdup("SLA_GGCSYLV_DUAL_REFINE");
        }
    }
    /*
     *  GELYAP
     */
    if ( strcasecmp("DLA_GELYAP", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gelyap,DLA_GELYAP)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1},
                NULL /*X*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GELYAP", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gelyap,SLA_GELYAP)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1},
                NULL /*X*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * GESTEIN
     */
    if ( strcasecmp("DLA_GESTEIN", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gestein,DLA_GESTEIN)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1},
                NULL /*X*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GESTEIN", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gestein,SLA_GESTEIN)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1},
                NULL /*X*/, &(Int) {M+1}, NULL /*SCALE*/, NULL /* WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }


    /*
     * GESYLV
     */
    if ( strcasecmp("DLA_GESYLV", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gesylv,DLA_GESYLV)(factA, factB, "N", "N", &(double){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M+1},  NULL /* B */, &(Int){N+1}, NULL /*QA*/, &(Int){M+1}, NULL /*QB*/, &(Int){N+1}, NULL /* X */, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GESYLV", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gesylv,SLA_GESYLV)(factA, factB, "N", "N", &(float){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M+1},  NULL /* B */, &(Int){N+1}, NULL /*QA*/, &(Int){M+1}, NULL /*QB*/, &(Int){N+1}, NULL /* X */, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * GESYLV2
     */
    if ( strcasecmp("DLA_GESYLV2", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gesylv2,DLA_GESYLV2)(factA, factB, "N", "N", &(double){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M+1},  NULL /* B */, &(Int){N+1}, NULL /*QA*/, &(Int){M+1}, NULL /*QB*/, &(Int){N+1}, NULL /* X */, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GESYLV2", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gesylv2,SLA_GESYLV2)(factA, factB, "N", "N", &(float){1},&(Int) {M}, &(Int) {N},
                NULL /*A*/, &(Int) {M+1},  NULL /* B */, &(Int){N+1}, NULL /*QA*/, &(Int){M+1}, NULL /*QB*/, &(Int){N+1}, NULL /* X */, &(Int){M+1},
                NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }


    /*
     * CGSYLV
     */
    if ( strcasecmp("DLA_GGCSYLV", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv,DLA_GGCSYLV)(factA, factB, "N", "N",
                &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1},  NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*QA*/, &(Int){M+1}, NULL /*ZA*/, &(Int){M+1},  NULL /*QB*/, &(Int){N+1}, NULL /*ZB*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GGCSYLV", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv,SLA_GGCSYLV)(factA, factB, "N", "N",
                &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1},  NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*QA*/, &(Int){M+1}, NULL /*ZA*/, &(Int){M+1},  NULL /*QB*/, &(Int){N+1}, NULL /*ZB*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * CGSYLV DUAL
     */
    if ( strcasecmp("DLA_GGCSYLV_DUAL", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv_dual,DLA_GGCSYLV_DUAL)(factA, factB, "N", "N",
                &(double) {1}, &(double){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1},  NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*QA*/, &(Int){M+1}, NULL /*ZA*/, &(Int){M+1},  NULL /*QB*/, &(Int){N+1}, NULL /*ZB*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GGCSYLV_DUAL", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv_dual,SLA_GGCSYLV_DUAL)(factA, factB, "N", "N",
                &(float) {1}, &(float){1}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1},  NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*QA*/, &(Int){M+1}, NULL /*ZA*/, &(Int){M+1},  NULL /*QB*/, &(Int){N+1}, NULL /*ZB*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }


    /*
     * GGSYLV
     */
    if ( strcasecmp("DLA_GGSYLV", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggsylv,DLA_GGSYLV)(factA, factB, "N", "N",
                &(double) { 1 },  &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*QA*/, &(Int){M+1}, NULL /*ZA*/, &(Int){M+1}, NULL /*QB*/, &(Int){N+1}, NULL /*ZB*/, &(Int){N+1},
                NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GGSYLV", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggsylv,SLA_GGSYLV)(factA, factB, "N", "N",
                &(float) { 1 },  &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*QA*/, &(Int){M+1}, NULL /*ZA*/, &(Int){M+1}, NULL /*QB*/, &(Int){N+1}, NULL /*ZB*/, &(Int){N+1},
                NULL /*X*/, &(Int){M+1}, NULL /*SCALE*/,
                NULL /*WORK*/, &ldwork, &info, 1, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * GGLYAP
     */
    if ( strcasecmp("DLA_GGLYAP", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gglyap,DLA_GGLYAP)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*C*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int) {M+1},
                NULL /*X*/,  &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GGLYAP", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gglyap,SLA_GGLYAP)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*C*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int) {M+1},
                NULL /*X*/,  &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * GGSTEIN
     */
    if ( strcasecmp("DLA_GGSTEIN", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggstein,DLA_GGSTEIN)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*C*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int) {M+1},
                NULL /*X*/,  &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }
    if ( strcasecmp("SLA_GGSTEIN", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggstein,SLA_GGSTEIN)(factA, "N", &(Int){M},
                NULL /*A*/, &(Int){M+1}, NULL /*C*/, &(Int){M+1}, NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int) {M+1},
                NULL /*X*/,  &(Int){M+1}, NULL /*SCALE*/, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }



    /*
     * Refine GELYAP
     */
    if ( strcasecmp("DLA_GELYAP_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gelyap_refine, DLA_GELYAP_REFINE) ("N", "N", &(Int){M},
                NULL /*A*/, &(Int) {M+1}, NULL /*X*/, &(Int){M+1},
                NULL /*Y*/, &(Int) {M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*Q*/, &(Int){M+1},
                &(Int) {2}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GELYAP_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gelyap_refine, SLA_GELYAP_REFINE) ("N", "N", &(Int){M},
                NULL /*A*/, &(Int) {M+1}, NULL /*X*/, &(Int){M+1},
                NULL /*Y*/, &(Int) {M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*Q*/, &(Int){M+1},
                &(Int) {2}, &(float){1.0},NULL,  NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * Refine GESTEIN
     */
    if ( strcasecmp("DLA_GESTEIN_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gestein_refine, DLA_GESTEIN_REFINE) ("N", "N", &(Int){M},
                NULL /*A*/, &(Int) {M+1}, NULL /*X*/, &(Int){M+1}, NULL, &(Int){M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*Q*/, &(Int){M+1},
                &(Int) {2}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GESTEIN_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gestein_refine, sla_gestein_refine) ("N","N", &(Int){M},
                NULL /*A*/, &(Int) {M+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*Q*/, &(Int){M+1},
                &(Int) {2}, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }


    /*
     * Refine GESYLV
     */
    if ( strcasecmp("DLA_GESYLV_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gesylv_refine, DLA_GESYLV_REFINE) ("N", "N","N", &(double){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int){M+1}, NULL /*BS*/, &(Int){N+1}, NULL /*Q*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1},
                &(Int){10}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GESYLV_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gesylv_refine, SLA_GESYLV_REFINE) ("N", "N","N", &(float){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int){M+1}, NULL /*BS*/, &(Int){N+1}, NULL /*Q*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1},
                &(Int){10}, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * Refine GESYLV2
     */
    if ( strcasecmp("DLA_GESYLV2_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gesylv2_refine, DLA_GESYLV2_REFINE) ("N", "N","N", &(double){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int){M+1}, NULL /*BS*/, &(Int){N+1}, NULL /*Q*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1},
                &(Int){10}, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GESYLV2_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gesylv2_refine, SLA_GESYLV2_REFINE) ("N", "N","N", &(float){1.0}, &(Int){M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int){M+1}, NULL /*BS*/, &(Int){N+1}, NULL /*Q*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1},
                &(Int){10}, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * GGLYAP
     */
    if ( strcasecmp("DLA_GGLYAP_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_gglyap_refine, DLA_GGLYAP_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*BS*/, &(Int) {M+1},NULL /*Q*/, &(Int) {M+1},NULL /*Z*/, &(Int) {M+1},
                &(Int) {10}, &(double) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GGLYAP_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_gglyap_refine, SLA_GGLYAP_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*BS*/, &(Int) {M+1},NULL /*Q*/, &(Int) {M+1},NULL /*Z*/, &(Int) {M+1},
                &(Int) {10}, &(float) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * GGSTEIN
     */
    if ( strcasecmp("DLA_GGSTEIN_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggstein_refine, DLA_GGSTEIN_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*BS*/, &(Int) {M+1},NULL /*Q*/, &(Int) {M+1},NULL /*Z*/, &(Int) {M+1},
                &(Int) {10}, &(double) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GGSTEIN_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggstein_refine, SLA_GGSTEIN_REFINE) ( "N", "N",&(Int) {M},
                NULL /*A*/, &(Int) {M+1}, NULL /*B*/, &(Int){M+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*AS*/, &(Int) {M+1}, NULL /*BS*/, &(Int) {M+1},NULL /*Q*/, &(Int) {M+1},NULL /*Z*/, &(Int) {M+1},
                &(Int) {10}, &(float) {10}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * GGSYLV
     */
    if ( strcasecmp("DLA_GGSYLV_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggsylv_refine, DLA_GGSYLV_REFINE) ("N", "N", "N",&(double) {1.0}, &(Int){M}, &(Int) {N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*As*/, &(Int){M+1}, NULL /*Bs*/, &(Int){N+1}, NULL /*Cs*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1}, NULL /*V*/, &(Int){N+1},
                &(Int){ 10}, &(double) {1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GGSYLV_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggsylv_refine, SLA_GGSYLV_REFINE) ("N", "N", "N",&(float) {1.0}, &(Int){M}, &(Int) {N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1}, NULL /*X*/, &(Int){M+1},NULL, &(Int){M+1},
                NULL /*As*/, &(Int){M+1}, NULL /*Bs*/, &(Int){N+1}, NULL /*Cs*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1}, NULL /*V*/, &(Int){N+1},
                &(Int){ 10}, &(float) {1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * CGSYLV
     */
    if ( strcasecmp("DLA_GGCSYLV_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv_refine, DLA_GGCSYLV_REFINE) ("N", "N", "N",&(double) { 1 }, &(double) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1},NULL, &(Int){M+1}, NULL, &(Int){M+1},
                NULL /*As*/, &(Int){M+1}, NULL /*Bs*/, &(Int){N+1}, NULL /*Cs*/, &(Int){M+1}, NULL /*Ds*/, &(Int){N+1},
                NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1}, NULL /*V*/, &(Int){N+1},
                &(Int){ 10 }, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GGCSYLV_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv_refine, SLA_GGCSYLV_REFINE) ("N", "N", "N",&(float) { 1 }, &(float) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1},NULL, &(Int){M+1}, NULL, &(Int){M+1},
                NULL /*As*/, &(Int){M+1}, NULL /*Bs*/, &(Int){N+1}, NULL /*Cs*/, &(Int){M+1}, NULL /*Ds*/, &(Int){N+1},
                NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1}, NULL /*V*/, &(Int){N+1},
                &(Int){ 10 }, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    /*
     * CGSYLV DUAL
     */
    if ( strcasecmp("DLA_GGCSYLV_DUAL_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(dla_ggcsylv_dual_refine, DLA_GGCSYLV_DUAL_REFINE) ("N", "N", "N",&(double) { 1 }, &(double) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1},NULL, &(Int){M+1}, NULL, &(Int){M+1},
                NULL /*As*/, &(Int){M+1}, NULL /*Bs*/, &(Int){N+1}, NULL /*Cs*/, &(Int){M+1}, NULL /*Ds*/, &(Int){N+1},
                NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1}, NULL /*V*/, &(Int){N+1},
                &(Int){ 10 }, &(double){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);
        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }

    if ( strcasecmp("SLA_GGCSYLV_DUAL_REFINE", function_name_internal) == 0 ) {
        FC_GLOBAL_(sla_ggcsylv_dual_refine, SLA_GGCSYLV_DUAL_REFINE) ("N", "N", "N",&(float) { 1 }, &(float) {1}, &(Int) {M}, &(Int){N},
                NULL /*A*/, &(Int){M+1}, NULL /*B*/, &(Int){N+1}, NULL /*C*/, &(Int){M+1}, NULL /*D*/, &(Int){N+1},
                NULL /*E*/, &(Int){M+1}, NULL /*F*/, &(Int){M+1},NULL, &(Int){M+1}, NULL, &(Int){M+1},
                NULL /*As*/, &(Int){M+1}, NULL /*Bs*/, &(Int){N+1}, NULL /*Cs*/, &(Int){M+1}, NULL /*Ds*/, &(Int){N+1},
                NULL /*Q*/, &(Int){M+1}, NULL /*Z*/, &(Int){M+1}, NULL /*U*/, &(Int){N+1}, NULL /*V*/, &(Int){N+1},
                &(Int){ 10 }, &(float){1.0}, NULL, NULL /*WORK*/, &ldwork, &info, 1, 1, 1);

        if ( info < 0 ) {
            return_value = -1;
            goto end;
        }
        return_value = (ssize_t) ldwork;
        goto end;
    }


end:
    free(function_name_internal);
    return return_value;

}


/**
 * @}
 */

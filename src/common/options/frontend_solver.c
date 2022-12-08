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
 * @addtogroup frontend_solver
 * @{
 */

/**
 * @brief Set the frontend solver for TRSYLV (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_trsylv_frontend_solver_set function sets the frontend solver for the
 * TRSYLV equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 * @see MEPACK_FRONTEND_SOLVER_LAPACK
 *
 */
void mepack_trsylv_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,trsylv_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRSYLV_FRONTEND_SOLVER_SET)( & _FS );
}

/**
 * @brief Set the frontend solver for TRSYLV2 (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_trsylv2_frontend_solver_set function sets the frontend solver for the
 * TRSYLV2 equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 *
 */
void mepack_trsylv2_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,trsylv2_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRSYLV2_FRONTEND_SOLVER_SET)( & _FS );
}

/**
 * @brief Set the frontend solver for TRLYAP (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_trlyap_frontend_solver_set function sets the frontend solver for the
 * TRLYAP equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 *
 */
void mepack_trlyap_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,trlyap_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRLYAP_FRONTEND_SOLVER_SET)( & _FS );
}

/**
 * @brief Set the frontend solver for TGLYAP (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_tglyap_frontend_solver_set function sets the frontend solver for the
 * TGLYAP equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 *
 */
void mepack_tglyap_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,tglyap_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGLYAP_FRONTEND_SOLVER_SET)( & _FS );
}

/**
 * @brief Set the frontend solver for TRSTEIN (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_trstein_frontend_solver_set function sets the frontend solver for the
 * TRSTEIN equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 *
 */
void mepack_trstein_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,trstein_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TRSTEIN_FRONTEND_SOLVER_SET)( & _FS );
}

/**
 * @brief Set the frontend solver for TGSTEIN (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_tgstein_frontend_solver_set function sets the frontend solver for the
 * TGSTEIN equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 *
 */
void mepack_tgstein_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,tgstein_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGSTEIN_FRONTEND_SOLVER_SET)( & _FS );
}


/**
 * @brief Set the frontend solver for TGSYLV (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_tgsylv_frontend_solver_set function sets the frontend solver for the
 * TGSYLV equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 * @see MEPACK_FRONTEND_SOLVER_GARDINER_LAUB
 *
 */
void mepack_tgsylv_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,tgsylv_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGSYLV_FRONTEND_SOLVER_SET)( & _FS );
}


/**
 * @brief Set the frontend solver for TGCSYLV (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_tgcsylv_frontend_solver_set function sets the frontend solver for the
 * TGCSYLV equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 *
 */
void mepack_tgcsylv_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,tgcsylv_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGCSYLV_FRONTEND_SOLVER_SET)( & _FS );
}

/**
 * @brief Set the frontend solver for TGCSYLV_DUAL (C Interface)
 * @param[in]   FS      inner solver to set
 *
 * The mepack_tgcsylv_dual_frontend_solver_set function sets the frontend solver for the
 * TGCSYLV_DUAL equation.
 *
 * @see MEPACK_FRONTEND_SOLVER_LEVEL3
 * @see MEPACK_FRONTEND_SOLVER_LEVEL2
 * @see MEPACK_FRONTEND_SOLVER_DAG
 * @see MEPACK_FRONTEND_SOLVER_2STAGE
 * @see MEPACK_FRONTEND_SOLVER_RECURSIVE
 *
 */
void mepack_tgcsylv_dual_frontend_solver_set(int FS)
{
    Int _FS = FS;
    FC_MODULE_(mepack_options_frontend_solver,tgcsylv_dual_frontend_solver_set,MEPACK_OPTIONS_FRONTEND_SOLVER,TGCSYLV_DUAL_FRONTEND_SOLVER_SET)( & _FS );
}




/**
 * @}
 */

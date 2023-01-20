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
 * @addtogroup isolver
 * @{
 */

/**
 * @brief Set the inner solver for TRSYLV (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_trsylv_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRSYLV equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_trsylv_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,trsylv_isolver_set,MEPACK_OPTIONS_ISOLVER,TRSYLV_ISOLVER_SET)( & _IS );
}

/**
 * @brief Set the inner solver for TRSYLV2 (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_trsylv2_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRSYLV2 equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_trsylv2_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,trsylv2_isolver_set,MEPACK_OPTIONS_ISOLVER,TRSYLV2_ISOLVER_SET)( & _IS );
}

/**
 * @brief Set the inner solver for TGSYLV (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_tgsylv_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRSYLV equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_tgsylv_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,tgsylv_isolver_set,MEPACK_OPTIONS_ISOLVER,TGSYLV_ISOLVER_SET)( & _IS );
}

/**
 * @brief Set the inner solver for TGCSYLV (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_tgcsylv_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRSYLV equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_tgcsylv_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,tgcsylv_isolver_set,MEPACK_OPTIONS_ISOLVER,TGCSYLV_ISOLVER_SET)( & _IS );
}

/**
 * @brief Set the inner solver for TGCSYLV_DUAL (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_tgcsylv_dual_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRSYLV equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_tgcsylv_dual_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,tgcsylv_dual_isolver_set,MEPACK_OPTIONS_ISOLVER,TGCSYLV_DUAL_ISOLVER_SET)( & _IS );
}


/**
 * @brief Set the inner solver for TRLYAP (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_trlyap_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRLYAP equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_trlyap_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,trlyap_isolver_set,MEPACK_OPTIONS_ISOLVER,TRLYAP_ISOLVER_SET)( & _IS );
}

/**
 * @brief Set the inner solver for TGLYAP (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_trlyap_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRLYAP equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_tglyap_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,tglyap_isolver_set,MEPACK_OPTIONS_ISOLVER,TGLYAP_ISOLVER_SET)( & _IS );
}

/**
 * @brief Set the inner solver for TRSTEIN (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_trstein_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRSTEIN equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_trstein_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,trstein_isolver_set,MEPACK_OPTIONS_ISOLVER,TRSTEIN_ISOLVER_SET)( & _IS );
}

/**
 * @brief Set the inner solver for TGSTEIN (C Interface)
 * @param[in]   IS      inner solver to set
 *
 * The mepack_trstein_isolver_set function set the inner solver for the
 * level 3 and DAG solvers for the TRSTEIN equation.
 *
 * @see MEPACK_ISOLVER_DEFAULT
 * @see MEPACK_ISOLVER_LOCAL_COPY_ALIGN
 * @see MEPACK_ISOLVER_LOCAL_COPY
 * @see MEPACK_ISOLVER_REORDER
 * @see MEPACK_ISOLVER_LEVEL2
 * @see MEPACK_ISOLVER_LEVEL2_NAIV
 * @see MEPACK_ISOLVER_RECURSIVE
 *
 */
void mepack_tgstein_isolver_set(int IS)
{
    Int _IS = IS;
    FC_MODULE_(mepack_options_isolver,tgstein_isolver_set,MEPACK_OPTIONS_ISOLVER,TGSTEIN_ISOLVER_SET)( & _IS );
}


/**
 * @}
 */

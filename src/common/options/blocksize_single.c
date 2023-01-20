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
 * @addtogroup blocksize
 * @{
 */



/**
 * @brief Set the row blocksize for TRSYLV
 * @param[in]  mb       New row blocksize
 *
 * The mepack_single_trsylv_blocksize_mb_set function sets the row blocksize for the TRSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trsylv_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trsylv_blocksize_mb_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TRSYLV
 * @param[in]  nb       New column blocksize
 *
 * The mepack_single_trsylv_blocksize_mb_set function sets the column blocksize for the TRSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trsylv_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_single,trsylv_blocksize_nb_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TRSYLV2
 * @param[in]  mb       New row blocksize
 *
 * The mepack_single_trsylv2_blocksize_mb_set function sets the row blocksize for the TRSYLV2 level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trsylv2_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trsylv2_blocksize_mb_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV2_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TRSYLV2
 * @param[in]  nb       New column blocksize
 *
 * The mepack_single_trsylv2_blocksize_mb_set function sets the column blocksize for the TRSYLV2 level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trsylv2_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_single,trsylv2_blocksize_nb_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV2_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TGSYLV
 * @param[in]  mb       New row blocksize
 *
 * The mepack_single_tgsylv_blocksize_mb_set function sets the row blocksize for the TGSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgsylv_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgsylv_blocksize_mb_set,  MEPACK_BLOCKSIZE_SINGLE,TGSYLV_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TGSYLV
 * @param[in]  nb       New column blocksize
 *
 * The mepack_single_tgsylv_blocksize_mb_set function sets the column blocksize for the TGSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgsylv_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_single,tgsylv_blocksize_nb_set,  MEPACK_BLOCKSIZE_SINGLE,TGSYLV_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TGCSYLV
 * @param[in]  mb       New row blocksize
 *
 * The mepack_single_tgcsylv_blocksize_mb_set function sets the row blocksize for the TGCSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgcsylv_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgcsylv_blocksize_mb_set,  MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TGCSYLV
 * @param[in]  nb       New column blocksize
 *
 * The mepack_single_tgcsylv_blocksize_mb_set function sets the column blocksize for the TGCSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgcsylv_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_single,tgcsylv_blocksize_nb_set,  MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TGCSYLV_DUAL
 * @param[in]  mb       New row blocksize
 *
 * The mepack_single_tgcsylv_dual_blocksize_mb_set function sets the row blocksize for the TGCSYLV_DUAL level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgcsylv_dual_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgcsylv_dual_blocksize_mb_set,  MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_DUAL_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TGCSYLV_DUAL
 * @param[in]  nb       New column blocksize
 *
 * The mepack_single_tgcsylv_dual_blocksize_nb_set function sets the column blocksize for the TGCSYLV_DUAL level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgcsylv_dual_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_single,tgcsylv_dual_blocksize_nb_set,  MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_DUAL_BLOCKSIZE_NB_SET)( &_NB);
    return;
}


/**
 * @brief Set the blocksize for TRLYAP
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_trlyap_blocksize_set function sets the blocksize for the TRLYAP level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trlyap_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trlyap_blocksize_set,  MEPACK_BLOCKSIZE_SINGLE,TRLYAP_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGLYAP
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_tglyap_blocksize_set function sets the blocksize for the TGLYAP level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tglyap_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tglyap_blocksize_set,  MEPACK_BLOCKSIZE_SINGLE,TGLYAP_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRSTEIN
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_trstein_blocksize_set function sets the blocksize for the TRSTEIN level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trstein_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trstein_blocksize_set,  MEPACK_BLOCKSIZE_SINGLE,TRSTEIN_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGSTEIN
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_tgstein_blocksize_set function sets the blocksize for the TGSTEIN level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgstein_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgstein_blocksize_set,  MEPACK_BLOCKSIZE_SINGLE,TGSTEIN_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRLYAP level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_trlyap_blocksize_2stage_set function sets the blocksize for the TRLYAP level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trlyap_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trlyap_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TRLYAP_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGLYAP level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_tglyap_blocksize_2stage_set function sets the blocksize for the TGLYAP level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tglyap_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tglyap_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TGLYAP_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRSTEIN level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_trstein_blocksize_2stage_set function sets the blocksize for the TRSTEIN level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trstein_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trstein_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TRSTEIN_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGSTEIN level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_tgstein_blocksize_2stage_set function sets the blocksize for the TGSTEIN level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgstein_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgstein_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TGSTEIN_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}


/**
 * @brief Set the blocksize for TRSYLV level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_trsylv_blocksize_2stage_set function sets the blocksize for the TRSYLV level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trsylv_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trsylv_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRSYLV2 level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_trsylv2_blocksize_2stage_set function sets the blocksize for the TRSYLV2 level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_trsylv2_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,trsylv2_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TRSYLV2_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGSYLV level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_tgsylv_blocksize_2stage_set function sets the blocksize for the TGSYLV level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgsylv_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgsylv_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TGSYLV_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGCSYLV level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_tgcsylv_blocksize_2stage_set function sets the blocksize for the TGCSYLV level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgcsylv_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgcsylv_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGCSYLV_DUAL level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_single_tgcsylv_dual_blocksize_2stage_set function sets the blocksize for the TGCSYLV_DUAL level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_single_tgcsylv_dual_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_single,tgcsylv_dual_blocksize_2stage_set,  MEPACK_BLOCKSIZE_SINGLE,TGCSYLV_DUAL_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}



/**
 * @}
 */

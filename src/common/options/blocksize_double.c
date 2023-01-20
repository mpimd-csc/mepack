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
 * The mepack_double_trsylv_blocksize_mb_set function sets the row blocksize for the TRSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trsylv_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trsylv_blocksize_mb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TRSYLV
 * @param[in]  nb       New column blocksize
 *
 * The mepack_double_trsylv_blocksize_mb_set function sets the column blocksize for the TRSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trsylv_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_double,trsylv_blocksize_nb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TRSYLV2
 * @param[in]  mb       New row blocksize
 *
 * The mepack_double_trsylv2_blocksize_mb_set function sets the row blocksize for the TRSYLV2 level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trsylv2_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trsylv2_blocksize_mb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV2_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TRSYLV2
 * @param[in]  nb       New column blocksize
 *
 * The mepack_double_trsylv2_blocksize_mb_set function sets the column blocksize for the TRSYLV2 level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trsylv2_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_double,trsylv2_blocksize_nb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV2_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TGSYLV
 * @param[in]  mb       New row blocksize
 *
 * The mepack_double_tgsylv_blocksize_mb_set function sets the row blocksize for the TGSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgsylv_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgsylv_blocksize_mb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSYLV_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TGSYLV
 * @param[in]  nb       New column blocksize
 *
 * The mepack_double_tgsylv_blocksize_mb_set function sets the column blocksize for the TGSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgsylv_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_double,tgsylv_blocksize_nb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSYLV_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TGCSYLV
 * @param[in]  mb       New row blocksize
 *
 * The mepack_double_tgcsylv_blocksize_mb_set function sets the row blocksize for the TGCSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgcsylv_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgcsylv_blocksize_mb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TGCSYLV
 * @param[in]  nb       New column blocksize
 *
 * The mepack_double_tgcsylv_blocksize_nb_set function sets the column blocksize for the TGCSYLV level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgcsylv_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_double,tgcsylv_blocksize_nb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_BLOCKSIZE_NB_SET)( &_NB);
    return;
}

/**
 * @brief Set the row blocksize for TGCSYLV_DUAL
 * @param[in]  mb       New row blocksize
 *
 * The mepack_double_tgcsylv_dual_blocksize_mb_set function sets the row blocksize for the TGCSYLV_DUAL level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgcsylv_dual_blocksize_mb_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgcsylv_dual_blocksize_mb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_DUAL_BLOCKSIZE_MB_SET)( &_MB);
    return;
}

/**
 * @brief Set the column blocksize for TGCSYLV_DUAL
 * @param[in]  nb       New column blocksize
 *
 * The mepack_double_tgcsylv_dual_blocksize_nb_set function sets the column blocksize for the TGCSYLV_DUAL level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgcsylv_dual_blocksize_nb_set ( int nb ) {
    Int _NB = nb;
    FC_MODULE_(mepack_options_blocksize_double,tgcsylv_dual_blocksize_nb_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_DUAL_BLOCKSIZE_NB_SET)( &_NB);
    return;
}


/**
 * @brief Set the blocksize for TRLYAP
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_trlyap_blocksize_set function sets the blocksize for the TRLYAP level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trlyap_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trlyap_blocksize_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRLYAP_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGLYAP
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_tglyap_blocksize_set function sets the blocksize for the TGLYAP level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tglyap_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tglyap_blocksize_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGLYAP_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRSTEIN
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_trstein_blocksize_set function sets the blocksize for the TRSTEIN level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trstein_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trstein_blocksize_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSTEIN_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGSTEIN
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_tgstein_blocksize_set function sets the blocksize for the TGSTEIN level-3
 * and DAG algorithms. If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgstein_blocksize_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgstein_blocksize_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSTEIN_BLOCKSIZE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRLYAP level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_trlyap_blocksize_2stage_set function sets the blocksize for the TRLYAP level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trlyap_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trlyap_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRLYAP_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGLYAP level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_tglyap_blocksize_2stage_set function sets the blocksize for the TGLYAP level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tglyap_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tglyap_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGLYAP_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRSTEIN level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_trstein_blocksize_2stage_set function sets the blocksize for the TRSTEIN level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trstein_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trstein_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSTEIN_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGSTEIN level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_tgstein_blocksize_2stage_set function sets the blocksize for the TGSTEIN level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgstein_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgstein_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSTEIN_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}


/**
 * @brief Set the blocksize for TRSYLV level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_trsylv_blocksize_2stage_set function sets the blocksize for the TRSYLV level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trsylv_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trsylv_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TRSYLV2 level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_trsylv2_blocksize_2stage_set function sets the blocksize for the TRSYLV2 level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_trsylv2_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,trsylv2_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TRSYLV2_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGSYLV level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_tgsylv_blocksize_2stage_set function sets the blocksize for the TGSYLV level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgsylv_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgsylv_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGSYLV_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGCSYLV level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_tgcsylv_blocksize_2stage_set function sets the blocksize for the TGCSYLV level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgcsylv_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgcsylv_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}

/**
 * @brief Set the blocksize for TGCSYLV_DUAL level-3 two-stage solver.
 * @param[in]  mb       New  blocksize
 *
 * The mepack_double_tgcsylv_dual_blocksize_2stage_set function sets the blocksize for the TGCSYLV_DUAL level-3
 * two-stage solver.
 * If the value is set to 0 the automatic selection is used.
 *
 */
void mepack_double_tgcsylv_dual_blocksize_2stage_set ( int mb ) {
    Int _MB = mb;
    FC_MODULE_(mepack_options_blocksize_double,tgcsylv_dual_blocksize_2stage_set,  MEPACK_OPTIONS_BLOCKSIZE_DOUBLE,TGCSYLV_DUAL_BLOCKSIZE_2STAGE_SET)( &_MB);
    return;
}





/**
 * @}
 */

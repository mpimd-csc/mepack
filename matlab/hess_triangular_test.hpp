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


#pragma once

#include "triangular_test.hpp"

template<typename T>
static bool mepack_mexoct_check_hess_triangular(mexoct::ArrayType<T>& A, mexoct::ArrayType<T> &B)
{
    bool b_triag = mepack_mexoct_check_triangular(B);

    if (b_triag && A.rows <= 2 ) return true;

    ssize_t i, j;

    for ( j = 0; j < A.columns-2; j++) {
        for ( i = j+2; i < A.rows; i++) {
            if ( A(i,j) != T(0.0) )
                return false;
        }
    }
    return b_triag;
}


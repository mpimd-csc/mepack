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


#pragma once

template<typename T>
static bool mepack_mexoct_check_qtriangular(mexoct::ArrayType<T>& mat)
{
    ssize_t i, j;

    i = 0;
    j = 0;

    if ( mat.rows != mat.columns) return false;
    for ( j = 0; j < mat.columns; j++) {
        for ( i = j+2; i < mat.rows; i++) {
            if ( mat(i,j) != T(0.0)) return false;
        }
    }
    i = 0;
    for ( i = 0; i < mat.rows; i++) {
        if ( i < mat.rows-1 && mat(i+1,i) != T(0.0)){
            if ( i < mat.rows -2 && mat(i+2,i+1) != T(0.0)) return false;
        }
    }
    return true;
}


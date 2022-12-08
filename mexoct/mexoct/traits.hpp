/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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
 * Copyright (C) Martin Koehler, 2019
 */

#pragma once

#include <iostream>
#ifdef MEXOCT_OCTAVE
#define ROWS(a)  (a.rows())
#define COLS(a)  (a.columns())
#define ARG_T      const octave_value &
#else
#define ROWS(a)  (mxGetM(a))
#define COLS(a)  (mxGetN(a))
#define ARG_T    const mxArray *
#endif

namespace mexoct {
    namespace Traits {

        struct Vector : public Trait {
            static bool isValid(ARG_T  _array){
                if (( ROWS(_array)  == 1 && COLS(_array) >= 1)
                        || (COLS(_array)  == 1 && ROWS(_array) >= 1))
                    return true;
                else
                    return false;

            }
            static constexpr const char* NAME = "Vector";
        };

        struct RowVector : public Trait {
            static bool isValid(ARG_T  _array){
                if (( ROWS(_array)  == 1 && COLS(_array) >= 1))
                    return true;
                else
                    return false;

            }
            static constexpr const char* NAME = "Row Vector";
        };

        struct ColVector : public Trait {
            static bool isValid(ARG_T  _array){
                if (( COLS(_array)  == 1 && ROWS(_array) >= 1))
                    return true;
                else
                    return false;
            }
            static constexpr const char* NAME = "Column Vector";
        };

        struct SquareMatrix : public Trait {
            static bool isValid(ARG_T _array) {
                if ( COLS(_array) == ROWS(_array)) {
                    // std::cerr << " cols : " << COLS(_array) << "\n";
                    // std::cerr << " rows : " << ROWS(_array) << "\n";
                    return true;
                } else {
                    return false;
                }
            }
            static constexpr const char* NAME = "Square Matrix";
        };

        struct TallSkinnyMatrix : public Trait {
            static bool isValid(ARG_T _array) {
                return ROWS(_array) >= COLS(_array);
            }
            static constexpr const char* NAME = "Tall-Skinny Matrix";
        };
    }
}

#undef ROWS
#undef COLS
#undef ARG_T

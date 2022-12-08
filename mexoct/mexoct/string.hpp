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

namespace mexoct {


    template<>
        struct ParamType<std::string>
        {
#ifdef MEXOCT_MATLAB
            static bool isValid(const mxArray *_array)
            {
#if 0
                fprintf(stderr, "empty %d\n", mxIsEmpty(_array));
                fprintf(stderr, "ischar %d\n", mxIsChar(_array));
                fprintf(stderr, "M %d\n", mxGetM(_array));
                fprintf(stderr, "N %d\n", mxGetN(_array));
                fprintf(stderr, "ID %d\n", mxGetClassID(_array));
                fprintf(stderr, "String? %d\n", mxIsClass(_array, "string"));
#endif
                // check type through class_name
#ifdef MEXOCT_MATLAB_DQ_STRING
          return ! mxIsEmpty(_array)
                    && (((mxIsChar(_array))
                    && (mxGetM(_array) == 1 && mxGetN(_array) >= 1 )) || mxIsClass(_array, "string") );

#else
          return ! mxIsEmpty(_array)
                    && (((mxIsChar(_array))
                    && (mxGetM(_array) == 1 && mxGetN(_array) >= 1 )));
#endif

                  }

            static std::string copy_in(const mxArray *_array)
            {
                char * str;
#ifdef MEXOCT_MATLAB_DQ_STRING
                if ( mxIsClass(_array, "string")) {
                    // This is required since the C API does not allow to access the double quoted strings from MATLAB directly
                    mxArray *string_class[1];
                    mxArray *char_array[1];
                    string_class[0] = const_cast<mxArray*> (_array);
                    mexCallMATLAB(1, char_array, 1, string_class, "char");
                    str = mxArrayToString(char_array[0]);
                    fprintf(stderr, "str %s\n", str);
                } else {
                    str = mxArrayToString(_array);
                }
#else
                str = mxArrayToString(_array);
#endif
                std::string s(str);
                mxFree(str);
                return s;
            }

            static mxArray * copy_out(std::string _array)
            {
                mxArray * mx = mxCreateString(_array.c_str());
                return mx;
            }

            static constexpr const char* NAME = "string";

#else

            static bool isValid(const octave_value & _array)
            {
                // check type through class_name
                //
#if OCTAVE_MAJOR_VERSION == 4 && OCTAVE_MINOR_VERSION < 4
                return !_array.is_empty()
                    && (_array.is_sq_string() || _array.is_dq_string())
                    && ( (_array.rows() == 1 && _array.columns() >= 1 ));
#else
                // std::cout << "empty: " << _array.isempty() << "\n";
                // std::cout << "sq: " << _array.is_sq_string() << "\n";
                // std::cout << "dq:"  << _array.is_dq_string() << "\n";
                // std::cout << _array.rows() << " x "  << _array.columns() << "\n";
                return !_array.isempty()
                    && (_array.is_sq_string() || _array.is_dq_string())
                    && ( (_array.rows() == 1 && _array.columns() >= 1 ));
#endif

            }

            static std::string copy_in(const octave_value & _array)
            {
                return _array.string_value();
            }

            static octave_value copy_out(std::string _array)
            {
                return octave_value(_array);
            }

            static constexpr const char* NAME = "string";


#endif
        };


}

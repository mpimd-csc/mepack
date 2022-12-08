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

namespace mexoct {
    template<typename T>
        struct ScalarType
        {
            ScalarType() = default;
            ScalarType(const ScalarType& _ost)
                : value ( _ost.value), cmplx(_ost.cmplx)
            {

            }
            ScalarType(const T & v) : value(v), cmplx(false) {}
            operator T () const {
                return value;
            }
            T & operator=(T & other) {
                value = other.value;
                return value;
            }
            std::string type_name () {
                return std::string(ParamType<ScalarType<T>>::NAME);
            }

            T value;
            bool cmplx;

        };


    template<typename T>
        struct ParamType<ScalarType<T>>
        {
#ifdef MEXOCT_MATLAB
            static bool isValid(const mxArray* _array)
            {
                // check type through ClassId
                return !mxIsEmpty(_array)       // Empty is void.
                    && mxIsScalar(_array)       // Check if scalar
                    && mxGetClassID(_array) == MatlabClassId<T>::value // Valid base type
                    && !mxIsComplex(_array); // is not complex
            }

            static ScalarType<T> copy_in(const mxArray* _array)
            {
                ScalarType<T> result;
                result.value = *(reinterpret_cast<T*>(mxGetData(_array)));
                result.cmplx = false;
                return result;
            }

            static mxArray* copy_out(const ScalarType<T>& _array)
            {
                mxArray * ret;
                mwSize dims[1] = {1};
                ret = mxCreateNumericArray(1, dims, MatlabClassId<T>::value, mxREAL);
                *(reinterpret_cast<T*>(mxGetData(ret))) = _array.value;
                return ret;
            }

            static constexpr const char* NAME = MatlabClassId<T>::name;
#else

            static bool isValid(const octave_value & _array)
            {
                // octave_stdout << "class name: " << _array.class_name() << " type_id = " << _array.type_id() << "\n";
                // check type through class_name
#if OCTAVE_MAJOR_VERSION == 4 && OCTAVE_MINOR_VERSION < 4
                return !_array.is_empty()
                    && _array.is_real_scalar()
                    && _array.class_name() == OctaveType<T>::class_name;
#else
                return !_array.isempty()
                    && _array.is_real_scalar()
                    && _array.class_name() == OctaveType<T>::class_name;
#endif

            }

            static ScalarType<T> copy_in(const octave_value & _array)
            {
                ScalarType<T> result;
                result.value = OctaveType<T>::get(_array);
                result.cmplx = false;
                return result;
            }

            static octave_value copy_out(const ScalarType<T>& _array)
            {
                // octave_stdout << "copy out " << OctaveType<T>::type_name << " " <<_array.value << "\n";
                return octave_value(_array.value);
            }

            static constexpr const char* NAME = OctaveType<T>::type_name;


#endif
        };


    template<typename T>
        struct ParamType<ScalarType<std::complex<T>>>
        {
#ifdef MEXOCT_MATLAB
            static bool isValid(const mxArray* _array)
            {
                // check type through ClassId
                return !mxIsEmpty(_array)       // Empty is void.
                    && mxIsScalar(_array)       // Check if scalar
                    && mxGetClassID(_array) == MatlabClassId<std::complex<T>>::value // Valid base type
                    && mxIsComplex(_array);  // T complex;
            }

            static ScalarType<std::complex<T>> copy_in(const mxArray* _array)
            {
                ScalarType<std::complex<T>> result;
#if MX_HAS_INTERLEAVED_COMPLEX
                if ( std::is_same<T, double>::value ) {
                    result.value = *(reinterpret_cast<std::complex<T>*>(mxGetComplexDoubles(_array)));
                } else if (  std::is_same<T, float>::value ) {
                    result.value = *(reinterpret_cast<std::complex<T>*>(mxGetComplexSingles(_array)));
                }
#else
                std::complex<T> v(*(reinterpret_cast<T*>(mxGetPr(_array))), *(reinterpret_cast<T*>(mxGetPi(_array))));
                result.value = v;
#endif
                result.cmplx = true;

                return result;
            }

            static mxArray* copy_out(const ScalarType<std::complex<T>>& _array)
            {
                mxArray * ret;
                mwSize dims[1] = {1};
                ret = mxCreateNumericArray(1, dims, MatlabClassId<std::complex<T>>::value, mxCOMPLEX);
#if MX_HAS_INTERLEAVED_COMPLEX
                if ( std::is_same<T, double>::value ) {
                    *(reinterpret_cast<std::complex<T>*>(mxGetComplexDoubles(ret))) = _array.value;
                } else if (  std::is_same<T, float>::value ) {
                    *(reinterpret_cast<std::complex<T>*>(mxGetComplexSingles(ret))) = _array.value;
                }
#else
                T * re = reinterpret_cast<T*>(mxGetPr(ret));
                T * im = reinterpret_cast<T*>(mxGetPi(ret));
                *re = static_cast<T>( _array.value.real());
                *im = static_cast<T>( _array.value.imag());
#endif

                return ret;
            }

            static constexpr const char* NAME = MatlabClassId<std::complex<T>>::name;
#else
            static bool isValid(const octave_value & _array)
            {
                // octave_stdout << "class name: " << _array.class_name() << " type_id = " << _array.type_id() << "\n";
#if OCTAVE_MAJOR_VERSION == 4 && OCTAVE_MINOR_VERSION < 4
                return !_array.is_empty()
                    && _array.is_complex_scalar()
                    && _array.class_name() == OctaveType<T>::class_name;
#else
                return !_array.isempty()
                    && _array.is_complex_scalar()
                    && _array.class_name() == OctaveType<T>::class_name;

#endif

            }
            static ScalarType<std::complex<T>> copy_in(const octave_value & _array)
            {
                ScalarType<std::complex<T>> result;
                result.value = OctaveType<std::complex<T>>::get(_array);
                result.cmplx = true;
                return result;

            }

            static octave_value copy_out(const ScalarType<std::complex<T>>& _array)
            {
                return octave_value(_array.value);
            }
            static constexpr const char* NAME = OctaveType<std::complex<T>>::type_name;
#endif
        };

}

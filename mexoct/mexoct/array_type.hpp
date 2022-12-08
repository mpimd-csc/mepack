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

#include <memory>
#include <cstring>
#include <iostream>

namespace mexoct {

    template<typename T>
        struct ArrayType
        {
            std::shared_ptr<T> values;
            ssize_t rows;
            ssize_t columns;
            ssize_t ld;
            bool cmplx;

            ArrayType(const ArrayType& _oat)
                : values ( _oat.values),
                rows   ( _oat.rows),
                columns (_oat.columns),
                ld     ( _oat.ld),
                cmplx(_oat.cmplx)
            {
            }
            ArrayType(const ssize_t m, const ssize_t n, T * v, bool c = false)
                : rows(m), columns(n), ld(m), cmplx(c)
            {
                values.reset(v,  std::default_delete<T[]>());
            }

            ArrayType(const ssize_t m, const ssize_t n, bool c = false) : rows(m), columns(n), ld(m), cmplx(c)
            {
                values.reset(new T[m*n], std::default_delete<T[]>());
            }


            // ArrayType(ArrayType&& _oat)
            //     : values ( _oat.values),
            //     rows   ( _oat.rows),
            //     columns (_oat.columns),
            //     ld     ( _oat.ld),
            //     cmplx(_oat.cmplx)
            // {
            //         _oat.values.reset();
            // }
            //


            operator T* () const {
                return this->values.get();
            }
            T& operator[](ssize_t idx) {
                T *v = this->values.get();
                return v[idx];
            }
            const T operator[](ssize_t idx) const {
                T *v = this->values.get();
                return v[idx];
            }

            T& operator()(ssize_t i, ssize_t j) {
                T *v = this->values.get();
                return v[i+j*ld];

            }
            const T operator()(ssize_t i, ssize_t j) const {
                T *v = this->values.get();
                return v[i+j*ld];
            }
            T * ptr () const {
                return this->values.get();
            }

            // T *values;
        };


    template<typename T>
        struct ParamType<ArrayType<T>>
        {
#ifdef MEXOCT_MATLAB
            static bool isValid(const mxArray* _array)
            {
                return !mxIsEmpty(_array)       // Empty is void.
                    && mxGetClassID(_array) == MatlabClassId<T>::value // Valid base type
                    && static_cast<bool>(mxIsComplex(_array)) == MatlabClassId<T>::cmplx
                    && !mxIsSparse(_array) ;

            }

            static ArrayType<T> copy_in(const mxArray* _array)
            {
                //mwSize m = mxGetM(_array);
                //mwSize n = mxGetN(_array);
                ssize_t m = mxGetM(_array);
                ssize_t n = mxGetN(_array);
                bool cpx = mxIsComplex(_array);


                if ( cpx ) {
#if MX_HAS_INTERLEAVED_COMPLEX
                    T * _xval = new T [m * n];
                    if ( std::is_same<T,std::complex<double>>::value ) {
                        std::memcpy(_xval, reinterpret_cast<T*>(mxGetComplexDoubles(_array)), sizeof(T)* m * n);
                    } else if (  std::is_same<T, std::complex<float>>::value ) {
                        std::memcpy(_xval, reinterpret_cast<T*>(mxGetComplexSingles(_array)), sizeof(T)* m * n);
                    }
                    return ArrayType<T>(m , n , reinterpret_cast<T*>(_xval), cpx);

#else
                    if ( std::is_same<T,std::complex<double>>::value ) {
                        double *re, *im;
                        ssize_t k,l;
                        re = (double *) mxGetPr(_array);
                        im = (double *) mxGetPi(_array);
                        std::complex<double> * _xval = new std::complex<double> [m * n];
                        for (l = 0; l < n; l++) {
                            for (k = 0; k < m; k++) {
                                _xval[k+l*m] = std::complex<double>(re[k+l*m], im[k+l*m]);
                            }
                        }
                        return ArrayType<T>(m , n , reinterpret_cast<T*>(_xval), cpx);
                    } else if (  std::is_same<T, std::complex<float>>::value ) {
                        float *re, *im;
                        ssize_t k,l;
                        re = (float  *) mxGetPr(_array);
                        im = (float  *) mxGetPi(_array);
                        std::complex<float> * _xval = new std::complex<float> [m * n];

                        for (l = 0; l < n; l++) {
                            for (k = 0; k < m; k++) {
                                _xval[k+l*m] = std::complex<float>(re[k+l*m], im[k+l*m]);
                            }
                        }
                        return ArrayType<T>(m , n , reinterpret_cast<T*>(_xval), cpx);
                    } else {
                        return ArrayType<T>(0,0);
                    }
#endif
                } else {
                    T * _xval = new T [m * n];
                    std::memcpy ( _xval, mxGetData(_array), sizeof(T) * m * n );
                    return ArrayType<T>(m , n , reinterpret_cast<T*>(_xval), cpx);
                }
                return ArrayType<T>(0,0);
            }

            static mxArray* copy_out(const ArrayType<T>& _array)
            {
                mxComplexity cflag = ( MatlabClassId<T>::cmplx ) ? mxCOMPLEX: mxREAL;
                mxArray* result = mxCreateNumericMatrix(_array.rows, _array.columns, MatlabClassId<T>::value, cflag);

                if ( cflag == mxREAL ) {
                    if ( _array.rows == _array.ld) {
                        std::memcpy(reinterpret_cast<T*>(mxGetData(result)), _array.ptr(), sizeof(T) * _array.rows * _array.columns);
                    } else {
                        ssize_t k;
                        T * ptr = reinterpret_cast<T*>(mxGetData(result));
                        T * inptr = _array.ptr();
                        for (k = 0; k < _array.columns; k++) {
                            std::memcpy( ptr, _array.ptr(), sizeof(T)* _array.rows );
                            ptr += _array.ld;
                            inptr += _array.ld;
                        }

                    }
                } else {
                    ssize_t m = _array.rows;
                    ssize_t n = _array.columns;
                    ssize_t ld = _array.ld;
#if MX_HAS_INTERLEAVED_COMPLEX
                    if ( std::is_same<T,std::complex<double>>::value ) {
                        if ( m == ld ) {
                            std::memcpy(reinterpret_cast<T*>(mxGetComplexDoubles(result)), _array.ptr(), sizeof(T) * _array.rows * _array.columns);
                        } else {
                            ssize_t k;
                            T * ptr = reinterpret_cast<T*>(mxGetComplexDoubles(result));
                            T * inptr = _array.ptr();
                            for (k = 0; k < _array.columns; k++) {
                                std::memcpy( ptr, _array.ptr(), sizeof(T)* _array.rows );
                                ptr += _array.ld;
                                inptr += _array.ld;
                            }
                        }
                    } else if (  std::is_same<T, std::complex<float>>::value ) {
                        if ( m == ld ) {
                            std::memcpy(reinterpret_cast<T*>(mxGetComplexSingles(result)), _array.ptr(), sizeof(T) * _array.rows * _array.columns);
                        } else {
                            ssize_t k;
                            T * ptr = reinterpret_cast<T*>(mxGetComplexSingles(result));
                            T * inptr = _array.ptr();
                            for (k = 0; k < _array.columns; k++) {
                                std::memcpy( ptr, _array.ptr(), sizeof(T)* _array.rows );
                                ptr += _array.ld;
                                inptr += _array.ld;
                            }
                        }

                    }

#else
                    if ( std::is_same<T,std::complex<double>>::value ) {
                        ssize_t k, l;
                        double * re = (double *) mxGetPr(result);
                        double * im = (double *) mxGetPi(result);
                        std::complex<double> * _xval = reinterpret_cast<std::complex<double>*>(_array.ptr());

                        for (l = 0; l < n; l++) {
                            for (k = 0; k < m; k++) {
                                re[k+l*m] = _xval[k+l*ld].real();
                                im[k+l*m] = _xval[k+l*ld].imag();

                            }
                        }
                    } else if (  std::is_same<T, std::complex<float>>::value ) {
                        ssize_t k, l;
                        float * re = (float *) mxGetPr(result);
                        float * im = (float *) mxGetPi(result);
                        std::complex<float> * _xval = reinterpret_cast<std::complex<float>*>(_array.ptr());

                        for (l = 0; l < n; l++) {
                            for (k = 0; k < m; k++) {
                                re[k+l*m] = _xval[k+l*ld].real();
                                im[k+l*m] = _xval[k+l*ld].imag();

                            }
                        }

                    }
#endif
                }
                return result ;

            }

            static constexpr const char* NAME = MatlabClassId<T>::name;

#else
            /* Octave part */

            static bool isValid(const octave_value & _array)
            {
                // octave_stdout << "class name: " << _array.class_name() << " type_id = " << _array.type_id() << "\n";
                // check type through class_name
#if OCTAVE_MAJOR_VERSION == 4 && OCTAVE_MINOR_VERSION < 4
                return !_array.is_empty()
                    && _array.class_name() == OctaveType<T>::class_name
                    && _array.is_complex_type() == OctaveType<T>::cmplx
                    && !_array.is_sparse_type();
#else
                return !_array.isempty()
                    && _array.class_name() == OctaveType<T>::class_name
                    && _array.iscomplex() == OctaveType<T>::cmplx
                    && !_array.issparse();
#endif

            }

            static ArrayType<T> copy_in(const octave_value & _array)
            {
                OctaveArrayType<T> wrapper;
                typename OctaveArrayType<T>::type r = wrapper.get(_array);
                T * _xvals = new T [_array.rows() * _array.columns()];
                std::memcpy(_xvals, r.fortran_vec(), sizeof(T) * _array.rows() * _array.columns());
                ArrayType<T> result(r.rows(), r.columns(), _xvals, false);
                return result;
            }

            static octave_value copy_out(const ArrayType<T> & array)
            {
                dim_vector dv(array.rows, array.columns);
                typename OctaveArrayType<T>::type _val (dv);
                if ( array.rows != array.ld ) {
                    ssize_t k;
                    T * ptr = reinterpret_cast<T*>(_val.fortran_vec());
                    T * inptr = array.ptr();
                    for (k = 0; k < array.columns; k++) {
                        std::memcpy( ptr, array.ptr(), sizeof(T)* array.rows );
                        ptr += array.ld;
                        inptr += array.ld;
                    }
                } else {
                    std::memcpy(_val.fortran_vec(), array.ptr(), sizeof(T)* array.rows * array.columns);
                }
                return octave_value(_val);
            }

            static constexpr const char* NAME = OctaveType<T>::class_name;


#endif
        };



}

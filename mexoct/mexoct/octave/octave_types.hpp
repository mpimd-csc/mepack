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

#include <cstdint>
#include <type_traits>
#include <vector>
#include <map>
#include <complex>
#include <string>

namespace mexoct {

    struct OctaveTypeBase
    {
        static constexpr const char * class_name = "";
        static constexpr const char * type_name = "";

        static constexpr bool      cmplx = false;
    };

    template<typename T>
        struct OctaveType : public OctaveTypeBase
    {
        static T get(const octave_value & arg) {
            return T();
        }
    };

    template<>
        struct OctaveType<int8_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "int8";
            static constexpr const char * type_name = "int8_t";

            static int8_t get(const octave_value & arg) {
                return arg.int8_scalar_value();
            }
        };

    template<>
        struct OctaveType<uint8_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "uint8";
            static constexpr const char * type_name = "uint8_t";

            static int8_t get(const octave_value & arg) {
                return arg.uint8_scalar_value();
            }
        };

    template<>
        struct OctaveType<int16_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "int16";
            static constexpr const char * type_name = "int16_t";

            static int16_t get(const octave_value & arg) {
                return arg.int16_scalar_value();
            }
        };

    template<>
        struct OctaveType<uint16_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "uint16";
            static constexpr const char * type_name = "uint16_t";

            static int16_t get(const octave_value & arg) {
                return arg.uint16_scalar_value();
            }
        };

    template<>
        struct OctaveType<int32_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "int32";
            static constexpr const char * type_name = "int32_t";

            static int32_t get(const octave_value & arg) {
                return arg.int32_scalar_value();
            }
        };

    template<>
        struct OctaveType<uint32_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "uint32";
            static constexpr const char * type_name = "uint32_t";

            static int32_t get(const octave_value & arg) {
                return arg.uint32_scalar_value();
            }
        };

    template<>
        struct OctaveType<int64_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "int64";
            static constexpr const char * type_name = "int64_t";

            static int64_t get(const octave_value & arg) {
                return arg.int64_scalar_value();
            }
        };

    template<>
        struct OctaveType<uint64_t> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "uint64";
            static constexpr const char * type_name = "uint64_t";

            static int64_t get(const octave_value & arg) {
                return arg.uint64_scalar_value();
            }
        };

    template<>
        struct OctaveType<float> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "single";
            static constexpr const char * type_name = "float";
            static float get(const octave_value & arg) {
                return arg.float_scalar_value();
            }
        };

    template<>
        struct OctaveType<double> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "double";
            static constexpr const char * type_name = "double";

            static double get(const octave_value & arg) {
                return arg.scalar_value();
            }
        };

    template<>
        struct OctaveType<bool> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "logical";
            static constexpr const char * type_name = "bool";

            static bool  get(const octave_value & arg) {
                return arg.bool_value();
            }
        };

    template<>
        struct OctaveType<std::complex<float>> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "single";
            static constexpr const char * type_name = "float complex";
            static constexpr bool cmplx = true;

            static std::complex<float> get(const octave_value & arg) {
                return arg.float_complex_value();
            }
        };

    template<>
        struct OctaveType<std::complex<double>> : public OctaveTypeBase
        {
            static constexpr const char * class_name = "double";
            static constexpr const char * type_name = "double complex";
            static constexpr bool cmplx = true;

            static std::complex<double> get(const octave_value & arg) {
                return arg.complex_value();
            }
        };


    template<typename oct_type>
        struct OctaveArrayType{
            using type = void*;
            type get (const octave_value & octval) {
                return NULL;
            }
            OctaveArrayType() = default;

        };

    template<>
        struct OctaveArrayType<double>
        {
            OctaveArrayType() = default;
            using type = Matrix;
            type get(const octave_value & octval) {
                return octval.matrix_value();
            }
        };

    template<>
        struct OctaveArrayType<float>
        {
            OctaveArrayType() = default;
            using type = FloatMatrix;
            type get(const octave_value & octval) {
                return octval.float_matrix_value();
            }
        };

    template<>
        struct OctaveArrayType<int8_t>
        {
            OctaveArrayType() = default;
            using type = int8NDArray;
            type get(const octave_value & octval) {
                return octval.int8_array_value();
            }
        };

    template<>
        struct OctaveArrayType<uint8_t>
        {
            OctaveArrayType() = default;
            using type = uint8NDArray;
            type get(const octave_value & octval) {
                return octval.uint8_array_value();
            }
        };

    template<>
        struct OctaveArrayType<int16_t>
        {
            OctaveArrayType() = default;
            using type = int16NDArray;
            type get(const octave_value & octval) {
                return octval.int16_array_value();
            }
        };

    template<>
        struct OctaveArrayType<uint16_t>
        {
            OctaveArrayType() = default;
            using type = uint16NDArray;
            type get(const octave_value & octval) {
                return octval.uint16_array_value();
            }
        };

    template<>
        struct OctaveArrayType<int32_t>
        {
            OctaveArrayType() = default;
            using type = int32NDArray;
            type get(const octave_value & octval) {
                return octval.int32_array_value();
            }
        };

    template<>
        struct OctaveArrayType<uint32_t>
        {
            OctaveArrayType() = default;
            using type = uint32NDArray;
            type get(const octave_value & octval) {
                return octval.uint32_array_value();
            }
        };

    template<>
        struct OctaveArrayType<int64_t>
        {
            OctaveArrayType() = default;
            using type = int64NDArray;
            type get(const octave_value & octval) {
                return octval.int64_array_value();
            }
        };

    template<>
        struct OctaveArrayType<uint64_t>
        {
            OctaveArrayType() = default;
            using type = uint64NDArray;
            type get(const octave_value & octval) {
                return octval.uint64_array_value();
            }
        };



    template<>
        struct OctaveArrayType<std::complex<double>>
        {
            OctaveArrayType() = default;
            using type = ComplexMatrix;
            type get(const octave_value & octval) {
                return octval.complex_matrix_value();
            }
        };

    template<>
        struct OctaveArrayType<std::complex<float>>
        {
            OctaveArrayType() = default;
            using type = FloatComplexMatrix;
            type get(const octave_value & octval) {
                return octval.float_complex_matrix_value();
            }
        };


}

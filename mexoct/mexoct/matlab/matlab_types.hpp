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

#include "mexoct/paramtypes.hpp"
#include <cstdint>
#include <type_traits>
#include <vector>
#include <map>
#include <complex>

namespace mexoct {
    template<typename T>
        struct MatlabClassIdBase
        {
            static constexpr mxClassID value = mxUNKNOWN_CLASS;
            static constexpr bool      cmplx = false;
            static constexpr const char *    name  = "scalar";
        };

    /*
     * Possible mxClassIDs
     *
     */
    template<typename T>
        struct MatlabClassId : public MatlabClassIdBase<T>
    {
        // condition is dependent on T to only trigger on specialization.
        static_assert(std::is_void<T>::value, "Unknown type");

        static constexpr mxClassID value = mxUNKNOWN_CLASS;
    };

    /* Integers  */
    template<>
        struct MatlabClassId<int8_t> : public MatlabClassIdBase<int8_t>
        {
            static constexpr mxClassID value = mxINT8_CLASS;
            static constexpr const char *    name  = "int8_t";

        };

    template<>
        struct MatlabClassId<uint8_t> : public MatlabClassIdBase<uint8_t>
        {
            static constexpr mxClassID value = mxUINT8_CLASS;
            static constexpr const char *    name  = "uint8_t";

        };

    template<>
        struct MatlabClassId<int16_t> : public MatlabClassIdBase<int16_t>
        {
            static constexpr mxClassID value = mxINT16_CLASS;
            static constexpr const char *    name  = "int16_t";

        };

    template<>
        struct MatlabClassId<uint16_t> : public MatlabClassIdBase<uint16_t>
        {
            static constexpr mxClassID value = mxUINT16_CLASS;
            static constexpr const char *    name  = "uint16_t";

        };

    template<>
        struct MatlabClassId<int32_t> : public MatlabClassIdBase<int32_t>
        {
            static constexpr mxClassID value = mxINT32_CLASS;
            static constexpr const char *    name  = "int32_t";

        };

    template<>
        struct MatlabClassId<uint32_t> : public MatlabClassIdBase<uint32_t>
        {
            static constexpr mxClassID value = mxUINT32_CLASS;
            static constexpr const char *    name  = "uint32_t";

        };


    template<>
        struct MatlabClassId<int64_t> : public MatlabClassIdBase<int64_t>
        {
            static constexpr mxClassID value = mxINT64_CLASS;
            static constexpr const char *    name  = "int64_t";

        };

    template<>
        struct MatlabClassId<uint64_t> : public MatlabClassIdBase<uint64_t>
        {
            static constexpr mxClassID value = mxUINT64_CLASS;
            static constexpr const char *    name  = "uint64_t";

        };


    /* Float types  */
    template<>
        struct MatlabClassId<double> : public MatlabClassIdBase<double>
        {
            static constexpr mxClassID value = mxDOUBLE_CLASS;
            static constexpr const char *    name  = "double";

        };

    template<>
        struct MatlabClassId<float> : public MatlabClassIdBase<float>
        {
            static constexpr mxClassID value = mxSINGLE_CLASS;
            static constexpr const char *    name  = "float";

        };

    template<>
        struct MatlabClassId<std::complex<double>> : public MatlabClassIdBase<std::complex<double>>
        {
            static constexpr mxClassID value = mxDOUBLE_CLASS;
            static constexpr bool      cmplx = true;
            static constexpr const char *    name  = "double complex";


        };

    template<>
        struct MatlabClassId<std::complex<float>> : public MatlabClassIdBase<std::complex<float>>
        {
            static constexpr mxClassID value = mxSINGLE_CLASS;
            static constexpr bool      cmplx = true;
            static constexpr const char *    name  = "float complex";

        };


    template<>
        struct MatlabClassId<bool> : public MatlabClassIdBase<bool>
        {
            static constexpr mxClassID value = mxLOGICAL_CLASS;
            static constexpr const char *    name  = "bool";

        };

    template<>
        struct MatlabClassId<char> : public MatlabClassIdBase<char>
        {
            static constexpr mxClassID value = mxCHAR_CLASS;
            static constexpr const char *    name  = "char";

        };


}

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

#include "mexoct/interface.hpp"

using namespace mexoct;

template<typename T>
std::tuple<T> copy_val(T val)
{
    return std::make_tuple(val);
}



MEXOCT_ENTRY(copy_input,
"test_type_info  -- Check the type information returned by mexoct\n"\
"\n"\
"    str = test_type_info(var) \n"\
"\n"\
"    str contains a string with the type name of var\n")
{
    MEXOCT_INIT();
    auto paramInt8      = Parameter<ScalarType<int8_t>>("var");
    auto paramUInt8     = Parameter<ScalarType<uint8_t>>("var");
    auto paramInt16     = Parameter<ScalarType<int16_t>>("var");
    auto paramUInt16    = Parameter<ScalarType<uint16_t>>("var");
    auto paramInt32     = Parameter<ScalarType<int32_t>>("var");
    auto paramUInt32    = Parameter<ScalarType<uint32_t>>("var");
    auto paramInt64     = Parameter<ScalarType<int64_t>>("var");
    auto paramUInt64    = Parameter<ScalarType<uint64_t>>("var");
    auto paramFloat     = Parameter<ScalarType<float>>("var");
    auto paramDouble    = Parameter<ScalarType<double>>("var");
    auto paramFloatCpx  = Parameter<ScalarType<std::complex<float>>>("var");
    auto paramDoubleCpx = Parameter<ScalarType<std::complex<double>>>("var");

    auto paramArrInt8   = Parameter<ArrayType<int8_t>>("var");
    auto paramArrUInt8  = Parameter<ArrayType<uint8_t>>("var");
    auto paramArrInt16  = Parameter<ArrayType<int16_t>>("var");
    auto paramArrUInt16 = Parameter<ArrayType<uint16_t>>("var");
    auto paramArrInt32  = Parameter<ArrayType<int32_t>>("var");
    auto paramArrUInt32 = Parameter<ArrayType<uint32_t>>("var");
    auto paramArrInt64  = Parameter<ArrayType<int64_t>>("var");
    auto paramArrUInt64 = Parameter<ArrayType<uint64_t>>("var");
    auto paramArrFloat  = Parameter<ArrayType<float>>("var");
    auto paramArrDouble = Parameter<ArrayType<double>>("var");
    auto paramArrFloatCpx  = Parameter<ArrayType<std::complex<float>>>("var");
    auto paramArrDoubleCpx = Parameter<ArrayType<std::complex<double>>>("var");


    auto fInt8   = makeFunctionExt<1>(copy_val<ScalarType<int8_t>>, paramInt8);
    auto fUInt8  = makeFunctionExt<1>(copy_val<ScalarType<uint8_t>>, paramUInt8);
    auto fInt16  = makeFunctionExt<1>(copy_val<ScalarType<int16_t>>, paramInt16);
    auto fUInt16 = makeFunctionExt<1>(copy_val<ScalarType<uint16_t>>, paramUInt16);
    auto fInt32  = makeFunctionExt<1>(copy_val<ScalarType<int32_t>>, paramInt32);
    auto fUInt32 = makeFunctionExt<1>(copy_val<ScalarType<uint32_t>>, paramUInt32);
    auto fInt64  = makeFunctionExt<1>(copy_val<ScalarType<int64_t>>, paramInt64);
    auto fUInt64 = makeFunctionExt<1>(copy_val<ScalarType<uint64_t>>,paramUInt64);

    auto fDouble = makeFunctionExt<1>(copy_val<ScalarType<double>>, paramDouble);
    auto fFloat  = makeFunctionExt<1>(copy_val<ScalarType<float>>, paramFloat);
    auto fDoubleCpx = makeFunctionExt<1>(copy_val<ScalarType<std::complex<double>>>, paramDoubleCpx);
    auto fFloatCpx  = makeFunctionExt<1>(copy_val<ScalarType<std::complex<float>>>, paramFloatCpx);

    auto fArrInt8 = makeFunctionExt<1>(copy_val<ArrayType<int8_t>>, paramArrInt8);
    auto fArrUInt8 = makeFunctionExt<1>(copy_val<ArrayType<uint8_t>>, paramArrUInt8);
    auto fArrInt16 = makeFunctionExt<1>(copy_val<ArrayType<int16_t>>, paramArrInt16);
    auto fArrUInt16 = makeFunctionExt<1>(copy_val<ArrayType<uint16_t>>, paramArrUInt16);
    auto fArrInt32 = makeFunctionExt<1>(copy_val<ArrayType<int32_t>>, paramArrInt32);
    auto fArrUInt32 = makeFunctionExt<1>(copy_val<ArrayType<uint32_t>>, paramArrUInt32);
    auto fArrInt64 = makeFunctionExt<1>(copy_val<ArrayType<int64_t>>, paramArrInt64);
    auto fArrUInt64 = makeFunctionExt<1>(copy_val<ArrayType<uint64_t>>, paramArrUInt64);

    auto fArrDouble = makeFunctionExt<1>(copy_val<ArrayType<double>>, paramArrDouble);
    auto fArrFloat  = makeFunctionExt<1>(copy_val<ArrayType<float>>, paramArrFloat);
    auto fArrDoubleCpx = makeFunctionExt<1>(copy_val<ArrayType<std::complex<double>>>, paramArrDoubleCpx);
    auto fArrFloatCpx  = makeFunctionExt<1>(copy_val<ArrayType<std::complex<float>>>, paramArrFloatCpx);


    auto parser = makeArgParser("copy_input", fInt8,fUInt8, fInt16, fUInt16, fInt32, fUInt32, fInt64, fUInt64,
                                              fDouble, fDoubleCpx, fFloat, fFloatCpx,
                                              fArrInt8,fArrUInt8, fArrInt16, fArrUInt16, fArrInt32, fArrUInt32, fArrInt64, fArrUInt64,
                                              fArrDouble, fArrDoubleCpx, fArrFloat, fArrFloatCpx
                                              );


    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;


}


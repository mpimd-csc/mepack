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
std::tuple<std::string> type_name(ScalarType<T> val)
{
    std::string retval(val.type_name());
    return std::make_tuple(retval);
}


MEXOCT_ENTRY(type_info,
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


    auto fInt8   = makeFunctionExt<1>(type_name<int8_t>,  paramInt8);
    auto fUInt8  = makeFunctionExt<1>(type_name<uint8_t>,  paramUInt8);
    auto fInt16  = makeFunctionExt<1>(type_name<int16_t>,  paramInt16);
    auto fUInt16 = makeFunctionExt<1>(type_name<uint16_t>, paramUInt16);
    auto fInt32  = makeFunctionExt<1>(type_name<int32_t>,  paramInt32);
    auto fUInt32 = makeFunctionExt<1>(type_name<uint32_t>, paramUInt32);
    auto fInt64  = makeFunctionExt<1>(type_name<int64_t>,  paramInt64);
    auto fUInt64 = makeFunctionExt<1>(type_name<uint64_t>, paramUInt64);

    auto fDouble    = makeFunctionExt<1>(type_name<double>, paramDouble);
    auto fDoubleCpx = makeFunctionExt<1>(type_name<std::complex<double>>, paramDoubleCpx);
    auto fFloat     = makeFunctionExt<1>(type_name<float>, paramFloat);
    auto fFloatCpx  = makeFunctionExt<1>(type_name<std::complex<float>>, paramFloatCpx);

    auto parser = makeArgParser("type_info",fInt8,fUInt8, fInt16, fUInt16, fInt32, fUInt32, fInt64, fUInt64,
                                            fDouble, fDoubleCpx, fFloat, fFloatCpx);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;


}

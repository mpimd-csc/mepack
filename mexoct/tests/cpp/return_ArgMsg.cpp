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

std::tuple<>
out0_and_in0()
{
    return std::make_tuple();
}

std::tuple<ScalarType<int8_t>>
out1_and_in0()
{
    return std::make_tuple(ScalarType<int8_t>(0));
}

std::tuple<>
out0_and_in1(ScalarType<int8_t> A)
{
    return std::make_tuple();
}

std::tuple<ScalarType<int8_t>>
out1_and_in1(ScalarType<int8_t> A)
{
    return std::make_tuple(ScalarType<int8_t>(0));
}

std::tuple<ScalarType<int8_t>,ScalarType<int8_t>>
out2_and_in0()
{
    return std::make_tuple(ScalarType<int8_t>(0), ScalarType<int8_t>(0));
}

std::tuple<>
out0_and_in2(ScalarType<int8_t> A, ScalarType<int8_t> B)
{
    return std::make_tuple();
}

std::tuple<ScalarType<int8_t>, ScalarType<int8_t>>
out2_and_in1(ScalarType<int8_t> A)
{
    return std::make_tuple(ScalarType<int8_t>(0), ScalarType<int8_t>(0));
}

std::tuple<ScalarType<int8_t>>
out1_and_in2(ScalarType<int8_t> A, ScalarType<int8_t> B)
{
    return std::make_tuple(ScalarType<int8_t>(0));
}

std::tuple<ScalarType<int8_t>, ScalarType<int8_t>>
out2_and_in2(ScalarType<int8_t> A, ScalarType<int8_t> B)
{
    return std::make_tuple(ScalarType<int8_t>(0), ScalarType<int8_t>(0));
}

MEXOCT_ENTRY(return_ArgMsg,
"test_return_ArgMsg  -- Check the user defined return argument message on function mismatch by mexoct\n"\
"\n"\
"    out = test_return_ArgMsg(in) \n"\
"\n")
{
    MEXOCT_INIT();

    /*-----------------------------------------------------------------------------
     *  parameters
     *-----------------------------------------------------------------------------*/
    auto paramInt8 = Parameter<ScalarType<int8_t>>("var");

    /*-----------------------------------------------------------------------------
     *  MATLAB/OCTAVE calls
     *-----------------------------------------------------------------------------*/
    //auto call_out0_in0_customArgMsg = makeFunctionExtArgMsg<0>(out0_and_in0, false, Utils::string("NOTHING_FROM_NOTHING"));
    auto call_out0_in0_customArgMsg = makeFunctionExt<0>(out0_and_in0, Utils::string("NOTHING_FROM_NOTHING"));
    auto call_out0_in0_defaultArgMsg = makeFunctionExt<0>(out0_and_in0);

    //auto call_out1_in0_customArgMsg = makeFunctionExtArgMsg<1>(out1_and_in0, true, Utils::string("ONE_FROM_NOTHING"));
    auto call_out1_in0_customArgMsg = makeFunctionExt<1>(out1_and_in0, Utils::string("ONE_FROM_NOTHING"));
    auto call_out1_in0_defaultArgMsg = makeFunctionExt<1>(out1_and_in0);

    //auto call_out0_in1_customArgMsg = makeFunctionExtArgMsg<0>(out0_and_in1, false, Utils::string("NOTHING_FROM ONE"), paramInt8);
    auto call_out0_in1_customArgMsg = makeFunctionExt<0>(out0_and_in1, Utils::string("NOTHING_FROM ONE"), paramInt8);
    auto call_out0_in1_defaultArgMsg = makeFunctionExt<0>(out0_and_in1, paramInt8);

    //auto call_out1_in1_customArgMsg = makeFunctionExtArgMsg<1>(out1_and_in1, true, Utils::string("ONE_FROM_ONE"), paramInt8);
    auto call_out1_in1_customArgMsg = makeFunctionExt<1>(out1_and_in1, Utils::string("ONE_FROM_ONE"), paramInt8);
    auto call_out1_in1_defaultArgMsg = makeFunctionExt<1>(out1_and_in1, paramInt8);

    //auto call_out2_in0_customArgMsg = makeFunctionExtArgMsg<2>(out2_and_in0, true, Utils::string("TWO_FROM_NOTHING"));
    auto call_out2_in0_customArgMsg = makeFunctionExt<2>(out2_and_in0, Utils::string("TWO_FROM_NOTHING"));
    auto call_out2_in0_defaultArgMsg = makeFunctionExt<2>(out2_and_in0);

    //auto call_out0_in2_customArgMsg = makeFunctionExtArgMsg<0>(out0_and_in2, false, Utils::string("NOTHING_FROM_TWO"), paramInt8, paramInt8);
    auto call_out0_in2_customArgMsg = makeFunctionExt<0>(out0_and_in2, Utils::string("NOTHING_FROM_TWO"), paramInt8, paramInt8);
    auto call_out0_in2_defaultArgMsg = makeFunctionExt<0>(out0_and_in2,paramInt8, paramInt8);

    //auto call_out2_in1_customArgMsg = makeFunctionExtArgMsg<2>(out2_and_in1, true, Utils::string("TWO_FROM_ONE"), paramInt8);
    auto call_out2_in1_customArgMsg = makeFunctionExt<2>(out2_and_in1, Utils::string("TWO_FROM_ONE"), paramInt8);
    auto call_out2_in1_defaultArgMsg = makeFunctionExt<2>(out2_and_in1,paramInt8);

    //auto call_out1_in2_customArgMsg = makeFunctionExtArgMsg<1>(out1_and_in2, true, Utils::string("ONE_FROM_TWO"), paramInt8, paramInt8);
    auto call_out1_in2_customArgMsg = makeFunctionExt<1>(out1_and_in2, Utils::string("ONE_FROM_TWO"), paramInt8, paramInt8);
    auto call_out1_in2_defaultArgMsg = makeFunctionExt<1>(out1_and_in2, paramInt8, paramInt8);

    //auto call_out2_in2_customArgMsg = makeFunctionExtArgMsg<2>(out2_and_in2, true, Utils::string("TWO_FROM_TWO"), paramInt8, paramInt8);
    auto call_out2_in2_customArgMsg = makeFunctionExt<2>(out2_and_in2,  Utils::string("TWO_FROM_TWO"), paramInt8, paramInt8);
    auto call_out2_in2_defaultArgMsg = makeFunctionExt<2>(out2_and_in2, paramInt8, paramInt8);

    /*-----------------------------------------------------------------------------
     *  parse calls
     *-----------------------------------------------------------------------------*/

    auto parser = makeArgParser("return_ArgMsg",
        call_out0_in0_customArgMsg, call_out0_in0_defaultArgMsg,
        call_out1_in0_customArgMsg, call_out1_in0_defaultArgMsg,
        call_out0_in1_customArgMsg, call_out0_in1_defaultArgMsg,
        call_out1_in1_customArgMsg, call_out1_in1_defaultArgMsg,
        call_out2_in1_customArgMsg, call_out2_in1_defaultArgMsg,
        call_out1_in2_customArgMsg, call_out1_in2_defaultArgMsg,
        call_out2_in0_customArgMsg, call_out2_in0_defaultArgMsg,
        call_out0_in2_customArgMsg, call_out0_in2_defaultArgMsg,
        call_out2_in2_customArgMsg, call_out2_in2_defaultArgMsg
    );

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}


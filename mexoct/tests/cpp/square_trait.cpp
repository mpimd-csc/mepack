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

#include <string>

using namespace mexoct;


MEXOCT_ENTRY(square_trait,
"square_trait  -- Check if the square trait works\n"\
"\n"\
"    str = vector_trait(trait, v) \n")
{

    MEXOCT_INIT();
    std::string s_trait;
    auto paramSquare = Parameter<ArrayType<double>, Traits::SquareMatrix>("square vector");


    auto varSquare = makeFunction<1>(paramSquare);


    auto parser = makeArgParser("square_trait", varSquare);
    MEXOCT_PARSE(parser);

    if ( parser.isSelected(varSquare)){
        parser.setReturn(std::string("square"));
    } else {
        parser.setReturn(std::string("unknown"));

    }
    MEXOCT_RETURN;
}



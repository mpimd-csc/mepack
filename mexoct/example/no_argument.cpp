/**
 * EXAMPLE: std::sort
 * This example shows how a generic algorithm can be exposed to matlab
 * for relevant data types.
 */

#include "mexoct/interface.hpp"

using namespace mexoct;


MEXOCT_ENTRY(no_argument, "\
no_argument -- Example code with no argument\n\
\n")
{
    MEXOCT_INIT();

    auto fn = makeFunction<1>();

	auto parser = makeArgParser("no_argument",fn);
    MEXOCT_PARSE(parser);

    parser.setReturn(std::string("Hello World."));

    MEXOCT_RETURN;
}

/**
 * EXAMPLE: std::sort
 * This example shows how a generic algorithm can be exposed to matlab
 * for relevant data types.
 */

#include "mexoct/interface.hpp"

using namespace mexoct;

void print(ScalarType<double> v)
{
    std::cout << "double XXX " << v.value << "\n";
    mexoct_stdout << "double " <<  v.value << "\n";
    return;
}

std::tuple<ScalarType<double>> copy_2nd(ScalarType<double> a, ScalarType<double> b)
{
    return std::make_tuple(b);
}


MEXOCT_ENTRY(example_function_variants, "\
example_function_variants  Test example of MEXOCT\n\
 \n\
    [L,U] = LU(A) stores an upper triangular matrix in U and a\n\
    \"psychologically lower triangular matrix\" (i.e. a product of lower\n\
    triangular and permutation matrices) in L, so that A = L*U. A can be\n\
    rectangular.\n")
{
    MEXOCT_INIT();
	// descriptions of valid signatures
	// in this case each type that should be available needs to be instantiated
    auto ParamScalarDouble = Parameter<ScalarType<double>>("double");

    auto variant1 = makeFunctionExt<0>(print,   ParamScalarDouble);
    auto variant2 = makeFunctionExt<1>(copy_2nd,  ParamScalarDouble, ParamScalarDouble);

    auto parser = makeArgParser("example_function_variants", variant1, variant2);
//
    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
   	// check inputs and call function
}

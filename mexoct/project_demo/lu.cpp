
#include "mexoct/interface.hpp"
using namespace mexoct;


MEXOCT_ENTRY(lu, "\
lu\n\
 \n\
    [L,U] = LU(A) stores an upper triangular matrix in U and a\n\
    \"psychologically lower triangular matrix\" (i.e. a product of lower\n\
    triangular and permutation matrices) in L, so that A = L*U. A can be\n\
    rectangular.\n")
{
    MEXOCT_INIT();

    auto fn = makeFunction<1>();

	auto parser = makeArgParser("lu",fn);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
   	// check inputs and call function
}

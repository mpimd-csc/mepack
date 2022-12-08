/**
 * EXAMPLE: std::sort
 * This example shows how a generic algorithm can be exposed to matlab
 * for relevant data types.
 */

#include "mexoct/interface.hpp"

using namespace mexoct;

static void fail_function()
{
#ifdef MEXOCT_MATLAB
    mexErrMsgTxt("Error!");
#else
    error("Error!");
#endif


}

MEXOCT_ENTRY(allocate_memory, "\
allocate_memory  -- Example code to check the memory cleanup\n\
\n")
{
    MEXOCT_INIT();

    ScalarType<uint64_t> val;
    auto ParamScalarUInt64 = Parameter<ScalarType<uint64_t>>("size", &val);
    auto fn = makeFunction<0>(ParamScalarUInt64);

	auto parser = makeArgParser("allocate_memory",fn);
    try {
        MEXOCT_PARSE_NOEXCEPT(parser);
    } catch (MexOctException & e) {
        mexoct_error(e);
    }

    mexoct_stdout << "Allocate : " << val/(1024*1024) << " MB\n";

    // char * mem = new char[val];
    // std::shared_ptr<char> mem(new char[val], std::default_delete<char []>());
    MEXOCT_BUFFER(char, mem, val);
    fail_function();

    mexoct_stdout << "Hier\n";

    MEXOCT_RETURN;
   	// check inputs and call function
}

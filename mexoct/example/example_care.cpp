/**
 * EXAMPLE: care
 * This example shows basic usage of mex-interface and customization 
 * points for arbitrary use cases and a c-style backend.
 */ 

#include "paramtypes_care.hpp"
#include "argparser.hpp"

enum struct Method : int
{
	Foo,
	Bar
};

// Implementations for different variants that share the same name in
// matlab.
Matrix care(Matrix A, Matrix B, Matrix C);
Matrix care_sparse(Matrix A, Matrix B, Matrix C, Method method);
Matrix gcare(Matrix A, Matrix B, Matrix C, Matrix E);

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	// destination for the arguments
	Matrix A,B,C,E;
	Method method; 

	// Define possible parameters.
	// The first template argument determines the expected type.
	// Parameter<Matrix> needs to be implemented to enable Matrix in this
	// situation.
	// The name should be consistent with the public documentation.
	// By giving pointers to Matrix-objects the arguments are directly
	// written on success.
	auto paramA = Parameter<Matrix, Traits::Dense>("A", &A);
	auto paramB = Parameter<Matrix, Traits::Dense>("B", &B);
	auto paramC = Parameter<Matrix, Traits::Dense>("C", &C);
	auto paramE = Parameter<Matrix, Traits::Dense>("E", &E);

	// The Matrix type handles both sparse and dense matrices so an 
	// additional trait is used to differentiate.
	auto paramAs = Parameter<Matrix, Traits::Sparse>("A", &A);
	auto paramBs = Parameter<Matrix, Traits::Sparse>("B", &B);
	auto paramCs = Parameter<Matrix, Traits::Sparse>("C", &C);

	// The destination pointer is optional
	auto paramMethod = Parameter<Utils::string>("method");

	// Describe function signatures from the created parameters.
	// All variants have one return value.
	auto variant1 = makeFunction<1>(paramA, paramB, paramC);
	// They may share some of the parameters.
	auto variant2 = makeFunction<1>(paramA, paramB, paramC, paramE);
	// The makeFunctionExt variant takes either a function pointer or 
	// lambda as first argument. This can be used to either automaticaly
	// call the implementation on match or to handle some special argument.
	// As first argument this function is given the return value array
	// from mexFunction.
	// The rest of the signature is defined by the parameter descriptions.
	auto variant3 = makeFunctionExt<1>([&method]
	(const Matrix& A, const Matrix& B, const Matrix& C, 
		const Utils::string& methodStr)
	{
		if(methodStr == "foo")
			method = Method::Foo;
		else 
			method = Method::Bar; 
	},paramAs, paramBs, paramCs, paramMethod);

	// Create an ArgParser from all possible functions.
	auto parser = makeArgParser("example_care", variant1, variant2, variant3);
	// Check inputs, assign arguments and call the function.
	// Leave the function if the arguments are invalid.
	if(!parser.parse(nlhs, plhs, nrhs, prhs))
		return;

	Matrix result{};

	// After parsing determine which variant was choosen and call the 
	// appropriate implementation.
	if(parser.isSelected(variant1))
		result = care(A,B,C);
	else if(parser.isSelected(variant2))
		result = gcare(A,B,C,E);
	else
		result = care_sparse(A,B,C, method);

	// Convert value back and set the return array.
	parser.setReturn(result);
}

Matrix care(Matrix A, Matrix B, Matrix C)
{
	mexPrintf("called care\n");
	return A;
}

Matrix care_sparse(Matrix A, Matrix B, Matrix C, Method method)
{
	mexPrintf("called care_sparse  with alg=%d\n", static_cast<int>(method));
	return A;
}

Matrix gcare(Matrix A, Matrix B, Matrix C, Matrix E)
{
	mexPrintf("called gcare\n");
	return A;
}
#pragma once

#include "paramtypes.hpp"
#include "utils.hpp"

// pod dynamic size matrix
struct Matrix
{
    int rows, columns;
    double* elements;
    bool isSparse;
};

// Specialize ParamType<> for Matrix so that this type can be used in a 
// mex interface. 
template<>
struct ParamType<Matrix>
{
	static bool isValid(const mxArray* _array) 
    {
        // Matrix only supports doubles 
        return !mxIsEmpty(_array)  
			&& mxIsDouble(_array)
            && !mxIsComplex(_array);
    }

    // This function creates a Matrix object from an mxArray
    static Matrix fromMx(const mxArray* _array)
    {
        return Matrix{};
    }

	static mxArray* toMx(const Matrix& _matrix)
	{
		return mxCreateDoubleMatrix(_matrix.rows, _matrix.columns, mxREAL);
	}

    // This name is displayed when a type missmatch is detected.
    // It should make sense for the matlab user!
	static constexpr const char* NAME = "real matrix";
};

template<>
struct ParamType<Utils::string>
{
	static bool isValid(const mxArray* _array) 
    {
        return mxIsChar(_array) && mxGetM(_array);
    }

    static Utils::string fromMx(const mxArray* _array)
    {
        char* mxs = mxArrayToString(_array);
        Utils::string s(mxs ? mxs : "");
        // the string creates a copy of the buffer
        mxFree(mxs);
        return s;
    }

	static constexpr const char* NAME = "string";
};


// the namespace is not necessary but allows for shorter names without
// collisions
namespace Traits{

// A trait is an additional property that a argument can be required 
// to possess and that is not determined by the base type.
struct Dense
{
    // verify property off the original input
    // returns true when _array is dense
	static bool isValid(const mxArray* _array)
	{
		return !mxIsSparse(_array);
	}

    // name as displayed in a error message
	static constexpr const char* NAME = "dense";
};

struct Sparse
{
	static bool isValid(const mxArray* _array)
	{
		return mxIsSparse(_array);
	}
	static constexpr const char* NAME = "sparse";
};

}
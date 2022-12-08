/**
 * EXAMPLE: std::sort
 * This example shows how a generic algorithm can be exposed to matlab
 * for relevant data types.
 */

#include "mexoct/interface.hpp"

using namespace mexoct;

#if 0
// wrapper around the actual sort routine to handle the return value
template<typename T>
std::tuple<FixedArray<T>> sort(FixedArray<T> _array)
{
	FixedArray<T> result(_array);

	std::sort(result.elements, result.elements + result.size);

	return std::make_tuple(result);
}

template<typename T>
std::tuple<FixedArray<T>, FixedArray<T>> sort2(FixedArray<T> _array1, FixedArray<T> _array2)
{
	FixedArray<T> result(_array1);
	FixedArray<T> result2(_array2);

	std::sort(result.elements, result.elements + result.size);
	std::sort(result2.elements, result2.elements + result2.size);

	// the tuple allows to return multiple values (of different types)
	return std::make_tuple(result, result2);
}

#endif

class TestClass {
    private:
        int m_x;
    public:

        TestClass (int x ){
            m_x = x;
        }

        ~TestClass() {
            mexoct_stdout << "That's it \n";
        }
};

template<typename T>
std::tuple<> print(T s)
{
    // error("Unkown type.\n");
    return std::make_tuple();

}

template<>
std::tuple<> print(std::string s) {
    mexoct_stdout << "string " << s << "\n";
    return std::make_tuple();
}


template<>
std::tuple<> print(ScalarType<int8_t> v){
    mexoct_stdout << "int8_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<uint8_t> v){
    mexoct_stdout << "uint8_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}
template<>
std::tuple<> print(ScalarType<int16_t> v){
    mexoct_stdout << "int16_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<uint16_t> v){
    mexoct_stdout << "uint16_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}
template<>
std::tuple<> print(ScalarType<int32_t> v){
    mexoct_stdout << "int32_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<uint32_t> v){
    mexoct_stdout << "uint32_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}
template<>
std::tuple<> print(ScalarType<int64_t> v){
    mexoct_stdout << "int64_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<uint64_t> v){
    mexoct_stdout << "uint64_t " << static_cast<int>( v.value) << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<double> v){
    mexoct_stdout << "double " <<  v.value << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<float> v){
    mexoct_stdout << "float " <<  v.value << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<std::complex<double>> v){
    mexoct_stdout << "double " <<  v.value << "\n";
    return std::make_tuple();
}

template<>
std::tuple<> print(ScalarType<std::complex<float>> v){
    mexoct_stdout << "float " <<  v.value << "\n";
    return std::make_tuple();
}



template<typename T>
std::tuple<ScalarType<T>> copy_val ( ScalarType<T> v ) {
    mexoct_stdout << "copy" << "\n";
    T val = v.value;
    return std::make_tuple(ScalarType<T>(val));
}

template<typename T>
std::tuple<> printa(ArrayType<T> a) {
    mexoct_stdout << "Out array:" << "\n";
    mexoct_stdout << a.rows << " " << a.columns << " " << a.ld << "\n";
    mexoct_stdout << "Ptr: " << a.values << "\n";
    ssize_t i, j;

    for (i = 0; i < a.rows; i++) {
        for (j = 0; j < a.columns; j++) {
            mexoct_stdout << a[i+j*a.rows] << " ";
        }
        mexoct_stdout << "\n";
    }
    mexoct_stdout << "Count :" << a.values.use_count() << "\n";

    return std::make_tuple();
}

template<typename T>
auto copy_array(ArrayType<T> a) -> std::tuple<ArrayType<T>> {
    return std::make_tuple(a);
}


MEXOCT_ENTRY(example_info, "\
example_info  Test example of MEXOCT\n\
 \n\
    [L,U] = LU(A) stores an upper triangular matrix in U and a\n\
    \"psychologically lower triangular matrix\" (i.e. a product of lower\n\
    triangular and permutation matrices) in L, so that A = L*U. A can be\n\
    rectangular.\n")
{
    MEXOCT_INIT();
	// descriptions of valid signatures
	// in this case each type that should be available needs to be instantiated
    auto ParamScalarInt8 = Parameter<ScalarType<int8_t>>("i8", "signed integer 8 value");
    auto ParamScalarUInt8 = Parameter<ScalarType<uint8_t>>("u8", "bloeder unsigned 8 wert");
    auto ParamScalarInt16 = Parameter<ScalarType<int16_t>>("int16_t");
    auto ParamScalarUInt16 = Parameter<ScalarType<uint16_t>>("uint16_t");
    auto ParamScalarInt32 = Parameter<ScalarType<int32_t>>("int32_t");
    auto ParamScalarUInt32 = Parameter<ScalarType<uint32_t>>("uint32_t");
    auto ParamScalarInt64 = Parameter<ScalarType<int64_t>>("int64_t");
    auto ParamScalarUInt64 = Parameter<ScalarType<uint64_t>>("uint64_t");
    auto ParamScalarFloat  = Parameter<ScalarType<float>>("float");
    auto ParamScalarDouble = Parameter<ScalarType<double>>("double");
    auto ParamScalarFloatCpx  = Parameter<ScalarType<std::complex<float>>>("float complex");
    auto ParamScalarDoubleCpx = Parameter<ScalarType<std::complex<double>>>("double complex");
    auto ParamString          = Parameter<std::string>("string");
    auto ParamMatrix          = Parameter<ArrayType<double>>("double array");
    auto ParamComplexMatrix          = Parameter<ArrayType<std::complex<double>>>("double complex array");


    auto variant1 = makeFunctionExt<0>(print<ScalarType<int8_t>>,   ParamScalarInt8);
    auto variant2 = makeFunctionExt<0>(print<ScalarType<uint8_t>>,  ParamScalarUInt8);
    auto variant3 = makeFunctionExt<0>(print<ScalarType<int16_t>>,  ParamScalarInt16);
    auto variant4 = makeFunctionExt<0>(print<ScalarType<uint16_t>>, ParamScalarUInt16);
    auto variant5 = makeFunctionExt<0>(print<ScalarType<int32_t>>,  ParamScalarInt32);
    auto variant6 = makeFunctionExt<0>(print<ScalarType<uint32_t>>, ParamScalarUInt32);
    auto variant7 = makeFunctionExt<0>(print<ScalarType<int64_t>>,  ParamScalarInt64);
    auto variant8 = makeFunctionExt<0>(print<ScalarType<uint64_t>>, ParamScalarUInt64);
    auto variant9 = makeFunctionExt<0>(print<ScalarType<float>>,  ParamScalarFloat);
    auto variant10 = makeFunctionExt<0>(print<ScalarType<double>>, ParamScalarDouble);
    auto variant11 = makeFunctionExt<0>(print<ScalarType<std::complex<float>>>,  ParamScalarFloatCpx);
    auto variant12 = makeFunctionExt<0>(print<ScalarType<std::complex<double>>>, ParamScalarDoubleCpx);

    auto variant13 = makeFunctionExt<1>(copy_val<int8_t>, ParamScalarInt8);

    auto variant14 = makeFunctionExt<0>(print<std::string>, ParamString);
    auto variant15 = makeFunctionExt<0>(printa<double>, ParamMatrix);
    auto variant15a = makeFunctionExt<0>(printa<std::complex<double>>, ParamComplexMatrix);

    auto variant16 = makeFunctionExt<1>(copy_array<double>, ParamMatrix);
    auto variant17 = makeFunctionExt<1>(copy_array<std::complex<double>>, ParamComplexMatrix);
    auto variant18 = makeFunction<0>(ParamScalarUInt8, ParamScalarInt8);

	auto parser = makeArgParser("example_info", variant1, variant2, variant3, variant4, variant5, variant6,
                                                variant7, variant8, variant9, variant10, variant11, variant12,
                                                variant13, variant14, variant15, variant15a, variant16,variant17, variant18);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
   	// check inputs and call function
}

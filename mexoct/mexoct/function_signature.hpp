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

#pragma once

#ifdef MEXOCT_OCTAVE
#define MEXOCT_CONST_VALUE const octave_value &
#define MEXOCT_VALUE octave_value &
#define MEXOCT_CONST_VALUE_LIST const octave_value_list &
#define MEXOCT_VALUE_LIST octave_value_list &

#else
#define MEXOCT_CONST_VALUE const mxArray *
#define MEXOCT_VALUE mxArray *
#define MEXOCT_CONST_VALUE_LIST const mxArray **
#define MEXOCT_VALUE_LIST mxArray **

#endif

namespace mexoct {

    /**
     * Specifies the signature of a function callable from matlab.
     * Use makeFunction or makeFunctionExt to create an object of this type!
     */
    template<int Results, typename TargetFn, typename... Params>
        class FunctionSignature
        {
            public:
                static constexpr int NumResults = Results;
                static constexpr int NumArgs = sizeof...(Params);

                FunctionSignature(TargetFn _targetFn, Utils::string return_ArgMsg,  const Params&... _params)
                    : m_errorType(ErrorType::None),
                    m_errorBuildFn(nullptr),
                    m_targetFunction(std::move(_targetFn)),
                    m_parameters(_params...),
                    m_return_ArgMsg(return_ArgMsg)  {}

                // Checks whether the given arguments satisfy this signature.
                // Tested are the number of arguments, the number of return values
                // and the conditions for the individual parameters.
                // @return true for a valid match
                bool check(int nlhs, MEXOCT_VALUE_LIST plhs, int nrhs, MEXOCT_CONST_VALUE_LIST prhs)
                {
                    if(nrhs != NumArgs)
                    {
                        m_errorType = ErrorType::ArgNum;
                        return false;
                    }

                    /* Required to allow the typical as quick & dirty call of functions */
                    if ( nlhs == 0 && NumResults <= 1)
                    {
                        return checkTypes<NumArgs, Params...>(0, prhs);
                    }
                    if ( nlhs == 0 && NumResults >= 2) {
                        m_errorType = ErrorType::ResultNum;
                        return false;
                    }

                    if( nlhs > 0 && nlhs != NumResults)
                    {
                        m_errorType = ErrorType::ResultNum;
                        return false;
                    }
                    return checkTypes<NumArgs, Params...>(0, prhs);
                }

                // Applies the ParamType conversions and calls the linked function.
                // Does not verify the arguments.
                void operator()(int nlhs, MEXOCT_VALUE_LIST plhs, int nrhs, MEXOCT_CONST_VALUE_LIST prhs ) const
                {
                    call(plhs, prhs, Details::make_index_sequence<NumArgs>());
                }

                // Create a string describing this Signature.
                Utils::string toString() const
                {
                    return "(" + printImpl<0>() + ")";
                }

                Utils::string printArgument(int i) const
                {
                    return printArgumentImpl<0>(i);
                }

                Utils::string printReturn() const
                {
                    if (!m_return_ArgMsg.empty()){
                        return m_return_ArgMsg + " = ";
                    }

                    if ( NumResults == 0 ) {
                        return Utils::string("");
                    }
                    if ( NumResults == 1 ) {
                        return Utils::string("ret1 = ");
                    }
                    Utils::string ret("[");
                    for ( int i = 0; i < NumResults; i++) {
                        if ( i > 0 ) ret += ", ";
                        ret += "ret" + std::to_string(i+1);
                    }
                    ret +="] = ";
                    return ret;
                }


                // Get the message of the last error that occurred.
                const char* getErrorMsg() const
                {
                    switch (m_errorType)
                    {
                        case ErrorType::None:
                            return "Unknown error occurred.";
                        case ErrorType::ArgNum:
                            return "Incorrect number of arguments.";
                        case ErrorType::ResultNum:
                            return "Incorrect number of return values.";
                        case ErrorType::ArgType:
                            m_errorMsg = m_errorBuildFn(*this);
                            break;
                    };
                    return m_errorMsg.c_str();
                }
            private:
                template<int ArgCount, typename Param, typename... ParamsLeft>
                    bool checkTypes(int pos, MEXOCT_CONST_VALUE_LIST prhs)
                    {
                        #ifdef MEXOCT_MATLAB
                        if(!Param::isValid(prhs[pos]))
                        #else
                        if(!Param::isValid(prhs(pos)))
                        #endif
                        {
                            m_errorType = ErrorType::ArgType;
                            m_errorBuildFn = buildTypeErrorMessage<Param, NumArgs - ArgCount>;
                            return false;
                        }

                        return checkTypes<ArgCount-1, ParamsLeft...>(pos+1, prhs);
                    }

                // recursion stop
                template<int ArgCount>
                    bool checkTypes(int pos, MEXOCT_CONST_VALUE_LIST prhs) const
                    {
                        return true;
                    }

                /* Call the function if m_targetFunction returns void (Using Tag dispatch) */
                template<typename... TS, int... Indices>
                    void callImpl(std::true_type, MEXOCT_VALUE_LIST plhs, std::tuple<TS...> const & tup, Details::index_sequence<Indices...>) const {
                        // std::cout << "Function returns void. \n ";
                        m_targetFunction(assignParameter<Indices>(std::get<Indices>(tup))...);
                        return;
                    }

                /* Call the function of m_targetFunction does not return void.  */
                template<typename... TS, int... Indices>
                    void callImpl(std::false_type, MEXOCT_VALUE_LIST plhs, std::tuple<TS...> const & tup, Details::index_sequence<Indices...>) const {
                        // std::cout << "Function returns something\n";
                        auto results = m_targetFunction(assignParameter<Indices>(std::get<Indices>(tup))...);
                        // static_assert(std::tuple_size<decltype(results)>::value == NumResults, "Incorrect number of return values provided.");
                        if ( std::tuple_size<decltype(results)>::value >= NumResults) {
                            expandResults(plhs, results, Details::make_index_sequence<NumResults>());
                        }
                        return;
                    }

                template<int... Indices>
                    void call(MEXOCT_VALUE_LIST plhs, MEXOCT_CONST_VALUE_LIST prhs, Details::index_sequence<Indices...> num_args_seq) const
                    {
                        /*  Setup the argument list  */
                        #ifdef MEXOCT_MATLAB
                        auto args = std::make_tuple(Params::TypeInfo::copy_in(prhs[Indices])...);
                        #else
                        auto args = std::make_tuple(Params::TypeInfo::copy_in(prhs(Indices))...);
                        #endif

                        /* Determine the return type  */
                        using return_type = decltype(m_targetFunction(assignParameter<Indices>(std::get<Indices>(args))...));

                        /* Call the function depending on the return type.  */
                        callImpl(std::integral_constant<bool,std::is_same<return_type, void>::value>{}, plhs, args,  num_args_seq  );

                    }

                // Assigns the given argument to the target referenced in the
                // associated Parameter<>.
                // Returns a reference to either the original value or the Parameter
                // target if available.
                template<int Ind, typename T>
                    T& assignParameter(T& value) const
                    {
                        auto& param = std::get<Ind>(m_parameters);
                        // constness is not transitive and param::target is conceptionally
                        // not part of the FunctionSignature
                        if(param.target)
                        {
                            // original value is not needed anymore
                            *param.target = std::move(value);
                            return *param.target;
                        }
                        return value;
                    }

                template<typename TT, int... Indices>
                    void expandResults(MEXOCT_VALUE_LIST plhs, const TT& _tuple, Details::index_sequence<Indices...>) const
                    {
                        #ifdef MEXOCT_MATLAB
                        std::make_tuple(plhs[Indices] = ParamType<typename std::tuple_element<Indices, TT>::type>::copy_out(std::get<Indices>(_tuple))...);
                        #else
                        std::make_tuple(plhs(Indices) = ParamType<typename std::tuple_element<Indices, TT>::type>::copy_out(std::get<Indices>(_tuple))...);
                        #endif
                    }

                /* The function has more than one argument */
                template<int Ind, typename std::enable_if<Ind < NumArgs, int>::type = 0>
                    Utils::string printImpl() const
                    {
                        const auto& param = std::get<Ind>(m_parameters);
                        if ( Ind > 0 ) {
                            return Utils::string(", ") + param.name + printImpl<Ind+1>();
                        } else {
                            return Utils::string(param.name) + printImpl<Ind+1>();
                        }
                    }

                /* Recursion stop, the function has one argument */
                template<int Ind, typename std::enable_if<Ind >= NumArgs, int>::type = 0>
                    Utils::string printImpl() const
                    {
                        return Utils::string("");
                    }


                /*  The function has no arguments */
                template<int Ind, typename std::enable_if< Ind >= NumArgs, int>::type = 0>
                    Utils::string printArgumentImpl(int i) const
                    {
                        return Utils::string("");
                    }


                /* The function has more than one argument */
                template<int Ind, typename std::enable_if<Ind < NumArgs, int>::type = 0>
                    Utils::string printArgumentImpl(int i) const
                    {
                        if ( Ind == i ) {
                            const auto& param = std::get<Ind>(m_parameters);
                            if( strlen(param.desc) > 0 ) {
                                return Utils::string(param.name) + " -- " + param.desc + " (expected type: " + param.getTypeName() + ")";
                            } else {
                                return Utils::string(param.name) + " (expected type: " + param.getTypeName() + ")";
                            }
                        } else {

                        // return param.name + Utils::string(", ") + printImpl<Ind+1>();
                            return printArgumentImpl<Ind+1>(i);
                        }
                    }

                template<typename Param, size_t ParamInd>
                    static Utils::string buildTypeErrorMessage(const FunctionSignature& _this)
                    {
                        const auto& param = std::get<ParamInd>(_this.m_parameters);
                        return Utils::string("Argument ")
                            + param.name
                            + " is of the wrong type, expected: "
                            + Param::getTraitNames() + " "
                            + Param::TypeInfo::NAME
                            + ".";
                    }

                enum struct ErrorType{
                    None,
                    ArgNum,
                    ResultNum,
                    ArgType
                };
                ErrorType m_errorType;
                using ErrorBuildFn = Utils::string(*)(const FunctionSignature&);
                ErrorBuildFn m_errorBuildFn;
                mutable Utils::string m_errorMsg;

                TargetFn m_targetFunction;
                std::tuple<Params...> m_parameters;
                Utils::string m_return_ArgMsg;
        };

    /**
     * Make a FunctionSignature with template deduction pre C++17.
     * @param NumResults - The number of return values this function has.
     * @param _fn - A function with the signature void(Params::Type...)
     *   that is called when the signature is matched. The first argument is the
     *   array of return values for mexFunction.
     *   Both a lambda or a regular pointer are valid.
     *   If USE_EXT_RETURN is defined, _fn should return an std::tuple<...>
     *   instead, which contents are converted and returned to mex.
     * @param _params - Any number of Parameter<...> objects that describe the
     *   arguments this function variant takes.
     * @return The complete FunctionSignature object.
     */
    // Variant with custom return Argument Message
    template<int NumResults, typename... Params, typename TargetFn>
        auto makeFunctionExt(TargetFn _fn, Utils::string return_ArgMsg, const Params&... _params)
        -> FunctionSignature<NumResults, TargetFn, Params...>
        {
            using namespace Details;

            static_assert(all_true<is_parameter_spec<Params>::value...>::value,
                    "All params are expected to be specializations of Parameter<...>.");

            return FunctionSignature<NumResults, TargetFn, Params...>(_fn,return_ArgMsg, _params...);
        }
    template<int NumResults, typename... Params, typename TargetFn>
        auto makeFunctionExt(TargetFn _fn, const char *return_ArgMsg, const Params&... _params)
        -> FunctionSignature<NumResults, TargetFn, Params...>
        {
            using namespace Details;

            static_assert(all_true<is_parameter_spec<Params>::value...>::value,
                    "All params are expected to be specializations of Parameter<...>.");

            return FunctionSignature<NumResults, TargetFn, Params...>(_fn,Utils::string(return_ArgMsg), _params...);
        }


    // Variant with default return Argument Message
    template<int NumResults, typename... Params, typename TargetFn>
        auto makeFunctionExt(TargetFn _fn, const Params&... _params)
        -> FunctionSignature<NumResults, TargetFn, Params...>
        {
            return makeFunctionExt<NumResults, Params...>(_fn, Utils::string(""), _params...);
        }

    // Simpler variant if no custom function is needed.
    // See makeFunctionExt for more information.
    template<int NumResults, typename... Params>
        auto makeFunction(const Params&... _params)
        -> FunctionSignature<NumResults, void(*)(typename Params::Type...), Params...>
        {
            // lambda that does nothing but has the right signature
            auto emptyFn = [](typename Params::Type...){};
            // cast to pointer
            using FnPtr = void(*)(typename Params::Type...);
            FnPtr fnPtr = emptyFn;
            return makeFunctionExt<NumResults, Params...>(fnPtr, Utils::string(""), _params...);
        }

}

#undef MEXOCT_CONST_VALUE
#undef MEXOCT_VALUE
#undef MEXOCT_CONST_VALUE_LIST
#undef MEXOCT_VALUE_LIST

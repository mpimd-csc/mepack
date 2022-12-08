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

#include <utility>
#include <string>
#include <tuple>
#include <cassert>

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
     * Handles arguments given by matlab.
     * Use makeArgParser() to create objects of this type!
     */
    template<typename... FnSignatures>
        class ArgParser
        {
            static constexpr size_t NumSignatures = sizeof...(FnSignatures);
            public:
            ArgParser(const char* _name, FnSignatures... _functions)
                : m_functions(_functions...),
                m_name(_name),
                m_plhs()
            {}

            // Checks the given arguments against known signatures.
            // If a match was found argument conversion takes place and
            // the associated function is invoked.
            // @params The arguments as given to mexFunction().
            // @return true if a matching signature was found.
#ifdef MEXOCT_MATLAB
            bool parse(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
#else
            bool parse(int nlhs, octave_value_list ** plhs, int nrhs, const octave_value_list & prhs)
#endif
            {

#ifdef MEXOCT_MATLAB
                if (!iterateChecks<0, FnSignatures...>(nlhs, plhs, nrhs, prhs))
#else
                *plhs = &(this->m_plhs);
                if (!iterateChecks<0, FnSignatures...>(nlhs, m_plhs, nrhs, prhs))
#endif
                {
                    Utils::string errorMsg;
                    if(NumSignatures > 1)
                        errorMsg = "No overload matches the given argument or return value combination:\n";
                    printFunctions<0>(errorMsg);

                    ParserException err(errorMsg);
                    throw err;

                    return false;
                }
#ifdef MEXOCT_MATLAB
                m_plhs = plhs;
#endif

                return true;
            }

            // Checks whether the given FunctionSignature has been choosen.
            // When parse() has not been executed or false was returned by it,
            // this function always returns false.
            template<typename T>
                bool isSelected(const T&) const
                {
                    return Details::type_index<0,T, std::tuple<FnSignatures...>>::index == m_choice;
                }

            template<typename... Args>
                void setReturn(const Args&... _args)
                {
                    // assert(sizeof...(_args) == m_numResults);
                    setReturnImpl<0>(_args...);
                }

            private:
            template<int Count, typename Fn, typename... Fns>
                bool iterateChecks(int nlhs, MEXOCT_VALUE_LIST plhs, int nrhs, MEXOCT_CONST_VALUE_LIST prhs)
                {
                    auto& fn = std::get<Count>(m_functions);

                    if(fn.check(nlhs, plhs, nrhs, prhs))
                    {
                        // save choice to allow for quick isSelected() checks after
                        m_choice = Count;
                        m_numResults = fn.NumResults;
                        // call matching function
                        fn(nlhs,plhs,nrhs,prhs);
                        return true;
                    }

                    return iterateChecks<Count+1, Fns...>(nlhs, plhs, nrhs, prhs);
                }

            // recursion end, no match found
            template<int Count>
                bool iterateChecks(int nlhs, MEXOCT_VALUE_LIST plhs, int nrhs, MEXOCT_CONST_VALUE_LIST prhs)
                {
                    return false;
                }


            template<int Ind, typename std::enable_if<Ind < NumSignatures, int>::type = 0>
                void printFunctions(Utils::string& _errorMessage) const
                {
                    const auto& fn = std::get<Ind>(m_functions);
                    if (Ind > 0) _errorMessage +="\n";
                    _errorMessage += fn.printReturn() + m_name + fn.toString() + "\n";
                    for (int i = 0; i < fn.NumArgs; i++) {
                        _errorMessage += Utils::string("    ") + fn.printArgument(i) +"\n";
                    }
                    _errorMessage += Utils::string("    -- Error -- ") +fn.getErrorMsg() + "\n";
                    printFunctions<Ind+1>(_errorMessage);
                }

            template<int Ind, typename std::enable_if<Ind >= NumSignatures, int>::type = 0>
                void printFunctions(Utils::string&) const
                {
                }

            template<int Ind, typename Arg, typename... Args>
                void setReturnImpl(const Arg& _arg, const Args&... _args)
                {
                    if ( Ind >= m_numResults ) return;
#ifdef MEXOCT_MATLAB
                    m_plhs[Ind] = ParamType<Arg>::copy_out(_arg);
#else
                    m_plhs(Ind) = ParamType<Arg>::copy_out(_arg);
#endif

                    setReturnImpl<Ind+1>(_args...);
                }

            template<int Ind>
                void setReturnImpl()
                {
                }

            std::tuple<FnSignatures...> m_functions;
            const char* m_name;
            int m_choice = -1;
            int m_numResults = 0;
#ifdef MEXOCT_MATLAB
            mxArray** m_plhs;
#else
            octave_value_list m_plhs;
#endif
        };

    /**
     * \brief Create an argument parser based on a set of possible function calls.
     * \param _name - Name of the function to be displayed on error.
     * \param _functions - Any amount of FunctionSignature<...> objects that are considered when parsing the input.
     *
     * The makeArgParser function creates and \ref mexoct::ArgParser object from a given set of possible function
     * signatures. The _name argument should be provide the same name as the function is named from a MATLAB or an
     * Octave point of view. The name is required to provide a meaningful error message if no function could
     * match the acutal set of provide arguments. The function signatures are create using
     * \ref mexoct::makeFunction and \ref mexoct::makeFunctionEx functions.
     *
     */
    template<typename... Functions>
        auto makeArgParser(const char* _name, Functions... _functions)
        -> ArgParser<Functions...>
        {
            return ArgParser<Functions...>(_name, _functions...);
        }

}

#undef MEXOCT_CONST_VALUE
#undef MEXOCT_VALUE
#undef MEXOCT_CONST_VALUE_LIST
#undef MEXOCT_VALUE_LIST

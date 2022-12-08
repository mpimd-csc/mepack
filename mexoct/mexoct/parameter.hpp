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
#else
#define MEXOCT_CONST_VALUE const mxArray *
#endif

namespace mexoct {
    /**
     * Holds information on how a single argument should be handled.
     * @param T - The type that should be converted to. It is required that
     *   T is copy assignable and a specialization of ParamType is available.
     * @param Traits - Any number of types that implement a isValid function
     *   and have a static member NAME. For details see paramtypes.hpp.
     */
    template<typename T, typename... Traits>
        struct Parameter
        {
            using TypeInfo = ParamType<T>;
            using Type = T;
            static constexpr int NumTraits = sizeof...(Traits);

            // _name used in error messages
            // _target variable in which the result is written into
            Parameter(const char* _name, T* _target = nullptr)
                : name(_name), desc("") ,target(_target)
            {}
            Parameter(const char* _name, const char *_desc, T* _target = nullptr)
                : name(_name), desc(_desc) ,target(_target)
            {}

            // Verifies that _array can be converted to type T and satisfies
            // all traits.
            static bool isValid(MEXOCT_CONST_VALUE _array)
            {
                return TypeInfo::isValid(_array) && isValid<0, Traits...>(_array);
            }

            // Returns a whitespace seperated list of trait names.
            static Utils::string getTraitNames()
            {
                Utils::string str;
                appendTraitNames<0, Traits...>(str);

                return str;
            }

            Utils::string getTypeName() const {
                return Utils::string(TypeInfo::NAME);
            }

            const char* name;
            const char* desc;
            T* target;

            private:
            template<int, typename Trait, typename... TraitsLeft>
                static bool isValid(MEXOCT_CONST_VALUE _array)
                {
                    return Trait::isValid(_array) && isValid<0,TraitsLeft...>(_array);
                }

            template<int>
                static bool isValid(MEXOCT_CONST_VALUE _array)
                {
                    return true;
                }

            template<int, typename Trait, typename... TraitsLeft>
                static void appendTraitNames(Utils::string& _str)
                {
                    _str += Trait::NAME;
                    appendTraitNames<0,TraitsLeft...>(_str);
                }

            template<int>
                static void appendTraitNames(Utils::string& _str)
                {
                }
        };

    namespace Details{
        template<typename T>
            struct is_parameter_spec : std::false_type {};

        template<typename T, typename... Traits>
            struct is_parameter_spec<Parameter<T,Traits...>> : std::true_type {};
    }



}
#undef MEXOCT_CONST_VALUE

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

/**
 * Defines the requirements and conversion operation for a parameter type
 * on the mex side.
 */
namespace mexoct {


    template<typename T>
        struct ParamType
        {
            static_assert(sizeof(T) == -1, "A type is missing a specialization for ParamType.");

#ifdef MEXOCT_OCTAVE
            static bool isValid(const octave_value & _array);
            static T copy_in(const octave_value & _array);
            static octave_value copy_out(const T& _array);


#else
            // @return true when the given mxArray can be converted to T.
            static bool isValid(const mxArray* _array);

            // Create a T object from the given mxArray.
            static T copy_in(const mxArray* _array);
            // Convert a T object to an mxArray.
            static mxArray* copy_out(const T& _object);
#endif
            // Name that is displayed as part of the error message
            // when a type mismatch is registered.
            static constexpr const char* NAME = "undefined";
        };

    /**
     * A Trait defines an additional property that a type should possess.
     */
    struct Trait
    {

#ifdef MEXOCT_OCTAVE
        static bool isValid(const octave_value & _array);
#else
        //@return true when the given mxArray satisfies this trait.
        static bool isValid(const mxArray* _array);
#endif
        // Name that is displayed in an error message.
        static constexpr const char* NAME = "undefined";
    };

    class Validator {
        public:
#ifdef MEXOCT_OCTAVE
        virtual bool isValid(const octave_value & _array) const  = 0;
#else
        //@return true when the given mxArray satisfies this trait.
        virtual bool isValid(const mxArray* _array) const = 0;
#endif

    };
}

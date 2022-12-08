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

#include <iostream>
#include <ostream>

namespace mexoct {

    class OutputStream : public std::streambuf
    {
        protected:
            virtual int_type overflow(int_type ch) override
            {
#ifdef MEXOCT_MATLAB
                mexPrintf("%c", ch);
#else
                char c = static_cast<char>(ch);
                octave_stdout << c;
#endif
                return 1;
            }
    };

    class Output {
        private:
            OutputStream os;
            std::ostream *ox;

            Output() : os()  {
                ox = new std::ostream(&os);
            }

            ~Output() {
                delete ox;
            }
        public:
            Output(const Output &) = delete;
            Output& operator=(const Output &) = delete;

            static Output &instance () {
                static Output o;
                return o;
            }
            static std::ostream & stream() {
                return *instance().ox;
            }
    };
}

#define mexoct_stdout (mexoct::Output::stream())


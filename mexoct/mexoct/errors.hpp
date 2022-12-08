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

#include <exception>

namespace mexoct {
    class MexOctException : public std::exception
    {
        private:
            static constexpr const char * m_id ="MEXOCT:ERROR";
        public:

            const char *id() const noexcept
            {
                return m_id;
            }

    };
    class ParserException : public MexOctException
    {
        private:
            Utils::string m_msg;
            static constexpr const char * m_id = "MEXOCT:ParserException";


        public:
            ParserException(const char * _msg) : m_msg(_msg) {}
            ParserException(const Utils::string & _msg) : m_msg(_msg) {}
            ParserException() : m_msg("Unknown Exception") {};


            const char* what() const noexcept {
                return m_msg.c_str();
            }

            const char *id() const noexcept {
                return m_id;
            }
    };

    void mexoct_error(const char * s) {
        #ifdef MEXOCT_MATLAB
            mexErrMsgTxt(s);
        #else
            error("%s",s);
        #endif
    }

    void mexoct_error(const Utils::string & s) {
        mexoct_error(s.c_str());
    }

    void mexoct_error(const char *id, const char * s) {
        #ifdef MEXOCT_MATLAB
            mexErrMsgIdAndTxt(id, s);
        #else
            error_with_id(id, "%s", s);
        #endif
    }

    template<typename... Args>
    void mexoct_error(const char *id, const char * s, Args... args) {
        #ifdef MEXOCT_MATLAB
            mexErrMsgIdAndTxt(id, s, args...);
        #else
            error_with_id(id, s, args...);
        #endif
    }

    void mexoct_error(const char *id, const Utils::string & s) {
        mexoct_error(id, s.c_str());
    }

    void mexoct_error(const MexOctException &e) {
        mexoct_error(e.id(), e.what());
    }

    void mexoct_error(const std::exception &e) {
        mexoct_error("std::exception", e.what());
    }

    void mexoct_warning(const char * s) {
        #ifdef MEXOCT_MATLAB
            mexWarnMsgTxt(s);
        #else
            warning("%s", s);
        #endif
    }

    void mexoct_warning(const Utils::string & s) {
        mexoct_warning(s.c_str());
    }

    void mexoct_warning(const char *id, const char * s) {
        #ifdef MEXOCT_MATLAB
            mexWarnMsgIdAndTxt(id, s);
        #else
            warning_with_id(id, "%s", s);
        #endif
    }

    template<typename... Args>
    void mexoct_warning(const char *id, const char * s, Args... args) {
        #ifdef MEXOCT_MATLAB
            mexWarnMsgIdAndTxt(id, s, args...);
        #else
            warning_with_id(id, s, args...);
        #endif
    }

    void mexoct_warning(const char *id, const Utils::string & s) {
        mexoct_warning(id, s.c_str());
    }


    void mexoct_warning(const MexOctException &e) {
        mexoct_warning(e.id(), e.what());
    }

    void mexoct_warning(const std::exception &e) {
        mexoct_warning("std::exception", e.what());
    }

}



/*
 * Copyright (C) Martin Koehler, 2017-2023
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include <list>
#include <string>

#include "mepack.h"

namespace mexoct {

    class MepackOptions {
        public:
            MepackOptions() : m_mb(-1), m_nb(-1), m_bignb(-1), m_isolver(-1), m_fsolver(-1), m_openmp(1),
                              m_sign(1), m_sign2(1), m_maxit(15), m_tau(1.0) {
#ifdef MEXOCT_MATLAB
                                  mepack_openmp_enable(0);
                                  m_openmp = 0;
#else
                                  mepack_openmp_enable(1);
                                  m_openmp = 1;
#endif
                              }
            ~MepackOptions() {}

            void setMB(int mb) { m_mb = mb; }
            void setNB(int nb) { m_nb = nb; }
            void setBigNB(int bigmb) { m_bignb = bigmb; }
            void setISolver(int isolver) { m_isolver = isolver; }
            void setFSolver(int fsolver) { m_fsolver = fsolver; }
            void setOpenMP(int openmp) {
                mepack_openmp_enable(openmp);
                m_openmp = openmp;

            }
            void setSign(int sign) { m_sign = sign; }
            void setSign2(int sign) { m_sign2 = sign; }
            void setMaxit(size_t maxit) { m_maxit = maxit;}
            void setTau(double tau) { m_tau = tau;}


            int getMB() const { return m_mb; }
            int getNB() const { return m_nb; }
            int getBigNB() const { return m_bignb; }
            int getISolver() const { return m_isolver; }
            int getFSolver() const { return m_fsolver; }
            int getOpenMP() const { return m_openmp; }
            int getSign() const { return m_sign; }
            int getSign2() const { return m_sign2; }
            size_t getMaxit() const { return m_maxit; }
            double getTau() const {return m_tau;}

            static std::list<std::string> getMembers() {
                return std::list<std::string>{"mb", "nb", "bignb", "isolver", "fsolver", "openmp", "sign", "sign2", "maxit", "tau"};
            }
        private:
            int m_mb;
            int m_nb;
            int m_bignb;
            int m_isolver;
            int m_fsolver;
            int m_openmp;
            int m_sign;
            int m_sign2;
            size_t m_maxit;
            double m_tau;
    };

    template<>
        struct ParamType<MepackOptions> {
#ifdef MEXOCT_MATLAB
            static bool isValid(const mxArray *_array)
            {
                if ( ! mxIsStruct(_array) ) return false;

                auto lst = MepackOptions::getMembers();
                int nfields = mxGetNumberOfFields(_array);
                for ( int i = 0; i < nfields; i++) {
                    std::string key(mxGetFieldNameByNumber(_array, i));
                    bool found = false;
                    for (auto memb: lst) {
                        if ( memb == key) found = true;
                    }
                    if (!found) {
                        mexoct_error("MEPACK:MepackOptions", "option structure contains unknown key: %s", key.c_str());
                        return false;
                    }

                    mxArray *val = mxGetFieldByNumber(_array, 0, i);
                    if ( ! (mxIsNumeric(val) && ! mxIsComplex(val))) {
                        mexoct_error("MEPACK:MepackOptions", "structure component %s needs to be a real or integer", key.c_str());
                        return false;
                    }
                    if ( !(mxIsScalar(val))) {
                        mexoct_error("MEPACK:MepackOptions", "structure component %s needs to be a scalar", key.c_str());
                        return false;
                    }
                }

                return true;
            }

            static MepackOptions copy_in(const mxArray *_array)
            {
                MepackOptions optset;
                mxArray * val;
                if ( (val = mxGetField(_array, 0, "mb"))) {
                    optset.setMB((int) mxGetScalar(val));
                }
                if ( (val = mxGetField(_array, 0, "nb"))) {
                    optset.setNB((int) mxGetScalar(val));
                }
                if ( (val = mxGetField(_array, 0, "bignb"))) {
                    optset.setBigNB((int) mxGetScalar(val));
                }
                if ( (val = mxGetField(_array, 0, "isolver"))) {
                    optset.setISolver((int) mxGetScalar(val));
                }
                if ( (val = mxGetField(_array, 0, "fsolver"))) {
                    if ( mxGetScalar(val) == 3 || mxGetScalar(val) == 4 ) {
                        mexoct_warning("MEPACK:MatlabOpenMP", "Due to the handling of OpenMP code in MATLAB, using OpenMP with MEPACK/MATLAB lead to crashes. If you know what you are do, this warning can be deactivated by warning('OFF','MEPACK:MatlabOpenMP').");
                    }
                    optset.setFSolver((int) mxGetScalar(val));
                }
                if ( (val = mxGetField(_array, 0, "openmp"))) {
                    if ( mxGetScalar(val) != 0 ) {
                        mexoct_warning("MEPACK:MatlabOpenMP", "Due to the handling of OpenMP code in MATLAB, using OpenMP with MEPACK/MATLAB lead to crashes. If you know what you are do, this warning can be deactivated by warning('OFF','MEPACK:MatlabOpenMP').");
                    }
                    optset.setOpenMP((int) mxGetScalar(val));
                }

                if ( (val = mxGetField(_array, 0, "sign"))) {
                    optset.setSign((int) mxGetScalar(val));
                }

                if ( (val = mxGetField(_array, 0, "sign2"))) {
                    optset.setSign2((int) mxGetScalar(val));
                }

                if ( ( val = mxGetField(_array, 0, "maxit"))) {
                    optset.setMaxit((size_t) mxGetScalar(val));
                }

                if ( (val = mxGetField(_array, 0, "tol") )) {
                    optset.setTau((double) mxGetScalar(val));
                }
                return optset;
            }

            static mxArray * copy_out(MepackOptions & _array)
            {
                return NULL;
            }

            static constexpr const char* NAME = "MepackOptions";

#else

            static bool isValid(const octave_value & _array)
            {
                if ( ! _array.isstruct()) return false;
                octave_scalar_map map = _array.scalar_map_value();
                auto lst = MepackOptions::getMembers();

                for ( auto mapiter = map.begin(); mapiter != map.end(); mapiter ++) {
                    std::string key = map.key(mapiter);
                    bool found = false;
                    for (auto memb: lst) {
                        if ( memb == key) found = true;
                    }
                    if (!found) {
                        mexoct_error("MEPACK:MepackOptions", "option structure contains unknown key: %s", key.c_str());
                        return false;
                    }


                    octave_value& val = map.contents(key);
                    if ( !val.is_scalar_type()) {
                        mexoct_error("MEPACK:MepackOptions", "structure component %s needs to be a scalar", key.c_str());
                        return false;
                    }
                    if ( !(val.isreal() || val.isinteger())) {
                        mexoct_error("MEPACK:MepackOptions", "structure component %s needs to be a real or integer", key.c_str());
                        return false;
                    }
                }
                return true;
            }

            static MepackOptions copy_in(const octave_value & _array)
            {
                MepackOptions optset;
                octave_scalar_map map = _array.scalar_map_value();
                if ( map.contains("mb")) {
                    optset.setMB(map.contents("mb").int_value());
                }
                if ( map.contains("nb")) {
                    optset.setNB(map.contents("nb").int_value());
                }
                if ( map.contains("bignb")) {
                    optset.setBigNB(map.contents("bignb").int_value());
                }
                if ( map.contains("isolver")) {
                    optset.setISolver(map.contents("isolver").int_value());
                }
                if ( map.contains("fsolver")) {
                    optset.setFSolver(map.contents("fsolver").int_value());
                }
                if ( map.contains("openmp")) {
                    optset.setOpenMP(map.contents("openmp").int_value());
                }
                if ( map.contains("sign")) {
                    optset.setSign(map.contents("sign").int_value());
                }
                if ( map.contains("sign2")) {
                    optset.setSign2(map.contents("sign2").int_value());
                }
                if ( map.contains("maxit")) {
                    optset.setMaxit((size_t) map.contents("maxit").int_value());
                }
                if ( map.contains("tau")) {
                    optset.setTau(map.contents("tau").double_value());
                }

                return optset;
            }

            static octave_value copy_out(std::string _array)
            {
                return octave_value();
            }

            static constexpr const char* NAME = "MepackOptions";


#endif
        };
}


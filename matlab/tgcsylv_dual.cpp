/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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
 * Copyright (C) Martin Koehler, 2017-2022
 */

#define MEXOCT_MATLAB_DQ_STRING

#include "mepack.h"
#include "mexoct/interface.hpp"
#include "mepackoptions.hpp"

using namespace mexoct;

#ifdef MEXOCT_OCTAVE
#define ROWS(a)  (a.rows())
#define COLS(a)  (a.columns())
#define ARG_T      const octave_value &
#else
#define ROWS(a)  (mxGetM(a))
#define COLS(a)  (mxGetN(a))
#define ARG_T    const mxArray *
#endif

template<typename T>
struct tcgd_solver {
    constexpr static const char *fn = NULL;
    constexpr static const char *fn_openmp = NULL;

    virtual void solve(const char *opA, const char * opB, T sign1, T sign2, int M, int N,       T * Aptr, int lda, T* Bptr, int ldb, T* Cptr, int ldc, T* Dptr, int ldd, T *Yptr, int ldy, T *Xptr, int ldx, T * scale, T * buffer, int *info) = 0;
    virtual void solve_openmp(const char *opA, const char *opB, T sign1, T sign2, int M, int N, T * Aptr, int lda, T* Bptr, int ldb, T* Cptr, int ldc, T* Dptr, int ldd, T *Yptr, int ldy, T *Xptr, int ldx, T * scale, T * buffer, int *info) = 0;
    virtual void blocksize_mb_set(int mb) = 0;
    virtual void blocksize_nb_set(int nb) = 0;
    virtual void isolver_set(int) = 0;
};

template<>
struct tcgd_solver<double> {
    constexpr static const char *fn = "DLA_TGCSYLV_DUAL_L3";
    constexpr static const char *fn_openmp = "DLA_TGCSYLV_DUAL_DAG";

    static void solve(const char *opA, const char *opB, double sign1, double sign2, int M, int N, double * Aptr, int lda, double *Bptr, int ldb, double *Cptr, int ldc, double *Dptr, int ldd,
            double *Yptr, int ldy, double *Xptr, int ldx, double * scale, double * buffer, int *info) {
        mepack_double_tgcsylv_dual_level3(opA, opB, sign1, sign2, M, N, Aptr, lda, Bptr, ldb, Cptr, ldc, Dptr, ldd, Yptr, ldy, Xptr, ldx, scale, buffer, info);
    };

    static void solve_openmp(const char *opA, const char *opB, double sign1, double sign2, int M, int N, double * Aptr, int lda, double *Bptr, int ldb, double *Cptr, int ldc, double *Dptr, int ldd,
            double *Yptr, int ldy, double *Xptr, int ldx, double * scale, double * buffer, int *info) {
        mepack_double_tgcsylv_dual_dag(opA, opB, sign1, sign2, M, N, Aptr, lda, Bptr, ldb, Cptr, ldc, Dptr, ldd, Yptr, ldy, Xptr, ldx, scale, buffer, info);
    };

    static void blocksize_mb_set(int mb) {
        mepack_double_tgcsylv_dual_blocksize_mb_set(mb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_double_tgcsylv_dual_blocksize_nb_set(nb);
    };

    static void isolver_set(int isolver) {
        mepack_tgcsylv_dual_isolver_set(isolver);
    };
};

template<>
struct tcgd_solver<float> {
    constexpr static const char *fn = "SLA_TGCSYLV_DUAL_L3";
    constexpr static const char *fn_openmp = "SLA_TGCSYLV_DUAL_DAG";

    static void solve(const char *opA, const char *opB, float sign1, float sign2, int M, int N, float * Aptr, int lda, float *Bptr, int ldb, float *Cptr, int ldc, float *Dptr, int ldd,
            float *Yptr, int ldy, float *Xptr, int ldx, float * scale, float * buffer, int *info) {
        mepack_single_tgcsylv_dual_level3(opA, opB, sign1, sign2, M, N, Aptr, lda, Bptr, ldb, Cptr, ldc, Dptr, ldd, Yptr, ldy, Xptr, ldx, scale, buffer, info);
    };

    static void solve_openmp(const char *opA, const char *opB, float sign1, float sign2, int M, int N, float * Aptr, int lda, float *Bptr, int ldb, float *Cptr, int ldc, float *Dptr, int ldd,
            float *Yptr, int ldy, float *Xptr, int ldx, float * scale, float * buffer, int *info) {
        mepack_single_tgcsylv_dual_dag(opA, opB, sign1, sign2, M, N, Aptr, lda, Bptr, ldb, Cptr, ldc, Dptr, ldd, Yptr, ldy, Xptr, ldx, scale, buffer, info);
    };

    static void blocksize_mb_set(int mb) {
        mepack_single_tgcsylv_dual_blocksize_mb_set(mb);
    };
    static void blocksize_nb_set(int nb) {
        mepack_single_tgcsylv_dual_blocksize_nb_set(nb);
    };

    static void isolver_set(int isolver) {
        mepack_tgcsylv_dual_isolver_set(isolver);
    };
};



template <typename T>
static std::tuple<ArrayType<T>, ArrayType<T>>
call_solver_op_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D, ArrayType<T> Y, ArrayType<T> X,
                      std::string opA = std::string("N"), std::string opB = std::string("N"),
                      MepackOptions optset = MepackOptions())
{
    if ( A.rows != Y.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments one and five need to have the same number of rows size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) A.rows, (int) A.columns, (int) Y.rows, (int) Y.columns);
    }
    if ( A.rows != X.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments one and six need to have the same number of rows size(A) = [ %d %d ] size (X) = [ %d %d ]", (int) A.rows, (int) A.columns, (int) X.rows, (int) X.columns);
    }
    if ( C.rows != Y.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments three and five need to have the same number of rows size(A) = [ %d %d ] size (Y) = [ %d %d ]", (int) C.rows, (int) C.columns, (int) Y.rows, (int) Y.columns);
    }
    if ( C.rows != X.rows) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments three and six need to have the same number of rows size(A) = [ %d %d ] size (X) = [ %d %d ]", (int) C.rows, (int) C.columns, (int) X.rows, (int) X.columns);
    }


    if ( B.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments two and five need to have the same number of columns size(B) = [ %d %d ] size (Y) = [ %d %d ]", (int) B.rows, (int) B.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( B.columns != X.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments two and six need to have the same number of columns size(B) = [ %d %d ] size (X) = [ %d %d ]", (int) B.rows, (int) B.columns, (int) X.rows, (int) X.columns);
    }


    if ( D.columns != Y.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments four and five to have the same number of columns size(B) = [ %d %d ] size (Y) = [ %d %d ]", (int) D.rows, (int) D.columns, (int) Y.rows, (int) Y.columns);
    }

    if ( D.columns != X.columns) {
        mexoct_error("MEPACK:DimMissMatch",
                "The input arguments four and six to have the same number of columns size(B) = [ %d %d ] size (X) = [ %d %d ]", (int) D.rows, (int) D.columns, (int) X.rows, (int) X.columns);
    }



    if ( optset.getMB() >= 0 ) {
        tcgd_solver<T>::blocksize_mb_set(optset.getMB());
    }
    if ( optset.getNB() >= 0 ) {
        tcgd_solver<T>::blocksize_nb_set(optset.getNB());
    }

    if ( optset.getISolver() >= 0 ) {
        tcgd_solver<T>::isolver_set(optset.getISolver());
    }
    ssize_t auxmem;
    if ( optset.getOpenMP()) {
        auxmem = mepack_memory(tcgd_solver<T>::fn_openmp, A.rows, B.rows);
    } else {
        auxmem = mepack_memory(tcgd_solver<T>::fn, A.rows, B.rows);
    }

    if ( auxmem < 0 ) {
        mexoct_error("MEPACK:mepack_memory", "Failed to obtain the size of the auxiliary memory.");
    }

    MEXOCT_BUFFER(T, auxbuf, auxmem);
    T scale = 1;
    T sign1 = 1;
    T sign2 = 1;
    int info = 0;

    if ( optset.getSign() < 0 ) {
        sign1 = -1;
    }
    if ( optset.getSign2() < 0 ) {
        sign2 = -1;
    }


    if ( optset.getOpenMP()) {
        tcgd_solver<T>::solve_openmp(opA.c_str(), opB.c_str(), sign1, sign2, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld, C.ptr(), C.ld, D.ptr(), D.ld,
                Y.ptr(), Y.ld, X.ptr(), X.ld, &scale, auxbuf, &info);
    } else {
        tcgd_solver<T>::solve(opA.c_str(), opB.c_str(), sign1, sign2, A.rows, B.rows, A.ptr(), A.ld, B.ptr(), B.ld, C.ptr(), C.ld, D.ptr(), D.ld,
                Y.ptr(), Y.ld, X.ptr(), X.ld, &scale, auxbuf, &info);
    }
    if ( info != 0) {
        mexoct_error("MEPACK:TGCSYLV_LEVEL3", "TGCSYLV_LEVEL3 failed with info = %d", (int) info);
    }

    return std::make_tuple(Y, X);


}

template<typename T>
static std::tuple<ArrayType<T>,ArrayType<T>>
call_solver(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D, ArrayType<T> Y, ArrayType<T> X)
{
    return call_solver_op_optset(A, B, C, D, Y, X);
}

template<typename T>
static std::tuple<ArrayType<T>,ArrayType<T>>
call_solver_op(ArrayType<T>  A,ArrayType<T> B,  ArrayType<T> C, ArrayType<T> D, ArrayType<T> Y, ArrayType<T> X, std::string opA, std::string opB)
{
    return call_solver_op_optset(A, B, C, D, Y, X, opA, opB);
}

template<typename T>
static std::tuple<ArrayType<T>, ArrayType<T>>
call_solver_optset(ArrayType<T>  A, ArrayType<T> B, ArrayType<T> C, ArrayType<T> D, ArrayType<T> Y, ArrayType<T> X, MepackOptions optset)
{
// #define DEBUG_ALL
#ifdef DEBUG_ALL
    std::cout << "Optset Function \n";
    std::cout << " - MB " << optset.getMB() << "\n";
    std::cout << " - NB " << optset.getNB() << "\n";
    std::cout << " - BIGNB " << optset.getBigNB() << "\n";
    std::cout << " - ISOLVER " << optset.getISolver() << "\n";
    std::cout << " - FSOLVER " << optset.getFSolver() << "\n";
#endif
    return call_solver_op_optset(A, B, C, D, Y, X, "N", "N",optset);
}



MEXOCT_ENTRY(mepack_tgcsylv_dual,
R"===(
[R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F)
[R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F, OpA, OpB)
[R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F, OptSet)
[R, L] = mepack_tgcsylv_dual(A, B, C, D, E, F, OpA, OpB, OptSet)
    A      -- double precision (quasi) upper triangular coefficient matrix
    B      -- double precision (quasi) upper triangular coefficient matrix
    C      -- double precision upper triangular coefficient matrix
    D      -- double precision upper triangular coefficient matrix,
              can be quasi upper triangular if B is upper triangular
    E, F   -- double precision right hand sides
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

[Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs)
[Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs, Op)
[Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs, OptSet)
[Rs, Ls] = mepack_tgcsylv_dual(As, Bs, Cs, Ds, Es, Fs, Op, OptSet)
    As     -- single precision (quasi) upper triangular coefficient matrix
    Bs     -- single precision (quasi) upper triangular coefficient matrix
    Cs     -- single precision upper triangular coefficient matrix
    Ds     -- single precision upper triangular coefficient matrix,
              can be quasi upper triangular if Bs is upper triangular
    Es, Fs -- single precision right hand sides
    OpA    -- Transpose operation on the left coefficient matrix
    OpB    -- Transpose operation on the right coefficient matrix
    OptSet -- Options structure to override the default settings
              of MEPACK.

The mepack_tgcsylv_dual function solves the dual coupled generalized Sylvester equation

    opA(A)' * R            + opA(C)' * R     = E
          +/- R * opB(B)'  +/- L * opB(D)'   = F

where (A, C) and (B,D) (or (D,B)) are (quasi) upper triangular matrix pairs,
both coming from a generalized Schur decomposition. The operation op
either returns the matrix or if Op == 'T' its transpose. In this way, the
transposed equation can be solved without transposing the coefficient
matrix before.

The function uses the level-3 solver D/SLA_TGCSYLV_DUAL_* of MEPACK.

If the OptSet argument is given, the default settings of MEPACK can
be overwritten this structure can have the following members for
mepack_tgcsylv_dual:

   - mb      -- row block size of the level-3 algorithm (mb >=2)
   - nb      -- column block size of the level-3 algorithm (nb >=2)
   - isolver -- internal level-2 solver in the algorithm
       = 0   -- Default/Automatic selection
       = 1   -- Local Copies, with alignments (if supported by the compiler)
       = 2   -- Local Copies, without alignments
       = 3   -- Reordering and Fortran Vectorization
       = 4   -- Classic level-2 solver
       = 5   -- Naive level-2 solver
       = 6   -- Recursive Blocking
   - sign    -- swap the sign between both terms in the first half of the second equation
       = 1   -- solve the equation with '+' in the first half of the second equation
       = -1  -- solve the equation with '-' in the first half of the second equation
   - sign2   -- swap the sign between both terms in the second half of the second equation
       = 1   -- solve the equation with '+' in the second half of the second equation
       = -1  -- solve the equation with '-' in the second half of the second equation
   - openmp  -- OpenMP selection
       = 1   -- OpenMP - DAG solver enabled and OpenMP accelerated operations
       = 0   -- Standard level-3 solver without OpenMP acceleration

This file file is part of MEPACK, the Matrix Equation Package
by Martin Koehler.
)===")
{

    MEXOCT_INIT();
    mepack_init();

    auto ParamOpAString = Parameter<std::string>("OpA", "Transpose operation on the left coefficient matrix");
    auto ParamOpBString = Parameter<std::string>("OpB", "Transpose operation on the right coefficient matrix");
    auto ParamOptSet   = Parameter<MepackOptions>("optset", "Options for MEPACK");

    auto ParamDoubleA = Parameter<ArrayType<double>, Traits::SquareMatrix>("A", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleB = Parameter<ArrayType<double>, Traits::SquareMatrix>("B", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleC = Parameter<ArrayType<double>, Traits::SquareMatrix>("C", "double precision upper triangular coefficient matrix");
    auto ParamDoubleD = Parameter<ArrayType<double>, Traits::SquareMatrix>("D", "double precision (quasi) upper triangular coefficient matrix");
    auto ParamDoubleE = Parameter<ArrayType<double>>("E", "double precision right hand side");
    auto ParamDoubleF = Parameter<ArrayType<double>>("F", "double precision right hand side");



    auto ParamFloatA = Parameter<ArrayType<float>, Traits::SquareMatrix>("As", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatB = Parameter<ArrayType<float>, Traits::SquareMatrix>("Bs", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatC = Parameter<ArrayType<float>, Traits::SquareMatrix>("Cs", "single precision upper triangular coefficient matrix");
    auto ParamFloatD = Parameter<ArrayType<float>, Traits::SquareMatrix>("Ds", "single precision (quasi) upper triangular coefficient matrix");
    auto ParamFloatE = Parameter<ArrayType<float>>("Es", "single precision right hand side");
    auto ParamFloatF = Parameter<ArrayType<float>>("Fs", "single precision right hand side");



    auto tgcsylv_double            = makeFunctionExt<2>(call_solver<double>, "[R, L]",           ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF);
    auto tgcsylv_double_optset     = makeFunctionExt<2>(call_solver_optset<double>, "[R, L]",    ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOptSet);
    auto tgcsylv_double_op         = makeFunctionExt<2>(call_solver_op<double>, "[R, L]",        ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString);
    auto tgcsylv_double_op_optset  = makeFunctionExt<2>(call_solver_op_optset<double>, "[R, L]", ParamDoubleA, ParamDoubleB, ParamDoubleC, ParamDoubleD, ParamDoubleE, ParamDoubleF, ParamOpAString, ParamOpBString,ParamOptSet);



    auto tgcsylv_float            = makeFunctionExt<2>(call_solver<float>, "[Rs, Ls]",             ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF);
    auto tgcsylv_float_optset     = makeFunctionExt<2>(call_solver_optset<float>, "[Rs, Ls]",      ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOptSet);
    auto tgcsylv_float_op         = makeFunctionExt<2>(call_solver_op<float>, "[Rs, Ls]",          ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString);
    auto tgcsylv_float_op_optset  = makeFunctionExt<2>(call_solver_op_optset<float>, "[Rs, Ls]",   ParamFloatA, ParamFloatB, ParamFloatC, ParamFloatD, ParamFloatE, ParamFloatF, ParamOpAString, ParamOpBString, ParamOptSet);


	auto parser = makeArgParser("mepack_tgcsylv_dual", tgcsylv_double,  tgcsylv_double_optset, tgcsylv_double_op, tgcsylv_double_op_optset,
                                                 tgcsylv_float, tgcsylv_float_op, tgcsylv_float_optset, tgcsylv_float_op_optset);

    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}



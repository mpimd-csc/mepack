/*
 * MEXOCT - QR examples
 * Copyright © 2021 Martin Koehler
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

#include "mexoct/interface.hpp"
using namespace mexoct;


/*
 * Definition of the Fortran interfaces of the QR Decomposition
 */

extern "C" void dgeqrf_( mexoct_fortran_int *M, mexoct_fortran_int *N, double *A, mexoct_fortran_int *LDA,
                        double *TAU, double *WORK, mexoct_fortran_int *LWORK, mexoct_fortran_int *INFO);
extern "C" void sgeqrf_( mexoct_fortran_int *M, mexoct_fortran_int *N, float *A, mexoct_fortran_int *LDA,
                        float *TAU, float *WORK, mexoct_fortran_int *LWORK, mexoct_fortran_int *INFO);


/*
 * Double precision routine
 */
std::tuple<ArrayType<double>, ArrayType<double>>
getrf_interface_double(ArrayType<double>  A)
{
    mexoct_fortran_int lwork;
    mexoct_fortran_int info;
    mexoct_fortran_int M = A.rows;
    mexoct_fortran_int N = A.columns;
    mexoct_fortran_int LD = A.ld;

    mexoct_fortran_int mn = (M < N ) ? M:N;

    lwork = A.columns * 64 * sizeof(double);

    ArrayType<double> TAU(mn, 1, new double[mn], false);

    MEXOCT_BUFFER(double, work, lwork);

    dgeqrf_(&M, &N, (double*) A, &LD, (double *) TAU, work, &lwork, &info);

    if (info != 0) {
        mexoct_error("GEQRF returned with an error.");
    }

    return std::make_tuple(A, TAU);
}

/*
 * Single precision routine
 */
std::tuple<ArrayType<float>, ArrayType<float>>
getrf_interface_float(ArrayType<float>  A)
{
    mexoct_fortran_int lwork;
    mexoct_fortran_int info;
    mexoct_fortran_int M = A.rows;
    mexoct_fortran_int N = A.columns;
    mexoct_fortran_int LD = A.ld;

    mexoct_fortran_int mn = (M < N ) ? M:N;

    lwork = A.columns * 64 * sizeof(float);

    ArrayType<float> TAU(mn, 1, new float[mn], false);

    MEXOCT_BUFFER(float, work, lwork);

    sgeqrf_(&M, &N, (float*) A, &LD, (float *) TAU, work, &lwork, &info);

    if (info != 0) {
        mexoct_error("GEQRF returned with an error.");
    }

    return std::make_tuple(A, TAU);

}


/*
 * Interface function
 */

MEXOCT_ENTRY(geqrf, "\
geqrf(A) -- Interface to GEQRF from LAPACK\n\
 \n\
    [AQR,TAU] = GEQRF(A) computes the QR decompostion of a matrix using\n\
    the QR decomposition from LAPACK. \n\
    \n\
    Inputs:\n\
        A    -  Input matrix A.\n\
        \n\
    Outputs:\n\
        AQR  -  QR decomposotion of A as computed by LAPACK GEQRF.\n\
                triu(AQR) contains the matrix R and tril(A,-1) contains\n\
                the Housevectors without ones on the diagonal.\n\
        TAU  -  Vector containing the scaling factors for the Householder\n\
                transformations.\n\
    \n")
{
    /* Initialize MexOct */
    MEXOCT_INIT();


    /* Define the possible parameters */
    auto ParamDoubleMatrix = Parameter<ArrayType<double>>("A", "double precision input matrix.");
    auto ParamSingleMatrix = Parameter<ArrayType<float>>("A", "single precision input matrix.");

    /* Define the possible function calls */
    auto dgetrf_double = makeFunctionExt<2>(getrf_interface_double,Utils::string("[AQR, TAU]"), ParamDoubleMatrix);
    auto dgetrf_float  = makeFunctionExt<2>(getrf_interface_float,"[AQR, TAU]", ParamSingleMatrix);

    /* Create a parser out of the possible function calls */
	auto parser = makeArgParser("geqrf", dgetrf_double, dgetrf_float);

    /* Execute the parser */
    MEXOCT_PARSE(parser);

    MEXOCT_RETURN;
}

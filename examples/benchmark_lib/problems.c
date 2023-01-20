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
 * Copyright (C) Martin Koehler, 2017-2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "cscutils/hdf.h"

#include "benchmark.h"

void benchmark_random_evp_double(Int M, Int *iseed,double *Ar, double *Q, double *A, Int sort)
{
        Int N2 = M *M ;
    	Int IDIST = 2;
        char groupname[128];
        char groupname_outer[128], all[512];
        hid_t h5file, group_outer, group;
        double *iA;
        double ts, te;
        double tsc, tse;
        Int NB = sort;
        char hdffile[PATH_LEN];


        snprintf(groupname_outer, 128, "N=%d", (int) M);
        snprintf(groupname, 128, "%d-%d-%d-%d-%d", (int) IDIST,(int) iseed[0], (int) iseed[1], (int) iseed[2], (int) iseed[3]);

        iA  = (double *) malloc(sizeof(double) * (M*M));
        FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, iA);

        snprintf(hdffile, PATH_LEN, "%s/qr/qr_%06d.h5", hdfstore_path, (int) M);
        h5file = csc_hdf5_open(hdffile, "rw");

        snprintf(all, 512, "%s/%s", groupname_outer, groupname);

        if ( csc_hdf5_group_exist(h5file, all) ) {
            /* printf("%s - exists. load and return.\n", all); */
            free(iA);
            group_outer = csc_hdf5_group_open(h5file, groupname_outer);
            group = csc_hdf5_group_open(group_outer, groupname);

            if ( Ar != NULL ) {
                csc_hdf5_matrix_read_real(group, "Aschur", M, M, M, Ar);
            }

            if ( A != NULL ) {
                csc_hdf5_matrix_read_real(group, "A", M, M, M, A);
            }
            if ( Q != NULL ) {
                csc_hdf5_matrix_read_real(group, "Q", M, M, M, Q);
            }


            if ( sort ) {
                Int info;
                Int freeQ = 0;
                Int ldwork = 256*M+16;
                double *work = (double *) malloc(sizeof(double) * (ldwork));

                if ( Q == NULL ) {
                    Q  = (double *) malloc(sizeof(double) * (M*M));
                    freeQ = 1;
                    csc_hdf5_matrix_read_real(group, "Q", M, M, M, Q);
                }
                FC_GLOBAL(dla_sort_ev,DLA_SORT_EV)( &M,  Ar, &M, Q, &M, &NB, work, &ldwork, &info );
                if ( freeQ ) free(Q);
                free(work);
            }

            csc_hdf5_group_close(group);
            csc_hdf5_group_close(group_outer);
            csc_hdf5_close(h5file);

            return;
        }



        group_outer = csc_hdf5_group_open(h5file, groupname_outer);
        group = csc_hdf5_group_open(group_outer, groupname);


        if (sizeof(Int) == 4 ) {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "idist", 1, &IDIST);

        } else {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "idist", 1, &IDIST);
        }

        csc_hdf5_matrix_write_real(group, "A", M, M, M, iA);
        if ( A != NULL ) {
            FC_GLOBAL(dlacpy,DLACPY)("All", &M, &M, iA, &M, A, &M, 1);
        }


        /*-----------------------------------------------------------------------------
         *  Solve EVP
         *-----------------------------------------------------------------------------*/
        double * alphar = malloc(sizeof(double) *3*M);
        double * alphai = alphar+M;
        Int ldwork = 256*M+16;
        double *work = (double *) malloc(sizeof(double) * (ldwork));
        Int *bwork = (Int *) malloc(sizeof(Int) * (M));
        Int sdim;
        Int freeQ = 0;
        Int info;

        if ( Q == NULL ) {
            Q = (double *) malloc(sizeof(double) * (M*M));
            freeQ = 1;
        }


        ts = get_wtime();
        tsc = get_ctime();
        FC_GLOBAL(dgees,DGEES)("V","N", NULL, &M, iA, &M, &sdim, alphar, alphai, Q, &M, work, &ldwork, bwork, &info);
        tse = get_ctime();
        te = get_wtime();
        csc_hdf5_matrix_write_real(group, "Aschur", M, M, M, iA);
        csc_hdf5_matrix_write_real(group, "Q", M, M, M, Q);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphar", M, alphar);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphai", M, alphai);
        te = te - ts;
        tse = tse - tsc;
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "walltime", 1, &te);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "cputime", 1, &tse);

        if ( sort ) {
            FC_GLOBAL(dla_sort_ev,DLA_SORT_EV)( &M,  iA, &M, Q, &M, &NB, work, &ldwork, &info );
        }

        if ( Ar != NULL )
            FC_GLOBAL(dlacpy,DLACPY)("All", &M, &M, iA, &M, Ar, &M, 1);

        if ( freeQ ) free(Q);
        free(alphar);
        free(work);
        free(bwork);
        free(iA);

        csc_hdf5_group_close(group);
        csc_hdf5_group_close(group_outer);
        csc_hdf5_close(h5file);

        return ;

}



/*-----------------------------------------------------------------------------
 *  Single Precision
 *-----------------------------------------------------------------------------*/
void benchmark_random_evp_float(Int M, Int *iseed,float *Ar, float *Q, float *A, Int sort)
{
        Int N2 = M *M ;
    	Int IDIST = 2;
        char groupname[128];
        char groupname_outer[128], all[512];
        hid_t h5file, group_outer, group;
        float *iA;
        float ts, te;
        float tsc, tse;
        Int NB = sort;
        char hdffile[PATH_LEN];


        snprintf(groupname_outer, 128, "N=%d", (int) M);
        snprintf(groupname, 128, "%d-%d-%d-%d-%d", (int) IDIST, (int) iseed[0], (int) iseed[1], (int) iseed[2],(int)  iseed[3]);

        iA  = (float *) malloc(sizeof(float) * (M*M));
        FC_GLOBAL(slarnv,SLARNV)(&IDIST, iseed, &N2, iA);

        snprintf(hdffile, PATH_LEN, "%s/qr/qr_float_%06d.h5", hdfstore_path, (int) M);
        h5file = csc_hdf5_open(hdffile, "rw");

        snprintf(all, 512, "%s/%s", groupname_outer, groupname);

        if ( csc_hdf5_group_exist(h5file, all) ) {
            // printf("%s - exists. load and return.\n", all);
            free(iA);
            group_outer = csc_hdf5_group_open(h5file, groupname_outer);
            group = csc_hdf5_group_open(group_outer, groupname);

            if ( Ar != NULL ) {
                csc_hdf5_matrix_read_real_single(group, "Aschur", M, M, M, Ar);
            }

            if ( A != NULL ) {
                csc_hdf5_matrix_read_real_single(group, "A", M, M, M, A);
            }
            if ( Q != NULL ) {
                csc_hdf5_matrix_read_real_single(group, "Q", M, M, M, Q);
            }


            if ( sort ) {
                Int info;
                Int freeQ = 0;
                Int ldwork = 256*M+16;
                float *work = (float *) malloc(sizeof(float) * (ldwork));

                if ( Q == NULL ) {
                    Q  = (float *) malloc(sizeof(float) * (M*M));
                    freeQ = 1;
                    csc_hdf5_matrix_read_real_single(group, "Q", M, M, M, Q);
                }
                FC_GLOBAL(sla_sort_ev,SLA_SORT_EV)( &M,  Ar, &M, Q, &M, &NB, work, &ldwork, &info );
                if ( freeQ ) free(Q);
                free(work);
            }

            csc_hdf5_group_close(group);
            csc_hdf5_group_close(group_outer);
            csc_hdf5_close(h5file);

            return;
        }



        group_outer = csc_hdf5_group_open(h5file, groupname_outer);
        group = csc_hdf5_group_open(group_outer, groupname);


        if (sizeof(Int) == 4 ) {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "idist", 1, &IDIST);

        } else {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "idist", 1, &IDIST);
        }

        csc_hdf5_matrix_write_real_single(group, "A", M, M, M, iA);
        if ( A != NULL ) {
            FC_GLOBAL(slacpy,SLACPY)("All", &M, &M, iA, &M, A, &M, 1);
        }


        /*-----------------------------------------------------------------------------
         *  Solve EVP
         *-----------------------------------------------------------------------------*/
        float * alphar = malloc(sizeof(float) *3*M);
        float * alphai = alphar+M;
        Int ldwork = 256*M+16;
        float *work = (float *) malloc(sizeof(float) * (ldwork));
        Int *bwork = (Int *) malloc(sizeof(Int) * (M));
        Int sdim;
        Int freeQ = 0;
        Int info;

        if ( Q == NULL ) {
            Q = (float *) malloc(sizeof(float) * (M*M));
            freeQ = 1;
        }


        ts = get_wtime();
        tsc = get_ctime();
        FC_GLOBAL(sgees,SGEES)("V","N", NULL, &M, iA, &M, &sdim, alphar, alphai, Q, &M, work, &ldwork, bwork, &info);
        tse = get_ctime();
        te = get_wtime();
        csc_hdf5_matrix_write_real_single(group, "Aschur", M, M, M, iA);
        csc_hdf5_matrix_write_real_single(group, "Q", M, M, M, Q);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphar", M, alphar);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphai", M, alphai);
        te = te - ts;
        tse = tse - tsc;
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "walltime", 1, &te);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "cputime", 1, &tse);

        if ( sort ) {
            FC_GLOBAL(sla_sort_ev,SLA_SORT_EV)( &M,  iA, &M, Q, &M, &NB, work, &ldwork, &info );
        }

        if ( Ar != NULL )
            FC_GLOBAL(slacpy,SLACPY)("All", &M, &M, iA, &M, Ar, &M, 1);

        if ( freeQ ) free(Q);
        free(alphar);
        free(work);
        free(bwork);
        free(iA);

        csc_hdf5_group_close(group);
        csc_hdf5_group_close(group_outer);
        csc_hdf5_close(h5file);

        return ;

}


void benchmark_random_gevp_double(Int M, Int *iseed, double *Ar, double *Cr, double *Q, double *Z, double *Ao, double *Co,  Int sort) {
        Int N2 = M *M ;
    	Int IDIST = 2;
        char groupname[128];
        char groupname_outer[128], all[512];
        hid_t h5file, group_outer, group;
        double *A, *C;
        double ts, te;
        double tsc, tse;
        Int NB = sort;
        char hdffile[PATH_LEN];

        snprintf(groupname_outer, 128, "N=%d", (int) M);
        snprintf(groupname, 128, "%d-%d-%d-%d-%d", (int) IDIST, (int) iseed[0], (int) iseed[1], (int) iseed[2], (int) iseed[3]);

        A  = (double *) malloc(sizeof(double) * (M*M));
        C  = (double *) malloc(sizeof(double) * (M*M));
        FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, A);
    	FC_GLOBAL(dlarnv,DLARNV)(&IDIST, iseed, &N2, C);

        if ( Ao != NULL ) {
            FC_GLOBAL(dlacpy,DLACPY)("All", &M, &M, A, &M, Ao, &M, 1);
        }
        if ( Co != NULL ) {
            FC_GLOBAL(dlacpy,DLACPY)("All", &M, &M, C, &M, Co, &M, 1);
        }


        snprintf(hdffile, PATH_LEN, "%s/qz/qz_%06d.h5", hdfstore_path, (int) M);
        h5file = csc_hdf5_open(hdffile, "rw");


        snprintf(all, 512, "%s/%s", groupname_outer, groupname);
        if ( csc_hdf5_group_exist(h5file, all) ) {
            // printf("%s - exists. load and return.\n", all);
            free(A);
            free(C);
            group_outer = csc_hdf5_group_open(h5file, groupname_outer);
            group = csc_hdf5_group_open(group_outer, groupname);
            csc_hdf5_matrix_read_real(group, "Aschur", M, M, M, Ar);
            csc_hdf5_matrix_read_real(group, "Bschur", M, M, M, Cr);

            if ( Q != NULL ) {
                csc_hdf5_matrix_read_real(group, "Q", M, M, M, Q);
            }
            if ( Z != NULL ) {
                csc_hdf5_matrix_read_real(group, "Z", M, M, M, Z);
            }

            if ( Ao != NULL ) {
                csc_hdf5_matrix_read_real(group, "A", M, M, M, Ao);
            }
            if ( Co != NULL ) {
                csc_hdf5_matrix_read_real(group, "B", M, M, M, Co);
            }


            if ( sort ) {
                Int info;
                Int freeQ = 0, freeZ = 0;
                Int ldwork = 256*M+16;
                double *work = (double *) malloc(sizeof(double) * (ldwork));

                if ( Q == NULL) {
                    Q  = (double *) malloc(sizeof(double) * (M*M));
                    csc_hdf5_matrix_read_real(group, "Q", M, M, M, Q);
                    freeQ = 1;
                }
                if ( Z == NULL ) {
                    Z  = (double *) malloc(sizeof(double) * (M*M));
                    csc_hdf5_matrix_read_real(group, "Z", M, M, M, Z);
                    freeZ = 1;
                }

                FC_GLOBAL(dla_sort_gev,DLA_SORT_GEV)( &M,  Ar, &M, Cr, &M, Q, &M, Z, &M, &NB, work, &ldwork, &info );
                if ( freeQ) free(Q);
                if ( freeZ) free(Z);
                free(work);
            }

            csc_hdf5_group_close(group);
            csc_hdf5_group_close(group_outer);
            csc_hdf5_close(h5file);

            return;
        }



        group_outer = csc_hdf5_group_open(h5file, groupname_outer);
        group = csc_hdf5_group_open(group_outer, groupname);


        if (sizeof(Int) == 4 ) {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "idist", 1, &IDIST);

        } else {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "idist", 1, &IDIST);
        }

        csc_hdf5_matrix_write_real(group, "A", M, M, M, A);
        csc_hdf5_matrix_write_real(group, "B", M, M, M, C);


        /*-----------------------------------------------------------------------------
         *  Solve GEVP
         *-----------------------------------------------------------------------------*/
        double * alphar = malloc(sizeof(double) *3*M);
        double * alphai = alphar+M;
        double * betar  = alphar+2*M;
        Int ldwork = 256*M+16;
        double *work = (double *) malloc(sizeof(double) * (ldwork));
        Int *bwork = (Int *) malloc(sizeof(Int) * (M));
        Int sdim;
        Int freeQ = 0, freeZ = 0;
        if ( Q == NULL ) {
            Q = (double *) malloc(sizeof(double) * (M*M));
            freeQ = 1;
        }
        if ( Z == NULL ) {
            Z = (double *) malloc(sizeof(double) * (M*M));
            freeZ = 1;
        }
        Int info;

        csc_hdf5_matrix_write_real(group, "A", M, M, M, A);
        csc_hdf5_matrix_write_real(group, "B", M, M, M, C);

        ts = get_wtime();
        tsc = get_ctime();
        FC_GLOBAL(dgges,DGGES)("V","V","N", NULL, &M, A, &M, C,&M, &sdim, alphar, alphai, betar, Q, &M, Z, &M, work, &ldwork, bwork, &info);
        tse = get_ctime();
        te = get_wtime();
        csc_hdf5_matrix_write_real(group, "Aschur", M, M, M, A);
        csc_hdf5_matrix_write_real(group, "Bschur", M, M, M, C);
        csc_hdf5_matrix_write_real(group, "Q", M, M, M, Q);
        csc_hdf5_matrix_write_real(group, "Z", M, M, M, Z);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphar", M, alphar);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphai", M, alphai);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "beta", M, betar);
        te = te - ts;
        tse = tse - tsc;
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "walltime", 1, &te);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "cputime", 1, &tse);

        if ( sort ) {
            FC_GLOBAL(dla_sort_gev,DLA_SORT_GEV)( &M,  A, &M, C, &M, Q, &M, Z, &M, &NB, work, &ldwork, &info );
        }

        FC_GLOBAL(dlacpy,DLACPY)("All", &M, &M, A, &M, Ar, &M, 1);
        FC_GLOBAL(dlacpy,DLACPY)("All", &M, &M, C, &M, Cr, &M, 1);




        if (freeQ) free(Q);
        if (freeZ) free(Z);
        free(alphar);
        free(work);
        free(bwork);
        free(A);
        free(C);

        csc_hdf5_group_close(group);
        csc_hdf5_group_close(group_outer);
        csc_hdf5_close(h5file);

        return ;

}

void FC_GLOBAL(benchmark_random_gevp_double_f,BENCHMARK_RANDOM_GEVP_DOUBLE_F)
    (Int *M, Int *iseed, double *Ar, double *Cr, Int *sort) {
    benchmark_random_gevp_double(*M,  iseed, Ar, Cr, NULL, NULL, NULL, NULL, *sort);

}

void FC_GLOBAL(benchmark_random_gevp_double_full_f, BENCHMARK_RANDOM_GEVP_DOUBLE_FULL_F)
        (Int *M, Int *iseed, double *Ar, double *Cr, double *Q, double *Z, double *Ao, double *Co,  Int *sort) {
        benchmark_random_gevp_double(*M, iseed, Ar, Cr, Q,  Z, Ao, Co, *sort);
}


void benchmark_random_gevp_float(Int M, Int *iseed, float *Ar, float *Cr, float *Q, float *Z, float *Ao, float *Co,  Int sort) {
        Int N2 = M *M ;
    	Int IDIST = 2;
        char groupname[128];
        char groupname_outer[128], all[512];
        hid_t h5file, group_outer, group;
        float *A, *C;
        float ts, te;
        float tsc, tse;
        Int NB = sort;
        char hdffile[PATH_LEN];

        snprintf(groupname_outer, 128, "N=%d", (int) M);
        snprintf(groupname, 128, "%d-%d-%d-%d-%d", (int) IDIST, (int) iseed[0], (int) iseed[1], (int) iseed[2], (int) iseed[3]);

        A  = (float *) malloc(sizeof(float) * (M*M));
        C  = (float *) malloc(sizeof(float) * (M*M));
        FC_GLOBAL(slarnv,SLARNV)(&IDIST, iseed, &N2, A);
    	FC_GLOBAL(slarnv,SLARNV)(&IDIST, iseed, &N2, C);

        if ( Ao != NULL ) {
            FC_GLOBAL(slacpy,SLACPY)("All", &M, &M, A, &M, Ao, &M, 1);
        }
        if ( Co != NULL ) {
            FC_GLOBAL(slacpy,SLACPY)("All", &M, &M, C, &M, Co, &M, 1);
        }


        snprintf(hdffile, PATH_LEN, "%s/qz/qz_float_%06d.h5", hdfstore_path, (int) M);
        h5file = csc_hdf5_open(hdffile, "rw");


        snprintf(all, 512, "%s/%s", groupname_outer, groupname);
        if ( csc_hdf5_group_exist(h5file, all) ) {
            /* printf("%s - exists. load and return.\n", all); */
            free(A);
            free(C);
            group_outer = csc_hdf5_group_open(h5file, groupname_outer);
            group = csc_hdf5_group_open(group_outer, groupname);
            csc_hdf5_matrix_read_real_single(group, "Aschur", M, M, M, Ar);
            csc_hdf5_matrix_read_real_single(group, "Bschur", M, M, M, Cr);

            if ( Q != NULL ) {
                csc_hdf5_matrix_read_real_single(group, "Q", M, M, M, Q);
            }
            if ( Z != NULL ) {
                csc_hdf5_matrix_read_real_single(group, "Z", M, M, M, Z);
            }

            if ( sort ) {
                Int info;
                Int freeQ = 0, freeZ = 0;
                Int ldwork = 256*M+16;
                float *work = (float *) malloc(sizeof(float) * (ldwork));

                if ( Q == NULL) {
                    Q  = (float *) malloc(sizeof(float) * (M*M));
                    csc_hdf5_matrix_read_real_single(group, "Q", M, M, M, Q);
                    freeQ = 1;
                }
                if ( Z == NULL ) {
                    Z  = (float *) malloc(sizeof(float) * (M*M));
                    csc_hdf5_matrix_read_real_single(group, "Z", M, M, M, Z);
                    freeZ = 1;
                }

                FC_GLOBAL(sla_sort_gev,SLA_SORT_GEV)( &M,  Ar, &M, Cr, &M, Q, &M, Z, &M, &NB, work, &ldwork, &info );
                if ( freeQ) free(Q);
                if ( freeZ) free(Z);
                free(work);
            }

            csc_hdf5_group_close(group);
            csc_hdf5_group_close(group_outer);
            csc_hdf5_close(h5file);

            return;
        }



        group_outer = csc_hdf5_group_open(h5file, groupname_outer);
        group = csc_hdf5_group_open(group_outer, groupname);


        if (sizeof(Int) == 4 ) {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER32, group, "idist", 1, &IDIST);

        } else {
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "seed", 4, iseed);
            csc_hdf5_vector_write(CSC_HDF5_INTEGER64, group, "idist", 1, &IDIST);
        }

        csc_hdf5_matrix_write_real_single(group, "A", M, M, M, A);
        csc_hdf5_matrix_write_real_single(group, "B", M, M, M, C);


        /*-----------------------------------------------------------------------------
         *  Solve GEVP
         *-----------------------------------------------------------------------------*/
        float * alphar = malloc(sizeof(float) *3*M);
        float * alphai = alphar+M;
        float * betar  = alphar+2*M;
        Int ldwork = 256*M+16;
        float *work = (float *) malloc(sizeof(float) * (ldwork));
        Int *bwork = (Int *) malloc(sizeof(Int) * (M));
        Int sdim;
        Int freeQ = 0, freeZ = 0;
        if ( Q == NULL ) {
            Q = (float *) malloc(sizeof(float) * (M*M));
            freeQ = 1;
        }
        if ( Z == NULL ) {
            Z = (float *) malloc(sizeof(float) * (M*M));
            freeZ = 1;
        }
        Int info;


        ts = get_wtime();
        tsc = get_ctime();
        FC_GLOBAL(sgges,SGGES)("V","V","N", NULL, &M, A, &M, C,&M, &sdim, alphar, alphai, betar, Q, &M, Z, &M, work, &ldwork, bwork, &info);
        tse = get_ctime();
        te = get_wtime();
        csc_hdf5_matrix_write_real_single(group, "Aschur", M, M, M, A);
        csc_hdf5_matrix_write_real_single(group, "Bschur", M, M, M, C);
        csc_hdf5_matrix_write_real_single(group, "Q", M, M, M, Q);
        csc_hdf5_matrix_write_real_single(group, "Z", M, M, M, Z);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphar", M, alphar);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "alphai", M, alphai);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "beta", M, betar);
        te = te - ts;
        tse = tse - tsc;
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "walltime", 1, &te);
        csc_hdf5_vector_write(CSC_HDF5_REAL, group, "cputime", 1, &tse);

        if ( sort ) {
            FC_GLOBAL(sla_sort_gev,SLA_SORT_GEV)( &M,  A, &M, C, &M, Q, &M, Z, &M, &NB, work, &ldwork, &info );
        }

        FC_GLOBAL(slacpy,SLACPY)("All", &M, &M, A, &M, Ar, &M, 1);
        FC_GLOBAL(slacpy,SLACPY)("All", &M, &M, C, &M, Cr, &M, 1);




        if (freeQ) free(Q);
        if (freeZ) free(Z);
        free(alphar);
        free(work);
        free(bwork);
        free(A);
        free(C);

        csc_hdf5_group_close(group);
        csc_hdf5_group_close(group_outer);
        csc_hdf5_close(h5file);

        return ;

}



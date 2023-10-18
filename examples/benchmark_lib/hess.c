#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "benchmark.h"

#define IDXF2C(k,l,ld) ((k-1)+(l-1)*(ld))
#define IDXF1C(k) (k-1)

//for benchmark problems, LDA=LDB=N
void benchmark_ghess_double(Int N, double *A, double *Q, double *B, double *Z, double *Work, Int ldwork){
    Int M = N;
    Int infoh = 0;
    Int ileft = 1;
    Int iright = M + 1;
    Int iwrk = iright + M;
    Int ilo, ihi;

    FC_GLOBAL(dggbal, DGGBAL)("P", &M, A, &M, B, &M, &ilo, &ihi, &Work[IDXF1C(ileft)], &Work[IDXF1C(iright)], &Work[IDXF1C(iwrk)], &infoh, 1);

    Int irows = ihi + 1 - ilo;
    Int icols = M + 1 - ilo;
    Int itau = iwrk;
    iwrk = itau + irows;
    Int dgeqrf_lwork = ldwork + 1 - iwrk;
    infoh = 0;

    FC_GLOBAL(dgeqrf, DGEQRF)(&irows, &icols, &B[IDXF2C(ilo,ilo,M)], &M,  &Work[IDXF1C(itau)], &Work[IDXF1C(iwrk)], &dgeqrf_lwork, &infoh);

    Int dormqr_lwork = ldwork + 1 - iwrk;

    FC_GLOBAL(dormqr, DORMQR)("L", "T", &irows, &icols, &irows, &B[IDXF2C(ilo,ilo,M)], &M, &Work[IDXF1C(itau)], &A[IDXF2C(ilo,ilo,M)], &M, &Work[IDXF1C(iwrk)], &dormqr_lwork, &infoh, 1, 1);

    double zero = 0.0, one = 1.0;
    FC_GLOBAL(dlaset,DLASET)("Full", &M, &M, &zero, &one, Q, &M, 1);

    if(irows > 1){
        irows = irows - 1;
        FC_GLOBAL_(dlacpy,DLACPY)("L", &irows, &irows, &B[IDXF2C(ilo+1,ilo,M)], &M, &Q[IDXF2C(ilo+1,ilo,M)], &M, 1);
    }
    irows += 1;
    Int dorgqr_lwork = ldwork + 1 - iwrk;
    FC_GLOBAL(dorgqr, DORGQR)(&irows, &irows, &irows, &Q[IDXF2C(ilo,ilo,M)], &M, &Work[IDXF1C(itau)], &Work[IDXF1C(iwrk)], &dorgqr_lwork, &infoh);

    FC_GLOBAL(dlaset,DLASET)("Full", &M, &M, &zero, &one, Z, &M, 1);

    FC_GLOBAL(dgghrd, DGGHRD)("V", "V", &M, &ilo, &ihi, A, &M, B, &M, Q, &M, Z, &M, &infoh, 1, 1);

    //remove all additional information from Hessenberg matrix A
    int i;
    for ( i = 0; i < M*M; i++ ) {
        if ( i > (i / M) * M + (i/M) + 1 && i < ((i/M) + 1) * M) {
            A[i] = 0.0;
	}
    }

    //clear out lower triangle of B
    for( i =  0; i < M*M; i++) {
        if(i > (i / M) * M + (i / M) && i < ((i / M) + 1) * M ) {
            B[i] = 0.0;
        }
    }
}

void benchmark_ghess_single(Int N, float *A, float *Q, float *B, float *Z, float *Work, Int ldwork){
    Int M = N;
    Int infoh = 0;
    Int ileft = 1;
    Int iright = M + 1;
    Int iwrk = iright + M;
    Int ilo, ihi;

    FC_GLOBAL(sggbal,SGGBAL)("P", &M, A, &M, B, &M, &ilo, &ihi, &Work[IDXF1C(ileft)], &Work[IDXF1C(iright)], &Work[IDXF1C(iwrk)], &infoh, 1);

    Int irows = ihi + 1 - ilo;
    Int icols = M + 1 - ilo;
    Int itau = iwrk;
    iwrk = itau + irows;
    Int sgeqrf_lwork = ldwork + 1 - iwrk;
    infoh = 0;

    FC_GLOBAL(sgeqrf,SGEQRF)(&irows, &icols, &B[IDXF2C(ilo,ilo,M)], &M,  &Work[IDXF1C(itau)], &Work[IDXF1C(iwrk)], &sgeqrf_lwork, &infoh);

    Int sormqr_lwork = ldwork + 1 - iwrk;

    FC_GLOBAL(sormqr,SORMQR)("L", "T", &irows, &icols, &irows, &B[IDXF2C(ilo,ilo,M)], &M, &Work[IDXF1C(itau)], &A[IDXF2C(ilo,ilo,M)], &M, &Work[IDXF1C(iwrk)], &sormqr_lwork, &infoh, 1, 1);

    float zero = 0.0f, one = 1.0f;
    FC_GLOBAL(slaset,SLASET)("Full", &M, &M, &zero, &one, Q, &M, 1);

    if(irows > 1){
        irows = irows - 1;
        FC_GLOBAL_(slacpy,SLACPY)("L", &irows, &irows, &B[IDXF2C(ilo+1,ilo,M)], &M, &Q[IDXF2C(ilo+1,ilo,M)], &M, 1);
    }
    irows += 1;
    Int sorgqr_lwork = ldwork + 1 - iwrk;
    FC_GLOBAL(sorgqr, SORGQR)(&irows, &irows, &irows, &Q[IDXF2C(ilo,ilo,M)], &M, &Work[IDXF1C(itau)], &Work[IDXF1C(iwrk)], &sorgqr_lwork, &infoh);

    FC_GLOBAL(slaset,SLASET)("Full", &M, &M, &zero, &one, Z, &M, 1);

    FC_GLOBAL(sgghrd, SGGHRD)("V", "V", &M, &ilo, &ihi, A, &M, B, &M, Q, &M, Z, &M, &infoh, 1, 1);

    //remove all additional information from Hessenberg matrix A
    int i;
    for ( i = 0; i < M*M; i++ ) {
        if ( i > (i / M) * M + (i/M) + 1 && i < ((i/M) + 1) * M) {
            A[i] = 0.0;
        }
    }

    //clear out lower triangle of B
    for( i =  0; i < M*M; i++) {
        if(i > (i / M) * M + (i / M) && i < ((i / M) + 1) * M ) {
            B[i] = 0.0;
        }
    }
}

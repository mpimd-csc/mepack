!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) Martin Koehler, 2017-2023
!

!> \brief Iterative Refinement for the  Coupled Generalized Sylvester Equations.
!
!  Definition:
!  ===========
!
!   SUBROUTINE SLA_GGCSYLV_REFINE(TRANSA, TRANSB, GUESS, SGN1, SGN2, M , N,  &
!        & A, LDA, B, LDB, C, LDC, D, LDD, R, LDR, L, LDL, E, LDE, F, LDF, &
!        & AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,  &
!        & MAXIT, TAU, CONVLOG, &
!        & WORK, LDWORK, INFO)
!
!       CHARACTER, DIMENSION(1) :: TRANSA, TRANSB, GUESS
!       REAL :: SGN1, SGN2
!       INTEGER     :: M, N, LDA, LDB, LDC, LDD, LDE, LDF, LDR, LDL, LDAS, LDBS, LDCS, LDDS, LDQ, LDZ, LDU, LDV, LDWORK, INFO
!       INTEGER :: MAXIT
!       REAL :: TAU
!       REAL :: A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), AS(LDAS, *), BS(LDBS,*), CS(LDCS, *), DS(LDDS, *)
!       REAL :: R(LDR, *), L(LDL, *)
!       REAL :: Q(LDQ, *), Z(LDZ, *), U(LDU, *), V(LDV, *), E ( LDE , * ), F( LDF, *)
!       REAL :: CONVLOG(*), WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SLA_GGCSYLV_REFINE solves a coupled generalized Sylvester equation of the following forms
!>
!>    op1(A) * R  + SGN1 * L  * op2(B) =  E                              (1)
!>    op1(C) * R  + SGN2 * L  * op2(D) =  F
!>
!> with iterative refinement, Thereby  (A,C) is a M-by-M matrix pencil and
!> (B,D) is a N-by-N matrix pencil.
!> The right hand side (E,F) and the solution (R,L) are M-by-N matrices.
!> The pencils (A,C) and (B,D) need to be given in the original form as well
!> as in their generalized Schur decomposition since both are required in the
!> iterative refinement procedure.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER
!>          Specifies the form of the system of equations with respect to A and C :
!>          == 'N':  op1(A) = A
!>          == 'T':  op1(A) = A**T
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER
!>          Specifies the form of the system of equations with respect to B and D:
!>          == 'N':  op2(B) = B,
!>          == 'T':  op2(B) = B**T
!> \endverbatim
!>
!> \param[in] GUESS
!> \verbatim
!>          GUESS is CHARACTER
!>          Specifies whether (R,L) contains an initial guess or nor not.
!>          =  'I': (R, L) contains an initial guess
!>          =  'N': No initial guess, (R,L) is set to zero at the begin of the iteration.
!> \endverbatim
!>
!>
!> \param[in] SGN1
!> \verbatim
!>          SGN1 is REAL, allowed values: +/-1
!>          Specifies the sign between in the first equation.
!> \endverbatim
!>
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is REAL, allowed values: +/-1
!>          Specifies the sign between in the second equation.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices B and D.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,M)
!>          The array A contains the original matrix A defining the equation.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          The array B contains the original matrix B defining the equation.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array A.  LDB >= max(1,N).
!> \endverbatim
!
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (LDC,M)
!>          The array C contains the original matrix C defining the equation.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (LDD,N)
!>          The array D contains the original matrix D defining the equation.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D.  LDD >= max(1,N).
!> \endverbatim
!
!> \param[in,out] R
!> \verbatim
!>          R is REAL array, dimension (LDR,N)
!>          On input, the array R contains the initial guess R0 for the first solution matrix.
!>          On output, the array R contains the refine solution matrix R.
!> \endverbatim
!>
!> \param[in] LDR
!> \verbatim
!>          LDR is INTEGER
!>          The leading dimension of the array R.  LDR >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] L
!> \verbatim
!>          L is REAL array, dimension (LDL,N)
!>          On input, the array L contains the initial guess for the second solution matrix.
!>          On output, the array L contains the solution L.
!> \endverbatim
!>
!> \param[in] LDL
!> \verbatim
!>          LDL is INTEGER
!>          The leading dimension of the array L.  LDF >= max(1,M).
!> \endverbatim

!
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (LDE,N)
!>          On input, the array E contains the right hand side E.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of the array E.  LDE >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is REAL array, dimension (LDF,N)
!>          On input, the array F contains the right hand side F.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F.  LDF >= max(1,M).
!> \endverbatim
!>
!> \param[in] AS
!> \verbatim
!>          AS is REAL array, dimension (LDAS,M)
!>          The array AS contains the generalized Schur decomposition of the
!>          A.
!> \endverbatim
!
!>
!> \param[in] LDAS
!> \verbatim
!>          LDAS is INTEGER
!>          The leading dimension of the array AS.  LDAS >= max(1,M).
!> \endverbatim
!
!>
!> \param[in] BS
!> \verbatim
!>          BS is REAL array, dimension (LDBS,N)
!>          The array BS contains the generalized Schur decomposition of the
!>          B.
!> \endverbatim
!
!>
!> \param[in] LDBS
!> \verbatim
!>          LDBS is INTEGER
!>          The leading dimension of the array BS.  LDBS >= max(1,N).
!> \endverbatim
!
!>
!> \param[in] CS
!> \verbatim
!>          CS is REAL array, dimension (LDCS,M)
!>          The array CS contains the generalized Schur decomposition of the
!>          C.
!> \endverbatim
!
!>
!> \param[in] LDCS
!> \verbatim
!>          LDCS is INTEGER
!>          The leading dimension of the array CS.  LDCS >= max(1,M).
!> \endverbatim
!
!>
!> \param[in] DS
!> \verbatim
!>          DS is REAL array, dimension (LDDS,N)
!>          The array DS contains the generalized Schur decomposition of the
!>          D.
!> \endverbatim
!
!>
!> \param[in] LDDS
!> \verbatim
!>          LDDS is INTEGER
!>          The leading dimension of the array DS.  LDDS >= max(1,N).
!> \endverbatim
!
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,M)
!>          The array Q contains the left generalized Schur vectors for (A,C) as returned by SGGES.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,M).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ,M)
!>          The array Z contains the right generalized Schur vectors for (A,C) as returned by SGGES.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= max(1,M).
!> \endverbatim
!>
!>
!> \param[in,out] U
!> \verbatim
!>          U is REAL array, dimension (LDU,N)
!>          The array U contains the left generalized Schur vectors for (B,D) as returned by SGGES.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,N).
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is REAL array, dimension (LDV,N)
!>          The array V contains the right generalized Schur vectors for (B,D) as returned by SGGES.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.  LDV >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] MAXIT
!> \verbatim
!>          MAXIT is INTEGER
!>          On input, MAXIT contains the maximum number of iteration that are performed, 2 <= MAXIT <= 100
!>          On exit, MAXIT contains the number of iteration steps taken by the algorithm.
!> \endverbatim
!
!> \param[in,out] TAU
!> \verbatim
!>          TAU is REAL
!>          On input, TAU contains the additional security factor for the stopping criterion, typical values are 0.1
!>          On exit, TAU contains the last relative residual when the stopping criterion got valid.
!> \endverbatim
!
!> \param[in,out] CONVLOG
!> \verbatim
!>          CONVLOG is REAL array, dimension (MAXIT)
!>          The CONVLOG array contains the convergence history of the iterative refinement. CONVLOG(I) contains the maximum
!>          relative residual of both equations before it is solved for the I-th time.
!> \endverbatim
!
!> \param[in,out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LDWORK))
!>          Workspace for the algorithm. The optmimal workspace is returned in LDWORK, if LDWORK == -1 on input. In this
!>          case no computations are performed.
!> \endverbatim
!>
!> \param[in,out] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          If LDWORK == -1 the subroutine will return the required size of the workspace in LDWORK on exit. No computations are
!>          performed and none of the arrays are referenced.
!> \endverbatim
!>
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          == 0:  Success
!>          > 0:  Iteration failed in step INFO
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          = -50: Some of the internal settings like NB,... are incorrect.
!> \endverbatim
!>
!> \see SLA_TGCSYLV_DAG
!> \see DLA_TGCSYLV_LEVEL3
!> \see SLA_TGCSYLV_L3_2S
!> \see SLA_TGCSYLV_L2_UNOPT
!> \see SLA_TGCSYLV_L2
!> \see SLA_TGCSYLV_L2_REORDER
!> \see SLA_TGCSYLV_L2_LOCAL_COPY_32
!> \see SLA_TGCSYLV_L2_LOCAL_COPY_64
!> \see SLA_TGCSYLV_L2_LOCAL_COPY_96
!> \see SLA_TGCSYLV_L2_LOCAL_COPY_128
!> \see SLA_TGCSYLV_L2_LOCAL_COPY
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup sglggsylv
!
SUBROUTINE SLA_GGCSYLV_REFINE(TRANSA, TRANSB, GUESS, SGN1, SGN2, M , N,  &
        & A, LDA, B, LDB, C, LDC, D, LDD, R, LDR, L, LDL, E, LDE, F, LDF, &
        & AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,  &
        & MAXIT, TAU, CONVLOG, &
        & WORK, LDWORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_SINGLE
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_FRONTEND_SOLVER

    IMPLICIT NONE

    ! Scalar
    CHARACTER, DIMENSION(1) :: TRANSA, TRANSB, GUESS
    REAL :: SGN1, SGN2
    INTEGER     :: M, N, LDA, LDB, LDC, LDD, LDE, LDF, LDAS, LDBS, LDCS, LDDS, LDQ, LDZ, LDU, LDV, LDWORK, INFO
    INTEGER     :: LDR, LDL
    INTEGER :: MAXIT
    REAL :: TAU
    REAL :: A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), AS(LDAS, *), BS(LDBS,*), CS(LDCS, *), DS(LDDS, *)
    REAL :: Q(LDQ, *), Z(LDZ, *), U(LDU, *), V(LDV, *), E ( LDE , * ), F( LDF, *)
    REAL :: R(LDR, *), L(LDL, *)
    REAL :: CONVLOG(*), WORK(*)

    ! Locals
    INTEGER :: IT, LL
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: MONE = -1.0

    REAL :: NRMER,NRMFR, SCALE
    REAL :: NRMRHSE, NRMRHSF, NRMA, NRMB, NRMC, NRMD, EPS , TOL
    INTEGER :: IINFO, LDWORKI, MB, NB, ISOLVER, WS, ERS, FRS
    INTEGER :: L3SOLVER, MAXMN
    LOGICAL :: LASTSTEP
    LOGICAL :: BTRANSA, BTRANSB, BGUESS

    EXTERNAL SLACPY
    EXTERNAL SGEMM
    EXTERNAL SLA_TRANSFORM_GENERAL
    EXTERNAL SLA_TGCSYLV_DAG
    EXTERNAL SLA_TGCSYLV_L3
    EXTERNAL SLA_TGCSYLV_L3_2S
    EXTERNAL SLA_TGCSYLV_L2
    EXTERNAL SLA_TGCSYLV_L2_UNOPT
    EXTERNAL SLA_TGCSYLV_L2_REORDER
    EXTERNAL SLA_TGCSYLV_RECURSIVE
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_32
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_64
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_96
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY_128
    EXTERNAL SLA_TGCSYLV_L2_LOCAL_COPY
    EXTERNAL LSAME
    EXTERNAL SLANGE
    REAL SLANGE
    LOGICAL LSAME

    EXTERNAL XERROR_HANDLER

    INFO = 0
    BTRANSA = LSAME(TRANSA, 'N')
    BTRANSB = LSAME(TRANSB, 'N')
    BGUESS  = LSAME(GUESS, 'I')
    MB = TGCSYLV_BLOCKSIZE_MB(M,N)
    NB = TGCSYLV_BLOCKSIZE_NB(M,N)
    ISOLVER = TGCSYLV_ISOLVER()
    L3SOLVER = TGCSYLV_FRONTEND_SOLVER()
    MAXMN = MAX(M,N)

    !
    ! Argument Checks
    !
    IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'T')) THEN
        INFO = -1
    ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'T')) THEN
        INFO = -2
    ELSE IF ( .NOT. BGUESS .AND. .NOT. LSAME(GUESS, 'N')) THEN
        INFO = -3
    ELSE IF ( SGN1 .NE. -ONE .AND. SGN1 .NE. ONE ) THEN
        INFO = -4
    ELSE IF ( SGN2 .NE. -ONE .AND. SGN2 .NE. ONE ) THEN
        INFO = -5
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -6
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -7
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -9
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -11
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -13
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -15
    ELSE IF ( LDR .LT. MAX(1, M)) THEN
        INFO = -17
    ELSE IF ( LDL .LT. MAX(1, M)) THEN
        INFO = -19
    ELSE IF ( LDE .LT. MAX(1, M)) THEN
        INFO = -21
    ELSE IF ( LDF .LT. MAX(1, M)) THEN
        INFO = -23
    ELSE IF ( LDAS .LT. MAX(1, M)) THEN
        INFO = -25
    ELSE IF ( LDBS .LT. MAX(1, N)) THEN
        INFO = -27
    ELSE IF ( LDCS .LT. MAX(1, M)) THEN
        INFO = -29
    ELSE IF ( LDDS .LT. MAX(1, N)) THEN
        INFO = -31
    ELSE IF ( LDQ .LT. MAX(1, M)) THEN
        INFO = -33
    ELSE IF ( LDZ .LT. MAX(1, M)) THEN
        INFO = -35
    ELSE IF ( LDU .LT. MAX(1, N)) THEN
        INFO = -37
    ELSE IF ( LDV .LT. MAX(1, N)) THEN
        INFO = -39
    ELSE IF ( MAXIT .LE. 1 .OR. MAXIT .GT. 100 ) THEN
        INFO = -40
    ELSE IF ( TAU .LE. ZERO) THEN
        INFO = -41
    ELSE IF ( MB .LT. 2 .OR. NB .LT. 2 .OR. ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 )) THEN
        INFO = -50
    ELSE IF ( L3SOLVER .LT. 1 .OR. L3SOLVER .GT. 5) THEN
        INFO = -50
    END IF

    IF (INFO .NE. 0 ) THEN
        CALL XERROR_HANDLER ('SLA_GGCSYLV_REFINE', -INFO)
        RETURN
    END IF

    !
    ! Workspace Query
    !
    LDWORKI = 2 * M * N
    IINFO = -1
    SELECT CASE(L3SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL SLA_TGCSYLV_L3(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXMN .LE. 32 ) THEN
                        CALL SLA_TGCSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 64) THEN
                        CALL SLA_TGCSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 96) THEN
                        CALL SLA_TGCSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 128) THEN
                        CALL SLA_TGCSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2,M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE
                        CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    END IF
                CASE(2)
                    IF ( M .GT. 128 .OR. N .GT. 128) THEN
                        CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE
                        CALL SLA_TGCSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    END IF
                CASE(3)
                    CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(4)
                    CALL SLA_TGCSYLV_L2(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(5)
                    CALL SLA_TGCSYLV_L2_UNOPT(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(6)
                    CALL SLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL SLA_TGCSYLV_DAG(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_2STAGE)
            CALL SLA_TGCSYLV_L3_2S(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_RECURSIVE)
           CALL SLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)

    END SELECT

    IF ( IINFO .LT. 0 ) THEN
        INFO = -43
        CALL XERROR_HANDLER("SLA_GGCSYLV_REFINE", -INFO)
        RETURN
    END IF

    LDWORKI = LDWORKI + MAX(M*N,IINFO)

    IF (LDWORK .LT. 0) THEN
        LDWORK = LDWORKI
        RETURN
    END IF
    IF ( LDWORK .LT. LDWORKI) THEN
        INFO = -43
    END IF

    IF (INFO .NE. 0 ) THEN
        CALL XERROR_HANDLER ('SLA_GGCSYLV_REFINE', -INFO)
        RETURN
    END IF

    !
    ! Quick return
    !
    IF ( M .EQ. 0 .OR. N.EQ.0 ) THEN
        RETURN
    END IF

    ! Setup Workspace
    EPS = MACHINE_PRECISION
    ERS = 1
    FRS = M*N + ERS
    WS = M*N + FRS


    ! Compute the stopping criterion
    NRMA = SLANGE("F", M, M, A, LDA, WORK(WS))
    NRMC = SLANGE("F", M, M, C, LDC, WORK(WS))
    NRMB = SLANGE("F", N, N, B, LDB, WORK(WS))
    NRMD = SLANGE("F", N, N, D, LDD, WORK(WS))
    NRMER = ONE
    NRMFR = ONE
    TOL = SQRT(REAL(M*N)) * ( NRMA + NRMB + NRMC + NRMD) * EPS * TAU

    NRMRHSE = SLANGE("F", M, N, E(1,1), LDE, WORK(WS))
    NRMRHSF = SLANGE("F", M, N, F(1,1), LDF, WORK(WS))
    IT = 1
    CONVLOG(1:MAXIT) = ZERO
    LASTSTEP = .FALSE.

    ! Main Loop
    DO WHILE ( IT .LE. MAXIT .AND. .NOT. LASTSTEP)


        ! Compute the right hand side
        IF ( IT .EQ. 1 .AND. .NOT. BGUESS ) THEN
            CALL SLACPY("All", M, N, E(1,1), LDE, WORK(ERS), M)
            CALL SLACPY("All", M, N, F(1,1), LDF, WORK(FRS), M)
        ELSE
            CALL SLACPY("All", M, N, E(1,1), LDE, WORK(ERS), M)
            CALL SLACPY("All", M, N, F(1,1), LDF, WORK(FRS), M)

            CALL SGEMM(TRANSA, "N", M, N, M, MONE, A, LDA, R, LDR, ONE, WORK(ERS), M)
            CALL SGEMM("N", TRANSB, M, N, N, SGN1*MONE, L, LDL, B, LDB, ONE, WORK(ERS), M)
            CALL SGEMM(TRANSA, "N", M, N, M, MONE, C, LDC, R, LDR, ONE, WORK(FRS), M)
            CALL SGEMM("N", TRANSB, M, N, N, SGN2*MONE, L, LDL, D, LDD, ONE, WORK(FRS), M)
        END IF

        ! Check the residual norm.
        NRMER = SLANGE("F", M, N, WORK(ERS), M, WORK(WS))
        NRMFR = SLANGE("F", M, N, WORK(FRS), M, WORK(WS))


        CONVLOG(IT) =  MAX(NRMER / NRMRHSE, NRMFR / NRMRHSF)

        IF ( (NRMER .LT. NRMRHSE * TOL) .AND. (NRMFR .LT. NRMRHSF * TOL) .AND. IT.GT.1) THEN
            LASTSTEP = .TRUE.
        END IF

        ! Transform Y
        !
        ! TRANS A      TRANS B         E <-        F <-
        !   N            N             QA**T*E*ZB  QA**T*F*ZB
        !   N            T             QA**T*E*QB  QA**T*F*QB
        !   T            N             ZA**T*E*ZB  ZA**T*F*ZB
        !   T            T             ZA**T*E*QB  ZA**T*F*QB
        IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(ERS), M, Q, LDQ, V, LDV, WORK(WS))
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(FRS), M, Q, LDQ, V, LDV, WORK(WS))

        ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(ERS), M, Q, LDQ, U, LDU, WORK(WS))
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(FRS), M, Q, LDQ, U, LDU, WORK(WS))

        ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(ERS), M, Z, LDZ, V, LDV, WORK(WS))
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(FRS), M, Z, LDZ, V, LDV, WORK(WS))
        ELSE !(T,T)
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(ERS), M, Z, LDZ, U, LDU, WORK(WS))
            CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(FRS), M, Z, LDZ, U, LDU, WORK(WS))
        END IF

        !
        ! Solve the Equation
        !
        IINFO = 0
        SELECT CASE(L3SOLVER)
            CASE (FRONTEND_SOLVER_LEVEL3)
                CALL SLA_TGCSYLV_L3(TRANSA, TRANSB, SGN1, SGN2, M, N, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, &
                    & WORK(ERS), M, WORK(FRS), M, &
                    & SCALE, WORK(WS), IINFO )
            CASE (FRONTEND_SOLVER_LEVEL2)
                SELECT CASE(ISOLVER)
                    CASE(1)
                        IF ( MAXMN .LE. 32 ) THEN
                            CALL SLA_TGCSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                        ELSE IF ( MAXMN .LE. 64) THEN
                            CALL SLA_TGCSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                        ELSE IF ( MAXMN .LE. 96) THEN
                            CALL SLA_TGCSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                        ELSE IF ( MAXMN .LE. 128) THEN
                            CALL SLA_TGCSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2,M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                        ELSE
                            CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                        END IF
                    CASE(2)
                        IF ( M .GT. 128 .OR. N .GT. 128) THEN
                            CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                        ELSE
                            CALL SLA_TGCSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                        END IF
                    CASE(3)
                        CALL SLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                    CASE(4)
                        CALL SLA_TGCSYLV_L2(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                    CASE(5)
                        CALL SLA_TGCSYLV_L2_UNOPT(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                    CASE(6)
                        CALL SLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)
                END SELECT
            CASE (FRONTEND_SOLVER_DAG)
                CALL SLA_TGCSYLV_DAG(TRANSA, TRANSB, SGN1, SGN2, M, N, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, &
                    & WORK(ERS), M, WORK(FRS), M, &
                    & SCALE, WORK(WS), IINFO )

            CASE (FRONTEND_SOLVER_2STAGE)
                CALL SLA_TGCSYLV_L3_2S(TRANSA, TRANSB, SGN1, SGN2, M, N, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, &
                    & WORK(ERS), M, WORK(FRS), M, &
                    & SCALE, WORK(WS), IINFO )

            CASE (FRONTEND_SOLVER_RECURSIVE)
               CALL SLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  AS, LDAS, BS, LDBS, &
                    & CS, LDCS, DS, LDDS,  WORK(ERS), M, WORK(FRS), M, SCALE, WORK(WS), IINFO)

        END SELECT


        IF ( IINFO .NE. 0 ) THEN
            INFO = IT
            EXIT
        END IF

        ! Transform X
        !
        ! TRANS A      TRANS B         R <-        L <-
        !   N            N             ZA*R*ZB**T  QA*L*QB**T
        !   N            T             ZA*R*QB**T  QA*L*ZB**T
        !   T            N             QA*R*ZB**T  ZA*L*QB**T
        !   T            T             QA*R*QB**T  ZA*L*ZB**T
        IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(ERS), M, Z, LDZ, V, LDV, WORK(WS) )
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(FRS), M, Q, LDQ, U, LDU, WORK(WS) )
        ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(ERS), M, Z, LDZ, U, LDU, WORK(WS) )
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(FRS), M, Q, LDQ, V, LDV, WORK(WS) )

        ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(ERS), M, Q, LDQ, V, LDV, WORK(WS) )
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(FRS), M, Z, LDZ, U, LDU, WORK(WS) )
        ELSE !(T,T)
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(ERS), M, Q, LDQ, U, LDU, WORK(WS) )
            CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(FRS), M, Z, LDZ, V, LDV, WORK(WS) )
        END IF


        !
        ! Update January 2021
        !
        IF ( IT .GT. 1 ) THEN
            !$omp parallel do
            DO LL = 1, N
                R(1:M, LL) = R(1:M, LL) + WORK(ERS+(LL-1)*M:ERS+LL*M-1)
            END DO
            !$omp end parallel do
            !$omp parallel do
            DO LL = 1, N
                L(1:M, LL) = L(1:M, LL) + WORK(FRS+(LL-1)*M:FRS+LL*M-1)
            END DO
            !$omp end parallel do

        ELSE
            CALL SLACPY("All", M, N, WORK(ERS), M, R, LDR)
            CALL SLACPY("All", M, N, WORK(FRS), M, L, LDL)
        END IF

        IT = IT + 1

    END DO

    IF (LDWORK .GT. MAXIT) THEN
        WORK(1:MAXIT) = CONVLOG(1:MAXIT)
    END IF

    IF ( IT .GT. MAXIT ) THEN
        INFO = MAXIT+1
    END IF
    MAXIT = IT-1
    TAU = MAX(NRMER/NRMRHSE, NRMFR / NRMRHSF)


    RETURN
END SUBROUTINE SLA_GGCSYLV_REFINE

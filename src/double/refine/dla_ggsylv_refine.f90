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

!> \brief Iterative Refinement for the Generalized Sylvester Equations.
!
!  Definition:
!  ===========
!
!   SUBROUTINE DLA_GGSYLV_REFINE(TRANSA, TRANSB, GUESS, SGN, M , N,  &
!        & A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, Y, LDY, &
!        & AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,  &
!        & MAXIT, TAU, CONVLOG, &
!        & WORK, LDWORK, INFO)
!
!       CHARACTER, DIMENSION(1) :: TRANSA, TRANSB, GUESS
!       DOUBLE PRECISION :: SGN, TAU
!       INTEGER     :: M, N, LDA, LDB, LDC, LDD, LDX, LDY, LDAS, LDBS, LDCS, LDDS, LDQ, LDZ, LDU, LDV, LDWORK, INFO
!       INTEGER :: MAXIT
!       DOUBLE PRECISION :: A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), AS(LDAS, *), BS(LDBS,*), CS(LDCS, *), DS(LDDS, *)
!       DOUBLE PRECISION :: Q(LDQ, *), Z(LDZ, *), U(LDU, *), V(LDV, *), X ( LDX , * ), Y (LDY, *)
!       DOUBLE PRECISION :: CONVLOG(*), WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DLA_GGSYLV_REFINE solves a generalized Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) + SGN * op1(C) * X * op2(D) = Y                    (1)
!>
!> with iterative refinement, Thereby  (A,C) is a M-by-M matrix pencil and
!> (B,D) is a N-by-N matrix pencil.
!> The right hand side Y and the solution X are M-by-N matrices.
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
!>          Specifies whether X contains an initial guess or nor not.
!>          =  'I': X contains an initial guess
!>          =  'N': No initial guess, X is set to zero at the begin of the iteration.
!> \endverbatim
!>
!>
!> \param[in] SGN
!> \verbatim
!>          SGN is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between in the first equation.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          The array A contains the original matrix A defining the eqaution.
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
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The array B contains the original matrix B defining the eqaution.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,M)
!>          The array C contains the original matrix C defining the eqaution.
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
!>          D is DOUBLE PRECISION array, dimension (LDD,N)
!>          The array D contains the original matrix D defining the eqaution.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D.  LDD >= max(1,N).
!> \endverbatim
!
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,N)
!>          On input, the array X contains the initial guess, if GUESS = 'I'.
!>          On output, the array X contains the solution X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,M).
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension (LDY,N)
!>          On input, the array Y contains the right hand side Y.
!>          The array stays unchanged during the iteration.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of the array Y.  LDY >= max(1,M).
!> \endverbatim
!>
!> \param[in] AS
!> \verbatim
!>          AS is DOUBLE PRECISION array, dimension (LDAS,M)
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
!>          BS is DOUBLE PRECISION array, dimension (LDBS,N)
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
!>          CS is DOUBLE PRECISION array, dimension (LDCS,M)
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
!>          DS is DOUBLE PRECISION array, dimension (LDDS,N)
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
!>          Q is DOUBLE PRECISION array, dimension (LDQ,M)
!>          The array Q contains the left generalized Schur vectors for (A,C) as returned by DGGES.
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
!>          Z is DOUBLE PRECISION array, dimension (LDZ,M)
!>          The array Z contains the right generalized Schur vectors for (A,C) as returned by DGGES.
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
!>          U is DOUBLE PRECISION array, dimension (LDU,N)
!>          The array U contains the left generalized Schur vectors for (B,D) as returned by DGGES.
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
!>          V is DOUBLE PRECISION array, dimension (LDV,N)
!>          The array V contains the right generalized Schur vectors for (B,D) as returned by DGGES.
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
!>          TAU is DOUBLE PRECISION
!>          On input, TAU contains the additional security factor for the stopping criterion, typical values are 0.1
!>          On exit, TAU contains the last relative residual when the stopping criterion got valid.
!> \endverbatim
!
!> \param[in,out] CONVLOG
!> \verbatim
!>          CONVLOG is DOUBLE PRECISION array, dimension (MAXIT)
!>          The CONVLOG array contains the convergence history of the iterative refinement. CONVLOG(I) contains the maximum
!>          relative residual before it is solved for the I-th time.
!> \endverbatim
!
!>
!> \param[in,out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LDWORK))
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
!> \see DLA_TGSYLV_DAG
!> \see DLA_TGSYLV_LEVEL3
!> \see DLA_TGSYLV_L3_2S
!> \see DLA_TGSYLV_L2_UNOPT
!> \see DLA_TGSYLV_L2
!> \see DLA_TGSYLV_L2_REORDER
!> \see DLA_TGSYLV_L2_LOCAL_COPY_32
!> \see DLA_TGSYLV_L2_LOCAL_COPY_64
!> \see DLA_TGSYLV_L2_LOCAL_COPY_96
!> \see DLA_TGSYLV_L2_LOCAL_COPY_128
!> \see DLA_TGSYLV_L2_LOCAL_COPY
!> \see DLA_TGSYLV_GARDINER_LAUB
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Januar 2023
!> \ingroup dblggsylv
!
SUBROUTINE DLA_GGSYLV_REFINE(TRANSA, TRANSB, GUESS, SGN, M , N,  &
        & A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, Y, LDY, &
        & AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, Q, LDQ, Z, LDZ, U, LDU, V, LDV,  &
        & MAXIT, TAU, CONVLOG, &
        & WORK, LDWORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_FRONTEND_SOLVER
    IMPLICIT NONE

    ! Scalar
    CHARACTER, DIMENSION(1):: TRANSA, TRANSB, GUESS
    DOUBLE PRECISION :: SGN
    INTEGER     :: M, N, LDA, LDB, LDC, LDD, LDX, LDY, LDAS, LDBS, LDCS, LDDS, LDQ, LDZ, LDU, LDV, LDWORK, INFO
    INTEGER :: MAXIT
    DOUBLE PRECISION :: TAU
    DOUBLE PRECISION :: A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), AS(LDAS, *), BS(LDBS,*), CS(LDCS, *), DS(LDDS, *)
    DOUBLE PRECISION :: Q(LDQ, *), Z(LDZ, *), U(LDU, *), V(LDV, *), X ( LDX , * ), Y (LDY, *)
    DOUBLE PRECISION :: CONVLOG(*), WORK(*)

    ! Locals



    INTEGER :: IT, L
    DOUBLE PRECISION, PARAMETER :: ONE  = 1.0D0
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
    DOUBLE PRECISION, PARAMETER :: MONE = -1.0D0
    DOUBLE PRECISION :: NRMR, SCALE
    DOUBLE PRECISION :: NRMRHS, NRMA, NRMB, NRMC, NRMD, EPS , TOL
    INTEGER :: IINFO, LDWORKI, MB, NB, ISOLVER, WS, RS
    INTEGER :: L3SOLVER, MAXMN
    LOGICAL :: LASTSTEP
    LOGICAL :: BTRANSA, BTRANSB, BGUESS

    EXTERNAL DLACPY
    EXTERNAL DGEMM
    EXTERNAL DLA_TRANSFORM_GENERAL
    EXTERNAL DLA_TGSYLV_DAG
    EXTERNAL DLA_TGSYLV_L3
    EXTERNAL DLA_TGSYLV_L3_2S
    EXTERNAL DLA_TGSYLV_L2
    EXTERNAL DLA_TGSYLV_L2_UNOPT
    EXTERNAL DLA_TGSYLV_L2_REORDER
    EXTERNAL DLA_TGSYLV_RECURSIVE
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGSYLV_L2_LOCAL_COPY
    EXTERNAL DLA_TGSYLV_GARDINER_LAUB
    EXTERNAL LSAME
    EXTERNAL DLANGE
    DOUBLE PRECISION DLANGE
    LOGICAL LSAME

    EXTERNAL XERROR_HANDLER

    INFO = 0
    BTRANSA = LSAME(TRANSA, 'N')
    BTRANSB = LSAME(TRANSB, 'N')
    BGUESS  = LSAME(GUESS, 'I')
    MB = TGSYLV_BLOCKSIZE_MB(M,N)
    NB = TGSYLV_BLOCKSIZE_NB(M,N)
    ISOLVER = TGSYLV_ISOLVER()
    L3SOLVER = TGSYLV_FRONTEND_SOLVER()
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
    ELSE IF ( SGN .NE. -ONE .AND. SGN .NE. ONE ) THEN
        INFO = -4
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -5
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -6
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -8
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -10
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -12
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -14
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -16
    ELSE IF ( LDY .LT. MAX(1, M)) THEN
        INFO = -18
    ELSE IF ( LDAS .LT. MAX(1, M)) THEN
        INFO = -20
    ELSE IF ( LDBS .LT. MAX(1, N)) THEN
        INFO = -22
    ELSE IF ( LDCS .LT. MAX(1, M)) THEN
        INFO = -24
    ELSE IF ( LDDS .LT. MAX(1, N)) THEN
        INFO = -26
    ELSE IF ( LDQ .LT. MAX(1, M)) THEN
        INFO = -28
    ELSE IF ( LDZ .LT. MAX(1, M)) THEN
        INFO = -30
    ELSE IF ( LDU .LT. MAX(1, N)) THEN
        INFO = -32
    ELSE IF ( LDV .LT. MAX(1, N)) THEN
        INFO = -34
    ELSE IF ( MAXIT .LE. 1 .OR. MAXIT .GT. 100 ) THEN
        INFO = -35
    ELSE IF ( TAU .LE. ZERO) THEN
        INFO = -36
    ELSE IF ( MB .LT. 2 .OR. NB .LT. 2 .OR. ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 )) THEN
        INFO = -50
    ELSE IF ( L3SOLVER .LT. 1 .OR. L3SOLVER .GT. 6) THEN
        INFO = -50
    END IF

    IF (INFO .NE. 0 ) THEN
        CALL XERROR_HANDLER ('DLA_GGSYLV_REFINE', -INFO)
        RETURN
    END IF


    !
    ! Workspace Query
    !
    LDWORKI = M * N
    IINFO = -1
    SELECT CASE(L3SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL DLA_TGSYLV_L3(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, &
                & SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXMN .LE. 32 ) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 64) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 96) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 128) THEN
                        CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                    ELSE
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                    END IF
                CASE(2)
                    IF ( M .GT. 128 .OR. N .GT. 128) THEN
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                    ELSE
                        CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                    END IF
                CASE(3)
                    CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                CASE(4)
                    CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                CASE(5)
                    CALL DLA_TGSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
                CASE(6)
                    CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL DLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_2STAGE)
            CALL DLA_TGSYLV_L3_2S(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, C, LDC, D, LDD, X, LDX, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_RECURSIVE)
           CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_GARDINER_LAUB)
           CALL DLA_TGSYLV_GARDINER_LAUB(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  X, LDX, SCALE, WORK, IINFO)


    END SELECT

    IF ( IINFO .LT. 0 ) THEN
        INFO = -38
        CALL XERROR_HANDLER("DLA_GGSYLV_REFINE", -INFO)
        RETURN
    END IF

    LDWORKI = LDWORKI + MAX(M*N,IINFO)

    IF (LDWORK .LT. 0) THEN
        LDWORK = LDWORKI
        RETURN
    END IF
    IF ( LDWORK .LT. LDWORKI) THEN
        INFO = -38
    END IF

    IF (INFO .NE. 0 ) THEN
        CALL XERROR_HANDLER ('DLA_GGSYLV_REFINE', -INFO)
        RETURN
    END IF




    INFO = 0
    ! Quick return
    IF ( M .EQ. 0 .OR. N.EQ.0 ) THEN
        RETURN
    END IF


    EPS = MACHINE_PRECISION
    RS = 1
    WS = M*N + RS


    NRMA = DLANGE("F", M, M, A, LDA, WORK(WS))
    NRMC = DLANGE("F", M, M, C, LDC, WORK(WS))
    NRMB = DLANGE("F", N, N, B, LDB, WORK(WS))
    NRMD = DLANGE("F", N, N, D, LDD, WORK(WS))
    NRMRHS = DLANGE("F", M, N, Y(1,1), LDY, WORK(WS))
    NRMR = ONE

    TOL = DSQRT(DBLE(M*N)) * ( NRMA * NRMB + NRMC * NRMD) * EPS * TAU

    IT = 1
    CONVLOG(1:MAXIT) = ZERO
    LASTSTEP = .FALSE.

    !
    ! Main loop
    !
    DO WHILE ( IT .LE. MAXIT .AND. .NOT. LASTSTEP)
        IF ( IT .EQ. 1 .AND. .NOT. BGUESS) THEN
            CALL DLACPY("All", M, N, Y(1,1), LDY, WORK(RS), M)
        ELSE
            CALL DLACPY("All", M, N, Y(1,1), LDY, WORK(RS), M)
            CALL DGEMM(TRANSA, "N", M, N, M, ONE, A, LDA, X, LDX, ZERO, WORK(WS), M)
            CALL DGEMM("N", TRANSB, M, N, N, MONE, WORK(WS), M, B, LDB, ONE, WORK(RS), M)
            CALL DGEMM(TRANSA, "N", M, N, M, ONE, C, LDC, X, LDX, ZERO, WORK(WS), M)
            CALL DGEMM("N", TRANSB, M, N, N, SGN*MONE, WORK(WS), M , D, LDD, ONE, WORK(RS), M)
        END IF

        NRMR = DLANGE("F", M, N, WORK(RS), M, WORK(WS))

        CONVLOG(IT) =  NRMR / NRMRHS

        IF ( NRMR .LT. NRMRHS * TOL .AND. IT .GT. 1) THEN
            LASTSTEP = .TRUE.
        END IF

        ! Transform the right hand side
        !
        ! TRANS A      TRANS B         Y <-
        !   N            N             QA**T*Y*ZB
        !   N            T             QA**T*Y*QB
        !   T            N             ZA**T*Y*ZB
        !   T            T             ZA**T*Y*QB
        IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
            CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(RS), M, Q, LDQ, V, LDV, WORK(WS))
        ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
            CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(RS), M, Q, LDQ, U, LDU, WORK(WS))
        ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
            CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(RS), M, Z, LDZ, V, LDV, WORK(WS))
        ELSE !(T,T)
            CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, WORK(RS), M, Z, LDZ, U, LDU, WORK(WS))
        END IF

        !
        ! Solve the Equation
        !
        IINFO = 0
        SELECT CASE(L3SOLVER)
            CASE (FRONTEND_SOLVER_LEVEL3)
                CALL DLA_TGSYLV_L3(TRANSA, TRANSB, SGN, M, N, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, &
                    & WORK(RS), M, &
                    & SCALE, WORK(WS), IINFO )
            CASE (FRONTEND_SOLVER_LEVEL2)
                SELECT CASE(ISOLVER)
                    CASE(1)
                        IF ( MAXMN .LE. 32 ) THEN
                            CALL DLA_TGSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                        ELSE IF ( MAXMN .LE. 64) THEN
                            CALL DLA_TGSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                        ELSE IF ( MAXMN .LE. 96) THEN
                            CALL DLA_TGSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                        ELSE IF ( MAXMN .LE. 128) THEN
                            CALL DLA_TGSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN,M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                        ELSE
                            CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                        END IF
                    CASE(2)
                        IF ( M .GT. 128 .OR. N .GT. 128) THEN
                            CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                        ELSE
                            CALL DLA_TGSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                                & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                        END IF
                    CASE(3)
                        CALL DLA_TGSYLV_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                    CASE(4)
                        CALL DLA_TGSYLV_L2(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                    CASE(5)
                        CALL DLA_TGSYLV_L2_UNOPT(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                    CASE(6)
                        CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                            & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
                END SELECT
            CASE (FRONTEND_SOLVER_DAG)
                CALL DLA_TGSYLV_DAG(TRANSA, TRANSB, SGN, M, N, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, &
                    & WORK(RS), M, &
                    & SCALE, WORK(WS), IINFO )

            CASE (FRONTEND_SOLVER_2STAGE)
                CALL DLA_TGSYLV_L3_2S(TRANSA, TRANSB, SGN, M, N, AS, LDAS, BS, LDBS, CS, LDCS, DS, LDDS, &
                    & WORK(RS), M, &
                    & SCALE, WORK(WS), IINFO )

            CASE (FRONTEND_SOLVER_RECURSIVE)
               CALL DLA_TGSYLV_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                    & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
            CASE (FRONTEND_SOLVER_GARDINER_LAUB)
               CALL DLA_TGSYLV_GARDINER_LAUB(TRANSA, TRANSB, SGN, M, N,  AS, LDAS, BS, LDBS, &
                    & CS, LDCS, DS, LDDS,  WORK(RS), M, SCALE, WORK(WS), IINFO)
        END SELECT

        IF ( IINFO .NE. 0 ) THEN
            INFO = IT
            EXIT
        END IF

        !
        ! Back transform
        !
        IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
            CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(RS), M, Z, LDZ, U, LDU, WORK(WS))
        ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
            CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(RS), M, Z, LDZ, V, LDV, WORK(WS))
        ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
            CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(RS), M, Q, LDQ, U, LDU, WORK(WS))
        ELSE !(T,T)
            CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, WORK(RS), M, Q, LDQ, V, LDV, WORK(WS))
        END IF

        !
        ! Update January 2021
        !
        IF ( IT .GT. 1 ) THEN
            !$omp parallel do
            DO L = 1, N
                X(1:M, L) = X(1:M, L) + WORK(RS+(L-1)*M:RS+L*M-1)
            END DO
            !$omp end parallel do
        ELSE
            CALL DLACPY("All", M, N, WORK(RS), M, X, LDX)
        END IF

        IT = IT + 1

    END DO

    IF ( IT .GT. MAXIT ) THEN
        INFO = MAXIT+1
    END IF
    MAXIT = IT - 1
    TAU = NRMR/NRMRHS

    RETURN
END SUBROUTINE DLA_GGSYLV_REFINE

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

!> \brief Iterative Refinement for the Standard Stein Equations.
!
!  Definition:
!  ===========
!
!   SUBROUTINE DLA_GESTEIN_REFINE(TRANS, GUESS, M, A, LDA, X, LDX, Y, LDY, &
!                            & AS, LDAS, BS, LDBS, Q, LDQ, &
!                            & MAXIT, TAU, CONVLOG, &
!                            & WORK, LDWORK, INFO)
!           CHARACTER TRANS(1), GUESS(1)
!           INTEGER M , LDA, LDY, LDAS, LDQ, LDX, MAXIT, INFO, LDWORK
!           DOUBLE PRECISION TAU
!           DOUBLE PRECISION A(LDA,*), Q(LDQ, *), X(LDX, *), WORK(*), CONVLOG(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DLA_GESTEIN_REFINE solves a standard Stein equation of the following forms
!>
!>    A * X * A^T  -  X  = SCALE * Y                                              (1)
!>
!> or
!>
!>    A^T * X * A  -  X  = SCALE * Y                                             (2)
!>
!> where A is a M-by-M matrix using iterative refinement.
!> The right hand side Y and the solution X are M-by-M matrices.
!> The matrix A needs to be provided as the original data
!> as well as in Schur decomposition since both are required in the
!> iterative refinement process.
!>
!> \endverbatim
!>
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          Specifies the form of the system of equations with respect to A :
!>          == 'N':  Equation (1) is solved
!>          == 'T':  Equation (2) is solved
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
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and B.  M >= 0.
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
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,M)
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
!>          Y is DOUBLE PRECISION array, dimension (LDY,M)
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
!>
!> \param[in] AS
!> \verbatim
!>          AS is DOUBLE PRECISION array, dimension (LDAS,M)
!>          The array AS contains the Schur decomposition of A.
!> \endverbatim
!>
!> \param[in] LDAS
!> \verbatim
!>          LDAS is INTEGER
!>          The leading dimension of the array AS.  LDAS >= max(1,M).
!> \endverbatim
!
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,M)
!>          The array Q contains the Schur vectors for A as returned by DGEES.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,M).
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
!>
!> \see DLA_TRSTEIN_L3
!> \see DLA_TRSTEIN_L2
!> \see DLA_TRSTEIN_L3_2S
!> \see DLA_TRSTEIN_DAG
!> \see DLA_TRSTEIN_RECURSIVE
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup dblgelyap
!
SUBROUTINE DLA_GESTEIN_REFINE(TRANS, GUESS, M, A, LDA, X, LDX, Y, LDY, &
                            & AS, LDAS, Q, LDQ, &
                            & MAXIT, TAU, CONVLOG, &
                            & WORK, LDWORK, INFO)
    USE MEPACK_OPTIONS_MACHINE_DOUBLE
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_FRONTEND_SOLVER
    IMPLICIT NONE

    ! Scalar
    CHARACTER, DIMENSION(1) :: TRANS, GUESS
    INTEGER     :: M, LDA, LDX, LDY, LDAS, LDQ, LDWORK, INFO
    INTEGER :: MAXIT
    DOUBLE PRECISION :: TAU
    DOUBLE PRECISION :: A(LDA, *), AS(LDAS, *)
    DOUBLE PRECISION :: Q(LDQ, *), X( LDX, *), Y(LDY, *)
    DOUBLE PRECISION :: WORK(*), CONVLOG(*)

    ! Locals
    INTEGER :: IT, L
    DOUBLE PRECISION, PARAMETER :: ONE  = 1.0D0
    DOUBLE PRECISION, PARAMETER :: TWO  = 2.0D0
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
    DOUBLE PRECISION, PARAMETER :: MONE = -1.0D0

    DOUBLE PRECISION :: SCALE
    DOUBLE PRECISION :: NRMRHS, NRMR, NRMA, EPS , TOL
    INTEGER :: IINFO, LDWORKI, MB, ISOLVER, WS, RS
    INTEGER :: L3SOLVER
    LOGICAL :: LASTSTEP
    LOGICAL :: BTRANS, BGUESS

    EXTERNAL DLACPY
    EXTERNAL DGEMM
    EXTERNAL DLA_TRANSFORM_STANDARD
    EXTERNAL DLA_TRSTEIN_L3
    EXTERNAL DLA_TRSTEIN_DAG
    EXTERNAL DLA_TRSTEIN_L3_2S
    EXTERNAL DLA_TRSTEIN_RECURSIVE
    EXTERNAL DLA_TRSTEIN_L2
    EXTERNAL LSAME
    EXTERNAL DLANGE
    DOUBLE PRECISION DLANGE
    LOGICAL LSAME

    EXTERNAL XERROR_HANDLER

    INFO = 0
    BTRANS = LSAME(TRANS, 'T')
    BGUESS = LSAME(GUESS, 'I')
    MB = TRSTEIN_BLOCKSIZE(M)
    ISOLVER = TRSTEIN_ISOLVER()
    L3SOLVER = TRSTEIN_FRONTEND_SOLVER()

    !
    ! Argument Checks
    !
    IF ( .NOT. BTRANS .AND. .NOT. LSAME(TRANS, 'N')) THEN
        INFO = -1
    ELSE IF (.NOT. BGUESS .AND. .NOT. LSAME(GUESS, 'N')) THEN
        INFO = -2
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -3
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -5
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -7
    ELSE IF ( LDY .LT. MAX(1, M)) THEN
        INFO = -9
    ELSE IF ( LDAS .LT. MAX(1, M)) THEN
        INFO = -11
    ELSE IF ( LDQ .LT. MAX(1, M)) THEN
        INFO = -13
    ELSE IF ( MAXIT .LE. 1 .OR. MAXIT .GT. 100) THEN
        INFO = -14
    ELSE IF ( TAU .LE. ZERO) THEN
        INFO = -15
    ELSE IF ( MB .LT. 0 .OR. ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 )) THEN
        INFO = -50
    ELSE IF ( L3SOLVER .LT. 1 .OR. L3SOLVER .GT. 5) THEN
        INFO = -50
    END IF

    IF (INFO .NE. 0 ) THEN
        CALL XERROR_HANDLER ('DLA_GESTEIN_REFINE', -INFO)
        RETURN
    END IF

    !
    ! Workspace Query
    !
    LDWORKI =  M * M
    IINFO = -1
    SELECT CASE(L3SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL DLA_TRSTEIN_L3(TRANS, M, A, LDA, X, LDX, SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_LEVEL2)
            CALL DLA_TRSTEIN_L2(TRANS, M, A, LDA, X, LDX, SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_DAG)
            CALL DLA_TRSTEIN_DAG(TRANS, M, A, LDA, X, LDX, SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_2STAGE)
            CALL DLA_TRSTEIN_L3_2S(TRANS, M, A, LDA, X, LDX, SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_RECURSIVE)
            CALL DLA_TRSTEIN_RECURSIVE(TRANS, M, A, LDA, X, LDX, SCALE, WORK, IINFO )
    END SELECT

    IF ( IINFO .LT. 0 ) THEN
        INFO = -17
        CALL XERROR_HANDLER("DLA_GESTEIN_REFINE", -INFO)
        RETURN
    END IF

    LDWORKI = LDWORKI + MAX(M**2, IINFO)

    IF (LDWORK .LT. 0) THEN
        LDWORK = LDWORKI
        RETURN
    END IF
    IF ( LDWORK .LT. LDWORKI) THEN
        INFO = -17
    END IF

    IF (INFO .NE. 0 ) THEN
        CALL XERROR_HANDLER ('DLA_GESTEIN_REFINE', -INFO)
        RETURN
    END IF

    INFO = 0
    !
    ! Quick return
    !
    IF ( M .EQ. 0 ) THEN
        RETURN
    END IF

    !
    ! Workspace setup
    !
    EPS = MACHINE_PRECISION
    RS = 1
    WS = M*M + RS


    NRMA = DLANGE("F", M, M, A, LDA, WORK(WS))
    NRMRHS = DLANGE("F", M, M, Y(1,1), LDY, WORK(WS))
    NRMR = ONE

    TOL = DBLE(M) * ( NRMA**TWO ) * EPS * TAU

    IT = 1
    CONVLOG(1:MAXIT) = ZERO
    LASTSTEP = .FALSE.


    !
    ! Main Loop
    !
    DO WHILE ( IT .LE. MAXIT .AND. .NOT. LASTSTEP)
        IF ( IT .EQ. 1 .AND. .NOT. BGUESS) THEN
            CALL DLACPY("All", M, M, Y(1,1), LDY, WORK(RS), M)
        ELSE
            CALL DLACPY("All", M, M, Y(1,1), LDY, WORK(RS), M)
            IF ( BTRANS ) THEN  ! Transposed case
                CALL DGEMM("T", "N", M, M, M,  ONE, A, LDA, X, LDX, ZERO, WORK(WS), M)
                CALL DGEMM("N", "N", M, M, M, MONE, WORK(WS), M, A, LDA, ONE, WORK(RS), M)
                !$omp parallel do
                DO L = 1, M
                    WORK(RS+(L-1)*M:RS+L*M-1) = X(1:M, L) + WORK(RS+(L-1)*M:RS+L*M-1)
                END DO
                !$omp end parallel do
            ELSE
                CALL DGEMM("N", "N", M, M, M,  ONE, A, LDA, X, LDX, ZERO, WORK(WS), M)
                CALL DGEMM("N", "T", M, M, M, MONE, WORK(WS), M, A, LDA, ONE, WORK(RS), M)
                !$omp parallel do
                DO L = 1, M
                    WORK(RS+(L-1)*M:RS+L*M-1) = X(1:M, L) + WORK(RS+(L-1)*M:RS+L*M-1)
                END DO
                !$omp end parallel do
            END IF

        END IF

        NRMR = DLANGE("F", M, M, WORK(RS), M, WORK(WS))


        CONVLOG(IT) =  NRMR / NRMRHS

        IF ( IT.GT.1 .AND. (NRMR .LT. NRMRHS * TOL) .AND. IT.GT.1 ) THEN
            LASTSTEP = .TRUE.
        END IF

        ! Transform Y <- Q**T * Y * Q
        CALL DLA_TRANSFORM_STANDARD("Transpose", M, WORK(RS), M, Q, LDQ, WORK(WS))

        !
        ! Solve the Equation
        !
        IINFO = 0
        SELECT CASE(L3SOLVER)
            CASE (FRONTEND_SOLVER_LEVEL3)
                CALL DLA_TRSTEIN_L3(TRANS, M, AS, LDAS, WORK(RS), M, SCALE, WORK(WS), IINFO)
            CASE (FRONTEND_SOLVER_LEVEL2)
                CALL DLA_TRSTEIN_L2(TRANS, M, AS, LDAS, WORK(RS), M, SCALE, WORK(WS), IINFO)
            CASE (FRONTEND_SOLVER_DAG)
                CALL DLA_TRSTEIN_DAG(TRANS, M, AS, LDAS, WORK(RS), M, SCALE, WORK(WS), IINFO)
            CASE (FRONTEND_SOLVER_2STAGE)
                CALL DLA_TRSTEIN_L3_2S(TRANS, M, AS, LDAS, WORK(RS), M, SCALE, WORK(WS), IINFO)
            CASE (FRONTEND_SOLVER_RECURSIVE)
                CALL DLA_TRSTEIN_RECURSIVE(TRANS, M, AS, LDAS, WORK(RS), M, SCALE, WORK(WS), IINFO)
        END SELECT

        IF ( IINFO .NE. 0 ) THEN
            INFO = IT
            EXIT
        END IF

        ! Transform X <- Q * X * Q**T
        CALL DLA_TRANSFORM_STANDARD("NonTranspose", M, WORK(RS), M, Q, LDQ, WORK(WS))

        !
        ! Update January 2021
        !
        IF ( IT .GT. 1 ) THEN
            !$omp parallel do
            DO L = 1, M
                X(1:M, L) = X(1:M, L) + WORK(RS+(L-1)*M:RS+L*M-1)
            END DO
            !$omp end parallel do
        ELSE
            CALL DLACPY("All", M, M, WORK(RS), M, X, LDX)
        END IF
       IT = IT + 1

    END DO

    IF ( IT .GT. MAXIT ) THEN
        INFO = MAXIT+1
    END IF
    MAXIT = IT
    TAU = NRMR/NRMRHS


    RETURN
END SUBROUTINE DLA_GESTEIN_REFINE


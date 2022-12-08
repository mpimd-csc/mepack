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
! Copyright (C) Martin Koehler, 2017-2022
!

!> \brief Frontend for the solution of Standard Stein Equations.
!
!  Definition:
!  ===========
!
!   SUBROUTINE SLA_GESTEIN(FACT, TRANS, M , A, LDA, Q, LDQ, X, LDX, SCALE, WORK, LDWORK, INFO)
!           CHARACTER FACT(1), TRANS(1)
!           INTEGER M , LDA, LDQ, LDX, INFO, LDWORK
!           REAL SCALE
!           REAL A(LDA,*), Q(LDQ, *), X(LDX, *), WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SLA_GESTEIN solves a standard Stein equation of the following forms
!>
!>    A * X * A^T - X  = SCALE * Y                                              (2)
!>
!> or
!>
!>    A^T * X * A - X  =  SCALE * Y                                             (1)
!>
!> where A is a M-by-M general matrix.
!> The right hand side Y and the solution X are M-by-M matrices.  The matrix A
!> can be either supplied as general matrix or factorized in terms of its
!> Schur decomposition.
!>
!> \endverbatim
!>
!
!  Arguments:
!  ==========
!
!> \param[in] FACT
!> \verbatim
!>          FACT is CHARACTER
!>          Specifies how the matrix A is given.
!>          == 'N':  The matrix A is given as a general matrix and its Schur decomposition
!>                  A = Q*S*Q**T will be computed.
!>          == 'F':  The matrix A is given as its Schur decomposition in terms of S and Q
!>                  form A = Q*S*Q**T
!> \endverbatim
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          Specifies the form of the system of equations with respect to A:
!>          == 'N':  Equation (1) is solved.
!>          == 'T':  Equation (2) is solved.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrices A and C.  M >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,M)
!>          If FACT == "N", the matrix A is a general matrix and it is overwritten with its
!>          schur decomposition S.
!>          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S being the
!>          Schur decomposition of A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDA,M)
!>          If FACT == "N", the matrix Q is an empty M-by-M matrix on input and contains the
!>          Schur vectors of A on output.
!>          If FACT == "F", the matrix Q contains the Schur vectors of A.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,M).
!> \endverbatim
!>
!
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1) or (2)
!>          Right hand side Y and the solution X are symmetric M-by-M matrices.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
!>
!> \param[in,out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LDWORK))
!>          Workspace for the algorithm. The optmimal workspace is given either by \ref mepack_memory_frontend
!>          or a previous call to the this routine with LDWORK === -1.
!> \endverbatim
!>
!> \param[in,out] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          Size of the workspace for the algorithm. This can be determined by a call \ref mepack_memory_frontend .
!>          Alternatively, if LDWORK == -1 on input the subroutine will return the required size of the workspace in LDWORK again
!>          without performing any computations.
!> \endverbatim
!>
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          == 0:  successful exit
!>          = 1:  SGEES failed
!>          = 2:  SLA_SORT_EV failed
!>          = 3:  Inner solver failed
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!>
!
!> \see SLA_TRSTEIN_L3
!> \see SLA_TRSTEIN_L3_2S
!> \see SLA_TRSTEIN_DAG
!> \see SLA_TRSTEIN_L2
!> \see SLA_TRSTEIN_RECURSIVE
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Dezember 2022
!> \ingroup sglgelyap
!
SUBROUTINE SLA_GESTEIN(FACT, TRANS, M , A, LDA, Q, LDQ, X, LDX, SCALE, WORK, LDWORK, INFO)
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_SORTEV
    USE MEPACK_OPTIONS_FRONTEND_SOLVER
    IMPLICIT NONE

    ! Arguments
    CHARACTER FACT(1), TRANS(1)
    INTEGER M , LDA, LDQ, LDX, INFO, LDWORK
    REAL SCALE
    REAL A(LDA,*), Q(LDQ, *), X(LDX, *), WORK(*)

    ! Local Variables
    INTEGER ININFO
    INTEGER LDWORKI, IINFO
    LOGICAL BTRANS, BFACT, DUMMY
    INTEGER SDIM, MB, ISOLVER
    INTEGER SORTEV
    INTEGER SOLVER
    INTEGER BIGNB

    REAL ONE, ZERO
    PARAMETER(ONE = 1.0, ZERO = 0.0)


    ! External Functions

    EXTERNAL SLA_TRANSFORM_STANDARD
    EXTERNAL SLA_TRSTEIN_DAG
    EXTERNAL SLA_TRSTEIN_L3
    EXTERNAL SLA_TRSTEIN_L3_2S
    EXTERNAL SLA_TRSTEIN_L2
    EXTERNAL SLA_TRSTEIN_RECURSIVE
    EXTERNAL SGEMM
    EXTERNAL SGEES
    EXTERNAL SLA_SORT_EV
    EXTERNAL LSAME
    EXTERNAL XERROR_HANDLER
    LOGICAL LSAME


    ! Check Input
    BTRANS = LSAME(TRANS, 'N')
    BFACT  = LSAME(FACT, 'F')
    MB = TRSTEIN_BLOCKSIZE(M)
    ISOLVER = TRSTEIN_ISOLVER()
    SOLVER  = TRSTEIN_FRONTEND_SOLVER()
    SORTEV = SORTEV_ENABLED()
    BIGNB = TRSTEIN_BLOCKSIZE_2STAGE(M)

    ININFO = INFO
    INFO = 0
    IF ( .NOT. BFACT .AND. .NOT. LSAME(FACT, 'N')) THEN
        INFO = -1
    ELSE IF ( .NOT. BTRANS .AND. .NOT. LSAME(TRANS, 'T')) THEN
        INFO = -2
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -3
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -5
    ELSE IF ( LDQ .LT. MAX(1, M)) THEN
        INFO = -7
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -9
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    ELSE IF ( SOLVER .LT. 1 .OR. SOLVER .GT. 5 .OR. (SOLVER.EQ.FRONTEND_SOLVER_2STAGE .AND.  BIGNB .LT. MB)) THEN
        INFO = -33
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_GESTEIN', -INFO)
        RETURN
    END IF

    !
    ! Compute Workspace
    !
    LDWORKI = 1
    IF ( .NOT. BFACT) THEN
        LDWORKI = MAX(LDWORKI, 3*M)
    END IF


    IINFO = -1
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL SLA_TRSTEIN_L3(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_LEVEL2)
            CALL SLA_TRSTEIN_L2(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_DAG)
            CALL SLA_TRSTEIN_DAG(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_2STAGE)
            CALL SLA_TRSTEIN_L3_2S(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_RECURSIVE)
            CALL SLA_TRSTEIN_RECURSIVE(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
    END SELECT

    IF ( IINFO .LT. 0 ) THEN
        INFO = -12
        CALL XERROR_HANDLER("SLA_GESTEIN", -INFO)
        RETURN
    END IF
    LDWORKI = MAX(LDWORKI, IINFO, M*M)

    IF(LDWORK .EQ. -1) THEN
        LDWORK = LDWORKI
        RETURN
    ELSE IF (LDWORKI .GT. LDWORK) THEN
        INFO = -12
        CALL XERROR_HANDLER("SLA_GESTEIN", -INFO)
        RETURN
    END IF


    ! Quick Return
    INFO = 0
    IF ( M .EQ. 0 ) RETURN

    IF (.NOT. BFACT) THEN
        CALL SGEES("Vectors", "NoSort", DUMMY, M, A, LDA, SDIM, WORK(1), WORK(M+1), Q, LDQ,  &
            & WORK(2*M+1), LDWORKI-2*M, DUMMY, IINFO)
        IF ( IINFO .NE. 0 ) THEN
            INFO = 1
            RETURN
        END IF

        IF ( SORTEV .NE. 0 ) THEN
            CALL SLA_SORT_EV(M, A, LDA, Q, LDQ, MB, WORK, IINFO)
            IF ( IINFO .NE. 0 ) THEN
                INFO = 2
                RETURN
            END IF
        END IF
    END IF

    ! Transform Y <- Q**T * Y * Q
    CALL SLA_TRANSFORM_STANDARD("Transpose", M, X, LDX, Q, LDQ, WORK)

    ! Solve
    IINFO = 0
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL SLA_TRSTEIN_L3(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_LEVEL2)
            CALL SLA_TRSTEIN_L2(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_DAG)
            CALL SLA_TRSTEIN_DAG(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_2STAGE)
            CALL SLA_TRSTEIN_L3_2S(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_RECURSIVE)
            CALL SLA_TRSTEIN_RECURSIVE(TRANS, M , A, LDA, X, LDX, SCALE, WORK, IINFO)
    END SELECT

    IF ( IINFO .NE. 0) THEN
        INFO = 3
        RETURN
    END IF

    ! Transform X <- Q * X * Q**T
    CALL SLA_TRANSFORM_STANDARD("NonTranspose", M, X, LDX, Q, LDQ, WORK)

    RETURN
END SUBROUTINE


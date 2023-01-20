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

!> \brief Frontend for the solution of Standard Sylvester Equations.
!
!  Definition:
!  ===========
!
!   SUBROUTINE SLA_GESYLV2(FACTA, FACTB, TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, QA, LDQA, QB, LDQB,
!                            X, LDX, SCALE, WORK, LDWORK, INFO)
!           CHARACTER FACTA(1), FACTB(1), TRANSA(1), TRANSB(1)
!           INTEGER M , LDA, LDB, LDQA, LDQB, LDX, INFO, LDWORK
!           REAL SCALE, SGN
!           REAL A(LDA,*), B(LDB, *), QA(LDQA, *), QB(LDQB, *), X(LDX, *), WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> SLA_GESYLV2 solves a Sylvester equation of the following forms
!>
!>    op1(A) * X * op2(B) +  X  = SCALE * Y                              (1)
!>
!> or
!>
!>    op1(A) * X * op2(B) -  X  = SCALE * Y                              (2)
!>
!> where A is a M-by-M matrix and B is a N-by-N matrix. The right hand
!> side Y and the solution X are M-by-N matrices. The matrices A and B can be
!> either a general unreduced matrix or a (quasi-) upper triangular factor.
!> In the later case QA and QB provide the Schur-vectors of the matrices A
!> and B.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FACTA
!> \verbatim
!>          FACTA is CHARACTER
!>          Specifies how the matrix A is given.
!>          == 'N':  The matrix A is given as a general matrix and its Schur decomposition
!>                  A = QA*S*QA**T will be computed.
!>          == 'F':  The matrix A is given as its Schur decomposition in terms of S and QA
!>                  form A = QA*S*QA**T
!> \endverbatim
!
!> \param[in] FACTB
!> \verbatim
!>          FACTB is CHARACTER
!>          Specifies how the matrix B is given.
!>          == 'N':  The matrix B is given as a general matrix and its Schur decomposition
!>                  B = QB*R*QB**T will be computed.
!>          == 'F':  The matrix A is given as its Schur decomposition in terms of R and QB
!>                  form A = QB*R*QB**T
!> \endverbatim
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER
!>          Specifies the form of the system of equations with respect to A:
!>          == 'N':  op1(A) = A
!>          == 'T':  op1(A) = A**T
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER
!>          Specifies the form of the system of equations with respect to B:
!>          == 'N':  op2(B) = B,
!>          == 'T':  op2(B) = B**T
!> \endverbatim
!>
!> \param[in] SGN
!> \verbatim
!>          SGN is REAL, allowed values: +/-1
!>          Specifies the sign between the two parts of the Sylvester equation.
!>          = 1 :  Solve Equation (1)
!>          == -1:  Solve Equation (2)
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
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          If FACT == "N", the matrix B is a general matrix and it is overwritten with its
!>          schur decomposition R.
!>          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix R being the
!>          Schur decomposition of B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] QA
!> \verbatim
!>          QA is REAL array, dimension (LDQA,M)
!>          If FACT == "N", the matrix QA is an empty M-by-M matrix on input and contains the
!>          Schur vectors of A on output.
!>          If FACT == "F", the matrix QA contains the Schur vectors of A.
!> \endverbatim
!>
!> \param[in] LDQA
!> \verbatim
!>          LDQA is INTEGER
!>          The leading dimension of the array QA.  LDQA >= max(1,M).
!> \endverbatim
!>
!>
!> \param[in,out] QB
!> \verbatim
!>          QB is REAL array, dimension (LDQB,N)
!>          If FACT == "N", the matrix QB is an empty N-by-N matrix on input and contains the
!>          Schur vectors of B on output.
!>          If FACT == "F", the matrix QB contains the Schur vectors of B.
!> \endverbatim
!>
!> \param[in] LDQB
!> \verbatim
!>          LDQB is INTEGER
!>          The leading dimension of the array QB.  LDQB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,N)
!>          On input, the matrix X contains the right hand side Y.
!>          On output, the matrix X contains the solution of Equation (1) or (2)
!>          Right hand side Y and the solution X are M-by-N matrices.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,M).
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
!>          Alternatively, if LDWORK == -1 on input the subroutine will return the required size of the workspace in LDWORK
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
!>          = 3:  SLA_TRLYAP_DAG failed
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!>
!
!> \see SLA_TRSYLV2_L3
!> \see SLA_TRSYLV2_L3_2S
!> \see SLA_TRSYLV2_DAG
!> \see SLA_TRSYLV2_L2_UNOPT
!> \see SLA_TRSYLV2_L2
!> \see SLA_TRSYLV2_L2_REORDER
!> \see SLA_TRSYLV2_L2_LOCAL_COPY
!> \see SLA_TRSYLV2_L2_LOCAL_COPY_32
!> \see SLA_TRSYLV2_L2_LOCAL_COPY_64
!> \see SLA_TRSYLV2_L2_LOCAL_COPY_96
!> \see SLA_TRSYLV2_L2_LOCAL_COPY_128
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date Januar 2023
!> \ingroup sglgesylv
!
SUBROUTINE SLA_GESYLV2(FACTA, FACTB, TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, QA, LDQA, QB, LDQB, &
        & X, LDX, SCALE, WORK, LDWORK, INFO)
    USE MEPACK_OPTIONS_BLOCKSIZE_SINGLE
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_SORTEV
    USE MEPACK_OPTIONS_FRONTEND_SOLVER
    IMPLICIT NONE

    ! Arguments
    CHARACTER FACTA(1), FACTB(1), TRANSA(1), TRANSB(1)
    INTEGER M, N , LDA, LDB, LDQA, LDQB, LDX, INFO, LDWORK
    REAL SCALE, SGN
    REAL A(LDA,*), B(LDB, *), QA(LDQA, *), QB(LDQB, *), X(LDX, *), WORK(*)



    ! Local Variables
    INTEGER ININFO
    INTEGER LDWORKI, IINFO
    LOGICAL BTRANSA, BTRANSB
    LOGICAL BFACTA, BFACTB, DUMMY
    INTEGER SDIM, MB, ISOLVER, NB
    INTEGER SORTEV
    INTEGER BIGNB
    INTEGER SOLVER
    INTEGER MAXMN
    REAL ONE, ZERO
    PARAMETER(ONE = 1.0, ZERO = 0.0)


    ! External Functions

    EXTERNAL SLA_TRANSFORM_GENERAL
    EXTERNAL SLA_TRSYLV2_L3
    EXTERNAL SLA_TRSYLV2_L2
    EXTERNAL SLA_TRSYLV2_L2_REORDER
    EXTERNAL SLA_TRSYLV2_L2_UNOPT
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_32
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_64
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_96
    EXTERNAL SLA_TRSYLV2_L2_LOCAL_COPY_128
    EXTERNAL SLA_TRSYLV2_L3_2S
    EXTERNAL SLA_TRSYLV2_DAG
    EXTERNAL SLA_TRSYLV2_RECURSIVE

    EXTERNAL SGEMM
    EXTERNAL SGEES
    EXTERNAL SLA_SORT_EV
    EXTERNAL LSAME
    EXTERNAL XERROR_HANDLER
    LOGICAL LSAME


    ! Check Input
    BTRANSA = LSAME(TRANSA, 'N')
    BTRANSB = LSAME(TRANSB, 'N')
    BFACTA  = LSAME(FACTA, 'F')
    BFACTB  = LSAME(FACTB, 'F')

    MB = TRSYLV2_BLOCKSIZE_MB(M,N)
    NB = TRSYLV2_BLOCKSIZE_NB(M,N)
    ISOLVER = TRSYLV2_ISOLVER()
    SOLVER = TRSYLV2_FRONTEND_SOLVER()
    SORTEV = SORTEV_ENABLED()
    BIGNB = TRSYLV2_BLOCKSIZE_2STAGE(M,N)

    ININFO = INFO
    INFO = 0
    IF ( .NOT. BFACTA .AND. .NOT. LSAME(FACTA, 'N')) THEN
        INFO = -1
    ELSE IF ( .NOT. BFACTB .AND. .NOT. LSAME(FACTB, 'N')) THEN
        INFO = -2
    ELSE IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'T')) THEN
        INFO = -3
    ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'T')) THEN
        INFO = -4
    ELSE IF ( SGN .NE. -ONE .AND. SGN .NE. ONE ) THEN
        INFO = -5
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -6
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -7
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -9
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -11
    ELSE IF ( LDQA .LT. MAX(1, M)) THEN
        INFO = -13
    ELSE IF ( LDQB .LT. MAX(1, N)) THEN
        INFO = -15
    ELSE IF ( LDX .LT. MAX(1, M)) THEN
        INFO = -17
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -30
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -31
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -32
    ELSE IF ( SOLVER .LT. 1 .OR. SOLVER .GT. 5 .OR. (SOLVER .EQ. FRONTEND_SOLVER_2STAGE .AND. BIGNB .LT. MIN(MB,NB))) THEN
        INFO = -33
    END IF


    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('SLA_GESYLV2', -INFO)
        RETURN
    END IF

    MAXMN = MAX (M, N)
    !
    ! Compute Workspace
    !
    LDWORKI = 1
    IF ( .NOT. BFACTA) THEN
        LDWORKI = MAX(LDWORKI, 3*M )
    END IF
    IF ( .NOT. BFACTB) THEN
        LDWORKI = MAX(LDWORKI, 3*N )
    END IF
    IINFO = -1
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL SLA_TRSYLV2_L3(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXMN .LE. 32 ) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 64) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 96) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 128) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB,X, LDX, SCALE, WORK, IINFO)
                    ELSE
                        CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    END IF
                CASE(2)
                    IF ( M .GT. 128 .OR. N .GT. 128) THEN
                        CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    END IF
                CASE(3)
                    CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                CASE(4)
                    CALL SLA_TRSYLV2_L2(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                CASE(5)
                    CALL SLA_TRSYLV2_L2_UNOPT(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                CASE(6)
                    CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB,  X, LDX, SCALE, WORK, IINFO)

            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL SLA_TRSYLV2_DAG(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_2STAGE)
            CALL SLA_TRSYLV2_L3_2S(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_RECURSIVE)
            CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
    END SELECT

    IF ( IINFO .LT. 0 ) THEN
        INFO = -20
        CALL XERROR_HANDLER("SLA_GESYLV2", -INFO)
        RETURN
    END IF
    LDWORKI = MAX(LDWORKI, IINFO, M*N)

    IF(LDWORK .EQ. -1) THEN
        LDWORK = LDWORKI
        RETURN
    ELSE IF (LDWORKI .GT. LDWORK) THEN
        INFO = -20
        CALL XERROR_HANDLER("SLA_GESYLV2", -INFO)
        RETURN
    END IF



    !
    ! Quick Return
    !
    INFO = 0
    IF ( M .EQ. 0 ) RETURN
    IF ( N .EQ. 0 ) RETURN


    !
    ! Factorize A
    !
    IF (.NOT. BFACTA) THEN
        CALL SGEES("Vectors", "NoSort", DUMMY, M, A, LDA, SDIM, WORK(1), WORK(M+1), QA, LDQA,  &
            & WORK(2*M+1), LDWORKI-2*M, DUMMY, IINFO)
        IF ( IINFO .NE. 0 ) THEN
            INFO = 1
            RETURN
        END IF

        IF ( SORTEV .NE. 0 ) THEN
            CALL SLA_SORT_EV(M, A, LDA, QA, LDQA, MB, WORK, IINFO)
            IF ( IINFO .NE. 0 ) THEN
                INFO = 2
                RETURN
            END IF
        END IF
    END IF

    !
    ! Factorize B
    !
    IF (.NOT. BFACTB) THEN
        CALL SGEES("Vectors", "NoSort", DUMMY, N, B, LDB, SDIM, WORK(1), WORK(N+1), QB, LDQB,  &
            & WORK(2*N+1), LDWORKI-2*N, DUMMY, IINFO)
        IF ( IINFO .NE. 0 ) THEN
            INFO = 1
            RETURN
        END IF

        IF ( SORTEV .NE. 0 ) THEN
            CALL SLA_SORT_EV(N, B, LDB, QB, LDQB, NB, WORK, IINFO)
            IF ( IINFO .NE. 0 ) THEN
                INFO = 2
                RETURN
            END IF
        END IF
    END IF

    ! Transform Y <- QA**T * Y * QB
    CALL SLA_TRANSFORM_GENERAL("Transpose", M, N, X, LDX, QA, LDQA, QB, LDQB, WORK)

    ! Solve
    IINFO = 0
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL SLA_TRSYLV2_L3(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXMN .LE. 32 ) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 64) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 96) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 128) THEN
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB,X, LDX, SCALE, WORK, IINFO)
                    ELSE
                        CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    END IF
                CASE(2)
                    IF ( M .GT. 128 .OR. N .GT. 128) THEN
                        CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    ELSE
                        CALL SLA_TRSYLV2_L2_LOCAL_COPY(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                    END IF
                CASE(3)
                    CALL SLA_TRSYLV2_L2_REORDER(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                CASE(4)
                    CALL SLA_TRSYLV2_L2(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                CASE(5)
                    CALL SLA_TRSYLV2_L2_UNOPT(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
                CASE(6)
                    CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, M, N,  A, LDA, B, LDB,  X, LDX, SCALE, WORK, IINFO)

            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL SLA_TRSYLV2_DAG(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_2STAGE)
            CALL SLA_TRSYLV2_L3_2S(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
        CASE (FRONTEND_SOLVER_RECURSIVE)
            CALL SLA_TRSYLV2_RECURSIVE(TRANSA, TRANSB, SGN, M, N, A, LDA, B, LDB, X, LDX, SCALE, WORK, IINFO)
    END SELECT
    IF ( IINFO .NE. 0) THEN
        INFO = 3
        RETURN
    END IF

    ! Transform X <- QA * X * QB**T
    CALL SLA_TRANSFORM_GENERAL("NonTranspose", M, N, X, LDX, QA, LDQA, QB, LDQB, WORK)

    RETURN
END SUBROUTINE


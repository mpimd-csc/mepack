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

!> \brief Frontend for the solution of Coupled Generalized Sylvester Equations.
!
!  Definition:
!  ===========
!
!   SUBROUTINE DLA_GGCSYLV(FACTA, FACTB, TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD,
!                            QA, LDQA, ZA, LDZA, QB, LDQB, ZB, LDZB,
!                            E, LDE, F, LDF, SCALE, WORK, LDWORK, INFO)
!           CHARACTER FACTA(1), FACTB(1), TRANSA(1), TRANSB(1)
!           INTEGER M , LDA, LDB, LDC, LDD, LDQA, LDQB, LDZA, LDZB, LDE, LDF, INFO, LDWORK
!           DOUBLE PRECISION SCALE, SGN1, SGN2
!           DOUBLE PRECISION A(LDA,*), B(LDB, *), C(LDC, *), D(LDD,*)
!           DOUBLE PRECISION QA(LDQA, *), ZA(LDZA, *), QB(LDQB, *), ZB(LDZB, *), E(LDE, *), F(LDF,*),
!           DOUBLE PRECISION WORK(*)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DLA_GGCSYLV solves a coupled generalized Sylvester equation of the following forms
!>
!>    op1(A) * R  + SGN1 * L  * op2(B) = SCALE * E                              (1)
!>    op1(C) * R  + SGN2 * L  * op2(D) = SCALE * F
!>
!> where (A,C) is a M-by-M matrix pencil and (B,D) is a N-by-N matrix pencil.
!> The right hand side (E,F) and the solution (R,L) are M-by-N matrix pencils. The pencils (A,C)
!> and (B,D) can be either given as general unreduced matrices, as generalized
!> Hessenberg form, or in terms of their generalized Schur decomposition.
!> If they are given as general matrices or as a generalized Hessenberg form
!> their generalized Schur decomposition will be computed.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FACTA
!> \verbatim
!>          FACTA is CHARACTER
!>          Specifies how the matrix pencil (A,C) is given.
!>          == 'N':  The matrix pencil (A,C) is given as a general matrices and its Schur decomposition
!>                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
!>          == 'F':  The matrix pencil (A,C) is already in generalized Schur form and S, R, QA, and ZA
!>                  are given.
!>          == 'H': The matrix pencil (A,C) is given in generalized Hessenberg form and its Schur decomposition
!>                  A = QA*S*ZA**T, C = QA*R*ZA**T will be computed.
!> \endverbatim
!
!> \param[in] FACTB
!> \verbatim
!>          FACTB is CHARACTER
!>          Specifies how the matrix pencil (B,D) is given.
!>          == 'N':  The matrix pencil (B,D) is given as a general matrices and its Schur decomposition
!>                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
!>          == 'F':  The matrix pencil (B,D) is already in generalized Schur form and U, V, QB, and ZB
!>                  are given.
!>          == 'H': The matrix pencil (B,D) is given in generalized Hessenberg form and its Schur decomposition
!>                  B = QB*U*ZB**T, D = QB*V*ZB**T will be computed.
!> \endverbatim
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
!> \param[in] SGN1
!> \verbatim
!>          SGN1 is DOUBLE PRECISION, allowed values: +/-1
!>          Specifies the sign between in the first equation.
!> \endverbatim
!>
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is DOUBLE PRECISION, allowed values: +/-1
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
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          If FACT == "N", the matrix A is a general matrix and it is overwritten with the
!>          (quasi-) upper triangular factor S of the Schur decomposition of (A,C).
!>          If FACT == "F", the matrix A contains its (quasi-) upper triangular matrix S of
!>          the Schur decomposition of (A,C).
!>          If FACT == "H", the matrix A is an upper Hessenberg matrix of the generalized
!>          Hessenberg form (A,C) and it is overwritten with the (quasi-) upper triangular
!>          factor S of the Schur decomposition of (A,C).
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
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          If FACT == "N",  the matrix B is a general matrix and it is overwritten with the
!>          (quasi-) upper triangular factor U of the Schur decomposition of (B,D).
!>          If FACT == "F", the matrix B contains its (quasi-) upper triangular matrix U of
!>          the Schur decomposition of (B,D).
!>          If FACT == "H", the matrix B is an upper Hessenberg matrix of the generalized
!>          Hessenberg form (B,D) and it is overwritten with the (quasi-) upper triangular
!>          factor U of the Schur decomposition of (B,D).
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,M)
!>          If FACT == "N", the matrix C is a general matrix and it is overwritten with the
!>          upper triangular factor R of the Schur decomposition of (A,C).
!>          If FACT == "F", the matrix C contains its upper triangular matrix R of
!>          the Schur decomposition of (A,C).
!>          If FACT == "H", the matrix C is the upper triangular matrix of the generalized Hessenberg form
!>          (A,C) and it is overwritten with the upper triangular factor R of the Schur decomposition of (A,C).
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (LDD,N)
!>          If FACT == "N",  the matrix D is a general matrix and it is overwritten with the
!>          upper triangular factor V of the Schur decomposition of (B,D).
!>          If FACT == "F", the matrix D contains its upper triangular matrix V of
!>          the Schur decomposition of (B,D).
!>          If FACT == "H", the matrix D is the upper triangular matrix of the generalized Hessenberg form
!>          (B,D) and it is overwritten with the upper triangular factor V of the Schur decomposition of (B,D).
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array D.  LDD >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] QA
!> \verbatim
!>          QA is DOUBLE PRECISION array, dimension (LDQA,M)
!>          If FACT == "N", the matrix QA is an empty M-by-M matrix on input and contains the
!>          left Schur vectors of (A,C) on output.
!>          If FACT == "F", the matrix QA contains the left Schur vectors of (A,C).
!>          If FACT == "H", the matrix QA is an empty M-by-M matrix on input and contains the
!>          left Schur vectors of (A,C) on output.
!> \endverbatim
!>
!> \param[in] LDQA
!> \verbatim
!>          LDQA is INTEGER
!>          The leading dimension of the array QA.  LDQA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] ZA
!> \verbatim
!>          ZA is DOUBLE PRECISION array, dimension (LDZA,M)
!>          If FACT == "N", the matrix ZA is an empty M-by-M matrix on input and contains the
!>          right Schur vectors of (A,C) on output.
!>          If FACT == "F", the matrix ZA contains the right Schur vectors of (A,C).
!>          If FACT == "H", the matrix ZA is an empty M-by-M matrix on input and contains the
!>          right Schur vectors of (A,C) on output.
!> \endverbatim
!>
!> \param[in] LDZA
!> \verbatim
!>          LDZA is INTEGER
!>          The leading dimension of the array ZA.  LDZA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] QB
!> \verbatim
!>          QB is DOUBLE PRECISION array, dimension (LDQB,N)
!>          If FACT == "N", the matrix QB is an empty N-by-N matrix on input and contains the
!>          left Schur vectors of (B,D) on output.
!>          If FACT == "F", the matrix QB contains the left Schur vectors of (B,D).
!>          If FACT == "H", the matrix QB is an empty M-by-M matrix on input and contains the
!>          left Schur vectors of (B,D) on output.
!> \endverbatim
!>
!> \param[in] LDQB
!> \verbatim
!>          LDQB is INTEGER
!>          The leading dimension of the array QB.  LDQB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] ZB
!> \verbatim
!>          ZB is DOUBLE PRECISION array, dimension (LDZB,N)
!>          If FACT == "N", the matrix ZB is an empty N-by-N matrix on input and contains the
!>          right Schur vectors of (B,D) on output.
!>          If FACT == "F", the matrix ZB contains the right Schur vectors of (B,D).
!>          If FACT == "H", the matrix ZB is an empty M-by-M matrix on input and contains the
!>          right Schur vectors of (B,D) on output.
!> \endverbatim
!>
!> \param[in] LDZB
!> \verbatim
!>          LDZB is INTEGER
!>          The leading dimension of the array ZB.  LDZB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (LDE,N)
!>          On input, the matrix E contains the right hand side E.
!>          On output, the matrix E contains the solution R.
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
!>          F is DOUBLE PRECISION array, dimension (LDF,N)
!>          On input, the matrix F contains the right hand side F.
!>          On output, the matrix F contains the solution L.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F.  LDF >= max(1,M).
!> \endverbatim
!
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          SCALE is a scaling factor to prevent the overflow in the result.
!>          If INFO == 0 then SCALE is 1.0D0 otherwise if one of the inner systems
!>          could not be solved correctly, 0 < SCALE <= 1 holds true.
!> \endverbatim
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
!>          If LDWORK == -1 on input the subroutine will return the required size of the workspace in LDWORK on exit.
!>          No computations are performed and none of the arrays are referenced.
!> \endverbatim
!>
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          == 0:  successful exit
!>          = 1:  DHGGES failed
!>          = 2:  DLA_SORT_GEV failed
!>          = 3:  Inner solver failed
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!>
!> \see DLA_TGCSYLV_DAG
!> \see DLA_TGCSYLV_LEVEL3
!> \see DLA_TGCSYLV_L3_2S
!> \see DLA_TGCSYLV_L2_UNOPT
!> \see DLA_TGCSYLV_L2
!> \see DLA_TGCSYLV_L2_REORDER
!> \see DLA_TGCSYLV_L2_LOCAL_COPY_32
!> \see DLA_TGCSYLV_L2_LOCAL_COPY_64
!> \see DLA_TGCSYLV_L2_LOCAL_COPY_96
!> \see DLA_TGCSYLV_L2_LOCAL_COPY_128
!> \see DLA_TGCSYLV_L2_LOCAL_COPY
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date October 2023
!> \ingroup dblggsylv
!
SUBROUTINE DLA_GGCSYLV(FACTA, FACTB, TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, &
        & QA, LDQA, ZA, LDZA, QB, LDQB, ZB, LDZB, &
        & E, LDE, F, LDF, SCALE, WORK, LDWORK, INFO)
    USE MEPACK_OPTIONS_BLOCKSIZE_DOUBLE
    USE MEPACK_OPTIONS_ISOLVER
    USE MEPACK_OPTIONS_SORTEV
    USE MEPACK_OPTIONS_FRONTEND_SOLVER
    IMPLICIT NONE

    ! Arguments
    CHARACTER FACTA(1), FACTB(1), TRANSA(1), TRANSB(1)
    INTEGER M, N , LDA, LDB, LDC, LDD, LDQA, LDQB, LDZA, LDZB, LDE, LDF, INFO, LDWORK
    DOUBLE PRECISION SCALE, SGN1, SGN2
    DOUBLE PRECISION A(LDA,*), B(LDB, *), C(LDC, *), D(LDD,*)
    DOUBLE PRECISION QA(LDQA, *), QB(LDQB, *), ZA(LDZA, *), ZB(LDZB, *), E(LDE, *),F(LDF,*),  WORK(*)



    ! Local Variables
    INTEGER ININFO
    INTEGER LDWORKI, IINFO
    LOGICAL BTRANSA, BTRANSB
    LOGICAL BFACTA, BFACTB, DUMMY, AHESS, BHESS
    INTEGER SDIM, MB, ISOLVER, NB, BIGNB
    INTEGER SORTEV
    INTEGER SOLVER
    INTEGER MAXMN

    DOUBLE PRECISION ONE, ZERO
    PARAMETER(ONE = 1.0D0, ZERO = 0.0D0)
    CHARACTER ASHAPE, BSHAPE


    ! External Functions

    EXTERNAL DLA_TRANSFORM_GENERAL
    EXTERNAL DLA_TGCSYLV_DAG
    EXTERNAL DLA_TGCSYLV_L3
    EXTERNAL DLA_TGCSYLV_L3_2S
    EXTERNAL DLA_TGCSYLV_L2
    EXTERNAL DLA_TGCSYLV_L2_UNOPT
    EXTERNAL DLA_TGCSYLV_L2_REORDER
    EXTERNAL DLA_TGCSYLV_RECURSIVE
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_32
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_64
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_96
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY_128
    EXTERNAL DLA_TGCSYLV_L2_LOCAL_COPY

    EXTERNAL DHGGES
    EXTERNAL DLA_SORT_GEV
    EXTERNAL LSAME
    EXTERNAL XERROR_HANDLER
    LOGICAL LSAME


    ! Check Input
    BTRANSA = LSAME(TRANSA, 'N')
    BTRANSB = LSAME(TRANSB, 'N')
    BFACTA  = LSAME(FACTA, 'F')
    BFACTB  = LSAME(FACTB, 'F')
    AHESS  = LSAME(FACTA, 'H')
    BHESS  = LSAME(FACTB, 'H')

    MB = TGCSYLV_BLOCKSIZE_MB(M,N)
    NB = TGCSYLV_BLOCKSIZE_NB(M,N)
    ISOLVER = TGCSYLV_ISOLVER()
    SOLVER  = TGCSYLV_FRONTEND_SOLVER()
    SORTEV = SORTEV_ENABLED()
    MAXMN = MAX(M, N)
    BIGNB = TGCSYLV_BLOCKSIZE_2STAGE(M,N)

    ININFO = INFO
    INFO = 0
    IF ( .NOT. BFACTA .AND. .NOT. AHESS .AND. .NOT. LSAME(FACTA, 'N')) THEN
        INFO = -1
    ELSE IF ( .NOT. BFACTB .AND. .NOT. BHESS .AND. .NOT. LSAME(FACTB, 'N')) THEN
        INFO = -2
    ELSE IF ( .NOT. BTRANSA .AND. .NOT. LSAME(TRANSA, 'T')) THEN
        INFO = -3
    ELSE IF ( .NOT. BTRANSB .AND. .NOT. LSAME(TRANSB, 'T')) THEN
        INFO = -4
    ELSE IF ( SGN1 .NE. -ONE .AND. SGN1 .NE. ONE ) THEN
        INFO = -5
    ELSE IF ( SGN2 .NE. -ONE .AND. SGN2 .NE. ONE ) THEN
        INFO = -6
    ELSE IF ( M .LT. 0 ) THEN
        INFO = -7
    ELSE IF ( N .LT. 0 ) THEN
        INFO = -8
    ELSE IF ( LDA .LT. MAX(1, M)) THEN
        INFO = -10
    ELSE IF ( LDB .LT. MAX(1, N)) THEN
        INFO = -12
    ELSE IF ( LDC .LT. MAX(1, M)) THEN
        INFO = -14
    ELSE IF ( LDD .LT. MAX(1, N)) THEN
        INFO = -16
    ELSE IF ( LDQA .LT. MAX(1, M)) THEN
        INFO = -18
    ELSE IF ( LDZA .LT. MAX(1, M)) THEN
        INFO = -20
    ELSE IF ( LDQB .LT. MAX(1, N)) THEN
        INFO = -22
    ELSE IF ( LDZB .LT. MAX(1, N)) THEN
        INFO = -24
    ELSE IF ( LDE .LT. MAX(1, M)) THEN
        INFO = -26
    ELSE IF ( LDF .LT. MAX(1, M)) THEN
        INFO = -28
    ELSE IF ( MB .LT. 2 ) THEN
        INFO = -40
    ELSE IF ( NB .LT. 2 ) THEN
        INFO = -41
    ELSE IF ( ISOLVER .LT. 1 .OR. ISOLVER .GT. 6 ) THEN
        INFO = -42
    ELSE IF ( SOLVER .LT. 1 .OR. SOLVER .GT. 5 .OR. (SOLVER .EQ. FRONTEND_SOLVER_2STAGE .AND. BIGNB .LT. MIN(MB,NB)) ) THEN
        INFO = -43
    END IF

    IF ( INFO .LT. 0 ) THEN
        CALL XERROR_HANDLER('DLA_GGCSYLV', -INFO)
        RETURN
    END IF

    !
    ! Compute the required workspace
    !
    LDWORKI = 1
    IF ( .NOT. BFACTA) THEN
        LDWORKI = MAX(LDWORKI, 3*M + 8*M+16)
    END IF
    IF ( .NOT. BFACTB) THEN
        LDWORKI = MAX(LDWORKI, 3*N + 8*N+16)
    END IF
    IINFO = -1
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL DLA_TGCSYLV_L3(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXMN .LE. 32 ) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 64) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 96) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 128) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2,M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE
                        CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    END IF
                CASE(2)
                    IF ( M .GT. 128 .OR. N .GT. 128) THEN
                        CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    END IF
                CASE(3)
                    CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(4)
                    CALL DLA_TGCSYLV_L2(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(5)
                    CALL DLA_TGCSYLV_L2_UNOPT(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(6)
                    CALL DLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL DLA_TGCSYLV_DAG(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_2STAGE)
            CALL DLA_TGCSYLV_L3_2S(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_RECURSIVE)
           CALL DLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)

    END SELECT

    IF ( IINFO .LT. 0 ) THEN
        INFO = -31
        CALL XERROR_HANDLER("DLA_GGCSYLV", -INFO)
        RETURN
    END IF
    LDWORKI = MAX(LDWORKI, IINFO, M*N)

    IF(LDWORK .EQ. -1) THEN
        LDWORK = LDWORKI
        RETURN
    ELSE IF (LDWORKI .GT. LDWORK) THEN
        INFO = -31
        CALL XERROR_HANDLER("DLA_GGCSYLV", -INFO)
        RETURN
    END IF


    INFO = 0
    IF ( M .EQ. 0 ) RETURN
    IF ( N .EQ. 0 ) RETURN
    !
    ! Factorize A
    !
    IF (.NOT. BFACTA) THEN
        ASHAPE = 'G'
        IF ( AHESS ) THEN
            ASHAPE = 'H'
        END IF
        CALL DHGGES(ASHAPE, "Vectors", "Vectors", "NoSort", DUMMY, M, A, LDA, C, LDC, SDIM, WORK(1), WORK(M+1), WORK(2*M+1), &
            & QA, LDQA,  ZA, LDZA,  WORK(3*M+1), LDWORKI-3*M, DUMMY, IINFO)
        IF ( IINFO .NE. 0 ) THEN
            INFO = 1
            RETURN
        END IF

        IF ( SORTEV .NE. 0 ) THEN
            CALL DLA_SORT_GEV(M, A, LDA, C, LDC, QA, LDQA, ZA, LDZA, MB, WORK, LDWORKI, IINFO)
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
        BSHAPE = 'G'
        IF ( BHESS ) THEN
            BSHAPE = 'H'
        END IF
         CALL DHGGES(BSHAPE, "Vectors", "Vectors", "NoSort", DUMMY, N, B, LDB, D, LDD, SDIM, WORK(1), WORK(N+1), WORK(2*N+1), &
            & QB, LDQB,  ZB, LDZB,  WORK(3*N+1), LDWORKI-3*M, DUMMY, IINFO)
        IF ( IINFO .NE. 0 ) THEN
            INFO = 1
            RETURN
        END IF

        IF ( SORTEV .NE. 0 ) THEN
            CALL DLA_SORT_GEV(N, B, LDB, D, LDD, QB, LDQB, ZB, LDZB, NB, WORK, LDWORKI, IINFO)
            IF ( IINFO .NE. 0 ) THEN
                INFO = 2
                RETURN
            END IF
        END IF

    END IF

    ! Transform Y
    !
    ! TRANS A      TRANS B         E <-        F <-
    !   N            N             QA**T*E*ZB  QA**T*F*ZB
    !   N            T             QA**T*E*QB  QA**T*F*QB
    !   T            N             ZA**T*E*ZB  ZA**T*F*ZB
    !   T            T             ZA**T*E*QB  ZA**T*F*QB
    IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, E, LDE, QA, LDQA, ZB, LDZB, WORK)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, F, LDF, QA, LDQA, ZB, LDZB, WORK)

    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, E, LDE, QA, LDQA, QB, LDQB, WORK)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, F, LDF, QA, LDQA, QB, LDQB, WORK)

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, E, LDE, ZA, LDZA, ZB, LDZB, WORK)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, F, LDF, ZA, LDZA, ZB, LDZB, WORK)
    ELSE !(T,T)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, E, LDE, ZA, LDZA, QB, LDQB, WORK)
        CALL DLA_TRANSFORM_GENERAL("Transpose", M, N, F, LDF, ZA, LDZA, QB, LDQB, WORK)
    END IF

    !
    ! Solve the triangular equation
    !
    IINFO = 0
    SELECT CASE(SOLVER)
        CASE (FRONTEND_SOLVER_LEVEL3)
            CALL DLA_TGCSYLV_L3(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )
        CASE (FRONTEND_SOLVER_LEVEL2)
            SELECT CASE(ISOLVER)
                CASE(1)
                    IF ( MAXMN .LE. 32 ) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_32(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 64) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_64(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 96) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_96(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE IF ( MAXMN .LE. 128) THEN
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY_128(TRANSA, TRANSB, SGN1, SGN2,M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE
                        CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    END IF
                CASE(2)
                    IF ( M .GT. 128 .OR. N .GT. 128) THEN
                        CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    ELSE
                        CALL DLA_TGCSYLV_L2_LOCAL_COPY(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                            & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                    END IF
                CASE(3)
                    CALL DLA_TGCSYLV_L2_REORDER(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(4)
                    CALL DLA_TGCSYLV_L2(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(5)
                    CALL DLA_TGCSYLV_L2_UNOPT(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
                CASE(6)
                    CALL DLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                        & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)
            END SELECT
        CASE (FRONTEND_SOLVER_DAG)
            CALL DLA_TGCSYLV_DAG(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_2STAGE)
            CALL DLA_TGCSYLV_L3_2S(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, &
                & SCALE, WORK, IINFO )

        CASE (FRONTEND_SOLVER_RECURSIVE)
           CALL DLA_TGCSYLV_RECURSIVE(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, &
                & C, LDC, D, LDD,  E, LDE, F, LDF, SCALE, WORK, IINFO)

    END SELECT

    IF ( IINFO .NE. 0) THEN
        INFO = 3
        RETURN
    END IF

    ! Transform X
    !
    ! TRANS A      TRANS B         R <-        L <-
    !   N            N             ZA*R*ZB**T  QA*L*QB**T
    !   N            T             ZA*R*QB**T  QA*L*ZB**T
    !   T            N             QA*R*ZB**T  ZA*L*QB**T
    !   T            T             QA*R*QB**T  ZA*L*ZB**T
    IF ( BTRANSA .AND. BTRANSB ) THEN   !(N,N)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, E, LDE, ZA, LDZA, ZB, LDZB, WORK)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, F, LDF, QA, LDQA, QB, LDQB, WORK)
    ELSE IF ( BTRANSA .AND. .NOT. BTRANSB ) THEN   !(N,T)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, E, LDE, ZA, LDZA, QB, LDQB, WORK)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, F, LDF, QA, LDQA, ZB, LDZB, WORK)

    ELSE IF ( .NOT. BTRANSA .AND. BTRANSB ) THEN   !(T,N)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, E, LDE, QA, LDQA, ZB, LDZB, WORK)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, F, LDF, ZA, LDZA, QB, LDQB, WORK)
    ELSE !(T,T)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, E, LDE, QA, LDQA, QB, LDQB, WORK)
        CALL DLA_TRANSFORM_GENERAL("NonTranspose", M, N, F, LDF, ZA, LDZA, ZB, LDZB, WORK)
    END IF

    RETURN
END SUBROUTINE


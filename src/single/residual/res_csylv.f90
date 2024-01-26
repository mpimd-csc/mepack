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

!> \brief Compute the relative residual for the coupled Sylvester equation.
!
!  Definition:
!  ===========
!       FUNCTION DLA_RESIDUAL CSYLV(TRANSA, TRANSB, SGN1, SGN2, M, N,  A, LDA, B, LDB, C, LDC, D, LDD,
!                          R, LDR, L, LDL, E, LDE, F, LDF, SCALE)
!           CHARACTER TRANSA(1), TRANSB(1)
!           REAL SCALE, SGN1, SGN2
!           INTEGER M, N, LDA, LDB. LDC, LDD, LDR, LDL, LDE, LDF
!           REAL A(LDA, *), B(LDB, *), C(LDC, *), D(LDD, *), R(LDR, *), L(LDL, *), E(LDE, *), F(LDF, *)
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Compute the relative residual of the coupled Sylvester equation
!>
!>  RelRes = max( || SCALE*E - opA(A)*R - SGN1 * L *opB(B) || / ( (||A||*||R||+||B||*||L||) + SCALE * ||E|| ),
!>                || SCALE*F - opA(C)*R - SGN2 * L *opB(D) || / ( (||C||*||R||+||D||*||L||) + SCALE * ||F|| ))    (1)
!>
!> The norms are evaluated in terms of the Frobenius norm.
!> \endverbatim
!> \remark Auxiliary memory is managed by the function itself.
!>
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to A :
!>          == 'N':  opA(A) = A (No transpose for A)
!>          == 'T':  opA(A) = A**T (Transpose A)
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER(1)
!>          Specifies the form of the system of equations with respect to B :
!>          == 'N':  opB(B) = B (No transpose for B)
!>          == 'T':  opB(B) = B**T (Transpose B)
!> \endverbatim
!
!> \param[in] SGN1
!> \verbatim
!>          SGN1 is REAL, allowed values: +/-1
!>          Specifies the sign between the two parts of the first equation.
!> \endverbatim
!
!> \param[in] SGN2
!> \verbatim
!>          SGN2 is REAL, allowed values: +/-1
!>          Specifies the sign between the two parts of the second equation.
!> \endverbatim
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The order of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,M)
!>          The coefficient matrix A.
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
!>          The coefficient matrix B.
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
!>          The coefficient matrix C.
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
!>          The coefficient matrix D.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of the array A.  LDD >= max(1,N).
!> \endverbatim
!
!> \param[in] R
!> \verbatim
!>          R is REAL array, dimension (LDR,N)
!>          First part of the solution of the Sylvester Equation.
!> \endverbatim
!>
!> \param[in] LDR
!> \verbatim
!>          LDR is INTEGER
!>          The leading dimension of the array R.  LDR >= max(1,M).
!> \endverbatim
!
!> \param[in] L
!> \verbatim
!>          L is REAL array, dimension (LDL,N)
!>          Second part of the solution of the Sylvester Equation.
!> \endverbatim
!>
!> \param[in] LDL
!> \verbatim
!>          LDL is INTEGER
!>          The leading dimension of the array L.  LDL >= max(1,M).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (LDE,N)
!>          The first right hand side of the Sylvester Equation.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of the array E.  LDE >= max(1,M).
!> \endverbatim
!>
!> \param[in] F
!> \verbatim
!>          F is REAL array, dimension (LDF,N)
!>          The second right hand side of the Sylvester Equation.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F.  LDF >= max(1,M).
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is REAL
!>          SCALE is a scaling factor to prevent the overflow in the result.
!> \endverbatim
!>
!> \return
!> \verbatim
!>          The function returns the relative residual as defined in (1). If an error
!>          occurs a negative integer I is returned, signalizing that the parameter -I
!>          has a wrong/illegal value. If I = -1000, the auxiliary memory could not be allocated.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Martin Koehler, MPI Magdeburg
!
!> \date January 2024
!> \ingroup residual
!
FUNCTION SLA_RESIDUAL_CSYLV(TRANSA, TRANSB, SGN1, SGN2, M, N, A, LDA, B, LDB, C, LDC, D, LDD,  &
        & R, LDR, L, LDL, E, LDE, F, LDF, SCALE)
    IMPLICIT NONE
    CHARACTER :: TRANSA, TRANSB
    INTEGER :: M, N, LDA, LDB, LDC, LDD, LDE, LDF, LDR, LDL
    REAL :: A(LDA, *), B(LDB, *), C(LDC,*), D(LDD, *), E(LDE, *), F(LDF, *), R(LDR, *), L(LDL, *)
    REAL :: SCALE, SGN1, SGN2, SLA_RESIDUAL_CSYLV

    REAL :: SLANGE
    EXTERNAL :: SGEMM, SLANGE, LSAME, SLACPY, XERROR_HANDLER
    LOGICAL :: LSAME

    INTEGER :: INFO, ALLOC_STATUS
    REAL :: NRMA, NRMR, NRML, NRMR1, NRMR2, NRMB, NRMC, NRMD, NRME, NRMF
    REAL, ALLOCATABLE :: WORK(:,:)
    REAL, PARAMETER :: DONE = 1.0
    REAL, PARAMETER :: DZERO = 0.0
    REAL, PARAMETER :: MONE = -1.0
    REAL :: DUMMY(1)

    INFO = 0
    IF ( .NOT. LSAME(TRANSA, "N") .AND. .NOT. LSAME(TRANSA, "T" )) THEN
        INFO = -1
    ELSE IF ( .NOT. LSAME(TRANSB, "N") .AND. .NOT. LSAME(TRANSB, "T" )) THEN
        INFO = -2
    ELSE IF ( SGN1 .NE. DONE .AND. SGN1 .NE. MONE) THEN
        INFO = -3
    ELSE IF ( SGN2 .NE. DONE .AND. SGN2 .NE. MONE) THEN
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
    ELSE IF ( LDR .LT. MAX(1, M)) THEN
        INFO = -16
    ELSE IF ( LDL .LT. MAX(1, M)) THEN
        INFO = -18
    ELSE IF ( LDE .LT. MAX(1, M)) THEN
        INFO = -20
    ELSE IF ( LDF .LT. MAX(1, M)) THEN
        INFO = -22
    END IF

    IF (INFO .NE. 0 ) THEN
        SLA_RESIDUAL_CSYLV = REAL(INFO)
        CALL XERROR_HANDLER("SLA_RESIDUAL_CSYLV", -INFO)
        RETURN
    END IF

    IF ( M .EQ. 0 .OR. N .EQ. 0) THEN
        SLA_RESIDUAL_CSYLV = 0.0
        RETURN
    END IF

    ALLOCATE(WORK(M,N), STAT = ALLOC_STATUS)
    IF ( ALLOC_STATUS .NE. 0) THEN
        SLA_RESIDUAL_CSYLV = REAL(-1000)
        RETURN
    END IF

    NRMA = SLANGE("F", M, M, A, LDA, WORK)
    NRMB = SLANGE("F", N, N, B, LDB, WORK)
    NRMC = SLANGE("F", M, M, C, LDC, WORK)
    NRMD = SLANGE("F", N, N, D, LDD, WORK)
    NRMR = SLANGE("F", M, N, R, LDR, WORK)
    NRML = SLANGE("F", M, N, L, LDL, WORK)
    NRME = SLANGE("F", M, N, E, LDE, WORK)
    NRMF = SLANGE("F", M, N, F, LDF, WORK)

    CALL SLACPY("A", M, N, E, LDE, WORK(1,1), M)
    CALL SGEMM(TRANSA, "N", M, N, M, MONE, A, LDA, R, LDR, SCALE, WORK, M)
    CALL SGEMM("N", TRANSB, M, N, N, MONE*SGN1 , L, LDL, B, LDB, DONE, WORK, M)
    NRMR1 = SLANGE("F", M, N, WORK(1,1), M, DUMMY)

    CALL SLACPY("A", M, N, F, LDF, WORK(1,1), M)
    CALL SGEMM(TRANSA, "N", M, N, M, MONE, C, LDC, R, LDR, SCALE, WORK, M)
    CALL SGEMM("N", TRANSB, M, N, N, MONE*SGN2 , L, LDL, D, LDD, DONE, WORK, M)
    NRMR2 = SLANGE("F", M, N, WORK(1,1), M, DUMMY)


    SLA_RESIDUAL_CSYLV = MAX(NRMR1/ ( (NRMA  * NRMR + NRMB * NRML ) + SCALE * NRME), &
                        &  NRMR2 / ( (NRMC  * NRMR + NRMD * NRML ) + SCALE * NRMF))


    DEALLOCATE(WORK)
    RETURN

END FUNCTION SLA_RESIDUAL_CSYLV

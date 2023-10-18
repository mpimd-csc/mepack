*> \brief <b> SGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* This is a modified version of LAPACK SGGES to allow a
* generalized Hessenberg form on input
*
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGGES + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgges.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgges.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgges.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SHGGES( ASHAPE, JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB,
*                         SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR,
*                         LDVSR, WORK, LWORK, BWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          ASHAPE, JOBVSL, JOBVSR, SORT
*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
*       ..
*       .. Array Arguments ..
*       LOGICAL            BWORK( * )
*       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
*      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ),
*      $                   VSR( LDVSR, * ), WORK( * )
*       ..
*       .. Function Arguments ..
*       LOGICAL            SELCTG
*       EXTERNAL           SELCTG
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SHGGES computes for a pair of N-by-N real matrices (A,B) which are
*> either nonsymmetric or generalized Hessenberg form,
*> the generalized eigenvalues, the generalized real Schur form (S,T),
*> optionally, the left and/or right matrices of Schur vectors (VSL and
*> VSR). This gives the generalized Schur factorization
*>
*>          (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
*>
*> Optionally, it also orders the eigenvalues so that a selected cluster
*> of eigenvalues appears in the leading diagonal blocks of the upper
*> quasi-triangular matrix S and the upper triangular matrix T.The
*> leading columns of VSL and VSR then form an orthonormal basis for the
*> corresponding left and right eigenspaces (deflating subspaces).
*>
*> (If only the generalized eigenvalues are needed, use the driver
*> SGGEV instead, which is faster.)
*>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
*> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
*> usually represented as the pair (alpha,beta), as there is a
*> reasonable interpretation for beta=0 or both being zero.
*>
*> A pair of matrices (S,T) is in generalized real Schur form if T is
*> upper triangular with non-negative diagonal and S is block upper
*> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
*> to real generalized eigenvalues, while 2-by-2 blocks of S will be
*> "standardized" by making the corresponding elements of T have the
*> form:
*>         [  a  0  ]
*>         [  0  b  ]
*>
*> and the pair of corresponding 2-by-2 blocks in S and T will have a
*> complex conjugate pair of generalized eigenvalues.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ASHAPE
*> \verbatim
*>          ASHAPE is CHARACTER*1
*>          = 'G': (A,B) is a pair of real nonsymmetric matrices;
*>          = 'H': (A,B) is in generalized Hessenberg form.
*>     
*> \param[in] JOBVSL
*> \verbatim
*>          JOBVSL is CHARACTER*1
*>          = 'N':  do not compute the left Schur vectors;
*>          = 'V':  compute the left Schur vectors.
*> \endverbatim
*>
*> \param[in] JOBVSR
*> \verbatim
*>          JOBVSR is CHARACTER*1
*>          = 'N':  do not compute the right Schur vectors;
*>          = 'V':  compute the right Schur vectors.
*> \endverbatim
*>
*> \param[in] SORT
*> \verbatim
*>          SORT is CHARACTER*1
*>          Specifies whether or not to order the eigenvalues on the
*>          diagonal of the generalized Schur form.
*>          = 'N':  Eigenvalues are not ordered;
*>          = 'S':  Eigenvalues are ordered (see SELCTG);
*> \endverbatim
*>
*> \param[in] SELCTG
*> \verbatim
*>          SELCTG is a LOGICAL FUNCTION of three REAL arguments
*>          SELCTG must be declared EXTERNAL in the calling subroutine.
*>          If SORT = 'N', SELCTG is not referenced.
*>          If SORT = 'S', SELCTG is used to select eigenvalues to sort
*>          to the top left of the Schur form.
*>          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
*>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
*>          one of a complex conjugate pair of eigenvalues is selected,
*>          then both complex eigenvalues are selected.
*>
*>          Note that in the ill-conditioned case, a selected complex
*>          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),
*>          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2
*>          in this case.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VSL, and VSR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the first of the pair of matrices.
*>          On exit, A has been overwritten by its generalized Schur
*>          form S.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB, N)
*>          On entry, the second of the pair of matrices.
*>          On exit, B has been overwritten by its generalized Schur
*>          form T.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] SDIM
*> \verbatim
*>          SDIM is INTEGER
*>          If SORT = 'N', SDIM = 0.
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*>          for which SELCTG is true.  (Complex conjugate pairs for which
*>          SELCTG is true for either eigenvalue count as 2.)
*> \endverbatim
*>
*> \param[out] ALPHAR
*> \verbatim
*>          ALPHAR is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] ALPHAI
*> \verbatim
*>          ALPHAI is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is REAL array, dimension (N)
*>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
*>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
*>          and  BETA(j),j=1,...,N are the diagonals of the complex Schur
*>          form (S,T) that would result if the 2-by-2 diagonal blocks of
*>          the real Schur form of (A,B) were further reduced to
*>          triangular form using 2-by-2 complex unitary transformations.
*>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
*>          positive, then the j-th and (j+1)-st eigenvalues are a
*>          complex conjugate pair, with ALPHAI(j+1) negative.
*>
*>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
*>          may easily over- or underflow, and BETA(j) may even be zero.
*>          Thus, the user should avoid naively computing the ratio.
*>          However, ALPHAR and ALPHAI will be always less than and
*>          usually comparable with norm(A) in magnitude, and BETA always
*>          less than and usually comparable with norm(B).
*> \endverbatim
*>
*> \param[out] VSL
*> \verbatim
*>          VSL is REAL array, dimension (LDVSL,N)
*>          If JOBVSL = 'V', VSL will contain the left Schur vectors.
*>          Not referenced if JOBVSL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSL
*> \verbatim
*>          LDVSL is INTEGER
*>          The leading dimension of the matrix VSL. LDVSL >=1, and
*>          if JOBVSL = 'V', LDVSL >= N.
*> \endverbatim
*>
*> \param[out] VSR
*> \verbatim
*>          VSR is REAL array, dimension (LDVSR,N)
*>          If JOBVSR = 'V', VSR will contain the right Schur vectors.
*>          Not referenced if JOBVSR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSR
*> \verbatim
*>          LDVSR is INTEGER
*>          The leading dimension of the matrix VSR. LDVSR >= 1, and
*>          if JOBVSR = 'V', LDVSR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If N = 0, LWORK >= 1, else LWORK >= max(8*N,6*N+16).
*>          For good performance , LWORK must generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] BWORK
*> \verbatim
*>          BWORK is LOGICAL array, dimension (N)
*>          Not referenced if SORT = 'N'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          = 1,...,N:
*>                The QZ iteration failed.  (A,B) are not in Schur
*>                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
*>                be correct for j=INFO+1,...,N.
*>          > N:  =N+1: other than QZ iteration failed in SHGEQZ.
*>                =N+2: after reordering, roundoff changed values of
*>                      some complex eigenvalues so that leading
*>                      eigenvalues in the Generalized Schur form no
*>                      longer satisfy SELCTG=.TRUE.  This could also
*>                      be caused due to scaling.
*>                =N+3: reordering failed in STGSEN.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup realGEeigen
*
*  =====================================================================
      SUBROUTINE shgges( ASHAPE, JOBVSL, JOBVSR, SORT, SELCTG, N, A,
     $                  LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL,
     $                  LDVSL, VSR, LDVSR, WORK, LWORK, BWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          ASHAPE, JOBVSL, JOBVSR, SORT
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   b( ldb, * ), beta( * ), vsl( ldvsl, * ),
     $                   vsr( ldvsr, * ), work( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0e+0, one = 1.0e+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL,
     $                   LQUERY, LST2SL, WANTST, AHESS
      INTEGER            I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT,
     $                   ILO, IP, IRIGHT, IROWS, ITAU, IWRK, MAXWRK,
     $                   minwrk
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL,
     $                   PVSR, SAFMAX, SAFMIN, SMLNUM
      CHARACTER          COMPQ, COMPZ
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      REAL               DIF( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           sgeqrf, sggbak, sggbal, sgghrd, shgeqz, slabad,
     $                   slacpy, slascl, slaset, sorgqr, sormqr, stgsen,
     $                   xerbla
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANGE
      EXTERNAL           lsame, ilaenv, slamch, slange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( lsame( jobvsl, 'N' ) ) THEN
         ijobvl = 1
         ilvsl = .false.
      ELSE IF( lsame( jobvsl, 'V' ) ) THEN
         ijobvl = 2
         ilvsl = .true.
      ELSE
         ijobvl = -1
         ilvsl = .false.
      END IF
*
      IF( lsame( jobvsr, 'N' ) ) THEN
         ijobvr = 1
         ilvsr = .false.
      ELSE IF( lsame( jobvsr, 'V' ) ) THEN
         ijobvr = 2
         ilvsr = .true.
      ELSE
         ijobvr = -1
         ilvsr = .false.
      END IF
*
      wantst = lsame( sort, 'S' )
      ahess = lsame( ashape, 'H' )
*
*     Test the input arguments
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      IF( ijobvl.LE.0 ) THEN
         info = -2
      ELSE IF( ijobvr.LE.0 ) THEN
         info = -3
      ELSE IF( ( .NOT.wantst ) .AND. ( .NOT.lsame( sort, 'N' ) ) ) THEN
         info = -4
      ELSE IF( n.LT.0 ) THEN
         info = -6
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -8
      ELSE IF( ldb.LT.max( 1, n ) ) THEN
         info = -10
      ELSE IF( ldvsl.LT.1 .OR. ( ilvsl .AND. ldvsl.LT.n ) ) THEN
         info = -16
      ELSE IF( ldvsr.LT.1 .OR. ( ilvsr .AND. ldvsr.LT.n ) ) THEN
         info = -18
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.)
*
      IF( info.EQ.0 ) THEN
         IF( n.GT.0 )THEN
            minwrk = max( 8*n, 6*n + 16 )
            maxwrk = minwrk - n +
     $               n*ilaenv( 1, 'SGEQRF', ' ', n, 1, n, 0 )
            maxwrk = max( maxwrk, minwrk - n +
     $                    n*ilaenv( 1, 'SORMQR', ' ', n, 1, n, -1 ) )
            IF( ilvsl ) THEN
               maxwrk = max( maxwrk, minwrk - n +
     $                       n*ilaenv( 1, 'SORGQR', ' ', n, 1, n, -1 ) )
            END IF
         ELSE
            minwrk = 1
            maxwrk = 1
         END IF
         work( 1 ) = maxwrk
*
         IF( lwork.LT.minwrk .AND. .NOT.lquery )
     $      info = -20
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'SGGES ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 ) THEN
         sdim = 0
         RETURN
      END IF
*
*     Get machine constants
*
      eps = slamch( 'P' )
      safmin = slamch( 'S' )
      safmax = one / safmin
      CALL slabad( safmin, safmax )
      smlnum = sqrt( safmin ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = slange( 'M', n, n, a, lda, work )
      ilascl = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         anrmto = smlnum
         ilascl = .true.
      ELSE IF( anrm.GT.bignum ) THEN
         anrmto = bignum
         ilascl = .true.
      END IF
      IF( ilascl )
     $   CALL slascl( 'G', 0, 0, anrm, anrmto, n, n, a, lda, ierr )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      bnrm = slange( 'M', n, n, b, ldb, work )
      ilbscl = .false.
      IF( bnrm.GT.zero .AND. bnrm.LT.smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .true.
      ELSE IF( bnrm.GT.bignum ) THEN
         bnrmto = bignum
         ilbscl = .true.
      END IF
      IF( ilbscl )
     $   CALL slascl( 'G', 0, 0, bnrm, bnrmto, n, n, b, ldb, ierr )
*
*     Permute the matrix to make it more nearly triangular
*     (Workspace: need 6*N + 2*N space for storing balancing factors)
*
      ileft = 1
      iright = n + 1
      iwrk = iright + n
      CALL sggbal( 'P', n, a, lda, b, ldb, ilo, ihi, work( ileft ),
     $             work( iright ), work( iwrk ), ierr )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Workspace: need N, prefer N*NB)
*
      irows = ihi + 1 - ilo
      icols = n + 1 - ilo
      itau = iwrk
      iwrk = itau + irows
      IF ( .NOT.ahess ) THEN
         CALL sgeqrf( irows, icols, b( ilo, ilo ), ldb, work( itau ),
     $             work( iwrk ), lwork+1-iwrk, ierr )
      END IF
*
*     Apply the orthogonal transformation to matrix A
*     (Workspace: need N, prefer N*NB)
*
      IF ( .NOT.ahess ) THEN
         CALL sormqr( 'L', 'T', irows, icols, irows, b( ilo, ilo ),
     $             ldb, work( itau ), a( ilo, ilo ), lda, work( iwrk ),
     $             lwork+1-iwrk, ierr )
      END IF
*
*     Initialize VSL
*     (Workspace: need N, prefer N*NB)
*
      IF( ilvsl ) THEN
         CALL slaset( 'Full', n, n, zero, one, vsl, ldvsl )
         IF( irows.GT.1 ) THEN
            CALL slacpy( 'L', irows-1, irows-1, b( ilo+1, ilo ), ldb,
     $                   vsl( ilo+1, ilo ), ldvsl )
         END IF
         IF ( .NOT.ahess ) THEN
            CALL sorgqr( irows, irows, irows, vsl( ilo, ilo ), ldvsl,
     $                work( itau ), work( iwrk ), lwork+1-iwrk, ierr )
         END IF
      END IF
*
*     Initialize VSR
*
      IF( ilvsr )
     $   CALL slaset( 'Full', n, n, zero, one, vsr, ldvsr )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      IF ( .NOT.ahess ) THEN
      CALL sgghrd( jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl,
     $             ldvsl, vsr, ldvsr, ierr )
      END IF
*
*     Perform QZ algorithm, computing Schur vectors if desired
*     (Workspace: need N)
*
      compq = jobvsl
      compz = jobvsr
      iwrk = itau
      IF( ahess .AND. ilvsl) THEN
         compq = 'I'
      END IF
      IF( ahess .AND. ilvsr) THEN
         compz = 'I'
      END IF
      CALL shgeqz( 'S', compq, compz, n, ilo, ihi, a, lda, b, ldb,
     $             alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr,
     $             work( iwrk ), lwork+1-iwrk, ierr )
      IF( ierr.NE.0 ) THEN
         IF( ierr.GT.0 .AND. ierr.LE.n ) THEN
            info = ierr
         ELSE IF( ierr.GT.n .AND. ierr.LE.2*n ) THEN
            info = ierr - n
         ELSE
            info = n + 1
         END IF
         GO TO 40
      END IF
*
*     Sort eigenvalues ALPHA/BETA if desired
*     (Workspace: need 4*N+16 )
*
      sdim = 0
      IF( wantst ) THEN
*
*        Undo scaling on eigenvalues before SELCTGing
*
         IF( ilascl ) THEN
            CALL slascl( 'G', 0, 0, anrmto, anrm, n, 1, alphar, n,
     $                   ierr )
            CALL slascl( 'G', 0, 0, anrmto, anrm, n, 1, alphai, n,
     $                   ierr )
         END IF
         IF( ilbscl )
     $      CALL slascl( 'G', 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr )
*
*        Select eigenvalues
*
         DO 10 i = 1, n
            bwork( i ) = selctg( alphar( i ), alphai( i ), beta( i ) )
   10    CONTINUE
*
         CALL stgsen( 0, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alphar,
     $                alphai, beta, vsl, ldvsl, vsr, ldvsr, sdim, pvsl,
     $                pvsr, dif, work( iwrk ), lwork-iwrk+1, idum, 1,
     $                ierr )
         IF( ierr.EQ.1 )
     $      info = n + 3
*
      END IF
*
*     Apply back-permutation to VSL and VSR
*     (Workspace: none needed)
*
      IF( ilvsl )
     $   CALL sggbak( 'P', 'L', n, ilo, ihi, work( ileft ),
     $                work( iright ), n, vsl, ldvsl, ierr )
*
      IF( ilvsr )
     $   CALL sggbak( 'P', 'R', n, ilo, ihi, work( ileft ),
     $                work( iright ), n, vsr, ldvsr, ierr )
*
*     Check if unscaling would cause over/underflow, if so, rescale
*     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
*     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)
*
      IF( ilascl )THEN
         DO 50 i = 1, n
            IF( alphai( i ).NE.zero ) THEN
               IF( ( alphar( i )/safmax ).GT.( anrmto/anrm ) .OR.
     $             ( safmin/alphar( i ) ).GT.( anrm/anrmto ) ) THEN
                  work( 1 ) = abs( a( i, i )/alphar( i ) )
                  beta( i ) = beta( i )*work( 1 )
                  alphar( i ) = alphar( i )*work( 1 )
                  alphai( i ) = alphai( i )*work( 1 )
               ELSE IF( ( alphai( i )/safmax ).GT.( anrmto/anrm ) .OR.
     $             ( safmin/alphai( i ) ).GT.( anrm/anrmto ) ) THEN
                  work( 1 ) = abs( a( i, i+1 )/alphai( i ) )
                  beta( i ) = beta( i )*work( 1 )
                  alphar( i ) = alphar( i )*work( 1 )
                  alphai( i ) = alphai( i )*work( 1 )
               END IF
            END IF
   50    CONTINUE
      END IF
*
      IF( ilbscl )THEN
         DO 60 i = 1, n
            IF( alphai( i ).NE.zero ) THEN
                IF( ( beta( i )/safmax ).GT.( bnrmto/bnrm ) .OR.
     $              ( safmin/beta( i ) ).GT.( bnrm/bnrmto ) ) THEN
                   work( 1 ) = abs(b( i, i )/beta( i ))
                   beta( i ) = beta( i )*work( 1 )
                   alphar( i ) = alphar( i )*work( 1 )
                   alphai( i ) = alphai( i )*work( 1 )
                END IF
             END IF
   60    CONTINUE
      END IF
*
*     Undo scaling
*
      IF( ilascl ) THEN
         CALL slascl( 'H', 0, 0, anrmto, anrm, n, n, a, lda, ierr )
         CALL slascl( 'G', 0, 0, anrmto, anrm, n, 1, alphar, n, ierr )
         CALL slascl( 'G', 0, 0, anrmto, anrm, n, 1, alphai, n, ierr )
      END IF
*
      IF( ilbscl ) THEN
         CALL slascl( 'U', 0, 0, bnrmto, bnrm, n, n, b, ldb, ierr )
         CALL slascl( 'G', 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr )
      END IF
*
      IF( wantst ) THEN
*
*        Check if reordering is correct
*
         lastsl = .true.
         lst2sl = .true.
         sdim = 0
         ip = 0
         DO 30 i = 1, n
            cursl = selctg( alphar( i ), alphai( i ), beta( i ) )
            IF( alphai( i ).EQ.zero ) THEN
               IF( cursl )
     $            sdim = sdim + 1
               ip = 0
               IF( cursl .AND. .NOT.lastsl )
     $            info = n + 2
            ELSE
               IF( ip.EQ.1 ) THEN
*
*                 Last eigenvalue of conjugate pair
*
                  cursl = cursl .OR. lastsl
                  lastsl = cursl
                  IF( cursl )
     $               sdim = sdim + 2
                  ip = -1
                  IF( cursl .AND. .NOT.lst2sl )
     $               info = n + 2
               ELSE
*
*                 First eigenvalue of conjugate pair
*
                  ip = 1
               END IF
            END IF
            lst2sl = lastsl
            lastsl = cursl
   30    CONTINUE
*
      END IF
*
   40 CONTINUE
*
      work( 1 ) = maxwrk
*
      RETURN
*
*     End of SGGES
*
      END

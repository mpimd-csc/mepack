*> \brief <b> SHGEES computes the eigenvalues, the Schur form, and,
* optionally, the matrix of Schur vectors for GE or HESS matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* This is a modified version of LAPACK SGEES to allow an
* upper-Hessenberg form on input
*
* Online html documentation of SGEES available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SGEES + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgees.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgees.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgees.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SHGEES( ASHAPE, JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,
*                         VS, LDVS, WORK, LWORK, BWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVS, SORT
*       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
*       ..
*       .. Array Arguments ..
*       LOGICAL            BWORK( * )
*       REAL               A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
*      $                   WR( * )
*       ..
*       .. Function Arguments ..
*       LOGICAL            SELECT
*       EXTERNAL           SELECT
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SHGEES computes for a matrix A,
*> which is either an N-by-N real nonsymmetric matrix or
*> an upper Hessenberg matrix, the eigenvalues, the real Schur form T,
*> and, optionally, the matrix of Schur vectors Z.
*> This gives the Schur factorization A = Z*T*(Z**T).
*>
*> Optionally, it also orders the eigenvalues on the diagonal of the
*> real Schur form so that selected eigenvalues are at the top left.
*> The leading columns of Z then form an orthonormal basis for the
*> invariant subspace corresponding to the selected eigenvalues.
*>
*> A matrix is in real Schur form if it is upper quasi-triangular with
*> 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
*> form
*>         [  a  b  ]
*>         [  c  a  ]
*>
*> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ASHAPE
*> \verbatim
*>          ASHAPE is CHARACTER*1
*>          = 'G': A is a full matrix;
*>          = 'H': A is an upper Hessenberg matrix.
*> \endverbatim
*>
*> \param[in] JOBVS
*> \verbatim
*>          JOBVS is CHARACTER*1
*>          = 'N': Schur vectors are not computed;
*>          = 'V': Schur vectors are computed.
*> \endverbatim
*>
*> \param[in] SORT
*> \verbatim
*>          SORT is CHARACTER*1
*>          Specifies whether or not to order the eigenvalues on the
*>          diagonal of the Schur form.
*>          = 'N': Eigenvalues are not ordered;
*>          = 'S': Eigenvalues are ordered (see SELECT).
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is a LOGICAL FUNCTION of two REAL arguments
*>          SELECT must be declared EXTERNAL in the calling subroutine.
*>          If SORT = 'S', SELECT is used to select eigenvalues to sort
*>          to the top left of the Schur form.
*>          If SORT = 'N', SELECT is not referenced.
*>          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
*>          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
*>          conjugate pair of eigenvalues is selected, then both complex
*>          eigenvalues are selected.
*>          Note that a selected complex eigenvalue may no longer
*>          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
*>          ordering may change the value of complex eigenvalues
*>          (especially if the eigenvalue is ill-conditioned); in this
*>          case INFO is set to N+2 (see INFO below).
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the N-by-N matrix A.
*>          On exit, A has been overwritten by its real Schur form T.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] SDIM
*> \verbatim
*>          SDIM is INTEGER
*>          If SORT = 'N', SDIM = 0.
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*>                         for which SELECT is true. (Complex conjugate
*>                         pairs for which SELECT is true for either
*>                         eigenvalue count as 2.)
*> \endverbatim
*>
*> \param[out] WR
*> \verbatim
*>          WR is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] WI
*> \verbatim
*>          WI is REAL array, dimension (N)
*>          WR and WI contain the real and imaginary parts,
*>          respectively, of the computed eigenvalues in the same order
*>          that they appear on the diagonal of the output Schur form T.
*>          Complex conjugate pairs of eigenvalues will appear
*>          consecutively with the eigenvalue having the positive
*>          imaginary part first.
*> \endverbatim
*>
*> \param[out] VS
*> \verbatim
*>          VS is REAL array, dimension (LDVS,N)
*>          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
*>          vectors.
*>          If JOBVS = 'N', VS is not referenced.
*> \endverbatim
*>
*> \param[in] LDVS
*> \verbatim
*>          LDVS is INTEGER
*>          The leading dimension of the array VS.  LDVS >= 1; if
*>          JOBVS = 'V', LDVS >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,3*N).
*>          For good performance, LWORK must generally be larger.
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
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value.
*>          > 0: if INFO = i, and i is
*>             <= N: the QR algorithm failed to compute all the
*>                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
*>                   contain those eigenvalues which have converged; if
*>                   JOBVS = 'V', VS contains the matrix which reduces A
*>                   to its partially converged Schur form.
*>             = N+1: the eigenvalues could not be reordered because some
*>                   eigenvalues were too close to separate (the problem
*>                   is very ill-conditioned);
*>             = N+2: after reordering, roundoff changed values of some
*>                   complex eigenvalues so that leading eigenvalues in
*>                   the Schur form no longer satisfy SELECT=.TRUE.  This
*>                   could also be caused by underflow due to scaling.
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE shgees( ASHAPE, JOBVS, SORT, SELECT, N, A, LDA, SDIM,
     $                  WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          ASHAPE, JOBVS, SORT
      INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      REAL               A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
     $                   wr( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      parameter( zero = 0.0e0, one = 1.0e0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            CURSL, LASTSL, LQUERY, LST2SL, SCALEA, WANTST,
     $                   wantvs, AHESS
      INTEGER            HSWORK, I, I1, I2, IBAL, ICOND, IERR, IEVAL,
     $                   ihi, ilo, inxt, ip, itau, iwrk, maxwrk, minwrk
      REAL               ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
      CHARACTER          COMPZ
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      REAL               DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           scopy, sgebak, sgebal, sgehrd, shseqr, slabad,
     $                   slacpy, slascl, sorghr, sswap, strsen, xerbla
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANGE
      EXTERNAL           lsame, ilaenv, slamch, slange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      wantvs = lsame( jobvs, 'V' )
      wantst = lsame( sort, 'S' )
      ahess = lsame( ashape, 'H' )
      IF( ( .NOT.wantvs ) .AND. ( .NOT.lsame( jobvs, 'N' ) ) ) THEN
         info = -2
      ELSE IF( ( .NOT.wantst ) .AND. ( .NOT.lsame( sort, 'N' ) ) ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldvs.LT.1 .OR. ( wantvs .AND. ldvs.LT.n ) ) THEN
         info = -12
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by SHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
*
*     Workspace for ASHAPE = 'G'
      IF ( .NOT.AHESS ) THEN
         IF( info.EQ.0 ) THEN
            IF( n.EQ.0 ) THEN
               minwrk = 1
               maxwrk = 1
            ELSE
               maxwrk = 2*n + n*ilaenv( 1, 'SGEHRD', ' ', n, 1, n, 0 )
               minwrk = 3*n
*
               CALL shseqr( 'S', jobvs, n, 1, n, a, lda, wr, wi, vs,
     $             ldvs, work, -1, ieval )
               hswork = int( work( 1 ) )
*
               IF( .NOT.wantvs ) THEN
                  maxwrk = max( maxwrk, n + hswork )
               ELSE
                  maxwrk = max( maxwrk, 2*n + ( n - 1 )*ilaenv( 1,
     $                       'SORGHR', ' ', n, 1, n, -1 ) )
                  maxwrk = max( maxwrk, n + hswork )
               END IF
            END IF
            work( 1 ) = maxwrk
*
            IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
               info = -14
            END IF
         END IF
*
*     Workspace for ASHAPE = 'H'
      ELSE
         IF( info.EQ.0 ) THEN
            IF( n.EQ.0 ) THEN
               minwrk = 1
               maxwrk = 1
            ELSE
               minwrk = 1+n
*
               CALL shseqr( 'S', jobvs, n, 1, n, a, lda, wr, wi, vs,
     $             ldvs, work, -1, ieval )
               hswork = int( work( 1 ) )
*
               maxwrk = n+hswork
            END IF
            work( 1 ) = maxwrk
*
            IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
               info = -14
            END IF
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'SHGEES ', -info )
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
      smlnum = slamch( 'S' )
      bignum = one / smlnum
      CALL slabad( smlnum, bignum )
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = slange( 'M', n, n, a, lda, dum )
      scalea = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         scalea = .true.
         cscale = smlnum
      ELSE IF( anrm.GT.bignum ) THEN
         scalea = .true.
         cscale = bignum
      END IF
      IF( scalea )
     $   CALL slascl( 'G', 0, 0, anrm, cscale, n, n, a, lda, ierr )
*
*     Permute the matrix to make it more nearly triangular
*     (Workspace: need N)
*
      ibal = 1
      CALL sgebal( 'P', n, a, lda, ilo, ihi, work( ibal ), ierr )
*
*     Reduce to upper Hessenberg form
*     (Workspace: need 3*N, prefer 2*N+N*NB)
*
      itau = n + ibal
      iwrk = n + itau
      IF ( .NOT.ahess ) THEN
         CALL sgehrd( n, ilo, ihi, a, lda, work( itau ), work( iwrk ),
     $             lwork-iwrk+1, ierr )
*
         IF( wantvs ) THEN
*
*           Copy Householder vectors to VS
*
            CALL slacpy( 'L', n, n, a, lda, vs, ldvs )
*
*           Generate orthogonal matrix in VS
*           (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*
            CALL sorghr( n, ilo, ihi, vs, ldvs, work( itau ),
     $                work( iwrk ), lwork-iwrk+1, ierr )
         END IF
      END IF
*
      sdim = 0
*
*     Perform QR iteration, accumulating Schur vectors in VS if desired
*     (Workspace: need N+1, prefer N+HSWORK (see comments) )
*
      iwrk = itau
      COMPZ = jobvs
      IF( ahess .AND. wantvs ) THEN
         COMPZ = 'I'
      END IF
*
      CALL shseqr( 'S', COMPZ, n, ilo, ihi, a, lda, wr, wi, vs, ldvs,
     $             work( iwrk ), lwork-iwrk+1, ieval )
      IF( ieval.GT.0 )
     $   info = ieval
*
*     Sort eigenvalues if desired
*
      IF( wantst .AND. info.EQ.0 ) THEN
         IF( scalea ) THEN
            CALL slascl( 'G', 0, 0, cscale, anrm, n, 1, wr, n, ierr )
            CALL slascl( 'G', 0, 0, cscale, anrm, n, 1, wi, n, ierr )
         END IF
         DO 10 i = 1, n
            bwork( i ) = SELECT( wr( i ), wi( i ) )
   10    CONTINUE
*
*        Reorder eigenvalues and transform Schur vectors
*        (Workspace: none needed)
*
         CALL strsen( 'N', jobvs, bwork, n, a, lda, vs, ldvs, wr, wi,
     $                sdim, s, sep, work( iwrk ), lwork-iwrk+1, idum, 1,
     $                icond )
         IF( icond.GT.0 )
     $      info = n + icond
      END IF
*
      IF( wantvs ) THEN
*
*        Undo balancing
*        (Workspace: need N)
*
         CALL sgebak( 'P', 'R', n, ilo, ihi, work( ibal ), n, vs, ldvs,
     $                ierr )
      END IF
*
      IF( scalea ) THEN
*
*        Undo scaling for the Schur form of A
*
         CALL slascl( 'H', 0, 0, cscale, anrm, n, n, a, lda, ierr )
         CALL scopy( n, a, lda+1, wr, 1 )
         IF( cscale.EQ.smlnum ) THEN
*
*           If scaling back towards underflow, adjust WI if an
*           offdiagonal element of a 2-by-2 block in the Schur form
*           underflows.
*
            IF( ieval.GT.0 ) THEN
               i1 = ieval + 1
               i2 = ihi - 1
               CALL slascl( 'G', 0, 0, cscale, anrm, ilo-1, 1, wi,
     $                      max( ilo-1, 1 ), ierr )
            ELSE IF( wantst ) THEN
               i1 = 1
               i2 = n - 1
            ELSE
               i1 = ilo
               i2 = ihi - 1
            END IF
            inxt = i1 - 1
            DO 20 i = i1, i2
               IF( i.LT.inxt )
     $            GO TO 20
               IF( wi( i ).EQ.zero ) THEN
                  inxt = i + 1
               ELSE
                  IF( a( i+1, i ).EQ.zero ) THEN
                     wi( i ) = zero
                     wi( i+1 ) = zero
                  ELSE IF( a( i+1, i ).NE.zero .AND. a( i, i+1 ).EQ.
     $                     zero ) THEN
                     wi( i ) = zero
                     wi( i+1 ) = zero
                     IF( i.GT.1 )
     $                  CALL sswap( i-1, a( 1, i ), 1, a( 1, i+1 ), 1 )
                     IF( n.GT.i+1 )
     $                  CALL sswap( n-i-1, a( i, i+2 ), lda,
     $                              a( i+1, i+2 ), lda )
                     IF( wantvs ) THEN
                        CALL sswap( n, vs( 1, i ), 1, vs( 1, i+1 ), 1 )
                     END IF
                     a( i, i+1 ) = a( i+1, i )
                     a( i+1, i ) = zero
                  END IF
                  inxt = i + 2
               END IF
   20       CONTINUE
         END IF
*
*        Undo scaling for the imaginary part of the eigenvalues
*
         CALL slascl( 'G', 0, 0, cscale, anrm, n-ieval, 1,
     $                wi( ieval+1 ), max( n-ieval, 1 ), ierr )
      END IF
*
      IF( wantst .AND. info.EQ.0 ) THEN
*
*        Check if reordering successful
*
         lastsl = .true.
         lst2sl = .true.
         sdim = 0
         ip = 0
         DO 30 i = 1, n
            cursl = SELECT( wr( i ), wi( i ) )
            IF( wi( i ).EQ.zero ) THEN
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
      END IF
*
      work( 1 ) = maxwrk
      RETURN
*
*     End of SGEES
*
      END

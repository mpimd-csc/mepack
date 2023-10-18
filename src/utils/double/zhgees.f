*> \brief <b> ZHGEES computes the eigenvalues, the Schur form, and,
* optionally, the matrix of Schur vectors for GE or HESS matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* This is a modified version of LAPACK ZGEES to allow an
* upper-Hessenberg form on input
*
* Online html documentation of ZGEES available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEES + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgees.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgees.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgees.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHGEES( ASHAPE, JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,
*                         LDVS, WORK, LWORK, RWORK, BWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          ASHAPE, JOBVS, SORT
*       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
*       ..
*       .. Array Arguments ..
*       LOGICAL            BWORK( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
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
*> ZHGEES computes for a matrix A,
*> which is either an N-by-N complex nonsymmetric matrix A
*> or an upper Hessenberg matrix, the eigenvalues,
*> the Schur form T, and, optionally, the matrix of Schur
*> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
*>
*> Optionally, it also orders the eigenvalues on the diagonal of the
*> Schur form so that selected eigenvalues are at the top left.
*> The leading columns of Z then form an orthonormal basis for the
*> invariant subspace corresponding to the selected eigenvalues.
*>
*> A complex matrix is in Schur form if it is upper triangular.
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
*>          = 'N': Eigenvalues are not ordered:
*>          = 'S': Eigenvalues are ordered (see SELECT).
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is a LOGICAL FUNCTION of one COMPLEX*16 argument
*>          SELECT must be declared EXTERNAL in the calling subroutine.
*>          If SORT = 'S', SELECT is used to select eigenvalues to order
*>          to the top left of the Schur form.
*>          IF SORT = 'N', SELECT is not referenced.
*>          The eigenvalue W(j) is selected if SELECT(W(j)) is true.
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
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the N-by-N matrix A.
*>          On exit, A has been overwritten by its Schur form T.
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
*>          If SORT = 'S', SDIM = number of eigenvalues for which
*>                         SELECT is true.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (N)
*>          W contains the computed eigenvalues, in the same order that
*>          they appear on the diagonal of the output Schur form T.
*> \endverbatim
*>
*> \param[out] VS
*> \verbatim
*>          VS is COMPLEX*16 array, dimension (LDVS,N)
*>          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
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
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,2*N).
*>          For good performance, LWORK must generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (N)
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
*>               <= N:  the QR algorithm failed to compute all the
*>                      eigenvalues; elements 1:ILO-1 and i+1:N of W
*>                      contain those eigenvalues which have converged;
*>                      if JOBVS = 'V', VS contains the matrix which
*>                      reduces A to its partially converged Schur form.
*>               = N+1: the eigenvalues could not be reordered because
*>                      some eigenvalues were too close to separate (the
*>                      problem is very ill-conditioned);
*>               = N+2: after reordering, roundoff changed values of
*>                      some complex eigenvalues so that leading
*>                      eigenvalues in the Schur form no longer satisfy
*>                      SELECT = .TRUE..  This could also be caused by
*>                      underflow due to scaling.
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE zhgees( ASHAPE, JOBVS, SORT, SELECT, N, A, LDA, SDIM,
     $                  W, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO )
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
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELECT
      EXTERNAL           SELECT
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTST, WANTVS, AHESS
      INTEGER            HSWORK, I, IBAL, ICOND, IERR, IEVAL, IHI, ILO,
     $                   itau, iwrk, maxwrk, minwrk
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
      CHARACTER          COMPZ
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, xerbla, zcopy, zgebak, zgebal, zgehrd,
     $                   zhseqr, zlacpy, zlascl, ztrsen, zunghr
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           lsame, ilaenv, dlamch, zlange
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
         info = -11
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       CWorkspace refers to complex workspace, and RWorkspace to real
*       workspace. NB refers to the optimal block size for the
*       immediately following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by ZHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
*     Workspace for ASHAPE = 'G'
*
      IF( .NOT.ahess ) THEN
         IF( info.EQ.0 ) THEN
            IF( n.EQ.0 ) THEN
               minwrk = 1
               maxwrk = 1
            ELSE
               maxwrk = n + n*ilaenv( 1, 'ZGEHRD', ' ', n, 1, n, 0 )
               minwrk = 2*n
*
               CALL zhseqr( 'S', jobvs, n, 1, n, a, lda, w, vs, ldvs,
     $             work, -1, ieval )
               hswork = int( work( 1 ) )
*
               IF( .NOT.wantvs ) THEN
                  maxwrk = max( maxwrk, hswork )
               ELSE
                  maxwrk = max( maxwrk, n + ( n - 1 )*ilaenv( 1,
     $                      'ZUNGHR', ' ', n, 1, n, -1 ) )
                  maxwrk = max( maxwrk, hswork )
               END IF
            END IF
            work( 1 ) = maxwrk
*
            IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
               info = -13
            END IF
         END IF
*
*     Workspace for ASHAPE = 'H'
*
      ELSE
         IF( info.EQ.0 ) THEN
            IF( n.EQ.0 ) THEN
               minwrk = 1
               maxwrk = 1
            ELSE
               minwrk = n
*
               CALL zhseqr( 'S', jobvs, n, 1, n, a, lda, w, vs, ldvs,
     $             work, -1, ieval )
               hswork = int( work( 1 ) )
*
               maxwrk = hswork
            END IF
            work( 1 ) = maxwrk
*
            IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
               info = -13
            END IF
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEES ', -info )
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
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
      CALL dlabad( smlnum, bignum )
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = zlange( 'M', n, n, a, lda, dum )
      scalea = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         scalea = .true.
         cscale = smlnum
      ELSE IF( anrm.GT.bignum ) THEN
         scalea = .true.
         cscale = bignum
      END IF
      IF( scalea )
     $   CALL zlascl( 'G', 0, 0, anrm, cscale, n, n, a, lda, ierr )
*
*     Permute the matrix to make it more nearly triangular
*     (CWorkspace: none)
*     (RWorkspace: need N)
*
      ibal = 1
      CALL zgebal( 'P', n, a, lda, ilo, ihi, rwork( ibal ), ierr )
*
*     Reduce to upper Hessenberg form
*     (CWorkspace: need 2*N, prefer N+N*NB)
*     (RWorkspace: none)
*
      itau = 1
      iwrk = n + itau
      IF ( .NOT.ahess ) THEN
         CALL zgehrd( n, ilo, ihi, a, lda, work( itau ), work( iwrk ),
     $             lwork-iwrk+1, ierr )
*
         IF( wantvs ) THEN
*
*           Copy Householder vectors to VS
*
            CALL zlacpy( 'L', n, n, a, lda, vs, ldvs )
*
*           Generate unitary matrix in VS
*           (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*           (RWorkspace: none)
*
            CALL zunghr( n, ilo, ihi, vs, ldvs, work( itau ),
     $                work( iwrk ), lwork-iwrk+1, ierr )
         END IF
      END IF
*
      sdim = 0
*
*     Perform QR iteration, accumulating Schur vectors in VS if desired
*     (CWorkspace: need 1, prefer HSWORK (see comments) )
*     (RWorkspace: none)
*
      iwrk = itau
      compz = jobvs
      IF( ahess .AND. wantvs ) THEN
         compz = 'I'
      END IF
*
      CALL zhseqr( 'S', compz, n, ilo, ihi, a, lda, w, vs, ldvs,
     $             work( iwrk ), lwork-iwrk+1, ieval )
      IF( ieval.GT.0 )
     $   info = ieval
*
*     Sort eigenvalues if desired
*
      IF( wantst .AND. info.EQ.0 ) THEN
         IF( scalea )
     $      CALL zlascl( 'G', 0, 0, cscale, anrm, n, 1, w, n, ierr )
         DO 10 i = 1, n
            bwork( i ) = SELECT( w( i ) )
   10    CONTINUE
*
*        Reorder eigenvalues and transform Schur vectors
*        (CWorkspace: none)
*        (RWorkspace: none)
*
         CALL ztrsen( 'N', jobvs, bwork, n, a, lda, vs, ldvs, w, sdim,
     $                s, sep, work( iwrk ), lwork-iwrk+1, icond )
      END IF
*
      IF( wantvs ) THEN
*
*        Undo balancing
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL zgebak( 'P', 'R', n, ilo, ihi, rwork( ibal ), n, vs, ldvs,
     $                ierr )
      END IF
*
      IF( scalea ) THEN
*
*        Undo scaling for the Schur form of A
*
         CALL zlascl( 'U', 0, 0, cscale, anrm, n, n, a, lda, ierr )
         CALL zcopy( n, a, lda+1, w, 1 )
      END IF
*
      work( 1 ) = maxwrk
      RETURN
*
*     End of ZGEES
*
      END

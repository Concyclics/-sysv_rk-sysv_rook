* spot02 dpot02 csyt02 zsyt02

      SUBROUTINE spot02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, LDX, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), RWORK( * ),
     $                   x( ldx, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      parameter( zero = 0.0e+0, one = 1.0e+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      REAL               ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANSY
      EXTERNAL           sasum, slamch, slansy
*     ..
*     .. External Subroutines ..
      EXTERNAL           ssymm
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( n.LE.0 .OR. nrhs.LE.0 ) THEN
      resid = zero
      RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      eps = slamch( 'Epsilon' )
      anorm = slansy( '1', uplo, n, a, lda, rwork )
      IF( anorm.LE.zero ) THEN
      resid = one / eps
      RETURN
      END IF
*
*     Compute  B - A*X
*
      CALL ssymm( 'Left', uplo, n, nrhs, -one, a, lda, x, ldx, one, b,
     $            ldb )
*
*     Compute the maximum over the number of right hand sides of
*        norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .
*
      resid = zero
      DO 10 j = 1, nrhs
      bnorm = sasum( n, b( 1, j ), 1 )
      xnorm = sasum( n, x( 1, j ), 1 )
      IF( xnorm.LE.zero ) THEN
            resid = one / eps
      ELSE
            resid = max( resid, ( ( bnorm / anorm ) / xnorm ) / eps )
      END IF
10    CONTINUE
*
      RETURN
*
*     End of SPOT02
*
      END

      SUBROUTINE DPOT02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), RWORK( * ),
     $                   X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DASUM, DLAMCH, DLANSY
      EXTERNAL           DASUM, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSYMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute  B - A*X
*
      CALL DSYMM( 'Left', UPLO, N, NRHS, -ONE, A, LDA, X, LDX, ONE, B,
     $            LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = DASUM( N, B( 1, J ), 1 )
         XNORM = DASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of DPOT02
*
      END

*> \brief \b CSYT02
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CSYT02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK,
*                          RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, LDB, LDX, N, NRHS
*       REAL               RESID
*       ..
*       .. Array Arguments ..
*       REAL               RWORK( * )
*       COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CSYT02 computes the residual for a solution to a complex symmetric
*> system of linear equations  A*x = b:
*>
*>    RESID = norm(B - A*X) / ( norm(A) * norm(X) * EPS ),
*>
*> where EPS is the machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows and columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of columns of B, the matrix of right hand sides.
*>          NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The original complex symmetric matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N)
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX array, dimension (LDX,NRHS)
*>          The computed solution vectors for the system of linear
*>          equations.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.   LDX >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,NRHS)
*>          On entry, the right hand side vectors for the system of
*>          linear equations.
*>          On exit, B is overwritten with the difference B - A*X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
*>          The maximum over the number of right hand sides of
*>          norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
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
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CSYT02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, LDX, N, NRHS
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      REAL               ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      REAL               CLANSY, SCASUM, SLAMCH
      EXTERNAL           CLANSY, SCASUM, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSYMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANSY( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute  B - A*X  (or  B - A'*X ) and store in B .
*
      CALL CSYMM( 'Left', UPLO, N, NRHS, -CONE, A, LDA, X, LDX, CONE, B,
     $            LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = SCASUM( N, B( 1, J ), 1 )
         XNORM = SCASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM/ANORM )/XNORM )/EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of CSYT02
*
      END

*> \brief \b ZSYT02
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSYT02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK,
*                          RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, LDB, LDX, N, NRHS
*       DOUBLE PRECISION   RESID
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZSYT02 computes the residual for a solution to a complex symmetric
*> system of linear equations  A*x = b:
*>
*>    RESID = norm(B - A*X) / ( norm(A) * norm(X) * EPS ),
*>
*> where EPS is the machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows and columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of columns of B, the matrix of right hand sides.
*>          NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The original complex symmetric matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N)
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (LDX,NRHS)
*>          The computed solution vectors for the system of linear
*>          equations.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.   LDX >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,NRHS)
*>          On entry, the right hand side vectors for the system of
*>          linear equations.
*>          On exit, B is overwritten with the difference B - A*X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is DOUBLE PRECISION
*>          The maximum over the number of right hand sides of
*>          norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
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
*> \ingroup complex16_lin
*
*  =====================================================================
      SUBROUTINE ZSYT02( UPLO, N, NRHS, A, LDA, X, LDX, B, LDB, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANSY
      EXTERNAL           DLAMCH, DZASUM, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZSYMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if ANORM = 0.
*
      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF
*
*     Compute  B - A*X  (or  B - A'*X ) and store in B .
*
      CALL ZSYMM( 'Left', UPLO, N, NRHS, -CONE, A, LDA, X, LDX, CONE, B,
     $            LDB )
*
*     Compute the maximum over the number of right hand sides of
*        norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .
*
      RESID = ZERO
      DO 10 J = 1, NRHS
         BNORM = DZASUM( N, B( 1, J ), 1 )
         XNORM = DZASUM( N, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   10 CONTINUE
*
      RETURN
*
*     End of ZSYT02
*
      END

      SUBROUTINE sget04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDXACT, N, NRHS
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      REAL               X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      parameter( zero = 0.0e+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, J
      REAL               DIFFNM, EPS, XNORM
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH
      EXTERNAL           isamax, slamch
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( n.LE.0 .OR. nrhs.LE.0 ) THEN
         resid = zero
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if RCOND is invalid.
*
      eps = slamch( 'Epsilon' )
      IF( rcond.LT.zero ) THEN
         resid = 1.0 / eps
         RETURN
      END IF
*
*     Compute the maximum of
*        norm(X - XACT) / ( norm(XACT) * EPS )
*     over all the vectors X and XACT .
*
      resid = zero
      DO 20 j = 1, nrhs
         ix = isamax( n, xact( 1, j ), 1 )
         xnorm = abs( xact( ix, j ) )
         diffnm = zero
         DO 10 i = 1, n
            diffnm = max( diffnm, abs( x( i, j )-xact( i, j ) ) )
   10    CONTINUE
         IF( xnorm.LE.zero ) THEN
            IF( diffnm.GT.zero )
     $         resid = 1.0 / eps
         ELSE
            resid = max( resid, ( diffnm / xnorm )*rcond )
         END IF
   20 CONTINUE
      IF( resid*eps.LT.1.0 )
     $   resid = resid / eps
*
      RETURN
*
*     End of SGET04
*
      END

      SUBROUTINE dget04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDXACT, N, NRHS
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, J
      DOUBLE PRECISION   DIFFNM, EPS, XNORM
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           idamax, dlamch
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( n.LE.0 .OR. nrhs.LE.0 ) THEN
         resid = zero
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if RCOND is invalid.
*
      eps = dlamch( 'Epsilon' )
      IF( rcond.LT.zero ) THEN
         resid = 1.0d0 / eps
         RETURN
      END IF
*
*     Compute the maximum of
*        norm(X - XACT) / ( norm(XACT) * EPS )
*     over all the vectors X and XACT .
*
      resid = zero
      DO 20 j = 1, nrhs
         ix = idamax( n, xact( 1, j ), 1 )
         xnorm = abs( xact( ix, j ) )
         diffnm = zero
         DO 10 i = 1, n
            diffnm = max( diffnm, abs( x( i, j )-xact( i, j ) ) )
   10    CONTINUE
         IF( xnorm.LE.zero ) THEN
            IF( diffnm.GT.zero )
     $         resid = 1.0d0 / eps
         ELSE
            resid = max( resid, ( diffnm / xnorm )*rcond )
         END IF
   20 CONTINUE
      IF( resid*eps.LT.1.0d0 )
     $   resid = resid / eps
*
      RETURN
*
*     End of DGET04
*
      END

      SUBROUTINE cget04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDXACT, N, NRHS
      REAL               RCOND, RESID
*     ..
*     .. Array Arguments ..
      COMPLEX            X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      parameter( zero = 0.0e+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, J
      REAL               DIFFNM, EPS, XNORM
      COMPLEX            ZDUM
*     ..
*     .. External Functions ..
      INTEGER            ICAMAX
      REAL               SLAMCH
      EXTERNAL           icamax, slamch
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, aimag, max, real
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( zdum ) = abs( real( zdum ) ) + abs( aimag( zdum ) )
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( n.LE.0 .OR. nrhs.LE.0 ) THEN
         resid = zero
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if RCOND is invalid.
*
      eps = slamch( 'Epsilon' )
      IF( rcond.LT.zero ) THEN
         resid = 1.0 / eps
         RETURN
      END IF
*
*     Compute the maximum of
*        norm(X - XACT) / ( norm(XACT) * EPS )
*     over all the vectors X and XACT .
*
      resid = zero
      DO 20 j = 1, nrhs
         ix = icamax( n, xact( 1, j ), 1 )
         xnorm = cabs1( xact( ix, j ) )
         diffnm = zero
         DO 10 i = 1, n
            diffnm = max( diffnm, cabs1( x( i, j )-xact( i, j ) ) )
   10    CONTINUE
         IF( xnorm.LE.zero ) THEN
            IF( diffnm.GT.zero )
     $         resid = 1.0 / eps
         ELSE
            resid = max( resid, ( diffnm / xnorm )*rcond )
         END IF
   20 CONTINUE
      IF( resid*eps.LT.1.0 )
     $   resid = resid / eps
*
      RETURN
*
*     End of CGET04
*
      END

      SUBROUTINE zget04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDXACT, N, NRHS
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, J
      DOUBLE PRECISION   DIFFNM, EPS, XNORM
      COMPLEX*16         ZDUM
*     ..
*     .. External Functions ..
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           izamax, dlamch
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dimag, max
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( n.LE.0 .OR. nrhs.LE.0 ) THEN
         resid = zero
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if RCOND is invalid.
*
      eps = dlamch( 'Epsilon' )
      IF( rcond.LT.zero ) THEN
         resid = 1.0d0 / eps
         RETURN
      END IF
*
*     Compute the maximum of
*        norm(X - XACT) / ( norm(XACT) * EPS )
*     over all the vectors X and XACT .
*
      resid = zero
      DO 20 j = 1, nrhs
         ix = izamax( n, xact( 1, j ), 1 )
         xnorm = cabs1( xact( ix, j ) )
         diffnm = zero
         DO 10 i = 1, n
            diffnm = max( diffnm, cabs1( x( i, j )-xact( i, j ) ) )
   10    CONTINUE
         IF( xnorm.LE.zero ) THEN
            IF( diffnm.GT.zero )
     $         resid = 1.0d0 / eps
         ELSE
            resid = max( resid, ( diffnm / xnorm )*rcond )
         END IF
   20 CONTINUE
      IF( resid*eps.LT.1.0d0 )
     $   resid = resid / eps
*
      RETURN
*
*     End of ZGET04
*
      END

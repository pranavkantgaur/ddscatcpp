#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"

//!----------- cgcommon package ------------------------------------------
//! these routines are required for use by the conjugate gradient machiner
//! history
//! 03.01.28 (BTD) moved routines CAXPY and CSWAP to blas.f package of
//!                basic linear algebra subroutines
//! end history
//!-----------------------------------------------------------------------

// ChB: Scsetrhsstop moved to special Cgcommon2.cpp in Pimlib
// ChB: Smachcons, Stopcrit, Progress too

void Pcsum(int isize, Complex *x)
{
// ChB: it was empty
}


/* **
    SUBROUTINE CVPROD(N,CX,INCX,CY,INCY)
      USE DDPRECISION, ONLY : WP
! part of PIM
      IMPLICIT NONE

!     Modified from saxpy level 1 BLAS
!     element-wise vector multiplication, y<-x*y
!     Rudnei Dias da Cunha, 16/6/93

!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.


!     .. Scalar Arguments ..
      INTEGER :: INCX, INCY, N
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: CX(*), CY(*)
!     ..
!     .. Local Scalars ..
      INTEGER :: I, IX, IY, M, MP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX + 1
      IF (INCY<0) IY = (-N+1)*INCY + 1
      DO I = 1, N
        CY(IY) = CY(IY)*CX(IX)
        IX = IX + INCX
        IY = IY + INCY
      END DO
      RETURN

!        code for both increments equal to 1


!        clean-up loop

20    M = MOD(N,4)
      IF (M==0) GO TO 40
      DO I = 1, M
        CY(I) = CY(I)*CX(I)
      END DO
      IF (N<4) RETURN
40    MP1 = M + 1
      DO I = MP1, N, 4
        CY(I) = CY(I)*CX(I)
        CY(I+1) = CY(I+1)*CX(I+1)
        CY(I+2) = CY(I+2)*CX(I+2)
        CY(I+3) = CY(I+3)*CX(I+3)
      END DO
      RETURN

    END SUBROUTINE CVPROD
    SUBROUTINE PRINTV(N,U)
      USE DDPRECISION, ONLY : WP
! part of PIM
      IMPLICIT NONE
!     .. Scalar Arguments ..
      INTEGER :: N
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: U(*)
!     ..
!     .. Local Scalars ..
      INTEGER :: I
!     ..
      WRITE (6,FMT=9000) (U(I),I=1,N)
      RETURN
9000  FORMAT (8(E14.8,1X))
    END SUBROUTINE PRINTV
** */

/* **
    FUNCTION CSIGN(X)
      USE DDPRECISION, ONLY : WP
! part of PIM
      IMPLICIT NONE
      COMPLEX (WP) :: CSIGN
!     .. Scalar Arguments ..
      COMPLEX (WP) :: X
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
      CSIGN = X/ABS(X)
      RETURN

    END FUNCTION CSIGN

    SUBROUTINE DECODE(RHO,C,S)
      USE DDPRECISION, ONLY : WP
! part of PIM
      IMPLICIT NONE
!     .. Scalar Arguments ..
      COMPLEX (WP) :: C, RHO, S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
!     ..
!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0_WP,0.0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0_WP,0.0_WP))
!     ..
      IF (RHO==CONE) THEN
        C = CZERO
        S = CONE

      ELSE IF (ABS(RHO)<ONE) THEN
        S = 2.0_WP*RHO
        C = SQRT(CONE-S**2)

      ELSE
        C = 2.0_WP/RHO
        S = SQRT(CONE-C**2)
      END IF

      RETURN

    END SUBROUTINE DECODE

    SUBROUTINE ENCODE(RHO,C,S)
      USE DDPRECISION, ONLY : WP
! part of PIM
      IMPLICIT NONE
!     .. Scalar Arguments ..
      COMPLEX (WP) :: C, RHO, S
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CSIGN
      EXTERNAL CSIGN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Parameters ..
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0_WP,0.0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0_WP,0.0_WP))
!     ..
      IF (C==CZERO) THEN
        RHO = CONE

      ELSE IF (ABS(S)<ABS(C)) THEN
        RHO = CSIGN(C)*S/2.0_WP

      ELSE
        RHO = 2.0_WP*CSIGN(S)/C
      END IF

      RETURN

    END SUBROUTINE ENCODE

    SUBROUTINE GIVENS(A,B,C,S)
      USE DDPRECISION, ONLY : WP
! part of PIM
      IMPLICIT NONE
!     .. Scalar Arguments ..
      COMPLEX (WP) :: A, B, C, S
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: TAU
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
!     ..
!     .. Parameters ..
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0_WP,0.0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0_WP,0.0_WP))
!     ..
      IF (B==CZERO) THEN
        C = CONE
        S = CZERO

      ELSE IF (ABS(B)>ABS(A)) THEN
        TAU = -A/B
        S = CONE/SQRT(CONE+TAU**2)
        C = S*TAU

      ELSE
        TAU = -B/A
        C = CONE/SQRT(CONE+TAU**2)
        S = C*TAU
      END IF

      RETURN

    END SUBROUTINE GIVENS
** */

/* **

    SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)
      USE DDPRECISION, ONLY : WP

!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 4/11/78.

      COMPLEX (WP) :: CX(1), CY(1)
      INTEGER :: I, INCX, INCY, IX, IY, N

      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX + 1
      IF (INCY<0) IY = (-N+1)*INCY + 1
      DO I = 1, N
        CY(IY) = CX(IX)
        IX = IX + INCX
        IY = IY + INCY
      END DO
      RETURN

!        code for both increments equal to 1

20    DO I = 1, N
        CY(I) = CX(I)
      END DO
      RETURN
    END SUBROUTINE CCOPY
    FUNCTION CDOTC(N,CX,INCX,CY,INCY)
      USE DDPRECISION, ONLY : WP
      COMPLEX (WP) :: CDOTC

!     forms the dot product of a vector.
!     jack dongarra, 3/11/78.

      COMPLEX (WP) :: CX(1), CY(1), CTEMP
      INTEGER :: I, IX, IY, N, INCX, INCY

      CTEMP = (0.0E0_WP,0.0E0_WP)
      CDOTC = (0.0E0_WP,0.0E0_WP)
      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX + 1
      IF (INCY<0) IY = (-N+1)*INCY + 1
      DO I = 1, N
        CTEMP = CTEMP + CONJG(CX(IX))*CY(IY)
        IX = IX + INCX
        IY = IY + INCY
      END DO
      CDOTC = CTEMP
      RETURN

!        code for both increments equal to 1

20    DO I = 1, N
        CTEMP = CTEMP + CONJG(CX(I))*CY(I)
      END DO
      CDOTC = CTEMP
      RETURN
    END FUNCTION CDOTC
    FUNCTION CDOTU(N,CX,INCX,CY,INCY)
      USE DDPRECISION, ONLY : WP
      COMPLEX (WP) :: CDOTU

!     forms the dot product of a vector.
!     jack dongarra, 3/11/78.

      COMPLEX (WP) :: CX(1), CY(1), CTEMP
      INTEGER :: I, IX, IY, N, INCX, INCY

      CTEMP = (0.0E0_WP,0.0E0_WP)
      CDOTU = (0.0E0_WP,0.0E0_WP)
      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX + 1
      IF (INCY<0) IY = (-N+1)*INCY + 1
      DO I = 1, N
        CTEMP = CTEMP + CX(IX)*CY(IY)
        IX = IX + INCX
        IY = IY + INCY
      END DO
      CDOTU = CTEMP
      RETURN

!        code for both increments equal to 1

20    DO I = 1, N
        CTEMP = CTEMP + CX(I)*CY(I)
      END DO
      CDOTU = CTEMP
      RETURN
    END FUNCTION CDOTU
    SUBROUTINE CROTG(CA,CB,C,S)
      USE DDPRECISION, ONLY : WP
      COMPLEX (WP) :: CA, CB, S
      REAL (WP) :: C
      REAL (WP) :: NORM, SCALE
      COMPLEX (WP) :: ALPHA

      IF (ABS(CA)/=0.0E0_WP) GO TO 10
      C = 0.0E0_WP
      S = (1.0E0_WP,0.0E0_WP)
      CA = CB
      GO TO 20
10    CONTINUE
      SCALE = ABS(CA) + ABS(CB)
      NORM = SCALE*SQRT((ABS(CA/CMPLX(SCALE,0.0E0_WP,KIND=WP)))**2+(ABS(CB/ &
        CMPLX(SCALE,0.0E0_WP,KIND=WP)))**2)
      ALPHA = CA/ABS(CA)
      C = ABS(CA)/NORM
      S = ALPHA*CONJG(CB)/NORM
      CA = ALPHA*NORM
20    CONTINUE
      RETURN
    END SUBROUTINE CROTG
    SUBROUTINE CSROT(N,CX,INCX,CY,INCY,C,S)
      USE DDPRECISION, ONLY : WP

!     applies a plane rotation, where the cos and sin (c and s) are
!     real and the vectors cx and cy are Complex.
!     jack dongarra, linpack, 3/11/78.

      COMPLEX (WP) :: CX(1), CY(1), CTEMP
      REAL (WP) :: C, S
      INTEGER :: I, INCX, INCY, IX, IY, N

      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!       code for unequal increments or equal increments not equal
!         to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX + 1
      IF (INCY<0) IY = (-N+1)*INCY + 1
      DO I = 1, N
        CTEMP = C*CX(IX) + S*CY(IY)
        CY(IY) = C*CY(IY) - S*CX(IX)
        CX(IX) = CTEMP
        IX = IX + INCX
        IY = IY + INCY
      END DO
      RETURN

!       code for both increments equal to 1

20    DO I = 1, N
        CTEMP = C*CX(I) + S*CY(I)
        CY(I) = C*CY(I) - S*CX(I)
        CX(I) = CTEMP
      END DO
      RETURN
    END SUBROUTINE CSROT
** */

/* **
!-----------------------------------------------------------------------
! Following code due to C.L. Lawson has been replaced by above module
! because it generated bad code when compiled with g77.
! For use by DDSCAT the overflow/underflow avoidance strategies
! used by this routine do not appear to be necessary, so they can
! be omitted in interests of speed and clean code.
! 98.10.07 BTD

!      real function scnrm2( n, cx, incx)
!      logical imag, scale
!      int i, incx, ix, n, next
!      real cutlo, cuthi, hitest, sum, xmax, absx, zero, one
!      Complex      cx(1)
!      real real,aimag
!      Complex zdumr,zdumi
!      real(zdumr) = zdumr
!      aimag(zdumi) = (0.0e0,-1.0e0)*zdumi
!      data         zero, one /0.0e0, 1.0e0/
!c
!c     unitary norm of the Complex n-vector stored in cx() with storage
!c     increment incx .
!c     if    n .le. 0 return with result = 0.
!c     if n .ge. 1 then incx must be .ge. 1
!c
!c           c.l.lawson , 1978 jan 08
!c     modified to correct problem with negative increment, 8/21/90.
!c
!c     four phase method     using two built-in constants that are
!c     hopefully applicable to all machines.
!c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
!c         cuthi = minimum of  sqrt(v)      over all known machines.
!c     where
!c         eps = smallest no. such that eps + 1. .gt. 1.
!c         u   = smallest positive no.   (underflow limit)
!c         v   = largest  no.            (overflow  limit)
!c
!c     brief outline of algorithm..
!c
!c     phase 1    scans zero components.
!c     move to phase 2 when a component is nonzero and .le. cutlo
!c     move to phase 3 when a component is .gt. cutlo
!c     move to phase 4 when a component is .ge. cuthi/m
!c     where m = n for x() real and m = 2*n for Complex.
!c
!c     values for cutlo and cuthi..
!c     from the environmental parameters listed in the imsl converter
!c     document the limiting values are as follows..
!c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds ar
!c                   univac and dec at 2**(-103)
!c                   thus cutlo = 2**(-51) = 4.44089e-16
!c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
!c                   thus cuthi = 2**(63.5) = 1.30438e19
!c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
!c                   thus cutlo = 2**(-33.5) = 8.23181d-11
!c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
!c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
!c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
!      data cutlo, cuthi / 8.232d-11,  1.304d19 /
!c
!      if(n .gt. 0) go to 10
!         scnrm2  = zero
!         go to 300
!c
!   10 assign 30 to next
!      sum = zero
!      i = 1
!      if( incx .lt. 0 )i = (-n+1)*incx + 1
!c                                                 begin main loop
!      do 220 ix = 1,n
!         absx = abs(real(cx(i)))
!         imag = .false.
!         go to next,(30, 50, 70, 90, 110)
!   30 if( absx .gt. cutlo) go to 85
!      assign 50 to next
!      scale = .false.
!c
!c                        phase 1.  sum is zero
!c
!   50 if( absx .eq. zero) go to 200
!      if( absx .gt. cutlo) go to 85
!c
!c                                prepare for phase 2.
!      assign 70 to next
!      go to 105
!c
!c                                prepare for phase 4.
!c
!  100 assign 110 to next
!      sum = (sum / absx) / absx
!  105 scale = .true.
!      xmax = absx
!      go to 115
!c
!c                   phase 2.  sum is small.
!c                             scale to avoid destructive underflow.
!c
!   70 if( absx .gt. cutlo ) go to 75
!c
!c                     common code for phases 2 and 4.
!c                     in phase 4 sum is large.  scale to avoid overflow
!c
!  110 if( absx .le. xmax ) go to 115
!         sum = one + sum * (xmax / absx)**2
!         xmax = absx
!         go to 200
!c
!  115 sum = sum + (absx/xmax)**2
!      go to 200
!c
!c
!c                  prepare for phase 3.
!c
!   75 sum = (sum * xmax) * xmax
!c
!   85 assign 90 to next
!      scale = .false.
!c
!c     for real or d.p. set hitest = cuthi/n
!c     for Complex      set hitest = cuthi/(2*n)
!c
!      hitest = cuthi/float( 2*n )
!c
!c                   phase 3.  sum is mid-range.  no scaling.
!c
!   90 if(absx .ge. hitest) go to 100
!         sum = sum + absx**2
!  200 continue
!c                  control selection of real and imaginary parts.
!c
!      if(imag) go to 210
!         absx = abs(aimag(cx(i)))
!         imag = .true.
!      go to next,(  50, 70, 90, 110 )
!c
!  210 continue
!      i = i + incx
!  220 continue
!c
!c              end of main loop.
!c              compute square root and adjust for scaling.
!c
!      scnrm2 = sqrt(sum)
!      if(scale) scnrm2 = scnrm2 * xmax
!  300 continue
!      return
!      end
!-----------------------------------------------------------------------
    SUBROUTINE SROTG(DA,DB,C,S)
      USE DDPRECISION, ONLY : WP

!     construct givens plane rotation.
!     jack dongarra, linpack, 3/11/78.
!                    modified 9/27/86.

      REAL (WP) :: DA, DB, C, S, ROE, SCALE, R, Z

      ROE = DB
      IF (ABS(DA)>ABS(DB)) ROE = DA
      SCALE = ABS(DA) + ABS(DB)
      IF (SCALE/=0.0E0_WP) GO TO 10
      C = 1.0E0_WP
      S = 0.0E0_WP
      R = 0.0E0_WP
      GO TO 20
10    R = SCALE*SQRT((DA/SCALE)**2+(DB/SCALE)**2)
      R = SIGN(1.0E0_WP,ROE)*R
      C = DA/R
      S = DB/R
20    Z = S
      IF (ABS(C)>0.0E0_WP .AND. ABS(C)<=S) Z = 1.0E0_WP/C
      DA = R
      DB = Z
      RETURN
    END SUBROUTINE SROTG

    SUBROUTINE DMACHCONS(WHAT,RESULT)
      USE DDPRECISION, ONLY : WP

! These values are for IEEE-754 arithmetic
! 07.08.05 (BTD) Modified for compatibility with both single or
!                double precision -- eliminate PARAMETER statements,
!                add executable statements.
!                Will these slow execution?? 
!                How often is DMACHCONS called?
! 08.08.05 (BTD) Modified to use preprocessing to modify file
!                to use appropriate precision MACHEPS,OVERFlOW,
!                end UNDERFLOW
! end history --------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: MACHEPS,OVERFLOW,UNDERFLOW
#ifdef sp
      PARAMETER(MACHEPS=1.192093E-7_WP)
      PARAMETER(OVERFLOW=3.402823E+38_WP)
      PARAMETER(UNDERFLOW=1.17549435E-38_WP)
#endif
#ifdef dp
      PARAMETER(MACHEPS=2.2204460492503E-16_WP)
      PARAMETER(OVERFLOW=1.7976313E+308_WP)
      PARAMETER(UNDERFLOW=2.2250739E-308_WP)
#endif

!     .. Scalar Arguments ..
      REAL (WP) :: RESULT
      CHARACTER :: WHAT

      IF((WHAT=='M').OR.(WHAT=='m'))THEN
         RESULT=MACHEPS
      ELSEIF((WHAT=='U').OR.(WHAT=='u'))THEN
         RESULT=UNDERFLOW
      ELSEIF((WHAT=='O').OR.(WHAT=='o'))THEN
         RESULT=OVERFLOW
      ENDIF

      RETURN

    END SUBROUTINE DMACHCONS
    SUBROUTINE PIMDGETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)
      USE DDPRECISION, ONLY : WP

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!               Centro de Processamento de Dados,
!         Universidade Federal do Rio Grande do Sul, Brasil
!                              and
!     Computing Laboratory, University of Kent at Canterbury, U.K.

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!  Description of parameter arrays
!   IPAR (INPUT)  : int
!     ipar( 1): lda    (Leading dimension of a)
!           2 : n      (Number of rows/columns of a)
!           3 : blksz  (Size of block of data; used when data is
!                       partitioned using cyclic mode)
!           4 : loclen (Number of elements stored locally;
!                       *PARALLEL: Equal to at least m/nprocs or
!                                  n/procs depending if row or
!                                  column partitioning is used or,
!                                  in the case of cyclic partitioning,
!                                  it is a multiple of either
!                                  m/(blksz*nprocs) or n/(blksz*nprocs).
!                       *SEQUENTIAL: equal to n)
!           5 : basisdim (Dimension of orthogonal basis, used in
!                       GMRES)
!           6 : nprocs (Number of processors)
!           7 : procid (Processor identification)
!           8 : precontype (Type of preconditioning; one of
!                           0 : No preconditioning,
!                           1 : Left preconditioning,
!                           2 : Right preconditioning,
!                           3 : Symmetric preconditioning)
!           9 : stoptype (Type of stopping criteria used)
!          10 : maxit  (Maximum number of iterations allowed)

!   IPAR (OUTPUT) : int
!     ipar(11): itno   (Number of iterations executed)
!          12 : status (On exit of iterative method, one of
!                        0: converged to solution
!                       -1: no convergence has been achieved
!                       -2: "soft"-breakdown, solution may have
!                           been found
!                       -3: "hard"-breakdown, no solution)
!                       -4: conflict in preconditioner and stopping
!                           criterion selected
!                       -5: error in stopping criterion 3, r^{T}z<0)
!          13 : steperr (If status is either -2 or -3, it gives
!                         the step number of the respective algorithm
!                         where a breakdown has occurred. Refer to the
!                         User's Guide for further information)

!   RPAR/DPAR (INPUT)  : real/double precision
!     dpar( 1): epsilon (The value of epsilon for use in the
!                        stopping criterion)

!   RPAR/DPAR (OUTPUT) : real/double precision
!     dpar( 2): exitnorm (The value of a norm of either the residual
!                         vector or the difference between two
!                         successive solution estimates according to
!                         the value of stoptype)
!           3,
!           4 : smallest and largest eigenvalues of Q1AQ2 (in the
!               symmetric case) OR smallest and largest real parts
!               (in the nonsymmetric case)
!           5,
!           6 : smallest and largest imaginary parts (only in the
!               nonsymmetric case)


!     .. Parameters ..
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: DPARSIZ
      PARAMETER (DPARSIZ=6)
!     ..
!     .. Scalar Arguments ..
      REAL (WP) :: EPSILON, EXITNORM
      INTEGER :: BASISDIM, BLKSZ, ITNO, LDA, LOCLEN, MAXIT, N, NPROCS, &
        PRECONTYPE, PROCID, STATUS, STEPERR, STOPTYPE
!     ..
!     .. Array Arguments ..
      REAL (WP) :: DPAR(DPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
      LDA = IPAR(1)
      N = IPAR(2)
      BLKSZ = IPAR(3)
      LOCLEN = IPAR(4)
      BASISDIM = IPAR(5)
      NPROCS = IPAR(6)
      PROCID = IPAR(7)
      PRECONTYPE = IPAR(8)
      STOPTYPE = IPAR(9)
      MAXIT = IPAR(10)
      ITNO = IPAR(11)
      STATUS = IPAR(12)
      STEPERR = IPAR(13)

      EPSILON = DPAR(1)
      EXITNORM = DPAR(2)

      RETURN

    END SUBROUTINE PIMDGETPAR
    SUBROUTINE PIMDPRTPAR(IPAR,DPAR)
      USE DDPRECISION, ONLY : WP

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!               Centro de Processamento de Dados,
!         Universidade Federal do Rio Grande do Sul, Brasil
!                              and
!     Computing Laboratory, University of Kent at Canterbury, U.K.

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!  Description of parameter arrays
!   IPAR (INPUT)  : int
!     ipar( 1): lda    (Leading dimension of a)
!           2 : n      (Number of rows/columns of a)
!           3 : blksz  (Size of block of data; used when data is
!                       partitioned using cyclic mode)
!           4 : loclen (Number of elements stored locally;
!                       *PARALLEL: Equal to at least m/nprocs or
!                                  n/procs depending if row or
!                                  column partitioning is used or,
!                                  in the case of cyclic partitioning,
!                                  it is a multiple of either
!                                  m/(blksz*nprocs) or n/(blksz*nprocs).
!                       *SEQUENTIAL: equal to n)
!           5 : basisdim (Dimension of orthogonal basis, used in
!                       GMRES)
!           6 : nprocs (Number of processors)
!           7 : procid (Processor identification)
!           8 : precontype (Type of preconditioning; one of
!                           0 : No preconditioning,
!                           1 : Left preconditioning,
!                           2 : Right preconditioning,
!                           3 : Symmetric preconditioning)
!           9 : stoptype (Type of stopping criteria used)
!          10 : maxit  (Maximum number of iterations allowed)

!   IPAR (OUTPUT) : int
!     ipar(11): itno   (Number of iterations executed)
!          12 : status (On exit of iterative method, one of
!                        0: converged to solution
!                       -1: no convergence has been achieved
!                       -2: "soft"-breakdown, solution may have
!                           been found
!                       -3: "hard"-breakdown, no solution)
!                       -4: conflict in preconditioner and stopping
!                           criterion selected
!                       -5: error in stopping criterion 3, r^{T}z<0)
!          13 : steperr (If status is either -2 or -3, it gives
!                         the step number of the respective algorithm
!                         where a breakdown has occurred. Refer to the
!                         User's Guide for further information)

!   RPAR/DPAR (INPUT)  : real/double precision
!     dpar( 1): epsilon (The value of epsilon for use in the
!                        stopping criterion)

!   RPAR/DPAR (OUTPUT) : real/double precision
!     dpar( 2): exitnorm (The value of a norm of either the residual
!                         vector or the difference between two
!                         successive solution estimates according to
!                         the value of stoptype)
!           3,
!           4 : smallest and largest eigenvalues of Q1AQ2 (in the
!               symmetric case) OR smallest and largest real parts
!               (in the nonsymmetric case)
!           5,
!           6 : smallest and largest imaginary parts (only in the
!               nonsymmetric case)


!     .. Parameters ..
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: DPARSIZ
      PARAMETER (DPARSIZ=6)
!     ..
!     .. Array Arguments ..
      REAL (WP) :: DPAR(DPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Local Scalars ..
      INTEGER :: I
!     ..
      WRITE (6,FMT=10) (IPAR(I),I=1,IPARSIZ)

10    FORMAT ('lda=',I6/'n=',I6/'blksz=',I6/'loclen=',I4/'basisdim=', &
        I4/'nprocs=',I4/'procid=',I4/'precontype=',I4/'stoptype=',I4/'maxit=', &
        I4/'itno=',I4/'status=',I4/'steperr=',I4)

      WRITE (6,FMT=20) (DPAR(I),I=1,DPARSIZ)

20    FORMAT ('epsilon=',D20.12/'exitnorm=',D20.12/'eigenvalues region:'/4( &
        D20.12,1X))

      RETURN

    END SUBROUTINE PIMDPRTPAR
    SUBROUTINE PIMDSETPAR(IPAR,DPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,EPSILON)
      USE DDPRECISION, ONLY : WP

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!               Centro de Processamento de Dados,
!         Universidade Federal do Rio Grande do Sul, Brasil
!                              and
!     Computing Laboratory, University of Kent at Canterbury, U.K.

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!  Description of parameter arrays
!   IPAR (INPUT)  : int
!     ipar( 1): lda    (Leading dimension of a)
!           2 : n      (Number of rows/columns of a)
!           3 : blksz  (Size of block of data; used when data is
!                       partitioned using cyclic mode)
!           4 : loclen (Number of elements stored locally;
!                       *PARALLEL: Equal to at least m/nprocs or
!                                  n/procs depending if row or
!                                  column partitioning is used or,
!                                  in the case of cyclic partitioning,
!                                  it is a multiple of either
!                                  m/(blksz*nprocs) or n/(blksz*nprocs).
!                       *SEQUENTIAL: equal to n)
!           5 : basisdim (Dimension of orthogonal basis, used in
!                       GMRES)
!           6 : nprocs (Number of processors)
!           7 : procid (Processor identification)
!           8 : precontype (Type of preconditioning; one of
!                           0 : No preconditioning,
!                           1 : Left preconditioning,
!                           2 : Right preconditioning,
!                           3 : Symmetric preconditioning)
!           9 : stoptype (Type of stopping criteria used)
!          10 : maxit  (Maximum number of iterations allowed)

!   IPAR (OUTPUT) : int
!     ipar(11): itno   (Number of iterations executed)
!          12 : status (On exit of iterative method, one of
!                        0: converged to solution
!                       -1: no convergence has been achieved
!                       -2: "soft"-breakdown, solution may have
!                           been found
!                       -3: "hard"-breakdown, no solution)
!                       -4: conflict in preconditioner and stopping
!                           criterion selected
!                       -5: error in stopping criterion 3, r^{T}z<0)
!          13 : steperr (If status is either -2 or -3, it gives
!                         the step number of the respective algorithm
!                         where a breakdown has occurred. Refer to the
!                         User's Guide for further information)

!   RPAR/DPAR (INPUT)  : real/double precision
!     dpar( 1): epsilon (The value of epsilon for use in the
!                        stopping criterion)

!   RPAR/DPAR (OUTPUT) : real/double precision
!     dpar( 2): exitnorm (The value of a norm of either the residual
!                         vector or the difference between two
!                         successive solution estimates according to
!                         the value of stoptype)
!           3,
!           4 : smallest and largest eigenvalues of Q1AQ2 (in the
!               symmetric case) OR smallest and largest real parts
!               (in the nonsymmetric case)
!           5,
!           6 : smallest and largest imaginary parts (only in the
!               nonsymmetric case)


!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: DPARSIZ
      PARAMETER (DPARSIZ=6)
!     ..
!     .. Scalar Arguments ..
      REAL (WP) :: EPSILON
      INTEGER :: BASISDIM, BLKSZ, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, &
        PROCID, STOPTYPE
!     ..
!     .. Array Arguments ..
      REAL (WP) :: DPAR(DPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Local Scalars ..
      REAL (WP) :: EXITNORM
      INTEGER :: ITNO, STATUS, STEPERR
!     ..
      IPAR(1) = LDA
      IPAR(2) = N
      IPAR(3) = BLKSZ
      IPAR(4) = LOCLEN
      IPAR(5) = BASISDIM
      IPAR(6) = NPROCS
      IPAR(7) = PROCID
      IPAR(8) = PRECONTYPE
      IPAR(9) = STOPTYPE
      IPAR(10) = MAXIT
      IPAR(11) = -1
      IPAR(12) = -1
      IPAR(13) = -1

      DPAR(1) = EPSILON
      DPAR(2) = -ONE

      RETURN

    END SUBROUTINE PIMDSETPAR
    SUBROUTINE PIMSPRTPAR(IPAR,SPAR)
      USE DDPRECISION, ONLY : WP

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!               Centro de Processamento de Dados,
!         Universidade Federal do Rio Grande do Sul, Brasil
!                              and
!     Computing Laboratory, University of Kent at Canterbury, U.K.

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!  Description of parameter arrays
!   IPAR (INPUT)  : int
!     ipar( 1): lda    (Leading dimension of a)
!           2 : n      (Number of rows/columns of a)
!           3 : blksz  (Size of block of data; used when data is
!                       partitioned using cyclic mode)
!           4 : loclen (Number of elements stored locally;
!                       *PARALLEL: Equal to at least m/nprocs or
!                                  n/procs depending if row or
!                                  column partitioning is used or,
!                                  in the case of cyclic partitioning,
!                                  it is a multiple of either
!                                  m/(blksz*nprocs) or n/(blksz*nprocs).
!                       *SEQUENTIAL: equal to n)
!           5 : basisdim (Dimension of orthogonal basis, used in
!                       GMRES)
!           6 : nprocs (Number of processors)
!           7 : procid (Processor identification)
!           8 : precontype (Type of preconditioning; one of
!                           0 : No preconditioning,
!                           1 : Left preconditioning,
!                           2 : Right preconditioning,
!                           3 : Symmetric preconditioning)
!           9 : stoptype (Type of stopping criteria used)
!          10 : maxit  (Maximum number of iterations allowed)

!   IPAR (OUTPUT) : int
!     ipar(11): itno   (Number of iterations executed)
!          12 : status (On exit of iterative method, one of
!                        0: converged to solution
!                       -1: no convergence has been achieved
!                       -2: "soft"-breakdown, solution may have
!                           been found
!                       -3: "hard"-breakdown, no solution)
!                       -4: conflict in preconditioner and stopping
!                           criterion selected
!                       -5: error in stopping criterion 3, r^{T}z<0)
!          13 : steperr (If status is either -2 or -3, it gives
!                         the step number of the respective algorithm
!                         where a breakdown has occurred. Refer to the
!                         User's Guide for further information)

!   RPAR/DPAR (INPUT)  : real/real
!     spar( 1): epsilon (The value of epsilon for use in the
!                        stopping criterion)

!   RPAR/DPAR (OUTPUT) : real/real
!     spar( 2): exitnorm (The value of a norm of either the residual
!                         vector or the difference between two
!                         successive solution estimates according to
!                         the value of stoptype)
!           3,
!           4 : smallest and largest eigenvalues of Q1AQ2 (in the
!               symmetric case) OR smallest and largest real parts
!               (in the nonsymmetric case)
!           5,
!           6 : smallest and largest imaginary parts (only in the
!               nonsymmetric case)


!     .. Parameters ..
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Local Scalars ..
      INTEGER :: I
!     ..
      WRITE (6,FMT=10) (IPAR(I),I=1,IPARSIZ)

10    FORMAT ('lda=',I6/'n=',I6/'blksz=',I6/'loclen=',I4/'basisdim=', &
        I4/'nprocs=',I4/'procid=',I4/'precontype=',I4/'stoptype=',I4/'maxit=', &
        I4/'itno=',I4/'status=',I4/'steperr=',I4)

      WRITE (6,FMT=20) (SPAR(I),I=1,SPARSIZ)

20    FORMAT ('epsilon=',E20.12/'exitnorm=',E20.12,'eigenvalues region:'/4( &
        E20.12,1X))


      RETURN

    END SUBROUTINE PIMSPRTPAR
** */

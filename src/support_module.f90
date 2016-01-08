!*******************************************************************************
    module support_module

    implicit none

    contains
!*******************************************************************************

!## Following a selection from BLAS Level 1

      SUBROUTINE DAXPY(N,Da,Dx,Incx,Dy,Incy)
      IMPLICIT NONE

!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      DOUBLE PRECISION Dx(*) , Dy(*) , Da
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N

      IF ( N<=0 ) RETURN
      IF ( Da==0.0D0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1

!        CLEAN-UP LOOP

         m = MOD(N,4)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Dy(i) = Dy(i) + Da*Dx(i)
            ENDDO
            IF ( N<4 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 4
            Dy(i) = Dy(i) + Da*Dx(i)
            Dy(i+1) = Dy(i+1) + Da*Dx(i+1)
            Dy(i+2) = Dy(i+2) + Da*Dx(i+2)
            Dy(i+3) = Dy(i+3) + Da*Dx(i+3)
         ENDDO
      ELSE

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!        NOT EQUAL TO 1

         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            Dy(iy) = Dy(iy) + Da*Dx(ix)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
         RETURN
      ENDIF
      END SUBROUTINE DAXPY

      SUBROUTINE DCOPY(N,Dx,Incx,Dy,Incy)
      IMPLICIT NONE

!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      DOUBLE PRECISION Dx(*) , Dy(*)
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N

      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1

!        CLEAN-UP LOOP

         m = MOD(N,7)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Dy(i) = Dx(i)
            ENDDO
            IF ( N<7 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 7
            Dy(i) = Dx(i)
            Dy(i+1) = Dx(i+1)
            Dy(i+2) = Dx(i+2)
            Dy(i+3) = Dx(i+3)
            Dy(i+4) = Dx(i+4)
            Dy(i+5) = Dx(i+5)
            Dy(i+6) = Dx(i+6)
         ENDDO
      ELSE

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!        NOT EQUAL TO 1

         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            Dy(iy) = Dx(ix)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
         RETURN
      ENDIF
      END SUBROUTINE DCOPY

      DOUBLE PRECISION FUNCTION DDOT(N,Dx,Incx,Dy,Incy)
      IMPLICIT NONE

!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      DOUBLE PRECISION Dx(*) , Dy(*) , dtemp
      INTEGER i , Incx , Incy , ix , iy , m , mp1 , N

      DDOT = 0.0D0
      dtemp = 0.0D0
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1

!        CLEAN-UP LOOP

         m = MOD(N,5)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               dtemp = dtemp + Dx(i)*Dy(i)
            ENDDO
            IF ( N<5 ) THEN
               DDOT = dtemp
               return
            ENDIF
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 5
            dtemp = dtemp + Dx(i)*Dy(i) + Dx(i+1)*Dy(i+1) + Dx(i+2)     &
                    *Dy(i+2) + Dx(i+3)*Dy(i+3) + Dx(i+4)*Dy(i+4)
         ENDDO
         DDOT = dtemp
      ELSE

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            dtemp = dtemp + Dx(ix)*Dy(iy)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
         DDOT = dtemp
         RETURN
      ENDIF
      END FUNCTION DDOT

      DOUBLE PRECISION FUNCTION DNRM1(N,X,I,J)
      IMPLICIT NONE
      INTEGER N , I , J , k
      DOUBLE PRECISION snormx , sum , X(N) , zero , one , scale , temp
      DATA zero/0.0D0/ , one/1.0D0/

!      DNRM1 - COMPUTES THE I-NORM OF A VECTOR
!              BETWEEN THE ITH AND THE JTH ELEMENTS

!      INPUT -
!      N       LENGTH OF VECTOR
!      X       VECTOR OF LENGTH N
!      I       INITIAL ELEMENT OF VECTOR TO BE USED
!      J       FINAL ELEMENT TO USE

!      OUTPUT -
!      DNRM1   NORM

      snormx = zero
      DO k = I , J
         snormx = MAX(snormx,ABS(X(k)))
      ENDDO
      DNRM1 = snormx
      IF ( snormx==zero ) RETURN
      scale = snormx
      IF ( snormx>=one ) scale = SQRT(snormx)
      sum = zero
      DO k = I , J
         temp = zero
         IF ( ABS(X(k))+scale/=scale ) temp = X(k)/snormx
         IF ( one+temp/=one ) sum = sum + temp*temp
      ENDDO
      sum = SQRT(sum)
      DNRM1 = snormx*sum
      END FUNCTION DNRM1

!replaced original routine with this one from BLAS:
!  http://netlib.sandia.gov/blas/dnrm2.f

        DOUBLE PRECISION FUNCTION DNRM2(N,X,Incx)
        IMPLICIT NONE
  !     .. Scalar Arguments ..
        INTEGER Incx , N
  !     ..
  !     .. Array Arguments ..
        DOUBLE PRECISION X(*)
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  DNRM2 returns the euclidean norm of a vector via the function
  !  name, so that
  !
  !     DNRM2 := sqrt( x'*x )
  !
  !  Further Details
  !  ===============
  !
  !  -- This version written on 25-October-1982.
  !     Modified on 14-October-1993 to inline the call to DLASSQ.
  !     Sven Hammarling, Nag Ltd.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
        DOUBLE PRECISION ONE , ZERO
        PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
  !     ..
  !     .. Local Scalars ..
        DOUBLE PRECISION absxi , norm , scale , ssq
        INTEGER ix
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC ABS , SQRT
  !     ..
        IF ( N<1 .OR. Incx<1 ) THEN
           norm = ZERO
        ELSEIF ( N==1 ) THEN
           norm = ABS(X(1))
        ELSE
           scale = ZERO
           ssq = ONE
  !        The following loop is equivalent to this call to the LAPACK
  !        auxiliary routine:
  !        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
  !
           DO ix = 1 , 1 + (N-1)*Incx , Incx
              IF ( X(ix)/=ZERO ) THEN
                 absxi = ABS(X(ix))
                 IF ( scale<absxi ) THEN
                    ssq = ONE + ssq*(scale/absxi)**2
                    scale = absxi
                 ELSE
                    ssq = ssq + (absxi/scale)**2
                 ENDIF
              ENDIF
           ENDDO
           norm = scale*SQRT(ssq)
        ENDIF
  !
        DNRM2 = norm
  !
  !     End of DNRM2.
  !
        END FUNCTION DNRM2

!      DOUBLE PRECISION FUNCTION DNRM2(N,Dx,Incx)
!      IMPLICIT NONE
!      INTEGER N , i , j , nn , next , Incx
!      DOUBLE PRECISION Dx(*) , cutlo , cuthi , hitest , sum , xmax ,    &
!                       zero , one
!      DATA zero , one/0.0D0 , 1.0D0/
!
!!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
!!     INCREMENT INCX .
!!     IF    N <= 0 RETURN WITH RESULT = 0.
!!     IF N >= 1 THEN INCX MUST BE >= 1
!
!!           C.L.LAWSON, 1978 JAN 08
!
!!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!!         CUTLO = MAXIMUM OF  SQRT(U/EPS)   OVER ALL KNOWN MACHINES.
!!         CUTHI = MINIMUM OF  SQRT(V)       OVER ALL KNOWN MACHINES.
!!     WHERE
!!         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!!         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
!!         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
!
!!     BRIEF OUTLINE OF ALGORITHM..
!
!!     PHASE 1    SCANS ZERO COMPONENTS.
!!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND <= CUTLO
!!     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!!     MOVE TO PHASE 4 WHEN A COMPONENT IS >= CUTHI/M
!!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!
!!     VALUES FOR CUTLO AND CUTHI..
!!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!!                   UNIVAC AND DEC AT 2**(-103)
!!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
!      DATA cutlo , cuthi/8.232D-11 , 1.304D19/
!
!      IF ( N>0 ) THEN
!
!         ASSIGN 200 TO next
!         sum = zero
!         nn = N*Incx
!!                       BEGIN MAIN LOOP
!         i = 1
!      ELSE
!         DNRM2 = zero
!         return
!      ENDIF
! 100  GOTO next
! 200  IF ( ABS(Dx(i))>cutlo ) GOTO 800
!      ASSIGN 300 TO next
!      xmax = zero
!
!!                        PHASE 1.  SUM IS ZERO
!
! 300  IF ( Dx(i)==zero ) GOTO 900
!      IF ( ABS(Dx(i))>cutlo ) GOTO 800
!
!!                        PREPARE FOR PHASE 2.
!
!      ASSIGN 600 TO next
!      GOTO 500
!
!!                        PREPARE FOR PHASE 4.
!
! 400  i = j
!      ASSIGN 700 TO next
!      sum = (sum/Dx(i))/Dx(i)
! 500  xmax = ABS(Dx(i))
!
!      sum = sum + (Dx(i)/xmax)**2
!      GOTO 900
!
!!                   PHASE 2.  SUM IS SMALL.
!!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
! 600  IF ( ABS(Dx(i))>cutlo ) THEN
!
!!                  PREPARE FOR PHASE 3.
!
!         sum = (sum*xmax)*xmax
!         GOTO 800
!      ENDIF
!
!!                   COMMON CODE FOR PHASES 2 AND 4.
!!                   IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
! 700  IF ( ABS(Dx(i))<=xmax ) THEN
!         sum = sum + (Dx(i)/xmax)**2
!      ELSE
!         sum = one + sum*(xmax/Dx(i))**2
!         xmax = ABS(Dx(i))
!      ENDIF
!      GOTO 900
!
!!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
! 800  hitest = cuthi/FLOAT(N)
!
!!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
!      DO j = i , nn , Incx
!         IF ( ABS(Dx(j))>=hitest ) GOTO 400
!         sum = sum + Dx(j)**2
!      ENDDO
!      DNRM2 = SQRT(sum)
!      return
!
! 900  i = i + Incx
!      IF ( i<=nn ) GOTO 100
!
!!              END OF MAIN LOOP.
!
!!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
!      DNRM2 = xmax*SQRT(sum)
!      END FUNCTION DNRM2

      SUBROUTINE DSROT(N,Dx,Incx,Dy,Incy,C,S)
      IMPLICIT NONE

!     APPLIES A PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.

      DOUBLE PRECISION Dx(*) , Dy(*) , dtemp , C , S
      INTEGER i , Incx , Incy , ix , iy , N

      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN

!       CODE FOR BOTH INCREMENTS EQUAL TO 1

         DO i = 1 , N
            dtemp = C*Dx(i) + S*Dy(i)
            Dy(i) = C*Dy(i) - S*Dx(i)
            Dx(i) = dtemp
         ENDDO
         return
      ENDIF

!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1

      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         dtemp = C*Dx(ix) + S*Dy(iy)
         Dy(iy) = C*Dy(iy) - S*Dx(ix)
         Dx(ix) = dtemp
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
      END SUBROUTINE DSROT

      SUBROUTINE DSROTG(Da,Db,C,S)
      IMPLICIT NONE

!     CONSTRUCT GIVENS PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.
!                    MODIFIED 9/27/86.

      DOUBLE PRECISION Da , Db , C , S , roe , scale , r , z , one ,    &
                       zero
      DATA one , zero/1.0D+00 , 0.0D+00/

      roe = Db
      IF ( ABS(Da)>ABS(Db) ) roe = Da
      scale = ABS(Da) + ABS(Db)
      IF ( scale/=zero ) THEN
         r = scale*SQRT((Da/scale)**2+(Db/scale)**2)
         r = SIGN(one,roe)*r
         C = Da/r
         S = Db/r
      ELSE
         C = one
         S = zero
         r = zero
      ENDIF
      z = S
      IF ( ABS(C)>zero .AND. ABS(C)<=S ) z = one/C
      Da = r
      Db = z
      END SUBROUTINE DSROTG

      SUBROUTINE DSCAL(N,Da,Dx,Incx)
      IMPLICIT NONE

!     SCALES A VECTOR BY A CONSTANT.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      DOUBLE PRECISION Da , Dx(*)
      INTEGER i , Incx , m , mp1 , N , nincx

      IF ( N<=0 ) RETURN
      IF ( Incx==1 ) THEN

!        CODE FOR INCREMENT EQUAL TO 1

!        CLEAN-UP LOOP

         m = MOD(N,5)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Dx(i) = Da*Dx(i)
            ENDDO
            IF ( N<5 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 5
            Dx(i) = Da*Dx(i)
            Dx(i+1) = Da*Dx(i+1)
            Dx(i+2) = Da*Dx(i+2)
            Dx(i+3) = Da*Dx(i+3)
            Dx(i+4) = Da*Dx(i+4)
         ENDDO
      ELSE

!        CODE FOR INCREMENT NOT EQUAL TO 1

         nincx = N*Incx
         DO i = 1 , nincx , Incx
            Dx(i) = Da*Dx(i)
         ENDDO

      ENDIF

      END SUBROUTINE DSCAL

!*******************************************************************************
    end module support_module
!*******************************************************************************

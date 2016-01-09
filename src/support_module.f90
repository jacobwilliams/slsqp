!*******************************************************************************
!>
!  Support routines for SLSQP. For example, a selection from BLAS level 1.
!  These have also been refactored into modern Fortran.

    module support_module

    use iso_fortran_env, only: wp => real64

    implicit none

    real(wp),parameter :: epmach = epsilon(1.0_wp)
    real(wp),parameter :: zero   = 0.0_wp
    real(wp),parameter :: one    = 1.0_wp
    real(wp),parameter :: two    = 2.0_wp
    real(wp),parameter :: four   = 4.0_wp
    real(wp),parameter :: ten    = 10.0_wp
    real(wp),parameter :: hun    = 100.0_wp
    real(wp),parameter :: alfmin = 0.1_wp

    contains
!*******************************************************************************

      subroutine daxpy(n,da,dx,incx,dy,incy)
      implicit none

      !! constant times a vector plus a vector.
      !! uses unrolled loops for increments equal to one.
      !! jack dongarra, linpack, 3/11/78.

      real(wp) :: dx(*) , dy(*) , da
      integer :: i , incx , incy , ix , iy , m , mp1 , n

      if ( n<=0 ) return
      if ( da==zero ) return
      if ( incx==1 .and. incy==1 ) then

        ! code for both increments equal to 1

        ! clean-up loop

         m = mod(n,4)
         if ( m/=0 ) then
            do i = 1 , m
               dy(i) = dy(i) + da*dx(i)
            enddo
            if ( n<4 ) return
         endif
         mp1 = m + 1
         do i = mp1 , n , 4
            dy(i) = dy(i) + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
         enddo
      else

         ! code for unequal increments or equal increments
         ! not equal to 1

         ix = 1
         iy = 1
         if ( incx<0 ) ix = (-n+1)*incx + 1
         if ( incy<0 ) iy = (-n+1)*incy + 1
         do i = 1 , n
            dy(iy) = dy(iy) + da*dx(ix)
            ix = ix + incx
            iy = iy + incy
         enddo
         return
      endif
      end subroutine daxpy

      subroutine dcopy(n,dx,incx,dy,incy)
      implicit none

!! copies a vector, x, to a vector, y.
!! uses unrolled loops for increments equal to one.
!! jack dongarra, linpack, 3/11/78.

      real(wp) :: dx(*) , dy(*)
      integer :: i , incx , incy , ix , iy , m , mp1 , n

      if ( n<=0 ) return
      if ( incx==1 .and. incy==1 ) then

!        code for both increments equal to 1

!        clean-up loop

         m = mod(n,7)
         if ( m/=0 ) then
            do i = 1 , m
               dy(i) = dx(i)
            enddo
            if ( n<7 ) return
         endif
         mp1 = m + 1
         do i = mp1 , n , 7
            dy(i) = dx(i)
            dy(i+1) = dx(i+1)
            dy(i+2) = dx(i+2)
            dy(i+3) = dx(i+3)
            dy(i+4) = dx(i+4)
            dy(i+5) = dx(i+5)
            dy(i+6) = dx(i+6)
         enddo
      else

!        code for unequal increments or equal increments
!        not equal to 1

         ix = 1
         iy = 1
         if ( incx<0 ) ix = (-n+1)*incx + 1
         if ( incy<0 ) iy = (-n+1)*incy + 1
         do i = 1 , n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
         enddo
         return
      endif
      end subroutine dcopy

      real(wp) function ddot(n,dx,incx,dy,incy)
      implicit none

!!  forms the dot product of two vectors.
!!  uses unrolled loops for increments equal to one.
!!  jack dongarra, linpack, 3/11/78.

      real(wp) :: dx(*) , dy(*) , dtemp
      integer :: i , incx , incy , ix , iy , m , mp1 , n

      ddot = zero
      dtemp = zero
      if ( n<=0 ) return
      if ( incx==1 .and. incy==1 ) then

!        code for both increments equal to 1

!        clean-up loop

         m = mod(n,5)
         if ( m/=0 ) then
            do i = 1 , m
               dtemp = dtemp + dx(i)*dy(i)
            enddo
            if ( n<5 ) then
               ddot = dtemp
               return
            endif
         endif
         mp1 = m + 1
         do i = mp1 , n , 5
            dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2) &
                    *dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
         enddo
         ddot = dtemp
      else

!        code for unequal increments or equal increments
!          not equal to 1

         ix = 1
         iy = 1
         if ( incx<0 ) ix = (-n+1)*incx + 1
         if ( incy<0 ) iy = (-n+1)*incy + 1
         do i = 1 , n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         enddo
         ddot = dtemp
         return
      endif
      end function ddot

      real(wp) function dnrm1(n,x,i,j)
      implicit none
      integer :: n , i , j , k
      real(wp) :: snormx , sum , x(n) , scale , temp

!!  dnrm1 - computes the i-norm of a vector
!!  between the ith and the jth elements

!      input -
!      n       length of vector
!      x       vector of length n
!      i       initial element of vector to be used
!      j       final element to use

!      output -
!      dnrm1   norm

      snormx = zero
      do k = i , j
         snormx = max(snormx,abs(x(k)))
      enddo
      dnrm1 = snormx
      if ( snormx==zero ) return
      scale = snormx
      if ( snormx>=one ) scale = sqrt(snormx)
      sum = zero
      do k = i , j
         temp = zero
         if ( abs(x(k))+scale/=scale ) temp = x(k)/snormx
         if ( one+temp/=one ) sum = sum + temp*temp
      enddo
      sum = sqrt(sum)
      dnrm1 = snormx*sum
      end function dnrm1

!>
!  replaced original routine with this one from blas:
!  http://netlib.sandia.gov/blas/dnrm2.f

        real(wp) function dnrm2(n,x,incx)
        implicit none

        integer,intent(in) :: incx
        integer,intent(in) :: n
        real(wp),intent(in) :: x(*)

  !  purpose
  !  =======
  !
  !  dnrm2 returns the euclidean norm of a vector via the function
  !  name, so that
  !
  !     dnrm2 := sqrt( x'*x )
  !
  !  further details
  !  ===============
  !
  !  -- this version written on 25-october-1982.
  !     modified on 14-october-1993 to inline the call to dlassq.
  !     sven hammarling, nag ltd.
  !
  !  =====================================================================

        real(wp) :: absxi , norm , scale , ssq
        integer :: ix

        if ( n<1 .or. incx<1 ) then
           norm = zero
        elseif ( n==1 ) then
           norm = abs(x(1))
        else
           scale = zero
           ssq = one
  !        the following loop is equivalent to this call to the lapack
  !        auxiliary routine:
  !        call dlassq( n, x, incx, scale, ssq )
           do ix = 1 , 1 + (n-1)*incx , incx
              if ( x(ix)/=zero ) then
                 absxi = abs(x(ix))
                 if ( scale<absxi ) then
                    ssq = one + ssq*(scale/absxi)**2
                    scale = absxi
                 else
                    ssq = ssq + (absxi/scale)**2
                 endif
              endif
           enddo
           norm = scale*sqrt(ssq)
        endif

        dnrm2 = norm

    end function dnrm2

      subroutine dsrot(n,dx,incx,dy,incy,c,s)
      implicit none

    !!  applies a plane rotation.
    !!  jack dongarra, linpack, 3/11/78.

      real(wp) :: dx(*) , dy(*) , dtemp , c , s
      integer :: i , incx , incy , ix , iy , n

      if ( n<=0 ) return
      if ( incx==1 .and. incy==1 ) then

         !code for both increments equal to 1

         do i = 1 , n
            dtemp = c*dx(i) + s*dy(i)
            dy(i) = c*dy(i) - s*dx(i)
            dx(i) = dtemp
         enddo
         return
      endif

      ! code for unequal increments or equal increments not equal to 1

      ix = 1
      iy = 1
      if ( incx<0 ) ix = (-n+1)*incx + 1
      if ( incy<0 ) iy = (-n+1)*incy + 1
      do i = 1 , n
         dtemp = c*dx(ix) + s*dy(iy)
         dy(iy) = c*dy(iy) - s*dx(ix)
         dx(ix) = dtemp
         ix = ix + incx
         iy = iy + incy
      enddo
      return
      end subroutine dsrot

      subroutine dsrotg(da,db,c,s)
      implicit none

    !!  construct givens plane rotation.
    !!  jack dongarra, linpack, 3/11/78.
    !!                 modified 9/27/86.

      real(wp) :: da , db , c , s , roe , scale , r , z

      roe = db
      if ( abs(da)>abs(db) ) roe = da
      scale = abs(da) + abs(db)
      if ( scale/=zero ) then
         r = scale*sqrt((da/scale)**2+(db/scale)**2)
         r = sign(one,roe)*r
         c = da/r
         s = db/r
      else
         c = one
         s = zero
         r = zero
      endif
      z = s
      if ( abs(c)>zero .and. abs(c)<=s ) z = one/c
      da = r
      db = z
      end subroutine dsrotg

      subroutine dscal(n,da,dx,incx)
      implicit none

    !!  scales a vector by a constant.
    !!  uses unrolled loops for increment equal to one.
    !!  jack dongarra, linpack, 3/11/78.

      real(wp) :: da , dx(*)
      integer :: i , incx , m , mp1 , n , nincx

      if ( n<=0 .or. incx<=0 ) return
      if ( incx==1 ) then

         ! code for increment equal to 1

         ! clean-up loop

         m = mod(n,5)
         if ( m/=0 ) then
            do i = 1 , m
               dx(i) = da*dx(i)
            enddo
            if ( n<5 ) return
         endif
         mp1 = m + 1
         do i = mp1 , n , 5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         enddo
      else

         ! code for increment not equal to 1

         nincx = n*incx
         do i = 1 , nincx , incx
            dx(i) = da*dx(i)
         enddo

      endif

      end subroutine dscal

!*******************************************************************************
    end module support_module
!*******************************************************************************

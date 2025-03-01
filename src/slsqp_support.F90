!*******************************************************************************
!> license: BSD
!
!  Support routines for SLSQP. For example, routines from
!  [BLAS](http://www.netlib.org/blas/) and [LINPACK](http://www.netlib.org/linpack/).
!  These have also been refactored into modern Fortran.

    module slsqp_support

    use slsqp_kinds

    implicit none

    private

    real(wp),parameter,public :: epmach = epsilon(1.0_wp)
    real(wp),parameter,public :: zero   = 0.0_wp
    real(wp),parameter,public :: one    = 1.0_wp
    real(wp),parameter,public :: two    = 2.0_wp
    real(wp),parameter,public :: four   = 4.0_wp
    real(wp),parameter,public :: ten    = 10.0_wp
    real(wp),parameter,public :: hun    = 100.0_wp

    public :: daxpy,dcopy,ddot,dnrm2,dscal

#ifdef HAS_BLAS

    ! linking with an external BLAS library.
    ! Define the interfaces here. Note that this
    ! will only work if the `wp` is the same kind
    ! as used in the BLAS library.

    interface
        pure subroutine daxpy(n,da,dx,incx,dy,incy)
            import :: wp
            implicit none
            integer,intent(in)                  :: n
            real(wp),intent(in)                 :: da
            real(wp),dimension(*),intent(in)    :: dx
            integer,intent(in)                  :: incx
            real(wp),dimension(*),intent(inout) :: dy
            integer,intent(in)                  :: incy
        end subroutine daxpy
        pure subroutine dcopy(n,dx,incx,dy,incy)
            import :: wp
            implicit none
            integer,intent(in)                :: n
            real(wp),dimension(*),intent(in)  :: dx
            integer,intent(in)                :: incx
            real(wp),dimension(*),intent(out) :: dy
            integer,intent(in)                :: incy
        end subroutine dcopy
        pure real(wp) function ddot(n,dx,incx,dy,incy)
            import :: wp
            implicit none
            integer,intent(in)               :: n
            real(wp),dimension(*),intent(in) :: dx
            integer,intent(in)               :: incx
            real(wp),dimension(*),intent(in) :: dy
            integer,intent(in)               :: incy
        end function ddot
        pure function dnrm2(n,x,incx) result(norm)
            import :: wp
            implicit none
            integer,intent(in)               :: incx
            integer,intent(in)               :: n
            real(wp),dimension(*),intent(in) :: x
            real(wp)                         :: norm
        end function dnrm2
        pure subroutine dscal(n,da,dx,incx)
            import :: wp
            implicit none
            integer,intent(in)                  :: n
            real(wp),intent(in)                 :: da
            real(wp),dimension(*),intent(inout) :: dx
            integer,intent(in)                  :: incx
        end subroutine dscal
    end interface

#else

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  constant times a vector plus a vector.
!  uses unrolled loops for increments equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

    pure subroutine daxpy(n,da,dx,incx,dy,incy)

    implicit none

    integer,intent(in)                  :: n
    real(wp),intent(in)                 :: da
    real(wp),dimension(*),intent(in)    :: dx
    integer,intent(in)                  :: incx
    real(wp),dimension(*),intent(inout) :: dy
    integer,intent(in)                  :: incy

    integer :: i , ix , iy , m , mp1

    if ( n<=0 ) return
    if ( abs(da)<=zero ) return
    if ( incx==1 .and. incy==1 ) then

        ! code for both increments equal to 1

        ! clean-up loop

        m = mod(n,4)
        if ( m/=0 ) then
            do i = 1 , m
                dy(i) = dy(i) + da*dx(i)
            end do
            if ( n<4 ) return
        end if
        mp1 = m + 1
        do i = mp1 , n , 4
            dy(i) = dy(i) + da*dx(i)
            dy(i+1) = dy(i+1) + da*dx(i+1)
            dy(i+2) = dy(i+2) + da*dx(i+2)
            dy(i+3) = dy(i+3) + da*dx(i+3)
        end do

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
        end do

    end if

    end subroutine daxpy
!*******************************************************************************

!*******************************************************************************
!>
!  copies a vector, x, to a vector, y.
!  uses unrolled loops for increments equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

    pure subroutine dcopy(n,dx,incx,dy,incy)

    implicit none

    integer,intent(in)                :: n
    real(wp),dimension(*),intent(in)  :: dx
    integer,intent(in)                :: incx
    real(wp),dimension(*),intent(out) :: dy
    integer,intent(in)                :: incy

    integer :: i , ix , iy , m , mp1

    if ( n<=0 ) return
    if ( incx==1 .and. incy==1 ) then

        ! code for both increments equal to 1

        ! clean-up loop

        m = mod(n,7)
        if ( m/=0 ) then
            do i = 1 , m
                dy(i) = dx(i)
            end do
            if ( n<7 ) return
        end if
        mp1 = m + 1
        do i = mp1 , n , 7
            dy(i) = dx(i)
            dy(i+1) = dx(i+1)
            dy(i+2) = dx(i+2)
            dy(i+3) = dx(i+3)
            dy(i+4) = dx(i+4)
            dy(i+5) = dx(i+5)
            dy(i+6) = dx(i+6)
        end do

    else

        ! code for unequal increments or equal increments
        ! not equal to 1

        ix = 1
        iy = 1
        if ( incx<0 ) ix = (-n+1)*incx + 1
        if ( incy<0 ) iy = (-n+1)*incy + 1
        do i = 1 , n
            dy(iy) = dx(ix)
            ix = ix + incx
            iy = iy + incy
        end do

    end if

    end subroutine dcopy
!*******************************************************************************

!*******************************************************************************
!>
!  forms the dot product of two vectors.
!  uses unrolled loops for increments equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

    pure real(wp) function ddot(n,dx,incx,dy,incy)

    implicit none

    integer,intent(in)               :: n
    real(wp),dimension(*),intent(in) :: dx
    integer,intent(in)               :: incx
    real(wp),dimension(*),intent(in) :: dy
    integer,intent(in)               :: incy

    real(wp) :: dtemp
    integer :: i , ix , iy , m , mp1

    ddot = zero
    dtemp = zero
    if ( n<=0 ) return
    if ( incx==1 .and. incy==1 ) then

        ! code for both increments equal to 1

        ! clean-up loop

        m = mod(n,5)
        if ( m/=0 ) then
            do i = 1 , m
                dtemp = dtemp + dx(i)*dy(i)
            end do
            if ( n<5 ) then
                ddot = dtemp
                return
            end if
        end if
        mp1 = m + 1
        do i = mp1 , n , 5
            dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + &
                    dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
        end do
        ddot = dtemp

    else

        ! code for unequal increments or equal increments
        ! not equal to 1

        ix = 1
        iy = 1
        if ( incx<0 ) ix = (-n+1)*incx + 1
        if ( incy<0 ) iy = (-n+1)*incy + 1
        do i = 1 , n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
        end do
        ddot = dtemp

    end if

    end function ddot
!*******************************************************************************

!*******************************************************************************
!>
!  Function that returns the Euclidean norm
!  \( \sqrt{ \mathbf{x}^T \mathbf{x} } \) of a vector \( \mathbf{x} \).
!
!### Further details
!
!  * this version written on 25-october-1982.
!  * modified on 14-october-1993 to inline the call to dlassq.
!    sven hammarling, nag ltd.
!  * Converted to modern Fortran, Jacob Williams, Jan. 2016.
!
!@note Replaced original SLSQP routine with this one from
!      [BLAS](http://netlib.sandia.gov/blas/dnrm2.f).

    pure function dnrm2(n,x,incx) result(norm)

    implicit none

    integer,intent(in)               :: incx
    integer,intent(in)               :: n
    real(wp),dimension(*),intent(in) :: x
    real(wp)                         :: norm

    real(wp) :: absxi , scale , ssq
    integer :: ix

    if ( n<1 .or. incx<1 ) then
        norm = zero
    elseif ( n==1 ) then
        norm = abs(x(1))
    else
        scale = zero
        ssq = one
        ! the following loop is equivalent to this call to the lapack
        ! auxiliary routine:
        ! call dlassq( n, x, incx, scale, ssq )
        do ix = 1 , 1 + (n-1)*incx , incx
            if ( abs(x(ix))>zero ) then
                absxi = abs(x(ix))
                if ( scale<absxi ) then
                    ssq = one + ssq*(scale/absxi)**2
                    scale = absxi
                else
                    ssq = ssq + (absxi/scale)**2
                end if
            end if
        end do
        norm = scale*sqrt(ssq)
    end if

    end function dnrm2
!*******************************************************************************

!*******************************************************************************
!>
!  scales a vector by a constant.
!  uses unrolled loops for increment equal to one.
!
!### Author
!  jack dongarra, linpack, 3/11/78.

    pure subroutine dscal(n,da,dx,incx)

    implicit none

    integer,intent(in)                  :: n
    real(wp),intent(in)                 :: da
    real(wp),dimension(*),intent(inout) :: dx
    integer,intent(in)                  :: incx

    integer :: i , m , mp1 , nincx

    if ( n<=0 .or. incx<=0 ) return
    if ( incx==1 ) then

        ! code for increment equal to 1

        ! clean-up loop

        m = mod(n,5)
        if ( m/=0 ) then
            do i = 1 , m
                dx(i) = da*dx(i)
            end do
            if ( n<5 ) return
        end if
        mp1 = m + 1
        do i = mp1 , n , 5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
        end do

    else

        ! code for increment not equal to 1

        nincx = n*incx
        do i = 1 , nincx , incx
            dx(i) = da*dx(i)
        end do

    end if

    end subroutine dscal
!*******************************************************************************

#endif

!*******************************************************************************
    end module slsqp_support
!*******************************************************************************

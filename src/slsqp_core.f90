!*******************************************************************************
!> license: BSD
!
!  Core subroutines for the SLSQP optimization method.
!  These are refactoried versions of the original routines.

    module slsqp_core

    use slsqp_kinds
    use slsqp_support
    use bvls_module,     only: bvls_wrapper
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan

    implicit none

    private

    real(wp),parameter :: eps = epsilon(1.0_wp)  !! machine precision

    type,public :: linmin_data
        !! data formerly saved in [[linmin]] routine.
        real(wp) :: a    = zero
        real(wp) :: b    = zero
        real(wp) :: d    = zero
        real(wp) :: e    = zero
        real(wp) :: p    = zero
        real(wp) :: q    = zero
        real(wp) :: r    = zero
        real(wp) :: u    = zero
        real(wp) :: v    = zero
        real(wp) :: w    = zero
        real(wp) :: x    = zero
        real(wp) :: m    = zero
        real(wp) :: fu   = zero
        real(wp) :: fv   = zero
        real(wp) :: fw   = zero
        real(wp) :: fx   = zero
        real(wp) :: tol1 = zero
        real(wp) :: tol2 = zero
    contains
        procedure :: destroy => destroy_linmin_data
    end type linmin_data

    type,public :: slsqpb_data
        !! data formerly saved in [[slsqpb]].
        real(wp) :: t       = zero
        real(wp) :: f0      = zero
        real(wp) :: h1      = zero
        real(wp) :: h2      = zero
        real(wp) :: h3      = zero
        real(wp) :: h4      = zero
        real(wp) :: t0      = zero
        real(wp) :: gs      = zero
        real(wp) :: tol     = zero
        real(wp) :: alpha   = zero
        integer  :: line    = 0
        integer  :: iexact  = 0
        integer  :: incons  = 0
        integer  :: ireset  = 0
        integer  :: itermx  = 0
        integer  :: n1      = 0
        integer  :: n2      = 0
        integer  :: n3      = 0
    contains
        procedure :: destroy => destroy_slsqpb_data
    end type slsqpb_data

    public :: slsqp

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  **slsqp**: **s**equential **l**east **sq**uares **p**rogramming
!  to solve general nonlinear optimization problems
!
!  a nonlinear programming method with quadratic programming subproblems
!  this subroutine solves the general nonlinear programming problem:
!
!  **minimize**
!
!  * \( f(x) \)
!
!  **subject to**
!
!  * \( c_j (x) = 0 \),           \( j = 1,...,meq   \)
!  * \( c_j (x) \ge 0 \),         \( j = meq+1,...,m \)
!  * \( xl_i \le x_i \le xu_i \), \( i = 1,...,n     \)
!
!  the algorithm implements the method of Han and Powell
!  with BFGS-update of the b-matrix and L1-test function
!  within the steplength algorithm.
!
!### Reference
!   * Dieter Kraft: "A software package for sequential quadratic programming",
!     DFVLR-FB 88-28, 1988
!
!### History
!   * implemented by: Dieter Kraft, DFVLR oberpfaffenhofen
!   * date: april - october, 1981.
!   * December, 31-st, 1984.
!   * March   , 21-st, 1987, revised to fortran 77
!   * March   , 20-th, 1989, revised to ms-fortran
!   * April   , 14-th, 1989, hesse   in-line coded
!   * February, 28-th, 1991, fortran/2 version 1.04 accepts statement functions
!   * March   ,  1-st, 1991, tested with salford ftn77/386 compiler vers 2.40 in protected mode
!   * January ,        2016, Refactored into modern Fortran by Jacob Williams
!
!### License
!  Original version copyright 1991: Dieter Kraft, FHM.
!  Released under a BSD license.
!
!@note `f`, `c`, `g`, `a` must all be set by the user before each call.

    subroutine slsqp(m,meq,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,w,l_w, &
                     sdat,ldat,alphamin,alphamax,tolf,toldf,toldx,&
                     max_iter_ls,nnls_mode,infinite_bound)

    implicit none

    integer,intent(in) :: m                     !! is the total number of constraints, \( m \ge 0 \)
    integer,intent(in) :: meq                   !! is the number of equality constraints, \( m_{eq} \ge 0 \)
    integer,intent(in) :: la                    !! see `a`, \( la \ge \max(m,1) \)
    integer,intent(in) :: n                     !! is the number of variables, \( n \ge 1 \)
    real(wp),dimension(n),intent(inout) :: x    !! `x()` stores the current iterate of the `n` vector `x`
                                                !! on entry `x()` must be initialized. on exit `x()`
                                                !! stores the solution vector `x` if `mode = 0`.
    real(wp),dimension(n),intent(in) :: xl      !! `xl()` stores an n vector of lower bounds `xl` to `x`.
    real(wp),dimension(n),intent(in) :: xu      !! `xu()` stores an n vector of upper bounds `xu` to `x`.
    real(wp),intent(in) :: f                    !! is the value of the objective function.
    real(wp),dimension(la),intent(in) :: c      !! `c()` stores the `m` vector `c` of constraints,
                                                !! equality constraints (if any) first.
                                                !! dimension of `c` must be greater or equal `la`,
                                                !! which must be greater or equal `max(1,m)`.
    real(wp),dimension(n+1),intent(in) :: g     !! `g()` stores the `n` vector `g` of partials of the
                                                !! objective function; dimension of `g` must be
                                                !! greater or equal `n+1`.
    real(wp),dimension(la,n+1),intent(in) ::  a !! the `la` by `n + 1` array `a()` stores
                                                !! the `m` by `n` matrix `a` of constraint normals.
                                                !! `a()` has first dimensioning parameter `la`,
                                                !! which must be greater or equal `max(1,m)`.
    real(wp),intent(inout) :: acc   !! `abs(acc)` controls the final accuracy.
                                    !! if `acc` < zero an exact linesearch is performed,
                                    !! otherwise an armijo-type linesearch is used.
    integer,intent(inout) :: iter   !! prescribes the maximum number of iterations.
                                    !! on exit `iter` indicates the number of iterations.
    integer,intent(inout) :: mode   !! mode controls calculation:
                                    !!
                                    !! reverse communication is used in the sense that
                                    !! the program is initialized by `mode = 0`; then it is
                                    !! to be called repeatedly by the user until a return
                                    !! with `mode /= abs(1)` takes place.
                                    !! if `mode = -1` gradients have to be calculated,
                                    !! while with `mode = 1` functions have to be calculated.
                                    !! mode must not be changed between subsequent calls of [[slsqp]].
                                    !!
                                    !! **evaluation modes**:
                                    !!
                                    !! * ** -1 **: gradient evaluation, (`g` & `a`)
                                    !! * **  0 **: *on entry*: initialization, (`f`, `g`, `c`, `a`),
                                    !!   *on exit*: required accuracy for solution obtained
                                    !! * **  1 **: function evaluation, (`f` & `c`)
                                    !!
                                    !! **failure modes**:
                                    !!
                                    !! * ** 2 **: number of equality constraints larger than `n`
                                    !! * ** 3 **: more than `3*n` iterations in [[lsq]] subproblem
                                    !! * ** 4 **: inequality constraints incompatible
                                    !! * ** 5 **: singular matrix `e` in [[lsq]] subproblem
                                    !! * ** 6 **: singular matrix `c` in [[lsq]] subproblem
                                    !! * ** 7 **: rank-deficient equality constraint subproblem [[hfti]]
                                    !! * ** 8 **: positive directional derivative for linesearch
                                    !! * ** 9 **: more than `iter` iterations in sqp
                                    !! * ** >=10 **: working space `w` too small,
                                    !!   `w` should be enlarged to `l_w=mode/1000`,
    integer,intent(in) :: l_w       !! the length of `w`, which should be at least:
                                    !!
                                    !! * `(3*n1+m)*(n1+1)`                     **for lsq**
                                    !! * `+(n1-meq+1)*(mineq+2) + 2*mineq`     **for lsi**
                                    !! * `+(n1+mineq)*(n1-meq) + 2*meq + n1`   **for lsei**
                                    !! * `+ n1*n/2 + 2*m + 3*n + 3*n1 + 1`     **for slsqpb**
                                    !!
                                    !! with `mineq = m - meq + 2*n1` & `n1 = n+1`
    real(wp),dimension(l_w),intent(inout) :: w  !! `w()` is a one dimensional working space.
                                                !! the first `m+n+n*n1/2` elements of `w` must not be
                                                !! changed between subsequent calls of [[slsqp]].
                                                !! on return `w(1) ... w(m)` contain the multipliers
                                                !! associated with the general constraints, while
                                                !! `w(m+1) ... w(m+n(n+1)/2)` store the cholesky factor
                                                !! `l*d*l(t)` of the approximate hessian of the
                                                !! lagrangian columnwise dense as lower triangular
                                                !! unit matrix `l` with `d` in its 'diagonal' and
                                                !! `w(m+n(n+1)/2+n+2 ... w(m+n(n+1)/2+n+2+m+2n)`
                                                !! contain the multipliers associated with all
                                                !! constraints of the quadratic program finding
                                                !! the search direction to the solution `x*`
    type(slsqpb_data),intent(inout) :: sdat  !! data for [[slsqpb]].
    type(linmin_data),intent(inout) :: ldat  !! data for [[linmin]].
    real(wp),intent(in) :: alphamin  !! min \( \alpha \) for line search
                                     !! \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)
    real(wp),intent(in) :: alphamax  !! max \( \alpha \) for line search
                                     !! \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)
    real(wp),intent(in) :: tolf      !! stopping criterion if \( |f| < tolf \) then stop.
    real(wp),intent(in) :: toldf     !! stopping criterion if \( |f_{n+1} - f_n| < toldf \) then stop.
    real(wp),intent(in) :: toldx     !! stopping criterion if \( ||x_{n+1} - x_n|| < toldx \) then stop.
    integer,intent(in)  :: max_iter_ls !! maximum number of iterations in the [[nnls]] problem
    integer,intent(in)  :: nnls_mode !! which NNLS method to use
    real(wp),intent(in) :: infinite_bound !! "infinity" for the upper and lower bounds.
                                          !! if `xl<=-infinite_bound` or `xu>=infinite_bound`
                                          !! then these bounds are considered nonexistant.
                                          !! If `infinite_bound=0` then `huge()` is used for this.

    integer :: il , im , ir , is , iu , iv , iw , ix , mineq, n1
    real(wp) :: infBnd !! local copy of `infinite_bound`

    if (infinite_bound==zero) then
        infBnd = huge(one)
    else
        infBnd = abs(infinite_bound)
    end if

    ! check length of working arrays
    n1 = n + 1
    mineq = m - meq + n1 + n1
    il = (3*n1+m)*(n1+1) + (n1-meq+1)*(mineq+2) + 2*mineq + (n1+mineq)&
         *(n1-meq) + 2*meq + n1*n/2 + 2*m + 3*n + 4*n1 + 1
    im = max(mineq,n1-meq)
    if ( l_w<il ) then
        mode = 1000*max(10,il)
        mode = mode + max(10,im)
        iter = 0
        return
    end if
    if (meq>n) then
        ! note: calling lsq when meq>n is corrupting the
        ! memory in some way, so just catch this here.
        mode = 2
        iter = 0
        return
    end if

    ! prepare data for calling sqpbdy  -  initial addresses in w

    im = 1
    il = im + max(1,m)
    il = im + la
    ix = il + n1*n/2 + 1
    ir = ix + n
    is = ir + n + n + max(1,m)
    is = ir + n + n + la
    iu = is + n1
    iv = iu + n1
    iw = iv + n1

    sdat%n1 = n1

    call slsqpb(m,meq,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,&
                w(ir),w(il),w(ix),w(im),w(is),w(iu),w(iv),w(iw),&
                sdat%t,sdat%f0,sdat%h1,sdat%h2,sdat%h3,sdat%h4,&
                sdat%n1,sdat%n2,sdat%n3,sdat%t0,sdat%gs,sdat%tol,sdat%line,&
                sdat%alpha,sdat%iexact,sdat%incons,sdat%ireset,sdat%itermx,&
                ldat,alphamin,alphamax,tolf,toldf,toldx,&
                max_iter_ls,nnls_mode,infBnd)

    end subroutine slsqp
!*******************************************************************************

!*******************************************************************************
!>
!  nonlinear programming by solving sequentially quadratic programs
!
!  l1 - line search, positive definite bfgs update

    subroutine slsqpb(m,meq,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,&
                      r,l,x0,mu,s,u,v,w,&
                      t,f0,h1,h2,h3,h4,n1,n2,n3,t0,gs,tol,line,&
                      alpha,iexact,incons,ireset,itermx,ldat,&
                      alphamin,alphamax,tolf,toldf,toldx,&
                      max_iter_ls,nnls_mode,infBnd)
    implicit none

    integer,intent(in)                  :: m
    integer,intent(in)                  :: meq
    integer,intent(in)                  :: la
    integer,intent(in)                  :: n
    real(wp),dimension(n)               :: x
    real(wp),dimension(n)               :: xl
    real(wp),dimension(n)               :: xu
    real(wp)                            :: f
    real(wp),dimension(la)              :: c
    real(wp),dimension(n+1)             :: g
    real(wp),dimension(la,n+1)          :: a
    real(wp)                            :: acc
    integer,intent(inout)               :: iter     !! **in:**  maximum number of iterations.
                                                    !! **out:** actual number of iterations.
    integer,intent(inout)               :: mode
    real(wp),dimension(m+n+n+2)         :: r
    real(wp),dimension((n+1)*(n+2)/2)   :: l
    real(wp),dimension(n)               :: x0
    real(wp),dimension(la)              :: mu
    real(wp),dimension(n+1)             :: s
    real(wp),dimension(n+1)             :: u
    real(wp),dimension(n+1)             :: v
    real(wp),dimension(*),intent(inout) :: w   !! `dim(w)` =
                                               !!
                                               !! * `n1*(n1+1) + meq*(n1+1) + mineq*(n1+1)`   for [[lsq]]
                                               !! * `+(n1-meq+1)*(mineq+2) + 2*mineq`         for [[lsi]]
                                               !! * `+(n1+mineq)*(n1-meq) + 2*meq + n1`       for [[lsei]]
                                               !!
                                               !! with `mineq = m - meq + 2*n1` & `n1 = n+1`
    real(wp),intent(inout)              :: t
    real(wp),intent(inout)              :: f0
    real(wp),intent(inout)              :: h1
    real(wp),intent(inout)              :: h2
    real(wp),intent(inout)              :: h3
    real(wp),intent(inout)              :: h4
    integer,intent(inout)               :: n1
    integer,intent(inout)               :: n2
    integer,intent(inout)               :: n3
    real(wp),intent(inout)              :: t0
    real(wp),intent(inout)              :: gs
    real(wp),intent(inout)              :: tol
    integer,intent(inout)               :: line
    real(wp),intent(inout)              :: alpha
    integer,intent(inout)               :: iexact
    integer,intent(inout)               :: incons
    integer,intent(inout)               :: ireset
    integer,intent(inout)               :: itermx
    type(linmin_data),intent(inout)     :: ldat      !! data for [[linmin]].
    real(wp),intent(in)                 :: alphamin  !! min \( \alpha \) for line search
                                                     !! \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)
    real(wp),intent(in)                 :: alphamax  !! max \( \alpha \) for line search
                                                     !! \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)
    real(wp),intent(in)                 :: tolf      !! stopping criterion if \( |f| < tolf \) then stop.
    real(wp),intent(in)                 :: toldf     !! stopping criterion if \( |f_{n+1} - f_n| < toldf \) then stop
    real(wp),intent(in)                 :: toldx     !! stopping criterion if \( ||x_{n+1} - x_n|| < toldx \) then stop
    integer,intent(in)                  :: max_iter_ls !! maximum number of iterations in the [[nnls]] problem
    integer,intent(in)                  :: nnls_mode !! which NNLS method to use
    real(wp),intent(in)                 :: infBnd !! "infinity" for the upper and lower bounds.

    integer :: i, j, k
    logical :: inconsistent_linearization !! if the SQP problem is inconsistent

    inconsistent_linearization = .false. ! initialize

    if ( mode<0 ) then

        ! call jacobian at current x

        ! update cholesky-factors of hessian matrix by modified bfgs formula

        do i = 1 , n
            u(i) = g(i) - ddot(m,a(1,i),1,r,1) - v(i)
        end do

        ! l'*s

        k = 0
        do i = 1 , n
            h1 = zero
            k = k + 1
            do j = i + 1 , n
                k = k + 1
                h1 = h1 + l(k)*s(j)
            end do
            v(i) = s(i) + h1
        end do

        ! d*l'*s

        k = 1
        do i = 1 , n
            v(i) = l(k)*v(i)
            k = k + n1 - i
        end do

        ! l*d*l'*s

        do i = n , 1 , -1
            h1 = zero
            k = i
            do j = 1 , i - 1
                h1 = h1 + l(k)*v(j)
                k = k + n - j
            end do
            v(i) = v(i) + h1
        end do

        h1 = ddot(n,s,1,u,1)
        h2 = ddot(n,s,1,v,1)
        h3 = 0.2_wp*h2
        if ( h1<h3 ) then
            h4 = (h2-h3)/(h2-h1)
            h1 = h3
            call dscal(n,h4,u,1)
            call daxpy(n,one-h4,v,1,u,1)
        end if
        if (h1==zero .or. h2==zero) then
            ! Singular update: reset hessian.
            ! [ JW : this is based on a SciPy update ]
            call reset_bfgs_matrix()
            if ( ireset>5 ) return
        else
            call ldl(n,l,u,+one/h1,v)
            call ldl(n,l,v,-one/h2,u)
        end if
        ! end of main iteration

    else if ( mode==0 ) then

        itermx = iter
        if ( acc>=zero ) then
            iexact = 0
        else
            iexact = 1
        end if
        acc = abs(acc)
        tol = ten*acc
        iter = 0
        ireset = 0
        n1 = n + 1
        n2 = n1*n/2
        n3 = n2 + 1
        s(1) = zero
        mu(1) = zero
        call dcopy(n,s(1),0,s,1)
        call dcopy(m,mu(1),0,mu,1)

        call reset_bfgs_matrix()
        if ( ireset>5 ) return

    else

        ! call functions at current x

        t = f
        do j = 1 , m
            if ( j<=meq ) then
                h1 = c(j)
            else
                h1 = zero
            end if
            t = t + mu(j)*max(-c(j),h1)
        end do
        h1 = t - t0
        select case (iexact)
        case(0)
            if ( h1<=h3/ten .or. line>10 ) then
                call convergence_check(acc,0,-1)
            else
                alpha = min(max(h3/(two*(h3-h1)),alphamin),alphamax)
                call inexact_linesearch()
            end if
        case(1)
            call exact_linesearch()
            if ( line==3 ) call convergence_check(acc,0,-1)
        end select

        return

    end if

    do
        ! main iteration : search direction, steplength, ldl'-update

        iter = iter + 1
        mode = 9
        if ( iter>itermx ) return

        ! search direction as solution of qp - subproblem

        call dcopy(n,xl,1,u,1)
        call dcopy(n,xu,1,v,1)
        call daxpy(n,-one,x,1,u,1)
        call daxpy(n,-one,x,1,v,1)
        h4 = one
        call lsq(m,meq,n,n3,la,l,g,a,c,u,v,s,r,w,mode,max_iter_ls,nnls_mode,infBnd)

        ! augmented problem for inconsistent linearization

        inconsistent_linearization = .false. ! initialize
        if ( mode==6 ) then
            if ( n==meq ) mode = 4
        end if
        if ( mode==4 ) then
            ! Will reject this iteration if the SQP problem is inconsistent,
            inconsistent_linearization = .true.
            do j = 1 , m
                if ( j<=meq ) then
                    a(j,n1) = -c(j)
                else
                    a(j,n1) = max(-c(j),zero)
                end if
            end do
            s(1) = zero
            call dcopy(n,s(1),0,s,1)
            h3 = zero
            g(n1) = zero
            l(n3) = hun
            s(n1) = one
            u(n1) = zero
            v(n1) = one
            incons = 0
            do
                call lsq(m,meq,n1,n3,la,l,g,a,c,u,v,s,r,w,mode,max_iter_ls,nnls_mode,infBnd)
                h4 = one - s(n1)
                if ( mode==4 ) then
                    l(n3) = ten*l(n3)
                    incons = incons + 1
                    if ( incons<=5 ) cycle
                    return
                else if ( mode/=1 ) then
                    return
                else
                    exit
                end if
            end do

        else if ( mode/=1 ) then
            return
        end if

        ! update multipliers for l1-test

        do i = 1 , n
            v(i) = g(i) - ddot(m,a(1,i),1,r,1)
        end do
        f0 = f
        call dcopy(n,x,1,x0,1)
        gs = ddot(n,g,1,s,1)
        h1 = abs(gs)
        h2 = zero
        do j = 1 , m
            if ( j<=meq ) then
                h3 = c(j)
            else
                h3 = zero
            end if
            h2 = h2 + max(-c(j),h3)
            h3 = abs(r(j))
            mu(j) = max(h3,(mu(j)+h3)/two)
            h1 = h1 + h3*abs(c(j))
        end do

        ! check convergence

        mode = 0
        if ( h1<acc .and. h2<acc .and. &
             .not. inconsistent_linearization .and. &
             .not. ieee_is_nan(f)) return
        h1 = zero
        do j = 1 , m
            if ( j<=meq ) then
                h3 = c(j)
            else
                h3 = zero
            end if
            h1 = h1 + mu(j)*max(-c(j),h3)
        end do
        t0 = f + h1
        h3 = gs - h1*h4
        mode = 8
        if ( h3>=zero ) then
            call reset_bfgs_matrix()
            if ( ireset>5 ) return
        else
            exit
        end if

    end do

    ! line search with an l1-testfunction

    line = 0
    alpha = alphamax
    if ( iexact==1 ) then
        call exact_linesearch()
        if ( line==3 ) call convergence_check(acc,0,-1)
    else
        call inexact_linesearch()
    end if

    contains

        subroutine reset_bfgs_matrix()  ! 100
            !! reset BFGS matrix
            ireset = ireset + 1
            if ( ireset>5 ) then
                ! check relaxed convergence in case of positive directional derivative
                ! [ JW: reuse this routine so that h3 is recomputed.
                !   this is based on a SciPy update to SLSQP ]
                call convergence_check(tol,0,8)
                ! the caller should return in this case
            else
                l(1) = zero
                call dcopy(n2,l(1),0,l,1)
                j = 1
                do i = 1 , n
                    l(j) = one
                    j = j + n1 - i
                end do
            end if
        end subroutine reset_bfgs_matrix

        subroutine inexact_linesearch()  ! 300
            line = line + 1
            h3 = alpha*h3
            call dscal(n,alpha,s,1)
            call dcopy(n,x0,1,x,1)
            call daxpy(n,one,s,1,x,1)
            call enforce_bounds(x,xl,xu,infBnd)  ! ensure that x doesn't violate bounds
            mode = 1
        end subroutine inexact_linesearch

        subroutine exact_linesearch() ! 400
            if ( line/=3 ) then
                alpha = linmin(line,alphamin,alphamax,t,tol, &
                               ldat%a, ldat%b, ldat%d, ldat%e, ldat%p,   ldat%q,   &
                               ldat%r, ldat%u, ldat%v, ldat%w, ldat%x,   ldat%m,   &
                               ldat%fu,ldat%fv,ldat%fw,ldat%fx,ldat%tol1,ldat%tol2 )
                call dcopy(n,x0,1,x,1)
                call daxpy(n,alpha,s,1,x,1)
                mode = 1
            else
                call dscal(n,alpha,s,1)
            end if
        end subroutine exact_linesearch

        subroutine convergence_check(tolerance,converged,not_converged)  ! 500
            real(wp),intent(in) :: tolerance !! tolerance
            integer,intent(in) :: converged     !! mode value if converged
            integer,intent(in) :: not_converged !! mode value if not converged

            h3 = zero
            do j = 1 , m
                if (j<=meq) then
                    h1 = c(j)
                else
                    h1 = zero
                end if
                h3 = h3 + max(-c(j),h1)
            end do
            mode = check_convergence(n,f,f0,x,x0,s,h3,tolerance,tolf,toldf,toldx,&
                                     converged,not_converged,inconsistent_linearization)
        end subroutine convergence_check

    end subroutine slsqpb
!*******************************************************************************

!*******************************************************************************
!>
!  Check for convergence.

    pure function check_convergence(n,f,f0,x,x0,s,h3,acc,tolf,toldf,toldx,&
                                    converged,not_converged,inconsistent_linearization) result(mode)

    implicit none

    integer,intent(in) :: n
    real(wp),intent(in) :: f
    real(wp),intent(in) :: f0
    real(wp),dimension(:),intent(in) :: x
    real(wp),dimension(:),intent(in) :: x0
    real(wp),dimension(:),intent(in) :: s
    real(wp),intent(in) :: h3
    real(wp),intent(in) :: acc
    real(wp),intent(in) :: tolf
    real(wp),intent(in) :: toldf
    real(wp),intent(in) :: toldx
    integer,intent(in) :: converged     !! mode value if converged
    integer,intent(in) :: not_converged !! mode value if not converged
    logical,intent(in) :: inconsistent_linearization !! if the SQP problem is inconsistent (will return `not_converged`)
    integer :: mode

    logical :: ok ! temp variable
    real(wp),dimension(n) :: xmx0

    if (h3>=acc .or. inconsistent_linearization .or. ieee_is_nan(f)) then
        mode = not_converged
    else

        ! if any are OK then it is converged
        ok = .false.
        if (.not. ok) ok = abs(f-f0)<acc
        if (.not. ok) ok = dnrm2(n,s,1)<acc
        ! note that these can be ignored if they are < 0:
        if (.not. ok .and. tolf>=zero)  ok = abs(f)<tolf
        if (.not. ok .and. toldf>=zero) ok = abs(f-f0)<toldf
        if (.not. ok .and. toldx>=zero) then
            xmx0 = x-x0 ! to avoid array temporary warning
            ok = dnrm2(n,xmx0,1)<toldx
        end if

        if (ok) then
            mode = converged
        else
            mode = not_converged
        end if

    end if

    end function check_convergence
!*******************************************************************************

!*******************************************************************************
!>
!  Minimize \( || e x - f || \) with respect to \(x\),
!  with upper triangular matrix \( e = + d ^{1/2} l^T \),
!  and vector \( f = -d^{-1/2} l^{-1} g \),
!  where the unit lower tridiangular matrix \(l\) is stored columnwise
!  dense in the \(n*(n+1)/2\) array \(l\) with vector \(d\) stored in its
!  'diagonal' thus substituting the one-elements of \(l\)
!
!  subject to:
!
!  * \( a(j)*x - b(j) = 0,              j=1,...,meq  \),
!  * \( a(j)*x - b(j) \ge 0,            j=meq+1,...,m\),
!  * \( x_l(i) \le x(i) \le x_u(i),     i=1,...,n    \),
!
!  On entry, the user has to provide the arrays `l`, `g`, `a`, `b`, `xl`, `xu`.
!  with dimensions: `l(n*(n+1)/2)`, `g(n)`, `a(la,n)`, `b(m)`, `xl(n)`, `xu(n)`.
!
!  The working array `w` must have at least the following dimension: dim(w) =
!
!  * `(3*n+m)*(n+1)`                    for [[lsq]]
!  * `+(n-meq+1)*(mineq+2) + 2*mineq`   for [[lsi]]
!  * `+(n+mineq)*(n-meq) + 2*meq + n`   for [[lsei]]
!
!  with `mineq = m - meq + 2*n`
!
!  On return, no array will be changed by the subroutine.
!
!### History
!  * coded dieter kraft, april 1987
!  * revised march 1989

    subroutine lsq(m,meq,n,nl,la,l,g,a,b,xl,xu,x,y,w,mode,max_iter_ls,nnls_mode,infBnd)

    implicit none

    integer,intent(in)        :: m
    integer,intent(in)        :: n
    integer,intent(in)        :: meq
    integer,intent(in)        :: nl
    integer,intent(in)        :: la
    real(wp),dimension(n)     :: x  !! stores the n-dimensional solution vector
    real(wp),dimension(m+n+n) :: y  !! stores the vector of lagrange multipliers of dimension
                                    !! m+n+n (constraints+lower+upper bounds)
    integer :: mode                 !! is a success-failure flag with the following meanings:
                                    !!
                                    !! * **1:** successful computation,
                                    !! * **2:** error return because of wrong dimensions (`n<1`),
                                    !! * **3:** iteration count exceeded by [[nnls]],
                                    !! * **4:** inequality constraints incompatible,
                                    !! * **5:** matrix `e` is not of full rank,
                                    !! * **6:** matrix `c` is not of full rank,
                                    !! * **7:** rank defect in [[hfti]]
    integer,intent(in) :: max_iter_ls !! maximum number of iterations in the [[nnls]] problem
    integer,intent(in) :: nnls_mode !! which NNLS method to use
    real(wp),intent(in) :: infbnd !! "infinity" for the upper and lower bounds.

    real(wp),dimension(nl)   :: l
    real(wp),dimension(n)    :: g
    real(wp),dimension(la,n) :: a
    real(wp),dimension(la)   :: b
    real(wp),dimension(*)    :: w
    real(wp),dimension(n)    :: xl
    real(wp),dimension(n)    :: xu
    real(wp) :: diag , xnorm
    integer :: i , ic , id , ie , if , ig , ih , il , im , ip , &
               iu , iw , i1 , i2 , i3 , i4 , mineq , &
               m1 , n1 , n2 , n3, num_unbounded, j

      n1 = n + 1
      mineq = m - meq
      m1 = mineq + n + n

      !  determine whether to solve problem
      !  with inconsistent linerarization (n2=1)
      !  or not (n2=0)

      n2 = n1*n/2 + 1
      if ( n2==nl ) then
         n2 = 0
      else
         n2 = 1
      end if
      n3 = n - n2

      !  recover matrix e and vector f from l and g

      i2 = 1
      i3 = 1
      i4 = 1
      ie = 1
      if = n*n + 1
      do i = 1 , n3
         i1 = n1 - i
         diag = sqrt(l(i2))
         w(i3) = zero
         call dcopy(i1,w(i3),0,w(i3),1)
         call dcopy(i1-n2,l(i2),1,w(i3),n)
         call dscal(i1-n2,diag,w(i3),n)
         w(i3) = diag
         w(if-1+i) = (g(i)-ddot(i-1,w(i4),1,w(if),1))/diag
         i2 = i2 + i1 - n2
         i3 = i3 + n1
         i4 = i4 + n
      end do
      if ( n2==1 ) then
         w(i3) = l(nl)
         w(i4) = zero
         call dcopy(n3,w(i4),0,w(i4),1)
         w(if-1+n) = zero
      end if
      call dscal(n,-one,w(if),1)

      ic = if + n
      id = ic + meq*n

      if ( meq>0 ) then

         !  recover matrix c from upper part of a

         do i = 1 , meq
            call dcopy(n,a(i,1),la,w(ic-1+i),meq)
         end do

         !  recover vector d from upper part of b

         call dcopy(meq,b(1),1,w(id),1)
         call dscal(meq,-one,w(id),1)

      end if

      ig = id + meq

      if ( mineq>0 ) then
         ! recover matrix g from lower part of a
         ! The matrix G(mineq+2*n,m1) is stored at w(ig)
         ! Not all rows will be filled if some of the upper/lower
         ! bounds are unbounded.
         do i = 1 , mineq
            call dcopy(n,a(meq+i,1),la,w(ig-1+i),m1)
         end do
      end if

      ih = ig + m1*n
      iw = ih + mineq + 2*n

      if ( mineq>0 ) then
         ! recover h from lower part of b
         ! The vector H(mineq+2*n) is stored at w(ih)
         call dcopy(mineq,b(meq+1),1,w(ih),1)
         call dscal(mineq,-one,w(ih),1)
      end if

      !  augment matrix g by +i and -i, and,
      !  augment vector h by xl and xu
      !  NaN or infBnd value indicates no bound

      ip = ig + mineq
      il = ih + mineq
      num_unbounded = 0

      do i=1,n
         if (ieee_is_nan(xl(i)) .or. xl(i)<=-infbnd) then
            num_unbounded = num_unbounded + 1
         else
            call update_w(xl(i), one)
         end if
      end do

      do i=1,n
         if (ieee_is_nan(xu(i)) .or. xu(i)>=infbnd) then
            num_unbounded = num_unbounded + 1
         else
            call update_w(xu(i), -one)
         end if
      end do

      call lsei(w(ic),w(id),w(ie),w(if),w(ig),w(ih),max(1,meq),meq,n,n, &
                m1,m1-num_unbounded,n,x,xnorm,w(iw),mode,max_iter_ls,nnls_mode)

      if ( mode==1 ) then
         ! restore lagrange multipliers (only for user-defined variables)
         call dcopy(m,w(iw),1,y(1),1)
         if (n3 > 0) then
            !set rest of the multipliers to nan (they are not used)
            y(m+1) = ieee_value(one, ieee_quiet_nan)
            do i=m+2,m+n3+n3
                y(i) = y(m+1)
            end do
         end if
         call enforce_bounds(x,xl,xu,infbnd)  ! to ensure that bounds are not violated
      end if

      contains

      subroutine update_w(val, fact)
        real(wp),intent(in) :: val  !! xu(i) or xl(i)
        real(wp),intent(in) :: fact !! -1 or 1
        w(il) = fact*val
        do j=1,n
           w(ip + m1*(j-1)) = zero
        end do
        w(ip + m1*(i-1)) = fact
        ip = ip + 1
        il = il + 1
      end subroutine update_w

      end subroutine lsq
!*******************************************************************************

!*******************************************************************************
!>
!  for `mode=1`, the subroutine returns the solution `x` of
!  equality & inequality constrained least squares problem lsei :
!
!  \( \underset{x}{\min} ||E x - f|| \)
!
!  s.t.  \( C x  = d  \) and \( G x \ge h \).
!
!  using QR decomposition & orthogonal basis of nullspace of \( C \).
!
!  The following dimensions of the arrays defining the problem
!  are necessary:
!````
!        dim(c) :   formal (lc,n),    actual (mc,n)
!        dim(d) :   formal (lc  ),    actual (mc  )
!        dim(e) :   formal (le,n),    actual (me,n)
!        dim(f) :   formal (le  ),    actual (me  )
!        dim(g) :   formal (lg,n),    actual (mg,n)
!        dim(h) :   formal (lg  ),    actual (mg  )
!        dim(x) :   formal (n   ),    actual (n   )
!        dim(w) :   2*mc+me+(me+mg)*(n-mc)  for lsei
!                 +(n-mc+1)*(mg+2)+2*mg     for lsi
!        dim(jw):   max(mg,l)
!````
!
!  On entry, the user has to provide the arrays C, d, E, f, G, and h.
!  On return, all arrays will be changed by the subroutine.
!
!### Reference
!  * Chapter 23.6 of Lawson & Hanson: Solving least squares problems.
!
!### History
!  * 18.5.1981, dieter kraft, dfvlr oberpfaffenhofen
!  * 20.3.1987, dieter kraft, dfvlr oberpfaffenhofen

    subroutine lsei(c,d,e,f,g,h,lc,mc,le,me,lg,mg,n,x,xnrm,w,mode,&
                    max_iter_ls,nnls_mode)

    implicit none

    integer,intent(in)                      :: lc
    integer,intent(in)                      :: mc
    integer,intent(in)                      :: le
    integer,intent(in)                      :: me
    integer,intent(in)                      :: lg
    integer,intent(in)                      :: mg
    integer,intent(in)                      :: n
    real(wp),dimension(lc,n),intent(inout)  :: c
    real(wp),dimension(lc)  ,intent(inout)  :: d
    real(wp),dimension(le,n),intent(inout)  :: e
    real(wp),dimension(le)  ,intent(inout)  :: f
    real(wp),dimension(lg,n),intent(inout)  :: g
    real(wp),dimension(lg)  ,intent(inout)  :: h
    real(wp),dimension(n)   ,intent(out)    :: x    !! stores the solution vector
    real(wp),intent(out)                    :: xnrm !! stores the residuum of the solution in euclidian norm
    real(wp),dimension(*)   ,intent(inout)  :: w    !! on return, stores the vector of lagrange multipliers
                                                    !! in its first `mc+mg` elements
    integer,intent(out)                     :: mode !! is a success-failure flag with the following meanings:
                                                    !!
                                                    !! * ***1:*** successful computation,
                                                    !! * ***2:*** error return because of wrong dimensions (`n<1`),
                                                    !! * ***3:*** iteration count exceeded by [[nnls]],
                                                    !! * ***4:*** inequality constraints incompatible,
                                                    !! * ***5:*** matrix `e` is not of full rank,
                                                    !! * ***6:*** matrix `c` is not of full rank,
                                                    !! * ***7:*** rank defect in [[hfti]]
    integer,intent(in) :: max_iter_ls !! maximum number of iterations in the [[nnls]] problem
    integer,intent(in) :: nnls_mode !! which NNLS method to use

    integer :: i , ie, if , ig , iw , j , k , krank , l , mc1
    real(wp) :: t , dum(1)

    mode = 2
    if ( mc<=n ) then
        l = n - mc
        mc1 = mc + 1
        iw = (l+1)*(mg+2) + 2*mg + mc
        ie = iw + mc + 1
        if = ie + me*l
        ig = if + me

        !  triangularize c and apply factors to e and g

        do i = 1 , mc
            j = min(i+1,lc)
            call h12(1,i,i+1,n,c(i,1),lc,w(iw+i),c(j,1),lc,1,mc-i)
            call h12(2,i,i+1,n,c(i,1),lc,w(iw+i),e,le,1,me)
            call h12(2,i,i+1,n,c(i,1),lc,w(iw+i),g,lg,1,mg)
        end do

        !  solve c*x=d and modify f

        mode = 6
        do i = 1 , mc
            if ( abs(c(i,i))<epmach ) return
            x(i) = (d(i)-ddot(i-1,c(i,1),lc,x,1))/c(i,i)
        end do
        mode = 1
        w(mc1) = zero
        !call dcopy(mg-mc,w(mc1),0,w(mc1),1)  ! original code
        call dcopy(mg,w(mc1),0,w(mc1),1)      ! bug fix for when meq = n

        if ( mc/=n ) then

            do i = 1 , me
                w(if-1+i) = f(i) - ddot(mc,e(i,1),le,x,1)
            end do

            !  store transformed e & g

            do i = 1 , me
                call dcopy(l,e(i,mc1),le,w(ie-1+i),me)
            end do
            do i = 1 , mg
                call dcopy(l,g(i,mc1),lg,w(ig-1+i),mg)
            end do

            if ( mg>0 ) then
                !  modify h and solve inequality constrained ls problem

                do i = 1 , mg
                    h(i) = h(i) - ddot(mc,g(i,1),lg,x,1)
                end do
                call lsi(w(ie),w(if),w(ig),h,me,me,mg,mg,l,x(mc1),xnrm,  &
                         w(mc1),mode,max_iter_ls,nnls_mode)

                if ( mc==0 ) return
                t = dnrm2(mc,x,1)
                xnrm = sqrt(xnrm*xnrm+t*t)
                if ( mode/=1 ) return
            else

                ! solve ls without inequality constraints

                mode = 7
                k = max(le,n)
                t = sqrt(epmach)
                call hfti(w(ie),me,me,l,w(if),k,1,t,krank,dum,w,w(l+1))
                xnrm = dum(1)
                call dcopy(l,w(if),1,x(mc1),1)
                if ( krank/=l ) return
                mode = 1
            end if
        end if

        !  solution of original problem and lagrange multipliers

        do i = 1 , me
            f(i) = ddot(n,e(i,1),le,x,1) - f(i)
        end do
        do i = 1 , mc
            d(i) = ddot(me,e(1,i),1,f,1) &
            - ddot(mg,g(1,i),1,w(mc1),1)
        end do

        do i = mc , 1 , -1
            call h12(2,i,i+1,n,c(i,1),lc,w(iw+i),x,1,1,1)
        end do

        do i = mc , 1 , -1
            j = min(i+1,lc)
            w(i) = (d(i)-ddot(mc-i,c(j,i),1,w(j),1))/c(i,i)
        end do
    end if

    end subroutine lsei
!*******************************************************************************

!*******************************************************************************
!>
!  for `mode=1`, the subroutine returns the solution `x` of
!  inequality constrained linear least squares problem:
!
!  \( \underset{x}{\min} ||E x - f|| \)
!
!  s.t. \( G x \ge h \).
!
!  the following dimensions of the arrays defining the problem
!  are necessary:
!````
!     dim(e) :   formal (le,n),    actual (me,n)
!     dim(f) :   formal (le  ),    actual (me  )
!     dim(g) :   formal (lg,n),    actual (mg,n)
!     dim(h) :   formal (lg  ),    actual (mg  )
!     dim(x) :   n
!     dim(w) :   (n+1)*(mg+2) + 2*mg
!     dim(jw):   lg
!````
!
!  on entry, the user has to provide the arrays `e`, `f`, `g`, and `h`.
!  on return, all arrays will be changed by the subroutine.
!
!### Reference
!  * Chapter 23.6 of Lawson & Hanson: Solving least squares problems.
!
!### History
!  * 03.01.1980, dieter kraft: coded
!  * 20.03.1987, dieter kraft: revised to fortran 77

    subroutine lsi(e,f,g,h,le,me,lg,mg,n,x,xnorm,w,mode,max_iter_ls,nnls_mode)

    implicit none

    integer,intent(in)                      :: le
    integer,intent(in)                      :: me
    integer,intent(in)                      :: lg
    integer,intent(in)                      :: mg
    integer,intent(in)                      :: n
    real(wp),dimension(le,n),intent(inout)  :: e
    real(wp),dimension(le)  ,intent(inout)  :: f
    real(wp),dimension(lg,n),intent(inout)  :: g
    real(wp),dimension(lg)  ,intent(inout)  :: h
    real(wp),dimension(n)   ,intent(out)    :: x     !! stores the solution vector
    real(wp),intent(out)                    :: xnorm !! stores the residuum of the solution in euclidian norm
    real(wp),dimension(*)   ,intent(inout)  :: w     !! stores the vector of lagrange multipliers in its first
                                                     !! `mg` elements
    integer,intent(out)                     :: mode  !! is a success-failure flag with the following meanings:
                                                     !!
                                                     !! * ***1:*** successful computation,
                                                     !! * ***2:*** error return because of wrong dimensions (`n<1`),
                                                     !! * ***3:*** iteration count exceeded by [[nnls]],
                                                     !! * ***4:*** inequality constraints incompatible,
                                                     !! * ***5:*** matrix `e` is not of full rank.
    integer,intent(in) :: max_iter_ls !! maximum number of iterations in the [[nnls]] problem
    integer,intent(in) :: nnls_mode !! which NNLS method to use

    integer :: i, j
    real(wp) :: t

    !  qr-factors of e and application to f

    do i = 1 , n
        j = min(i+1,n)
        call h12(1,i,i+1,me,e(1,i),1,t,e(1,j),1,le,n-i)
        call h12(2,i,i+1,me,e(1,i),1,t,f,1,1,1)
    end do

    !  transform g and h to get least distance problem

    mode = 5
    do i = 1 , mg
        do j = 1 , n
            if ( abs(e(j,j))<epmach .or. ieee_is_nan(e(j,j))) return
            g(i,j) = (g(i,j)-ddot(j-1,g(i,1),lg,e(1,j),1))/e(j,j)
        end do
        h(i) = h(i) - ddot(n,g(i,1),lg,f,1)
    end do

    !  solve least distance problem

    call ldp(g,lg,mg,n,h,x,xnorm,w,mode,max_iter_ls,nnls_mode)
    if ( mode==1 ) then

        !  solution of original problem

        call daxpy(n,one,f,1,x,1)
        do i = n , 1 , -1
            j = min(i+1,n)
            x(i) = (x(i)-ddot(n-i,e(i,j),le,x(j),1))/e(i,i)
        end do
        j = min(n+1,me)
        t = dnrm2(me-n,f(j),1)
        xnorm = sqrt(xnorm*xnorm+t*t)
    end if

    end subroutine lsi
!*******************************************************************************

!*******************************************************************************
!>
!  Least distance programming routine.
!  Minimize \( \frac{1}{2} \mathbf{x}^T \mathbf{x}  \) subject to
!  \( \mathbf{G} \mathbf{x} \ge \mathbf{h} \).
!
!  The declared dimension of `w` must be at least `(n+1)*(m+2)+2*m`:
!````
!       first (n+1)*m locs of w = matrix e for problem nnls.
!       next      n+1 locs of w = vector f for problem nnls.
!       next      n+1 locs of w = vector z for problem nnls.
!       next        m locs of w = vector y for problem nnls.
!       next        m locs of w = vector wdual for problem nnls.
!````
!
!### References
!  * C.L. Lawson, R.J. Hanson, 'Solving least squares problems'
!    Prentice Hall, 1974. (revised 1995 edition)
!  * [lawson-hanson](http://www.netlib.org/lawson-hanson/all) from Netlib.
!
!### History
!  * Jacob Williams, refactored into modern Fortran, Jan. 2016.
!
!@note The 1995 version of this routine may have some sort of problem.
!      Using a refactored version of the original routine.

    subroutine ldp(g,mg,m,n,h,x,xnorm,w,mode,max_iter_ls,nnls_mode)

    implicit none

    integer,intent(in)                   :: mg
    integer,intent(in)                   :: m
    integer,intent(in)                   :: n
    real(wp),dimension(mg,n),intent(in)  :: g     !! on entry `g` stores the `m` by `n` matrix of
                                                  !! linear inequality constraints. `g` has first
                                                  !! dimensioning parameter `mg`
    real(wp),dimension(m),intent(in)     :: h     !! the right side of the inequality system.
    real(wp),dimension(n),intent(out)    :: x     !! solution vector `x` if `mode=1`.
    real(wp),dimension(*),intent(inout)  :: w     !! `w` is a one dimensional working space, the length
                                                  !! of which should be at least `(m+2)*(n+1) + 2*m`.
                                                  !! on exit `w` stores the lagrange multipliers
                                                  !! associated with the constraints.
                                                  !! at the solution of problem `ldp`.
    real(wp),intent(out)                 :: xnorm !! euclidian norm of the solution vector
                                                  !! if computation is successful
    integer,intent(out)                  :: mode  !! success-failure flag with the following meanings:
                                                  !!
                                                  !! * ***1:*** successful computation,
                                                  !! * ***2:*** error return because of wrong dimensions (`n<=0`),
                                                  !! * ***3:*** iteration count exceeded by [[nnls]],
                                                  !! * ***4:*** inequality constraints incompatible.
    integer,intent(in) :: max_iter_ls !! maximum number of iterations in the [[nnls]] problem
    integer,intent(in) :: nnls_mode !! which NNLS method to use

    integer :: i , iw , iwdual , iy , iz , j , jf , n1
    real(wp) :: fac , rnorm

    if ( n<=0 ) then
       ! error return.
       mode = 2
    else
        ! state dual problem
        mode = 1
        x = zero
        xnorm = zero
        if (m/=0) then
            iw=0
            do j=1,m
                do i=1,n
                    iw=iw+1
                    w(iw)=g(j,i)
                end do
                iw=iw+1
                w(iw)=h(j)
            end do
            jf=iw+1
            do i=1,n
                iw=iw+1
                w(iw)=zero
            end do
            w(iw+1)=one
            n1=n+1
            iz=iw+2
            iy=iz+n1
            iwdual=iy+m
            ! solve dual problem
            select case (nnls_mode)
            case(1) ! original
                call nnls (w,n1,n1,m,w(jf),w(iy),rnorm,w(iwdual),w(iz),mode,max_iter_ls)
            case(2) ! new version
                call bvls_wrapper(w,n1,n1,m,w(jf),w(iy),rnorm,w(iwdual),w(iz),mode,max_iter_ls)
            case default
                error stop 'invalid nnls_mode'
            end select

            if (mode==1) then
                mode=4
                if (rnorm>zero) then
                    !  compute solution of primal problem
                    fac=one-ddot(m,h,1,w(iy),1)
                    if (ieee_is_nan(fac)) return
                    if (fac>=eps) then
                        mode=1
                        fac=one/fac
                        do j=1,n
                            x(j)=fac*ddot(m,g(1,j),1,w(iy),1)
                        end do
                        xnorm=dnrm2(n,x,1)
                        ! compute lagrange multipliers for primal problem
                        w(1)=zero
                        call dcopy(m,w(1),0,w,1)
                        call daxpy(m,fac,w(iy),1,w,1)
                    end if
                end if
            end if
        end if
    end if

    end subroutine ldp
!*******************************************************************************

!*******************************************************************************
!>
!  Nonnegative least squares algorithm.
!
!  Given an m by n matrix, \(\mathbf{A}\), and an m-vector, \(\mathbf{b}\),
!  compute an n-vector, \(\mathbf{x}\), that solves the least squares problem:
!
!  \( \mathbf{A} \mathbf{x} = \mathbf{b}\) subject to \( \mathbf{x} \ge 0 \)
!
!### References
!  * C.L. Lawson, R.J. Hanson, 'Solving least squares problems'
!    Prentice Hall, 1974. (revised 1995 edition)
!  * [lawson-hanson](http://www.netlib.org/lawson-hanson/all) from Netlib.
!
!### History
!  * Jacob Williams, refactored into modern Fortran, Jan. 2016.

    subroutine nnls(a,mda,m,n,b,x,rnorm,w,zz,mode,max_iter)

    implicit none

    integer,intent(in)                      :: mda     !! first dimensioning parameter for the array `a`.
    integer,intent(in)                      :: n
    real(wp),dimension(mda,n),intent(inout) :: a       !! on entry, contains the `m` by `n`
                                                       !! matrix, `a`. on exit, contains
                                                       !! the product matrix, `q*a`, where `q` is an
                                                       !! `m` by `m` orthogonal matrix generated implicitly by
                                                       !! this subroutine.
    integer,intent(in)                      :: m
    real(wp),dimension(m),intent(inout)     :: b       !! on entry, contains the m-vector `b`. on exit, contains `q*b`.
    real(wp),dimension(n),intent(out)       :: x       !! the solution vector.
    real(wp),intent(out)                    :: rnorm   !! euclidean norm of the residual vector.
    real(wp),dimension(n),intent(inout)     :: w       !! array of working space.  on exit `w` will contain
                                                       !! the dual solution vector. `w` will satisfy `w(i) = 0`
                                                       !! for all `i` in set `p` and `w(i) <= 0` for all `i` in set `z`.
    real(wp),dimension(m),intent(inout)     :: zz      !! an m-array of working space.
    integer,intent(out)                     :: mode    !! this is a success-failure flag with the following meanings:
                                                       !!
                                                       !! * ***1*** the solution has been computed successfully.
                                                       !! * ***2*** the dimensions of the problem are bad. either `m<=0` or `n<=0`.
                                                       !! * ***3*** iteration count exceeded. more than `3*n` iterations.
    integer,intent(in)                      :: max_iter !! maximum number of iterations (if <=0, then `3*n` is used)

    integer :: i,ii,ip,iter,itmax,iz,iz1,iz2,izmax,j,jj,jz,l,npp1,nsetp,rtnkey
    real(wp) :: alpha,asave,cc,sm,ss,t,temp,unorm,up,wmax,ztest
    real(wp),dimension(1) :: dummy
    integer,dimension(n) :: index   !! an integer working array.
                                    !! the contents of this array define the sets
                                    !! `p` and `z` as follows:
                                    !!
                                    !! * `index(1:nsetp) = set p`.
                                    !! * `index(iz1:iz2) = set z`.
                                    !!
                                    !! where: `iz1 = nsetp + 1 = npp1`, `iz2 = n`

    real(wp),parameter :: factor = 0.01_wp

    mode = 1
    if ( m<=0 .or. n<=0 ) then
        mode = 2
        return
    end if
    iter = 0

    if (max_iter<=0) then
        itmax = 3*n
    else
        itmax = max_iter
    end if

    ! initialize the arrays index(1:n) and x(1:n).
    x = zero
    index = [(i, i=1,n)]
    iz2 = n
    iz1 = 1
    nsetp = 0
    npp1 = 1

    main : do

        ! ******  main loop begins here  ******
        ! quit if all coefficients are already in the solution.
        ! or if m cols of a have been triangularized.

        if ( iz1<=iz2 .and. nsetp<m ) then

            ! compute components of the dual (negative gradient) vector w().

            do iz = iz1 , iz2
                j = index(iz)
                sm = zero
                do l = npp1 , m
                    sm = sm + a(l,j)*b(l)
                end do
                w(j) = sm
            end do

            do
                ! find largest positive w(j).
                wmax = zero
                do iz = iz1 , iz2
                    j = index(iz)
                    if ( w(j)>wmax ) then
                        wmax = w(j)
                        izmax = iz
                    end if
                end do

                ! if wmax <= 0. go to termination.
                ! this indicates satisfaction of the kuhn-tucker conditions.
                if ( wmax<=zero ) then
                    call termination()
                    return
                end if

                iz = izmax
                j = index(iz)

                ! the sign of w(j) is ok for j to be moved to set p.
                ! begin the transformation and check new diagonal element to avoid
                ! near linear dependence.

                asave = a(npp1,j)
                call h12(1,npp1,npp1+1,m,a(1,j),1,up,dummy,1,1,0)
                unorm = zero
                if ( nsetp/=0 ) then
                    do l = 1 , nsetp
                        unorm = unorm + a(l,j)**2
                    end do
                end if
                unorm = sqrt(unorm)
                if (abs(a(npp1,j))*factor >= unorm*eps) then

                    ! col j is sufficiently independent.  copy b into zz, update zz
                    ! and solve for ztest ( = proposed new value for x(j) ).

                    do l = 1 , m
                        zz(l) = b(l)
                    end do
                    call h12(2,npp1,npp1+1,m,a(1,j),1,up,zz,1,1,1)
                    ztest = zz(npp1)/a(npp1,j)

                    ! see if ztest is positive
                    if ( ztest>zero ) then

                        ! the index j=index(iz) has been selected to be moved from
                        ! set z to set p. update b, update indices, apply householder
                        ! transformations to cols in new set z, zero subdiagonal elts in
                        ! col j, set w(j)=0.

                        do l = 1 , m
                            b(l) = zz(l)
                        end do

                        index(iz) = index(iz1)
                        index(iz1) = j
                        iz1 = iz1 + 1
                        nsetp = npp1
                        npp1 = npp1 + 1

                        if ( iz1<=iz2 ) then
                            do jz = iz1 , iz2
                                jj = index(jz)
                                call h12(2,nsetp,npp1,m,a(1,j),1,up,a(1,jj),1,mda,1)
                            end do
                        end if

                        if ( nsetp/=m ) then
                            do l = npp1 , m
                                a(l,j) = zero
                            end do
                        end if

                        w(j) = zero
                        ! solve the triangular system.
                        ! store the solution temporarily in zz().
                        rtnkey = 1
                        exit
                    end if
                end if

                ! reject j as a candidate to be moved from set z to set p.
                ! restore a(npp1,j), set w(j)=0., and loop back to test dual
                ! coeffs again.

                a(npp1,j) = asave
                w(j) = zero

            end do

        else
            call termination()
            return
        end if

        ! ******  end of main loop  ******

        secondary : do

            ! the following block of code is used as an internal subroutine
            ! to solve the triangular system, putting the solution in zz().
            do l = 1 , nsetp
                ip = nsetp + 1 - l
                if ( l/=1 ) then
                    do ii = 1 , ip
                        zz(ii) = zz(ii) - a(ii,jj)*zz(ip+1)
                    end do
                end if
                jj = index(ip)
                zz(ip) = zz(ip)/a(ip,jj)
            end do
            if (rtnkey/=1 .and. rtnkey/=2) return

            ! ******  secondary loop begins here ******

            ! iteration counter.

            iter = iter + 1
            if ( iter>itmax ) then
                mode = 3
                !write (*,'(/a)') ' nnls quitting on iteration count.'
                call termination()
                return
            end if

            ! see if all new constrained coeffs are feasible.
            ! if not compute alpha.

            alpha = two
            do ip = 1 , nsetp
                l = index(ip)
                if ( zz(ip)<=zero ) then
                    t = -x(l)/(zz(ip)-x(l))
                    if ( alpha>t ) then
                        alpha = t
                        jj = ip
                    end if
                end if
            end do

            ! if all new constrained coeffs are feasible then alpha will
            ! still = 2.    if so exit from secondary loop to main loop.

            if ( abs(alpha-two)<=zero ) then
                ! ******  end of secondary loop  ******

                do ip = 1 , nsetp
                    i = index(ip)
                    x(i) = zz(ip)
                end do
                ! all new coeffs are positive.  loop back to beginning.
                cycle main
            end if

            ! otherwise use alpha which will be between 0. and 1. to
            ! interpolate between the old x and the new zz.

            do ip = 1 , nsetp
                l = index(ip)
                x(l) = x(l) + alpha*(zz(ip)-x(l))
            end do

            ! modify a and b and the index arrays to move coefficient i
            ! from set p to set z.

            i = index(jj)

            move_p : do
                x(i) = zero

                if ( jj/=nsetp ) then
                    jj = jj + 1
                    do j = jj , nsetp
                        ii = index(j)
                        index(j-1) = ii
                        call g1(a(j-1,ii),a(j,ii),cc,ss,a(j-1,ii))
                        a(j,ii) = zero
                        do l = 1 , n
                            if ( l/=ii ) then
                                ! apply procedure g2 (cc,ss,a(j-1,l),a(j,l))
                                temp = a(j-1,l)
                                a(j-1,l) = cc*temp + ss*a(j,l)
                                a(j,l) = -ss*temp + cc*a(j,l)
                            end if
                        end do
                        ! apply procedure g2 (cc,ss,b(j-1),b(j))
                        temp = b(j-1)
                        b(j-1) = cc*temp + ss*b(j)
                        b(j) = -ss*temp + cc*b(j)
                    end do
                end if

                npp1 = nsetp
                nsetp = nsetp - 1
                iz1 = iz1 - 1
                index(iz1) = i

                ! see if the remaining coeffs in set p are feasible.  they should
                ! be because of the way alpha was determined.
                ! if any are infeasible it is due to round-off error.  any
                ! that are nonpositive will be set to zero
                ! and moved from set p to set z.

                do jj = 1 , nsetp
                    i = index(jj)
                    if ( x(i)<=zero ) cycle move_p
                end do
                exit move_p

            end do move_p

            ! copy b( ) into zz( ).  then solve again and loop back.

            do i = 1 , m
                zz(i) = b(i)
            end do
            rtnkey = 2

        end do secondary

    end do main

    contains

    subroutine termination()
        !! come to here for termination.
        !! compute the norm of the final residual vector.

        sm = zero
        if ( npp1<=m ) then
            do i = npp1 , m
                sm = sm + b(i)**2
            end do
        else
            do j = 1 , n
                w(j) = zero
            end do
        end if
        rnorm = sqrt(sm)
    end subroutine termination

    end subroutine nnls
!*******************************************************************************

!*******************************************************************************
!>
!  Rank-deficient least squares algorithm using
!  householder forward triangulation with column interchanges.
!
!### References
!  * C.L. Lawson, R.J. Hanson, 'Solving least squares problems'
!    Prentice Hall, 1974. (revised 1995 edition)
!  * [lawson-hanson](http://www.netlib.org/lawson-hanson/all) from Netlib.
!
!### History
!  * Jacob Williams, refactored into modern Fortran, Jan. 2016.

    subroutine hfti(a,mda,m,n,b,mdb,nb,tau,krank,rnorm,h,g)

    implicit none

    integer,intent(in)                       :: mda   !! the first dimensioning parameter of matrix `a` (mda >= m).
    integer,intent(in)                       :: m
    integer,intent(in)                       :: n
    integer,intent(in)                       :: mdb   !! first dimensioning parameter of matrix `b` (mdb>=max(m,n))
    integer,intent(in)                       :: nb
    real(wp),dimension(mda,n),intent(inout)  :: a     !! the array `a` initially contains the \( m \times n \) matrix \(\mathbf{A}\)
                                                      !! of the least squares problem \( \mathbf{A} \mathbf{x} = \mathbf{b} \).
                                                      !! either `m >= n` or `m < n` is permitted.
                                                      !! there is no restriction on the rank of `a`.
                                                      !! the matrix `a` will be modified by the subroutine.
    real(wp),intent(in)                      :: tau   !! absolute tolerance parameter for pseudorank
                                                      !! determination, provided by the user.
    integer,intent(out)                      :: krank !! pseudorank of `a`, set by the subroutine.
    real(wp),dimension(nb),intent(out)       :: rnorm !! on exit, `rnorm(j)` will contain the euclidian
                                                      !! norm of the residual vector for the problem
                                                      !! defined by the `j-th` column vector of the array `b`.
    real(wp),dimension(n),intent(inout)      :: h     !! array of working space
    real(wp),dimension(n),intent(inout)      :: g     !! array of working space
    real(wp),dimension(mdb,nb),intent(inout) :: b     !! if `nb = 0` the subroutine will make no reference
                                                      !! to the array `b`. if `nb > 0` the array `b` must
                                                      !! initially contain the `m x nb` matrix `b` of the
                                                      !! the least squares problem `ax = b` and on return
                                                      !! the array `b` will contain the `n x nb` solution `x`.

    integer  :: i, ii, ip1, j, jb, jj, k, kp1, l, ldiag, lmax
    real(wp) :: hmax, sm, tmp
    logical  :: need_lmax
    integer,dimension(n) :: ip    !! integer array of working space
                                  !! recording permutation indices of column vectors

    real(wp),parameter :: factor = 0.001_wp

    k = 0
    ldiag = min(m,n)
    if ( ldiag<=0 ) then

        ! the solution vectors, x, are now
        ! in the first  n  rows of the array b(,).
        krank = k
        return

    else

        do j = 1 , ldiag

            need_lmax = .true.
            if ( j/=1 ) then
                ! update squared column lengths and find lmax
                lmax = j
                do l = j , n
                    h(l) = h(l) - a(j-1,l)**2
                    if ( h(l)>h(lmax) ) lmax = l
                end do
                if (factor*h(lmax) >= hmax*eps) need_lmax = .false.

            end if

            if (need_lmax) then
                ! compute squared column lengths and find lmax
                lmax = j
                do l = j , n
                    h(l) = zero
                    do i = j , m
                        h(l) = h(l) + a(i,l)**2
                    end do
                    if ( h(l)>h(lmax) ) lmax = l
                end do
                hmax = h(lmax)
            end if

            ! lmax has been determined

            ! do column interchanges if needed.

            ip(j) = lmax
            if ( ip(j)/=j ) then
                do i = 1 , m
                    tmp = a(i,j)
                    a(i,j) = a(i,lmax)
                    a(i,lmax) = tmp
                end do
                h(lmax) = h(j)
            end if

            ! compute the j-th transformation and apply it to a and b.

            call h12(1,j,j+1,m,a(1,j),1,h(j),a(1,j+1),1,mda,n-j)
            call h12(2,j,j+1,m,a(1,j),1,h(j),b,1,mdb,nb)

        end do

        ! determine the pseudorank, k, using the tolerance, tau.
        do j=1,ldiag
            if (abs(a(j,j))<=tau) exit
        end do
        k=j-1
        kp1=j

    end if

    ! compute the norms of the residual vectors.

    if ( nb>0 ) then
        do jb = 1 , nb
            tmp = zero
            if ( kp1<=m ) then
                do i = kp1 , m
                    tmp = tmp + b(i,jb)**2
                end do
            end if
            rnorm(jb) = sqrt(tmp)
        end do
    end if
    ! special for pseudorank = 0
    if ( k>0 ) then

        ! if the pseudorank is less than n compute householder
        ! decomposition of first k rows.

        if ( k/=n ) then
            do ii = 1 , k
                i = kp1 - ii
                call h12(1,i,kp1,n,a(i,1),mda,g(i),a,mda,1,i-1)
            end do
        end if

        if ( nb>0 ) then
            do jb = 1 , nb

                ! solve the k by k triangular system.

                do l = 1 , k
                    sm = zero
                    i = kp1 - l
                    if ( i/=k ) then
                        ip1 = i + 1
                        do j = ip1 , k
                            sm = sm + a(i,j)*b(j,jb)
                        end do
                    end if
                    b(i,jb) = (b(i,jb)-sm)/a(i,i)
                end do

                ! complete computation of solution vector.

                if ( k/=n ) then
                    do j = kp1 , n
                        b(j,jb) = zero
                    end do
                    do i = 1 , k
                        call h12(2,i,kp1,n,a(i,1),mda,g(i),b(1,jb),1,mdb,1)
                    end do
                end if

                ! re-order the solution vector to compensate for the
                ! column interchanges.

                do jj = 1 , ldiag
                    j = ldiag + 1 - jj
                    if ( ip(j)/=j ) then
                        l = ip(j)
                        tmp = b(l,jb)
                        b(l,jb) = b(j,jb)
                        b(j,jb) = tmp
                    end if
                end do
            end do
        end if

    else if ( nb>0 ) then

        do jb = 1 , nb
            do i = 1 , n
                b(i,jb) = zero
            end do
        end do

    end if
    krank = k

    end subroutine hfti
!*******************************************************************************

!*******************************************************************************
!>
!  Construction and/or application of a single
!  householder transformation \( Q = I + u(u^t)/b \).
!
!### References
!  * C.L. Lawson, R.J. Hanson, 'Solving least squares problems'
!    Prentice Hall, 1974. (revised 1995 edition)
!  * [lawson-hanson](http://www.netlib.org/lawson-hanson/all) from Netlib.
!
!### History
!  * Jacob Williams, refactored into modern Fortran, Jan. 2016.

    subroutine h12(mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)

    implicit none

    integer,intent(in)                      :: mode   !! `1` or `2` -- selects algorithm ***h1*** to construct and apply a
                                                      !! householder transformation, or algorithm ***h2*** to apply a
                                                      !! previously constructed transformation.
    integer,intent(in)                      :: lpivot !! the index of the pivot element
    integer,intent(in)                      :: l1     !! if `l1 <= m` the transformation will be constructed to
                                                      !! zero elements indexed from `l1` through `m`.
                                                      !! if `l1 > m` the subroutine does an identity transformation.
    integer,intent(in)                      :: m      !! see `li`.
    integer,intent(in)                      :: iue    !! see `u`.
    real(wp),dimension(iue,*),intent(inout) :: u      !! on entry with `mode = 1`, `u` contains the pivot
                                                      !! vector.  `iue` is the storage increment between elements.
                                                      !! on exit when `mode = 1`, `u` and `up` contain quantities
                                                      !! defining the vector `u` of the householder transformation.
                                                      !! on entry with `mode = 2`, `u` and `up` should contain
                                                      !! quantities previously computed with `mode = 1`.  these will
                                                      !! not be modified during the entry with `mode = 2`.
                                                      !! `dimension[u(iue,m)]`
    real(wp),intent(inout)                  :: up     !! see `u`.
    real(wp),dimension(*),intent(inout)     :: c      !! on entry with `mode = 1 or 2`, `c` contains a matrix which
                                                      !! will be regarded as a set of vectors to which the
                                                      !! householder transformation is to be applied.
                                                      !! on exit `c` contains the set of transformed vectors.
    integer,intent(in)                      :: ice    !! storage increment between elements of vectors in `c`.
    integer,intent(in)                      :: icv    !! storage increment between vectors in `c`.
    integer,intent(in)                      :: ncv    !! number of vectors in `c` to be transformed. if `ncv <= 0`
                                                      !! no operations will be done on `c`.

    integer  :: i, i2, i3, i4, incr, j
    real(wp) :: b, cl, clinv, sm

    if ( 0>=lpivot .or. lpivot>=l1 .or. l1>m ) return
    cl = abs(u(1,lpivot))
    if ( mode/=2 ) then
        ! construct the transformation.
        do j = l1 , m
            cl = max(abs(u(1,j)),cl)
        end do
        if ( cl<=zero ) return
        clinv = one/cl
        sm = (u(1,lpivot)*clinv)**2
        do j = l1 , m
            sm = sm + (u(1,j)*clinv)**2
        end do
        cl = cl*sqrt(sm)
        if ( u(1,lpivot)>zero ) cl = -cl
        up = u(1,lpivot) - cl
        u(1,lpivot) = cl
    else if ( cl<=zero ) then
        return
    end if

    if ( ncv>0 ) then
        ! apply the transformation i+u*(u**t)/b to c.
        b = up*u(1,lpivot)
        ! b must be nonpositive here.
        if ( b<zero ) then
            b = one/b
            i2 = 1 - icv + ice*(lpivot-1)
            incr = ice*(l1-lpivot)
            do j = 1 , ncv
                i2 = i2 + icv
                i3 = i2 + incr
                i4 = i3
                sm = c(i2)*up
                do i = l1 , m
                    sm = sm + c(i3)*u(1,i)
                    i3 = i3 + ice
                end do
                if ( abs(sm)>zero ) then
                    sm = sm*b
                    c(i2) = c(i2) + sm*up
                    do i = l1 , m
                        c(i4) = c(i4) + sm*u(1,i)
                        i4 = i4 + ice
                    end do
                end if
            end do
        end if
    end if

    end subroutine h12
!*******************************************************************************

!*******************************************************************************
!>
!  Compute orthogonal rotation matrix.
!
!  Compute matrix $$ \left[ \begin{array}{cc} c & s \\ -s & c \end{array} \right] $$
!  so that
!
!  $$
!  \left[ \begin{array}{cc} c & s \\ -s & c \end{array} \right]
!  \left[ \begin{array}{c} a \\ b \end{array} \right]  =
!  \left[ \begin{array}{c} \sqrt{a^2+b^2} \\ 0 \end{array} \right]
!  $$
!
!  Compute \( \sigma = \sqrt{a^2+b^2} \)
!
!  \( \sigma \) is computed last to allow for the possibility that
!  `sig` may be in the same location as `a` or `b`.
!
!### References
!  * C.L. Lawson, R.J. Hanson, 'Solving least squares problems'
!    Prentice Hall, 1974. (revised 1995 edition)
!  * [lawson-hanson](http://www.netlib.org/lawson-hanson/all) from Netlib.
!
!### History
!  * Jacob Williams, refactored into modern Fortran, Jan. 2016.

    subroutine g1(a,b,c,s,sig)

    implicit none

    real(wp)             :: a
    real(wp)             :: b
    real(wp)             :: sig
    real(wp),intent(out) :: c
    real(wp),intent(out) :: s

    real(wp) :: xr, yr

    if ( abs(a)>abs(b) ) then
        xr = b/a
        yr = sqrt(one+xr**2)
        c = sign(one/yr,a)
        s = c*xr
        sig = abs(a)*yr
    else
        if ( abs(b)>zero ) then
            xr = a/b
            yr = sqrt(one+xr**2)
            s = sign(one/yr,b)
            c = s*xr
            sig = abs(b)*yr
        else
            sig = zero
            c = zero
            s = one
        end if
    end if

    end subroutine g1
!*******************************************************************************

!*******************************************************************************
!>
!  \(LDL^T\) - rank-one - update
!
!### Purpose:
!
!  Updates the \(LDL^T\) factors of matrix \(A\)
!  by rank-one matrix \(\sigma z z^T \).
!
!### Reference
!  * R. Fletcher, M.J.D. Powell,
!    "[On the modification of LDL' factorization](http://www.ams.org/journals/mcom/1974-28-128/S0025-5718-1974-0359297-1/S0025-5718-1974-0359297-1.pdf)".
!    Mathematics of Computation Vol. 28, No. 128, p. 1067-1087, October 1974.
!
!### History
!  * D. Kraft, DFVLR - institut fuer dynamik der flugsysteme
!    d-8031  oberpfaffenhofen
!  * Status: 15. january 1980

    subroutine ldl(n,a,z,sigma,w)

    implicit none

    integer,intent(in)                  :: n     !! order of the coefficient matrix `a`
    real(wp),intent(in)                 :: sigma !! scalar factor by which the modifying dyade \(z z^T\) is multiplied.
    real(wp),dimension(*),intent(inout) :: a     !! ***In:*** positive definite matrix of dimension `n`;
                                                 !! only the lower triangle is used and is stored column by
                                                 !! column as one dimensional array of dimension `n*(n+1)/2`.
                                                 !!
                                                 !! ***Out:*** updated \(LDL^T\) factors
    real(wp),dimension(*),intent(inout) :: w     !! working array of dimension `n` (used only if \( \sigma \lt 0 \) ).
    real(wp),dimension(*),intent(inout) :: z     !! vector of dimension `n` of updating elements.

    integer :: i , ij , j
    real(wp) :: t , v , u , tp , beta , alpha , delta , gamma

    if ( abs(sigma)>zero ) then
        ij = 1
        t = one/sigma
        if ( sigma<=zero ) then
            ! prepare negative update
            do i = 1 , n
                w(i) = z(i)
            end do
            do i = 1 , n
                v = w(i)
                t = t + v*v/a(ij)
                do j = i + 1 , n
                    ij = ij + 1
                    w(j) = w(j) - v*a(ij)
                end do
                ij = ij + 1
            end do
            if ( t>=zero ) t = epmach/sigma
            do i = 1 , n
                j = n + 1 - i
                ij = ij - i
                u = w(j)
                w(j) = t
                t = t - u*u/a(ij)
            end do
        end if
        ! here updating begins
        do i = 1 , n
            v = z(i)
            delta = v/a(ij)
            if ( sigma<zero ) tp = w(i)
            if ( sigma>zero ) tp = t + delta*v
            alpha = tp/t
            a(ij) = alpha*a(ij)
            if ( i==n ) return
            beta = delta/tp
            if ( alpha>four ) then
                gamma = t/tp
                do j = i + 1 , n
                    ij = ij + 1
                    u = a(ij)
                    a(ij) = gamma*u + beta*z(j)
                    z(j) = z(j) - v*u
                end do
            else
                do j = i + 1 , n
                    ij = ij + 1
                    z(j) = z(j) - v*a(ij)
                    a(ij) = a(ij) + beta*z(j)
                end do
            end if
            ij = ij + 1
            t = tp
        end do
    end if

    end subroutine ldl
!*******************************************************************************

!*******************************************************************************
!>
!  Linesearch without derivatives (used by [[slsqp]] if `linesearch_mode=2`).
!  Returns the abscissa approximating the point where `f` attains a minimum.
!
!### purpose:
!
!  to find the argument linmin where the function `f` takes it's minimum
!  on the interval `ax`, `bx`. It uses a combination of golden section
!  and successive quadratic interpolation.
!
!### Reference
!
!  This function subprogram is a slightly modified version of the
!  ALGOL 60 procedure `localmin` given in R.P. Brent:
!  "[Algorithms for minimization without derivatives](https://maths-people.anu.edu.au/~brent/pub/pub011.html)",
!  Prentice-Hall (1973).
!
!### History
!
!  * Kraft, D., DFVLR - institut fuer dynamik der flugsysteme
!    d-8031  oberpfaffenhofen
!  * status: 31. august 1984
!  * Jacob Williams, Jan 2016, Refactored into modern Fortran.
!    Added saved variables as `inout`s to make the routine thread-safe.

    real(wp) function linmin(mode,ax,bx,f,tol,&
                                a,b,d,e,p,q,r,u,v,&
                                w,x,m,fu,fv,fw,fx,tol1,tol2)

    implicit none

    integer,intent(inout) :: mode  !! controls reverse communication
                                   !! must be set to 0 initially, returns with intermediate
                                   !! values 1 and 2 which must not be changed by the user,
                                   !! ends with convergence with value 3.
    real(wp) :: f                  !! function value at `linmin` which is to be brought in by
                                   !! reverse communication controlled by `mode`
    real(wp),intent(in) :: tol     !! desired length of interval of uncertainty of final result
    real(wp),intent(in) :: ax      !! left endpoint of initial interval
    real(wp),intent(in) :: bx      !! right endpoint of initial interval
    real(wp),intent(inout) :: a,b,d,e,p,q,r,u,v,w,x,m,fu,fv,fw,fx,tol1,tol2

    real(wp),parameter :: c = (3.0_wp-sqrt(5.0_wp))/2.0_wp  !! golden section ratio = `0.381966011`
    real(wp),parameter :: sqrteps = sqrt(eps)  !! square root of machine precision

    select case (mode)
    case(1)

        ! main loop starts here

        fx = f
        fv = fx
        fw = fv

    case(2)

        fu = f
        ! update a, b, v, w, and x
        if ( fu>fx ) then
            if ( u<x ) a = u
            if ( u>=x ) b = u
            if ( fu<=fw .or. abs(w-x)<=zero ) then
                v = w
                fv = fw
                w = u
                fw = fu
            else if ( fu<=fv .or. abs(v-x)<=zero .or. abs(v-w)<=zero ) then
                v = u
                fv = fu
            end if
        else
            if ( u>=x ) a = x
            if ( u<x ) b = x
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
        end if

    case default

        ! initialization
        a = ax
        b = bx
        e = zero
        v = a + c*(b-a)
        w = v
        x = w
        linmin = x
        mode = 1
        return

    end select

    m = 0.5_wp*(a+b)
    tol1 = sqrteps*abs(x) + tol
    tol2 = tol1 + tol1

    ! test convergence

    if ( abs(x-m)<=tol2-0.5_wp*(b-a) ) then
        ! end of main loop
        linmin = x
        mode = 3
    else
        r = zero
        q = r
        p = q
        if ( abs(e)>tol1 ) then
            ! fit parabola
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q - (x-w)*r
            q = q - r
            q = q + q
            if ( q>zero ) p = -p
            if ( q<zero ) q = -q
            r = e
            e = d
        end if

        ! is parabola acceptable
        if ( abs(p)>=0.5_wp*abs(q*r) .or. p<=q*(a-x) .or. p>=q*(b-x) ) then
            ! golden section step
            if ( x>=m ) e = a - x
            if ( x<m ) e = b - x
            d = c*e
        else
            ! parabolic interpolation step
            d = p/q
            ! f must not be evaluated too close to a or b
            if ( u-a<tol2 ) d = sign(tol1,m-x)
            if ( b-u<tol2 ) d = sign(tol1,m-x)
        end if

        ! f must not be evaluated too close to x
        if ( abs(d)<tol1 ) d = sign(tol1,d)
        u = x + d
        linmin = u
        mode = 2

    end if

    end function linmin
!*******************************************************************************

!*******************************************************************************
!>
!  enforce the bound constraints on `x`.

    subroutine enforce_bounds(x,xl,xu,infbnd)

    implicit none

    real(wp),dimension(:),intent(inout) :: x   !! optimization variable vector
    real(wp),dimension(:),intent(in)    :: xl  !! lower bounds (must be same dimension as `x`)
    real(wp),dimension(:),intent(in)    :: xu  !! upper bounds (must be same dimension as `x`)
    real(wp),intent(in) :: infbnd !! "infinity" for the upper and lower bounds.
                                  !! Note that `NaN` may also be used to indicate no bound.

    where (x<xl .and. xl>-infbnd .and. .not. ieee_is_nan(xl))
        x = xl
    elsewhere (x>xu .and. xu<infbnd .and. .not. ieee_is_nan(xu))
        x = xu
    end where

    end subroutine enforce_bounds
!*******************************************************************************

!*******************************************************************************
!>
!  Destructor for [[slsqpb_data]] type.

    subroutine destroy_slsqpb_data(me)

    implicit none

    class(slsqpb_data),intent(out) :: me

    end subroutine destroy_slsqpb_data
!*******************************************************************************

!*******************************************************************************
!>
!  Destructor for [[linmin_data]] type.

    subroutine destroy_linmin_data(me)

    implicit none

    class(linmin_data),intent(out) :: me

    end subroutine destroy_linmin_data
!*******************************************************************************

!*******************************************************************************
    end module slsqp_core
!*******************************************************************************

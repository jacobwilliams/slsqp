!*******************************************************************************
!>
!  module for the SLSQP optimization method.

    module slsqp_module

    use support_module

    implicit none

    type :: linmin_data
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
    end type linmin_data

    type :: slsqpb_data
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
    end type slsqpb_data

    type,public :: slsqp_solver

        private

        integer :: n        = 0
        integer :: m        = 0
        integer :: meq      = 0
        integer :: max_iter = 0    !! maximum number of iterations

        real(wp) :: acc = 0.0_wp   !! accuracy tolerance

        real(wp),dimension(:),allocatable :: xl  !! lower bound on x
        real(wp),dimension(:),allocatable :: xu  !! upper bound on x

        integer :: l_w = 0 !! size of `work`
        real(wp),dimension(:),allocatable :: w !! real work array

        integer :: l_jw = 0 !! size of `jwork`
        integer,dimension(:),allocatable :: jw  !! integer work array

        procedure(func),pointer :: f => null()  !! problem function subroutine
        procedure(grad),pointer :: g => null()  !! gradient subroutine

        integer :: linesearch_mode = 1  !! linesearch mode.
                                        !! `1` = inexact (Armijo) linesearch,
                                        !! `2` = exact linesearch.
        type(linmin_data),allocatable :: linmin !! data formerly within [[linmin]].
                                                !! Only used when `linesearch_mode=2`
        type(slsqpb_data) :: slsqpb  !! data formerly within [[slsqpb]].

    contains

        private

        procedure,public :: initialize => initialize_slsqp
        procedure,public :: destroy    => destroy_slsqp
        procedure,public :: optimize   => slsqp_wrapper

    end type slsqp_solver

    abstract interface
        subroutine func(me,x,f,c)  !! function computation
            import :: wp,slsqp_solver
            implicit none
            class(slsqp_solver),intent(inout) :: me
            real(wp),dimension(:),intent(in)  :: x  !! optimization variable vector
            real(wp),intent(out)              :: f  !! value of the objective function
            real(wp),dimension(:),intent(out) :: c  !! the constraint vector `dimension(m)`,
                                                    !! equality constraints (if any) first.
        end subroutine func
        subroutine grad(me,x,g,a)
            import :: wp,slsqp_solver
            implicit none
            class(slsqp_solver),intent(inout)   :: me
            real(wp),dimension(:),intent(in)    :: x  !! optimization variable vector
            real(wp),dimension(:),intent(out)   :: g  !! objective function partials w.r.t x `dimension(n)`
            real(wp),dimension(:,:),intent(out) :: a  !! gradient matrix of constraints w.r.t. x `dimension(m,n)`
        end subroutine grad
    end interface

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the [[slsqp_solver]] class.  see [[slsqp]] for more details.

    subroutine initialize_slsqp(me,n,m,meq,max_iter,acc,f,g,xl,xu,linesearch_mode,status_ok)

    implicit none

    class(slsqp_solver),intent(inout) :: me
    integer,intent(in)                :: n               !! the number of varibles, n >= 1
    integer,intent(in)                :: m               !! total number of constraints, m >= 0
    integer,intent(in)                :: meq             !! number of equality constraints, meq >= 0
    integer,intent(in)                :: max_iter        !! maximum number of iterations
    procedure(func)                   :: f               !! problem function
    procedure(grad)                   :: g               !! function to compute gradients
    real(wp),dimension(n),intent(in)  :: xl              !! lower bound on `x`
    real(wp),dimension(n),intent(in)  :: xu              !! upper bound on `x`
    real(wp),intent(in)               :: acc             !! accuracy
    integer,intent(in)                :: linesearch_mode !! 1 = inexact, 2 = exact
    logical,intent(out)               :: status_ok       !! will be false if there were errors

    integer :: n1,mineq

    status_ok = .false.
    call me%destroy()

    if (size(xl)/=size(xu) .or. size(xl)/=n) then
        write(*,*) 'error: invalid upper or lower bound vector size'
    else if (meq<0 .or. meq>m) then
        write(*,*) 'error: invalid meq value:', meq
    else if (m<0) then
        write(*,*) 'error: invalid m value:', m
    else if (n<1) then
        write(*,*) 'error: invalid n value:', n
    else if (any(xl>xu)) then
        write(*,*) 'error: lower bounds must be <= upper bounds.'
    else

        !two linesearch modes:
        select case (linesearch_mode)
        case(1)     !inexact
            me%linesearch_mode = linesearch_mode
        case(2)     !exact
            me%linesearch_mode = linesearch_mode
            allocate(me%linmin)
        case default
            write(*,*) 'error: invalid linesearch_mode (must be 1 or 2): ',linesearch_mode
            call me%destroy()
            return
        end select

        status_ok = .true.
        me%n = n
        me%m = m
        me%meq = meq
        me%max_iter = max_iter
        me%acc = acc
        me%f => f
        me%g => g

        allocate(me%xl(n)); me%xl = xl
        allocate(me%xu(n)); me%xu = xu

        !work arrays:
        n1 = n+1
        mineq = m - meq + 2*n1
        me%l_w = n1*(n1+1) + meq*(n1+1) + mineq*(n1+1) + &   !for lsq
                 (n1-meq+1)*(mineq+2) + 2*mineq        + &   !for lsi
                 (n1+mineq)*(n1-meq) + 2*meq + n1      + &   !for lsei
                  n1*n/2 + 2*m + 3*n +3*n1 + 1               !for slsqpb
        allocate(me%w(me%l_w))
        me%w = 0.0_wp

        me%l_jw = mineq
        allocate(me%jw(me%l_jw))
        me%jw = 0

    end if

    end subroutine initialize_slsqp
!*******************************************************************************

!*******************************************************************************
!>
!  destructor for [[slsqp_solver]].

    subroutine destroy_slsqp(me)

    implicit none

    class(slsqp_solver),intent(out) :: me

    end subroutine destroy_slsqp
!*******************************************************************************

!*******************************************************************************
!>
!  main routine for calling [[slsqp]].

    subroutine slsqp_wrapper(me,x,istat)

    implicit none

    class(slsqp_solver),intent(inout)   :: me
    real(wp),dimension(:),intent(inout) :: x        !! in: initialize optimization variables,
                                                    !! out: solution.
    integer,intent(out)                 :: istat    !! status code

    real(wp)                               :: f        !! objective function
    real(wp),dimension(max(1,me%m))        :: c        !! constraint vector
    real(wp),dimension(max(1,me%m),me%n+1) :: a        !! a matrix for slsqp
    real(wp),dimension(me%n+1)             :: g        !! g matrix for slsqp
    real(wp),dimension(me%m)               :: cvec     !! constraint vector
    real(wp),dimension(me%n)               :: dfdx     !! objective function partials
    real(wp),dimension(me%m,me%n)          :: dcdx     !! constraint partials
    integer :: i,mode,la,iter
    real(wp) :: acc

    !check setup:
    if (size(x)/=me%n) error stop 'invalid size(x) in slsqp_wrapper'

    !initialize:
    i    = 0
    iter = me%max_iter
    la   = max(1,me%m)
    mode = 0
    a    = 0.0_wp
    g    = 0.0_wp
    c    = 0.0_wp

    !linesearch:
    select case(me%linesearch_mode)
    case(1)     !inexact (armijo-type linesearch)
        acc = abs(me%acc)
    case(2)     !exact
        acc = -abs(me%acc)
    case default
        error stop 'invalid linesearch_mode in slsqp_wrapper'
    end select

    !main solver loop:
    do

        if (mode==0 .or. mode==1) then  !function evaluation (f&c)
            call me%f(x,f,cvec)
            c(1:me%m)   = cvec

            !write(*,*) ''
            !write(*,*) 'func'       !........
            !write(*,*) 'x=',x
            !write(*,*) 'f=',f
            !write(*,*) 'c=',c

            write(*,*) i,x,f,norm2(c)
            i = i + 1

            !note: not really an iteration (can also be during exact linesearch)

        end if

        if (mode==0 .or. mode==-1) then  !gradient evaluation (g&a)
            call me%g(x,dfdx,dcdx)
            g(1:me%n)        = dfdx
            a(1:me%m,1:me%n) = dcdx

            !write(*,*) ''
            !write(*,*) 'grad'       !........
            !write(*,*) 'x=',x
            !write(*,*) 'g=',g
            !write(*,*) 'a=',a

            !write(*,*) i,x,f,norm2(c)
            !i = i + 1

        end if

        !main routine:
        call slsqp(me%m,me%meq,la,me%n,x,me%xl,me%xu,&
                    f,c,g,a,acc,iter,mode,&
                    me%w,me%l_w,me%jw,me%l_jw,&
                    me%slsqpb,me%linmin)

        select case (mode)
        case(0) !required accuracy for solution obtained
            write(*,*) ''
            write(*,*) 'solution: ',x
            write(*,*) ''
            exit
        case(1,-1)
            !continue to next call
        case(2);
            write(*,*) 'number of equality contraints larger than n'
            exit
        case(3);
            write(*,*) 'more than 3*n iterations in lsq subproblem'
            exit
        case(4);
            write(*,*) 'inequality constraints incompatible'
            exit
        case(5);
            write(*,*) 'singular matrix e in lsq subproblem'
            exit
        case(6);
            write(*,*) 'singular matrix c in lsq subproblem'
            exit
        case(7);
            write(*,*) 'rank-deficient equality constraint subproblem hfti'
            exit
        case(8);
            write(*,*) 'positive directional derivative for linesearch'
            exit
        case(9);
            write(*,*) 'more than iter iterations in sqp'
            exit
        case default
            write(*,*) 'unknown slsqp error'
            exit
        end select

    end do

    istat = mode

    end subroutine slsqp_wrapper
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
!  \( f(x) \)
!
!  **subject to**
!
!  \( c_j (x) = 0 \)        , \( j = 1,...,meq   \)
!
!  \(c_j (x) \ge 0 \)       , \( j = meq+1,...,m \)
!
!  \(xl_i \le x_i \le xu_i \) , \( i = 1,...,n     \)
!
!  the algorithm implements the method of Han and Powell
!  with BFGS-update of the b-matrix and L1-test function
!  within the steplength algorithm.
!
!# Reference
!   * Dieter Kraft: "a software package for sequential quadratic programming",
!     dfvlr-fb 88-28, 1988
!
!# History
!
!   * implemented by: Dieter Kraft, DFVLR oberpfaffenhofen
!   * date: april - october, 1981.
!   * december, 31-st, 1984.
!   * march   , 21-st, 1987, revised to fortan 77
!   * march   , 20-th, 1989, revised to ms-fortran
!   * april   , 14-th, 1989, hesse   in-line coded
!   * february, 28-th, 1991, fortran/2 version 1.04 accepts statement functions
!   * march   ,  1-st, 1991, tested with salford ftn77/386 compiler vers 2.40 in protected mode
!   * january, 2016 : refactoring into modern fortran by jacob williams
!
!# License
!
!  Copyright 1991: Dieter Kraft, FHM
!  BSD license
!
!@note: f,c,g,a must all be set by the user before each call.

    subroutine slsqp(m,meq,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,w,l_w, &
                     jw,l_jw,sdat,ldat)

    implicit none

    integer,intent(in) :: m                     !! is the total number of constraints, m >= 0
    integer,intent(in) :: meq                   !! is the number of equality constraints, meq >= 0
    integer,intent(in) :: la                    !! see a, la >= max(m,1)
    integer,intent(in) :: n                     !! is the number of varibles, n >= 1
    real(wp),dimension(n),intent(inout) :: x    !! x() stores the current iterate of the n vector x
                                                !! on entry x() must be initialized. on exit x()
                                                !! stores the solution vector x if mode = 0.
    real(wp),dimension(n),intent(in) :: xl      !! xl() stores an n vector of lower bounds xl to x.
    real(wp),dimension(n),intent(in) :: xu      !! xu() stores an n vector of upper bounds xu to x.
    real(wp),intent(in) :: f                    !! is the value of the objective function.
    real(wp),dimension(la),intent(in) :: c      !! c() stores the m vector c of constraints,
                                                !! equality constraints (if any) first.
                                                !! dimension of c must be greater or equal la,
                                                !! which must be greater or equal max(1,m).
    real(wp),dimension(n+1),intent(in) :: g     !! g() stores the n vector g of partials of the
                                                !! objective function; dimension of g must be
                                                !! greater or equal n+1.
    real(wp),dimension(la,n+1),intent(in) ::  a !! the la by n + 1 array a() stores
                                                !!  the m by n matrix a of constraint normals.
                                                !!  a() has first dimensioning parameter la,
                                                !!  which must be greater or equal max(1,m).
    real(wp),intent(inout) :: acc   !! abs(acc) controls the final accuracy.
                                    !! if acc < zero an exact linesearch is performed,
                                    !! otherwise an armijo-type linesearch is used.
    integer,intent(inout) :: iter   !! prescribes the maximum number of iterations.
                                    !! on exit iter indicates the number of iterations.
    integer,intent(inout) :: mode   !! mode controls calculation:
                                    !! reverse communication is used in the sense that
                                    !! the program is initialized by `mode = 0`; then it is
                                    !! to be called repeatedly by the user until a return
                                    !! with `mode /= abs(1)` takes place.
                                    !! if `mode = -1` gradients have to be calculated,
                                    !! while with `mode = 1` functions have to be calculated.
                                    !! mode must not be changed between subsequent calls of sqp.
                                    !! **evaluation modes**:
                                    !!    * -1 *: gradient evaluation, (g&a)
                                    !!    *  0 *: *on entry*: initialization, (f,g,c&a),
                                    !!            *on exit*: required accuracy for solution obtained
                                    !!    *  1 *: function evaluation, (f&c)
                                    !! **failure modes**:
                                    !!     * 2 *: number of equality contraints larger than n
                                    !!     * 3 *: more than 3*n iterations in lsq subproblem
                                    !!     * 4 *: inequality constraints incompatible
                                    !!     * 5 *: singular matrix e in lsq subproblem
                                    !!     * 6 *: singular matrix c in lsq subproblem
                                    !!     * 7 *: rank-deficient equality constraint subproblem hfti
                                    !!     * 8 *: positive directional derivative for linesearch
                                    !!     * 9 *: more than iter iterations in sqp
                                    !!  * >=10 *: working space w or jw too small,
                                    !!            w should be enlarged to l_w=mode/1000,
                                    !!            jw should be enlarged to l_jw=mode-1000*l_w
    integer,intent(in) :: l_w       !!   the length of w, which should be at least:
                                    !!   (3*n1+m)*(n1+1)                        *for lsq*
                                    !!  +(n1-meq+1)*(mineq+2) + 2*mineq         *for lsi*
                                    !!  +(n1+mineq)*(n1-meq) + 2*meq + n1       *for lsei*
                                    !!  + n1*n/2 + 2*m + 3*n + 3*n1 + 1         *for slsqpb*
                                    !!   with mineq = m - meq + 2*n1  &  n1 = n+1
    integer,intent(in) :: l_jw      !! the length of jw which should be at least
                                    !! `mineq = m - meq + 2*(n+1)`.
    real(wp),dimension(l_w),intent(inout) :: w  !! w() is a one dimensional working space.
                                                !! the first `m+n+n*n1/2` elements of w must not be
                                                !! changed between subsequent calls of slsqp.
                                                !! on return w(1) ... w(m) contain the multipliers
                                                !! associated with the general constraints, while
                                                !! w(m+1) ... w(m+n(n+1)/2) store the cholesky factor
                                                !! l*d*l(t) of the approximate hessian of the
                                                !! lagrangian columnwise dense as lower triangular
                                                !! unit matrix l with d in its 'diagonal' and
                                                !! w(m+n(n+1)/2+n+2 ... w(m+n(n+1)/2+n+2+m+2n)
                                                !! contain the multipliers associated with all
                                                !! all constraints of the quadratic program finding
                                                !! the search direction to the solution x*
    integer,dimension(l_jw),intent(inout) :: jw !! jw() is a one dimensional integer working space
    type(slsqpb_data),intent(inout) :: sdat  !! data for [[slsqpb]].
    type(linmin_data),intent(inout) :: ldat  !! data for [[linmin]].

    integer :: il , im , ir , is , iu , iv , iw , ix , mineq, n1

!.... note: there seems to be two slightly different specifications
!     of the appropriate length of w. are they equivalent???

     !         notice:    for proper dimensioning of w it is recommended to
     !                    copy the following statements into the head of
     !                    the calling program (and remove the comment c)
     !#######################################################################
     !     integer len_w, len_jw, m, n, n1, meq, mineq
     !     parameter (m=... , meq=... , n=...  )
     !     parameter (n1= n+1, mineq= m-meq+n1+n1)
     !     parameter (len_w=
     !    $           (3*n1+m)*(n1+1)
     !    $          +(n1-meq+1)*(mineq+2) + 2*mineq
     !    $          +(n1+mineq)*(n1-meq) + 2*meq + n1
     !    $          +(n+1)*n/2 + 2*m + 3*n + 3*n1 + 1,
     !    $           len_jw=mineq)
     !     double precision w(len_w)
     !     integer          jw(len_jw)
     !#######################################################################

!     dim(w) =         n1*(n1+1) + meq*(n1+1) + mineq*(n1+1)  for lsq
!                    +(n1-meq+1)*(mineq+2) + 2*mineq          for lsi
!                    +(n1+mineq)*(n1-meq) + 2*meq + n1        for lsei
!                    + n1*n/2 + 2*m + 3*n +3*n1 + 1           for slsqpb
!                      with mineq = m - meq + 2*n1  &  n1 = n+1

!   check length of working arrays

      n1 = n + 1
      mineq = m - meq + n1 + n1
      il = (3*n1+m)*(n1+1) + (n1-meq+1)*(mineq+2) + 2*mineq + (n1+mineq)&
           *(n1-meq) + 2*meq + n1*n/2 + 2*m + 3*n + 4*n1 + 1
      im = max(mineq,n1-meq)
      if ( l_w<il .or. l_jw<im ) then
         mode = 1000*max(10,il)
         mode = mode + max(10,im)
         return
      endif

!   prepare data for calling sqpbdy  -  initial addresses in w

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

      call slsqpb(m,meq,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,w(ir),w(il), &
                  w(ix),w(im),w(is),w(iu),w(iv),w(iw),jw,&
                  sdat%t,sdat%f0,sdat%h1,sdat%h2,sdat%h3,sdat%h4,&
                  sdat%n1,sdat%n2,sdat%n3,sdat%t0,sdat%gs,sdat%tol,sdat%line,&
                  sdat%alpha,sdat%iexact,sdat%incons,sdat%ireset,sdat%itermx,&
                  ldat)

      end subroutine slsqp
!*******************************************************************************

!*******************************************************************************
!>
!  nonlinear programming by solving sequentially quadratic programs
!
!  l1 - line search,  positive definite  bfgs update

    subroutine slsqpb(m,meq,la,n,x,xl,xu,f,c,g,a,acc,iter,mode,r,l,x0,&
                      mu,s,u,v,w,iw,&
                      t,f0,h1,h2,h3,h4,n1,n2,n3,t0,gs,tol,line,&
                      alpha,iexact,incons,ireset,itermx,ldat)
    implicit none

    real(wp) ,intent(inout) :: t
    real(wp) ,intent(inout) :: f0
    real(wp) ,intent(inout) :: h1
    real(wp) ,intent(inout) :: h2
    real(wp) ,intent(inout) :: h3
    real(wp) ,intent(inout) :: h4
    integer  ,intent(inout) :: n1
    integer  ,intent(inout) :: n2
    integer  ,intent(inout) :: n3
    real(wp) ,intent(inout) :: t0
    real(wp) ,intent(inout) :: gs
    real(wp) ,intent(inout) :: tol
    integer  ,intent(inout) :: line
    real(wp) ,intent(inout) :: alpha
    integer  ,intent(inout) :: iexact
    integer  ,intent(inout) :: incons
    integer  ,intent(inout) :: ireset
    integer  ,intent(inout) :: itermx
    type(linmin_data),intent(inout) :: ldat !! data for [[linmin]].

    integer :: iw(*), i, iter, k, j, la, m, meq, mode, n

    real(wp) :: a(la,n+1) , c(la) , g(n+1) , l((n+1)*(n+2)/2) , &
                  mu(la) , r(m+n+n+2) , s(n+1) , u(n+1) , v(n+1) , &
                  w(*) , x(n) , xl(n) , xu(n) , x0(n) , &
                  acc , f

!     dim(w) =         n1*(n1+1) + meq*(n1+1) + mineq*(n1+1)  for lsq
!                     +(n1-meq+1)*(mineq+2) + 2*mineq
!                     +(n1+mineq)*(n1-meq) + 2*meq + n1       for lsei
!                      with mineq = m - meq + 2*n1  &  n1 = n+1

      if ( mode<0 ) then

          !   call jacobian at current x

          !   update cholesky-factors of hessian matrix by modified bfgs formula

         do i = 1 , n
            u(i) = g(i) - ddot(m,a(1,i),1,r,1) - v(i)
         enddo

         !   l'*s

         k = 0
         do i = 1 , n
            h1 = zero
            k = k + 1
            do j = i + 1 , n
               k = k + 1
               h1 = h1 + l(k)*s(j)
            enddo
            v(i) = s(i) + h1
         enddo

         !   d*l'*s

         k = 1
         do i = 1 , n
            v(i) = l(k)*v(i)
            k = k + n1 - i
         enddo

         !   l*d*l'*s

         do i = n , 1 , -1
            h1 = zero
            k = i
            do j = 1 , i - 1
               h1 = h1 + l(k)*v(j)
               k = k + n - j
            enddo
            v(i) = v(i) + h1
         enddo

         h1 = ddot(n,s,1,u,1)
         h2 = ddot(n,s,1,v,1)
         h3 = 0.2_wp*h2
         if ( h1<h3 ) then
            h4 = (h2-h3)/(h2-h1)
            h1 = h3
            call dscal(n,h4,u,1)
            call daxpy(n,one-h4,v,1,u,1)
         endif
         call ldl(n,l,u,+one/h1,v)
         call ldl(n,l,v,-one/h2,u)

         !   end of main iteration

         goto 200

      elseif ( mode==0 ) then

         itermx = iter
         if ( acc>=zero ) then
            iexact = 0
         else
            iexact = 1
         endif
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

      else

         !   call functions at current x

         t = f
         do j = 1 , m
            if ( j<=meq ) then
               h1 = c(j)
            else
               h1 = zero
            endif
            t = t + mu(j)*max(-c(j),h1)
         enddo
         h1 = t - t0
         if ( iexact+1==1 ) then
            if ( h1<=h3/ten .or. line>10 ) goto 500
            alpha = max(h3/(two*(h3-h1)),alfmin)
            goto 300
         elseif ( iexact+1==2 ) then
            goto 400
         else
            goto 500
         endif

      endif

!   reset bfgs matrix

 100  ireset = ireset + 1
      if ( ireset>5 ) then
         ! check relaxed convergence in case of positive directional derivative
         if ( (abs(f-f0)<tol .or. dnrm2(n,s,1)<tol) .and. h3<tol ) then
            mode = 0
         else
            mode = 8
         endif
         return
      else
         l(1) = zero
         call dcopy(n2,l(1),0,l,1)
         j = 1
         do i = 1 , n
            l(j) = one
            j = j + n1 - i
         enddo
      endif

!   main iteration : search direction, steplength, ldl'-update

 200  iter = iter + 1
      mode = 9
      if ( iter>itermx ) return

!   search direction as solution of qp - subproblem

      call dcopy(n,xl,1,u,1)
      call dcopy(n,xu,1,v,1)
      call daxpy(n,-one,x,1,u,1)
      call daxpy(n,-one,x,1,v,1)
      h4 = one
      call lsq(m,meq,n,n3,la,l,g,a,c,u,v,s,r,w,iw,mode)

!   augmented problem for inconsistent linearization

      if ( mode==6 ) then
         if ( n==meq ) mode = 4
      endif
      if ( mode==4 ) then
         do j = 1 , m
            if ( j<=meq ) then
               a(j,n1) = -c(j)
            else
               a(j,n1) = max(-c(j),zero)
            endif
         enddo
         s(1) = zero
         call dcopy(n,s(1),0,s,1)
         h3 = zero
         g(n1) = zero
         l(n3) = hun
         s(n1) = one
         u(n1) = zero
         v(n1) = one
         incons = 0
 250     call lsq(m,meq,n1,n3,la,l,g,a,c,u,v,s,r,w,iw,mode)
         h4 = one - s(n1)
         if ( mode==4 ) then
            l(n3) = ten*l(n3)
            incons = incons + 1
            if ( incons<=5 ) goto 250
            return
         elseif ( mode/=1 ) then
            return
         endif
      elseif ( mode/=1 ) then
         return
      endif

!   update multipliers for l1-test

      do i = 1 , n
         v(i) = g(i) - ddot(m,a(1,i),1,r,1)
      enddo
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
         endif
         h2 = h2 + max(-c(j),h3)
         h3 = abs(r(j))
         mu(j) = max(h3,(mu(j)+h3)/two)
         h1 = h1 + h3*abs(c(j))
      enddo

!   check convergence

      mode = 0
      if ( h1<acc .and. h2<acc ) return
      h1 = zero
      do j = 1 , m
         if ( j<=meq ) then
            h3 = c(j)
         else
            h3 = zero
         endif
         h1 = h1 + mu(j)*max(-c(j),h3)
      enddo
      t0 = f + h1
      h3 = gs - h1*h4
      mode = 8
      if ( h3>=zero ) goto 100

!   line search with an l1-testfunction

      line = 0
      alpha = one
      if ( iexact==1 ) goto 400

!   inexact linesearch

 300  line = line + 1
      h3 = alpha*h3
      call dscal(n,alpha,s,1)
      call dcopy(n,x0,1,x,1)
      call daxpy(n,one,s,1,x,1)

      call enforce_bounds(x,xl,xu)  ! ensure that x doesn't violate bounds

      mode = 1
      return

!   exact linesearch

 400  if ( line/=3 ) then
         alpha = linmin(line,alfmin,one,t,tol, &
                         ldat%a, ldat%b, ldat%d, ldat%e, ldat%p,   ldat%q,   &
                         ldat%r, ldat%u, ldat%v, ldat%w, ldat%x,   ldat%m,   &
                         ldat%fu,ldat%fv,ldat%fw,ldat%fx,ldat%tol1,ldat%tol2 )
         call dcopy(n,x0,1,x,1)
         call daxpy(n,alpha,s,1,x,1)
         mode = 1
         return
      endif
      call dscal(n,alpha,s,1)

!   check convergence

 500  h3 = zero
      do j = 1 , m
         if ( j<=meq ) then
            h1 = c(j)
         else
            h1 = zero
         endif
         h3 = h3 + max(-c(j),h1)
      enddo
      if ( (abs(f-f0)<acc .or. dnrm2(n,s,1)<acc) .and. h3<acc ) then
         mode = 0
      else
         mode = -1
      endif

      end subroutine slsqpb
!*******************************************************************************

!*******************************************************************************
!>
!   minimize with respect to x
!
!             ||e*x - f||
!                                      1/2  t
!   with upper triangular matrix e = +d   *l ,
!
!                                      -1/2  -1
!                     and vector f = -d    *l  *g,
!
!  where the unit lower tridiangular matrix l is stored columnwise
!  dense in the n*(n+1)/2 array l with vector d stored in its
! 'diagonal' thus substituting the one-elements of l
!
!   subject to
!
!             a(j)*x - b(j) = 0 ,         j=1,...,meq,
!             a(j)*x - b(j) >=0,          j=meq+1,...,m,
!             xl(i) <= x(i) <= xu(i),     i=1,...,n,
!     on entry, the user has to provide the arrays l, g, a, b, xl, xu.
!     with dimensions: l(n*(n+1)/2), g(n), a(la,n), b(m), xl(n), xu(n)
!     the working array w must have at least the following dimension:
!     dim(w) =        (3*n+m)*(n+1)                        for lsq
!                    +(n-meq+1)*(mineq+2) + 2*mineq        for lsi
!                    +(n+mineq)*(n-meq) + 2*meq + n        for lsei
!                      with mineq = m - meq + 2*n
!     on return, no array will be changed by the subroutine.
!     x     stores the n-dimensional solution vector
!     y     stores the vector of lagrange multipliers of dimension
!           m+n+n (constraints+lower+upper bounds)
!     mode  is a success-failure flag with the following meanings:
!          mode=1: successful computation
!               2: error return because of wrong dimensions (n<1)
!               3: iteration count exceeded by nnls
!               4: inequality constraints incompatible
!               5: matrix e is not of full rank
!               6: matrix c is not of full rank
!               7: rank defect in hfti
!
!     coded            dieter kraft, april 1987
!     revised                        march 1989

    subroutine lsq(m,meq,n,nl,la,l,g,a,b,xl,xu,x,y,w,jw,mode)
    implicit none

      real(wp) :: l , g , a , b , w , xl , xu , x , y , diag , xnorm

      integer :: jw(*) , i , ic , id , ie , if , ig , ih , il , im , ip , &
                 iu , iw , i1 , i2 , i3 , i4 , la , m , meq , mineq , &
                 mode , m1 , n , nl , n1 , n2 , n3

      dimension a(la,n) , b(la) , g(n) , l(nl) , w(*) , x(n) , xl(n) , &
                xu(n) , y(m+n+n)

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
      endif
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
      enddo
      if ( n2==1 ) then
         w(i3) = l(nl)
         w(i4) = zero
         call dcopy(n3,w(i4),0,w(i4),1)
         w(if-1+n) = zero
      endif
      call dscal(n,-one,w(if),1)

      ic = if + n
      id = ic + meq*n

      if ( meq>0 ) then

         !  recover matrix c from upper part of a

         do i = 1 , meq
            call dcopy(n,a(i,1),la,w(ic-1+i),meq)
         enddo

         !  recover vector d from upper part of b

         call dcopy(meq,b(1),1,w(id),1)
         call dscal(meq,-one,w(id),1)

      endif

      ig = id + meq

      if ( mineq>0 ) then
         !  recover matrix g from lower part of a
         do i = 1 , mineq
            call dcopy(n,a(meq+i,1),la,w(ig-1+i),m1)
         enddo
      endif

      !  augment matrix g by +i and -i

      ip = ig + mineq
      do i = 1 , n
         w(ip-1+i) = zero
         call dcopy(n,w(ip-1+i),0,w(ip-1+i),m1)
      enddo
      w(ip) = one
      call dcopy(n,w(ip),0,w(ip),m1+1)

      im = ip + n
      do i = 1 , n
         w(im-1+i) = zero
         call dcopy(n,w(im-1+i),0,w(im-1+i),m1)
      enddo
      w(im) = -one
      call dcopy(n,w(im),0,w(im),m1+1)

      ih = ig + m1*n

      if ( mineq>0 ) then
         ! recover h from lower part of b
         call dcopy(mineq,b(meq+1),1,w(ih),1)
         call dscal(mineq,-one,w(ih),1)
      endif

      !  augment vector h by xl and xu

      il = ih + mineq
      call dcopy(n,xl,1,w(il),1)
      iu = il + n
      call dcopy(n,xu,1,w(iu),1)
      call dscal(n,-one,w(iu),1)

      iw = iu + n

      call lsei(w(ic),w(id),w(ie),w(if),w(ig),w(ih),max(1,meq),meq,n,n, &
                m1,m1,n,x,xnorm,w(iw),jw,mode)

      if ( mode==1 ) then
         ! restore lagrange multipliers
         call dcopy(m,w(iw),1,y(1),1)
         call dcopy(n3,w(iw+m),1,y(m+1),1)
         call dcopy(n3,w(iw+m+n),1,y(m+n3+1),1)
         call enforce_bounds(x,xl,xu)  ! to ensure that bounds are not violated
      endif

      end subroutine lsq
!*******************************************************************************

!*******************************************************************************
!>
!  for mode=1, the subroutine returns the solution x of
!  equality & inequality constrained least squares problem lsei :
!
!                min ||e*x - f||
!                 x
!
!                s.t.  c*x  = d,
!                      g*x >= h.
!
!     using qr decomposition & orthogonal basis of nullspace of c
!     chapter 23.6 of lawson & hanson: solving least squares problems.
!
!     the following dimensions of the arrays defining the problem
!     are necessary
!     dim(e) :   formal (le,n),    actual (me,n)
!     dim(f) :   formal (le  ),    actual (me  )
!     dim(c) :   formal (lc,n),    actual (mc,n)
!     dim(d) :   formal (lc  ),    actual (mc  )
!     dim(g) :   formal (lg,n),    actual (mg,n)
!     dim(h) :   formal (lg  ),    actual (mg  )
!     dim(x) :   formal (n   ),    actual (n   )
!     dim(w) :   2*mc+me+(me+mg)*(n-mc)  for lsei
!              +(n-mc+1)*(mg+2)+2*mg     for lsi
!     dim(jw):   max(mg,l)
!     on entry, the user has to provide the arrays c, d, e, f, g, and h.
!     on return, all arrays will be changed by the subroutine.
!     x     stores the solution vector
!     xnorm stores the residuum of the solution in euclidian norm
!     w     stores the vector of lagrange multipliers in its first
!           mc+mg elements
!     mode  is a success-failure flag with the following meanings:
!          mode=1: successful computation
!               2: error return because of wrong dimensions (n<1)
!               3: iteration count exceeded by nnls
!               4: inequality constraints incompatible
!               5: matrix e is not of full rank
!               6: matrix c is not of full rank
!               7: rank defect in hfti
!
!     18.5.1981, dieter kraft, dfvlr oberpfaffenhofen
!     20.3.1987, dieter kraft, dfvlr oberpfaffenhofen

    subroutine lsei(c,d,e,f,g,h,lc,mc,le,me,lg,mg,n,x,xnrm,w,jw,mode)
    implicit none

      integer :: jw(*) , i , ie , if , ig , iw , j , k , krank , l , lc , &
                 le , lg , mc , mc1 , me , mg , mode , n
      real(wp) :: c(lc,n) , e(le,n) , g(lg,n) , d(lc) , f(le) , &
                  h(lg) , x(n) , w(*) , t , xnrm , dum(1)

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
         enddo

!  solve c*x=d and modify f

         mode = 6
         do i = 1 , mc
            if ( abs(c(i,i))<epmach ) return
            x(i) = (d(i)-ddot(i-1,c(i,1),lc,x,1))/c(i,i)
         enddo
         mode = 1
         w(mc1) = zero
         !call dcopy(mg-mc,w(mc1),0,w(mc1),1)  ! original code
         call dcopy(mg,w(mc1),0,w(mc1),1)      ! bug fix for when meq = n

         if ( mc/=n ) then

            do i = 1 , me
               w(if-1+i) = f(i) - ddot(mc,e(i,1),le,x,1)
            enddo

!  store transformed e & g

            do i = 1 , me
               call dcopy(l,e(i,mc1),le,w(ie-1+i),me)
            enddo
            do i = 1 , mg
               call dcopy(l,g(i,mc1),lg,w(ig-1+i),mg)
            enddo

            if ( mg>0 ) then
!  modify h and solve inequality constrained ls problem

               do i = 1 , mg
                  h(i) = h(i) - ddot(mc,g(i,1),lg,x,1)
               enddo
               call lsi(w(ie),w(if),w(ig),h,me,me,mg,mg,l,x(mc1),xnrm,  &
                        w(mc1),jw,mode)
               if ( mc==0 ) return
               t = dnrm2(mc,x,1)
               xnrm = sqrt(xnrm*xnrm+t*t)
               if ( mode/=1 ) return
            else

!  solve ls without inequality constraints

               mode = 7
               k = max(le,n)
               t = sqrt(epmach)
               call hfti(w(ie),me,me,l,w(if),k,1,t,krank,dum,w,w(l+1),jw)
               xnrm = dum(1)
               call dcopy(l,w(if),1,x(mc1),1)
               if ( krank/=l ) return
               mode = 1
            endif
         endif

!  solution of original problem and lagrange multipliers

         do i = 1 , me
            f(i) = ddot(n,e(i,1),le,x,1) - f(i)
         enddo
         do i = 1 , mc
            d(i) = ddot(me,e(1,i),1,f,1) &
                   - ddot(mg,g(1,i),1,w(mc1),1)
         enddo

         do i = mc , 1 , -1
            call h12(2,i,i+1,n,c(i,1),lc,w(iw+i),x,1,1,1)
         enddo

         do i = mc , 1 , -1
            j = min(i+1,lc)
            w(i) = (d(i)-ddot(mc-i,c(j,i),1,w(j),1))/c(i,i)
         enddo
      endif

      end subroutine lsei
!*******************************************************************************

!*******************************************************************************
!>
!  for mode=1, the subroutine returns the solution x of
!  inequality constrained linear least squares problem:
!
!                    min ||e*x-f||
!                     x
!
!                    s.t.  g*x >= h
!
!     the algorithm is based on qr decomposition as described in
!     chapter 23.5 of lawson & hanson: solving least squares problems
!
!     the following dimensions of the arrays defining the problem
!     are necessary
!     dim(e) :   formal (le,n),    actual (me,n)
!     dim(f) :   formal (le  ),    actual (me  )
!     dim(g) :   formal (lg,n),    actual (mg,n)
!     dim(h) :   formal (lg  ),    actual (mg  )
!     dim(x) :   n
!     dim(w) :   (n+1)*(mg+2) + 2*mg
!     dim(jw):   lg
!     on entry, the user has to provide the arrays e, f, g, and h.
!     on return, all arrays will be changed by the subroutine.
!     x     stores the solution vector
!     xnorm stores the residuum of the solution in euclidian norm
!     w     stores the vector of lagrange multipliers in its first
!           mg elements
!     mode  is a success-failure flag with the following meanings:
!          mode=1: successful computation
!               2: error return because of wrong dimensions (n<1)
!               3: iteration count exceeded by nnls
!               4: inequality constraints incompatible
!               5: matrix e is not of full rank
!
!     03.01.1980, dieter kraft: coded
!     20.03.1987, dieter kraft: revised to fortran 77

    subroutine lsi(e,f,g,h,le,me,lg,mg,n,x,xnorm,w,jw,mode)
    implicit none

      integer :: i , j , le , lg , me , mg , mode , n , jw(lg)
      real(wp) :: e(le,n) , f(le) , g(lg,n) , h(lg) , x(n) , w(*) , &
                  xnorm , t

!  qr-factors of e and application to f

      do i = 1 , n
         j = min(i+1,n)
         call h12(1,i,i+1,me,e(1,i),1,t,e(1,j),1,le,n-i)
         call h12(2,i,i+1,me,e(1,i),1,t,f,1,1,1)
      enddo

!  transform g and h to get least distance problem

      mode = 5
      do i = 1 , mg
         do j = 1 , n
            if ( abs(e(j,j))<epmach ) return
            g(i,j) = (g(i,j)-ddot(j-1,g(i,1),lg,e(1,j),1))/e(j,j)
         enddo
         h(i) = h(i) - ddot(n,g(i,1),lg,f,1)
      enddo

!  solve least distance problem

      call ldp(g,lg,mg,n,h,x,xnorm,w,jw,mode)
      if ( mode==1 ) then

!  solution of original problem

         call daxpy(n,one,f,1,x,1)
         do i = n , 1 , -1
            j = min(i+1,n)
            x(i) = (x(i)-ddot(n-i,e(i,j),le,x(j),1))/e(i,i)
         enddo
         j = min(n+1,me)
         t = dnrm2(me-n,f(j),1)
         xnorm = sqrt(xnorm*xnorm+t*t)
      endif

      end subroutine lsi
!*******************************************************************************

!*******************************************************************************
!>
!
!                     t
!     minimize   1/2 x x    subject to   g * x >= h.
!
!       c.l. lawson, r.j. hanson: 'solving least squares problems'
!       prentice hall, englewood cliffs, new jersey, 1974.
!
!     parameter description:
!
!     g(),mg,m,n   on entry g() stores the m by n matrix of
!                  linear inequality constraints. g() has first
!                  dimensioning parameter mg
!     h()          on entry h() stores the m vector h representing
!                  the right side of the inequality system
!
!     remark: g(),h() will not be changed during calculations by ldp
!
!     x()          on entry x() need not be initialized.
!                  on exit x() stores the solution vector x if mode=1.
!     xnorm        on exit xnorm stores the euclidian norm of the
!                  solution vector if computation is successful
!     w()          w is a one dimensional working space, the length
!                  of which should be at least (m+2)*(n+1) + 2*m
!                  on exit w() stores the lagrange multipliers
!                  associated with the constraints
!                  at the solution of problem ldp
!     index()      index() is a one dimensional integer working space
!                  of length at least m
!     mode         mode is a success-failure flag with the following
!                  meanings:
!          mode=1: successful computation
!               2: error return because of wrong dimensions (n<=0)
!               3: iteration count exceeded by nnls
!               4: inequality constraints incompatible

    subroutine ldp(g,mg,m,n,h,x,xnorm,w,index,mode)
    implicit none

      real(wp) :: g , h , x , xnorm , w , u , v , fac , rnorm
      integer :: index , i , if , iw , iwdual , iy , iz , j , m , mg , &
                 mode , n , n1

      dimension g(mg,n) , h(m) , x(n) , w(*) , index(m)

      mode = 2
      if ( n>0 ) then

         ! state dual problem

         mode = 1
         x(1) = zero
         call dcopy(n,x(1),0,x,1)
         xnorm = zero
         if ( m/=0 ) then
            iw = 0
            do j = 1 , m
               do i = 1 , n
                  iw = iw + 1
                  w(iw) = g(j,i)
               enddo
               iw = iw + 1
               w(iw) = h(j)
            enddo
            if = iw + 1
            do i = 1 , n
               iw = iw + 1
               w(iw) = zero
            enddo
            w(iw+1) = one
            n1 = n + 1
            iz = iw + 2
            iy = iz + n1
            iwdual = iy + m

            ! solve dual problem

            call nnls(w,n1,n1,m,w(if),w(iy),rnorm,w(iwdual),w(iz),index,mode)

            if ( mode==1 ) then
               mode = 4
               if ( rnorm>zero ) then

                  ! compute solution of primal problem

                  fac = one - ddot(m,h,1,w(iy),1)
                  if ( diff(one+fac,one)>zero ) then
                     mode = 1
                     fac = one/fac
                     do j = 1 , n
                        x(j) = fac*ddot(m,g(1,j),1,w(iy),1)
                     enddo
                     xnorm = dnrm2(n,x,1)

                     ! compute lagrange multipliers for primal problem

                     w(1) = zero
                     call dcopy(m,w(1),0,w,1)
                     call daxpy(m,fac,w(iy),1,w,1)
                  endif
               endif
            endif
         endif
      endif

      end subroutine ldp
!*******************************************************************************

!*******************************************************************************
    pure elemental function diff(u,v) result(d)
        !! replaced statement function in original code
        implicit none
        real(wp),intent(in) :: u
        real(wp),intent(in) :: v
        real(wp) :: d
        d = u - v
    end function diff
!*******************************************************************************

!*******************************************************************************
!>
!  c.l.lawson and r.j.hanson, jet propulsion laboratory:
!  'solving least squares problems'. prentice-hall.1974
!
!      **********   nonnegative least squares   **********
!
!     given an m by n matrix, a, and an m-vector, b, compute an
!     n-vector, x, which solves the least squares problem
!
!                  a*x = b  subject to  x >= 0
!
!     a(),mda,m,n
!            mda is the first dimensioning parameter for the array,a().
!            on entry a()  contains the m by n matrix,a.
!            on exit a() contains the product q*a,
!            where q is an m by m orthogonal matrix generated
!            implicitly by this subroutine.
!            either m>=n or m<n is permissible.
!            there is no restriction on the rank of a.
!     b()    on entry b() contains the m-vector, b.
!            on exit b() contains q*b.
!     x()    on entry x() need not be initialized.
!            on exit x() will contain the solution vector.
!     rnorm  on exit rnorm contains the euclidean norm of the
!            residual vector.
!     w()    an n-array of working space.
!            on exit w() will contain the dual solution vector.
!            w will satisfy w(i)=0 for all i in set p
!            and w(i)<=0 for all i in set z
!     z()    an m-array of working space.
!     index()an integer working array of length at least n.
!            on exit the contents of this array define the sets
!            p and z as follows:
!            index(1)    thru index(nsetp) = set p.
!            index(iz1)  thru index (iz2)  = set z.
!            iz1=nsetp + 1 = npp1, iz2=n.
!     mode   this is a success-failure flag with the following meaning:
!            1    the solution has been computed successfully.
!            2    the dimensions of the problem are wrong,
!                 either m <= 0 or n <= 0.
!            3    iteration count exceeded, more than 3*n iterations.
!
!     revised          dieter kraft, march 1983

    subroutine nnls(a,mda,m,n,b,x,rnorm,w,z,index,mode)
    implicit none

      integer :: i , ii , ip , iter , itmax , iz , izmax , iz1 , iz2 , j , &
                 jj , jz , k , l , m , mda , mode , n , npp1 , nsetp , &
                 index(n)

      real(wp) :: a(mda,n) , b(m) , x(n) , w(n) , z(m) , asave , &
                  wmax , alpha , c , s , t , u , v , up , rnorm , unorm

      real(wp),parameter :: factor = 1.0e-2_wp

      mode = 2
      if ( m>0 .and. n>0 ) then
         mode = 1
         iter = 0
         itmax = 3*n

! step one (initialize)

         do i = 1 , n
            index(i) = i
         enddo
         iz1 = 1
         iz2 = n
         nsetp = 0
         npp1 = 1
         x(1) = zero
         call dcopy(n,x(1),0,x,1)

! step two (compute dual variables)
! .....entry loop a

 50      if ( iz1<=iz2 .and. nsetp<m ) then
            do iz = iz1 , iz2
               j = index(iz)
               w(j) = ddot(m-nsetp,a(npp1,j),1,b(npp1),1)
            enddo

! step three (test dual variables)

 60         wmax = zero
            do iz = iz1 , iz2
               j = index(iz)
               if ( w(j)>wmax ) then
                  wmax = w(j)
                  izmax = iz
               endif
            enddo

! .....exit loop a

            if ( wmax>zero ) then
               iz = izmax
               j = index(iz)

! step four (test index j for linear dependency)

               asave = a(npp1,j)
               call h12(1,npp1,npp1+1,m,a(1,j),1,up,z,1,1,0)
               unorm = dnrm2(nsetp,a(1,j),1)
               t = factor*abs(a(npp1,j))
               if ( diff(unorm+t,unorm)>zero ) then
                  call dcopy(m,b,1,z,1)
                  call h12(2,npp1,npp1+1,m,a(1,j),1,up,z,1,1,1)
                  if ( z(npp1)/a(npp1,j)>zero ) then
! step five (add column)

                     call dcopy(m,z,1,b,1)
                     index(iz) = index(iz1)
                     index(iz1) = j
                     iz1 = iz1 + 1
                     nsetp = npp1
                     npp1 = npp1 + 1
                     do jz = iz1 , iz2
                        jj = index(jz)
                        call h12(2,nsetp,npp1,m,a(1,j),1,up,a(1,jj),1,mda,1)
                     enddo
                     k = min(npp1,mda)
                     w(j) = zero
                     call dcopy(m-nsetp,w(j),0,a(k,j),1)

! step six (solve least squares sub-problem)
! .....entry loop b

 62                  do ip = nsetp , 1 , -1
                        if ( ip/=nsetp ) call daxpy(ip,-z(ip+1),a(1,jj),1,z,1)
                        jj = index(ip)
                        z(ip) = z(ip)/a(ip,jj)
                     enddo
                     iter = iter + 1
                     if ( iter<=itmax ) then
! step seven to ten (step length algorithm)

                        alpha = one
                        jj = 0
                        do ip = 1 , nsetp
                           if ( z(ip)<=zero ) then
                              l = index(ip)
                              t = -x(l)/(z(ip)-x(l))
                              if ( alpha>=t ) then
                                 alpha = t
                                 jj = ip
                              endif
                           endif
                        enddo
                        do ip = 1 , nsetp
                           l = index(ip)
                           x(l) = (one-alpha)*x(l) + alpha*z(ip)
                        enddo

! .....exit loop b

                        if ( jj==0 ) goto 50

! step eleven (delete column)

                        i = index(jj)
 64                     x(i) = zero
                        jj = jj + 1
                        do j = jj , nsetp
                           ii = index(j)
                           index(j-1) = ii
                           call dsrotg(a(j-1,ii),a(j,ii),c,s)
                           t = a(j-1,ii)
                           call dsrot(n,a(j-1,1),mda,a(j,1),mda,c,s)
                           a(j-1,ii) = t
                           a(j,ii) = zero
                           call dsrot(1,b(j-1),1,b(j),1,c,s)
                        enddo
                        npp1 = nsetp
                        nsetp = nsetp - 1
                        iz1 = iz1 - 1
                        index(iz1) = i
                        if ( nsetp<=0 ) then
                           mode = 3
                           goto 100
                        else
                           do jj = 1 , nsetp
                              i = index(jj)
                              if ( x(i)<=zero ) goto 64
                           enddo
                           call dcopy(m,b,1,z,1)
                           goto 62
                        endif
                     else
                        mode = 3
                        goto 100
                     endif
                  endif
               endif
               a(npp1,j) = asave
               w(j) = zero
               goto 60
            endif
         endif
! step twelve (solution)

 100     k = min(npp1,m)
         rnorm = dnrm2(m-nsetp,b(k),1)
         if ( npp1>m ) then
            w(1) = zero
            call dcopy(n,w(1),0,w,1)
         endif
      endif

      end subroutine nnls
!*******************************************************************************

!*******************************************************************************
!>
!     rank-deficient least squares algorithm as described in:
!     c.l.lawson and r.j.hanson, jet propulsion laboratory, 1973 jun 12
!     to appear in 'solving least squares problems', prentice-hall, 1974
!
!     a(*,*),mda,m,n   the array a initially contains the m x n matrix a
!                      of the least squares problem ax = b.
!                      the first dimensioning parameter mda must satisfy
!                      mda >= m. either m >= n or m < n is permitted.
!                      there is no restriction on the rank of a.
!                      the matrix a will be modified by the subroutine.
!     b(*,*),mdb,nb    if nb = 0 the subroutine will make no reference
!                      to the array b. if nb > 0 the array b() must
!                      initially contain the m x nb matrix b  of the
!                      the least squares problem ax = b and on return
!                      the array b() will contain the n x nb solution x.
!                      if nb>1 the array b() must be double subscripted
!                      with first dimensioning parameter mdb>=max(m,n),
!                      if nb=1 the array b() may be either single or
!                      double subscripted.
!     tau              absolute tolerance parameter for pseudorank
!                      determination, provided by the user.
!     krank            pseudorank of a, set by the subroutine.
!     rnorm            on exit, rnorm(j) will contain the euclidian
!                      norm of the residual vector for the problem
!                      defined by the j-th column vector of the array b.
!     h(), g()         arrays of working space of length >= n.
!     ip()             integer array of working space of length >= n
!                      recording permutation indices of column vectors

    subroutine hfti(a,mda,m,n,b,mdb,nb,tau,krank,rnorm,h,g,ip)
    implicit none

    integer :: i , j , jb , k , kp1 , krank , l , ldiag , lmax , m , &
               mda , mdb , n , nb , ip(n)
    real(wp) :: a(mda,n) , b(mdb,nb) , h(n) , g(n) , rnorm(nb) , &
                tau , hmax , tmp , &
                u , v

    real(wp),parameter :: factor = 1.0e-3_wp

    k = 0
    ldiag = min(m,n)
    if ( ldiag<=0 ) then
       krank = k
       return
    else

       ! compute lmax

       do j = 1 , ldiag
          if ( j/=1 ) then
             lmax = j
             do l = j , n
                h(l) = h(l) - a(j-1,l)**2
                if ( h(l)>h(lmax) ) lmax = l
             enddo
             if ( diff(hmax+factor*h(lmax),hmax)>zero ) goto 20
          endif
          lmax = j
          do l = j , n
             h(l) = zero
             do i = j , m
                h(l) = h(l) + a(i,l)**2
             enddo
             if ( h(l)>h(lmax) ) lmax = l
          enddo
          hmax = h(lmax)

          ! column interchanges if needed

20        ip(j) = lmax
          if ( ip(j)/=j ) then
             do i = 1 , m
                tmp = a(i,j)
                a(i,j) = a(i,lmax)
                a(i,lmax) = tmp
             enddo
             h(lmax) = h(j)
          endif

          ! j-th transformation and application to a and b

          i = min(j+1,n)
          call h12(1,j,j+1,m,a(1,j),1,h(j),a(1,i),1,mda,n-j)
          call h12(2,j,j+1,m,a(1,j),1,h(j),b,1,mdb,nb)
       enddo

       !determine pseudorank:

       do j=1,ldiag
          if (abs(a(j,j))<=tau) exit
       end do
       k=j-1
       kp1=j

    end if

    ! norm of residuals

    do jb = 1 , nb
       rnorm(jb) = dnrm2(m-k,b(kp1,jb),1)
    enddo
    if ( k>0 ) then
       if ( k/=n ) then
          ! householder decomposition of first k rows
          do i = k , 1 , -1
             call h12(1,i,kp1,n,a(i,1),mda,g(i),a,mda,1,i-1)
          enddo
       endif
       do jb = 1 , nb

          ! solve k*k triangular system

          do i = k , 1 , -1
             j = min(i+1,n)
             b(i,jb) = (b(i,jb)-ddot(k-i,a(i,j),mda,b(j,jb),1))/a(i,i)
          enddo

          ! complete solution vector

          if ( k/=n ) then
             do j = kp1 , n
                b(j,jb) = zero
             enddo
             do i = 1 , k
                call h12(2,i,kp1,n,a(i,1),mda,g(i),b(1,jb),1,mdb,1)
             enddo
          endif

          ! reorder solution according to previous column interchanges

          do j = ldiag , 1 , -1
             if ( ip(j)/=j ) then
                l = ip(j)
                tmp = b(l,jb)
                b(l,jb) = b(j,jb)
                b(j,jb) = tmp
             endif
          enddo
       enddo
    else
       do jb = 1 , nb
          do i = 1 , n
             b(i,jb) = zero
          enddo
       enddo
    endif
    krank = k
    end subroutine hfti
!*******************************************************************************

!*******************************************************************************
!>
!     c.l.lawson and r.j.hanson, jet propulsion laboratory, 1973 jun 12
!     to appear in 'solving least squares problems', prentice-hall, 1974
!
!     construction and/or application of a single
!     householder transformation  q = i + u*(u**t)/b
!
!     mode    = 1 or 2   to select algorithm  h1  or  h2 .
!     lpivot is the index of the pivot element.
!     l1,m   if l1 <= m   the transformation will be constructed to
!            zero elements indexed from l1 through m.
!            if l1 > m the subroutine does an identity transformation.
!     u(),iue,up
!            on entry to h1 u() stores the pivot vector.
!            iue is the storage increment between elements.
!            on exit from h1 u() and up store quantities defining
!            the vector u of the householder transformation.
!            on entry to h2 u() and up
!            should store quantities previously computed by h1.
!            these will not be modified by h2.
!     c()    on entry to h1 or h2 c() stores a matrix which will be
!            regarded as a set of vectors to which the householder
!            transformation is to be applied.
!            on exit c() stores the set of transformed vectors.
!     ice    storage increment between elements of vectors in c().
!     icv    storage increment between vectors in c().
!     ncv    number of vectors in c() to be transformed.
!            if ncv <= 0 no operations will be done on c().

    subroutine h12(mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)
    implicit none

      integer :: incr , ice , icv , iue , lpivot , l1 , mode , ncv
      integer :: i , i2 , i3 , i4 , j , m
      real(wp) :: u , up , c , cl , clinv , b , sm

      dimension u(iue,*) , c(*)

      if ( 0<lpivot .and. lpivot<l1 .and. l1<=m ) then
         cl = abs(u(1,lpivot))
         if ( mode/=2 ) then

             ! ****** construct the transformation ******

            do j = l1 , m
               sm = abs(u(1,j))
               cl = max(sm,cl)
            enddo
            if ( cl<=zero ) return
            clinv = one/cl
            sm = (u(1,lpivot)*clinv)**2
            do j = l1 , m
               sm = sm + (u(1,j)*clinv)**2
            enddo
            cl = cl*sqrt(sm)
            if ( u(1,lpivot)>zero ) cl = -cl
            up = u(1,lpivot) - cl
            u(1,lpivot) = cl

            ! ****** apply the transformation  i+u*(u**t)/b  to c ******

         elseif ( cl<=zero ) then
            return
         endif
         if ( ncv>0 ) then
            b = up*u(1,lpivot)
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
                  enddo
                  if ( sm/=zero ) then
                     sm = sm*b
                     c(i2) = c(i2) + sm*up
                     do i = l1 , m
                        c(i4) = c(i4) + sm*u(1,i)
                        i4 = i4 + ice
                     enddo
                  endif
               enddo
            endif
         endif
      endif

      end subroutine h12
!*******************************************************************************

!*******************************************************************************
!>
!   ldl     ldl' - rank-one - update
!
!   purpose:
!           updates the ldl' factors of matrix a by rank-one matrix
!           sigma*z*z'
!
!   input arguments: (* means parameters are changed during execution)
!     n     : order of the coefficient matrix a
!   * a     : positive definite matrix of dimension n;
!             only the lower triangle is used and is stored column by
!             column as one dimensional array of dimension n*(n+1)/2.
!   * z     : vector of dimension n of updating elements
!     sigma : scalar factor by which the modifying dyade z*z' is
!             multiplied
!
!   output arguments:
!     a     : updated ldl' factors
!
!   working array:
!     w     : vector op dimension n (used only if sigma < zero)
!
!   method:
!     that of fletcher and powell as described in :
!     fletcher,r.,(1974) on the modification of ldl' factorization.
!     powell,m.j.d.      math.computation 28, 1067-1078.
!
!   implemented by:
!     kraft,d., dfvlr - institut fuer dynamik der flugsysteme
!               d-8031  oberpfaffenhofen
!
!   status: 15. january 1980

    subroutine ldl(n,a,z,sigma,w)
    implicit none

      integer :: i , ij , j , n
      real(wp) :: a(*) , t , v , w(*) , z(*) , u , tp , &
                       beta , alpha , delta , gamma , sigma

      if ( sigma/=zero ) then
         ij = 1
         t = one/sigma
         if ( sigma<=zero ) then
            ! prepare negative update
            do i = 1 , n
               w(i) = z(i)
            enddo
            do i = 1 , n
               v = w(i)
               t = t + v*v/a(ij)
               do j = i + 1 , n
                  ij = ij + 1
                  w(j) = w(j) - v*a(ij)
               enddo
               ij = ij + 1
            enddo
            if ( t>=zero ) t = epmach/sigma
            do i = 1 , n
               j = n + 1 - i
               ij = ij - i
               u = w(j)
               w(j) = t
               t = t - u*u/a(ij)
            enddo
         endif
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
               enddo
            else
               do j = i + 1 , n
                  ij = ij + 1
                  z(j) = z(j) - v*a(ij)
                  a(ij) = a(ij) + beta*z(j)
               enddo
            endif
            ij = ij + 1
            t = tp
         enddo
      endif

      end subroutine ldl
!*******************************************************************************

!*******************************************************************************
!>
!   linmin  linesearch without derivatives
!   (used if exact = 1)
!
!   purpose:
!
!  to find the argument linmin where the function f takes it's minimum
!  on the interval ax, bx.
!  combination of golden section and successive quadratic interpolation.
!
!   input arguments: (* means parameters are changed during execution)
!
! *mode   see output arguments
!  ax     left endpoint of initial interval
!  bx     right endpoint of initial interval
!  f      function value at linmin which is to be brought in by
!         reverse communication controlled by mode
!  tol    desired length of interval of uncertainty of final result
!
!   output arguments:
!
!  linmin abscissa approximating the point where f attains a minimum
!  mode   controls reverse communication
!         must be set to 0 initially, returns with intermediate
!         values 1 and 2 which must not be changed by the user,
!         ends with convergence with value 3.
!
!   method:
!
!  this function subprogram is a slightly modified version of the
!  algol 60 procedure localmin given in
!  r.p. brent: algorithms for minimization without derivatives,
!              prentice-hall (1973).
!
!   implemented by:
!
!     kraft, d., dfvlr - institut fuer dynamik der flugsysteme
!                d-8031  oberpfaffenhofen
!
!   status: 31. august  1984

    real(wp) function linmin(mode,ax,bx,f,tol,&
                                a,b,d,e,p,q,r,u,v,&
                                w,x,m,fu,fv,fw,fx,tol1,tol2)

    implicit none

      integer :: mode
      real(wp) :: f,tol,ax,bx
      real(wp),intent(inout) :: a,b,d,e,p,q,r,u,v,w,x,m,fu,fv,fw,fx,tol1,tol2

      real(wp),parameter :: c    = (3.0_wp-sqrt(5.0_wp))/2.0_wp  !! golden section ratio = `0.381966011`
      real(wp),parameter :: eps  = sqrt(epsilon(1.0_wp))         !! square - root of machine precision

      if ( mode==1 ) then

          !  main loop starts here

         fx = f
         fv = fx
         fw = fv

      elseif ( mode==2 ) then

         fu = f
         !  update a, b, v, w, and x
         if ( fu>fx ) then
            if ( u<x ) a = u
            if ( u>=x ) b = u
            if ( fu<=fw .or. w==x ) then
               v = w
               fv = fw
               w = u
               fw = fu
            elseif ( fu<=fv .or. v==x .or. v==w ) then
               v = u
               fv = fu
            endif
         else
            if ( u>=x ) a = x
            if ( u<x ) b = x
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
         endif
      else
         !  initialization
         a = ax
         b = bx
         e = zero
         v = a + c*(b-a)
         w = v
         x = w
         linmin = x
         mode = 1
         return
      endif
      m = 0.5_wp*(a+b)
      tol1 = eps*abs(x) + tol
      tol2 = tol1 + tol1

      !  test convergence

      if ( abs(x-m)<=tol2-0.5_wp*(b-a) ) then
         !  end of main loop
         linmin = x
         mode = 3
      else
         r = zero
         q = r
         p = q
         if ( abs(e)>tol1 ) then
            !  fit parabola
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q - (x-w)*r
            q = q - r
            q = q + q
            if ( q>zero ) p = -p
            if ( q<zero ) q = -q
            r = e
            e = d
         endif

         !  is parabola acceptable
         if ( abs(p)>=0.5_wp*abs(q*r) .or. p<=q*(a-x) .or. p>=q*(b-x) ) then
            !  golden section step
            if ( x>=m ) e = a - x
            if ( x<m ) e = b - x
            d = c*e
         else
            !  parabolic interpolation step
            d = p/q
            !  f must not be evaluated too close to a or b
            if ( u-a<tol2 ) d = sign(tol1,m-x)
            if ( b-u<tol2 ) d = sign(tol1,m-x)
         endif

         !  f must not be evaluated too close to x
         if ( abs(d)<tol1 ) d = sign(tol1,d)
         u = x + d
         linmin = u
         mode = 2

      endif

      end function linmin
!*******************************************************************************

!*******************************************************************************
!>
!  enforce the bound constraints on x.

    subroutine enforce_bounds(x,xl,xu)

    implicit none

    real(wp),dimension(:),intent(inout) :: x   !! optimization variable vector
    real(wp),dimension(:),intent(in)    :: xl  !! lower bounds (must be same dimension as `x`)
    real(wp),dimension(:),intent(in)    :: xu  !! upper bounds (must be same dimension as `x`)

    where (x<xl)
        x = xl
    elsewhere (x>xu)
        x = xu
    end where

    end subroutine enforce_bounds
!*******************************************************************************

!*******************************************************************************
    end module slsqp_module
!*******************************************************************************

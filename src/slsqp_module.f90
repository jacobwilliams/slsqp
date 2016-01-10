!*******************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!  Module containing the object-oriented interface to the SLSQP method.
!  It is called using the [[slsqp_solver]] class, which
!  is the only public entity in this module.

    module slsqp_module

    use slsqp_kinds
    use slsqp_support
    use slsqp_core
    use iso_fortran_env, only: error_unit

    implicit none

    private

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
        write(error_unit,*) 'error: invalid upper or lower bound vector size'
    else if (meq<0 .or. meq>m) then
        write(error_unit,*) 'error: invalid meq value:', meq
    else if (m<0) then
        write(error_unit,*) 'error: invalid m value:', m
    else if (n<1) then
        write(error_unit,*) 'error: invalid n value:', n
    else if (any(xl>xu)) then
        write(error_unit,*) 'error: lower bounds must be <= upper bounds.'
    else

        !two linesearch modes:
        select case (linesearch_mode)
        case(1)     !inexact
            me%linesearch_mode = linesearch_mode
        case(2)     !exact
            me%linesearch_mode = linesearch_mode
            allocate(me%linmin)
        case default
            write(error_unit,*) 'error: invalid linesearch_mode (must be 1 or 2): ',linesearch_mode
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
    if (size(x)/=me%n) then
        write(error_unit,*) 'invalid size(x) in slsqp_wrapper'
        istat = -100
        return
    end if

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
        write(error_unit,*) 'invalid linesearch_mode in slsqp_wrapper'
        istat = -101
        return
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
            write(*,*) 'more than max_iter iterations in slsqp'
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
    end module slsqp_module
!*******************************************************************************

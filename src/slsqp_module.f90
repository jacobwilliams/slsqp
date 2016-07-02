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
    use iso_fortran_env, only: error_unit,output_unit

    implicit none

    private

    type,public :: slsqp_solver

        !! The main class used to interface with the SLSQP solver.

        private

        integer  :: n        = 0        !! number of optimization variables (\( n > 0 \))
        integer  :: m        = 0        !! number of constraints (\( m \ge 0 \))
        integer  :: meq      = 0        !! number of equality constraints (\( m \ge m_{eq} \ge 0 \))
        integer  :: max_iter = 0        !! maximum number of iterations
        real(wp) :: acc      = zero     !! accuracy tolerance

        !these two were not in the original code:
        real(wp) :: alphamin = 0.1_wp   !! min \( \alpha \) for line search \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)
        real(wp) :: alphamax = 1.0_wp   !! max \( \alpha \) for line search \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)

        integer :: iprint = output_unit !! unit number of status printing (0 for no printing)

        real(wp),dimension(:),allocatable :: xl  !! lower bound on x
        real(wp),dimension(:),allocatable :: xu  !! upper bound on x

        integer :: l_w = 0 !! size of `w`
        real(wp),dimension(:),allocatable :: w !! real work array

        integer :: l_jw = 0 !! size of `jw`
        integer,dimension(:),allocatable :: jw  !! integer work array

        procedure(func),pointer     :: f      => null()         !! problem function subroutine
        procedure(grad),pointer     :: g      => null()         !! gradient subroutine
        procedure(iterfunc),pointer :: report => null()         !! for reporting an iteration

        integer :: linesearch_mode = 1  !! linesearch mode.
                                        !! `1` = inexact (Armijo) linesearch,
                                        !! `2` = exact linesearch.
        type(linmin_data) :: linmin !! data formerly within [[linmin]].
                                                !! Only used when `linesearch_mode=2`
        type(slsqpb_data) :: slsqpb  !! data formerly within [[slsqpb]].

        logical :: user_triggered_stop = .false.    !! if the `abort` method has been called
                                                    !! to stop the iterations

    contains

        private

        procedure,public :: initialize => initialize_slsqp
        procedure,public :: destroy    => destroy_slsqp
        procedure,public :: optimize   => slsqp_wrapper
        procedure,public :: abort      => stop_iterations

        procedure :: report_message  !! for reporting messages to the user

    end type slsqp_solver

    abstract interface
        subroutine func(me,x,f,c)  !! for computing the function
            import :: wp,slsqp_solver
            implicit none
            class(slsqp_solver),intent(inout) :: me
            real(wp),dimension(:),intent(in)  :: x  !! optimization variable vector
            real(wp),intent(out)              :: f  !! value of the objective function
            real(wp),dimension(:),intent(out) :: c  !! the constraint vector `dimension(m)`,
                                                    !! equality constraints (if any) first.
        end subroutine func
        subroutine grad(me,x,g,a)  !! for computing the gradients
            import :: wp,slsqp_solver
            implicit none
            class(slsqp_solver),intent(inout)   :: me
            real(wp),dimension(:),intent(in)    :: x  !! optimization variable vector
            real(wp),dimension(:),intent(out)   :: g  !! objective function partials w.r.t x `dimension(n)`
            real(wp),dimension(:,:),intent(out) :: a  !! gradient matrix of constraints w.r.t. x `dimension(m,n)`
        end subroutine grad
        subroutine iterfunc(me,iter,x,f,c)  !! for reporting an interation
            import :: wp,slsqp_solver
            implicit none
            class(slsqp_solver),intent(inout) :: me
            integer,intent(in)                :: iter  !! iteration number
            real(wp),dimension(:),intent(in)  :: x     !! optimization variable vector
            real(wp),intent(in)               :: f     !! value of the objective function
            real(wp),dimension(:),intent(in)  :: c     !! the constraint vector `dimension(m)`,
                                                       !! equality constraints (if any) first.
        end subroutine iterfunc
    end interface

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  A method that the user can call to stop the iterations.
!  (it can be called in any of the functions).
!  SLSQP will stop at the end of the next iteration.

    subroutine stop_iterations(me)

    implicit none

    class(slsqp_solver),intent(inout) :: me

    me%user_triggered_stop = .true.

    end subroutine stop_iterations
!*******************************************************************************

!*******************************************************************************
!>
!  initialize the [[slsqp_solver]] class.  see [[slsqp]] for more details.

    subroutine initialize_slsqp(me,n,m,meq,max_iter,acc,f,g,xl,xu,status_ok,&
                                linesearch_mode,iprint,report,alphamin,alphamax)

    implicit none

    class(slsqp_solver),intent(inout) :: me
    integer,intent(in)                :: n               !! the number of varibles, \( n \ge 1 \)
    integer,intent(in)                :: m               !! total number of constraints, \( m \ge 0 \)
    integer,intent(in)                :: meq             !! number of equality constraints, \( m_{eq} \ge 0 \)
    integer,intent(in)                :: max_iter        !! maximum number of iterations
    procedure(func)                   :: f               !! problem function
    procedure(grad)                   :: g               !! function to compute gradients
    real(wp),dimension(n),intent(in)  :: xl              !! lower bounds on `x`
    real(wp),dimension(n),intent(in)  :: xu              !! upper bounds on `x`
    real(wp),intent(in)               :: acc             !! accuracy
    logical,intent(out)               :: status_ok       !! will be false if there were errors
    integer,intent(in),optional       :: linesearch_mode !! 1 = inexact (default), 2 = exact
    integer,intent(in),optional       :: iprint          !! unit number of status messages (default=output_unit)
    procedure(iterfunc),optional      :: report          !! user-defined procedure that will be called once per iteration
    real(wp),intent(in),optional      :: alphamin        !! minimum alpha for linesearch [default 0.1]
    real(wp),intent(in),optional      :: alphamax        !! maximum alpha for linesearch [default 1.0]

    integer :: n1,mineq,i

    status_ok = .false.
    call me%destroy()

    if (present(iprint)) me%iprint = iprint

    if (size(xl)/=size(xu) .or. size(xl)/=n) then
        call me%report_message('error: invalid upper or lower bound vector size')
        call me%report_message('  size(xl) =',ival=size(xl))
        call me%report_message('  size(xu) =',ival=size(xu))
        call me%report_message('  n        =',ival=n)
    else if (meq<0 .or. meq>m) then
        call me%report_message('error: invalid meq value:', ival=meq)
    else if (m<0) then
        call me%report_message('error: invalid m value:', ival=m)
    else if (n<1) then
        call me%report_message('error: invalid n value:', ival=n)
    else if (any(xl>xu)) then
        call me%report_message('error: lower bounds must be <= upper bounds.')
        do i=1,n
            if (xl(i)>xu(i)) then
                call me%report_message('  xl(i)>xu(i) for variable',ival=i)
            end if
        end do
    else

        if (present(linesearch_mode)) then
            !two linesearch modes:
            select case (linesearch_mode)
            case(1)     !inexact
                me%linesearch_mode = linesearch_mode
            case(2)     !exact
                me%linesearch_mode = linesearch_mode
            case default
                call me%report_message('error: invalid linesearch_mode (must be 1 or 2): ',&
                                        ival=linesearch_mode)
                call me%destroy()
                return
            end select
        end if

        !optional linesearch bounds:
        if (present(alphamin)) me%alphamin = alphamin
        if (present(alphamax)) me%alphamax = alphamax

        !verify valid values for alphamin and alphamax: 0<alphamin<alphamax<=1
        if (me%alphamin<=0.0_wp .or. me%alphamax<=0.0_wp .or. &
            me%alphamax<=me%alphamin .or. &
            me%alphamin>=1.0_wp .or. me%alphamax>1.0_wp) then

            call me%report_message('error: invalid values for alphamin or alphamax.')
            call me%report_message('  alphamin =',rval=me%alphamin)
            call me%report_message('  alphamax =',rval=me%alphamax)
            call me%destroy()
            return

        end if

        status_ok = .true.
        me%n = n
        me%m = m
        me%meq = meq
        me%max_iter = max_iter
        me%acc = acc
        me%f => f
        me%g => g
        if (present(report)) me%report => report

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
        me%w = zero

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

    subroutine slsqp_wrapper(me,x,istat,iterations,status_message)

    implicit none

    class(slsqp_solver),intent(inout)   :: me
    real(wp),dimension(:),intent(inout) :: x          !! **in:**  initial optimization variables,
                                                      !! **out:** solution.
    integer,intent(out)                 :: istat      !! status code (see `mode` in [[slsqp]]).
    integer,intent(out),optional        :: iterations !! number of iterations
    character(len=:),intent(out),allocatable,optional :: status_message !! string status message corresponding to `istat`

    !local variables:
    real(wp)                               :: f        !! objective function
    real(wp),dimension(max(1,me%m))        :: c        !! constraint vector
    real(wp),dimension(max(1,me%m),me%n+1) :: a        !! a matrix for slsqp
    real(wp),dimension(me%n+1)             :: g        !! g matrix for slsqp
    real(wp),dimension(me%m)               :: cvec     !! constraint vector
    real(wp),dimension(me%n)               :: dfdx     !! objective function partials
    real(wp),dimension(me%m,me%n)          :: dcdx     !! constraint partials
    integer                                :: i        !! iteration counter
    integer                                :: mode     !! reverse communication flag for [[slsqp]]
    integer                                :: la       !! input to [[slsqp]]
    integer                                :: iter     !! in/out for [[slsqp]]
    real(wp)                               :: acc      !! in/out for [[slsqp]]

    !initialize:
    i    = 0
    iter = me%max_iter
    la   = max(1,me%m)
    mode = 0
    a    = zero
    g    = zero
    c    = zero
    if (present(iterations)) iterations = 0
    call me%linmin%destroy()
    call me%slsqpb%destroy()

    !check setup:
    if (size(x)/=me%n) then
        istat = -100
        call me%report_message(mode_to_status_message(istat))
        if (present(status_message)) status_message = mode_to_status_message(istat)
        return
    end if

    !linesearch:
    select case(me%linesearch_mode)
    case(1)     !inexact (armijo-type linesearch)
        acc = abs(me%acc)
    case(2)     !exact
        acc = -abs(me%acc)
    case default
        istat = -101
        call me%report_message(mode_to_status_message(istat))
        if (present(status_message)) status_message = mode_to_status_message(istat)
        return
    end select

    !make sure the functions have been associated:
    if (.not. associated(me%f)) then
        istat = -102
        call me%report_message(mode_to_status_message(istat))
        if (present(status_message)) status_message = mode_to_status_message(istat)
        return
    end if
    if (.not. associated(me%g)) then
        istat = -103
        call me%report_message(mode_to_status_message(istat))
        if (present(status_message)) status_message = mode_to_status_message(istat)
        return
    end if

    !main solver loop:
    do

        if (mode==0 .or. mode==1) then  !function evaluation (f&c)
            call me%f(x,f,cvec)
            c(1:me%m)   = cvec
        end if

        if (mode==0 .or. mode==-1) then  !gradient evaluation (g&a)
            call me%g(x,dfdx,dcdx)
            g(1:me%n)        = dfdx
            a(1:me%m,1:me%n) = dcdx

            !this is an iteration:
            !note: the initial guess is reported as iteration 0:
            if (associated(me%report)) call me%report(i,x,f,c) !report iteration
            i = i + 1  ! iteration counter
        end if

        !main routine:
        call slsqp(me%m,me%meq,la,me%n,x,me%xl,me%xu,&
                    f,c,g,a,acc,iter,mode,&
                    me%w,me%l_w,me%jw,me%l_jw,&
                    me%slsqpb,me%linmin,me%alphamin,me%alphamax)

        if (mode==1 .or. mode==-1) then
            !continue to next call
        else
            if (mode==0 .and. associated(me%report)) &
                call me%report(i,x,f,c) !report solution
            call me%report_message(mode_to_status_message(mode))
            exit
        end if

        if (me%user_triggered_stop) then
            mode = -2
            call me%report_message(mode_to_status_message(mode))
            me%user_triggered_stop = .false.    !have to reset in case
                                                !method is called again.
            exit
        end if

    end do

    istat = mode
    if (present(iterations)) iterations = iter
    if (present(status_message)) status_message = mode_to_status_message(istat)

    end subroutine slsqp_wrapper
!*******************************************************************************

!*******************************************************************************
!>
!  Report a message from an [[slsqp_solver]] class. This uses the `iprint`
!  variable in the class as the unit number for printing. Note: for fatal errors,
!  if no unit is specified, the `error_unit` is used.

    subroutine report_message(me,str,ival,rval,fatal)

    implicit none

    class(slsqp_solver),intent(in) :: me
    character(len=*),intent(in)    :: str    !! the message to report.
    integer,intent(in),optional    :: ival   !! optional integer to print after the message.
    real(wp),intent(in),optional   :: rval   !! optional real to print after the message.
    logical,intent(in),optional    :: fatal  !! if True, then the program is stopped (default=False).

    logical :: stop_program  !! true if the program is to be stopped
    logical :: write_message  !! true if the message is to be printed
    character(len=10) :: istr  !! string version of `ival`
    character(len=30) :: rstr  !! string version of `rval`
    character(len=:),allocatable :: str_to_write  !! the actual message to the printed
    integer :: istat  !! iostat for integer to string conversion

    !fatal error check:
    if (present(fatal)) then
        stop_program = fatal
    else
        stop_program = .false.
    end if

    !note: if stopping program, then the message is always printed:
    write_message = me%iprint/=0 .or. stop_program

    if (write_message) then

        if (present(ival)) then
            write(istr,fmt='(I10)',iostat=istat) ival
            if (istat/=0) istr = '*****'
            str_to_write = str//' '//trim(adjustl(istr))
        elseif (present(rval)) then
            write(istr,fmt='(F30.16)',iostat=istat) rval
            if (istat/=0) rstr = '*****'
            str_to_write = str//' '//trim(adjustl(rstr))
        else
            str_to_write = str
        end if

        if (me%iprint==0) then
            write(error_unit,'(A)') str_to_write  !in this case, use the error unit
        else
            write(me%iprint,'(A)') str_to_write   !user specified unit number
        end if
        deallocate(str_to_write)

        if (stop_program) error stop 'Fatal Error'

    end if

    end subroutine report_message
!*******************************************************************************

!*******************************************************************************
!>
!  Convert the [[slsqp]] `mode` flag to a message string.

    pure function mode_to_status_message(imode) result(message)

    implicit none

    integer,intent(in) :: imode
    character(len=:),allocatable :: message

    select case (imode)
    case(0) !required accuracy for solution obtained
        message = 'required accuracy for solution obtained'
    case(-100)
        message = 'invalid size(x) in slsqp_wrapper'
    case(-101)
        message = 'invalid linesearch_mode in slsqp_wrapper'
    case(-102)
        message = 'function is not associated'
    case(-103)
        message = 'gradient function is not associated'
    case(-2)
        message = 'user-triggered stop of slsqp'
    case(1,-1)
        message = 'in progress'
    case(2)
        message = 'number of equality contraints larger than n'
    case(3)
        message = 'more than 3*n iterations in lsq subproblem'
    case(4)
        message = 'inequality constraints incompatible'
    case(5)
        message = 'singular matrix e in lsq subproblem'
    case(6)
        message = 'singular matrix c in lsq subproblem'
    case(7)
        message = 'rank-deficient equality constraint subproblem hfti'
    case(8)
        message = 'positive directional derivative for linesearch'
    case(9)
        message = 'more than max_iter iterations in slsqp'
    case default
        message = 'unknown slsqp error'
    end select

    end function mode_to_status_message
!*******************************************************************************

!*******************************************************************************
    end module slsqp_module
!*******************************************************************************

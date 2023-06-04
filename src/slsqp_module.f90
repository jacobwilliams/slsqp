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
    use, intrinsic :: iso_fortran_env, only: error_unit,output_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

    implicit none

    private

    type,public :: slsqp_solver

        !! The main class used to interface with the SLSQP solver.

        private

        integer  :: n           = 0        !! number of optimization variables (\( n > 0 \))
        integer  :: m           = 0        !! number of constraints (\( m \ge 0 \))
        integer  :: meq         = 0        !! number of equality constraints (\( m \ge m_{eq} \ge 0 \))
        integer  :: max_iter    = 0        !! maximum number of iterations
        real(wp) :: acc         = zero     !! accuracy tolerance
        real(wp) :: tolf        = -one     !! accuracy tolerance over f:  if \( |f| < tolf \) then stop
        real(wp) :: toldf       = -one     !! accuracy tolerance over df: if \( |f_{n+1} - f_n| < toldf \) then stop.
                                           !! It's different from `acc` in the case of positive derivative
        real(wp) :: toldx       = -one     !! accuracy tolerance over dx: if \( |x_{n+1} - x_n| < toldx \) then stop

        integer  :: gradient_mode = 0      !! how the gradients are computed:
                                           !!
                                           !! * 0 - use the user-supplied `g` subroutine. [default]
                                           !! * 1 - approximate by basic backward differences
                                           !! * 2 - approximate by basic forward differences
                                           !! * 3 - approximate by basic central differences
        real(wp) :: gradient_delta  = 1.0e8_wp !! perturbation step size to approximate gradients
                                               !! by finite differences (`gradient_mode` 1-3).

        !these two were not in the original code:
        real(wp) :: alphamin = 0.1_wp   !! min \( \alpha \) for line search \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)
        real(wp) :: alphamax = 1.0_wp   !! max \( \alpha \) for line search \( 0 < \alpha_{min} < \alpha_{max} \le 1 \)

        integer :: iprint = output_unit !! unit number of status printing (0 for no printing)

        real(wp),dimension(:),allocatable :: xl  !! lower bound on x
        real(wp),dimension(:),allocatable :: xu  !! upper bound on x

        integer :: l_w = 0 !! size of `w`
        real(wp),dimension(:),allocatable :: w !! real work array

        procedure(func),pointer     :: f      => null()         !! problem function subroutine
        procedure(grad),pointer     :: g      => null()         !! gradient subroutine
        procedure(iterfunc),pointer :: report => null()         !! for reporting an iteration

        integer :: linesearch_mode = 1  !! linesearch mode:
                                        !!
                                        !! * `1` = inexact (Armijo) linesearch,
                                        !! * `2` = exact linesearch.
        type(linmin_data) :: linmin !! data formerly within [[linmin]].
                                    !! Only used when `linesearch_mode=2`
        type(slsqpb_data) :: slsqpb  !! data formerly within [[slsqpb]].

        ! note: the following two maybe should be combined into a separate type
        ! along with the two methods...
        integer :: nnls_mode = 1 !! Which NNLS method to use:
                                 !!
                                 !! 1. Use the original [[nnls]]
                                 !! 2. Use the newer [[bvls]]
        integer :: max_iter_ls = 0  !! max iterations in the least squares problem.
                                    !! if `<=0`, defaults to `3*n`.
                                    !! (use by either [[nnls]] or [[bvls]])

        logical :: user_triggered_stop = .false.    !! if the `abort` method has been called
                                                    !! to stop the iterations

        real(wp) :: infinite_bound = huge(one)  !! "infinity" for the upper and lower bounds.
                                                !! if `xl<=-infinite_bound` or `xu>=infinite_bound`
                                                !! then these bounds are considered nonexistant.

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
        subroutine iterfunc(me,iter,x,f,c)  !! for reporting an iteration
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
                                linesearch_mode,iprint,report,alphamin,alphamax,&
                                gradient_mode,gradient_delta,tolf,toldf,toldx,&
                                max_iter_ls,nnls_mode,infinite_bound)

    implicit none

    class(slsqp_solver),intent(inout) :: me
    integer,intent(in)                :: n               !! the number of variables, \( n \ge 1 \)
    integer,intent(in)                :: m               !! total number of constraints, \( m \ge 0 \)
    integer,intent(in)                :: meq             !! number of equality constraints, \( m_{eq} \ge 0 \)
    integer,intent(in)                :: max_iter        !! maximum number of iterations
    procedure(func)                   :: f               !! problem function
    procedure(grad)                   :: g               !! function to compute gradients (must be
                                                         !! associated if `gradient_mode=0`)
    real(wp),dimension(n),intent(in)  :: xl              !! lower bounds on `x`.
                                                         !! `xl(i)=NaN` (or `xl(i)<=-infinite_bound`) indicates to ignore `i`th bound
    real(wp),dimension(n),intent(in)  :: xu              !! upper bounds on `x`.
                                                         !! `xu(i)=NaN` (or `xu(i)>=infinite_bound`) indicates to ignore `i`th bound
    real(wp),intent(in)               :: acc             !! accuracy
    logical,intent(out)               :: status_ok       !! will be false if there were errors
    integer,intent(in),optional       :: linesearch_mode !! 1 = inexact (default), 2 = exact
    integer,intent(in),optional       :: iprint          !! unit number of status messages (default=`output_unit`)
    procedure(iterfunc),optional      :: report          !! user-defined procedure that will be called once per iteration
    real(wp),intent(in),optional      :: alphamin        !! minimum alpha for linesearch [default 0.1]
    real(wp),intent(in),optional      :: alphamax        !! maximum alpha for linesearch [default 1.0]
    integer,intent(in),optional       :: gradient_mode   !! how the gradients are to be computed:
                                                         !!
                                                         !! * 0 - use the user-supplied `g` subroutine. [default]
                                                         !! * 1 - approximate by basic backward differences
                                                         !! * 2 - approximate by basic forward differences
                                                         !! * 3 - approximate by basic central differences
                                                         !!
                                                         !! Note that modes 1-3 do not respect the variable bounds.
    real(wp),intent(in),optional      :: gradient_delta  !! perturbation step size (>epsilon) to compute the approximated
                                                         !! gradient by finite differences (`gradient_mode` 1-3).
                                                         !! note that this is an absolute step that does not respect
                                                         !! the `xl` or `xu` variable bounds.
    real(wp),intent(in),optional      :: tolf            !! stopping criterion if \( |f| < tolf \) then stop.
    real(wp),intent(in),optional      :: toldf           !! stopping criterion if \( |f_{n+1} - f_n| < toldf \) then stop
    real(wp),intent(in),optional      :: toldx           !! stopping criterion if \( ||x_{n+1} - x_n|| < toldx \) then stop
    integer,intent(in),optional       :: max_iter_ls     !! maximum number of iterations in the [[nnls]] problem
    integer,intent(in),optional       :: nnls_mode       !! Which NNLS method to use:
                                                         !!
                                                         !! 1. Use the original [[nnls]]
                                                         !! 2. Use the newer [[bvls]]
    real(wp),intent(in),optional :: infinite_bound !! "infinity" for the upper and lower bounds.
                                                   !! if `xl<=-infinite_bound` or `xu>=infinite_bound`
                                                   !! then these bounds are considered nonexistant.
                                                   !! If not present then `huge()` is used for this.

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
    else if (any(xl>xu .and. .not. ieee_is_nan(xl) .and. .not. ieee_is_nan(xu))) then
        call me%report_message('error: lower bounds must be <= upper bounds.')
        do i=1,n
            if (xl(i)>xu(i) .and. .not. ieee_is_nan(xl(i)) .and. .not. ieee_is_nan(xu(i))) then
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
        if (me%alphamin<=zero .or. me%alphamax<=zero .or. &
            me%alphamax<=me%alphamin .or. &
            me%alphamin>=one .or. me%alphamax>one) then

            call me%report_message('error: invalid values for alphamin or alphamax.')
            call me%report_message('  alphamin =',rval=me%alphamin)
            call me%report_message('  alphamax =',rval=me%alphamax)
            call me%destroy()
            return

        end if

        if (present(tolf))  me%tolf  = tolf
        if (present(toldf)) me%toldf = toldf
        if (present(toldx)) me%toldx = toldx

        if (present(max_iter_ls)) me%max_iter_ls = max_iter_ls
        if (present(nnls_mode)) then
            select case (nnls_mode)
            case(1:2)
                me%nnls_mode = nnls_mode
            case default
                call me%report_message('error: invalid value for nnls_mode. defaulting to 1.')
                me%nnls_mode = 1
            end select
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

        if (present(gradient_mode)) then
            me%gradient_mode = gradient_mode
            if (present(gradient_delta)) then
                me%gradient_delta = gradient_delta
            end if
        end if

        if (present(infinite_bound)) then
            me%infinite_bound = abs(infinite_bound)
        else
            me%infinite_bound = huge(one)
        end if

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
    character(len=:),intent(out),allocatable,optional :: status_message
                                                      !! string status message
                                                      !! corresponding to `istat`

    ! local variables:
    real(wp),dimension(:),allocatable   :: c        !! constraint vector -- `dimension(max(1,me%m))`
    real(wp),dimension(:,:),allocatable :: a        !! a matrix for [[slsqp]] -- `dimension(max(1,me%m),me%n+1)`
    real(wp),dimension(:),allocatable   :: g        !! g matrix for [[slsqp]] -- `dimension(me%n+1)`
    real(wp),dimension(:),allocatable   :: cvec     !! constraint vector -- `dimension(me%m)`
    real(wp),dimension(:),allocatable   :: dfdx     !! objective function partials -- `dimension(me%n)`
    real(wp),dimension(:,:),allocatable :: dcdx     !! constraint partials -- `dimension(me%m,me%n)`
    real(wp),dimension(:),allocatable   :: delta    !! perturbation step size to approximate gradient -- `dimension(me%n)`
    real(wp),dimension(:),allocatable   :: cvecr    !! right function value to approximate constraints vector's gradient -- `dimension(me%m)`
    real(wp),dimension(:),allocatable   :: cvecl    !! left function value to approximate constraints vector's gradient -- `dimension(me%m)`
    real(wp)                            :: f        !! objective function
    integer                             :: i        !! iteration counter
    integer                             :: mode     !! reverse communication flag for [[slsqp]]
    integer                             :: la       !! input to [[slsqp]]
    integer                             :: iter     !! in/out for [[slsqp]]
    real(wp)                            :: acc      !! in/out for [[slsqp]]
    integer                             :: ig       !! loop index to approximate gradient
    real(wp)                            :: fr       !! right function value to approximate objective function's gradient
    real(wp)                            :: fl       !! left function value to approximate objective function's gradient
    real(wp)                            :: fact     !! denominator factor for finite difference approximation

    !initialize:
    allocate(c(max(1,me%m))       )
    allocate(a(max(1,me%m),me%n+1))
    allocate(g(me%n+1)            )
    allocate(cvec(me%m)           )
    allocate(dfdx(me%n)           )
    allocate(dcdx(me%m,me%n)      )
    allocate(delta(me%n)          )
    allocate(cvecr(me%m)          )
    allocate(cvecl(me%m)          )
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
    if ((me%gradient_mode==0).and.(.not. associated(me%g))) then
        istat = -103
        call me%report_message(mode_to_status_message(istat))
        if (present(status_message)) status_message = mode_to_status_message(istat)
        return
    end if
    if (me%gradient_mode<0 .or. me%gradient_mode>3) then
        istat = -104
        call me%report_message(mode_to_status_message(istat))
        if (present(status_message)) status_message = mode_to_status_message(istat)
        return
    end if
    if (me%gradient_mode/=0 .and. me%gradient_delta<=epmach) then
        istat = -105
        call me%report_message(mode_to_status_message(istat))
        if (present(status_message)) status_message = mode_to_status_message(istat)
        return
    end if

    !main solver loop:
    do

        if (mode==0 .or. mode==1) then  !function evaluation (f&c)
            call me%f(x,f,cvec)
            c(1:me%m) = cvec
        end if

        if (mode==0 .or. mode==-1) then  !gradient evaluation (g&a)
            select case (me%gradient_mode)
            case (0) ! user supplied gradients
                call me%g(x,dfdx,dcdx)
                g(1:me%n)        = dfdx
                a(1:me%m,1:me%n) = dcdx
            case default ! approximate using finite differences
                if (me%gradient_mode==3) then
                    fact = two  ! central differences
                else
                    fact = one  ! forward/backward differences
                end if
                do ig=1,me%n
                    !initialize a delta to perturb the objective
                    !function and the constraint vector
                    delta     = zero
                    delta(ig) = me%gradient_delta
                    !get the right and left value of the objective
                    !function and the constraint vector
                    select case (me%gradient_mode)
                    case (1) ! backward difference
                        call me%f(x,fr,cvecr)
                        call me%f(x-delta,fl,cvecl)
                    case (2) ! forward difference
                        call me%f(x+delta,fr,cvecr)
                        call me%f(x,fl,cvecl)
                    case (3) ! central difference
                        call me%f(x+delta,fr,cvecr)
                        call me%f(x-delta,fl,cvecl)
                    end select
                    !compute the gradients by first-order finite differences
                    g(ig) = (fr-fl) / ( fact*delta(ig) )
                    if (me%m>0) then
                        a(:,ig) = (cvecr-cvecl) / ( fact*delta(ig) )
                    end if
                end do
            end select
            !this is an iteration:
            !note: the initial guess is reported as iteration 0:
            if (associated(me%report)) call me%report(i,x,f,c) !report iteration
            i = i + 1  ! iteration counter
        end if

        !main routine:
        call slsqp(me%m,me%meq,la,me%n,x,me%xl,me%xu,&
                   f,c,g,a,acc,iter,mode,&
                   me%w,me%l_w, &
                   me%slsqpb,me%linmin,me%alphamin,me%alphamax,&
                   me%tolf,me%toldf,me%toldx,&
                   me%max_iter_ls,me%nnls_mode,&
                   me%infinite_bound)

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
        message = 'Required accuracy for solution obtained'
    case(-100)
        message = 'Invalid size(x) in slsqp_wrapper'
    case(-101)
        message = 'Invalid linesearch_mode in slsqp_wrapper'
    case(-102)
        message = 'Function is not associated'
    case(-103)
        message = 'Gradient function is not associated'
    case(-104)
        message = 'Invalid gradient mode'
    case(-105)
        message = 'Invalid perturbation step size for finite difference gradients'
    case(-2)
        message = 'User-triggered stop of slsqp'
    case(1,-1)
        message = 'In progress'
    case(2)
        message = 'Number of equality constraints larger than n'
    case(3)
        message = 'More than 3*n iterations in lsq subproblem'
    case(4)
        message = 'Inequality constraints incompatible'
    case(5)
        message = 'Singular matrix e in lsq subproblem'
    case(6)
        message = 'Singular matrix c in lsq subproblem'
    case(7)
        message = 'Rank-deficient equality constraint subproblem hfti'
    case(8)
        message = 'Positive directional derivative for linesearch'
    case(9)
        message = 'More than max_iter iterations in slsqp'
    case default
        message = 'Unknown slsqp error'
    end select

    end function mode_to_status_message
!*******************************************************************************

!*******************************************************************************
    end module slsqp_module
!*******************************************************************************

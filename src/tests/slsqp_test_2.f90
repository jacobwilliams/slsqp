!*******************************************************************************
!> author: Jacob Williams
!
!  Test for the [[slsqp_module]].

    program slsqp_test_2

    use slsqp_module
    use slsqp_kinds

    implicit none

    integer,parameter               :: n = 3                    !! number of optimization variables
    integer,parameter               :: m = 2                    !! total number of constraints
    integer,parameter               :: meq = 1                  !! number of equality constraints
    integer,parameter               :: max_iter = 100           !! maximum number of allowed iterations
    real(wp),dimension(n),parameter :: xl = [-10.0_wp, -10.0_wp, -10.0_wp]  !! lower bounds
    real(wp),dimension(n),parameter :: xu = [ 10.0_wp,  10.0_wp,  10.0_wp]  !! upper bounds
    real(wp),parameter              :: acc = 1.0e-7_wp          !! tolerance
    integer,parameter               :: linesearch_mode = 1      !! use inexact linesearch.

    type(slsqp_solver)    :: solver      !! instantiate an slsqp solver
    real(wp),dimension(n) :: x           !! optimization variable vector
    integer               :: istat       !! for solver status check
    logical               :: status_ok   !! for initialization status check
    integer               :: iterations  !! number of iterations by the solver

    x = [1.0_wp, 2.0_wp, 3.0_wp] ! initial guess

    call solver%initialize(n,m,meq,max_iter,acc,test_func,test_grad,&
                            xl,xu,linesearch_mode=linesearch_mode,status_ok=status_ok,&
                            report=report_iteration,&
                            alphamin=0.1_wp, alphamax=0.5_wp) !to limit search steps

    if (status_ok) then
        call solver%optimize(x,istat,iterations)
        write(*,*) ''
        write(*,*) 'solution   :', x
        write(*,*) 'istat      :', istat
        write(*,*) 'iterations :', iterations
        write(*,*) ''
    else
        error stop 'error calling slsqp.'
    end if

    ! Solution is: x = [1,1,1], f = 3

    contains

    subroutine test_func(me,x,f,c)

        !!  Compute the objective function and constraints
        !!
        !!  Minimize:
        !!
        !!   * \( f = x_1^2 + x_2^2 + x_3 \)
        !!
        !!  Subject to:
        !!
        !!   * \( c_1 = x_1 x_2 - x_3 = 0 \)
        !!   * \( c_2 = x_3 - 1 \ge 0 \)

        implicit none

        class(slsqp_solver),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x   !! optimization variable vector
        real(wp),intent(out)              :: f   !! value of the objective function
        real(wp),dimension(:),intent(out) :: c   !! the constraint vector `dimension(m)`,
                                                 !! equality constraints (if any) first.

        f = x(1)**2 + x(2)**2 + x(3)  !objective function

        c(1) = x(1)*x(2) - x(3)       !equality constraint (==0)
        c(2) = x(3) - 1.0_wp          !inequality constraint (>=0)

    end subroutine test_func

    subroutine test_grad(me,x,g,a)

        !! compute the gradients.

        implicit none

        class(slsqp_solver),intent(inout)   :: me
        real(wp),dimension(:),intent(in)    :: x    !! optimization variable vector
        real(wp),dimension(:),intent(out)   :: g    !! objective function partials w.r.t x `dimension(n)`
        real(wp),dimension(:,:),intent(out) :: a    !! gradient matrix of constraints w.r.t. x `dimension(m,n)`

        g(1) = 2.0_wp*x(1)
        g(2) = 2.0_wp*x(2)
        g(3) = 1.0_wp

        a(1,1) = x(2)
        a(1,2) = x(1)
        a(1,3) = -1.0_wp

        a(2,1) = 0.0_wp
        a(2,2) = 0.0_wp
        a(2,3) = 1.0_wp

    end subroutine test_grad

    subroutine report_iteration(me,iter,x,f,c)

        !! report an iteration (print to the console).

        use, intrinsic :: iso_fortran_env, only: output_unit

        implicit none

        class(slsqp_solver),intent(inout) :: me
        integer,intent(in)                :: iter
        real(wp),dimension(:),intent(in)  :: x
        real(wp),intent(in)               :: f
        real(wp),dimension(:),intent(in)  :: c

        !write a header:
        if (iter==0) then
            write(output_unit,'(*(A20,1X))') 'iteration', &
                                             'x(1)', 'x(2)', 'x(3)', &
                                             'f(1)', 'c(1)', 'c(2)'
        end if

        !write the iteration data:
        write(output_unit,'(I20,1X,(*(F20.16,1X)))') iter,x,f,c

    end subroutine report_iteration

    end program slsqp_test_2
!*******************************************************************************

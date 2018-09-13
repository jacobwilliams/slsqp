!*******************************************************************************
!> author: Jacob Williams
!
!  Test for the [[slsqp_module]].

    program slsqp_test_3

    use slsqp_module
    use slsqp_kinds

    implicit none

    integer,parameter               :: n = 2                    !! number of optimization variables
    integer,parameter               :: m = 1                    !! total number of constraints
    integer,parameter               :: meq = 0                  !! number of equality constraints
    integer,parameter               :: max_iter = 100           !! maximum number of allowed iterations
    real(wp),dimension(n),parameter :: xl = [-1.0_wp, -1.0_wp]  !! lower bounds
    real(wp),dimension(n),parameter :: xu = [ 1.0_wp,  1.0_wp]  !! upper bounds
    real(wp),parameter              :: acc = 1.0e-8_wp          !! tolerance
    integer,parameter               :: linesearch_mode = 1      !! use inexact linesearch.
    real(wp),parameter              :: grad_delta = 1.0e-5_wp   !! step to approximate gradient

    type(slsqp_solver)      :: solver        !! instantiate an slsqp solver
    real(wp),dimension(n)   :: x             !! optimization variable vector
    integer                 :: istat         !! for solver status check
    logical                 :: status_ok     !! for initialization status check
    integer                 :: iterations    !! number of iterations by the solver
    integer                 :: gradient_mode !! gradient computation mode
    procedure(func),pointer :: f             !! pointer to `rosenbrock_func`
    procedure(grad),pointer :: g             !! not used here since we are letting
                                             !! slsqp compute the gradients

    ! test each of the gradient modes (backward, forward, and central diffs)
    do gradient_mode = 1, 3

        write(*,*) ''
        write(*,*) '---------------------------'
        write(*,*) 'gradient mode:', gradient_mode
        write(*,*) '---------------------------'

        x = [0.1_wp, 0.1_wp] !initial guess

        f => rosenbrock_func
        g => null()

        call solver%initialize(n,m,meq,max_iter,acc,f,g,&
                                xl,xu,linesearch_mode=linesearch_mode,status_ok=status_ok,&
                                report=report_iteration,&
                                gradient_mode=gradient_mode,grad_delta=grad_delta)

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

    end do

    !solution:  x1 = 0.7864151509699762
    !           x2 = 0.6176983165977831
    !           f  = 0.0456748087191604
    !           c  = -0.0000000000028654

    contains

    subroutine rosenbrock_func(me,x,f,c)
        !! Rosenbrock function
        !!
        !! Minimize the Rosenbrock function: \( f(x) = 100 (x_2 - x_1)^2 + (1 - x_1)^2 \),
        !! subject to the inequality constraint: \( x_1^2 + x_2^2 \le 1 \).
        !!
        !! see: http://www.mathworks.com/help/optim/ug/example-nonlinear-constrained-minimization.html
        implicit none
        class(slsqp_solver),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x      !! optimization variable vector
        real(wp),intent(out)              :: f      !! value of the objective function
        real(wp),dimension(:),intent(out) :: c      !! the constraint vector `dimension(m)`,
                                                    !! equality constraints (if any) first.

        f = 100.0_wp*(x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2  !objective function
        c(1) = 1.0_wp - x(1)**2 - x(2)**2  !equality constraint (>=0)

    end subroutine rosenbrock_func

    subroutine report_iteration(me,iter,x,f,c)
        use, intrinsic :: iso_fortran_env, only: output_unit
        !! report an iteration (print to the console).
        implicit none
        class(slsqp_solver),intent(inout) :: me
        integer,intent(in)                :: iter
        real(wp),dimension(:),intent(in)  :: x
        real(wp),intent(in)               :: f
        real(wp),dimension(:),intent(in)  :: c

        !write a header:
        if (iter==0) then
            write(output_unit,'(*(A20,1X))') 'iteration', 'x(1)', 'x(2)', 'f(1)', 'c(1)'
        end if

        !write the iteration data:
        write(output_unit,'(I20,1X,(*(F20.16,1X)))') iter,x,f,c

    end subroutine report_iteration

    end program slsqp_test_3
!*******************************************************************************

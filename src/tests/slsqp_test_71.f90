!*******************************************************************************
!> author: Jacob Williams
!
!  Test for the [[slsqp_module]].
!
!  This is problem 71 from the Hock-Schittkowsky test suite:
!
!```
!  min   x1*x4*(x1 + x2 + x3)  +  x3
!  s.t.  x1*x2*x3*x4  -  x5             -  25  =  0
!        x1**2 + x2**2 + x3**2 + x4**2  -  40  =  0
!        1 <=  x1,x2,x3,x4  <= 5
!        0 <=  x5
!
!  Starting point:
!     x = (1, 5, 5, 1, -24)
!
!  Optimal solution:
!     x = (1.00000000, 4.74299963, 3.82114998, 1.37940829, 0)
!```

    program slsqp_test_71

    use slsqp_module
    use slsqp_kinds

    implicit none

    integer,parameter               :: n = 5                    !! number of optimization variables
    integer,parameter               :: m = 2                    !! total number of constraints
    integer,parameter               :: meq = 2                  !! number of equality constraints
    integer,parameter               :: max_iter = 100           !! maximum number of allowed iterations
    real(wp),dimension(n),parameter :: xl = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]     !! lower bounds
    real(wp),dimension(n),parameter :: xu = [5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 1.0e10_wp]  !! upper bounds
    real(wp),parameter              :: acc = 1.0e-8_wp          !! tolerance
    integer,parameter               :: linesearch_mode = 1      !! use inexact linesearch.

    type(slsqp_solver)    :: solver      !! instantiate an slsqp solver
    real(wp),dimension(n) :: x           !! optimization variable vector
    integer               :: istat       !! for solver status check
    logical               :: status_ok   !! for initialization status check
    integer               :: iterations  !! number of iterations by the solver

    x = [1.0_wp, 5.0_wp, 5.0_wp, 1.0_wp, -24.0_wp] !initial guess

    call solver%initialize(n,m,meq,max_iter,acc,func,grad,&
                            xl,xu,linesearch_mode=linesearch_mode,&
                            status_ok=status_ok,&
                            report=report_iteration)

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

    contains

    subroutine func(me,x,f,c)

        implicit none
        class(slsqp_solver),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x      !! optimization variable vector
        real(wp),intent(out)              :: f      !! value of the objective function
        real(wp),dimension(:),intent(out) :: c      !! the constraint vector `dimension(m)`,
                                                    !! equality constraints (if any) first.

        f = x(1)*x(4)*(x(1)+x(2)+x(3)) + x(3)

        c(1) = x(1)*x(2)*x(3)*x(4) - x(5) - 25.0_wp
        c(2) = x(1)**2 + x(2)**2 + x(3)**2 + x(4)**2 - 40.0_wp

    end subroutine func

    subroutine grad(me,x,g,a)
        !! gradients for [[func]].
        implicit none
        class(slsqp_solver),intent(inout)   :: me
        real(wp),dimension(:),intent(in)    :: x    !! optimization variable vector
        real(wp),dimension(:),intent(out)   :: g    !! objective function partials w.r.t x `dimension(n)`
        real(wp),dimension(:,:),intent(out) :: a    !! gradient matrix of constraints w.r.t. x `dimension(m,n)`

        g(1) = x(4)*(2.0_wp*x(1)+x(2)+x(3))
        g(2) = x(1)*x(4)
        g(3) = x(1)*x(4) + 1.0_wp
        g(4) = x(1)*(x(1)+x(2)+x(3))
        g(5) = 0.0_wp

        a = 0.0_wp
        a(1,1) = x(2)*x(3)*x(4)
        a(1,2) = x(1)*x(3)*x(4)
        a(1,3) = x(1)*x(2)*x(4)
        a(1,4) = x(1)*x(2)*x(3)
        a(1,5) = -1.0_wp
        a(2,1) = 2.0_wp*x(1)
        a(2,2) = 2.0_wp*x(2)
        a(2,3) = 2.0_wp*x(3)
        a(2,4) = 2.0_wp*x(4)

    end subroutine grad

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
            write(output_unit,'(*(A20,1X))') 'iteration', 'x(1)', 'x(2)', 'x(3)', 'x(4)', 'x(5)', 'f', 'c(1)', 'c(2)'
        end if

        !write the iteration data:
        write(output_unit,'(I20,1X,(*(F20.16,1X)))') iter,x,f,c

    end subroutine report_iteration

    end program slsqp_test_71
!*******************************************************************************
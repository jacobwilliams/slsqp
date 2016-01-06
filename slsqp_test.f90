    program slsqp_test

    use slsqp_module
    use iso_fortran_env,    only: wp => real64

    implicit none

    type(slsqp_solver)              :: solver
    integer,parameter               :: n = 2
    integer,parameter               :: m = 1
    integer,parameter               :: meq = 1
    integer,parameter               :: max_iter = 100
    real(wp),dimension(n),parameter :: xl = [-1.0_wp, -1.0_wp]
    real(wp),dimension(n),parameter :: xu = [ 1.0_wp,  1.0_wp]
    real(wp),parameter              :: acc = 1.0e-8_wp

    real(wp),dimension(n) :: x

    x = [0.1_wp, 0.1_wp] !initial guess

    call solver%initialize(n,m,meq,max_iter,acc,rosenbrock_func,rosenbrock_grad,xl,xu)
    call solver%optimize(x)

    !solution:  x1 = 0.78641515415520813
    !           x2 = 0.61769831254040897
    !           f  = 4.5674808719472458E-002
    !           c  = 2.2848389846785722E-013

    contains

    subroutine rosenbrock_func(me,x,f,c)
        !! rosenbrock function computation
        !! see: http://www.mathworks.com/help/optim/ug/example-nonlinear-constrained-minimization.html
        implicit none
        class(slsqp_solver),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x      !! optimization variable vector
        real(wp),intent(out)              :: f      !! value of the objective function
        real(wp),dimension(:),intent(out) :: c      !! the constraint vector `dimension(m)`,
                                                    !! equality constraints (if any) first.

        f = 100.0_wp*(x(2) - x(1)**2)**2 + (1.0_wp - x(1))**2  !objective functdion
        c(1) = x(1)**2 + x(2)**2 - 1.0_wp  !equality constraint

    end subroutine rosenbrock_func

    subroutine rosenbrock_grad(me,x,g,a)
        implicit none
        class(slsqp_solver),intent(inout)   :: me
        real(wp),dimension(:),intent(in)    :: x    !! optimization variable vector
        real(wp),dimension(:),intent(out)   :: g    !! objective function partials w.r.t x `dimension(n)`
        real(wp),dimension(:,:),intent(out) :: a    !! gradient matrix of constraints w.r.t. x `dimension(m,n)`

        g(1) = -400.0_wp*(x(2)-x(1)**2)*x(1) - 2.0_wp*(1.0_wp-x(1))  !df/x1
        g(2) = 200.0_wp*(x(2) - x(1)**2)                             !df/x2

        a(1,1) = 2.0_wp*x(1)     ! dc/dx1
        a(1,2) = 2.0_wp*x(2)     ! dc/dx2

    end subroutine rosenbrock_grad

    end program slsqp_test

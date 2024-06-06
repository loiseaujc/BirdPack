module test_gauss_seidel
    use BirdPack
    use stdlib_io_npy, only: save_npy
    use testdrive, only: new_unittest, unittest_type, error_type, check
    implicit none
    private

    public :: collect_gauss_seidel_testsuite

    contains

    !------------------------------------------------------------
    !-----     DEFINITIONS OF THE UNIT TESTS FOR JACOBI     -----
    !------------------------------------------------------------

    subroutine collect_gauss_seidel_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ new_unittest("Manufactured solution for GS", manufactured_solution_gauss_seidel) ]

        return
    end subroutine collect_gauss_seidel_testsuite

    subroutine manufactured_solution_gauss_seidel(error)
        ! Error type to be returned.
        type(error_type), allocatable, intent(out) :: error
        ! Test problem.
        integer, parameter :: nx = 64+2, ny = 64+2
        real(dp), parameter :: Lx = 1.0_dp, Ly = 1.0_dp
        type(grid2D) :: grid
        real(dp), dimension(nx, ny) :: sol, rhs
        real(dp) :: abserror
        real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

        ! Creates the grid for the problem.
        grid = set_grid(nx, ny, Lx, Ly)

        ! Generate rhs and corresponding solution.
        rhs = create_rhs(grid) ; sol = make_solution(grid)

        ! Jacobi iteration.
        block
        real(dp), dimension(nx, ny) :: u = 0.0_dp
        real(dp), allocatable :: residual(:)
        real(dp) :: tol
        integer :: maxiter

        ! Set options.
        tol = atol ; maxiter = 10**6

        ! Jacobi solver.
        u = 0.0_dp ; call redblack_gauss_seidel_solver(u, rhs, grid, residual, tol, maxiter)

        ! Error estimate.
        abserror = maxval(abs(u - sol)) * grid%dx*grid%dy
        end block

        call check(error, abserror <= 1.0e-6_dp)

    contains
        pure function create_rhs(mesh) result(rhs)
            type(grid2D), intent(in) :: mesh
            real(dp), dimension(mesh%nx, mesh%ny) :: rhs
            integer :: i, j
            real(dp), parameter :: pi = 4*atan(1.0_dp)

            associate(dx => mesh%dx, dy => mesh%dy)
            do concurrent (i=1:mesh%nx, j=1:mesh%ny)
                rhs(i, j) = -8*pi**2 * sin((i-1)*dx*(2*pi)) * sin((j-1)*dy*(2*pi))
            enddo
            end associate

            return
        end function create_rhs

        pure function make_solution(mesh) result(sol)
            type(grid2D), intent(in) :: mesh
            real(dp), dimension(mesh%nx, mesh%ny) :: sol
            integer :: i, j
            real(dp), parameter :: pi = 4*atan(1.0_dp)

            associate(dx => mesh%dx, dy => mesh%dy)
             do concurrent (i=1:mesh%nx, j=1:mesh%ny)
                sol(i, j) = sin((i-1)*dx*(2*pi))*sin((j-1)*dy*(2*pi))
            enddo
            end associate

            return
        end function make_solution
    end subroutine manufactured_solution_gauss_seidel

end module test_gauss_seidel

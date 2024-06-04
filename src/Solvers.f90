module birdpack_solvers
  use iso_fortran_env
  use stdlib_optval, only: optval
  use birdpack_constants
  use birdpack_grid
  implicit none
  private

  public :: jacobi_solver

contains

  subroutine jacobi_solver(u, rhs, grid, residual, tol, maxiter)
    real(dp), intent(inout) :: u(:, :)
    !! Initial guess/solution of the Poisson equation.
    real(dp), intent(in) :: rhs(:, :)
    !! Right-hand side of the Poisson equation.
    type(grid2D), intent(in) :: grid
    !! Derived-type containing the grid's information.
    real(dp), allocatable, intent(out) :: residual(:)
    !! Time-history of the residual.
    real(dp), optional, intent(in) :: tol
    !! Tolerance of the solver.
    integer, optional, intent(in) :: maxiter
    !! Maximum number of iterations for the solver.

    ! Optional arguments.
    real(dp) :: tolerance
    integer  :: maxiterations

    !--------------------------------------
    !-----     OPTIONAL ARGUMENTS     -----
    !--------------------------------------

    tolerance     = optval(tol, atol)
    maxiterations = optval(maxiter, 10000)
    residual = [1.0_dp]

    !---------------------------------
    !-----     JACOBI SOLVER     -----
    !---------------------------------

    block
    ! Temporary array to store the updated solution.
    real(dp) :: up(grid%nx, grid%ny)
    ! Interation indices.
    integer :: i, j, k, iteration
    ! Residual.
    real(dp) :: absolute_error = 1.0_dp
    ! Coefficients.
    real(dp) :: alpha

    associate( dx => grid%dx, dy => grid%dy)
    iteration = 1 ; up = 0.0_dp

    alpha = 2/dx**2 + 2/dy**2 ; alpha = 1.0_dp / alpha

    ! Jacobi iteration.
    jacobi: do while (absolute_error > tolerance)
      ! Reset the absolute error.
      absolute_error = 0.0_dp

      ! Jacobi iteration.
      do concurrent (i=2:grid%nx-1, j=2:grid%ny-1)
          ! Jacobi update.
          up(i, j) = alpha*(u(i+1, j) + u(i-1, j))/dx**2 & ! Derivative in the x-direction.
                   + alpha*(u(i, j+1) + u(i, j-1))/dy**2 & ! Derivative in the y-direction.
                   - alpha*rhs(i, j)                       ! Source term.
      enddo

      do concurrent (i=2:grid%nx-1, j=2:grid%ny-1)
          ! Jacobi update.
          u(i, j) = alpha*(up(i+1, j) + up(i-1, j))/dx**2 & ! Derivative in the x-direction.
                  + alpha*(up(i, j+1) + up(i, j-1))/dy**2 & ! Derivative in the y-direction.
                  - alpha*rhs(i, j)                         ! Source term.
          ! On-the-fly residual computation.
          absolute_error = max(absolute_error, abs(up(i, j) - u(i, j)))
      enddo

      ! Update residual time-history.
      residual = [residual, absolute_error]

      ! Update iteration counter.
      iteration = iteration + 2 ; if (iteration > maxiterations) exit jacobi
    enddo jacobi

    end associate
    end block

    return
  end subroutine jacobi_solver

end module birdpack_solvers

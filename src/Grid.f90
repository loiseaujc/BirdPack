module birdpack_grid
  use birdpack_constants
  implicit none
  private

  public :: set_grid

  type, public :: grid2D
    integer :: nx, ny
    !! Number of grid points in the grid.
    real(dp) :: Lx, Ly
    !! Dimension of the spatial domain.
    real(dp) :: dx, dy
    !! Grid spacing in both directions.
  contains
    private
  end type

contains

  pure function set_grid(nx, ny, Lx, Ly) result(grid)
    integer, intent(in) :: nx, ny
    !! Number of grid points in both directions.
    real(dp), intent(in) :: Lx, Ly
    !! Dimension of the spatial domain.
    real(dp) :: dx, dy
    !! Grid spacing in both directions.
    type(grid2D) :: grid
    !! Derived-type to store the grid information.

    ! Compute the grid spacings.
    dx = Lx / (nx-1) ; dy = Ly / (ny-1)

    ! Create the grid object.
    grid = grid2d(nx, ny, Lx, Ly, dx, dy)

    return
  end function set_grid

end module birdpack_grid

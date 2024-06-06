module BirdPack

  ! Definitions of all the constants.
  use birdpack_constants
  ! Definition of the derived-type for the grid.
  use birdpack_grid
  ! Definition of the various solvers.
  use birdpack_solvers

  implicit none
  private

  !----------------------------------------------------------
  !-----     DEFINITION OF THE BASE NUMERICAL TYPES     -----
  !----------------------------------------------------------

  ! Double precision kind.
  public :: dp
  ! Default absolute tolerance.
  public :: atol
  ! Default relative tolerance.
  public :: rtol

  !--------------------------------------------------------
  !-----     DEFINITION OF GRID-RELATED FUNCTIONS     -----
  !--------------------------------------------------------

  ! Derived-type to store the grid's information.
  public :: grid2D
  ! Function to setup the grid.
  public :: set_grid

  !-------------------------------------
  !-----     ITERATIVE SOLVERS     -----
  !-------------------------------------

  public :: jacobi_solver
  public :: redblack_gauss_seidel_solver

contains

end module BirdPack

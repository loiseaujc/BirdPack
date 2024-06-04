module birdpack_constants
  implicit none
  private

  integer, parameter, public :: dp = selected_real_kind(15, 307)
  !! Definition of the double precision data type.
  real(dp), parameter, public :: atol = 10.0_dp ** -precision(1.0_dp)
  !! Absolute tolerance used by the solvers.
  real(dp), parameter, public :: rtol = sqrt(atol)
  !! Relative tolerance used by the solvers.

end module birdpack_constants

program Tester
    ! Fortran best practice.
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    ! Unit-test utilities.
    use testdrive, only: run_testsuite, new_testsuite, testsuite_type
    use test_jacobi
    use test_gauss_seidel
    implicit none

    ! Unit-test related.
    integer :: status, is, num_tests
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("+", *(1x, a))'

    status = 0

    testsuites = [ &
                   new_testsuite("Jacobi Test Suite", collect_jacobi_testsuite), &
                   new_testsuite("Gauss-Seidel Test Suite", collect_gauss_seidel_testsuite) &
                 ]

   do is = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(is)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(is)%collect, error_unit, status)
      write (*, *)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
      error stop
   else if (status == 0) then
      write (*, *) "All tests successfully passed!"
      write (*, *)
   end if

end program Tester

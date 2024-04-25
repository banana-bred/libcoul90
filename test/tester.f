! ================================================================================================================================ !
program tester
  !! Driver for units testing

  use iso_fortran_env,  only : error_unit
  use testdrive,        only : run_testsuite, new_testsuite,  testsuite_type, select_suite, run_selected, get_argument
  use test_coul90,  only : collect_coul90

  implicit none

  integer :: stat
  integer :: is

  character(len=:), allocatable :: suite_name
  character(len=:), allocatable :: test_name

  type(testsuite_type), allocatable :: testsuites(:)

  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  ! -- the tests to run
  testsuites = [                                     &
    new_testsuite("COUL90",   collect_coul90)&
  ]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then

    is = select_suite(testsuites, suite_name)

    if (is > 0 .AND. is <= size(testsuites)) then

      if (allocated(test_name)) then

        write(error_unit, fmt) "Suite:", testsuites(is)%name
        call run_selected( testsuites(is)%collect, test_name, error_unit, stat)

        if (stat < 0) error stop 1

      else

        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)

      end if

    else

      write(error_unit, fmt) "Available testsuites"

      do is = 1, size(testsuites)
        write(error_unit, fmt) "-", testsuites(is)%name
      end do

      error stop 1

    end if

  else

    do is = 1, size(testsuites)
      write(error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

  end if

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

! ================================================================================================================================ !
end program tester
! ================================================================================================================================ !

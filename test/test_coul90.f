! =================================================================================================================================!
module test_coul90
! =================================================================================================================================!
  use testdrive,  only: new_unittest, unittest_type, error_type, check
  use libcoul90,  only: coulf, coulg

  implicit none

  private

  public collect_coul90

! =================================================================================================================================!
contains
! =================================================================================================================================!

  ! -------------------------------------------------------------------------------------------------------------------------------!
  subroutine collect_coul90(testsuite)
    !! Collect all exported unit tests

    type(unittest_type), allocatable, intent(out) :: testsuite(:)
      !! Collection of tests to run

    testsuite = [                                                             &
      new_unittest("Test coulf",  test_coulf)  &
    ]
      ! new_unittest("invalid", test_invalid, should_fail=.true.) &

  end subroutine collect_coul90

  ! -------------------------------------------------------------------------------------------------------------------------------!
  subroutine test_coulf(error)

    use libcoul90,       only: coulf
    use iso_fortran_env, only: dp => real64

    type(error_type), allocatable, intent(out) :: error

!     integer, parameter :: nx   = 20
!     integer, parameter :: neta = 20

!     integer :: ix, ieta

!     real(dp), allocatable :: xvals(:)
!     real(dp), allocatable :: etavals(:)

!     xvals = [ -50.0_dp,  -50.0_dp, -4.0_dp, -4.0_dp, 10.0_dp, 10.0_dp, 100.0_dp ]
!     etavals = [ 5, 50  ]

    integer :: l = 0
    real(dp) :: eta = -0.5_dp
    real(dp) :: x = 20.0_dp
    real(dp) :: f
    f = coulf(l, eta, x)

    call check(error, f, -0.10237230180742793_dp)

    if(allocated(error)) return

  end subroutine test_coulf


! =================================================================================================================================!
end module test_coul90
! =================================================================================================================================!

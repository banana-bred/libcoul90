! ================================================================================================================================ !
module libcoul90__constants
  !! Contains constants used throughout coulcc

  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: ZERO
  public :: DZERO
  public :: HALF
  public :: ONE
  public :: TWO
  public :: THREE
  public :: SIX
  public :: TEN
  public :: TEN2
  public :: RT2DPI
  public :: RL35
  public :: LOGE
  public :: SMALL
  public :: ACCUR

  real(dp), parameter :: ZERO   = 0.0_dp
  real(dp), parameter :: DZERO  = ZERO
  real(dp), parameter :: HALF   = 0.5_dp
  real(dp), parameter :: ONE    = 1.0_dp
  real(dp), parameter :: TWO    = 2.0_dp
  real(dp), parameter :: THREE  = 3.0_dp
  real(dp), parameter :: SIX    = 6.0_dp
  real(dp), parameter :: TEN    = 10.0_dp
  real(dp), parameter :: TEN2   = 100.0_dp
  real(dp), parameter :: RT2DPI = 0.797884560802865_dp ! -- sqrt(2/π)
  real(dp), parameter :: RL35   = 35.0_dp
  real(dp), parameter :: LOGE   = 0.4242945_dp
  real(dp), parameter :: SMALL  = 1e-150_dp
  ! real(dp), parameter :: ACCUR = 100 * epsilon(ZERO)
  real(dp), parameter :: ACCUR = 1e-14_dp
  ! real(qp), parameter :: RT2DPI = 0.79788456080286535587989211986876373_dp ! -- sqrt(2/π)

! ================================================================================================================================ !
end module libcoul90__constants
! ================================================================================================================================ !

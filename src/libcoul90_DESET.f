! ================================================================================================================================ !
module libcoul90_deset
  !! Replaces the COMMON DESET block

  use iso_fortran_env, only: dp => real64

  implicit none

  save

  private

  public :: CF1
  public :: P
  public :: Q
  public :: F
  public :: GAMM
  public :: WRONSK

  real(dp) :: CF1
  real(dp) :: P
  real(dp) :: Q
  real(dp) :: F
  real(dp) :: GAMM
  real(dp) :: WRONSK

! ================================================================================================================================ !
end module libcoul90_deset
! ================================================================================================================================ !

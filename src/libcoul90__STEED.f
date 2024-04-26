! ================================================================================================================================ !
module libcoul90__steed
  !! Replaces the COMMON STEED block

  use iso_fortran_env, only: dp => real64

  implicit none

  save

  private

  public :: NFP
  public :: NPQ
  public :: IEXP
  public :: MINL
  public :: PACCQ

  integer  :: NFP
  integer  :: NPQ
  integer  :: IEXP
  integer  :: MINL
  real(dp) :: PACCQ

! ================================================================================================================================ !
end module libcoul90__steed
! ================================================================================================================================ !

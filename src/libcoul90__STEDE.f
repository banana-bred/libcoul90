! ================================================================================================================================ !
module libcoul90__stede
  !! Replaces the COMMON /STEDE/ block in sbesjy
  use iso_fortran_env, only: dp => real64

  implicit none

  private

  public :: ERR
  public :: NFP
  public :: NF

  integer  :: NF
  integer  :: NFP
  real(dp) :: ERR

! ================================================================================================================================ !
end module libcoul90__stede
! ================================================================================================================================ !

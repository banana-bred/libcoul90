! ================================================================================================================================ !
module libcoul90
  !! Provides interfaces and wrappers for calling `coul90()` that can be used as they are or as an example for
  !! implenting something similar. `coul90()` is a program (contained in the `coul90_submodule` submodule) that
  !! calculates complex Coulomb functions using Steed's method. `coul90` can also calculate other fuctions,
  !! such as spherical and cylindrical Bessel. See the implementation of the procedure `coul90_wrapper()`
  !! in this module or the implementation of the procedure `coul90()` in coul90_submodule.f for more details.

  implicit none

  private

  public :: coulf
  public :: coulg
  public :: sphbessj
  public :: sphbessy
  public :: cylbessj
  public :: cylbessy
  public :: coul90
  public :: coul90_wrapper
  public :: ricbes
  public :: sbesjy

  interface
    module subroutine coul90(X, ETA_IN, XLMIN, LRANGE, FC, GC, FCP, GCP, KFN, IFAIL)
      use iso_fortran_env, only: dp => real64
      integer, intent(in)  :: KFN
      integer, intent(in)  :: LRANGE
      integer, intent(out) :: IFAIL
      real(dp), intent(in) :: X
      real(dp), intent(in) :: XLMIN
      real(dp), intent(in) :: ETA_IN
      real(dp), intent(out) :: FC(0:)
      real(dp), intent(out) :: GC(0:)
      real(dp), intent(out) :: FCP(0:)
      real(dp), intent(out) :: GCP(0:)
    end subroutine
  end interface

  interface
    module subroutine ricbes(x,lmax, psi,chi,psid,chid, ifail)
      use iso_fortran_env, only: dp => real64
      integer, parameter :: LIMIT = 20000
      integer, parameter :: MAXL = 1001
      integer, intent(in) :: LMAX
      integer, intent(out) :: IFAIL
      real(dp), intent(in) :: X
      real(dp), intent(out) :: PSI(0:MAXL)
      real(dp), intent(out) :: PSID(0:MAXL)
      real(dp), intent(out) :: CHI(0:MAXL)
      real(dp), intent(out) :: CHID(0:MAXL)
    end subroutine ricbes
  end interface

  interface
    module subroutine sbesjy(X,LMAX, J,Y,JP,YP, IFAIL )
      use iso_fortran_env, only: dp => real64
      integer, parameter :: LIMIT = 20000
      integer, parameter :: MAXL = 1001
      integer, intent(in) :: LMAX
      integer, intent(out) :: IFAIL
      real(dp), intent(in) :: X
      real(dp), intent(out) :: J(0:MAXL)
      real(dp), intent(out) :: JP(0:MAXL)
      real(dp), intent(out) :: Y(0:MAXL)
      real(dp), intent(out) :: YP(0:MAXL)
    end subroutine sbesjy
  end interface

  interface coulf
    !! Returns the regular Coulomb function F_λ(η, x)
    module procedure coulf_ilambda
    module procedure coulf_rlambda
  end interface coulf

  interface coulg
    !! Returns the irregular Coulomb function G_λ(η, x)
    module procedure coulg_ilambda
    module procedure coulg_rlambda
  end interface coulg

  interface sphbessj
    !! Returns the spherical Bessel function of the first kind j_λ(x)
    module procedure sphbessj_ilambda
    module procedure sphbessj_rlambda
  end interface sphbessj

  interface sphbessy
    !! Returns the spherical Bessel function of the second kind y_λ(x)
    module procedure sphbessy_ilambda
    module procedure sphbessy_rlambda
  end interface sphbessy

  interface cylbessj
    !! Returns the (cylindrical) Bessel function of the first kind J_λ(x)
    module procedure cylbessj_ilambda
    module procedure cylbessj_rlambda
  end interface cylbessj

  interface cylbessy
    !! Returns the (cylindrical) Bessel function of the second kind Y_λ(x)
    module procedure cylbessy_ilambda
    module procedure cylbessy_rlambda
  end interface cylbessy

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  !------------------------------------------------------------------------------------------------------------------------------- !
  subroutine coul90_wrapper(xlmin, nl, eta, x, f, fp, g, gp, kfn)
    !! A wrapper routine to call coul90. The integers MODE and KFN
    !! must be set as follows to obtain the desired output from
    !! coul90 :
    !!
    !! The contents of F, FP, G, and GP  are controlled by KFN.
    !! ________________________________
    !!| ARRAY |          KFN          |
    !!|       |   0   |   1   |   2   |
    !!|_______|_______|_______|_______|
    !!|   F   | F_λ   |  j_λ  |  J_λ  |
    !!|   G   | G_λ   |  y_λ  |  Y_λ  |
    !! --------------------------------
    !! FP and GP are the derivatives of the functions represented by F and G,
    !! respectively.
    use iso_fortran_env, only: dp => real64
    real(dp),  intent(in)   :: xlmin
      !! The minimum value of λ
    integer,  intent(in)    :: nl
      !! The number of λ values
    integer,  intent(in)    :: kfn
      !! Controls which functions are represented by the output arrays
    real(dp), intent(in) :: eta
      !! The Sommerfeld parameter η in F_λ(η, x)
    real(dp), intent(in) :: x
      !! The distance x in F_λ(η, x)
    real(dp), intent(inout) :: f(:)
      !! The array containing values of F.
    real(dp), intent(inout) :: g(:)
      !! The array containing values of G.
    real(dp), intent(inout) :: fp(:)
      !! The array containing values of FP.
    real(dp), intent(inout) :: gp(:)
      !! The array containing values of GP.
    integer :: ifail
    integer :: lrange
    lrange = nl - 1
    call coul90(x, eta, xlmin, lrange, f, g, fp, gp, kfn, ifail)
    select case(ifail)
    case(0)
      continue
    case(1)
      call die("CF 1 DID NOT CONVERGE AFTER LIMIT ITERATIONS")
    case(2)
      call die("CF 2 DID NOT CONVERGE AFTER LIMIT ITERATIONS")
    case(-1)
      call die("X < SQRT(ACCUR)")
    case(-2)
      call die("INCONSISTENCY IN ORDER VALUES (L-VALUES)")
    case default
      call die("unknown error code from COULCC in WCLBES")
    end select
  end subroutine coul90_wrapper

  !------------------------------------------------------------------------------------------------------------------------------- !
  !  COULF INTERFACE
  !------------------------------------------------------------------------------------------------------------------------------- !

  !------------------------------------------------------------------------------------------------------------------------------- !
  function coulf_ilambda(lambda, eta, x) result(res)
    !! Returns the regular Coulomb function F_λ(η, x) for integer λ
    use iso_fortran_env, only: dp => real64
    integer,  intent(in) :: lambda
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 0
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = real(lambda, kind = dp)
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = f(nl)
  end function coulf_ilambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  function coulf_rlambda(lambda, eta, x) result(res)
    !! Returns the regular Coulomb function F_λ(η, x) for real λ
    use iso_fortran_env, only: dp => real64
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 0
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = lambda
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = f(nl)
  end function coulf_rlambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  !  COULG INTERFACE
  !------------------------------------------------------------------------------------------------------------------------------- !

  !------------------------------------------------------------------------------------------------------------------------------- !
  function coulg_ilambda(lambda, eta, x) result(res)
    !! Returns the irregular Coulomb function G_λ(η, x) for integer λ
    use iso_fortran_env, only: dp => real64
    integer,  intent(in) :: lambda
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 0
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = real(lambda, kind = dp)
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = g(nl)
  end function coulg_ilambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  function coulg_rlambda(lambda, eta, x) result(res)
    !! Returns the irregular Coulomb function G_λ(η, x) for real λ
    use iso_fortran_env, only: dp => real64
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: eta
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 0
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = lambda
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = g(nl)
  end function coulg_rlambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  !  SPHBESSJ INTERFACE
  !------------------------------------------------------------------------------------------------------------------------------- !

  !------------------------------------------------------------------------------------------------------------------------------- !
  function sphbessj_ilambda(lambda, x) result(res)
    !! Returns the spherical Bessel function of the first kind j_λ(x) for integer λ
    use iso_fortran_env, only: dp => real64
    integer,  intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 1
    real(dp) :: xlambda
    real(dp) :: eta = 0.0_dp
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = real(lambda, kind = dp)
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = f(nl)
  end function sphbessj_ilambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  function sphbessj_rlambda(lambda, x) result(res)
    !! Returns the spherical Bessel function of the first kind j_λ(x) for real λ
    use iso_fortran_env, only: dp => real64
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 1
    real(dp) :: eta = 0.0_dp
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = lambda
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = f(nl)
  end function sphbessj_rlambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  !  SPHBESSY INTERFACE
  !------------------------------------------------------------------------------------------------------------------------------- !

  !------------------------------------------------------------------------------------------------------------------------------- !
  function sphbessy_ilambda(lambda, x) result(res)
    !! Returns the spherical Bessel function of the second kind y_λ(x) for integer λ
    use iso_fortran_env, only: dp => real64
    integer,  intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 1
    real(dp) :: eta = 0.0_dp
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = real(lambda, kind = dp)
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = g(nl)
  end function sphbessy_ilambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  function sphbessy_rlambda(lambda, x) result(res)
    !! Returns the spherical Bessel function of the second kind y_λ(x) for real λ
    use iso_fortran_env, only: dp => real64
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 1
    real(dp) :: eta = 0.0_dp
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = lambda
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = g(nl)
  end function sphbessy_rlambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  !  CYLBESSJ INTERFACE
  !------------------------------------------------------------------------------------------------------------------------------- !

  !------------------------------------------------------------------------------------------------------------------------------- !
  function cylbessj_ilambda(lambda, x) result(res)
    !! Returns the (cylindrical) Bessel function of the second kind J_λ(x) for integer λ
    use iso_fortran_env, only: dp => real64
    integer,  intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 2
    real(dp) :: eta = 0.0_dp
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = real(lambda, kind = dp)
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = f(nl)
  end function cylbessj_ilambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  function cylbessj_rlambda(lambda, x) result(res)
    !! Returns the (cylindrical) Bessel function of the second kind J_λ(x) for real λ
    use iso_fortran_env, only: dp => real64
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 2
    real(dp) :: eta = 0.0_dp
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = lambda
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = f(nl)
  end function cylbessj_rlambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  !  CYLBESSY INTERFACE
  !------------------------------------------------------------------------------------------------------------------------------- !

  !------------------------------------------------------------------------------------------------------------------------------- !
  function cylbessy_ilambda(lambda, x) result(res)
    !! Returns the (cylindrical) Bessel function of the second kind Y_λ(x) for integer λ
    use iso_fortran_env, only: dp => real64
    integer,  intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 2
    real(dp) :: eta = 0.0_dp
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = real(lambda, kind = dp)
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = g(nl)
  end function cylbessy_ilambda

  !------------------------------------------------------------------------------------------------------------------------------- !
  function cylbessy_rlambda(lambda, x) result(res)
    !! Returns the (cylindrical) Bessel function of the second kind Y_λ(x) for real λ
    use iso_fortran_env, only: dp => real64
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: x
    real(dp) :: res
    integer, parameter :: nl  = 1
    integer, parameter :: kfn = 2
    real(dp) :: eta = 0.0_dp
    real(dp) :: xlambda
    real(dp) :: f(nl)
    real(dp) :: g(nl)
    real(dp) :: fp(nl)
    real(dp) :: gp(nl)
    xlambda = lambda
    call coul90_wrapper(xlambda, nl, eta, x, f, fp, g, gp, kfn)
    res = g(nl)
  end function cylbessy_rlambda

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine die(message)
    !! Stop program execution with a message
    use iso_fortran_env, only: stderr => error_unit
    character(*), intent(in), optional :: message
    if( .NOT. present(message)) error stop
    write(stderr,'("STOP",X,"::",X,A)') message
    write(stderr,*)
    error stop
  end subroutine die

! ================================================================================================================================ !
end module libcoul90
! ================================================================================================================================ !

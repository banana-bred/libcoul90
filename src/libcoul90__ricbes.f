! ================================================================================================================================ !
submodule (libcoul90) libcoul90__ricbes

  implicit none

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module SUBROUTINE RICBES(X,LMAX, PSI,CHI,PSID,CHID, IFAIL )
    !!   REAL RICCATI-BESSEL FUNCTIONS AND X-DERIVATIVES :
    !!   PSI = X . J/L/(X) , CHI = X . Y/L/(X)    FROM L=0 TO L=LMAX
    !!        FOR REAL X > SQRT(ACCUR) (E.G. 1D-7)  AND INTEGER LMAX
    !!
    !! PSI (L)  =      PSI/L/(X) STORES   REGULAR RICCATI BESSEL FUNCTION:
    !! PSID(L)  = D/DX PSI/L/(X)         PSI(0) =  SIN(X)
    !! CHI (L)  =      CHI/L/(X) STORES IRREGULAR RICCATI BESSEL FUNCTION:
    !! CHID(L)  = D/DX CHI/L/(X)         CHI(0) = -COS(X) [HANDBOOK DEFN.]
    !!
    !!    IFAIL = -1 FOR ARGUMENTS OUT OF RANGE
    !!          =  0 FOR ALL RESULTS SATISFACTORY
    !!
    !!   USING LENTZ-THOMPSON EVALUATION OF CONTINUED FRACTION CF1,
    !!   AND TRIGONOMETRIC FORMS FOR L = 0 SOLUTIONS.
    !!   LMAX IS LARGEST L NEEDED AND MUST BE <= MAXL, THE ARRAY INDEX.
    !!   MAXL CAN BE DELETED AND ALL THE ARRAYS DIMENSIONED (0:*)
    !!   SMALL IS MACHINE DEPENDENT, ABOUT SQRT(MINIMUM REAL NUMBER),
    !!         SO 1D-150 FOR DOUBLE PRECISION ON VAX, PCS ETC.
    !!   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
    !!   IN OSCILLATING REGION X .GE.  [ SQRT{LMAX*(LMAX+1)} ]
    !!   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8
    !!   IS THE SMALLEST NUMBER WITH 1+ACC8.NE.1 FOR OUR WORKING PRECISION
    !!   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
    !!   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
    !!   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
    !!   THE VARIABLE ERR IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
    !!
    !!   NOTE: FOR X=1 AND L=100   PSI = 7.4 E-190  CHI = -6.7+E186  1.4.94
    !!---------------------------------------------------------------------
    !!   AUTHOR :   A.R.BARNETT       MANCHESTER    12 MARCH 1990/95
    !!                                AUCKLAND      12 MARCH 1991
    !!---------------------------------------------------------------------

    use iso_fortran_env,      only: dp => real64, stderr => error_unit
    use libcoul90__stede,     only: ERR, NFP
    use libcoul90__constants, only: ZERO, ONE, TWO, THREE, ACCUR, SMALL

    integer, parameter :: LIMIT = 20000, MAXL = 1001

    integer, intent(in) :: LMAX
    integer, intent(out) :: IFAIL
    real(dp), intent(in) :: X
    real(dp), intent(out) :: PSI(0:MAXL)
    real(dp), intent(out) :: PSID(0:MAXL)
    real(dp), intent(out) :: CHI(0:MAXL)
    real(dp), intent(out) :: CHID(0:MAXL)
    ! real(dp) :: PSI (0:MAXL), CHI (0:MAXL)
    ! real(dp) :: PSID(0:MAXL), CHID(0:MAXL)

    integer :: L

    real(dp) :: TK,SL
    real(dp) :: XINV, CF1,DCF1, DEN, C,D, OMEGA, TWOXI

    ! PARAMETER ( ZERO  = 0.0D0  , ONE   = 1.0D0 , TWO = 2.0D0 )
    ! PARAMETER ( SMALL = 1.D-150, THREE = 3.0D0 )
    !OMMON /STEDE/    ERR,NFP       ! not required in code
    !-------
    ! ACCUR = 1.D-14                  ! suitable for double precision

    IFAIL = -1                      ! user to check on exit
    IF (X .LT. DSQRT(ACCUR) )       GOTO 50
    !-------TAKES CARE OF NEGATIVE X ... USE REFLECTION FORMULA
    !-------BEGIN CALCULATION OF CF1 UNLESS LMAX = 0, WHEN SOLUTIONS BELOW
    XINV  = ONE / X
    IF (LMAX .GT. 0) THEN
        TWOXI =     XINV + XINV
        SL  =  REAL(LMAX)* XINV     ! used also in do loop 3
        TK  =  TWO * SL  + XINV * THREE
        CF1 =  SL  + XINV           ! initial value of CF1
        DEN =  ONE                  ! unnormalised psi(Lmax,x)
        IF ( ABS(CF1) .LT. SMALL ) CF1 = SMALL
        C   = CF1                   ! inverse ratio of A convergents
        D   = ZERO                  ! direct  ratio of B convergents
        DO L = 1,LIMIT
        ! DO 10 L = 1,LIMIT
            C   = TK - ONE / C
            D   = TK - D
            IF ( ABS(C) .LT. SMALL ) C = SMALL
            IF ( ABS(D) .LT. SMALL ) D = SMALL
            D   = ONE / D
            DCF1= D   * C
            CF1 = CF1 * DCF1
            IF ( D .LT. ZERO ) DEN = - DEN
            IF ( ABS(DCF1 - ONE) .LE. ACCUR ) GOTO 20
           TK   = TK + TWOXI
           NFP  = L                 ! ie number in loop
         END DO
    ! 10     CONTINUE
         GOTO 50             ! error exit, no convergence
    20   CONTINUE
         ERR = ACCUR * DSQRT(DBLE(NFP))    ! error estimate
         PSI (LMAX) = DEN
         PSID(LMAX) = CF1 * DEN
    !------ DOWNWARD RECURSION TO L=0  AS RICCATI BESSEL FUNCTIONS
         ! DO 30 L =  LMAX , 1, -1
         DO L =  LMAX , 1, -1
             PSI (L-1) = SL * PSI(L)   +  PSID(L)
             PSID(L-1) = SL * PSI(L-1) -  PSI (L)
             SL        = SL - XINV
         END DO
    ! 30       CONTINUE
         DEN =PSI(0)
     ENDIF                        ! end loop for Lmax > 0
    !------ CALCULATE THE L=0   RICCATI BESSEL FUNCTIONS DIRECTLY
    PSI (0)   =  DSIN(X)
    !HI (0)   = -DCOS(X)
    PSID(0)   = -CHI(0)
    !HID(0)   =  PSI(0)
    IF (LMAX .GT. 0) THEN
        OMEGA  =PSI(0) / DEN
        ! DO 40 L = 1 , LMAX
        DO L = 1 , LMAX
              PSI (L) = OMEGA * PSI (L)
              PSID(L) = OMEGA * PSID(L)
                  SL  = XINV  * REAL(L)
              CHI (L) = SL * CHI(L-1)   - CHID(L-1)
              CHID(L) = CHI(L-1)  -  SL * CHI (L)
    ! 40       CONTINUE
        END DO
    ENDIF
    IFAIL = 0                       ! calculations successful
    RETURN
    !---------------------------------------------------------------------
    !       ERROR TRAPS
    !---------------------------------------------------------------------
    50   IF (X .LT. ZERO) THEN
            WRITE(stderr,1000) X
          ELSEIF (X .EQ. ZERO) THEN
            IFAIL = 0
            CHI(0) = ONE      ! irregular function
            ! DO 60 L = 1, LMAX
            DO L = 1, LMAX
                  CHI(L) = ZERO     ! remaining arrays untouched
            END DO
    ! 60           CONTINUE
      ELSE                          ! x .le. sqrt(accur), e.g. 1D-7
            WRITE(stderr,1001) X
      ENDIF
    1000   FORMAT(' X NEGATIVE !',1PE15.5,'    USE REFLECTION FORMULA'/)
    1001   FORMAT(' WITH X = ',1PE15.5,'    TRY SMALL-X SOLUTIONS',/, '  PSI/L/(X)  ->   X**(L+1) / (2L+1)!!       AND',/, &
       '  CHI/L/(X)  ->  -(2L-1)!! / X**L'/)
  END subroutine ricbes
  !---------------------------------------------------------------------
  !       END OF SUBROUTINE RICBES
  !---------------------------------------------------------------------

! ================================================================================================================================ !
end submodule libcoul90__ricbes
! ================================================================================================================================ !

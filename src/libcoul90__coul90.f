! ================================================================================================================================ !
submodule (libcoul90) libcoul90__coul90
  !! Submodule containing the implementation of the COUL90 procedure

  implicit none

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module SUBROUTINE COUL90(X, ETA_IN, XLMIN, LRANGE, FC, GC, FCP, GCP, KFN, IFAIL)
    !!  COULOMB & BESSEL FUNCTION PROGRAM-- COUL90 -- USING STEED'S METHOD
    !!
    !!  COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = (D/DX) F, GCP = (D/DX) G
    !!   FOR REAL X .GT. 0. ,REAL ETA (INCLUDING 0.), AND REAL XLMIN .GT.-1.
    !!   FOR (LRANGE+1) INTEGER-SPACED LAMBDA VALUES.
    !!   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
    !!   EQUATION, TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
    !!   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
    !!   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
    !!----------------------------------------------------------------------
    !!   CALLING VARIABLES; ALL REALS ARE real(dp) (REAL*8)
    !!
    !!   X       - REAL ARGUMENT FOR COULOMB FUNCTIONS > 0.0
    !!             [ X > SQRT(ACCUR) : ACCUR IS TARGET ACCURACY 1.0D-14 ]
    !!   ETA  - REAL SOMMERFELD PARAMETER, UNRESTRICTED > = < 0.0
    !!   XLMIN   - REAL MINIMUM LAMBDA-VALUE (L-VALUE OR ORDER),
    !!             GENERALLY IN RANGE 0.0 - 1.0 AND MOST USUALLY 0.0
    !!   LRANGE  - INTEGER NUMBER OF ADDITIONAL L-VALUES : RESULTS RETURNED
    !!             FOR L-VALUES XLMIN TO XLMIN + LRANGE INCLUSIVE
    !!   FC ,GC  - REAL VECTORS F,G OF REGULAR, IRREGULAR COULOMB FUNCTIONS
    !!   FCP,GCP - REAL VECTORS FOR THE X-DERIVATIVES OF  F,G
    !!             THESE VECTORS TO BE OF LENGTH AT LEAST MINL + LRANGE
    !!             STARTING ELEMENT MINL = MAX0( IDINT(XLMIN+ACCUR),0 )
    !!   KFN     - INTEGER CHOICE OF FUNCTIONS TO BE COMPUTED :
    !!           = 0         REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
    !!           = 1    SPHERICAL BESSEL      "      "     "        j & y
    !!           = 2  CYLINDRICAL BESSEL      "      "     "        J & Y
    !!
    !!   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
    !!   IN OSCILLATING REGION X .GE. [ETA + SQRT{ETA**2 + XLM*(XLM+1)}]
    !!   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8 IS
    !!   THE SMALLEST NUMBER WITH 1.+ACC8.NE.1. FOR OUR WORKING PRECISION.
    !!   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
    !!   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
    !!   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
    !!   THE VARIABLE PACCQ IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
    !!----------------------------------------------------------------------
    !!   ERROR RETURNS                THE USER SHOULD TEST IFAIL ON EXIT
    !!
    !!   IFAIL ON INPUT IS SET TO 0                        LIMIT = 20000
    !!   IFAIL IN OUTPUT =  0 : CALCULATIONS SATISFACTORY
    !!                   =  1 : CF1 DID NOT CONVERGE AFTER LIMIT ITERATIONS
    !!                   =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
    !!                   = -1 : X < 1D-7 = SQRT(ACCUR)
    !!                   = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES)
    !!----------------------------------------------------------------------
    !!  MACHINE-DEPENDENT PARAMETERS:    ACCUR - SEE ABOVE
    !!           SMALL - OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
    !!           IE 1D-30 FOR IBM REAL*8,    1D-150 FOR real(dp)
    !!----------------------------------------------------------------------
    !!  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
    !!  ORIGINAL PROGRAM  RCWFN       IN    CPC  8 (1974) 377-395
    !!                 +  RCWFF       IN    CPC 11 (1976) 141-142
    !!  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
    !!  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
    !!  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188
    !!  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
    !!  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
    !!  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
    !!  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509     1.4.94
    !!----------------------------------------------------------------------
    !!  AUTHOR: A. R. BARNETT           MANCHESTER  MARCH   1981/95
    !!                                  AUCKLAND    MARCH   1991
    !!----------------------------------------------------------------------

    use libcoul90__steed,     only: NFP, NPQ, IEXP, MINL, PACCQ
    use libcoul90__deset,     only: CF1, P, Q, F, GAMM, WRONSK
    use iso_fortran_env,     only: sp => real32, dp => real64, stderr => error_unit
    use libcoul90__constants, only: TWO, TEN2, RT2DPI

    integer, intent(in) :: LRANGE
    integer, intent(in) :: KFN
    integer, intent(out) :: IFAIL
    real(dp), intent(in) :: X
    real(dp), intent(in) :: XLMIN
    real(dp), intent(in) :: ETA_IN
    ! real(dp) :: FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
    ! real(dp), intent(out) :: FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
    real(dp), intent(out) :: FC (0:), GC (0:), FCP(0:), GCP(0:)
      !! ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM MINL

    integer, parameter  :: LIMIT = 20000
    real(dp), parameter :: SMALL = 1.0e-150_dp

    ! COMMON  /STEED/  PACCQ,NFP,NPQ,IEXP,MINL    !not required in code
    ! COMMON  /DESET/  CF1,P,Q,F,GAMM,WRONSK     !information only
    !----------------------------------------------------------------------
    !     COUL90 HAS CALLS TO: DSQRT,DABS,MAX0,IDINT,DSIGN,DFLOAT,DMIN1
    !----------------------------------------------------------------------
    ! DATA ZERO,ONE,TWO,TEN2,HALF /0.0D0, 1.0D0, 2.0D0, 1.0D2, 0.5D0/
    ! DATA RT2DPI /0.79788 45608 02865  D0/
    !Q    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
    !-----THIS CONSTANT IS  DSQRT(TWO / PI):
    !-----USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND real(dp)
    !----------------CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
    ! ACCUR = 1.0D-14

    real(dp), parameter :: ACC8 = epsilon(1.0_dp)
    real(dp), parameter :: ACCUR = 100 * ACC8

    integer :: L, MAXL
    logical :: ETANE0, XLTURN

    real(sp), parameter :: ZERO  = 0.0_sp
    real(sp), parameter :: HALF  = 0.5_sp
    real(sp), parameter :: ONE   = 1.0_sp
    real(sp), parameter :: SIX   = 6.0_sp
    real(sp), parameter :: TEN   = 10.0_sp
    real(sp), parameter :: RL35  = 35.0_sp
    real(sp), parameter :: ALOGE = 0.4342945_sp
    real(sp) :: GH2, XLL1, HLL, HL, SL, RL2, GH, PHI, PHI10
    real(dp) :: ACCH
    real(dp) :: XINV, PK, C, D, PK1, ETAK, RK2, TK, DCF1, DEN, XLM, XLL
    real(dp) :: EL, XL, RL, SL, FCMAXL, FCMINL, GCMINL, OMEGA
    real(dp) :: WI, A, B, AR, AI, BR, BI, DR, DI, DPP, DQ, ALPHA, BETA
    real(dp) :: E2MM1, FJWKB, GJWKB, GAMMAI
    real(dp) :: ETA

    ! -- so that ETA_IN may be of INTENT(IN)
    ETA = ETA_IN

    IFAIL = 0
    IEXP  = 1
    NPQ   = 0
    GJWKB = ZERO
    PACCQ = ONE
    IF(KFN /= 0) ETA = ZERO
    ETANE0  = ETA /= ZERO
    ACCH  = DSQRT(ACCUR)

    !-----   TEST RANGE OF X, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
    IF( X <= ACCH )                GO TO 100
    IF( KFN == 2 )   THEN
        XLM = XLMIN - HALF
    ELSE
        XLM = XLMIN
    ENDIF
    IF( XLM <= -ONE .OR. LRANGE < 0 )         GO TO 105
    E2MM1  = XLM * XLM + XLM
    XLTURN = X * (X -  TWO * ETA) < E2MM1
    E2MM1  = E2MM1  +  ETA * ETA
    XLL    = XLM + DFLOAT(LRANGE)

    !-----  LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
    !-----  XLL  IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
    !-----  DETERMINE STARTING ARRAY ELEMENT (MINL) FROM XLMIN
    MINL  = MAX0( IDINT(XLMIN + ACCUR),0 )     ! index from 0
    MAXL  = MINL + LRANGE

    !-----   EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
    XINV = ONE / X
    DEN  = ONE                       ! unnormalised F(MAXL,ETA,X)
    PK   = XLL + ONE
    CF1  = ETA / PK  +  PK * XINV
    IF( DABS(CF1) < SMALL )    CF1 = SMALL
    RK2  = ONE
    D = ZERO
    C = CF1

    !----- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA + 1: LENTZ-THOMPSON
    DO L =  1 , LIMIT             ! abort if reach LIMIT (20000)
        PK1 = PK + ONE
        IF( ETANE0 ) THEN
            ETAK = ETA / PK
            RK2  = ONE + ETAK * ETAK
            TK  = (PK + PK1) * (XINV + ETAK / PK1)
        ELSE
            TK  = (PK + PK1) * XINV
        ENDIF
        D   =  TK - RK2 * D          ! direct  ratio of B convergents
        C   =  TK - RK2 / C          ! inverse ratio of A convergents
        IF( DABS(C) < SMALL ) C = SMALL
        IF( DABS(D) < SMALL ) D = SMALL
        D   = ONE / D
        DCF1=   D * C
        CF1 = CF1 * DCF1
        IF( D < ZERO )    DEN = -DEN
        PK  = PK1
        IF( DABS(DCF1-ONE) < ACCUR )     GO TO  20 ! proper exit
    END DO
    GO TO 110 ! error exit
    20 NFP = PK - XLL - 1                        ! number of steps
    F = CF1                                 ! need DEN later

    !----DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
    IF( LRANGE > 0 )       THEN
        FCMAXL    = SMALL  * DEN
        FCP(MAXL) = FCMAXL * CF1
        FC (MAXL) = FCMAXL
        XL = XLL
        RL = ONE
        DO L =  MAXL, MINL+1, -1
            IF( ETANE0 )  THEN
                EL = ETA / XL
                RL = DSQRT( ONE + EL * EL )
                SL = XL * XINV  + EL
                GC (L) = RL                  ! storage
                GCP(L) = SL
            ELSE
                SL = XL * XINV
            ENDIF
            FC (L-1)  = ( FC(L)   * SL  +  FCP(L) ) / RL
            FCP(L-1)  =   FC(L-1) * SL  -  FC (L) * RL
            XL    =  XL - ONE                   ! end value is XLM
        END DO
        IF( DABS(FC(MINL)) < ACCUR*SMALL )  FC(MINL) = ACCUR * SMALL
        F   = FCP(MINL) / FC(MINL)             ! F'/F at min L-value
        DEN = FC (MINL)                        ! normalisation
    ENDIF

    !---------------------------------------------------------------------
    !-----   NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
    !-----   EVALUATE CF2 = P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)
    !---------------------------------------------------------------------

    IF( XLTURN ) CALL JWKB( X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP )
    IF( IEXP > 1 .OR. GJWKB > (ONE / (ACCH*TEN2)) ) THEN
        OMEGA = FJWKB
        GAMM = GJWKB * OMEGA
        gammai = ONE/GAMM
        P     = F
        Q     = ONE
    ELSE                                     ! find cf2
        XLTURN = .FALSE.
        PK =  ZERO
        WI =  ETA + ETA
        P  =  ZERO
        Q  =  ONE - ETA * XINV
        AR = -E2MM1
        AI =  ETA
        BR =  TWO * (X - ETA)
        BI =  TWO
        DR =  BR / (BR * BR + BI * BI)
        DI = -BI / (BR * BR + BI * BI)
        DPP = -XINV * (AR * DI + AI * DR)
        DQ =  XINV * (AR * DR - AI * DI)
        DO L = 1, LIMIT
            P  = P  + DPP
            Q  = Q  + DQ
            PK = PK + TWO
            AR = AR + PK
            AI = AI + WI
            BI = BI + TWO
            D  = AR * DR - AI * DI + BR
            DI = AI * DR + AR * DI + BI
            C  = ONE / (D * D + DI * DI)
            DR =  C * D
            DI = -C * DI
            A  = BR * DR - BI * DI - ONE
            B  = BI * DR + BR * DI
            C  = DPP * A  - DQ * B
            DQ = DPP * B  + DQ * A
            DPP = C
            IF( DABS(DPP)+DABS(DQ) < (DABS(P)+DABS(Q)) * ACCUR ) GO TO 50
        END DO
        GO TO 120 ! error exit
        50 NPQ   = PK / TWO                              ! proper exit
        PACCQ = HALF * ACCUR / DMIN1( DABS(Q),ONE )
        IF( DABS(P) > DABS(Q) ) PACCQ = PACCQ * DABS(P)
        !---------------------------------------------------------------------
        !    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
        !---------------------------------------------------------------------
        GAMM   = (F - P) / Q
        GAMMAI  = ONE / GAMM
        IF( DABS(GAMM) <= ONE )  THEN
            OMEGA  = DSQRT( ONE  +  GAMM * GAMM )
        ELSE
            OMEGA  = DSQRT( ONE  +  GAMMAI* GAMMAI) * DABS(GAMM)
        ENDIF
        OMEGA  = ONE / ( OMEGA * DSQRT(Q) )
        WRONSK = OMEGA
    ENDIF
    !---------------------------------------------------------------------
    !    RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
    !---------------------------------------------------------------------
    IF( KFN == 1 )       THEN         !   spherical Bessel functions
        ALPHA = XINV
        BETA  = XINV
    ELSEIF( KFN == 2 ) THEN         ! cylindrical Bessel functions
        ALPHA = HALF * XINV
        BETA  = DSQRT( XINV ) * RT2DPI
    ELSE                            ! kfn = 0,   Coulomb functions
        ALPHA = ZERO
        BETA  = ONE
    ENDIF
    FCMINL = DSIGN( OMEGA,DEN ) * BETA
    IF( XLTURN )   THEN
        GCMINL =   GJWKB * BETA
    ELSE
        GCMINL =  FCMINL * GAMM
    ENDIF
    IF( KFN /= 0 )    GCMINL = -GCMINL         ! Bessel sign differs
    FC (MINL) = FCMINL
    GC (MINL) = GCMINL
    GCP(MINL) = GCMINL * (P - Q * GAMMAI - ALPHA)
    FCP(MINL) = FCMINL * (F - ALPHA)
    IF( LRANGE == 0 )                          RETURN

    !---------------------------------------------------------------------
    !    UPWARD RECURRENCE FROM GC(MINL),GCP(MINL) STORED VALUES ARE RL,SL
    !    RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
    !      XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
    !---------------------------------------------------------------------
    OMEGA = BETA * OMEGA / DABS(DEN)
    XL = XLM
    RL = ONE
    DO L = MINL+1 , MAXL                   ! indexed from 0
        XL = XL + ONE
        IF( ETANE0 ) THEN
            RL = GC (L)
            SL = GCP(L)
        ELSE
            SL =  XL * XINV
        ENDIF
        GC (L)  = ( (SL - ALPHA) * GC(L-1) - GCP(L-1) ) / RL
        GCP(L)  =    RL *  GC(L-1)  -  (SL + ALPHA) * GC(L)
        FCP(L)  = OMEGA * ( FCP(L)  -  ALPHA * FC(L) )
        FC (L)  = OMEGA *   FC (L)
    END DO

    RETURN

    !------------------   ERROR MESSAGES
    100 IFAIL = -1
    WRITE(stderr, 1000) X,ACCH
    1000 FORMAT(' FOR X = ',1PD12.3,'     TRY SMALL-X  SOLUTIONS, OR X IS NEGATIVE'/ ,' SQUARE ROOT (ACCURACY) =  ',D12.3/)
    RETURN
    105 IFAIL = -2
    WRITE (stderr, 1005) LRANGE,XLMIN,XLM
    1005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ', I10,1P2D15.6/)
    RETURN
    110 IFAIL =  1
    WRITE (stderr, 1010) LIMIT, CF1,DCF1, PK,ACCUR
    1010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',I10,' ITERATIONS CF1,DCF1,PK,ACCUR =  ',1P4D12.3/)
    RETURN
    120 IFAIL =  2
    WRITE (stderr, 1020) LIMIT,P,Q,DPP,DQ,ACCUR
    1020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',I7,' ITERATIONS P,Q,DPP,DQ,ACCUR =  ',1P4D17.7,D12.3/)


  END SUBROUTINE COUL90

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  SUBROUTINE  JWKB(X, ETA, XL, FJWKB, GJWKB, IEXP)
    !!     COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
    !!     AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
    !!     CALCULATED IN SINGLE, RETURNED IN real(dp) VARIABLES
    !!     CALLS DMAX1, SQRT, LOG, EXP, ATAN2, FLOAT, INT
    !!     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991

    use iso_fortran_env,     only: dp => real64
    use libcoul90__constants, only: DZERO, ZERO, HALF, ONE, SIX, TEN, RL35, LOGE

    real(dp), intent(in) :: X
    real(dp), intent(in) :: ETA
    real(dp), intent(in) :: XL
    real(dp), intent(out) :: FJWKB
    real(dp), intent(out) :: GJWKB

    INTEGER, parameter :: MAXEXP = 300

    integer :: IEXP
    REAL(dp) :: GH2, XLL1, HLL, HL, SL, RL2, GH, PHI, PHI10

    !----------------------------------------------------------------------
    ! REAL(dp) ::    ZERO,HALF,ONE,SIX,TEN,RL35,LOGE
    ! PARAMETER  ( MAXEXP = 300 )
    ! DATA  ZERO,HALF,ONE,SIX,TEN  /0.0E0, 0.5E0, 1.0E0, 6.0E0, 1.0E1/
    ! DATA DZERO,RL35,LOGE /0.0D0, 35.0E0, 0.43429 45 E0 /
    !----------------------------------------------------------------------
    ! OOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR real(dp)
    !----------------------------------------------------------------------

    GH2   =  X * (ETA + ETA - X)
    XLL1  = DMAX1( XL * XL + XL, DZERO )
    IF( GH2 + XLL1 <= ZERO )                 RETURN
    HLL  = XLL1 + SIX / RL35
    HL   = SQRT(HLL)
    SL   = ETA / HL + HL / X
    RL2  = ONE + ETA * ETA / HLL
    GH   = SQRT(GH2 + HLL) / X
    PHI  = X*GH - HALF*( HL*LOG((GH + SL)**2 / RL2) - LOG(GH) )
    IF ( ETA /= ZERO ) PHI = PHI - ETA * ATAN2(X*GH,X - ETA)
    PHI10 = -PHI * LOGE
    IEXP  =  INT(PHI10)
    IF ( IEXP > MAXEXP ) THEN
        GJWKB = TEN**(PHI10 - FLOAT(IEXP))
    ELSE
        GJWKB = EXP(-PHI)
        IEXP  = 0
    ENDIF
    FJWKB = HALF / (GH * GJWKB)

  END SUBROUTINE

! ================================================================================================================================ !
end submodule libcoul90__coul90
! ================================================================================================================================ !

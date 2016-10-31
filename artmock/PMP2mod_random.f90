!------------------------------------------------------------------- 
    Module LUXURY
      Integer*4, Allocatable :: i24A(:),j24A(:),in24A(:),        &
                                kountA(:), iFlagA(:),SeedsPage(:,:)
      Real*4 ,   Allocatable ::  carryA(:),gSetA(:)
    end Module LUXURY
Module Random
contains
 
!--------------------------------------                                 
!				normal random numbers                                              
      FUNCTION GAUSS (M) 
!--------------------------------------                                 
      X = 0. 
      DO I = 1, 5 
      X = X + RANDd (M) 
      EndDo 
      X2 = 1.5491933 * (X - 2.5) 
      GAUSS = X2 * (1. - 0.01 * (3. - X2 * X2) ) 
      RETURN 
      END FUNCTION GAUSS                            

!------------------------------------------------                       
!				                                       random number generator     
      FUNCTION RANDd (M) 
!------------------------------------------------                       
      DATA LC, AM, KI, K1, K2, K3, K4, L1, L2, L3, L4	 / 453815927,     &
      2147483648., 2147483647, 536870912, 131072, 256, 	16777216, 4,    &
      16384, 8388608, 128 /
      ML = M / K1 * K1 
      M1 = (M - ML) * L1 
      ML = M / K2 * K2 
      M2 = (M - ML) * L2 
      ML = M / K3 * K3 
      M3 = (M - ML) * L3 
      ML = M / K4 * K4 
      M4 = (M - ML) * L4 
      M5 = KI - M 
      IF (M1.GE.M5) M1 = M1 - KI - 1 
      ML = M + M1 
      M5 = KI - ML 
      IF (M2.GE.M5) M2 = M2 - KI - 1 
      ML = ML + M2 
      M5 = KI - ML 
      IF (M3.GE.M5) M3 = M3 - KI - 1 
      ML = ML + M3 
      M5 = KI - ML 
      IF (M4.GE.M5) M4 = M4 - KI - 1 
      ML = ML + M4 
      M5 = KI - ML 
      IF (LC.GE.M5) ML = ML - KI - 1 
      M = ML + LC 
      RANDd = M / AM 
      RETURN 
    END FUNCTION RANDd

    !--------------------------------------                                 
!		normal random numbers                                                
      FUNCTION GAUSS3 (gSet,iFlag) 
!                        Uses ranlux for homogenous rand numbers        
!                        Uses Box-Muller inversion for gaussian         
!                 Excellent quality and speed                           
!  N=       100000000  GAUSS3+ranlux ------                             
!  sigma= 0.99997     mean= 0.18692E-03                                 
!  sigma Frac(>sig) Frac(<-sig)   True       n(>)    n(<)               
!  1.00 0.1587     0.1586     0.1587     15869095 15857724              
!  1.50 0.6682E-01 0.6680E-01 0.6681E-01  6681652  6679768              
!  2.00 0.2277E-01 0.2273E-01 0.2275E-01  2276651  2273197              
!  2.50 0.6222E-02 0.6206E-02 0.6210E-02   622190   620610              
!  3.00 0.1356E-02 0.1347E-02 0.1350E-02   135628   134674              
!  4.00 0.3242E-04 0.3145E-04 0.3170E-04     3242     3145              
!                                                                       
!--------------------------------------                                 
      DIMENSION RanNum (2)
      If (iFlag.eq.0) Then 
    1    CALL ranlux (RanNum, 2 ) 
         x1 = 2. * RanNum (1) - 1. 
         x2 = 2. * RanNum (2) - 1. 
         R = x1**2 + x2**2 
         If (R.ge.1.) GoTo 1 
         F              = sqrt ( - 2. * LOG (R) / R) 
         gSet        = x1 * F 
         GAUSS3 = x2 * F 
         iFlag        = 1 
      Else 
         GAUSS3 = gSet 
         iFlag        = 0 
      EndIf 
      RETURN 
      END FUNCTION GAUSS3                           

 !--------------------------------------                                 
!		                                               
      FUNCTION GAUSS3fake (gSet,iFlag) 
!
!--------------------------------------                                 
      DIMENSION RanNum (2)
      If (iFlag.eq.0) Then 
    1    CALL ranlux (RanNum, 2 ) 
         x1 = 2. * RanNum (1) - 1. 
         x2 = 2. * RanNum (2) - 1. 
         R = x1**2 + x2**2 
         If (R.ge.1.) GoTo 1 
         iFlag        = 1 
      Else 
         iFlag        = 0 
      EndIf 
      GAUSS3fake =0.
      RETURN 
      END FUNCTION GAUSS3fake                           

!---------------------------------------
!                   generate a vector of random numbers    
  SUBROUTINE getRandom(Gg,jp,kp,N)
!
!---------------------------------------
  use LUXURY
  Integer*8 :: kp,jp
  Real*4 :: Gg(N)

    Ns =SeedsPage(jp,kp) 
    lux = 2
    Call rluxgo(lux,Ns,0,0)
    call ranlux(Gg,N)
  end SUBROUTINE getRandom
!
                                                                        
!--------------------------------------------------------------------   
                                                                        
!   LUXURY LEVELS.                                                      
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia      
!           and Zaman, very long period, but fails many tests.          
!  level 1  (p=48): considerable improvement in quality over level 0,   
!           now passes the gap test, but still fails spectral test.     
!  level 2  (p=97): passes all known tests, but theoretically still     
!           defective.                                                  
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible         
!           correlations have very small chance of being observed.      
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.      
!---------------------------------------------                          
      SUBROUTINE ranlux (rvec, lenv) 
!                    returns a vector RVEC of LEN                       
!                   32-bit random floating point numbers between        
!                   zero (not included) and one (also not incl.).       
!           Next call to ranlux gives next LEN of random numbers        
!---------------------------------------------                          
      INTEGER lenv 
      REAL rvec (lenv)
      INTEGER iseeds (24) 

INCLUDE 'luxuryp.h' 
!      DATA ndskip / 0, 24, 73, 199, 365 / 
!      DATA i24, j24, luxlev / 24, 10, lxdflt / 
!      DATA notyet / .true. / 
!      DATA in24, kount, mkount, carry / 0, 0, 0, 0. / 
                                                                        
                                                                        
!  NOTYET is .TRUE. if no initialization has been performed yet.        
!              Default Initialization by Multiplicative Congruential    
                                                                        
      IF (notyet) THEN 
         stop ' Access to runlux before initialization STOP'
      ENDIF 
                                                                        
!          The Generator proper: "Subtract-with-borrow",                
!          as proposed by Marsaglia and Zaman,                          
!          Florida State University, March, 1989                        
                                                                        
      DO ivec = 1, lenv 
      uni = seeds (j24) - seeds (i24) - carry 
      IF (uni.LT.0.) THEN 
         uni = uni + 1.0 
         carry = twom24 
      ELSE 
         carry = 0. 
      ENDIF 
      seeds (i24) = uni 
      i24 = next (i24) 
      j24 = next (j24) 
      rvec (ivec) = uni 
      !  small numbers (with less than 12 "significant" bits) are "padde
      IF (uni.LT.twom12) THEN 
         rvec (ivec) = rvec (ivec) + twom24 * seeds (j24) 
      !        and zero is forbidden in case someone takes a logarithm  
         IF (rvec (ivec) .EQ.0.) rvec (ivec) = twom24 * twom24 
      ENDIF 
      !        Skipping to luxury.  As proposed by Martin Luscher.      
      in24 = in24 + 1 
      IF (in24.EQ.24) THEN 
         in24 = 0 
         kount = kount + nskip 
         DO isk = 1, nskip 
         uni = seeds (j24) - seeds (i24) - carry 
         IF (uni.LT.0.) THEN 
            uni = uni + 1.0 
            carry = twom24 
         ELSE 
            carry = 0. 
         ENDIF 
         seeds (i24) = uni 
         i24 = next (i24) 
         j24 = next (j24) 
         ENDDO 
      ENDIF 
      ENDDO 
      kount = kount + lenv 
      IF (kount.GE.igiga) THEN 
         mkount = mkount + 1 
         kount = kount - igiga 
      ENDIF 
      RETURN 
                                                                        
          ! SUBROUTINE ranlux                                           
      END SUBROUTINE ranlux                         
                                                                        
!---------------------------------------------                          
!                    Subroutine to initialize from one or three integers
      SUBROUTINE rluxgo (lux, ins, k1, k2) 
!                initializes the generator from                         
!               one 32-bit integer INT and sets Luxury Level LUX        
!               which is integer between zero and MAXLEV, or if         
!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2         
!               should be set to zero unless restarting at a break      
!               point given by output of RLUXAT                         
!---------------------------------------------                          
INCLUDE 'luxuryp.h' 
      INTEGER iseeds(24)
                                                                        
      IF (lux.LT.0) THEN 
         luxlev = lxdflt 
      ELSEIF (lux.LE.maxlev) THEN 
         luxlev = lux 
      ELSEIF (lux.LT.24.OR.lux.GT.2000) THEN 
         luxlev = maxlev 
         Write (6, '(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ', lux 
      ELSE 
         luxlev = lux 
         DO ilx = 0, maxlev 
         IF (lux.EQ.ndskip (ilx) + 24) luxlev = ilx 
         ENDDO 
      ENDIF 
      IF (luxlev.LE.maxlev) THEN 
         nskip = ndskip (luxlev) 
      !Write (6, '(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :', &
      !                 luxlev, '     P=', nskip + 24                                      
      ELSE 
         nskip = luxlev - 24 
      !   Write (6, '(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',       &
      !   luxlev                                                         
      ENDIF 
      in24 = 0 
      IF (ins.LT.0) Write (6,'(A)')' Illegal initialization by RLUXGO ',&
                               'negative input seed'                                             
      IF (ins.GT.0) THEN 
         jseed = ins 
      !Write (6, '(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS', &
      !                                  jseed, k1, k2                                                     
      ELSE 
         jseed = jsdflt 
      Write (6, '(A)') ' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED' 
      ENDIF 
      inseed = jseed 
      notyet = .false. 
      twom24 = 1. 
      DO i = 1, 24 
      twom24 = twom24 * 0.5 
      k = jseed / 53668 
      jseed = 40014 * (jseed-k * 53668) - k * 12211 
      IF (jseed.LT.0) jseed = jseed+icons 
      iseeds (i) = MOD (jseed, itwo24) 
      ENDDO 
      twom12 = twom24 * 4096. 
      DO i = 1, 24 
      seeds (i) = REAL (iseeds (i) ) * twom24 
      next (i) = i - 1 
      ENDDO 
      next (1) = 24 
      i24 = 24 
      j24 = 10 
      carry = 0. 
      IF (seeds (24) .EQ.0.) carry = twom24 
      !        If restarting at a break point, skip K1 + IGIGA*K2       
      !        Note that this is the number of numbers delivered to     
      !        the user PLUS the number skipped (if luxury .GT. 0).     
      kount = k1 
      mkount = k2 
      IF (k1 + k2.NE.0) THEN 
         DO iouter = 1, k2 + 1 
         inner = igiga 
         IF (iouter.EQ.k2 + 1) inner = k1 
         DO isk = 1, inner 
         uni = seeds (j24) - seeds (i24) - carry 
         IF (uni.LT.0.) THEN 
            uni = uni + 1.0 
            carry = twom24 
         ELSE 
            carry = 0. 
         ENDIF 
         seeds (i24) = uni 
         i24 = next (i24) 
         j24 = next (j24) 
         ENDDO 
         ENDDO 
      !         Get the right value of IN24 by direct calculation       
         in24 = MOD (kount, nskip + 24) 
         IF (mkount.GT.0) THEN 
            izip = MOD (igiga, nskip + 24) 
            izip2 = mkount * izip + in24 
            in24 = MOD (izip2, nskip + 24) 
         ENDIF 
      !       Now IN24 had better be between zero and 23 inclusive      
         IF (in24.GT.23) THEN 
      Write (6, '(A/A,3I11,A,I5)') '  Error in RESTARTING with RLUXGO:',&
           '  The values', ins, k1, k2, ' cannot occur at luxury level', luxlev 
            in24 = 0 
         ENDIF 
      ENDIF 
      RETURN 
                                                                        
      END SUBROUTINE rluxgo                     
end Module Random

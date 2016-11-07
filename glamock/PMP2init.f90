!--------------------------------------------------
!
!
!
!
!
!--------------------------------------------------
Module Param
      Integer*4, parameter :: NtabM = 100000
      Real*4  :: Om,        &      ! cosmology
                 Omb,       &
                 OmL,       &
                 AEXPN,     &
                 Sn,        &      ! normalization
                 sigma8,    &
                 hubble
      Integer*4 :: NROW,    &      ! number of particles in 1D
                   NGRID,   &      ! number of grid points in 1D
                   Nout,    &      ! number of outputs/mocks
                   Nbiaspars       ! number of bias parameters
      Real*4   :: Box, zinit,da,zfinal
      Real*4   :: densThr,sigV
      Real*4   :: zout(1000),BiasPars(1000)
      Real*8   :: rf,OmLOm0
      Real*4   :: xkt(0:NtabM),Pkt(0:NtabM)  ! Power spectrum
      Real*4   :: StepK,alog0
      INteger*4:: Ntab

    Contains
!---------------------------------------
REAL*8 FUNCTION P(x)
!                     interpolate table with p(k)
!                       x is in real 1/Mpc
!---------------------------------------
real*8 :: x,dk
  If(x.ge.xkt(Ntab))Then  ! slope is ns =-3
         P =Pkt(Ntab)/(x/xkt(Ntab))**3
         Return
  EndIf
  If(x.lt.xkt(1))Then                ! slope is ns=1
         P =Pkt(1)*(x/xkt(1))
         Return
  EndIf
     ind = INT((log10(x)-alog0)/StepK) +1 
      dk  = xkt(ind+1)-xkt(ind)
      P   = (Pkt(ind)*(xkt(ind+1)-x)+Pkt(ind+1)*(x-xkt(ind)))/dk
      Return
End FUNCTION P

!-------------------------------------- P*k^2*Top-Hat Filter
REAL*8 FUNCTION Ptophat(wk)
!---------------------------------------
real*8 :: wk,X,TOPHAT
        IF (wk.lt.1.d-4) THEN
            Ptophat =wk*wk**2
        ELSE
            X      =wk*rf
	    TOPHAT =( (SIN(X)-x*COS(X))*3./X**3 )**2
            Ptophat=P(wk)*wk**2*TOPHAT
        ENDIF
      END FUNCTION Ptophat
!
!-----------------------------------  Pcold(k)*k^2
REAL*8 FUNCTION P2(WK)
!
  Real*8 :: WK  
    P2=WK**2*P(WK)
END FUNCTION P2
!---------------------------------------
REAL*8 FUNCTION Hnorm(x)
!
   Real*8 :: x
         Hnorm =sqrt(1.d+0 +OmLOm0*x**3)
End FUNCTION Hnorm
!---------------------------------------
      REAL*8 FUNCTION Hage(x)
!
    Real*8 :: x
         Hage =sqrt(x)/Hnorm(x) 
End FUNCTION Hage
!---------------------------------------
      REAL*8 FUNCTION Hgrow(x)
!
    Real*8 :: x
        Hgrow =(sqrt(x)/Hnorm(x))**3 
End FUNCTION Hgrow

!-------------------------------------------------
!                       Age of the Universe: t0 (z=0)
!                       Expansion parameter: a =1/(1+z)
!                       Growth_Rate_Density at a: GrowthDen
!                          GrowthDen =\delta\rho/\rho
!                          normalized to GrowthDen =a for Omega=1
!                          GrowthDen < 1 at a=1 
!                     
      SUBROUTINE AGE(t0,GrowthDen,a)
!
!---------------------------------------
      real*8 :: t0,GrowthDen,a,ww        
      Real*8, PARAMETER :: zero =1.d-12
      Real*8  INTG

      OmLOm0 =OmL/Om
      t0     =9.766/hubble/sqrt(Om)*INTG(Hage,zero,a)
      !Hubble = hubble*sqrt(Om/a**3)*Hnorm(a)
      ww     = INTG(Hgrow,zero,a)
      GrowthDen =2.5*Hnorm(a)/sqrt(a**3)*ww
      
End SUBROUTINE AGE

End Module Param
!
!--------------------------------------------------
Program Initialize
  use Param
      CHARACTER  ::    Name*120,Header*45
      REAL*8 ::        INTG,wk,Uklow,Ukup,Sig8,sigma,a,t0,Growthden,Growthden_0
      Real*8, parameter :: PI     =3.1415926535d0
      logical    :: exst
      EXTERNAL         INTG

      Call CheckInit !---- check whether  all input files exist
      Call ReadInit  !---- read all input information
      
      OPEN(10,file='Setup.dat')
      OPEN(1,file='lcdm.dat')
      
      H=100.*hubble     ! Hubble constant
      a =1.
      CALL AGE(t0,GrowthDen_0,a)
      write (*,'(3(a,es12.4))')' T_universe= ',t0,' a=',a,' GrowthRate= ',GrowthDen_0
      write (1,'(3(a,es12.4))')' T_universe= ',t0,' a=',a,' GrowthRate= ',GrowthDen_0

      rf    =8.        ! top-hat for sigma_8
      Sn    = (Sigma8)**2 / &
               INTG(Ptophat,1.d-5,5.d+1) ! normalization of P(k)
                                   !          bias parameter      
      Sig8 = sqrt(Sn*(INTG(Ptophat,1.d-5,5.d+1)))
      write (*,10) hubble,Sigma8  ! check if all is self-consistent
      write (1,10) hubble,Sigma8  ! check if all is self-consistent
      If(abs(Sig8-sigma8)>0.005*sigma8)write(*,'(a,2es14.4)') &
         ' Difference betweeen requested and real Sigma8 is:',Sig8,sigma8
 10   Format(' Hubble=',f7.3,' Sigma_8 =',G11.3)

      Uklow =2.*pi/Box            ! frequency range for integrals
      Ukup  =2.*pi/Box*(NROW/2.)
      AEXPN =1./(1.+zinit)            ! expansion parameter
      a     =AEXPN
      CALL AGE(t0,GrowthDen,a)
      sigma =(GrowthDen/GrowthDen_0)* &
                sqrt(Sn*INTG(P2,Uklow/sqrt(2.),Ukup))
      WRITE (*,40) zinit,sigma,NROW,Box
      WRITE (1,40) zinit,sigma,NROW,Box
 40   format('  z=',f8.3,' delta\rho/rho in box=',f9.5,/ &
            5X,'Particles=',i4,' Box=',f8.2)
      write (1,*) ' k/h        Pk*h^3   Power Spectrum at z=',1/AEXPN-1.
      Do i=1,1000
         wk = 1.d-3*10.**((i-1)/20.)
         if(wk> 100.d0)exit
         pp = P(wk)
         write(1,'(2es14.5)') wk,pp
      end Do

      write (10,*) '--------------- Setup file for PMP2 simulations ------------------'
      write(HEADER,'(3(a,i4.4),a,i3.3,a,es8.2,a,f5.3)') &
       'N=',NROW,'x',NGRID,'L=',Int(Box),'zi',INT(zinit),'da=',da,'f',da/AEXPN
      write(10,'(a)')TRIM(HEADER)
      write (10,50) AEXPN,     'Expansion Parameter'
      write (10,'(es12.5,T20,a,es12.4)') da,         'Step in dAEXPN    da/a = ',da/AEXPN
      write (10,50) sigma,     'DRho/rho in box    '
      write (10,50) Box,       'Box in  Mpc/h   '
      write (10,50) sigma8,    'sigma8    '
      write (10,50) hubble,    'Hubble    '
      write (10,50) Om,        'Omega Matter'
      write (10,50) OmL,       'Omega Lambda'
      write (10,50) Omb,       'Omega Baryons'
      write (10,60) NROW,      'NROW  Number Particles   '
      write (10,60) NGRID,     'NGRID Number grid points '
      write (10,60) 0,         'Random seed'
      write (10,50) Box/Ngrid, 'Cell Size   '
      write (10,50) 2.5841e11*Om*hubble**2*(Box/1000./(NROW/1024.))**3, 'Particle Mass'   
      write (10,50) zinit,      'Initial redshift       '
      write (10,50) zfinal,     'Final redshift       '
      write (10,50) DensThr,    'Density Threshold for V correction '
      write (10,50) sigV,       'rms V correction factor'
      write (10,60) Nout,       'Number of redshifts for analysis'
      write (10,'(100f8.3)')(zout(i),i=1,Nout)
      write (10,60) Nbiaspars,  'Number of bias parameters'
      Do i=1,Nbiaspars
         write (10,50) BiasPars(i), 'Bias'
      End Do     
50    format(es12.5,T20,a)
60    format(i5,T20,a)
         write (*,*) ' Results were written to Setup.dat'
      CLOSE (10)

      a_init  = AEXPN
      da_init = da
      Call TestStepping(a_init,da_init)
      
    end Program Initialize
!
!---------------------------------------------------             
!
   Subroutine ReadInit
!---------------------------------------------------
     use Param
            !-- Read PkTable. Assign Omegas and hubble
      Omb = ParseLine(10)
      Omc = ParseLine(10)
      OmL = ParseLine(10)
      Om  = ParseLine(10)
      sigma8 = ParseLine(10)
      hubble =0.678
      write(*,'(3(a,ES12.3))')' Model from PkTable file: Ombar =', &
                       Omb,' Om_matter =',Om,' Sigma8 =',sigma8
      Ntab = 0
12   READ (10, *, end = 32, err = 32) xx, pp 
        Ntab = Ntab + 1 
        xkt (Ntab) = xx  !*hubble 
        Pkt (Ntab) = pp                !  Pk 
      GOTO 12 
32    Write (*,*) ' Read ', Ntab, ' lines from P(k) table '
          close(10)
      If (Ntab.le.1) stop 'wrong table for p(k)' 
      StepK = log10 (xkt (Ntab) / xkt (1) )/(Ntab-1) 
      alog0 = log10 (xkt (1) ) 
                                                                        
                            ! test that spacing is the same             
      Do k = 2, Ntab - 1 
        ss = log10 (xkt (k + 1) / xkt (k) ) 
        If (abs (ss / StepK - 1.) .gt.2.e-2) Then 
        Write (*,*) ' error in K spacing. k=', k, xkt (k + 1),xkt (k)
        STOP 
        EndIf 
      EndDo 

      !------------- Read input parameters from Init.dat
      Box    = ParseLine(11)
      Nrow   = iParseLine(11)
      Ngrid  = iParseLine(11)
      sigma8 = ParseLine(11)
      zinit  = ParseLine(11)
      da     = ParseLine(11)
      zfinal = ParseLine(11)
      Nout   = iParseLine(11)
      read(11,*)(zout(i),i=1,Nout)
      densThr = ParseLine(11)
      sigV    = ParseLine(11)
      Nbiaspars = iParseLine(11)
      Do i=1,Nbiaspars
         BiasPars(i) = ParseLine(11)
      End Do
                   !--- make new da
      fr = da*(1.+zinit)*100.     ! = da/a*100
      frnew = INT(fr*10.)/10./100.
      da  = frnew/(1.+zinit)
      write(*,*) ' Done reading input from Init.dat'
      close(11)
      write(*,*) 'Box   = ',Box
      write(*,*) 'Nrow  =',Nrow
      write(*,*) 'Ngrid =',Ngrid
      write(*,*) 'Nout  =',Nout
      write(*,*) 'Nbias =',Nbiaspars
    end Subroutine ReadInit
!
!---------------------------------------------------
!
!             check if all init files are present
!
   Subroutine CheckInit
!---------------------------------------------------
     use Param
     logical :: exst
      Inquire(file='PkTable.dat',exist = exst)
      if(.not.exst)Stop ' File PkTable with the power spectrum not found'
      open(10,file='PkTable.dat')
      
      Inquire(file='Init.dat',exist = exst)
      if(.not.exst)Then
         write(*,*)' File Init.dat with initial parameters not found.'
         write(*,*)' I create a new one. Edit it and restart the code'
         open(11,file='Init.dat')
         write(11,'(a,f9.3)')   'Box      = ',1000.
         write(11,'(a,i9)')     'Nrow     = ',1000
         write(11,'(a,i9)')     'Ngrid    = ',2000
         write(11,'(a,f9.3)')   'sig8     = ',0.828
         write(11,'(a,f9.3)')   'z_init   = ',100.
         write(11,'(a,es11.4)') 'step da  = ',4e-4
         write(11,'(a,f9.3)')   'z_final  = ',0.
         write(11,'(a,i9)')     '#outputs = ',10
         write(11,'(20f5.2)')  2.5,1.5,1.,0.8,0.7,0.5,0.3,0.2,0.1,0.
         write(11,'(a,f9.3)')   'dens_thr = ',30.
         write(11,'(a,f9.3)')   'Vrms     = ',sigv
         write(11,'(a,i9)')     '#Params  = ',10
         do i=1,10
            write(11,'(a,f9.3)') 'Parametr = ',0.
         endDo
         stop
      end if
         open(11,file='Init.dat')

      Inquire(file='TableSeeds.dat',exist = exst)
      if(.not.exst)Call SetSeeds
 
end Subroutine CheckInit
!------------------------------
!
!             Generate a table with random seeds
!
!------------------------------
Subroutine SetSeeds
  integer*8 :: is,ij,iM
  integer*4 :: Nseed
  is = 1232_8**3
  
    Nseed  = 1298302
    nslip = 137
    Noff  = 2357
    Ntable =5000
    NCount =0
    open(1,file='TableSeeds.dat')
    write(1,*) 'Seeds:',Nseed,Nslip
   Do ij=1,is
      x  =RANDd(Nseed)
      Nn =INT(x*nslip)+1
      Do jj =1,Noff+Nn
         x  =RANDd(Nseed)
      End Do
      Ncount =Ncount +1
      write(1,*) Nseed,Ncount
      If(Ncount>Ntable)exit
   endDo
   close(1)
 end Subroutine SetSeeds

   
!------------------------------------------------                       
!				          random number generator     
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

 
!--------------------------------------------------
!        read line from  input file iFile
!                real format
      Function ParseLine(iFile)
      Character :: Line*120,Line2*120,Line3(120)

      Read(iFile,'(a)')Line
      Ieq =INDEX(Line,'=',BACK=.TRUE.)
               !write(*,*) '  Ieq =',Ieq
      backspace (iFile)                  !--- go to line start
      write(Line2,'(a1,i2,a)') '(',Ieq,'a1,g12.5)' ! make format
           !write(*,'(a)') Line2
      Read(iFile,Line2)(Line3(i),i=1,Ieq),dummy    ! read
      ParseLine = dummy
           !write(*,'(a,ES12.3)') ' Result =',ParseLine
      end Function ParseLine
!--------------------------------------------------
!        read line from  input file iFile
!                          integer format
      Function iParseLine(iFile)
      Character :: Line*120,Line2*120,Line3(120)

      Read(iFile,'(a)')Line
      Ieq =INDEX(Line,'=',BACK=.TRUE.)
      !write(*,*) '  Ieq =',Ieq
           backspace (iFile)
           write(Line2,'(a1,i2,a)') '(',Ieq,'a1,i10)'
           !write(*,'(a)') Line2
      Read(iFile,Line2)(Line3(i),i=1,Ieq),idummy
      iParseLine = idummy
      !write(*,'(a,ES12.3)') ' Result =',ParseLine
      end Function iParseLine
                         

!-------------------------------------------- Simpson integration
      REAL*8 FUNCTION INTG(FUNC,A,B)
!---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (EPS=3.0d-5, JMAX=22)
      EXTERNAL FUNC
      OST=-1.D30
      OS= -1.D30
      ST =0.
      DO  J=1,JMAX
        CALL TRAPZD(FUNC,A,B,ST,J)
        INTG=(4.0d0*ST-OST)/3.0d0
        IF (ABS(INTG-OS).Le.EPS*ABS(OS)) RETURN
        OS=INTG
        OST=ST
     end DO
      WRITE (*,*)'Integration did not converge'
    END FUNCTION INTG
!----------------------------------------------
      SUBROUTINE TRAPZD(FUNCC,A,B,S,N)
!---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
        SAVE IT
        EXTERNAL FUNCC
      IF (N.EQ.1) THEN
        S=0.5d0*(B-A)*(FUNCC(A)+FUNCC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5D0*DEL
        SUM=0.0D0
        DO  J=1,IT
          SUM=SUM+FUNCC(X)
          X=X+DEL
        end DO
        S=0.5D0*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
    END SUBROUTINE TRAPZD
!----------------------------------------------
SUBROUTINE TestStepping(a_init,da_init)
!----------------------------------------------
     StepFactor = da_init/a_init
     a = a_init
         da= da_init
         i = 0
         write(*,'(2(a,f9.4))') ' a_init =',a, ' z_init=',1./a-1.
         write(*,'(2(a,f9.4))') ' da     =',da,' da/a  =',da/a
         write(*,'(2(a,f9.4))') ' stFact =',stepFactor
         write(*,'(a)') ' Step      a         da        da/a    step change'   
         iCount = 0
         rMax   = da/a
         Nsteps2= 0
         Do 
            ind = 0
            If(da<StepFactor/1.25*a.and.a<0.333)Then
               write(*,*) '                  Increase step:',da,1.25*StepFactor*a
               da =1.5*da        ! increase step
               ind =1
               iCount = iCount +1
            EndIf
            a = a +da
            i = i +1
            rmax = max(rMax,da/a)
            if(a>0.333)Nsteps2 = Nsteps2 +1
            write(*,'(i5,3es12.4,i3)') i,a,da,da/a,ind
            if(a>1.)exit
         end Do
         write(*,'(2(a,i5))')   ' Number of steps   = ',i, ' n_steps(z<2) =',Nsteps2 
         write(*,'(a,i5)')   ' Number of changes = ',iCount 
         write(*,'(a,f8.4)') ' Maximum da/a      = ',rmax
end SUBROUTINE TestStepping


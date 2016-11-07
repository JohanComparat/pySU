!-------------------------------------------------
!
!    Routines to estimate Power Spectrum
!
!-------------------------------------------------
Module Power
  use Tools
  use fft5
  Real*4, Allocatable, Dimension(:,:)   :: Pk,dNharm,kharm
  Integer*4                             :: iCIC=1
Contains
!----------------------------------------------------
!           power spectrum   P(k)
!
!               iFlag =0 - only real space
!                     =1 - real + 3 projections P0 + 3 projections P2
!
SUBROUTINE GetPower(iFlag)
!
!----------------------------------------------------
  use Density
  Integer*4   :: moment
  Real*4      :: kNyq
  Character*120 :: Name
      Mem_current = Memory(INT(21*NGRID,8))
    Allocate(Pk(NGRID,0:6),dNharm(NGRID,0:6),kharm(NGRID,0:6))
    moment = INT(100.*(1./AEXPN-1.)+0.5)  ! make integer out of Z_moment
    Pk(:,:)     = 0.
    dNharm(:,:) = 0.
    kharm(:,:)  = 0.
    boxsize = Box
    kNyq = 2.*3.1415926/Box*NGRID/2.
    If(iFlag==0)Then              ! --- only real space Pk
       write(Name,'(2(a,i4.4),3(a,i3.3))')'PowerDM.R.',moment,   &
                                        '.',Nrealization,'.dat'
       OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(18,'(a,f7.4,a,i4,a,f8.3,a,i4,a,f8.2)')' Aexpn =',AEXPN,' Step=',ISTEP,    &
            ' Redshift= ',1./AEXPN-1.,' Ngrid= ',Ngrid,' kNyq= ',kNyq
       CALL POWERfft5(0,0)
       WRITE (18,'(3x,a,T9,a,T19,a,T30,a,T43,a)') 'bin','k/Mpch',      &
                 'N','Preal'
       DO I=1,NGRID/2
         WRITE (18,'(i6,f9.5,es11.4,7es13.5)')I,kharm(I,0),dNharm(I,0),Pk(I,0)
      EndDo
      
   Else                         ! --- Real + redshift space
      !--------------  change to smaller grid

      Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
      Deallocate (FI)
      NGRID_old = NGRID                  ! store old value of NGRID
       NGRID = NGRID/2
      Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
      ALLOCATE(FI(NGRID,NGRID,NGRID))
      write(Name,'(2(a,i4.4),3(a,i3.3))')'PowerDM.Z.',moment,'.',Nrealization,   &
           '.V',INT(sigv),'.d',INT(DensThr),'.dat'
      OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(18,'(a,f7.4,a,i4,a,f8.3,a,i4,a,f8.2)')' Aexpn =',AEXPN,' Step=',ISTEP,    &
            ' Redshift= ',1./AEXPN-1.,' Ngrid= ',Ngrid,' kNyq= ',kNyq

        Do iSwitch = 0, 3
          Call DENSITrsd(iSwitch,NGRID_old)    ! find density for all particles in PMfiles
          CALL POWERfft5(0,iSwitch)               ! Power spectrum monopole
          If(iSwitch/=0)Then
            Call DENSITrsd(iSwitch,NGRID_old)    ! find density for all particles in PMfiles
            CALL POWERfft5(1,iSwitch) ! Power spectrum quadrupole
         End If
        EndDo

   WRITE (18,'(3x,a,T9,a,T19,a,T30,a,T43,a,T56,a)') 'bin','k/Mpch',      &
                 'N','Preal','P0_redshift','P2_redshift'
   
   DO I=1,NGRID/2
      P2 = (Pk(I,1)+ Pk(I,3)+Pk(I,5))/3.
      P0 = (Pk(I,2)+ Pk(I,4)+Pk(I,6))/3.
      !P2 = (Pk(I,1)+ Pk(I,3))/2.
      !P0 = (Pk(I,2)+ Pk(I,4))/2.
      WRITE (18,110)  I,kharm(I,0),dNharm(I,0),Pk(I,0),P0,P2 !(Pk(I,j),j=0,6)
   EndDo
   110 format(i6,f9.5,es11.4,7es13.5)
      Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
      Deallocate (FI)
      NGRID =  NGRID_old                  ! restore the old value of NGRID
      Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
      ALLOCATE(FI(NGRID,NGRID,NGRID))     
 End If
 
       Mem_current = Memory(-INT(21*NGRID,8))
    DeAllocate(Pk,dNharm,kharm)
    close(18)
    
      end SUBROUTINE GetPower
!----------------------------------------------------
!           power spectrum   P(k) for tracers
!
!               iFlag =0 - only real space
!                     =1 - real + 3 projections P0 + 3 projections P2
!
SUBROUTINE GetPowerGal(iFlag)
!
!----------------------------------------------------
  use Density
  Integer*4   :: moment
  Real*4      :: kNyq
  Character*120 :: Name
      Mem_current = Memory(INT(21*NGRID,8))
    Allocate(Pk(NGRID,0:6),dNharm(NGRID,0:6),kharm(NGRID,0:6))
    moment = INT(100.*(1./AEXPN-1.)+0.5)  ! make integer out of Z_moment
    Pk(:,:)     = 0.
    dNharm(:,:) = 0.
    kharm(:,:)  = 0.
    boxsize = Box
    kNyq = 2.*3.1415926/Box*NGRID/2.
    If(iFlag==0)Then              ! --- only real space Pk
       write(Name,'(2(a,i4.4),3(a,i3.3))')'PowerGal.R.',moment,   &
                                        '.',Nrealization,'.dat'
       OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(18,'(a,f7.4,a,i4,a,f8.3,a,i4,a,f8.2,a,i10)')' Aexpn =',AEXPN,' Step=',ISTEP,    &
            ' Redshift= ',1./AEXPN-1.,' Ngrid= ',Ngrid,' kNyq= ',kNyq
       write(18,'(a,i10,a,f8.3,a,es13.4,a,f8.3)'),' Ngalaxies=',Ngalaxies, &
            ' DensThresh= ',BiasPars(1),    &
            ' Normalize = ',BiasPars(2),   &
            ' Slope     = ',BiasPars(3)
       CALL POWERfft5(0,0)
       WRITE (18,'(3x,a,T9,a,T19,a,T30,a,T43,a)') 'bin','k/Mpch',      &
                 'N','Preal'
       DO I=1,NGRID/2
         WRITE (18,'(i6,f9.5,es11.4,7es13.5)')I,kharm(I,0),dNharm(I,0),Pk(I,0)
      EndDo
      
   Else                         ! --- Real + redshift space
      !--------------  change to smaller grid

      Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
      Deallocate (FI)
      NGRID_old = NGRID                  ! store old value of NGRID
       NGRID = NGRID/2
      Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
      ALLOCATE(FI(NGRID,NGRID,NGRID))
      write(Name,'(2(a,i4.4),3(a,i3.3))')'PowerGal.Z.',moment,'.',Nrealization,   &
           '.dat'
      OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(18,'(a,f7.4,a,i4,a,f8.3,a,i4,a,f8.2,a,i10)')' Aexpn =',AEXPN,' Step=',ISTEP,    &
            ' Redshift= ',1./AEXPN-1.,' Ngrid= ',Ngrid,' kNyq= ',kNyq
       write(18,'(a,i10,a,f8.3,a,es13.4,a,f8.3)'),' Ngalaxies=',Ngalaxies, &
            ' DensThresh= ',BiasPars(1),    &
            ' Normalize = ',BiasPars(2),   &
            ' Slope     = ',BiasPars(3)
        Do iSwitch = 0, 3
          Call DensGal(iSwitch)    ! find density for all particles in PMfiles
          CALL POWERfft5(0,iSwitch)               ! Power spectrum monopole
          If(iSwitch/=0)Then
            Call DensGal(iSwitch)    ! find density for all particles in PMfiles
            CALL POWERfft5(1,iSwitch) ! Power spectrum quadrupole
         End If
        EndDo

   WRITE (18,'(3x,a,T9,a,T19,a,T30,a,T43,a,T56,a)') 'bin','k/Mpch',      &
                 'N','Preal','P0_redshift','P2_redshift'
   
   DO I=1,NGRID/2
      P2 = (Pk(I,1)+ Pk(I,3)+Pk(I,5))/3.
      P0 = (Pk(I,2)+ Pk(I,4)+Pk(I,6))/3.
      !P2 = (Pk(I,1)+ Pk(I,3))/2.
      !P0 = (Pk(I,2)+ Pk(I,4))/2.
      WRITE (18,110)  I,kharm(I,0),dNharm(I,0),Pk(I,0),P0,P2 !(Pk(I,j),j=0,6)
   EndDo
   110 format(i6,f9.5,es11.4,7es13.5)
      Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
      Deallocate (FI)
      NGRID =  NGRID_old                  ! restore the old value of NGRID
      Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
      ALLOCATE(FI(NGRID,NGRID,NGRID))     
 End If
 
       Mem_current = Memory(-INT(21*NGRID,8))
    DeAllocate(Pk,dNharm,kharm)
    close(18)
    
  end SUBROUTINE GetPowerGal
!----------------------------------------------------
!                                 power spectrum   P(k)
!                                 for given density field FI
SUBROUTINE POWERfft5(iPol,iSwitch)
!----------------------------------------------------
      PARAMETER ( Pi=3.141592653)
      Real*8,    Allocatable,DIMENSION(:) :: DPOWER,dk
      Real*8, Allocatable,DIMENSION(:)    :: iVolume
      Real*8,    Allocatable,DIMENSION(:,:) :: DPW,dkk
      Real*8, Allocatable,DIMENSION(:,:)  :: iVol
      integer*4, parameter :: Nlensav = 16384
      integer*4, parameter :: Nlenwrk = 16384
      real*8,SAVE          :: work(1:Nlenwrk)
      real*8,SAVE          :: wsave(1:Nlensav)
      real*8               :: r(Nlenwrk)
      REAL*8               :: SIGMA,k_Ny,PkCorrec,Pkk,dens
      Integer*4            :: OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM
!$OMP THREADPRIVATE(work,wsave)
      Call Timing(4,-1)
      lensav = Ngrid+int(log(real(Ngrid,kind = 4))/log(2.0E+00))+4
      lenwrk = Ngrid
      NPOWER =NGRID
      boxsize = Box
      iThreads = OMP_GET_MAX_THREADS()  
      ALLOCATE( DPOWER(NPOWER),iVolume(NPOWER),dk(NPOWER))
      ALLOCATE( DPW(NPOWER,iThreads),iVol(NPOWER,iThreads), &
               dkk(NPOWER,iThreads))

      k_Ny = 2.*Pi/boxsize*NGrid/2.
      dens = Nparticles/boxsize**3
      NSPEC = NGRID/2
      If(Boxsize.le.0..or.Boxsize.gt.1.e+4)Then
         write (*,*) ' Wrong Box size=',Boxsize,' STOP'
         Return
      EndIf
      SIGMA =0.
      Scale = 2.*Pi/Boxsize
      DO M =1,NGRID
         iVolume(M) =0
         DPOWER(M)  =0.
         dk(M)      =0.
         Do i=1,iThreads
           iVol(M,i) =0
           DPW(M,i)     =0.
           dkk(M,i)     =0.
         EndDo 
      ENDDO

      call rfft1i ( Ngrid, wsave, lensav, ier ) !   Initialize FFT
      inc  = 1
      lenr = Ngrid

       write(*,*) ' time =',seconds()/60., ' go to first fft'
!$OMP PARALLEL DO DEFAULT(SHARED) copyin(wsave,work)  & 
!$OMP PRIVATE ( k,j,i ,r, ier)
       Do k=1,NGRID     ! ------ x-direction
       Do j=1,NGRID     
          Do i=1,NGRID
             r(i) = FI(i,j,k)
          EndDo
          call rfft1f ( Ngrid, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do i=1,NGRID
            FI(i,j,k) = r(i)
          EndDo
       EndDo
       EndDo
       Do k=1,NGRID
       Do j=1,NGRID
          FI(1,j,k) = FI(1,j,k)*sqrt(2.)
          FI(NGRID,j,k) = FI(NGRID,j,k)*sqrt(2.)
       end Do
       end Do
       write(*,*) ' time =',seconds()/60., ' go to second fft'
!$OMP PARALLEL DO DEFAULT(SHARED) copyin(wsave,work) & 
!$OMP PRIVATE (k, j,i ,r, ier )
       Do k=1,NGRID     ! ------ y-direction
       Do i=1,NGRID     
          Do j=1,NGRID
             r(j) = FI(i,j,k)
          EndDo
          call rfft1f ( Ngrid, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do j=1,NGRID
            FI(i,j,k) = r(j)
          EndDo
       EndDo
       EndDo
       Do k=1,NGRID
       Do i=1,NGRID
          FI(i,1,k) = FI(i,1,k)*sqrt(2.)
          FI(i,NGRID,k) = FI(i,NGRID,k)*sqrt(2.)
       end Do
       end Do
       write(*,*) ' time =',seconds()/60., ' go to transposition'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,Ngrid
      DO K=1,Ngrid-1
            DO I=K+1,Ngrid
               aa = FI(I,J,K)
               FI(I,J,K) =FI(K,J,I)
               FI(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO
       write(*,*) ' time =',seconds()/60., ' go to third fft'
!$OMP PARALLEL DO DEFAULT(SHARED) copyin(wsave,work)  & 
!$OMP PRIVATE ( k,j,i ,r, ier)
       Do j=1,NGRID     ! ------ z-direction
       Do i=1,NGRID     
          Do k=1,NGRID
             r(k) = FI(k,i,j)
          EndDo
          call rfft1f ( Ngrid, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do k=1,NGRID
            FI(k,i,j) = r(k)
          EndDo
       EndDo
       EndDo
       Do j=1,NGRID
       Do i=1,NGRID
          FI(1,i,j) = FI(1,i,j)*sqrt(2.)
          FI(NGRID,i,j) = FI(NGRID,i,j)*sqrt(2.)
       end Do
       end Do
       write(*,*) ' time =',seconds()/60.,'  Normalize '
!$OMP PARALLEL DO DEFAULT(SHARED)   & 
!$OMP PRIVATE ( k,j,i)
       Do k=1,NGRID
       Do j=1,NGRID
       Do i=1,NGRID
          FI(i,j,k) = FI(i,j,k)**2/8.
       end Do
       end Do
       end Do

      FI(1,1,1) =0.
!      kStep     =NGRID/iThreads
!      If(kstep*iThreads.ne.NGRID)STOP '-- Wrong number of iThreads'
      write (*,*) '    Power .... threads =',iThreads

!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE ( k,j,i,MI,MJ,MK,AMP,RK,MP,iOMP,MP2,d1,d2,w1,w2,cosmu2,P2) &
!$OMP REDUCTION(+:SIGMA)
      Do K=1,NGRID
         iOMP = OMP_GET_THREAD_NUM()+1
         !iOMP = 1
         MK = (K/2)**2
         !write(*,*) K,K/2,MK
      DO J=1,NGRID
         MJ = (J/2)**2
         DO I=1,NGRID
           MI = (I/2)**2
           AMP = FI(I,J,K)
           SIGMA = SIGMA +AMP
           RK =  MI +MJ+MK        ! k**2
!                                   Power Spectrum
           If(RK.lt.0.5)cycle     ! do not take k =0 component
           RK =Sqrt(RK)
           If(iCIC ==1)Then
              MP =max(INT(RK +0.0001),1)     ! left boundary
              MP2 = min(MP+1,NGRID)         ! right boundary
              d1  = MP
              d2  = d1 + 1.
              w1  =  MP - RK +1.
              w2  =  1. - w1              
              If(iPol==0)Then
                DPW(MP,iOMP)  =DPW(MP,iOMP) +AMP*w1
                DPW(MP2,iOMP)  =DPW(MP2,iOMP) +AMP*w2
              Else
                cosmu2 = MI/RK**2             ! square of cos()
                P2     = 1.5*cosmu2 -0.5   ! Lagrange Pol_2
                DPW(MP,iOMP)  =DPW(MP,iOMP) +AMP*w1*P2
                DPW(MP2,iOMP)  =DPW(MP2,iOMP) +AMP*w2*P2
             EndIf      ! iPol
           
              iVol(MP,iOMP) =iVol(MP,iOMP)+w1
              dkk(MP,iOMP)  =dkk(MP,iOMP) +d1*w1
              iVol(MP2,iOMP) =iVol(MP2,iOMP)+w2
              dkk(MP2,iOMP)  =dkk(MP2,iOMP) +d2*w2
            Else      ! NGP in k-space
               MP = max(INT(RK +0.5),1)     ! middle bin
               If(iPol==0)Then
                  DPW(MP,iOMP)  =DPW(MP,iOMP) +AMP
               Else
                  cosmu2 = MI/RK**2     ! square of cos()
                  P2     = 1.5*cosmu2 -0.5   ! Lagrange Pol_2
                  DPW(MP,iOMP)  =DPW(MP,iOMP) +AMP*P2                  
               EndIf
               iVol(MP,iOMP) =iVol(MP,iOMP)+1.
               dkk(MP,iOMP)  =dkk(MP,iOMP) +RK
            End If
            ENDDO
         ENDDO
      ENDDO
      Do MP=1,NGRID
      Do iP=1, iThreads
          DPOWER(MP)  =DPOWER(MP) +DPW(MP,iP)
          iVolume(MP) =iVolume(MP)+iVol(MP,iP)
          dk(MP)      =dk(MP)     +dkk(MP,iP)
      EndDo 
      EndDo
      If(iPol == 1)DPOWER(:) = 2.5*DPOWER(:) !! (2L+1)/2 for quadrupole
      WRITE (*,100) SQRT(SIGMA),2.*Pi*(float(Nparticles))**0.3333/boxsize/2., &
                                2.*Pi/boxsize*NGRID/2.
      !WRITE (17,100) SQRT(SIGMA),2.*Pi*(float(Nparticles))**0.3333/boxsize/2., &
      !                          2.*Pi/boxsize*NGRID/2.
100   FORMAT(8x,'RMS DRho/Rho            =',F7.3,   &
            /8x,'k_Nyquist for Particles =',f6.2, &
            /8x,'k_Nyquist for Grid      =',f6.2, &
            /4X,'n',2X,'k(h/Mpc)',1X, &
            'Harmonics ','Power Spectrum(h^3)  Pk_corrected')
110   FORMAT(I6,F11.6,1p,g12.4,3X,5G14.5)
            ii = 2*iSwitch-iPol
      !write(*,*) ' iPol/sw =',iPol,iSwitch,ii
      DO I=1,NGRID/2
         IF(iVolume(I).GT.0)THEN
            wk  = dk(I)/iVolume(I)*Scale
            Pkk  = Boxsize**3*DPOWER(I)/iVolume(I)
            PkCorrec = Pkk/(1.-0.6667*(sin(Pi*wk/k_Ny/2.))**2)-1./dens
            !if(I< 20)WRITE (*,110)  I,wk,iVolume(I),Pkk,PkCorrec
            !WRITE (17,110) I,wk,iVolume(I),Pkk,PkCorrec
            kharm(I,ii)  = wk
            dNharm(I,ii) = iVolume(I)
            Pk(I,ii)     = Pkk  !!!PkCorrec
         ENDIF
      ENDDO
      DEALLOCATE( DPOWER,iVolume,dk)
      DEALLOCATE( DPW,iVol,dkk)
      
      Call Timing(4,1)      
    END SUBROUTINE POWERfft5

end Module Power

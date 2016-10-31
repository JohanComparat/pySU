!---------------------------------------------------------
!
Module Density
!
!---------------------------------------------------------
!
  use Tools
  
Contains

!-------------------------------------------------------------
!
!         Make density contrast FI(NGRID,NGRID,NGRID)
!         for Nparticles XPAR,YPAR,ZPAR(Nparticles) in 1-Ngrid
!
SUBROUTINE DENSIT
    
!-------------------------------------------------------------
real*8 :: XN,YN, &
     X,Y,Z,D1,D2,D3,T1,T2,T3,T2W,D2W
integer*8 :: IN

    Call TimingMain(3,-1)
    XN   =FLOAT(NGRID)+1.-1.E-7
    YN   =FLOAT(NGRID)
    Wpar = YN**3/FLOAT(Nparticles)
				!       Subtract mean density
!write(*,*) ' Init density'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3)
    DO M3=1,NGRID
       DO M2=1,NGRID
          DO M1=1,NGRID
	        FI(M1,M2,M3) = -1.
	  END DO
       END DO
    END DO
    
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (IN,X,Y,Z,D1,D2,D3,T1,T2,T3,T2W,D2W) &
!$OMP PRIVATE (I,J,K,I1,J1,K1)    
     DO   IN=1,Nparticles         ! loop over particles 
	  X=XPAR(IN)
	  Y=YPAR(IN)
	  Z=ZPAR(IN)
	  I=INT(X)
	  J=INT(Y)
	  K=INT(Z)

	  D1=X-FLOAT(I)
	  D2=Y-FLOAT(J)
	  D3=Z-FLOAT(K)
	  T1=1.-D1
	  T2=1.-D2
	  T3=1.-D3
	  T2W =T2*WPAR
	  D2W =D2*WPAR
	  I1=I+1
	     IF(I1.GT.NGRID)I1=1
	  J1=J+1
	     IF(J1.GT.NGRID)J1=1
	  K1=K+1
          IF(K1.GT.NGRID)K1=1
!$OMP ATOMIC
          FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
!$OMP ATOMIC
	  FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
!$OMP ATOMIC
	  FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
!$OMP ATOMIC
	  FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
   
!$OMP ATOMIC
	  FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
!$OMP ATOMIC
	  FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
!$OMP ATOMIC
	  FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
!$OMP ATOMIC
	  FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W
	   
      ENDDO 
    Call TimingMain(3,1)
    D1 = 0. ; D2 =0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (M1,M2,M3) REDUCTION(+:D1,D2)
    DO M3=1,NGRID
       DO M2=1,NGRID
          DO M1=1,NGRID
	       D1 = D1 + FI(M1,M2,M3)
	       D2 = D2 + FI(M1,M2,M3)**2
          END DO
       END DO
    END DO
    D1 = D1/(float(NGRID))**3
    D2 = sqrt(D2/(float(NGRID))**3)
    write(*,*) ' Finished density: aver/rms=',D1,D2
    
  END SUBROUTINE DENSIT
!---------------------------------------
!                   generate a vector of gaussian numbers    
  SUBROUTINE getGauss(Gg,jp,kp,N)
!
!---------------------------------------
  use Random
  use LUXURY
  Integer*8 :: kp,jp
  Real*4 :: Gg(N)

    Ns =SeedsPage(jp,kp) 
    lux = 2
    Call rluxgo(lux,Ns,0,0)
    Do i=1,N
      Gg(i) = GAUSS3(gSet,iFlag)
    end Do
  end SUBROUTINE getGauss
!

!-----------------------------------------------------
!                              Density Field    
    SUBROUTINE DENSITrsd(iSwitch,NGRID_old)
!
!           iSwitch = 0 - real space
!                   = 1,2,3 - redshift space 
!----------------------------------------------------
use Random
Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8, parameter :: Nslip= 160481183, Nf=4538127
Integer*8 :: IN, ip,jp,kp,jpp,kpp
Integer*4, save :: Nseed =198239321
Real*4, allocatable :: Gg(:)
Call Timing(2,-1)      ! start reading time

      	XN   =FLOAT(NGRID)+1.-1.E-7
        YN   =FLOAT(NGRID)
        Nseed = 121071 + Nseed
        Xscale = Float(NGRID)/NGRID_old
        write(*,*) ' Inside Densit: Ngrid =',NGRID,iSwitch
        Nthreads = OMP_GET_MAX_THREADS()
!				       Subtract mean density
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3)
      DO M3=1,NGRID
       DO M2=1,NGRID
	      DO M1=1,NGRID
	        FI(M1,M2,M3) = -1.
	      END DO
       END DO
      END DO

      W     = FLOAT(NGRID)**3/float(Nparticles)
      ! factor = sqrt(AEXPN/(Om0+AEXPN**3*Oml0))/100. ! V(km/s)->space
      ! Vscale = (Box/Ngrid)100/Aexpn  ! km/s
      factor = sqrt(AEXPN/(Om+AEXPN**3*OmL))/AEXPN
      sigFactor = sigv/100.*AEXPN*Ngrid/Box
      PARTW = W  
      write(*,'(a,f8.3,a,i5,a,8es12.4)') ' DENSITrsd:Particle weight = ',PARTW, &
           ' Switch=',iSwitch, &
           ' Factors= ',factor,sigFactor
      xmin =1.e+5
      xmax =-1.e+5
      ymin =1.e+5
      ymax =-1.e+5
      zmin =1.e+5
      zmax =-1.e+5
      If(iSwitch/=0)Then
         !write(*,*) ' Allocate Gg: ',Nrow
         Allocate(Gg(Nrow))
         Gg = 0.
      end If
       
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (IN,X,Y,Z,I,J,K,D1,D2,D3,T1,T2,T3,T2W,D2W,I1,J1,K1, &
!$OMP       vv,ip,jpp,kpp,jp,kp,Gg)
      Do kp =1,NROW
         if(mod(kp,250)==0.and.iSwitch==1)write(*,*) ' k=',kp
         kpp = (kp-1)*NROW**2
      Do jp =1,Nrow
         jpp = (jp-1)*Nrow
         If(iSwitch/=0)Then
!$OMP critical            
            Call getGauss(Gg,jp,kp,Nrow)
!$OMP end critical
        end If
      Do ip =1,Nrow
         IN = ip + jpp +kpp
         
        If(iSwitch == 0)Then
            X = (XPAR(IN)-1.)*Xscale+1.
            Y = (YPAR(IN)-1.)*Xscale+1.
            Z = (ZPAR(IN)-1.)*Xscale+1.
         Else
            If(dens(IN)>DensThr)Then
               vv = sigFactor*Gg(ip)*(dens(IN)-DensThr)**0.3333
            Else
               vv  = 0.
            end If
               SELECT CASE (iSwitch)
                CASE (1)
                   Z = (XPAR(IN)-1.)*Xscale+1. + (Vx(IN)*Xscale+vv)*factor
                   Y = (YPAR(IN)-1.)*Xscale+1.
                   X = (ZPAR(IN)-1.)*Xscale+1.
                CASE (2)
                   X = (XPAR(IN)-1.)*Xscale+1.
                   Z = (YPAR(IN)-1.)*Xscale+1. + (Vy(IN)*Xscale+vv)*factor
                   Y = (ZPAR(IN)-1.)*Xscale+1.
                CASE (3)
                   X = (XPAR(IN)-1.)*Xscale+1.
                   Y = (YPAR(IN)-1.)*Xscale+1.
                   Z = (ZPAR(IN)-1.)*Xscale+1. + (Vz(IN)*Xscale+vv)*factor
                End SELECT
      end If

           IF(X.ge.NGRID+1.)X=X-NGRID
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID

	   I=INT(X)
	   J=INT(Y)
	   K=INT(Z)
           If(I.le.0)write (*,*) ' X:',X,Y,Z,' Irow=',IROW,IN
           If(J.le.0)write (*,*) ' Y:',X,Y,Z,' Irow=',IROW,IN
           If(K.le.0)write (*,*) ' Z:',X,Y,Z,' Irow=',IROW,IN
           If(I.gt.NGRID+1)write (*,*) ' X:',X,Y,Z,' Irow=',IROW,IN
           If(J.gt.NGRID+1)write (*,*) ' Y:',X,Y,Z,' Irow=',IROW,IN
           If(K.gt.NGRID+1)write (*,*) ' Z:',X,Y,Z,' Irow=',IROW,IN
                   !---------------------------------------- CIC
	        D1=X-FLOAT(I)
	        D2=Y-FLOAT(J)
	        D3=Z-FLOAT(K)
	        T1=1.-D1
	        T2=1.-D2
	        T3=1.-D3
	        T2W =T2*W
	        D2W =D2*W
	        I1=I+1
	           IF(I1.GT.NGRID)I1=1
	        J1=J+1
	           IF(J1.GT.NGRID)J1=1
	        K1=K+1
         IF(K1.GT.NGRID)K1=1
!$OMP Atomic         
	             FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W 
                !------------------------------------------ NGP
		!      FI(I ,J ,K ) =FI(I ,J ,K ) + W
           ENDDO
           ENDDO
           ENDDO

      If(iSwitch/=0)DEALLOCATE(Gg)

         Call Timing(2,1)      ! start reading time
       END SUBROUTINE DENSITrsd
!--------------------------------------------------                     
!                    : Store and retrieve seeds for parallelization 
!                      of random number generator
!--------------------------------------------------                     
SUBROUTINE SetRandom
   use Tools
   use LUXURY
   use Random
   write(*,*) ' Inside SetRandom'
   write(*,*) NROW,NGRID
   
   ALLOCATE(SeedsPage(NROW,NROW))

      Ns     = Nseed 
   Do j=1,NROW 
      Do i=1,NROW 
         SeedsPage(i,j) = Ns
         dummy = RANDd(Ns) 
      End Do
   End Do
   Nseed = Ns    ! new seed
   !write(*,'(a,6i12)') ' Initialized random seeds: ',(SeedsPage(i,1),i=1,6)
           End SUBROUTINE SetRandom        
!
! 
!---------------------------------------
!                           Density at Particle position    
    SUBROUTINE DENSPART
!
! 
!---------------------------------------
Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8 :: IN
Integer*4 :: Nseed
         Call Timing(3,-1)      ! start time
        write(*,*) ' Inside DensPart: Ngrid =',NGRID

      
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (IN,X,Y,Z,I,J,K,D1,D2,D3,T1,T2,T3,T2W,D2W,I1,J1,K1)
      Do IN =1,Nparticles
            X = XPAR(IN)
            Y = YPAR(IN)
            Z = ZPAR(IN)
           IF(X.ge.NGRID+1.)X=X-NGRID
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID

	   I=INT(X)
	   J=INT(Y)
	   K=INT(Z)
                   !---------------------------------------- CIC
	   D1=X-FLOAT(I)
	   D2=Y-FLOAT(J)
	   D3=Z-FLOAT(K)
	   T1=1.-D1
	   T2=1.-D2
	   T3=1.-D3
	   I1=I+1
	      IF(I1.GT.NGRID)I1=1
	   J1=J+1
	      IF(J1.GT.NGRID)J1=1
	   K1=K+1
              IF(K1.GT.NGRID)K1=1
           dens(IN)=  FI(I ,J ,K )*T3*T1*T2 + &
                      FI(I1,J ,K )*T3*D1*T2 + &
                      FI(I ,J1,K )*T3*T1*D2 + &
                      FI(I1,J1,K )*T3*D1*D2 + &
                      FI(I ,J ,K1)*D3*T1*T2 + &
                      FI(I1,J ,K1)*D3*D1*T2 + &
                      FI(I ,J1,K1)*D3*T1*D2 + &
                      FI(I1,J1,K1)*D3*D1*D2               
        ENDDO 

         Call Timing(3,1)    
       END SUBROUTINE DENSPART

!-----------------------------------------------------
!		                  Compute mean density and rms
      SUBROUTINE DENTES(DELR)
!-----------------------------------------------------
     real*8 :: SUM1,SUM2
     real*8, parameter :: dlog =0.05   
         Call Timing(3,-1)      ! start reading time
      SUM1 = 0.
      SUM2= 0.
      Nn  = 0
      Am  = 0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i) &
!$OMP REDUCTION(+:SUM1,SUM2)
      DO K=1,NGRID
      DO J=1,NGRID
      DO I=1,NGRID
       SUM1 = SUM1 + FI(I,J,K)
       SUM2=  SUM2 + (FI(I,J,K))**2
      ENDDO
      ENDDO
      ENDDO
      Total =(FLOAT(NGRID))**3
      DENM = SUM1/Total
      DELR = DSQRT(SUM2/Total-DENM**2)
      WRITE (*,150)  DELR,DENM
      WRITE (17,150) DELR,DENM
 150     format(20x,'Density is in units of average density', &
                    ' in the Box',/20x,   &
                    ' RMS Delta Rho/Rho   =',G11.4,/20x,  &
                    ' Mean Delta Rho/Rho  =',G11.4)
         Call Timing(3,1)  
     END SUBROUTINE DENTES

!-----------------------------------------------------
!		          Compute statistics of Density
      SUBROUTINE DensDistr
!-----------------------------------------------------
      real*8, parameter :: dlog =0.025   
      real*8 :: SUM1,SUM2
      Integer*8,   Allocatable,DIMENSION(:)   :: nCells
      Real*8,      Allocatable,DIMENSION(:)   :: den
      Integer*8,   Allocatable,DIMENSION(:,:) :: nCellsT
      Real*8,      Allocatable,DIMENSION(:,:) :: denT
      Integer*4  :: OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM

      Call Timing(3,-1)      ! start reading time
         
      DensMax = 1.e5    !-- assumed density maximum
      DensMin = 0.001    !-- assumed density minimum
      iMin = log10(DensMin)/dlog-1
      iMax = log10(DensMax)/dlog+1
      iThreads = OMP_GET_MAX_THREADS() 
      allocate(nCells(iMin:iMax),den(iMin:iMax))
      allocate(nCellsT(iMin:iMax,iThreads),denT(iMin:iMax,iThreads))

      nCells(:)    =0
      den(:)       =0.
      nCellsT(:,:) =0
      denT(:,:)    =0.
      
      SUM1 = 0.
      SUM2= 0.
      Nn  = 0
      Am  = 0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i) &
!$OMP REDUCTION(+:SUM1,SUM2)
      DO K=1,NGRID
      DO J=1,NGRID
      DO I=1,NGRID
       SUM1 = SUM1 + FI(I,J,K)+1.
       SUM2=  SUM2 + (FI(I,J,K)+1.)**2
      ENDDO
      ENDDO
      ENDDO
      Total =(FLOAT(NGRID))**3
      DENM = SUM1/Total
      DELR = DSQRT(SUM2/Total-DENM**2)
      WRITE (*,150)  DELR,DENM,NGRID
      WRITE (18,150) DELR,DENM,NGRID
 150     format(20x,'Density is in units of average density', &
                    ' in the Box',/20x,   &
                    ' RMS Delta Rho/Rho   =',G11.4,/20x,  &
                    ' Mean Delta Rho/Rho  =',G11.4,/20x,  &
                    ' Number grid points  =',i5)

!$OMP PARALLEL DO DEFAULT(SHARED)  &
!$OMP PRIVATE ( ind,k,j,i,iOMP)
      Do K=1,NGRID
         iOMP = OMP_GET_THREAD_NUM()+1
      DO J=1,NGRID
         DO I=1,NGRID
           ind = Floor(log10(FI(I,J,K)+1.)/dlog)
           ind = MIN(MAX(ind,iMin),iMax)
           nCellsT(ind,iOMP) = nCellsT(ind,iOMP) +1
           denT(ind,iOMP)    = denT(ind,iOMP) +FI(I,J,K)+1.
        END DO
     END DO
  END Do
  DO i=iMin,iMax            ! sum results of all threads
     DO iP=1,iThreads
        nCells(i) = nCells(i) + nCellsT(i,iP)
        den(i)    = den(i)    + denT(i,iP)
     end DO
  END DO
  dNorm = float(Ngrid)**3
  write(18,*)' dens_left   dens_right  density   dens*dN/d(dens)/Ncells  cells'
  DO i=iMin,iMax           ! print results
     d1 = 10.**(i*dlog)
     if(i==iMin)d1 =0.
    d2 = 10.**((i+1)*dlog) 
     den(i) = den(i)/max(nCells(i),1)
     write(18,'(3es12.4,3x,es13.5,i12)') d1,d2,den(i),   &
          nCells(i)/(d2-d1)/dnorm*den(i),nCells(i)     
  END DO
  
      deallocate(nCells,den)
      deallocate(nCellsT,denT)
      
         Call Timing(3,1)  
END SUBROUTINE DensDistr

!--------------------------------------------
!
!          Make density field of galaxies    
!
subroutine DensGal(iSwitch)
!
!--------------------------------------------
use Tools
character*120 :: Name

Real*8 :: ss,X,Y,Z,D1,D2,D3,T1,T2,T3
Integer*8 ::  ip
Call Timing(2,-1)      ! start reading time

      	XN   =FLOAT(NGRID)+1.-1.E-7
        YN   =FLOAT(NGRID)
!				       Subtract mean density
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (M1,M2,M3)
      DO M3=1,NGRID
       DO M2=1,NGRID
	      DO M1=1,NGRID
	        FI(M1,M2,M3) = -1.
	      END DO
       END DO
      END DO

      W     = FLOAT(NGRID)**3/float(Ngalaxies)

      Xscale = NGRID/Box
      Vscale = Xscale/100.*sqrt(AEXPN/(Om+AEXPN**3*OmL))
      PARTW = W  
      write(*,'(/a,es12.4,a,i5,a,8es12.4)') ' DensGal: Particle weight = ',PARTW, &
           ' Switch=',iSwitch, &
           ' Factors= ',Xscale,Vscale
      xmin =1.e+5
      xmax =-1.e+5
      ymin =1.e+5
      ymax =-1.e+5
      zmin =1.e+5
      zmax =-1.e+5

       
!$OMP PARALLEL DO DEFAULT(SHARED) & 
!$OMP PRIVATE (ip,X,Y,Z,I,J,K,D1,D2,D3,T1,T2,T3,T2W,D2W,I1,J1,K1)
      Do ip = 1,Ngalaxies
         !if(mod(ip,50000)==0)write(*,*) ' galaxy=',ip
         If(iSwitch == 0)Then
            X = Xb(ip)*Xscale+1.
            Y = Yb(ip)*Xscale+1.
            Z = Zb(ip)*Xscale+1.
         Else
               SELECT CASE (iSwitch)
               CASE (1)
                    Z = Xb(ip)*Xscale+1. +VXb(ip)*Vscale
                    Y = Yb(ip)*Xscale+1.
                    X = Zb(ip)*Xscale+1.
                 CASE (2)
                    X = Xb(ip)*Xscale+1. 
                    Z = Yb(ip)*Xscale+1. +VYb(ip)*Vscale
                    Y = Zb(ip)*Xscale+1.
                CASE (3)
                    X = Xb(ip)*Xscale+1. 
                    Y = Yb(ip)*Xscale+1.
                    Z = Zb(ip)*Xscale+1. +VZb(ip)*Vscale
                End SELECT
      end If

           IF(X.ge.NGRID+1.)X=X-NGRID
           IF(Y.ge.NGRID+1.)Y=Y-NGRID
           IF(Z.ge.NGRID+1.)Z=Z-NGRID

           IF(X.lt.1.)X=X+NGRID
           IF(Y.lt.1.)Y=Y+NGRID
           IF(Z.lt.1.)Z=Z+NGRID

	   I=INT(X)
	   J=INT(Y)
	   K=INT(Z)
           If(I.le.0)write (*,'(a,3es15.5,a,2i12)') ' X:',X,Y,Z,' Irow=',ip
           If(J.le.0)write (*,'(a,3es15.5,a,2i12)') ' Y:',X,Y,Z,' Irow=',ip
           If(K.le.0)write (*,'(a,5es15.5,a,2i12)') ' Z:',Z,Xb(ip),VXb(ip),Xb(ip)*Xscale+1.,VXb(ip)*Vscale,' Irow=',ip
           If(I.gt.NGRID+1)write (*,'(a,3es15.5,a,2i12)') ' X:',X,Y,Z,' Irow=',ip
           If(J.gt.NGRID+1)write (*,'(a,3es15.5,a,2i12)') ' Y:',X,Y,Z,' Irow=',ip
           If(K.gt.NGRID+1)write (*,'(a,5es15.5,a,2i12)') ' Z:',Z,Xb(ip),VXb(ip),Xb(ip)*Xscale+1.,VXb(ip)*Vscale,' Irow=',ip
                   !---------------------------------------- CIC
	        D1=X-FLOAT(I)
	        D2=Y-FLOAT(J)
	        D3=Z-FLOAT(K)
	        T1=1.-D1
	        T2=1.-D2
	        T3=1.-D3
	        T2W =T2*W
	        D2W =D2*W
	        I1=I+1
	           IF(I1.GT.NGRID)I1=1
	        J1=J+1
	           IF(J1.GT.NGRID)J1=1
	        K1=K+1
         IF(K1.GT.NGRID)K1=1
!$OMP Atomic         
	             FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
!$OMP Atomic         
                     FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
!$OMP Atomic         
	             FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
!$OMP Atomic         
	             FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
!$OMP Atomic         
	             FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
!$OMP Atomic         
	             FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W 
                !------------------------------------------ NGP
		!      FI(I ,J ,K ) =FI(I ,J ,K ) + W
           ENDDO

         Call Timing(2,1)      ! start reading time

end subroutine DensGal



!------------------------------------------------------------
   end Module Density
   

!
! _______________________ START 3-D SIMULATIONS                         
!                                                                       
!                         Klypin, August 2015                         
!
!   Uses FFT5 :  Nrow = 4*3**n*5**m
!

!-------------------------------------------------------------
Module setInitialConditions
  use Tools
  
    Integer*4, parameter  :: NtabM = 100000
    REAL*8,    PARAMETER  :: PI=3.1415926535
      Real*4   :: xkt(0:NtabM),Pkt(0:NtabM)  ! Power spectrum
      Real*4   :: StepK,alog0
      Integer*4:: Ntab
      Character (LEN=12), DIMENSION(4) :: NamesDir=(/ &
                                         'PMstrtXX.DAT','PMstrtVX.DAT',&
                                         'PMstrtYY.DAT','PMstrtVY.DAT'/)

      INTEGER*4,parameter ::   NPAGE = 1024**2  ! # particles in a record
      INTEGER*4 ::   NMAX,  &  !  = NGRID/2,     & 
                     NRECL, &  !  = NPAGE*6,   & ! # particles in a record
                     NSPEC     ! = NROW/2             ! No waves shorter than Ny
      REAL*4                ::   alpha  =0.,   &
                                 Qscale
      Integer*4, PARAMETER :: nbyteword = 4      ! defines length of direct-access:1or4
      Integer*4, PARAMETER :: LevelParticles = 0 


      REAL,   ALLOCATABLE, DIMENSION(:,:,:) ::  GR
      Real*4, ALLOCATABLE, DIMENSION(:)     ::  XPt,VXt
      Real*4,  DIMENSION(NPAGE)              ::  xbuf,vxbuf
      Real*4, dimension(NPAGE) ::  XPp,YPp,ZPp, &
     		                   VPx,VPy,VPz

!$OMP THREADPRIVATE(xbuf,vxbuf)
!$OMP THREADPRIVATE(XPp,YPp,ZPp,VPx,VPy,VPz)

    Contains
!--------------------------------------------------                     
!                              sqrt(Power spectrum)
!                                       k = (2pi/L) 
      FUNCTION TRUNF (WK) 
!-------------------------------------------------
        real*4 :: k
      IF (WK.GE.FLOAT (NSPEC) ) THEN 
         TRUNF = 0. 
         RETURN 
      ENDIF 
          k = QSCALE * wk 
          TRUNF = sqrt (Ppk (k) )
        END FUNCTION TRUNF

!---------------------------------------                                
!            interpolate table with p(k)                       
!                                             
        FUNCTION Ppk (x) 
! 
!---------------------------------------                                                        
                              ! slope is ns =-3                         
      If (x.ge.xkt (Ntab) ) Then 
         Ppk = Pkt (Ntab) / (x / xkt (Ntab) ) **3 
         Return 
      EndIf 
                                         
      If (x.lt.xkt (1) ) Then    ! slope is ns=1                
         Ppk = Pkt (1) * (x / xkt (1) ) 
         Return 
      EndIf 
      ind = INT ( (log10 (x) - alog0) / StepK) + 1 
      dk = xkt (ind+1) - xkt (ind) 
      Ppk = (Pkt (ind) * (xkt (ind+1) - x) + Pkt (ind+1) * (x - xkt (   &
                   ind) ) ) / dk                                                     

      END FUNCTION Ppk    

!--------------------------------------------------                     
!               Store  seeds for parallelization 
!                      of random number generator
!
      SUBROUTINE SetRandom
!                   
!-------------------------------------------------                                                     
Use Random
INCLUDE 'luxuryp.h' 

      gSet  = 0.
      Ns    = Nseed   
      lux   = 2
      Do k=1,NROW                               ! luxury
          CALL rluxgo (lux, Ns, 0, 0)      ! initialize luxury 
            !if(k/10*10==k)write(*,*) ' page number=',k
          i24A(k)    = i24
          j24A(k)    = j24
          in24A(k)  = in24
          kountA(k) = kount
          carryA(k) = carry
          gsetA(k)  = gSet
          iFlagA(k) = iFlag
             Do m=1,24
                 SeedsPage(m,k) = seeds(m)
             EndDo    
              dummy = RANDd(Ns) 
        end Do
End SUBROUTINE SetRandom
                                                   
!-------------------------------------------------
!        
SUBROUTINE ReadPkTable
!                   
!-------------------------------------------------                                                     
            !-- Read PkTable. Assign Omegas
      OmbT = ParseLine(10)
      OmcT = ParseLine(10)
      OmLT = ParseLine(10)
      OmT  = ParseLine(10)
      sigma8T = ParseLine(10)
           !-- check parameters
      If(abs(OmbT-Omb)>1.e-4*OmbT)Stop 'Error: OmegaBaryons is different from PkTable value'
      If(abs(OmT-Om)>1.e-4*OmT)Stop 'Error: OmegaMatter is different from PkTable value'
      write(*,'(3(a,ES12.3))')'  Ombar =', &
                       OmbT,' Om_matter =',OmT
      Ntab = 0
12   READ (10, *, end = 32, err = 32) xx, pp    !---- read P(k)
        Ntab = Ntab + 1 
        xkt (Ntab) = xx  !*hubble 
        Pkt (Ntab) = pp                !  Pk 
      GOTO 12 
32    close(10)
      !Write (*,*) ' Read ', Ntab, ' lines from P(k) table '
      
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
   END SUBROUTINE ReadPkTable
  
!-------------------------------------------------------------
!                                                                       
SUBROUTINE Initialize
!     
!-------------------------------------------------------------
                                                
   write(*,'(a,$)') ' Enter realization number for this run = '
   Read(*,*) Nrealization
   write(*,*)
   open(1,file='../TableSeeds.dat')
   read(1,*)     ! skip the first line == header
   Do ij=1,Nrealization
      read(1,*) Nseed0,Ncount
   endDo
   Nseed = Nseed0
   close(1)

   extras (:)  = 0.
   extras(100) = Box
   NMAX        = NGRID/2
   NRECL       = NPAGE*6
   NSPEC       = NROW/2
   NBYTE = NPAGE * 6 * 4 
   Nparticles  = INT(NROW,8)**3 
   ISTEP = 0                                                        
   TINTG = 0.                                                       
   AU0   = 0.                                                       
   AEU0  = 0.                                                       
   EKIN  = 0.                                                       
   EKIN1 = 0.                                                       
   EKIN2 = 0.                                                       
   QSCALE = 2.*PI/Box
   myMemory= Memory(2*Nparticles)
   Allocate(XPt(Nparticles),VXt(Nparticles))
   write(*,*) ' NROW   =',NROW
   write(*,*) ' NGRID  =',NGRID
   write(*,*) ' Npart  =',Nparticles
   ASTEP = ASTEP0
   AEXPN = AEXPN0
END SUBROUTINE Initialize

!---------------------------------------------------
!
!             check if all init files are present
!
   Subroutine CheckInit
!---------------------------------------------------
     logical :: exst
      Inquire(file='../PkTable.dat',exist = exst)
      if(.not.exst)Stop ' File PkTable with the power spectrum not found'
      open(10,file='../PkTable.dat')
      
      Inquire(file='../Setup.dat',exist = exst)
      if(.not.exst)Then
         write(*,*)' Error: File ../Setup.dat not found. Run PMP2init.exe'
         stop
      end if
         open(11,file='../Setup.dat')

      Inquire(file='../TableSeeds.dat',exist = exst)
      if(.not.exst)Call SetSeeds
 
end Subroutine CheckInit
    
!------------------------------
!
!             Generate a table with random seeds
!
!------------------------------
Subroutine SetSeeds
  use Random
  integer*8 :: is,ij,iM
  integer*4 :: Nseed0
  is = 1232_8**3
  
    Nseed0  = 1298302
    nslip = 137
    noff  = 2357
    Ntable =5000
    NCount =0
    open(1,file='TableSeeds.dat')
    write(1,*) 'Seeds:',Nseed0,Nslip
   Do ij=1,is
      x  =RANDd(Nseed0)
      Nn =INT(x*nslip)+1
      Do jj =1,noff+Nn
         x  =RANDd(Nseed0)
      End Do
      Ncount =Ncount +1
      write(1,*) Nseed0,Ncount
      If(Ncount>Ntable)exit
   endDo
   close(1)
 end Subroutine SetSeeds
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
!------------------------------------------------                       
!                             Make a realization of spectrum of perturbations
!                             ALPHA  = normalization factor for displacement
!                             NRAND = seed for random numbers           
SUBROUTINE SPECTR (iDirection) 
!------------------------------------------------                       
use Random
INCLUDE 'luxuryp.h' 
      REAL(8)   :: SUMM=0., Wi3, Wj3, Wk3, WD, TS, TRX,ss 
      REAL      ::    t0=0,t1=0.        ! counters for timing
      INTEGER*4, allocatable, DIMENSION(:) :: mapz, map3

      CALL Timing(3,-1)                  ! initialize time counter
      allocate(mapz(NROW), map3(NROW))

      SUMM = 0.
      map3(1) = 0
      mapz(1) = 1
      map3(NROW) = 0
      mapz(NROW) = NROW
      DO j=2,NROW-1,2
          j3 = j/2
          map3(j)   = -j3   ! -k   for cos
          map3(j+1) =  j3   !  k       sin
          mapz(j)   =  j+1  ! flip cos <--> sin
          mapz(j+1) =  j
      end DO
                                                           
!$OMP PARALLEL DO DEFAULT(SHARED)  &                                     
!$OMP PRIVATE(Mk3,Mk2,Mk1)                                            
      DO Mk3 = 1, NROW 
      DO Mk2 = 1, NROW 
      DO Mk1 = 1, NROW 
          FI (Mk1, Mk2, Mk3) = 0. 
      ENDDO 
      ENDDO 
      ENDDO 
!                                      Set Spectrum                     
!$OMP PARALLEL DO DEFAULT(SHARED)  &                                     
!$OMP PRIVATE(k,ksign,kz,k3,Wk3,j,jsign,jz,j3,Wj3,i,isign,iz,i3,Wi3)  &
!$OMP PRIVATE(WD,Wk,TS,TRX,m,NRAND,gSet,iFlag)           &
!$OMP REDUCTION(+:SUMM)
      DO k = 1, NROW 
           i24      = i24A(k) 
           j24      = j24A(k) 
           in24    = in24A(k) 
           kount   = kountA(k)
           mkount  = 0 
           carry  = carryA(k)
           iflag  = iFlagA(k)
           gset   = gSetA(k)
              Do m=1,24
                  seeds(m) =SeedsPage(m,k) 
              EndDo    
         Wk3 = map3(k)**2 

         DO j = 1, NROW 
            Wj3 = map3(j)**2 
            DO i = 1, NROW 
               Wi3 = map3(i)**2 
                TS = GAUSS3 (gSet,iFlag) 
               WD  = Wi3 + Wj3 + Wk3 

               IF (WD< 1.e-3) THEN 
                  FI (1, 1, 1) = 0. 
               ELSE 
                  Wk = SQRT (WD) 
                  TS = TRUNF (Wk) * TS
                  TRX = TS / WD 
                                                           ! start switch
                  If (iDirection.eq.1) Then 
                     FI (mapz(i), j, k) = TRX * map3(i)
                  ElseIF (iDirection.eq.2) Then 
                     FI (i, mapz(j), k) = TRX * map3(j) 
                  Else 
                     FI (i, j, mapz(k)) = TRX * map3(k)
                  EndIf                                  ! end iDirection
                  SUMM = SUMM + TS**2 
               ENDIF        ! end kx=ky=kz=0 switch
            ENDDO         ! i
         ENDDO            ! j
      ENDDO               ! k
      
      IF (SUMM.LE.0.) Write (*,  * ) ' Error!!! Summ over spectrum = 0'
      ALPHA = AMPLT / SQRT (SUMM) * sqrt (8.) 
      Write (*,'(10x,a20,3g13.6)') 'SPECTR ==>  SUMM=', SUMM, ALPHA 
      
      Do i =1,Nrow             ! adjust zero-k values of FI
         Do j=1,NROW
            FI(i,j,1) = FI(i,j,1)/sqrt(2.)
            FI(i,1,j) = FI(i,1,j)/sqrt(2.)
            FI(1,i,j) = FI(1,i,j)/sqrt(2.)
         End Do
      End Do

      deallocate(mapz, map3)
      CALL Timing(3,1)                  

      END SUBROUTINE SPECTR                         
!------------------------------------------------                       
!		                        Define coordinates and velocities for        
!                                  particles with resolution = Indx     
!                                  Indx =1 = high resolution            
!                                  Indx =2 = medium resolution          
!                                  Icurrent = number of particles       
      SUBROUTINE BLOCKS (XCONS, VCONS, Indx, Icurrent, iDirection) 
!------------------------------------------------                       
      Integer*8        nLevel (10) 
      Real(8)       :: sDispl =0. 
      REAL          :: t0=0,t1=0.        ! counters for timing
      Integer*8     :: Icurrent,JROW
                                 
      CALL Timing(4,-1)                    ! initialize time counter   

      sDispl = 0. 
      xShift =  0.5 !/2**LevelParticles
      JROW   = NROW  ! make it long integer
      write(*,*) ' XCONS, VCONS =',XCONS, VCONS
      
      QFACT = FLOAT (NGRID) / FLOAT (NROW)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ip,jp,kp,DX,Icurrent,Q) &
!$OMP REDUCTION(+:sDispl)      
      Do kp = 1, NROW
      Do jp = 1, NROW
      Do ip = 1, NROW
         DX = FI (ip, jp, kp)              !  find displacements 
         Icurrent = ip+ (jp-1_8)*JROW +(kp-1_8)*JROW*JROW 
         If (iDirection.eq.1) Then 
            Q = QFACT * (ip - 1.) + 1. 
         ElseIF (iDirection.eq.2) Then 
            Q = QFACT * (jp - 1.) + 1. 
         Else 
            Q = QFACT * (kp - 1.) + 1. 
         EndIf 
         XPt (Icurrent) = Q - XCONS * DX + xShift
         VXt (Icurrent) = VCONS * DX
         sDispl = sDispl + DX**2
      EndDo        ! end ip
      EndDo        ! end jp
      EndDo        ! end kp
                                                                                           
      Write (16, '(10x,a20,g12.4,a,i3)') ' RMS 3d diplacement=',  &
           XCONS*sqrt (sDispl/max(iCurrent,1)),                   &
           ' In direction = ',iDirection   

      CALL Timing(4,1)                                                                           
      END  SUBROUTINE BLOCKS                                         
!------------------------------------------------                       
!                                              				   FFT of the spectru
      SUBROUTINE VECTOR  
!------------------------------------------------                       
use fft5
      integer*4, parameter :: Nlensav = 8192
      integer*4, parameter :: Nlenwrk = 8192
      real*8,    parameter :: P16 = 6.28318530718
      real*4,    parameter :: sq2 = 1. !1.41421356237
      real*8,  save        :: wsave(1:Nlensav)
      real*8,  save        :: work(1:Nlenwrk)
      REAL*8               :: XX,D1,D2,A1,A2,A3,wi,wj,wk
      Integer*4            :: OMP_GET_MAX_THREADS,OMP_GET_THREAD_NUM
      Integer*4            :: Ng,ier,lensav,lenwrk,lenr,inc
      real*8               :: r(Nlenwrk)
      REAL                 :: t0=0,t1=0.        ! counters for timing
      Real*8               :: ss
!$OMP THREADPRIVATE(work,wsave)

      CALL Timing(5,-1)

      If(NROW>Nlenwrk)Stop ' Incresase Nlenwrk in POTENTfft5'
      Ng = NROW
      lensav = NROW+int(log(real(NROW,kind = 4))/log(2.0E+00))+4
      lenwrk = NROW

      call rfft1i ( Ng, wsave, lensav, ier ) !   Initialize FFT
      inc  = 1
      lenr = Ngrid

             write(*,*)  ' XY fft'
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k=1,NROW             ! fft for xy planes in x-dir
       Do j=1,NROW     
          Do i=1,NROW
             r(i) = FI(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do i=1,NROW
            FI(i,j,k) = r(i)
         EndDo
       EndDo
       
       Do i=1,NROW        ! fft xy planes in y-dir
          Do j=1,NROW
             r(j) = FI(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do j=1,NROW
            FI(i,j,k) = r(j)
          EndDo
       EndDo
    EndDo

      write (*,'(10x,a)') 'Swap k<->i'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = FI(I,J,K)
               FI(I,J,K) =FI(K,J,I)
               FI(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO

!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r, ier)
     Do j=1,NROW     ! ------ z-direction
       Do i=1,NROW     
           Do k=1,NROW
             r(k) = FI(k,j,i)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
           Do k=1,NROW
             FI(k,j,i) = r(k) 
          EndDo
       EndDo   
    end Do
      
      write (*,'(10x,a)') 'Swap i<->k'
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE ( k,j,i,aa)
      DO J=1,NROW
      DO K=1,NROW-1
            DO I=K+1,NROW
               aa = FI(I,J,K)
               FI(I,J,K) =FI(K,J,I)
               FI(K,J,I) =aa
            ENDDO
         ENDDO
      ENDDO

      Write (*,'(10x,a)') ' FFT is done' 
      CALL Timing(5,1)

      RETURN 
      END SUBROUTINE VECTOR                         
!---------------------------------------------                          
!               1)    Update SKINE and Wtotal
!               2)    Write current data to  disk/tape                
!                   do not write files for 3rd direction
      SUBROUTINE WriteData (Vscale, AexpV, Wtotal,SKINE,iDirection) 
!----------------------------------------------                         
Real*8             ::    Wtotal,SKINE,vvx,vx2
REAL               ::    t0=0,t1=0.        ! counters for timing
Integer*8          ::  i,jfirst,jlast,npp,JPAGE,j0

      JPAGE = NPAGE
      CALL Timing(6,-1)                    ! initialize time counter
      jfirst = 1
      vvx  = 0. 
      vx2  = 0. 
      W    =  (float(NROW)/NGRID)**3
      Wtotal = 0.
      vvx    = 0.
      vx2    = 0.
      SKINE  = 0.
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP REDUCTION(+:vvx,vx2,SKINE,Wtotal)
            Do i = 1,Nparticles 
               vvx = vvx + VXt (i) 
               vx2 = vx2 + VXt (i) **2 
               SKINE = SKINE+W * VXt (i) **2
               Wtotal = Wtotal + w 
            EndDo                             ! end current species
         v2    = sqrt ( vx2 / Nparticles) / AexpV * Vscale 
         vvx  = vvx / Nparticles / AexpV * Vscale 
         write(*,*) ' WriteData: ',Nparticles,AexpV,Vscale
         Write ( *,'(22x,"V(km/s)=",2g11.3," Npart=",i12)')   & 
                                        vvx, v2, Nparticles
If(iDirection == 3)Then
      CALL Timing(6,1)
      Return                   ! do not write the last page
EndIf

                                 ! write data page by page    
   Npages  = (Nparticles-1)/NPAGE+1
   write(*,*) ' Npages=',Npages
   If(Npages.eq.0)stop ' Zero number of pages requested. Stop'


Do i=1,Npages
   If(i==Npages)Then
      NinPage = Nparticles -(i-1)*JPAGE  ! # particles in the current page
   Else
      NinPage = JPAGE
   EndIf
   jfirst  = (i-1)*JPAGE +1
   jlast   = jfirst + NinPage-1
   If(mod(i,100)==0.or.i==Npages)write(*,'(10x,a,i5,a,4i11)')    &
        'Write page=',i,' Particles=',NinPage,jfirst,jlast
   Do j0 = jfirst,jlast
      xbuf(j0-jfirst+1)  = XPt(j0)
      vxbuf(j0-jfirst+1) = Vxt(j0)
   EndDo
   write(1,rec=i)xbuf
   write(2,rec=i)vxbuf
EndDo
   
      CALL Timing(6,1)

      END  SUBROUTINE WriteData                                         
                                 
                                                                        
!---------------------------------------------------                    
!                                  Read  current data from disk/tape,   
!                                  Open files                           
!                                  Nrecl is the number of values in a re
!                                  Npage is the number of particles in a
      SUBROUTINE FilesOpen 

Character (LEN=50) :: FileName
Logical            :: FileExists 

      !write(*,'(10x,a/,2(30x,a/))')'Create temporary Files====>',NamesDir

!                                     Open control file               
      OPEN(9, FILE = 'PMcrd.DAT', form = 'unformatted') 
                            ! this clears old header file               
      Write(9) 
      CLOSE (9) 
      OPEN(9, FILE = 'PMcrd.DAT', FORM = 'UNFORMATTED', STATUS =       &
         'UNKNOWN')                                                        
      OPEN(16, FILE = 'Results.log') 
      OPEN(25, file = 'pt.dat', form = 'unformatted') 
                              

      RETURN 
      END SUBROUTINE FilesOpen                      
!---------------------------------------------                          
!                       Write current data to  disk/tape                
!                                                                       
      SUBROUTINE FilesWrite 
!----------------------------------------------                         
!                                       write header and control data   
REAL               ::    t0=0,t1=0.        ! counters for timing
Integer*8          ::    j0,jfirst,jlast,NinPage,JPAGE,Npages
Character*80       ::    fname
      CALL Timing(7,-1)                   ! initialize time counter

      Write (9) HEADER, AEXPN, AEXP0, AMPLT, ASTEP, ISTEP, 0., TINTG,&
      EKIN, EKIN1, EKIN2, AU0, AEU0, NROW, NGRID, Nrealization, Nseed,    &
      Om, OmL, hubble, Nparticles, extras                             
      REWIND 9 

      NBYTE = NPAGE*4
      NACCES= NBYTE / nbyteword
      JPAGE = NPAGE              ! int8
      Nrecpage = 256            ! max number of records per file
      
      Npages  = (Nparticles-1)/JPAGE+1
      Nfiles  = (Npages-1)/Nrecpage +1    ! number of files
      write(*,*) ' Npages=',Npages
      write(*,*) ' Nparts=',Nparticles
      If(Npages.eq.0)stop ' Zero number of pages requested. Stop'      
      
       OPEN(21,FILE=NamesDir(1),ACCESS='DIRECT',form='unformatted', &
                         STATUS='UNKNOWN',RECL=NACCES)  ! x coord and veloc
       OPEN(22,FILE=NamesDir(2),ACCESS='DIRECT',form='unformatted', &
                         STATUS='UNKNOWN',RECL=NACCES)  ! y
       OPEN(23,FILE=NamesDir(3),ACCESS='DIRECT',form='unformatted', &
                         STATUS='UNKNOWN',RECL=NACCES)  ! x coord and veloc
       OPEN(24,FILE=NamesDir(4),ACCESS='DIRECT',form='unformatted', &
                         STATUS='UNKNOWN',RECL=NACCES)  ! y
      XMAX = FLOAT (NGRID) + 1. 
      XSHF = FLOAT (NGRID) 
      NBYTE = NRECL*4
      NACCES= NBYTE / nbyteword

      Nfiles = (Npages-1)/Nrecpage +1

      Do i=1,Nfiles                       !-- open files
         ifile = (i-1)
         If(ifile<10)Then
            write(Fname,'(a,i1.1,a)')'PMcrs',ifile,'.DAT'
         Else
            write(Fname,'(a,i2.2,a)')'PMcrs',ifile,'.DAT'
         EndIf
         Open(30+ifile,file=TRIM(Fname),ACCESS='DIRECT', &
             FORM='unformatted',STATUS='UNKNOWN',RECL=NACCES)
         write(*,'(2i7,2a,3x,i9)') i,ifile,' Open file = ',TRIM(Fname)
      End Do
         
 Do i=1,Npages         !-------- read buffers and dump into final files
   If(i==Npages)Then
      NinPage = Nparticles -(i-1)*JPAGE  ! # particles in the current page
   Else
      NinPage = JPAGE
   EndIf
   
   jfirst  = (i-1)*JPAGE +1
   jlast   = jfirst + NinPage-1
   If(mod(i,100)==0.or.i==Npages)    &
        write(*,'(10x,a,i5,a,4i11)')'Write page=',i,' Particles=',NinPage,jfirst,jlast
   read(21,rec=i)xbuf
   read(22,rec=i)vxbuf
   XPp = xbuf
   VPx   = vxbuf
   read(23,rec=i)xbuf
   read(24,rec=i)vxbuf
   YPp = xbuf
   VPy   = vxbuf
   Do j0 = jfirst,jlast
      ZPp(j0-jfirst+1) = XPt(j0)
      VPz (j0-jfirst+1)  = Vxt(j0)
   EndDo

   Do j=1,NinPage     !-----                 Impose  Periodical boundary conditions 
        IF (XPp(j) .GT.XMAX) XPp(j) = XPp(j) - XSHF 
        IF (XPp(j) .LE.1.)   XPp(j) = XPp(j) + XSHF 
        If (XPp (j) .lt.1..or.XPp(j) .GT.XMAX)   &
             write(*,'(a,2i10,3es13.5)') ' Boundary Error X!!', i,j,XPp(j),YPp(j),ZPp(j)
        If (XPp (j) .ge.XMAX)Then
         !write (*,*)' Fixing Boundary:', j, XPp(j)  
            XPp (j) = XPp (j) - 1.e-4 
         ENDIF
         
        IF (YPp(j) .GT.XMAX) YPp(j) = YPp(j) - XSHF 
        IF (YPp(j) .LE.1.)   YPp(j) = YPp(j) + XSHF 
        If (YPp(j) .lt.1..or.YPp(j) .GT.XMAX)   &
             write(*,'(a,2i10,3es13.5)') ' Boundary Error Y!!', i,j,XPp(j),YPp(j),ZPp(j)
        If (YPp (j) .ge.XMAX)Then
            !write (*,*)' Fixing Boundary:', j, YPp(j)  
            YPp (j) = YPp (j) - 1.e-4 
         ENDIF
         
        IF (ZPp(j) .GT.XMAX) ZPp(j) = ZPp(j) - XSHF 
        IF (ZPp(j) .LE.1.)   ZPp(j) = ZPp(j) + XSHF 
        If (ZPp(j) .lt.1..or.ZPp(j) .GT.XMAX)   &
             write(*,'(a,2i10,3es13.5)') ' Boundary Error Z!!', i,j,XPp(j),YPp(j),ZPp(j)

        If (ZPp (j) .ge.XMAX)Then
            !write (*,*)' Fixing Boundary:', j, ZPp(j)  
            ZPp (j) = ZPp (j) - 1.e-4 
        ENDIF
   EndDo

     WRITE ((i-1)/Nrecpage+30,REC=mod(i-1,Nrecpage)+1) XPp,YPp,ZPp,VPx,VPy,VPz

 EndDo                            ! end Npages                             
      CALL Timing(7,1)
      close(21)
      close(22)
      close(23)
      close(24)
      Call SYSTEM('rm -f PMstrtXX.DAT')
      Call SYSTEM('rm -f PMstrtVX.DAT')
      Call SYSTEM('rm -f PMstrtYY.DAT')
      Call SYSTEM('rm -f PMstrtVY.DAT')
      
      END SUBROUTINE FilesWrite                          
    end MODULE setInitialConditions
                                                                          
!-------------------------------------------------------------------------
PROGRAM  PMstartMp
  use      setInitialConditions
  use      Random
  use      FFT5
   Real*8             ::    Wtotal=0., SKINE=0.
   Integer*8          :: Icurrent
   CALL Timing(0,-1)
   
   CALL Timing(1,-1)
   Call CheckInit     ! test setup, open input files
   Call ReadSetup
   Call Initialize
   Call ReadPkTable
   Call SetRandom
   CALL Timing(1,1)

      write(*,*) ' ---- Nrow= ',Nrow
      Wtotal = (float(NROW)**3*4+2.*Nparticles*4. &
                + 31.*NROW*4                    &
                + 6. *NPAGE*4)/1024.**3
      write (*,'(//10x,a,g13.4,a/)') ' Memory required=',Wtotal,'Gb'
      CALL FilesOpen                   ! this opens files on disk
         AEXP0 = AEXPN 
         AEXPV = AEXPN - ASTEP / 2./2**LevelParticles 
         Fact   = sqrt (Om + OmL * AEXPV**3) 
         QFACT  = FLOAT (NGRID) / FLOAT (NROW) 
         Vscale = Box * 100. / NGRID 

         SKINE   = 0.                              ! counter for  kinetic energy
         Wtotal  = 0.                             ! counter for  total weight
         myMemory = Memory(Nparticles)
      ALLOCATE(FI(NROW,NROW,NROW))     ! allocate memory for spectrum

      Do iDirection = 1, 3                  ! ----------------- loop over x,y,z directions 
          Write(*,'(//3x,a12,i3,15("-"))' )'- Direction=',iDirection
          If(iDirection < 3) Then         
             NBYTE = NPAGE*4
             NACCES= NBYTE / nbyteword
                 OPEN(1,FILE=NamesDir(2*iDirection-1),ACCESS='DIRECT',FORM='UNFORMATTED',        &
                        STATUS ='UNKNOWN',RECL=NACCES)                                  
                 OPEN(2,FILE=NamesDir(2*iDirection),ACCESS='DIRECT',FORM='UNFORMATTED',        &
                        STATUS ='UNKNOWN',RECL=NACCES)                                  
              EndIf
              write(*,*) ' Go to Spectrum '
          CALL SPECTR (iDirection) 
			       !   get the displacement vector by FFT 
          VCONS = - ALPHA/(2.*PI/NGRID)*(AEXPV/AEXP0)*SQRT(AEXPV)*Fact
          XCONS = 	ALPHA / (2. * PI / NGRID) * (AEXPN / AEXP0) 
            Write(*,'(3x,a12,g12.4,a,g12.4)' ) 'Scaling:(x)=', XCONS, ' (v)=', VCONS 

          CALL VECTOR
             Write (*,'(3x,a,ES12.4)' ) 'SPECTR done: Alpha=', ALPHA 
          Icurrent = 0                            ! current particle   
          CALL BLOCKS (XCONS, VCONS, kSpec, Icurrent, iDirection) 
          CALL WriteData (Vscale, AexpV, Wtotal,SKINE,iDirection)
          CLOSE(1) 
          CLOSE(2) 
      EndDo 
      
      DEALLOCATE (FI)                                   
      EKIN = 0.5 * SKINE / AEXPV**2 
      Write (*, '('' Ekin='',E12.4,'' Weight per cell='',g12.5,'' (must be 1)'')') EKIN,Wtotal/float(NGRID)**3            
      Write (16,'('' Ekin='',E12.4,'' Weight per cell='',g12.5,'' (must be 1)'')') EKIN,Wtotal/float(NGRID)**3
                       
      CALL FilesWrite           ! write header and control data
      Write (25) astep/2**LevelParticles        ! write pt.dat file: time-step for particles             
      close(25)
      CALL Timing(0,1)
            write(*,'(10x,a/,(10x,a15,f10.2))') &
                     'Wall cpu_time (seconds) ', &
                     'Init        =',CPU(1), &
                     'Spectr      =',CPU(3), &
                     'Coords      =',CPU(4), &
                     'FFT         =',CPU(5), &
                     'Write temp  =',CPU(6), &
                     'Write final =',CPU(7), &
                     'Total       =',CPU(0)

      STOP 
      END  PROGRAM  PMstartMp                                         

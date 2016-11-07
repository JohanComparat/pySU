!-----------------------------------------------------
!     Module Tools for PMP2 code
!	 - Main shared variable and arrays
!        - Routines: ReadSetup, ReadDataPM,
!                  seconds, Timing, TimingMain, Memory
!-------------------------------------------------
!
Module Tools
        Real*4  :: Om,        &      ! cosmology
                 Omb,       &
                 OmL,       &
                 AEXPN,     &
                 ASTEP,     &
                 sigma8,    &
                 Box,       &
                 hubble
      Integer*4 :: NROW,    &      ! number of particles in 1D
                   NGRID,   &      ! number of grid points in 1D
                   Nout,    &      ! number of outputs/mocks
                   Nbiaspars       ! number of bias parameters
      Integer*8 :: Nparticles, Ngalaxies
      Real*4   :: zinit,da,zfinal
      Real*4   :: densThr,sigV
      Real*4   :: zout(1000),BiasPars(1000)

      Real*4     :: AEXPN0,ASTEP0,AMPLT,EKIN,EKIN1,EKIN2,AEU0,  &
                    TINTG,                 &
                    extras(100),ENKIN,ENPOT
      Integer*4  :: iStep,Nrealization,Nseed
      Real*4, Allocatable, Dimension(:,:,:) :: FI
      Real*4, allocatable,  Dimension(:) ::  Xb,Yb,Zb,VXb,Vyb,Vzb
      Real*4, allocatable,  Dimension(:) ::  XPAR,YPAR,ZPAR,VX,VY,VZ
      Real*4, allocatable,  Dimension(:) ::  dens,RandP
      Real*4       :: StepFactor =3.e-2
      Character*45 :: HEADER
      Real*8       :: CPU(0:10) =0.      
    Contains
!--------------------------------------------
subroutine ReadSetup
!--------------------------------------------
character*80 :: Line  
      read(11,'(a)') Line
          write(*,*)TRIM(Line)
      read(11,*) HEADER
          write(*,*) TRIM(HEADER)
      read(11,*) AEXPN0     !'Initial Expansion Parameter'
      write(*,*) ' Ainit=',AEXPN0
      
      read(11,*) ASTEP0   !'Initial Step in dAEXPN    da/a = ',da/AEXPN
      read(11,*) sigma     !'DRho/rho in box    '
      read(11,*) Box       !'Box in  Mpc/h   '
      read(11,*) sigma8    !'sigma8    '
      read(11,*) hubble    !'Hubble    '
      read(11,*) Om        !'Omega Matter'
      read(11,*) OmL       !'Omega Lambda'
      read(11,*) Omb       !'Omega Baryons'
      read(11,*) NROW      !'NROW  Number Particles   '
      read(11,*) NGRID     !'NGRID Number grid points '
      read(11,*) Nseed     !'Random seed'
      read(11,*) Cell      !'Cell Size   '
      read(11,*) aMass     ! 'Particle Mass'   
      read(11,*) zinit     ! 'Initial redshift       '
      read(11,*) zfinal    ! 'Final redshift       '
      read(11,*) DensThr   ! 'Density Threshold for V correction '
      read(11,*) sigV      ! 'rms V correction factor'
      read(11,*) Nout      ! 'Number of redshifts for analysis'
      read(11,*)(zout(i),i=1,Nout)
      read(11,*) Nbiaspars !  'Number of bias parameters'
      Do i=1,Nbiaspars
         read(11,*) BiasPars(i)   ! 'Bias parameter'
      End Do     
         write (*,*) ' Results were read from Setup.dat'
      CLOSE (11)
      AMPLT = sigma
    end subroutine ReadSetup

!---------------------------------------
!        Read    PMfiles
!             moment <0    use PMcrd.DAT, PMcrs0.DAT ...    
!             moment >= 0  use PMcrd.xxxx.DAT, PMcrs0,XXXX.DAT ..
    SUBROUTINE ReadDataPM(moment)
!      
!---------------------------------------
       Character*80 :: Name
        Logical      :: exst
        Integer*8    :: iCount,ii,ioff,ip
        Integer*8    :: Ngal,Nrecord,Jpage,NinPage
        Integer*4    :: moment
!
        !			Read data and open files
      If(moment<0)Then
         Open (4,file ='PMcrd.DAT',form ='UNFORMATTED',status ='UNKNOWN')
      Else
         write(Name,'(a,i4.4,a)')'PMcrd.',moment,'.DAT'
         Open (4,file =TRIM(Name),form ='UNFORMATTED',status ='UNKNOWN')
      end If
         
      READ  (4) HEADER,                        &
                       AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,     &
                       NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble, &
                       Nparticles,extras
      WRITE (*,'(a,/10x,a,f8.4,4(a,i7))') HEADER,                 & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
      WRITE (17,'(a,/10x,a,f8.4,4(a,i7))') HEADER,                 & 
                       ' a=',AEXPN, ' step= ',ISTEP,         &
                       ' Nrow= ', NROW, ' Ngrid=',NGRID
      close(4)

      Nrecord = 1024**2
      Naccess = Nrecord*6 !*4
      xR      = NGRID +1
      boxsize = extras(100)
      Box     = boxsize
      Npages   = (Nparticles-1)/Nrecord+1 ! number of records
      Nlast   = Nparticles - (Npages-1)*Nrecord ! number of particles in last record
      Jpage   = Nrecord

      write(*,'(a,i10)') ' NROW   =',NROW
      write(*,'(a,i10)') ' Ngal   =',Nparticles
      write(*,'(a,i10)') ' Ngrid  =',Ngrid
      write(*,'(a,f10.1)') ' Box    =',Box
      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))

      myMemory =Memory(6_8*Nparticles)
      Allocate (XPAR(Nparticles),YPAR(Nparticles),ZPAR(Nparticles))
      Allocate (VX(Nparticles),VY(Nparticles),VZ(Nparticles))

      iCount = 0

      ifile = 0
      jj    = 0
      If(moment<0)Then
         write(Name,'(a,i1.1,a)')'PMcrs',ifile,'.DAT'
      Else
         write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
      End If
        INQUIRE(file=TRIM(Name),EXIST=exst)
        If(.not.exst)Stop ' File PMcrs0.DAT does not exist. Error'
        OPEN(UNIT=20,FILE=TRIM(Name),ACCESS='DIRECT', &
              FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
   Do ii =1,Npages
      If(ii==Npages)Then
         NinPage = Nparticles -(ii-1)*JPAGE  ! # particles in the current page
      Else
         NinPage = JPAGE
      EndIf
      jj = jj +1
      If(ii<10.or.ii==Npages)write(*,'(3(a,i9))') ' Reading page= ',ii,' record =',jj,' NinPage= ',NinPage
10      Read(20,REC=jj,iostat=ierr) Xb,Yb,Zb,VXb,VYb,VZb
      If(ierr /=0)Then
         close(20)
         ifile = ifile +1
         If(Moment<0)Then
           If(ifile<10)Then
             write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.DAT'
           Else
             write(Name,'(a,i2.2,a,i4.4,a)')'PMcrs',ifile,'.DAT'
          EndIf
       Else
          If(ifile<10)Then
             write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
          Else
             write(Name,'(a,i2.2,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
          EndIf
       end If

           jj = 1
           INQUIRE(file=TRIM(Name),EXIST=exst)
           If(.not.exst)Then
             write(*,'(a,3i5)')' Attempting to read file number, record = ',ifile,ii,Npages
             write(*,'(2a)')' Attempting to read non-existing file: ',TRIM(Name)
             Stop ' Error reading PMcrs files: Did not get all files'
           End If
           Open(20,file=TRIM(Name),ACCESS='DIRECT', &
              FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
          write(*,'(2i7,2a,3x,i9)') ii,ifile,' Open file = ',TRIM(Name),Ninpage
          go to 10
       end If

       ioff = (ii-1)*JPAGE
!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE (ip)       
           Do ip =1,NinPage
              If(ip+ioff > Nparticles)STOP 'Attempt to read too many particles '
              if(INT(Xb(ip))==Ngrid+1)Xb(ip)=Xb(ip)-1.e-3
              if(INT(Yb(ip))==Ngrid+1)Yb(ip)=Yb(ip)-1.e-3
              if(INT(Zb(ip))==Ngrid+1)Zb(ip)=Zb(ip)-1.e-3
              if(INT(Xb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Xb(ip)),Xb(ip)
              if(INT(Yb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Yb(ip)),Yb(ip)
              if(INT(Zb(ip))==Ngrid+1)write(*,*)'Error in boundary: ',INT(Zb(ip)),Zb(ip)
              
                Xpar(ip+ioff) = Xb(ip) 
                Ypar(ip+ioff) = Yb(ip) 
                Zpar(ip+ioff) = Zb(ip) 
                  VX(ip+ioff) = VXb(ip) 
                  VY(ip+ioff) = VYb(ip) 
                  VZ(ip+ioff) = VZb(ip) 
           end Do
        end DO
        close (20)
      
           xx = MAXVAL(Xpar)
           xm = MINVAL(Xpar)
           write(*,*)' x     min/max= ',xm,xx
           xx = MAXVAL(Ypar)
           xm = MINVAL(Ypar)
           write(*,*)' y     min/max= ',xm,xx
           xx = MAXVAL(Zpar)
           xm = MINVAL(Zpar)
           write(*,*)' z     min/max= ',xm,xx
           Do ii =1,Nparticles
              If(Xpar(ii).lt.1.0.or.Ypar(ii).lt.1.0.or.Zpar(ii).lt.1.0)&
                write(*,'(a,i10,1p,3g14.5)')' Error coord: ',ii,Xpar(ii),Ypar(ii),Zpar(ii)
           EndDo
     DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)
     
    end SUBROUTINE ReadDataPM
!
!----------------------------------------------------
function seconds ()
!----------------------------------------------------
!
!     purpose: returns elapsed time in seconds
      Integer*8, SAVE :: first=0,rate=0,i0=0
      Integer*8       :: i

      If(first==0)Then
         CALL SYSTEM_CLOCK(i,rate)
         first =1
         i0    = i
         seconds = 0.
      Else
         CALL SYSTEM_CLOCK(i)
         seconds = float(i-i0)/float(rate)
      EndIf

    end function seconds
!--------------------------------------------
subroutine Timing ( ielement , isign )
!--------------------------------------------
        ! use Timing(i,-1) to start clock
        !     Timing(i,1)  to suspend clock
        !     Timing(0,0)  to print results
    Character*80 :: FName='timing.log'
    Integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
    If(isign == 0)Then
       Open(30,file=TRIM(FName),position='append')
       write(30,'(i5,1p,10G13.4)')CPU(0:10)/60.
       close(30)
       CPU(:) = 0
    EndIf
     CPU(ielement) = CPU(ielement) + float(isign) * seconds()
   end subroutine Timing
!--------------------------------------------
subroutine TimingMain ( ielement , isign )
!--------------------------------------------
        ! 0 - total
        ! 1 - force, 2- move
        ! 3 - density, 4- IO
    Character*80 :: FName='timing.log'
    Integer OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
    If(isign == 0)Then
       Open(30,file=TRIM(FName),position='append')
       If(iStep==1)Then
          write(30,'(3(a,i11))') ' Ngrid    = ',Ngrid, &
                                ' Npart    = ',Nparticles, &
                                ' Nthreads = ',OMP_GET_MAX_THREADS()
          write(30,'(T3,a,T10,a,T23,a,T35,a,T49,a,T62,a)')&
               'Step','Tot/min','Force','Move','Density','IO'
       End If
       write(30,'(i5,1p,10G13.4)'),iStep,CPU(0:4)/60.
       close(30)
       CPU(:) = 0
       return
    EndIf
     CPU(ielement) = CPU(ielement) + float(isign) * seconds()
   end subroutine TimingMain
!--------------------------------------------
Real*4 Function Memory(i_add)
!--------------------------------------------

    Real*8, SAVE :: mem =0.001   ! initial memory in Gb
    Integer*8 :: i_add           ! number of 4byte words
    Real*8    :: add
    mem    = mem + i_add*4./1024.**3
    Memory = mem
    write(*,'(a,f8.3,a)') ' Current Memory = ',Memory,'Gb'
  end Function Memory

end Module Tools

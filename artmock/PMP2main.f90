!-------------------------------------------------
!            Parallel (OMP) PM code
!	     2015, 1997  Anatoly Klypin (aklypin@nmsu.edu)
!	                           Astronomy Department, NMSU
!
!
!-------------------------------------------------
Module LocalData
  Integer*4  :: Nsteps,   &   ! # steps for this run
                Ntotal        ! # steps for the whole simulation
  Real*4     :: Alist(2000)   ! List of expansion parameter steps     
  Integer*4  :: Nlist(2000)   ! List of steps for analysis     
end Module LocalData
!
!-------------------------------------------------
Program PMP
  use Tools
  use fft5
  use Density
  use LocalData
  use Analyze
  
  Ncheckpoint  = 50   !  save for re-starting every  Ncheckpoint step

  Call Initialize     
                         
      DO  i=1,Nsteps         !  Main loop
         Call TimingMain(0,-1)
         CALL DENSIT            ! Define density
         CALL POTENTfft5            ! Define potential
         FactV  = 1.                ! normal step
         If(ASTEP<StepFactor/1.25*AEXPN.and.AEXPN<0.3333)Then
            ASTEP =1.5*ASTEP        ! increase step
            FactV = 1.2
         EndIf
         CALL MOVE(FactV)              ! Move particles
         CALL ADDTIME           ! Advance time

         If(mod(ISTEP,Ncheckpoint)==0) CALL WriteDataPM(0)  ! checkpoint snapshot
         If(Nlist(ISTEP)==1)Call Analysis
         Call TimingMain(0,1)
         Call TimingMain(0,0)
         IF(AEXPN.GE.1.)exit  ! Do again if a < 1
   end DO
   CALL WriteDataPM(1)     ! make snapshot
 END Program PMP

!--------------------------------------------
subroutine Initialize
!--------------------------------------------
  use Tools
  use LocalData
     logical :: exst

      OPEN(17,FILE='Result.log',STATUS='UNKNOWN',position='append')
      WRITE (*,'(A,$)') ' Enter number of steps for this run => '
      READ (*,*) Nsteps	 ! Make this number of steps

      Inquire(file='../Setup.dat',exist = exst)
      if(.not.exst)Then
         write(*,*)' Error: File ../Setup.dat not found. Run PMP2init.exe'
         stop
      end if
         open(11,file='../Setup.dat')

      Inquire(file='../TableSeeds.dat',exist = exst)
      if(.not.exst)Call SetSeeds
 
      Call ReadSetup
      CALL ReadDataPM(-1)
      myMemory =Memory(1_8*NGRID*NGRID*NGRID)
      Allocate (FI(NGRID,NGRID,NGRID))
      write(*,*) ASTEP0,AEXPN0
      StepFactor = ASTEP0/AEXPN0

      IF(AEXPN.GE.1.)THEN         ! change this if you need to run
	      WRITE (*,*) ' Cannot run over a=1' !  beyond a=1
	      STOP
      ENDIF
      write(*,'(a,T20,a,T30,i12,T45,a,i4))')' Start running: ', &
                      'Nparticles:',Nparticles,' Ngrid= ',Ngrid
      write(*,'(T20,a,T30,i12,T45,a,es12.4)')  'Step:',ISTEP,' da=',ASTEP
      write(*,'(T20,a,T30,es12.4)')  'StepFactor :',StepFactor

      !---- make table for steps
      Alist(:) = 0.
      Nlist(:) = 0
      da = ASTEP0
      a  = AEXPN0
      Alist(1) = a
      i= 0
      Do
         If(da<StepFactor/1.25*a.and.a<0.333)Then
            da =1.5*da        ! increase step
         EndIf
         a = a +da
         i = i +1
         If(i.gt.2000)Stop 'Too many timesteps.Increase length of Alist'
         Alist(i) = a
         If(a.ge.1.)exit
      end Do
      Ntotal =i
      Do i=1,Nout           !-- check every zout moment
         a = 1./(1.+zout(i))
         Do j=2,Ntotal      !-- find closest moment in all steps 
            if(a.lt.Alist(j))exit
         end Do
         da1 = Alist(j)-a
         da0 = a -Alist(j-1)
         If(da1<da0)Then    !-- mark closest time step for analysis
            Nlist(j) = 1
         else
            Nlist(j-1) =1
         end If
      end Do
      write(*,*)'Step  a_expansion   redshift   Analyze '
      do i=1,Ntotal
         if(Nlist(i)==1)&
         write(*,'(i5,2es13.4,i3)') i,Alist(i),1./Alist(i)-1.,Nlist(i)
      endDo
         
    end subroutine Initialize
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
!---------------------------------------
!        Write     PMfiles: PMcrd/crs   
 SUBROUTINE WriteDataPM(iFlag)
!
!---------------------------------------
    use Tools
       Character*80 :: Name
        Logical      :: exst
        Integer*8    :: iCount,ii,Jpage,Nrecord
        Integer*8    :: Ngal,moment,i,jfirst,jlast,j0
     Call TimingMain(4,-1)

      Npage   = 1024**2       ! number of particles per record
      Nrecord = Npage
      Jpage   = Npage
      Naccess = Nrecord*6 !!*4 !!!
      Nrecpage = 256            ! max number of records per file
      moment = ISTEP
      Npages  = (Nparticles-1)/JPAGE+1
      Nfiles  = (Npages-1)/Nrecpage +1    ! number of files

      write(*,*) ' Npages      = ',Npages
      write(*,*) ' Write files = ',iFlag,moment
        !			open files
     If(iFlag == 0)Then
        Open (4,file ='PMcrd.DAT',form ='UNFORMATTED',status ='UNKNOWN')
          write(Name,'(a,i1,a)')'PMcrs',0,'.DAT'
          OPEN(UNIT=20,FILE=TRIM(Name),ACCESS='DIRECT', &
                 FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
     Else
        write(Name,'(a,i4.4,a)')'PMcrd.',moment,'.DAT'
        Open (4,file =TRIM(Name),form ='UNFORMATTED',status ='UNKNOWN')
          write(Name,'(a,i1,a,i4.4,a)')'PMcrs',0,'.',moment,'.DAT'
          write(*,'(a)') TRIM(Name)
          OPEN(UNIT=20,FILE=TRIM(Name),ACCESS='DIRECT', &
                 FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
     end If
     PARTW = float(NGRID)**3/FLOAT(Nparticles)
     AU0   = 0.
      write  (4) HEADER,                        &
                       AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW, &
                       TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,     &
                       NROW,NGRID,Nrealization,Nseed,Om,OmL,hubble, &
                       Nparticles,extras
      close(4)
      myMemory =Memory(6_8*Nrecord)
      Allocate (Xb(Nrecord),Yb(Nrecord),Zb(Nrecord))
      Allocate (VXb(Nrecord),VYb(Nrecord),VZb(Nrecord))
      
        jj = 1
Do i=1,Npages         !-------- dump into  files
   If(i==Npages)Then
      NinPage = Nparticles -(i-1)*JPAGE  ! # particles in the current page
   Else
      NinPage = JPAGE
   EndIf
    If(mod(i-1,Nrecpage)==0.and.i>1)Then   ! close old and open new file
       close(20)
       jj = 1
       ifile = (i-1)/Nrecpage       ! construct file name
       If(iFlag==0)Then
          If(ifile<10)Then
             write(Name,'(a,i1.1,a)')'PMcrs',ifile,'.DAT'
          Else
             write(Name,'(a,i2.2,a)')'PMcrs',ifile,'.DAT'
          EndIf
       Else
          If(ifile<10)Then
             write(Name,'(a,i1.1,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
          Else
             write(Name,'(a,i2.2,a,i4.4,a)')'PMcrs',ifile,'.',moment,'.DAT'
          EndIf
       end If
       Open(20,file=TRIM(Name),ACCESS='DIRECT', &
                 FORM='unformatted',STATUS='UNKNOWN',RECL=NACCESS)
       write(*,'(2i7,2a,3x,i9)') i,ifile,' Open file = ',TRIM(Name),Ninpage
    end If

   jfirst  = (i-1)*JPAGE +1        ! first and last particles in current record
   jlast   = jfirst + NinPage-1
   If(mod(i,10)==0.or.i==Npages)    &
        write(*,'(10x,a,i5,a,4i11)')'Write page=',i,' Particles=',NinPage,jfirst,jlast

   Do j0 = jfirst,jlast
      Xb(j0-jfirst+1)   = Xpar(j0)
      Yb(j0-jfirst+1)   = Ypar(j0)
      Zb(j0-jfirst+1)   = Zpar(j0)
      Vxb(j0-jfirst+1)  = VX(j0)
      VYb(j0-jfirst+1)  = VY(j0)
      VZb(j0-jfirst+1)  = VZ(j0)
   EndDo

   !write(*,'(i10,1p,6g13.5)') (k,XPAR(k),YPAR(k),ZPAR(k),VX(k),VY(k),VZ(k),k=1,1024)
   WRITE (20,REC=jj) Xb,Yb,Zb,VXb,VYb,VZb
   jj = jj +1
EndDo                            ! end lspecies loop   
        
     myMemory =Memory(-6_8*Nrecord)
     DEALLOCATE (Xb,Yb,Zb,VXb,VYb,VZb)
     Call TimingMain(4,1)
     
   end SUBROUTINE WriteDataPM
!--------------------------------------------------
      SUBROUTINE SetTest
!--------------------------------------------------
use Tools
     Integer*8 :: ii

     ISTEP =0
     AEXPN = 0.01
     ASTEP = 0.001
     Om   = 1.
     Oml   = 0.
     write(*,*) ' Ngrid =',Ngrid
     ii = 0
     DO M3=1,NGRID,2
       DO M2=1,NGRID,2
          DO M1=1,NGRID,2
             ii = ii +1
             If(ii>Nparticles)Stop ' Number of particles is to big in SetTest'
             XPAR(ii) =  M1+0.05
             YPAR(ii) =  M2+0.05
             ZPAR(ii) =  M3+0.05

             VX(ii)   = 0.
             VY(ii)   = 0.
             VZ(ii)   = 0.
	  END DO
       END DO
    END DO
    If(ii/= Nparticles)Stop ' Wrong number of particles in SetTest'
end SUBROUTINE SetTest
!--------------------------------------------------
!	   Advance Aexpn, Istep, tIntg,...
!	   AEXPN = currnet expansion parameter
!          ISTEP = current step
!          ASTEP = step in the expansion parameter	    
      SUBROUTINE ADDTIME
!--------------------------------------------------
use Tools

      ISTEP = ISTEP + 1
      AEXPN = AEXPN + ASTEP
		                 !    Energy conservation
      IF(ISTEP.EQ.1)THEN
        EKIN1 = EKIN
        EKIN2 = 0.
        EKIN  = ENKIN
        AU0   = AEXP0*ENPOT
        AEU0  = AEXP0*ENPOT + AEXP0*(EKIN+EKIN1)/2.
        TINTG = 0.
        WRITE (*,40) ISTEP,AEXPN,EKIN,ENPOT,AU0,AEU0
        WRITE (17,40) ISTEP,AEXPN,EKIN,ENPOT,AU0,AEU0
40      FORMAT('**** STEP=',I3,' A=',F10.4,' E KIN=',E12.4,  &
              ' E POT=',E12.4,/'      AU0,AEU0=',2E12.4)
      ELSE
        EKIN2 = EKIN1
        EKIN1 = EKIN
        EKIN  = ENKIN
        TINTG = TINTG +ASTEP*(EKIN1 +(EKIN -2.*EKIN1 +EKIN2)/24.)
        ERROR = ((AEXPN-ASTEP)*((EKIN+EKIN1)/2.+ENPOT)-AEU0+TINTG)/ &
                       ((AEXPN-ASTEP)*ENPOT)*100.
        WRITE (*,50)  ISTEP,AEXPN,ERROR,EKIN,ENPOT,TINTG
        WRITE (17,50) ISTEP,AEXPN,ERROR,EKIN,ENPOT,TINTG
      ENDIF
50    FORMAT('Step = ',I4,' A=',F7.4,' Error(%)=',f7.2, &
             ' Ekin=',E11.3,' Epot=',E11.3,' Intg=',E11.3)
    END SUBROUTINE ADDTIME
!------------------------------------------------------------

!     Advance each particle:	    dA	  AEXPN     dA
!	   by one step	     I______._______I_______.______I	  ->  A
!			    i-1     .	    i	    .	  i+1	step
!				    .	  {Fi}	    .	   .
!				  { vx }  { x }     .	   .
!				  { vy }  { y }     ^	   .
!				    ._______________^	   .
!					    .		   ^
!					    .______________^
!
!			  0.5
!		  dP = - A     * Grad(Fi) * dA ; A =AEXPN
!			  i		i
!				   3/2
!		  dX =	 P(new)/A	* dA ; A      =AEXPN+dA/2
!				 i+1/2		i+1/2
!
!------------------------------------------------
      SUBROUTINE MOVE(FactV)
!------------------------------------------------
use Tools
!             PCONST = factor to change velocities
!             XCONST = factor to change coordinates
!	         Note: 0.5 is needed in Pconst because
!	            Fi(i+1)-Fi(i-1) is used as gradient
!             FactV = 1 - normal constant step
!                   = 1.2 - increase step by 1.5
real*8 :: SVEL,SPHI, PCONST, XCONST, XN, YN, &
     D1,D2,D3,T1,T2,T3, T2T1,T2D1,D2T1,D2D1, &
     GX,GY,GZ,FP,VVx,VVY,VVZ,                &
     GX111,GX211,GX121,GX221,                &
     GX112,GX212,GX122,GX222,                &
     GY111,GY211,GY121,GY221,                &
     GY112,GY212,GY122,GY222,                &
     GZ111,GZ211,GZ121,GZ221,                &
     GZ112,GZ212,GZ122,GZ222,                &
     X,Y,Z
integer*8 :: IN
     Call TimingMain(2,-1)
      PCONST = - SQRT(AEXPN/(Om+OmL*AEXPN**3))*ASTEP*0.5/FactV
      Ahalf  =   AEXPN+ASTEP/2.
      XCONST =   ASTEP/SQRT(Ahalf*(Om+OmL*Ahalf**3))/Ahalf
	   SVEL = 0.                      ! counter for \Sum(v_i**2)
	   SPHI = 0.                      ! counter for \Sum(phi_i)
	   XN   = FLOAT(NGRID)+1.-1.E-8   ! N+1
	   YN   = FLOAT(NGRID)            ! N
    Wpar    = YN**3/FLOAT(Nparticles)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE (X,Y,Z,GX,GY,GZ,FP,VVx,VVY,VVZ,I,J,K) &
!$OMP PRIVATE (I1,J1,K1,K2,K3,I2,J2,I0,J0) &
!$OMP PRIVATE (D1,D2,D3,T1,T2,T3, T2T1,T2D1,D2T1,D2D1) &
!$OMP PRIVATE (F111,F211,F121,F221,F112,F212,F122,F222) &
!$OMP PRIVATE (F113,F213,F123,F223,F110,F210,F120,F220) &
!$OMP PRIVATE (F311,F321,F131,F231,F312,F322,F132,F232) &
!$OMP PRIVATE (F011,F021,F101,F201,F012,F022,F102,F202) &
!$OMP PRIVATE (GX111,GX211,GX121,GX221,GX112,GX212,GX122,GX222) &
!$OMP PRIVATE (GY111,GY211,GY121,GY221,GY112,GY212,GY122,GY222) &
!$OMP PRIVATE (GZ111,GZ211,GZ121,GZ221,GZ112,GZ212,GZ122,GZ222) &
!$OMP REDUCTION(+:SVEL,SPHI)    
	    DO  IN=1,Nparticles            ! Loop over particles
  	       X=XPAR(IN)
	       Y=YPAR(IN)
	       Z=ZPAR(IN)
	       VVX=VX(IN)
	       VVY=VY(IN)
	       VVZ=VZ(IN)
	       I=INT(X)
	       J=INT(Y)
	       K=INT(Z)
	       D1=X-FLOAT(I)
	       D2=Y-FLOAT(J)
	       D3=Z-FLOAT(K)
	       T1=1.-D1
	       T2=1.-D2
	       T3=1.-D3
	       T2T1 =T2*T1
	       T2D1 =T2*D1
	       D2T1 =D2*T1
	       D2D1 =D2*D1
	       I1=I+1
	          IF(I1.GT.NGRID)I1=1
	       J1=J+1
	          IF(J1.GT.NGRID)J1=1
	       K1=K+1
	          IF(K1.GT.NGRID)K1=1
	       K2=K+2
	          IF(K2.GT.NGRID)K2=K2-NGRID
	       K3=K-1
	          IF(K3.LT.1    )K3=NGRID
	       F111 =FI(I ,J ,K )  !  Read potential to Fij vars
	       F211 =FI(I1,J ,K )
	       F121 =FI(I ,J1,K )
	       F221 =FI(I1,J1,K )
   
	       F112 =FI(I ,J ,K1)
	       F212 =FI(I1,J ,K1)
	       F122 =FI(I ,J1,K1)
	       F222 =FI(I1,J1,K1)
   
	       F113 =FI(I ,J ,K2)
	       F213 =FI(I1,J ,K2)
	       F123 =FI(I ,J1,K2)
	       F223 =FI(I1,J1,K2)
   
	       F110 =FI(I ,J ,K3)
	       F210 =FI(I1,J ,K3)
	       F120 =FI(I ,J1,K3)
	       F220 =FI(I1,J1,K3)
   
	       I2=I+2
	          IF(I2.GT.NGRID)I2=I2-NGRID
	       J2=J+2
	          IF(J2.GT.NGRID)J2=J2-NGRID
	       F311 =FI(I2,J ,K )
	       F321 =FI(I2,J1,K )
	       F131 =FI(I ,J2,K )
	       F231 =FI(I1,J2,K )
   
	       F312 =FI(I2,J ,K1)
	       F322 =FI(I2,J1,K1)
	       F132 =FI(I ,J2,K1)
	       F232 =FI(I1,J2,K1)
   
	       I0=I-1
	          IF(I0.LT.1)I0=NGRID
	       J0=J-1
	          IF(J0.LT.1)J0=NGRID
	       F011 =FI(I0,J ,K )
	       F021 =FI(I0,J1,K )
	       F101 =FI(I ,J0,K )
	       F201 =FI(I1,J0,K )
   
	       F012 =FI(I0,J ,K1)
	       F022 =FI(I0,J1,K1)
	       F102 =FI(I ,J0,K1)
	       F202 =FI(I1,J0,K1)
			!	 Find {2*gradient} in nods
	       GX111 =F211 -F011
	       GX211 =F311 -F111
	       GX121 =F221 -F021
	       GX221 =F321 -F121
       
	       GX112 =F212 -F012
	       GX212 =F312 -F112
	       GX122 =F222 -F022
	       GX222 =F322 -F122
       
	       GY111 =F121 -F101
	       GY211 =F221 -F201
	       GY121 =F131 -F111
	       GY221 =F231 -F211
       
	       GY112 =F122 -F102
	       GY212 =F222 -F202
	       GY122 =F132 -F112
	       GY222 =F232 -F212
       
	       GZ111 =F112 -F110
	       GZ211 =F212 -F210
	       GZ121 =F122 -F120
	       GZ221 =F222 -F220
       
	       GZ112 =F113 -F111
	       GZ212 =F213 -F211
	       GZ122 =F123 -F121
	       GZ222 =F223 -F221
			    !	 Interpolate to the point
      GX=PCONST*(T3*(T2T1*GX111+T2D1*GX211 +D2T1*GX121+D2D1*GX221 )+ &
     		 D3*(T2T1*GX112+T2D1*GX212 +D2T1*GX122+D2D1*GX222 ))

      GY=PCONST*(T3*(T2T1*GY111+T2D1*GY211 +D2T1*GY121+D2D1*GY221 )+ &
     		 D3*(T2T1*GY112+T2D1*GY212 +D2T1*GY122+D2D1*GY222 ))

      GZ=PCONST*(T3*(T2T1*GZ111+T2D1*GZ211 +D2T1*GZ121+D2D1*GZ221 )+ &
      		 D3*(T2T1*GZ112+T2D1*GZ212 +D2T1*GZ122+D2D1*GZ222 ))

			!	 Find potential of the point
      FP=	 T3*(T2T1*F111+T2D1*F211 +D2T1*F121+D2D1*F221 )+  &
      		 D3*(T2T1*F112+T2D1*F212 +D2T1*F122+D2D1*F222 )
	 SPHI = SPHI + FP*WPAR

         VVX =VVX+GX             ! Move points
	 VVY =VVY+GY
	 VVZ =VVZ+GZ
	 X	=X  +VVX*XCONST
	 Y	=Y  +VVY*XCONST
	 Z	=Z  +VVZ*XCONST
	 IF(X.LT.1.d0)X=X+YN     ! Periodical conditions
	 IF(X.GE.XN)X=X-YN
	 IF(Y.LT.1.d0)Y=Y+YN
	 IF(Y.GE.XN)Y=Y-YN
	 IF(Z.LT.1.d0)Z=Z+YN
         IF(Z.GE.XN)Z=Z-YN

         
	 SVEL=SVEL+(VVX**2+VVY**2+VVZ**2)*WPAR
	     	                                        
	  XPAR(IN)=X            ! Write new coordinates
	  YPAR(IN)=Y
	  ZPAR(IN)=Z
	  VX(IN)=VVX
	  VY(IN)=VVY
	  VZ(IN)=VVZ
           if(INT(Xpar(IN))==Ngrid+1)Xpar(IN)=Xpar(IN)-1.e-3
           if(INT(Ypar(IN))==Ngrid+1)Ypar(IN)=Ypar(IN)-1.e-3
           if(INT(Zpar(IN))==Ngrid+1)Zpar(IN)=Zpar(IN)-1.e-3
      ENDDO
!			   Set energies:
!			   Kin energy now at A(i+1/2)
!			   Pot energy		at A(i)
      ENKIN = SVEL / 2. / (AEXPN+ASTEP/2.)**2
      ENPOT = SPHI /2.
      Call TimingMain(2,1)

  end SUBROUTINE MOVE
!-------------------------------------------------
!	    Find potential on Grid FI:	DENSITY    ->	POTENTIAL
!
!		   O 1		    ^ - Fourier component
!		   |
!	     1	   |-4	 1	^      ^	2Pi
!	     O-----O-----O     Fi    =	Rho	/ (2cos(---  (i-1))+
!		   |	   i,j		i,j	Ngrid
!		   |
!		   O 1			  2Pi
!				       2cos(---  (j-1))-4)
!		       ^			Ngrid
!		       Fi	= 1 (?)
!			 11
!		   2
!		NABLA  Fi = 3/2  /A * (Rho - <Rho>) ;
!		   X
!			      <Rho> = 1
!
SUBROUTINE POTENTfft5
!---------------------------------------------
  use Tools
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
      Real*8  :: GREENf(Nlenwrk)
      real*8  :: r(Nlenwrk)
!$OMP THREADPRIVATE(work,wsave)

      Call TimingMain(1,-1)
      If(Ngrid>Nlenwrk)Stop ' Incresase Nlenwrk in POTENTfft5'
      Ng = Ngrid
      lensav = Ngrid+int(log(real(Ngrid,kind = 4))/log(2.0E+00))+4
      lenwrk = Ngrid

      trfi = 1.5*Om/aexpn
			 ! Set Green function components
     GREENf(1)     =  2.
     GREENf(Ngrid) = -2.
      DO i=1,Ngrid/2-1
	  XX = 2.*COS(P16*i/Ngrid)
	  GREENf(2*i)   = XX
	  GREENf(2*i+1) = XX
      End DO

      call rfft1i ( Ng, wsave, lensav, ier ) !   Initialize FFT
      inc  = 1
      lenr = Ngrid

      !write(*,'(a,3i8,3x,1p,10G15.7)')' PotentFFT5: ',Ngrid,lensav,lenwrk,(GREENf(i),i=1,10)
       write(*,*) ' time =',seconds(), ' XY fft'
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k=1,NGRID             ! fft for xy planes in x-dir
       Do j=1,NGRID     
          Do i=1,NGRID
             r(i) = FI(i,j,k)
          EndDo
          call rfft1f ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do i=1,NGRID
            FI(i,j,k) = r(i)
          EndDo
       EndDo

       Do i=1,NGRID        ! fft xy planes in y-dir
          Do j=1,NGRID
             r(j) = FI(i,j,k)
          EndDo
          call rfft1f ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do j=1,NGRID
            FI(i,j,k) = r(j)
          EndDo
       EndDo

    EndDo

    
       write(*,*) ' time =',seconds(), '  transposition'
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
       write(*,*) ' time =',seconds(), ' Z forw/backw'
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r, ier,A1,A2,A3)
     Do j=1,NGRID     ! ------ z-direction
       Do i=1,NGRID     
           Do k=1,NGRID
             r(k) = FI(k,j,i)
          EndDo
          call rfft1f ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
           Do k=1,NGRID
             FI(k,j,i) = r(k)
          EndDo
       EndDo   
    end Do
    ww = (P16/Ngrid)**2
!$OMP PARALLEL DO DEFAULT(SHARED)  copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r, ier,A1,A2,A3,wi,wj,wk)
     Do j=1,NGRID     ! ------ z-direction
        A3 = GREENf(J) -6.
        !wj = j/2
       Do i=1,NGRID     
          A2 = GREENf(I) + A3
          !wi = i/2 
          Do k=1,NGRID
             A1 =A2 +GREENf(K)               !--- use this for descrete Poisson solver 
             ! wk = k/2
             !A1 = -ww*(wi**2+wj**2+wk**2)   !--- use this for k**2 Green functions
             !if(ww<5.)write(50,'(3i4,3x,3f7.3,f9.4,2g14.5)') i,j,k,wi,wj,wk,ww,A1,FI(k,j,i)
             IF(ABS(A1).LT.1.d-7) A1=1.
             r(k) = FI(k,j,i)*trfi/A1
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do k=1,NGRID
            FI(k,j,i) = r(k)
          EndDo
       EndDo
    end Do


       write(*,*) ' time =',seconds(),'  transpose '
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
       write(*,*) ' time =',seconds(), ' XY fft'
!$OMP PARALLEL DO DEFAULT(SHARED) copyin(wsave,work) & 
!$OMP PRIVATE ( k,j,i ,r,ier)
    Do k=1,NGRID             ! fft for xy planes in x-dir
       Do j=1,NGRID     
          Do i=1,NGRID
             r(i) = FI(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do i=1,NGRID
            FI(i,j,k) = r(i)
          EndDo
       EndDo

       Do i=1,NGRID        ! fft xy planes in y-dir
          Do j=1,NGRID
             r(j) = FI(i,j,k)
          EndDo
          call rfft1b ( Ng, inc, r, lenr, wsave, lensav, work, lenwrk, ier )
          Do j=1,NGRID
            FI(i,j,k) = r(j)
          EndDo
       EndDo

    EndDo
       write(*,*) ' time =',seconds(), ' Finished Potent'

      Call TimingMain(1,1)     
            
    end SUBROUTINE POTENTfft5





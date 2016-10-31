!--------------------------------------------------
!
!      Routines to analyze results on PM simulations
!
!--------------------------------------------------
Module Analyze

CONTAINS

!--------------------------------------------
!
subroutine Analysis
!
!--------------------------------------------
use Tools
use Density
use LUXURY
use Random
use Power
character*120 :: Name


           !-----  Fill the table for things to analyze
    iPower         = 1   !-- all particles Power spectrum 
    iPowerRSD      = 1   !       redshift distortions for DM
    iDensityDistr  = 1   !-- PDF for DM for different cell-sizes
    iBias          = 1   !-- Biasing model
    iSave          = 0   !-- Save snapshot
    
       moment = 100.*(1./AEXPN-1.) ! redshift*100.
       moment = max(moment,0)
       Nseed = 13979821
    
      Call DENSIT   ! make dinsity on original Ngrid mesh

               !------------ Density Distribution in real space
    If(iDensityDistr == 1)Then
       write(Name,'(2(a,i4.4),3(a,i3.3))')'DensDistrDM.',moment,   &
                                        '.',Nrealization,'.dat'
       OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(18,'(a,f7.4,a,i4,a,f8.3,a,i4,a,f8.2)')' Aexpn =',AEXPN,' Step=',ISTEP,    &
            ' Redshift= ',1./AEXPN-1.,' Ngrid= ',Ngrid
      Call DensDistr
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID_old = NGRID                  ! store old value of NGRID
       NGRID = NGRID/2                    ! increase cell-size twice
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      Call DENSITrsd(0,NGRID_old)         !- density in real space with new Ngrid
      Call DensDistr                      !  statistics of PDF
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID = NGRID_old/4
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      Call DENSITrsd(0,NGRID_old)
      Call DensDistr
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID = NGRID_old/8
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      Call DENSITrsd(0,NGRID_old)
      Call DensDistr
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       Deallocate (FI)
       NGRID = NGRID_old/16
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      Call DENSITrsd(0,NGRID_old)
      Call DensDistr
      close(18)
                                         ! restore NGRID and Density FI
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       DEALLOCATE(FI)
      NGRID = NGRID_old
        Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
      CALL DENSIT
    End If
                      !----------------- Power spectrum
    If(iPower == 1)Then
       If(iPowerRSD == 1)Then
         Mem_current = Memory(1_8*Nparticles)
         ALLOCATE(dens(Nparticles))
         CALL DENSPART     ! get density for each particle. Need it for corrections
       end If
       Call GetPower(0)  ! pk only real space, Ngrid mesh
       If(iPowerRSD == 1)Then
         CALL SetRandom    ! set seeds for random numbers
         Call GetPower(1)  ! pk in real and redshift space, Ngrid/2 mesh
       end If
    End If
                      !----------------- Save snapshot    
    If(iSave == 1)Call WriteDataPM(0)
    
                      !----------------- Biasing model    
    If(iBias == 1)Then
      If(iPower /= 1)Then     ! get density for each particle
        Mem_current = Memory(1_8*Nparticles)
        ALLOCATE(dens(Nparticles))
        CALL DENSPART             ! find density for each particle
        CALL SetRandom            ! set seeds for random numbers
      EndIf    
       Mem_current = Memory(-1_8*NGRID*NGRID*NGRID)
       DEALLOCATE(FI)             ! release FI to make room for other arrays
       
       Call BiasParticles         ! select galaxies using biasing model
                                  ! write galaxies into a file
       If(Ngalaxies>1)Then
          Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
          ALLOCATE(FI(NGRID,NGRID,NGRID))
          write(Name,'(2(a,i4.4),3(a,i3.3))')'DensDistrGal.',moment,   &
                                           '.',Nrealization,'.dat'
          OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
          write(18,'(a)')HEADER
          write(18,'(a,i10,a,f8.3,a,es13.4,a,f8.3,a,f8.3)')' Ngalaxies =',Ngalaxies, &
               ' DensThresh= ',BiasPars(1),   &
               ' Normalize = ',BiasPars(2),   &
               ' Slope     = ',BiasPars(3),   &
               ' DensMax   = ',BiasPars(4)
          Call  DensGal(0)
          Call  DensDistr
          close(18)
          Call GetPowerGal(0)  ! pk only real space, Ngrid mesh
          Call GetPowerGal(1)  ! pk real+z space, Ngrid/2 mesh

           myMemory =Memory(-7_8*Ngalaxies)
          DeAllocate (Xb,Yb,Zb)
          DeAllocate (VXb,VYb,VZb,RandP)
        end If  ! end Ngalaxies >1
    end If      ! end Biasing model
                         !--- deallocate arrays which are not needed
    If(ALLOCATED(dens))Then
            myMemory =Memory(-1_8*Nparticles)
            DeAllocate (dens)
    End If
    If(ALLOCATED(SeedsPage))Then
            DeAllocate (SeedsPage)
         End If
    If(.not.ALLOCATED(FI))Then 
       Mem_current = Memory(1_8*NGRID*NGRID*NGRID)
       ALLOCATE(FI(NGRID,NGRID,NGRID))
    end If
         
  end subroutine Analysis

!--------------------------------------------
!
!            create a set of biased particles
!
subroutine BiasParticles
!
!           BiasPars(1) = density threshold
!                  (2) = normalization
!                  (3) = power law slope
!--------------------------------------------
use Tools
use LUXURY
use Random
character*120 :: Name
Real*4, allocatable :: Gg(:),SelGal(:)
Integer*8 :: ip,jp,kp,jpp,kpp,IN,iGal

         write(*,*) ' Inside BiasParticles   Allocate Gg Nrow= ',Nrow
       moment = 100.*(1./AEXPN-1.) ! redshift*100.
       moment = max(moment,0)
      myMemory =Memory(1_8*Nparticles)
      Allocate (SelGal(Nparticles))
      Allocate(Gg(Nrow))
      
         Gg = 0.
         SelGal(:) = 0.
         iGal =0
!!$OMP PARALLEL DO DEFAULT(SHARED) & 
!!$OMP PRIVATE (IN,ip,jpp,kpp,jp,kp,Gg,dd,ff) &
!!$OMP REDUCTION(+:iGal)
      Do kp =1,NROW
         if(mod(kp,50)==0.and.iSwitch==1)write(*,*) ' k=',kp
         kpp = (kp-1)*NROW**2
      Do jp =1,NROW
         jpp = (jp-1)*NROW
!!$OMP critical            
            Call getRandom(Gg,jp,kp,NROW)
!!$OMP end critical
      Do ip =1,NROW
         IN = ip + jpp +kpp
                        !---- Biasing model
         dd = dens(IN)+1.-BiasPars(1)
         If(dd>0.0 .AND. dd<BiasPars(4))Then             ! take only above density threshold
            
            ff = BiasPars(2)*dd**BiasPars(3)
            
            If(ff>1.)write(*,'(2(a,es13.4))') ' Too big bias: probability= ',ff, &
                    ' DM density= ',dens(IN)
            If(Gg(ip)<ff)Then      ! mark particle as galaxy
               iGal = iGal+1
               !write(*,'(a,i12,a,4es13.4,a,8i10)') '    particle =',IN, &
               !     ' dens=',dens(IN),dd,ff,Gg(ip),' i=',ip,jp,kp,iGal
               SelGal(IN) = dens(IN)
            End If
         end If
      end Do
      end Do
      end Do

      Ngalaxies = iGal
      write(*,*)' Selected Galaxies=',Ngalaxies
      If(Ngalaxies == 0)Then
         write(*,*) '  No selected galaxies. Skip analysis.'
         DeAllocate (SelGal)
         deallocate(Gg)
         return
      End If
      
      myMemory =Memory(7_8*Ngalaxies)
      Allocate (Xb(Ngalaxies),Yb(Ngalaxies),Zb(Ngalaxies))
      Allocate (VXb(Ngalaxies),VYb(Ngalaxies),VZb(Ngalaxies))
      Allocate (RandP(Ngalaxies))
      
      ip = 0
      Xscale = Box/NGRID            ! scale coords to Mpch
      Vscale = Xscale*100./AEXPN    ! scale V to km/s
      Nseed = 1287922+3177*Nrealization
      ip = 0
      Do IN =1,Nparticles           !--- put selected into buffers
         If(SelGal(IN)> 0.)Then
            ip =ip +1
            dd = 1.
            If(dens(IN)+1.-DensThr>0.)dd = sigv*(dens(IN)+1.-DensThr)**0.3333 ! correct for resolution
            RandP(ip) = dens(IN)                 ! DM density at galaxy
            Xb(ip) = (XPAR(IN)-1.)*Xscale        ! coords in Mpch units
            Yb(ip) = (YPAR(IN)-1.)*Xscale
            Zb(ip) = (ZPAR(IN)-1.)*Xscale
            VXb(ip) = VX(IN)*Vscale + dd*GAUSS(Nseed) !velocity in km/s
            VYb(ip) = VY(IN)*Vscale + dd*GAUSS(Nseed)
            VZb(ip) = VZ(IN)*Vscale + dd*GAUSS(Nseed)
         End If
      EndDo

!!!!------------- random distribution !!!!!
!      Do IN =1,Nparticles         
!         Xb(IN) = RANDd(Nseed)*Box
!      End Do
!      Do IN =1,Nparticles         
!         Yb(IN) = RANDd(Nseed)*Box
!      End Do
!      Do IN =1,Nparticles         
!         Zb(IN) = RANDd(Nseed)*Box
!      End Do
     
      If(ip/=Ngalaxies)Stop ' Error in number of galaxies'

       write(Name,'(2(a,i4.4),3(a,i3.3))')'GalaxiesZ.',moment,   &
                                        '.',Nrealization,'.dat'
       OPEN(18,FILE=TRIM(Name),STATUS='UNKNOWN')
       write(18,'(a)')HEADER
       write(Name,'(2(a,i4.4),3(a,i3.3))')'DensityZ.',moment,   &
                                        '.',Nrealization,'.dat'
       write(18,'(a,i10,a,f8.3,a,es13.4,a,f8.3)')' Ngalaxies=',Ngalaxies, &
            ' DensThresh= ',BiasPars(1),    &
            ' Normalize= ',BiasPars(2),   &
            ' Slope= ',BiasPars(3)
       write(18,'(T10,a,T43,a,T67,a)')'XYZ(Mpch)','Vxyz(km/s)','DensDM rho/<rho>'
       Do i=1,Ngalaxies
          write(18,'(3f11.5,3x,3f10.3,es12.3)') Xb(i),Yb(i),Zb(i), &
                        VXb(i),VYb(i),VZb(i),RandP(i)+1.
       End Do
       close(18)

      myMemory =Memory(-1_8*Nparticles)
      
      DeAllocate (SelGal)
      deallocate(Gg)
      
end subroutine BiasParticles


end Module Analyze

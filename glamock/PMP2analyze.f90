!-------------------------------------------------
!            Parallel (OMP) PM code
!	     2015, 1997  Anatoly Klypin (aklypin@nmsu.edu)
!	                           Astronomy Department, NMSU
!
!
!-------------------------------------------------
Module LocalData
  Integer*4 :: moment     ! Snapshot number
end Module LocalData
!
!-------------------------------------------------
Program Bias
  use Tools
  use fft5
  use Density
  use LocalData
  use Analyze

  Call Initialize     
  Call Analysis
  
END Program Bias

!--------------------------------------------
subroutine Initialize
!--------------------------------------------
use Tools
use LocalData
     logical :: exst

      WRITE (*,'(A,$)') ' Enter snapshot number to analyze => '
      READ (*,*) moment	 ! read snapshot number
      Inquire(file='../Setup.dat',exist = exst)
      if(.not.exst)Then
         write(*,*)' Error: File ../Setup.dat not found. Run PMP2init.exe'
         stop
      end if
         open(11,file='../Setup.dat')

      
      Call ReadSetup
      CALL ReadDataPM(moment)
      myMemory =Memory(1_8*NGRID*NGRID*NGRID)
      Allocate (FI(NGRID,NGRID,NGRID))
         
end subroutine Initialize





      PARAMETER ( maxlev = 4, lxdflt = 3, jsdflt = 314159265 )
      PARAMETER ( igiga     = 1000000000)
      PARAMETER ( twop12 = 4096.)
      PARAMETER ( itwo24 = 2**24, icons = 2147483563 )
      PARAMETER (NROWmax = 8192)
      LOGICAL           notyet
!                            default
!  Luxury Level     0   1   2  *3*    4
!    ndskip        /0, 24, 73, 199, 365/
! Corresponds to p=24  48  97  223  389
!     time factor   1   2   3    6   10   on slow workstation
!                   1 1.5   2    3    5   on fast mainframe
!                   1 1.5 2.5    5  8.5   on PC using LF90
      COMMON /LUX1/   &
             i24, j24, in24, kount, mkount, carry, seeds(24) 
      COMMON /LUX2/ notyet,  twom24, twom12, &
             nskip, ndskip(0:maxlev), next(24), inseed, luxlev
      COMMON/AllSeeds/i24A(NROWmax),j24A(NROWmax),in24A(NROWmax),        &
                           kountA(NROWmax), carryA(NROWmax),          &
	                   gSetA(NROWmax),  iFlagA(NROWmax),          &
                           SeedsPage(24,NROWmax)
!$OMP threadprivate(/LUX1/)                         


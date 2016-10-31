!------------------------------------------------------------------
!
!               FFT5 routines
!
MODULE FFT5
Contains
!----------------------------
subroutine rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1B: real single precision backward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1B computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the backward
!    transform or Fourier synthesis, transforming the sequence from 
!    spectral to physical space.  This transform is normalized since a 
!    call to RFFT1B followed by a call to RFFT1F (or vice-versa) reproduces
!    the original array within roundoff error. 
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    25 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, 
!    in array R, of two consecutive elements within the sequence. 
!
!    Input/output, real ( kind = 4 ) R(LENR), on input, the data to be 
!    transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1) + 1. 
!
!    Input, real ( kind = 4 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFT1I before the first call to routine
!    RFFT1F or RFFT1B for a given transform length N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4. 
!
!    Workspace, real ( kind = 4 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least N. 
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = 8 ) r(lenr)
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if ( lenr < inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'rfft1b ', 6 )
    return
  end if

  if ( lensav < n + int ( log ( real ( n, kind = 4 ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfft1b ', 8 )
    return
  end if

  if ( lenwrk < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RFFT1B - Fatal error!'
    write ( *, '(a)' ) '  LENWRK < N:'
    write ( *, '(a,i6)' ) '  LENWRK = ', lenwrk
    write ( *, '(a,i6)' ) '  N = ', n
    ier = 3
    call xerfft ( 'rfft1b ', 10 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call rfftb1 ( n, inc, r, work, wsave, wsave(n+1) )

  return
end subroutine rfft1b
subroutine rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! RFFT1F: real single precision forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    RFFT1F computes the one-dimensional Fourier transform of a periodic
!    sequence within a real array.  This is referred to as the forward 
!    transform or Fourier analysis, transforming the sequence from physical
!    to spectral space.  This transform is normalized since a call to 
!    RFFT1F followed by a call to RFFT1B (or vice-versa) reproduces the 
!    original array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    25 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array R, of two consecutive elements within the sequence. 
!
!    Input/output, real ( kind = 8 ) R(LENR), on input, contains the sequence 
!    to be transformed, and on output, the transformed data.
!
!    Input, integer ( kind = 4 ) LENR, the dimension of the R array.  
!    LENR must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to RFFT1I before the first call to routine RFFT1F
!    or RFFT1B for a given transform length N. 
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK). 
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array. 
!    LENWRK must be at least N. 
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENR not big enough:
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough.
!
  implicit none

  integer ( kind = 4 ) lenr
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)
  real ( kind = 8 ) r(lenr)

  ier = 0

  if ( lenr < inc * ( n - 1 ) + 1 ) then
    ier = 1
    call xerfft ( 'rfft1f ', 6 )
    return
  end if

  if ( lensav < n + int ( log ( real ( n, kind = 4 ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfft1f ', 8 )
    return
  end if

  if ( lenwrk < n ) then
    ier = 3
    call xerfft ( 'rfft1f ', 10 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call rfftf1 ( n, inc, r, work, wsave, wsave(n+1) )

  return
end subroutine rfft1f
subroutine rfft1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! RFFT1I: initialization for RFFT1B and RFFT1F.
!
!  Discussion:
!
!    RFFT1I initializes array WSAVE for use in its companion routines 
!    RFFT1B and RFFT1F.  The prime factorization of N together with a
!    tabulation of the trigonometric functions are computed and stored
!    in array WSAVE.  Separate WSAVE arrays are required for different
!    values of N. 
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    25 March 2005
!
!  Author:
!
!    Paul Swarztrauber
!    Richard Valent
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors of
!    N and also containing certain trigonometric values which will be used in
!    routines RFFT1B or RFFT1F.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.
!
  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if ( lensav < n + int ( log ( real ( n, kind = 4 ) ) ) + 4 ) then
    ier = 2
    call xerfft ( 'rfft1i ', 3 )
    return
  end if

  if ( n == 1 ) then
    return
  end if

  call rffti1 ( n, wsave(1), wsave(n+1) )

  return
end subroutine rfft1i
subroutine xerfft ( srname, info )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERFFT - Fatal error!'

  if ( 1 <= info ) then
    write ( *, '(a,a,a,i3,a)') '  On entry to ', trim ( srname ), &
      ' parameter number ', info, ' had an illegal value.'
  else if ( info == -1 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameters LOT, JUMP, N and INC are inconsistent.'
  else if ( info == -2 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter L is greater than LDIM.'
  else if ( info == -3 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter M is greater than MDIM.'
  else if ( info == -5 ) then
    write( *, '(a,a,a,a)') '  Within ', trim ( srname ), &
      ' input error returned by lower level routine.'
  else if ( info == -6 ) then
    write( *, '(a,a,a,a)') '  On entry to ', trim ( srname ), &
      ' parameter LDIM is less than 2*(L/2+1).'
  end if

  stop
end subroutine xerfft


subroutine rfftb1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(in,*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(15)
  real ( kind = 8 ) half
  real ( kind = 8 ) halfm
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) wa(n)

  nf = int ( fac(2) )
  na = 0

  do k1 = 1, nf

    ip = int ( fac(k1+2) )
    na = 1 - na

    if ( 5 < ip ) then
      if ( k1 /= nf ) then
        na = 1 - na
      end if
    end if

  end do

  half = 0.5E+00
  halfm = -0.5E+00
  modn = mod ( n, 2 )
  nl = n - 2
  if ( modn /= 0 ) then
    nl = n - 1
  end if

  if ( na == 0 ) then

    do j = 2, nl, 2
      c(1,j) = half * c(1,j)
      c(1,j+1) = halfm * c(1,j+1)
    end do

  else

    ch(1) = c(1,1)
    ch(n) = c(1,n)

    do j = 2, nl, 2
      ch(j) = half*c(1,j)
      ch(j+1) = halfm*c(1,j+1)
    end do

  end if

  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = int ( fac(k1+2) )
    l2 = ip * l1
    ido = n / l2
    idl1 = ido * l1

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call r1f4kb ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3) )
      else
        call r1f4kb ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call r1f2kb ( ido, l1, c, in, ch, 1, wa(iw) )
      else
        call r1f2kb ( ido, l1, ch, 1, c, in, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call r1f3kb ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2) )
      else
        call r1f3kb ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call r1f5kb ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call r1f5kb ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call r1fgkb ( ido, ip, l1, idl1, c, c, c, in, ch, ch, 1, wa(iw) )
      else
        call r1fgkb ( ido, ip, l1, idl1, ch, ch, ch, 1, c, c, in, wa(iw) )
      end if

      if ( ido == 1 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ido

  end do

  return
end subroutine rfftb1
!
!---------------------------------------------
subroutine r1fgkb ( ido, ip, l1, idl1, cc, c1, c2, in1, ch, ch2, in2, wa )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(in1,ido,l1,ip)
  real ( kind = 8 ) c2(in1,idl1,ip)
  real ( kind = 8 ) cc(in1,ido,ip,l1)
  real ( kind = 8 ) ch(in2,ido,l1,ip)
  real ( kind = 8 ) ch2(in2,idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido)

  tpi = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
  arg = tpi / real ( ip, kind = 4 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  idp2 = ido + 2
  nbd = ( ido - 1 ) / 2
  ipp2 = ip + 2
  ipph = ( ip + 1 ) / 2

  if ( ido < l1 ) then
    do i = 1, ido
      do k = 1, l1
        ch(1,i,k,1) = cc(1,i,1,k)
      end do
    end do
  else
    do k = 1, l1
      do i = 1, ido
        ch(1,i,k,1) = cc(1,i,1,k)
      end do
    end do
  end if

  do j = 2, ipph
    jc = ipp2 - j
    j2 = j + j
    do k = 1, l1
      ch(1,1,k,j) = cc(1,ido,j2-2,k)+cc(1,ido,j2-2,k)
      ch(1,1,k,jc) = cc(1,1,j2-1,k)+cc(1,1,j2-1,k)
    end do
  end do

  if ( ido == 1 ) then

  else if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      do i = 3, ido, 2
        ic = idp2 - i
        do k = 1, l1
          ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
          ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
          ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
          ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 3, ido, 2
          ic = idp2 - i
          ch(1,i-1,k,j) = cc(1,i-1,2*j-1,k)+cc(1,ic-1,2*j-2,k)
          ch(1,i-1,k,jc) = cc(1,i-1,2*j-1,k)-cc(1,ic-1,2*j-2,k)
          ch(1,i,k,j) = cc(1,i,2*j-1,k)-cc(1,ic,2*j-2,k)
          ch(1,i,k,jc) = cc(1,i,2*j-1,k)+cc(1,ic,2*j-2,k)
        end do
      end do
    end do

  end if

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ipp2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 = dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      c2(1,ik,l) = ch2(1,ik,1)+ar1*ch2(1,ik,2)
      c2(1,ik,lc) = ai1*ch2(1,ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ipp2 - j
      ar2h = dc2*ar2-ds2*ai2
      ai2 = dc2*ai2+ds2*ar2
      ar2 = ar2h

      do ik = 1, idl1
        c2(1,ik,l) = c2(1,ik,l)+ar2*ch2(1,ik,j)
        c2(1,ik,lc) = c2(1,ik,lc)+ai2*ch2(1,ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    do ik = 1, idl1
      ch2(1,ik,1) = ch2(1,ik,1)+ch2(1,ik,j)
    end do
  end do

  do j = 2, ipph
    jc = ipp2 - j
    do k = 1, l1
      ch(1,1,k,j) = c1(1,1,k,j)-c1(1,1,k,jc)
      ch(1,1,k,jc) = c1(1,1,k,j)+c1(1,1,k,jc)
    end do
  end do

  if ( ido == 1 ) then

  else if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      do i = 3, ido, 2
        do k = 1, l1
          ch(1,i-1,k,j)  = c1(1,i-1,k,j) - c1(1,i,k,jc)
          ch(1,i-1,k,jc) = c1(1,i-1,k,j) + c1(1,i,k,jc)
          ch(1,i,k,j)    = c1(1,i,k,j)   + c1(1,i-1,k,jc)
          ch(1,i,k,jc)   = c1(1,i,k,j)   - c1(1,i-1,k,jc)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      do k = 1, l1
        do i = 3, ido, 2
          ch(1,i-1,k,j) = c1(1,i-1,k,j)-c1(1,i,k,jc)
          ch(1,i-1,k,jc) = c1(1,i-1,k,j)+c1(1,i,k,jc)
          ch(1,i,k,j) = c1(1,i,k,j)+c1(1,i-1,k,jc)
          ch(1,i,k,jc) = c1(1,i,k,j)-c1(1,i-1,k,jc)
        end do
      end do
    end do

  end if

  if ( ido == 1 ) then
    return
  end if

  do ik = 1, idl1
    c2(1,ik,1) = ch2(1,ik,1)
  end do

  do j = 2, ip
    do k = 1, l1
      c1(1,1,k,j) = ch(1,1,k,j)
    end do
  end do

  if ( l1 < nbd ) then

    is = -ido
    do j = 2, ip
       is = is + ido
       do k = 1, l1
         idij = is
         do i = 3, ido, 2
           idij = idij + 2
           c1(1,i-1,k,j) = wa(idij-1)*ch(1,i-1,k,j)-wa(idij)* ch(1,i,k,j)
           c1(1,i,k,j) = wa(idij-1)*ch(1,i,k,j)+wa(idij)* ch(1,i-1,k,j)
         end do
       end do
    end do

  else

    is = -ido

    do j = 2, ip
      is = is + ido
      idij = is
      do i = 3, ido, 2
        idij = idij + 2
        do k = 1, l1
           c1(1,i-1,k,j) = wa(idij-1) * ch(1,i-1,k,j) - wa(idij) * ch(1,i,k,j)
           c1(1,i,k,j)   = wa(idij-1) * ch(1,i,k,j)   + wa(idij) * ch(1,i-1,k,j)
        end do
      end do
    end do

  end if

  return
end subroutine r1fgkb
!
!---------------------------------------------
subroutine r1fgkf ( ido, ip, l1, idl1, cc, c1, c2, in1, ch, ch2, in2, wa )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ai1
  real ( kind = 8 ) ai2
  real ( kind = 8 ) ar1
  real ( kind = 8 ) ar1h
  real ( kind = 8 ) ar2
  real ( kind = 8 ) ar2h
  real ( kind = 8 ) arg
  real ( kind = 8 ) c1(in1,ido,l1,ip)
  real ( kind = 8 ) c2(in1,idl1,ip)
  real ( kind = 8 ) cc(in1,ido,ip,l1)
  real ( kind = 8 ) ch(in2,ido,l1,ip)
  real ( kind = 8 ) ch2(in2,idl1,ip)
  real ( kind = 8 ) dc2
  real ( kind = 8 ) dcp
  real ( kind = 8 ) ds2
  real ( kind = 8 ) dsp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idij
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) ik
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) nbd
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido)

  tpi = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 )
  arg = tpi / real ( ip, kind = 4 )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  ipp2 = ip + 2
  idp2 = ido + 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then

    do ik = 1, idl1
      c2(1,ik,1) = ch2(1,ik,1)
    end do
 
  else

    do ik = 1, idl1
      ch2(1,ik,1) = c2(1,ik,1)
    end do

    do j = 2, ip
      do k = 1, l1
        ch(1,1,k,j) = c1(1,1,k,j)
      end do
    end do

    if ( l1 < nbd ) then

      is = -ido

      do j = 2, ip
        is = is + ido
        do k = 1, l1
          idij = is
          do i = 3, ido, 2
            idij = idij + 2
            ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij) *c1(1,i,k,j)
            ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij) *c1(1,i-1,k,j)
          end do
        end do
      end do

    else

      is = -ido

      do j = 2, ip
        is = is + ido
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          do k = 1, l1
            ch(1,i-1,k,j) = wa(idij-1)*c1(1,i-1,k,j)+wa(idij) *c1(1,i,k,j)
            ch(1,i,k,j) = wa(idij-1)*c1(1,i,k,j)-wa(idij) *c1(1,i-1,k,j)
          end do
        end do
      end do

    end if

    if ( nbd < l1 ) then

      do j = 2, ipph
        jc = ipp2 - j
        do i = 3, ido, 2
          do k = 1, l1
            c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
            c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
            c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
            c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ipp2 - j
        do k = 1, l1
          do i = 3, ido, 2
            c1(1,i-1,k,j) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
            c1(1,i-1,k,jc) = ch(1,i,k,j)-ch(1,i,k,jc)
            c1(1,i,k,j) = ch(1,i,k,j)+ch(1,i,k,jc)
            c1(1,i,k,jc) = ch(1,i-1,k,jc)-ch(1,i-1,k,j)
          end do
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ipp2 - j
    do k = 1, l1
      c1(1,1,k,j) = ch(1,1,k,j)+ch(1,1,k,jc)
      c1(1,1,k,jc) = ch(1,1,k,jc)-ch(1,1,k,j)
    end do
  end do

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ipp2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 = dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      ch2(1,ik,l) = c2(1,ik,1)+ar1*c2(1,ik,2)
      ch2(1,ik,lc) = ai1*c2(1,ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph
      jc = ipp2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2 = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h
      do ik = 1, idl1
        ch2(1,ik,l) = ch2(1,ik,l)+ar2*c2(1,ik,j)
        ch2(1,ik,lc) = ch2(1,ik,lc)+ai2*c2(1,ik,jc)
      end do
    end do

  end do

  do j = 2, ipph
    do ik = 1, idl1
      ch2(1,ik,1) = ch2(1,ik,1)+c2(1,ik,j)
    end do
  end do

  if ( ido < l1 ) then

    do i = 1, ido
      do k = 1, l1
        cc(1,i,1,k) = ch(1,i,k,1)
      end do
    end do

  else

    do k = 1, l1
      do i = 1, ido
        cc(1,i,1,k) = ch(1,i,k,1)
      end do
    end do

  end if

  do j = 2, ipph
    jc = ipp2 - j
    j2 = j+j
    do k = 1, l1
      cc(1,ido,j2-2,k) = ch(1,1,k,j)
      cc(1,1,j2-1,k) = ch(1,1,k,jc)
    end do
  end do

  if ( ido == 1 ) then
    return
  end if

  if ( nbd < l1 ) then

    do j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      do i = 3, ido, 2
        ic = idp2 - i
        do k = 1, l1
          cc(1,i-1,j2-1,k) = ch(1,i-1,k,j)+ch(1,i-1,k,jc)
          cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j)-ch(1,i-1,k,jc)
          cc(1,i,j2-1,k) = ch(1,i,k,j)+ch(1,i,k,jc)
          cc(1,ic,j2-2,k) = ch(1,i,k,jc)-ch(1,i,k,j)
        end do
      end do
   end do

  else

    do j = 2, ipph
      jc = ipp2 - j
      j2 = j + j
      do k = 1, l1
        do i = 3, ido, 2
          ic = idp2 - i
          cc(1,i-1,j2-1,k)  = ch(1,i-1,k,j) + ch(1,i-1,k,jc)
          cc(1,ic-1,j2-2,k) = ch(1,i-1,k,j) - ch(1,i-1,k,jc)
          cc(1,i,j2-1,k)    = ch(1,i,k,j)   + ch(1,i,k,jc)
          cc(1,ic,j2-2,k)   = ch(1,i,k,jc)  - ch(1,i,k,j)
        end do
      end do
    end do

  end if

  return
end subroutine r1fgkf
!
!-------------------------------------------------------
subroutine r1f2kb ( ido, l1, cc, in1, ch, in2, wa1 )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,2,l1)
  real ( kind = 8 ) ch(in2,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) wa1(ido)

  do k = 1, l1
    ch(1,1,k,1) = cc(1,1,1,k) + cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k) - cc(1,ido,2,k)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2
    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i

        ch(1,i-1,k,1) = cc(1,i-1,1,k) + cc(1,ic-1,2,k)
        ch(1,i,k,1)   = cc(1,i,1,k)   - cc(1,ic,2,k)

        ch(1,i-1,k,2) = wa1(i-2) * ( cc(1,i-1,1,k) - cc(1,ic-1,2,k) ) &
                      - wa1(i-1) * ( cc(1,i,1,k)   + cc(1,ic,2,k) )
        ch(1,i,k,2)   = wa1(i-2) * ( cc(1,i,1,k)   + cc(1,ic,2,k) ) &
                      + wa1(i-1) * ( cc(1,i-1,1,k) - cc(1,ic-1,2,k) )

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,ido,k,1) =     cc(1,ido,1,k) + cc(1,ido,1,k)
    ch(1,ido,k,2) = - ( cc(1,1,2,k)   + cc(1,1,2,k) )
  end do

  return
end subroutine r1f2kb
!
!-------------------------------------------------------
subroutine r1f2kf ( ido, l1, cc, in1, ch, in2, wa1 )

!*****************************************************************************80
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) ch(in2,ido,2,l1)
  real ( kind = 8 ) cc(in1,ido,l1,2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) wa1(ido)

  do k = 1, l1
    ch(1,1,1,k)   = cc(1,1,k,1) + cc(1,1,k,2)
    ch(1,ido,2,k) = cc(1,1,k,1) - cc(1,1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        ch(1,i,1,k) = cc(1,i,k,1)   + ( wa1(i-2) * cc(1,i,k,2) &
                                      - wa1(i-1) * cc(1,i-1,k,2) )
        ch(1,ic,2,k) = -cc(1,i,k,1) + ( wa1(i-2) * cc(1,i,k,2) &
                                      - wa1(i-1) * cc(1,i-1,k,2) ) 
        ch(1,i-1,1,k) = cc(1,i-1,k,1)  + ( wa1(i-2) * cc(1,i-1,k,2) &
                                         + wa1(i-1) * cc(1,i,k,2))
        ch(1,ic-1,2,k) = cc(1,i-1,k,1) - ( wa1(i-2) * cc(1,i-1,k,2) &
                                         + wa1(i-1) * cc(1,i,k,2))
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,1,2,k) = -cc(1,ido,k,2)
    ch(1,ido,1,k) = cc(1,ido,k,1)
  end do

  return
end subroutine r1f2kf
!
!-------------------------------------------------------
subroutine r1f3kb ( ido, l1, cc, in1, ch, in2, wa1, wa2 )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,3,l1)
  real ( kind = 8 ) ch(in2,ido,l1,3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) taui
  real ( kind = 8 ) taur
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    ch(1,1,k,1) = cc(1,1,1,k) + 2.0E+00 * cc(1,ido,2,k)
    ch(1,1,k,2) = cc(1,1,1,k) + 2.0E+00 * taur * cc(1,ido,2,k) &
                              - 2.0E+00 * taui * cc(1,1,3,k)
    ch(1,1,k,3) = cc(1,1,1,k) + 2.0E+00 * taur * cc(1,ido,2,k) &
                              + 2.0E+00 * taui * cc(1,1,3,k)
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
      ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k))
      ch(1,i-1,k,2) = wa1(i-2)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) -wa1(i-1)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))))
      ch(1,i,k,2) = wa1(i-2)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))+ &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) +wa1(i-1)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))- &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))
      ch(1,i-1,k,3) = wa2(i-2)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k)))) -wa2(i-1)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))))
      ch(1,i,k,3) = wa2(i-2)* &
        ((cc(1,i,1,k)+taur*(cc(1,i,3,k)-cc(1,ic,2,k)))- &
        (taui*(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))) +wa2(i-1)* &
        ((cc(1,i-1,1,k)+taur*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))+ &
        (taui*(cc(1,i,3,k)+cc(1,ic,2,k))))
    end do
  end do

  return
end subroutine r1f3kb
!
!-------------------------------------------------------
subroutine r1f3kf ( ido, l1, cc, in1, ch, in2, wa1, wa2 )

!*****************************************************************************80
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,l1,3)
  real ( kind = 8 ) ch(in2,ido,3,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) taui
  real ( kind = 8 ) taur
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 3.0E+00
  taur = cos ( arg )
  taui = sin ( arg )

  do k = 1, l1
    ch(1,1,1,k) = cc(1,1,k,1)          + ( cc(1,1,k,2) + cc(1,1,k,3) )
    ch(1,1,3,k) =                 taui * ( cc(1,1,k,3) - cc(1,1,k,2) )
    ch(1,ido,2,k) = cc(1,1,k,1) + taur * ( cc(1,1,k,2) + cc(1,1,k,3) )
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
        wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3)))
      ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3)))
      ch(1,i-1,3,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))+(taui*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))
      ch(1,ic-1,2,k) = (cc(1,i-1,k,1)+taur*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))))-(taui*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))))
      ch(1,i,3,k) = (cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))))+(taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2))))
      ch(1,ic,2,k) = (taui*((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2))))-(cc(1,i,k,1)+taur*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))))
    end do
  end do

  return
end subroutine r1f3kf
!
!-------------------------------------------------------
subroutine r1f4kb ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3 )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,4,l1)
  real ( kind = 8 ) ch(in2,ido,l1,4)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) sqrt2
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  sqrt2 = sqrt ( 2.0E+00 )

  do k = 1, l1
    ch(1,1,k,3) = ( cc(1,1,1,k)   + cc(1,ido,4,k) ) &
                - ( cc(1,ido,2,k) + cc(1,ido,2,k) )
    ch(1,1,k,1) = ( cc(1,1,1,k)   + cc(1,ido,4,k) ) &
                + ( cc(1,ido,2,k) + cc(1,ido,2,k) )
    ch(1,1,k,4) = ( cc(1,1,1,k)   - cc(1,ido,4,k) ) &
                + ( cc(1,1,3,k)   + cc(1,1,3,k) )
    ch(1,1,k,2) = ( cc(1,1,1,k)   - cc(1,ido,4,k) ) &
                - ( cc(1,1,3,k)   + cc(1,1,3,k) )
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        ch(1,i-1,k,1) = (cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
          +(cc(1,i-1,3,k)+cc(1,ic-1,2,k))
        ch(1,i,k,1) = (cc(1,i,1,k)-cc(1,ic,4,k)) &
          +(cc(1,i,3,k)-cc(1,ic,2,k))
        ch(1,i-1,k,2) = wa1(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
          -(cc(1,i,3,k)+cc(1,ic,2,k)))-wa1(i-1) &
          *((cc(1,i,1,k)+cc(1,ic,4,k))+(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,2) = wa1(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
          +(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa1(i-1) &
          *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))-(cc(1,i,3,k)+cc(1,ic,2,k)))
        ch(1,i-1,k,3) = wa2(i-2)*((cc(1,i-1,1,k)+cc(1,ic-1,4,k)) &
          -(cc(1,i-1,3,k)+cc(1,ic-1,2,k)))-wa2(i-1) &
          *((cc(1,i,1,k)-cc(1,ic,4,k))-(cc(1,i,3,k)-cc(1,ic,2,k)))
        ch(1,i,k,3) = wa2(i-2)*((cc(1,i,1,k)-cc(1,ic,4,k)) &
          -(cc(1,i,3,k)-cc(1,ic,2,k)))+wa2(i-1) &
          *((cc(1,i-1,1,k)+cc(1,ic-1,4,k))-(cc(1,i-1,3,k) &
          +cc(1,ic-1,2,k)))
        ch(1,i-1,k,4) = wa3(i-2)*((cc(1,i-1,1,k)-cc(1,ic-1,4,k)) &
          +(cc(1,i,3,k)+cc(1,ic,2,k)))-wa3(i-1) &
          *((cc(1,i,1,k)+cc(1,ic,4,k))-(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))
        ch(1,i,k,4) = wa3(i-2)*((cc(1,i,1,k)+cc(1,ic,4,k)) &
          -(cc(1,i-1,3,k)-cc(1,ic-1,2,k)))+wa3(i-1) &
          *((cc(1,i-1,1,k)-cc(1,ic-1,4,k))+(cc(1,i,3,k)+cc(1,ic,2,k)))
      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,ido,k,1) = ( cc(1,ido,1,k) + cc(1,ido,3,k) ) &
                  + ( cc(1,ido,1,k) + cc(1,ido,3,k))
    ch(1,ido,k,2) = sqrt2 * ( ( cc(1,ido,1,k) - cc(1,ido,3,k) ) &
                            - ( cc(1,1,2,k)   + cc(1,1,4,k) ) )
    ch(1,ido,k,3) = ( cc(1,1,4,k) - cc(1,1,2,k) ) &
                  + ( cc(1,1,4,k) - cc(1,1,2,k) )
    ch(1,ido,k,4) = -sqrt2 * ( ( cc(1,ido,1,k) - cc(1,ido,3,k) ) &
                             + ( cc(1,1,2,k) + cc(1,1,4,k) ) )
  end do

  return
end subroutine r1f4kb
!
!-------------------------------------------------------
subroutine r1f4kf ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3 )

!*****************************************************************************80
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,ido,l1,4)
  real ( kind = 8 ) ch(in2,ido,4,l1)
  real ( kind = 8 ) hsqt2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)

  hsqt2 = sqrt ( 2.0E+00 ) / 2.0E+00 

  do k = 1, l1
    ch(1,1,1,k)   = ( cc(1,1,k,2) + cc(1,1,k,4) ) &
                  + ( cc(1,1,k,1) + cc(1,1,k,3) )
    ch(1,ido,4,k) = ( cc(1,1,k,1) + cc(1,1,k,3) ) &
                  - ( cc(1,1,k,2) + cc(1,1,k,4) )
    ch(1,ido,2,k) = cc(1,1,k,1) - cc(1,1,k,3)
    ch(1,1,3,k)   = cc(1,1,k,4) - cc(1,1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( 2 < ido ) then

    idp2 = ido + 2

    do k = 1, l1
      do i = 3, ido, 2
        ic = idp2 - i
        ch(1,i-1,1,k) = ((wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
          cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
          cc(1,i,k,4)))+(cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))
        ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i-1,k,2)+ &
          wa1(i-1)*cc(1,i,k,2))+(wa3(i-2)*cc(1,i-1,k,4)+ &
          wa3(i-1)*cc(1,i,k,4)))
        ch(1,i,1,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))+(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
        ch(1,ic,4,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))-(cc(1,i,k,1)+(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
        ch(1,i-1,3,k) = ((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))+(cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))
        ch(1,ic-1,2,k) = (cc(1,i-1,k,1)-(wa2(i-2)*cc(1,i-1,k,3)+ &
          wa2(i-1)*cc(1,i,k,3)))-((wa1(i-2)*cc(1,i,k,2)-wa1(i-1)* &
          cc(1,i-1,k,2))-(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
          cc(1,i-1,k,4)))
        ch(1,i,3,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
          cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
          cc(1,i,k,2)))+(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
        ch(1,ic,2,k) = ((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
          cc(1,i,k,4))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
          cc(1,i,k,2)))-(cc(1,i,k,1)-(wa2(i-2)*cc(1,i,k,3)- &
          wa2(i-1)*cc(1,i-1,k,3)))
       end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1
    ch(1,ido,1,k) = (hsqt2*(cc(1,ido,k,2)-cc(1,ido,k,4)))+ cc(1,ido,k,1)
    ch(1,ido,3,k) = cc(1,ido,k,1)-(hsqt2*(cc(1,ido,k,2)- cc(1,ido,k,4)))
    ch(1,1,2,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))- cc(1,ido,k,3)
    ch(1,1,4,k) = (-hsqt2*(cc(1,ido,k,2)+cc(1,ido,k,4)))+ cc(1,ido,k,3)
  end do

  return
end subroutine r1f4kf
!
!-------------------------------------------------------
subroutine r1f5kb ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,5,l1)
  real ( kind = 8 ) ch(in2,ido,l1,5)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti11
  real ( kind = 8 ) ti12
  real ( kind = 8 ) tr11
  real ( kind = 8 ) tr12
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.0E+00
  tr11 = cos ( arg )
  ti11 = sin ( arg )
  tr12 = cos ( 2.0E+00 * arg )
  ti12 = sin ( 2.0E+00 * arg )

  do k = 1, l1

    ch(1,1,k,1) = cc(1,1,1,k) + 2.0E+00 * cc(1,ido,2,k) &
                              + 2.0E+00 * cc(1,ido,4,k)

    ch(1,1,k,2) = ( cc(1,1,1,k) &
      +   tr11 * 2.0E+00 * cc(1,ido,2,k) + tr12 * 2.0E+00 * cc(1,ido,4,k) ) &
      - ( ti11 * 2.0E+00 * cc(1,1,3,k)   + ti12 * 2.0E+00 * cc(1,1,5,k))

    ch(1,1,k,3) = ( cc(1,1,1,k) &
      +   tr12 * 2.0E+00 * cc(1,ido,2,k) + tr11 * 2.0E+00 * cc(1,ido,4,k) ) &
      - ( ti12 * 2.0E+00 * cc(1,1,3,k)   - ti11 * 2.0E+00 * cc(1,1,5,k))

    ch(1,1,k,4) = ( cc(1,1,1,k) &
      +   tr12 * 2.0E+00 * cc(1,ido,2,k) + tr11 * 2.0E+00 * cc(1,ido,4,k) ) &
      + ( ti12 * 2.0E+00 * cc(1,1,3,k)   - ti11 * 2.0E+00 * cc(1,1,5,k))

    ch(1,1,k,5) = ( cc(1,1,1,k) &
      +   tr11 * 2.0E+00 * cc(1,ido,2,k) + tr12 * 2.0E+00 * cc(1,ido,4,k) ) &
      + ( ti11 * 2.0E+00 * cc(1,1,3,k)   + ti12 * 2.0E+00 * cc(1,1,5,k) )

  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,k,1) = cc(1,i-1,1,k)+(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +(cc(1,i-1,5,k)+cc(1,ic-1,4,k))
      ch(1,i,k,1) = cc(1,i,1,k)+(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +(cc(1,i,5,k)-cc(1,ic,4,k))
      ch(1,i-1,k,2) = wa1(i-2)*((cc(1,i-1,1,k)+tr11* &
        (cc(1,i-1,3,k)+cc(1,ic-1,2,k))+tr12 &
        *(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa1(i-1)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))+(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,2) = wa1(i-2)*((cc(1,i,1,k)+tr11*(cc(1,i,3,k) &
        -cc(1,ic,2,k))+tr12*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti11*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))+ti12 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))+wa1(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k) &
        +cc(1,ic-1,2,k))+tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k))) &
        -(ti11*(cc(1,i,3,k)+cc(1,ic,2,k))+ti12 &
        *(cc(1,i,5,k)+cc(1,ic,4,k))))
      ch(1,i-1,k,3) = wa2(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa2(i-1) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
      cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,3) = wa2(i-2) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        +(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa2(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))-(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))
      ch(1,i-1,k,4) = wa3(i-2) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa3(i-1) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
      cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,4) = wa3(i-2) &
        *((cc(1,i,1,k)+tr12*(cc(1,i,3,k)- &
        cc(1,ic,2,k))+tr11*(cc(1,i,5,k)-cc(1,ic,4,k))) &
        -(ti12*(cc(1,i-1,3,k)-cc(1,ic-1,2,k))-ti11 &
        *(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa3(i-1) &
        *((cc(1,i-1,1,k)+tr12*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr11*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti12*(cc(1,i,3,k) &
        +cc(1,ic,2,k))-ti11*(cc(1,i,5,k)+cc(1,ic,4,k))))
      ch(1,i-1,k,5) = wa4(i-2) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k)))) &
        -wa4(i-1) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k))))
      ch(1,i,k,5) = wa4(i-2) &
        *((cc(1,i,1,k)+tr11*(cc(1,i,3,k)-cc(1,ic,2,k)) &
        +tr12*(cc(1,i,5,k)-cc(1,ic,4,k)))-(ti11*(cc(1,i-1,3,k) &
        -cc(1,ic-1,2,k))+ti12*(cc(1,i-1,5,k)-cc(1,ic-1,4,k)))) &
        +wa4(i-1) &
        *((cc(1,i-1,1,k)+tr11*(cc(1,i-1,3,k)+cc(1,ic-1,2,k)) &
        +tr12*(cc(1,i-1,5,k)+cc(1,ic-1,4,k)))+(ti11*(cc(1,i,3,k) &
        +cc(1,ic,2,k))+ti12*(cc(1,i,5,k)+cc(1,ic,4,k))))
    end do
  end do

  return
end subroutine r1f5kb
!
!-------------------------------------------------------
subroutine r1f5kf ( ido, l1, cc, in1, ch, in2, wa1, wa2, wa3, wa4 )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) arg
  real ( kind = 8 ) cc(in1,ido,l1,5)
  real ( kind = 8 ) ch(in2,ido,5,l1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ic
  integer ( kind = 4 ) idp2
  integer ( kind = 4 ) k
  real ( kind = 8 ) ti11
  real ( kind = 8 ) ti12
  real ( kind = 8 ) tr11
  real ( kind = 8 ) tr12
  real ( kind = 8 ) wa1(ido)
  real ( kind = 8 ) wa2(ido)
  real ( kind = 8 ) wa3(ido)
  real ( kind = 8 ) wa4(ido)

  arg = 2.0E+00 * 4.0E+00 * atan ( 1.0E+00 ) / 5.0E+00
  tr11 = cos ( arg )
  ti11 = sin ( arg )
  tr12 = cos ( 2.0E+00 * arg )
  ti12 = sin ( 2.0E+00 * arg )

  do k = 1, l1

    ch(1,1,1,k) = cc(1,1,k,1) + ( cc(1,1,k,5) + cc(1,1,k,2) ) &
                              + ( cc(1,1,k,4) + cc(1,1,k,3) )

    ch(1,ido,2,k) = cc(1,1,k,1) + tr11 * ( cc(1,1,k,5) + cc(1,1,k,2) ) &
                                + tr12 * ( cc(1,1,k,4) + cc(1,1,k,3) )

    ch(1,1,3,k) =                 ti11 * ( cc(1,1,k,5) - cc(1,1,k,2) ) &
                                + ti12 * ( cc(1,1,k,4) - cc(1,1,k,3) )

    ch(1,ido,4,k) = cc(1,1,k,1) + tr12 * ( cc(1,1,k,5) + cc(1,1,k,2) ) &
                                + tr11 * ( cc(1,1,k,4) + cc(1,1,k,3) )

    ch(1,1,5,k) =                 ti12 * ( cc(1,1,k,5) - cc(1,1,k,2) ) &
                                - ti11 * ( cc(1,1,k,4) - cc(1,1,k,3) )
  end do

  if ( ido == 1 ) then
    return
  end if

  idp2 = ido + 2

  do k = 1, l1
    do i = 3, ido, 2
      ic = idp2 - i
      ch(1,i-1,1,k) = cc(1,i-1,k,1)+((wa1(i-2)*cc(1,i-1,k,2)+ &
        wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
        cc(1,i,k,5)))+((wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))+(wa3(i-2)*cc(1,i-1,k,4)+ &
        wa3(i-1)*cc(1,i,k,4)))
      ch(1,i,1,k) = cc(1,i,k,1)+((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4)))
      ch(1,i-1,3,k) = cc(1,i-1,k,1)+tr11* &
        ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
        +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
        ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
        +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))+ti11* &
        ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
        -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
        ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
        -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4)))
      ch(1,ic-1,2,k) = cc(1,i-1,k,1)+tr11* &
        ( wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2) &
        +wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5))+tr12* &
        ( wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3) &
        +wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))-(ti11* &
        ( wa1(i-2)*cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2) &
        -(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))+ti12* &
        ( wa2(i-2)*cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3) &
        -(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))
      ch(1,i,3,k) = (cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))+(ti11*((wa4(i-2)*cc(1,i-1,k,5)+ &
        wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))
      ch(1,ic,2,k) = (ti11*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
        cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))+ti12*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))-(cc(1,i,k,1)+tr11*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr12*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))
      ch(1,i-1,5,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
        cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
        cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))+(ti12*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
        cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
        cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))
      ch(1,ic-1,4,k) = (cc(1,i-1,k,1)+tr12*((wa1(i-2)* &
        cc(1,i-1,k,2)+wa1(i-1)*cc(1,i,k,2))+(wa4(i-2)* &
        cc(1,i-1,k,5)+wa4(i-1)*cc(1,i,k,5)))+tr11*((wa2(i-2)* &
        cc(1,i-1,k,3)+wa2(i-1)*cc(1,i,k,3))+(wa3(i-2)* &
        cc(1,i-1,k,4)+wa3(i-1)*cc(1,i,k,4))))-(ti12*((wa1(i-2)* &
        cc(1,i,k,2)-wa1(i-1)*cc(1,i-1,k,2))-(wa4(i-2)* &
        cc(1,i,k,5)-wa4(i-1)*cc(1,i-1,k,5)))-ti11*((wa2(i-2)* &
        cc(1,i,k,3)-wa2(i-1)*cc(1,i-1,k,3))-(wa3(i-2)* &
        cc(1,i,k,4)-wa3(i-1)*cc(1,i-1,k,4))))
      ch(1,i,5,k) = (cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))+(ti12*((wa4(i-2)*cc(1,i-1,k,5)+ &
        wa4(i-1)*cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))
      ch(1,ic,4,k) = (ti12*((wa4(i-2)*cc(1,i-1,k,5)+wa4(i-1)* &
        cc(1,i,k,5))-(wa1(i-2)*cc(1,i-1,k,2)+wa1(i-1)* &
        cc(1,i,k,2)))-ti11*((wa3(i-2)*cc(1,i-1,k,4)+wa3(i-1)* &
        cc(1,i,k,4))-(wa2(i-2)*cc(1,i-1,k,3)+wa2(i-1)* &
        cc(1,i,k,3))))-(cc(1,i,k,1)+tr12*((wa1(i-2)*cc(1,i,k,2)- &
        wa1(i-1)*cc(1,i-1,k,2))+(wa4(i-2)*cc(1,i,k,5)-wa4(i-1)* &
        cc(1,i-1,k,5)))+tr11*((wa2(i-2)*cc(1,i,k,3)-wa2(i-1)* &
        cc(1,i-1,k,3))+(wa3(i-2)*cc(1,i,k,4)-wa3(i-1)* &
        cc(1,i-1,k,4))))
     end do
  end do

  return
end subroutine r1f5kf

!
!-------------------------------------------------------
subroutine r4_factor ( n, nf, fac )

!*****************************************************************************80
!
  implicit none

  real ( kind = 8 ) fac(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf) = real ( ntry, kind = 4 )
      nl = nq

    end do

  end do

  return
end subroutine r4_factor
!
!-----------------------------------------------------
subroutine rfftf1 ( n, in, c, ch, wa, fac )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) in
  integer ( kind = 4 ) n

  real ( kind = 8 ) c(in,*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(15)
  integer ( kind = 4 ) idl1
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) ix2
  integer ( kind = 4 ) ix3
  integer ( kind = 4 ) ix4
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) kh
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) modn
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  real ( kind = 8 ) sn
  real ( kind = 8 ) tsn
  real ( kind = 8 ) tsnm
  real ( kind = 8 ) wa(n)

  nf = int ( fac(2) )
  na = 1
  l2 = n
  iw = n

  do k1 = 1, nf

    kh = nf - k1
    ip = int ( fac(kh+3) )
    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call r1f4kf ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3) )
      else
        call r1f4kf ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call r1f2kf ( ido, l1, c, in, ch, 1, wa(iw) )
      else
        call r1f2kf ( ido, l1, ch, 1, c, in, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call r1f3kf ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2) )
      else
        call r1f3kf ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call r1f5kf ( ido, l1, c, in, ch, 1, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call r1f5kf ( ido, l1, ch, 1, c, in, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call r1fgkf ( ido, ip, l1, idl1, c, c, c, in, ch, ch, 1, wa(iw) )
        na = 1
      else
        call r1fgkf ( ido, ip, l1, idl1, ch, ch, ch, 1, c, c, in, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  sn = 1.0E+00 / real ( n, kind = 4 )
  tsn = 2.0E+00 / real ( n, kind = 4 )
  tsnm = -tsn
  modn = mod ( n, 2 )
  nl = n - 2

  if ( modn /= 0 ) then
    nl = n - 1
  end if

  if ( na == 0 ) then

    c(1,1) = sn * ch(1)
    do j = 2, nl, 2
      c(1,j) = tsn * ch(j)
      c(1,j+1) = tsnm * ch(j+1)
    end do

    if ( modn == 0 ) then
      c(1,n) = sn * ch(n)
    end if

  else

    c(1,1) = sn * c(1,1)

    do j = 2, nl, 2
      c(1,j) = tsn * c(1,j)
      c(1,j+1) = tsnm * c(1,j+1)
    end do

    if ( modn == 0 ) then
      c(1,n) = sn * c(1,n)
    end if

  end if

  return
end subroutine rfftf1
!
!------------------------------------------------------
subroutine rffti1 ( n, wa, fac )

!*****************************************************************************80
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) arg
  real ( kind = 8 ) argh
  real ( kind = 8 ) argld
  real ( kind = 8 ) fac(15)
  real ( kind = 8 ) fi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipm
  integer ( kind = 4 ) is
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) ld
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nfm1
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(n)

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf+2) = real ( ntry, kind = 4 )
      nl = nq
!
!  If 2 is a factor, make sure it appears first in the list of factors.
!
      if ( ntry == 2 ) then
        if ( nf /= 1 ) then
          do i = 2, nf
            ib = nf - i + 2
            fac(ib+2) = fac(ib+1)
          end do
          fac(3) = 2.0E+00
        end if
      end if

    end do

  end do

  fac(1) = real ( n, kind = 4 )
  fac(2) = real ( nf, kind = 4 )
  tpi = 8.0D+00 * atan ( 1.0D+00 )
  argh = tpi / real ( n, kind = 4 )
  is = 0
  nfm1 = nf-1
  l1 = 1

  do k1 = 1, nfm1
    ip = int ( fac(k1+2) )
    ld = 0
    l2 = l1 * ip
    ido = n / l2
    ipm = ip - 1
    do j = 1, ipm
      ld = ld + l1
      i = is
      argld = real ( ld, kind = 4 ) * argh
      fi = 0.0E+00
      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0E+00
        arg = fi * argld
        wa(i-1) = real ( cos ( arg ), kind = 8 )
        wa(i) = real ( sin ( arg ), kind = 8 )
      end do
      is = is + ido
    end do
    l1 = l2
  end do

  return
end subroutine rffti1

  
end MODULE FFT5

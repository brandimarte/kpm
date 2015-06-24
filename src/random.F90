!  *******************************************************************  !
!  KPM Fortran Code 2014                                                !
!                                                                       !
!  Written by Eric de Castro e Andrade (eandrade@ift.unesp.br) and      !
!             Pedro Brandimarte (brandimarte@gmail.com).                !
!                                                                       !
!  Copyright (c), All Rights Reserved                                   !
!                                                                       !
!  This program is free software. You can redistribute it and/or        !
!  modify it under the terms of the GNU General Public License          !
!  (version 3 or later) as published by the Free Software Foundation    !
!  <http://fsf.org/>.                                                   !
!                                                                       !
!  This program is distributed in the hope that it will be useful, but  !
!  WITHOUT ANY WARRANTY, without even the implied warranty of           !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     !
!  General Public License for more details (file 'LICENSE_GPL'          !
!  distributed along with this program or at                            !
!  <http://www.gnu.org/licenses/gpl.html>).                             !
!  *******************************************************************  !
!                             MODULE random                             !
!  *******************************************************************  !
!  Description: subroutines for generating and testing pseudorandom     !
!  numbers.                                                             !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE random

!
! Modules
!
  use precision,       only: dp

  implicit none

  PUBLIC  :: sgenrand, genrand
  PRIVATE ! default is private


CONTAINS


!  *******************************************************************  !
!                               sgenrand                                !
!  *******************************************************************  !
!  Description: Set initial values to the working area of 624 words.    !
!  ('seed' is any 32-bit integer except for 0).                         !
!                                                                       !
!  Coded by Takuji Nishimura, considering the suggestions by Topher     !
!  Cooper and Marc Rieffel in July-Aug. 1997.                           !
!                                                                       !
!  Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.            !
!  When you use this, send an email to: matumoto@math.keio.ac.jp        !
!  with an appropriate reference to your work.                          !
!  *******************************************************************  !
  subroutine sgenrand (seed)

    implicit integer(a-z)

!   Period parameters.
    parameter (N = 624)
    dimension mt(0:N-1) ! the array for the state vector
    common /block/ mti, mt
    save   /block/

!   Setting initial seeds to mt[N] using the generator
!   at Line 25 of Table 1 in [KNUTH 1981, The Art of
!   Computer Programming Vol. 2 (2nd Ed.), pp102].
    mt(0) = iand (seed, -1)
    do 100 mti = 1,N-1
       mt(mti) = iand (69069 * mt(mti-1), -1)
100 continue

  end subroutine sgenrand


!  *******************************************************************  !
!                                genrand                                !
!  *******************************************************************  !
!  Description: Generates one pseudorandom real number (double) which   !
!  is uniformly distributed on [0,1]-interval, for each call. Before    !
!  'genrand ()', 'sgenrand (seed)' must be called once.                 !
!                                                                       !
!  Integer generator is obtained by modifying two lines.                !
!                                                                       !
!  Coded by Takuji Nishimura, considering the suggestions by Topher     !
!  Cooper and Marc Rieffel in July-Aug. 1997.                           !
!                                                                       !
!  Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.            !
!  When you use this, send an email to: matumoto@math.keio.ac.jp        !
!  with an appropriate reference to your work.                          !
!  *******************************************************************  !
  real(dp) function genrand ()

    implicit integer(a-z)

!   Period parameters.
    parameter (N      =  624)
    parameter (N1     =  N+1)
    parameter (M      =  397)
    parameter (MATA   = -1727483681) ! constant vector a
    parameter (UMASK  = -2147483648) ! most significant w-r bits
    parameter (LMASK  =  2147483647) ! least significant r bits

!   Tempering parameters.
    parameter (TMASKB = -1658038656)
    parameter (TMASKC = -272236544)
    dimension mt(0:N-1) ! the array for the state vector
    common /block/ mti, mt
    save   /block/
    data   mti/N1/ ! 'mti==N+1' means 'mt[N]' is not initialized

    dimension mag01(0:1)
    data mag01/0, MATA/
    save mag01 ! mag01(x) = x * MATA for x=0,1

    TSHFTU(y) = ishft(y,-11)
    TSHFTS(y) = ishft(y,7)
    TSHFTT(y) = ishft(y,15)
    TSHFTL(y) = ishft(y,-18)

    if (mti >= N) then ! generate N words at one time

!      Check if 'sgenrand ()' has not been called. 
       if (mti == N+1) then
          call sgenrand (4357) ! a default initial seed is used
       endif

       do 100 kk = 0,N-M-1

          y = ior (iand(mt(kk),UMASK), iand(mt(kk+1),LMASK))
          mt(kk) = ieor (ieor(mt(kk+M),ishft(y,-1)), mag01(iand(y,1)))

100    continue

       do 110 kk = N-M,N-2
          y = ior (iand(mt(kk),UMASK), iand(mt(kk+1),LMASK))
          mt(kk) = ieor (ieor(mt(kk+(M-N)),ishft(y,-1)),                &
               mag01(iand(y,1)))

110    continue

       y = ior (iand(mt(N-1),UMASK), iand(mt(0),LMASK))
       mt(N-1) = ieor (ieor(mt(M-1),ishft(y,-1)), mag01(iand(y,1)))
       mti = 0

    endif

    y = mt(mti)
    mti = mti + 1
    y = ieor (y, TSHFTU(y))
    y = ieor (y, iand(TSHFTS(y),TMASKB))
    y = ieor (y, iand(TSHFTT(y),TMASKC))
    y = ieor (y, TSHFTL(y))

    if (y < 0) then
       genrand = (dble(y) + 2.0d0**32) / (2.0d0**32 - 1.0d0)
    else
       genrand = dble(y) / (2.0d0**32 - 1.0d0)
    endif


  end function genrand


!  *******************************************************************  !


END MODULE random


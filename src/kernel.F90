!  *******************************************************************  !
!  KPM Fortran Code 2014                                                !
!                                                                       !
!  By Pedro Brandimarte (brandimarte@gmail.com) and                     !
!     Eric de Castro e Andrade (eandrade@ift.unesp.br)                  !
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
!                             MODULE kernel                             !
!  *******************************************************************  !
!  Description: integral kernels implementations used to reduce the     !
!  effects due to truncating the infinite Chebyshev series (such as     !
!  poor precision and fluctuations, known as Gibbs oscilations).        !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *******************************************************************  !

MODULE kernel

!
! Modules
!
  use precision,       only: dp
  use parallel,        only: 
  use options,         only: 

  implicit none

  PUBLIC  :: KERcalc, KERfree, ker
  PRIVATE :: KERjackson, KERlorentz

  real(dp), allocatable, dimension (:) :: ker ! kernel coefficients


CONTAINS


!  *******************************************************************  !
!                                KERcalc                                !
!  *******************************************************************  !
!  Description: compute the kernel coefficients from all energy grid    !
!  according to the choosen kernel (Jackson or Lorentz).                !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  character(len=8) kerLabel   : Integral kernel to be used             !
!  *******************************************************************  !
  subroutine KERcalc

!
! Modules
!
    use parallel,        only: IOnode
    use options,         only: polDegree, kerLabel

!   Local variables.
    integer :: n
    logical, external :: leqi ! compare 2 strings (at 'fdf.f')

!   Allocatte kernel array.
    allocate (ker(polDegree))

    if (IOnode) write (6,'(a,/)') 'Computing kernel coefficients'

!   Compute kernel coefficients.
    if (leqi(kerLabel,'jackson')) then

       do n = 1,polDegree
          ker(n) = KERjackson (n-1)
       enddo

    else ! Lorentz kernel

       do n = 1,polDegree
          ker(n) = KERlorentz (n-1)
       enddo

    endif


  end subroutine KERcalc


!  *******************************************************************  !
!                              KERjackson                               !
!  *******************************************************************  !
!  Description: returns the 'n'-th coefficient ('n' in                  !
!  {0, ... , 'polDegree'}) from the Jackson kernel (the best option     !
!  for most applications).                      !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  ****************************** INPUT ******************************  !
!  integer n                   : Expansion index                        !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  *******************************************************************  !
  real(dp) function KERjackson (n)

!
! Modules
!
    use options,         only: polDegree

!   Input variables.
    integer, intent(in) :: n

!   Local variables.
    real(dp) :: foo
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp

    foo = (polDegree - n + 1) * DCOS(n * pi / (polDegree + 1))
    foo = foo + DSIN(n * pi / (polDegree + 1))                          &
         * DCOTAN(pi / (polDegree + 1))
    KERjackson = foo / (polDegree + 1)


  end function KERjackson


!  *******************************************************************  !
!                              KERlorentz                               !
!  *******************************************************************  !
!  Description: returns the 'n'-th coefficient ('n' in                  !
!  {0, ... , 'polDegree'}) from the Lorentz kernel (the best option     !
!  for Green's functions).                                              !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  ****************************** INPUT ******************************  !
!  integer n                   : Expansion index                        !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  real*8 Lresol               : "Resolution" of Lorentz kernel         !
!  *******************************************************************  !
  real(dp) function KERlorentz (n)

!
! Modules
!
    use options,         only: polDegree, Lresol

!   Input variables.
    integer, intent(in) :: n


    KERlorentz = DSINH(Lresol*(1 - n/polDegree)) / DSINH(Lresol)


  end function KERlorentz


!  *******************************************************************  !
!                                KERfree                                !
!  *******************************************************************  !
!  Description: free kernel array.                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *******************************************************************  !
  subroutine KERfree

!   Free memory.
    deallocate (ker)


  end subroutine KERfree



!  *******************************************************************  !


END MODULE kernel


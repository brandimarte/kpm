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
!                             MODULE moment                             !
!  *******************************************************************  !
!  Description: subroutines for computing the expectation values of     !
!  Chebyshev polynomials in the normalized Hamiltonian.                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE moment

!
! Modules
!
  use precision,       only: dp
  use options,         only: 
  use hsparse,         only: 
  use string,          only: 

  implicit none

  PUBLIC  :: Minit, MomentsH, MomentsH2, Mfree, muH
  PRIVATE ! default is private

  real(dp), allocatable, dimension (:) :: muH ! moments


CONTAINS


!  *******************************************************************  !
!                                 Minit                                 !
!  *******************************************************************  !
!  Description: allocate moments array.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  *******************************************************************  !
  subroutine Minit

!
! Modules
!
    use options,         only: polDegree

!   Allocatte moments array.
    allocate (muH(polDegree))
    muH(1) = 1.0_dp


  end subroutine Minit


!  *******************************************************************  !
!                               MomentsH                                !
!  *******************************************************************  !
!  Description: compute the moments, i.e. the expectation values of     !
!  Chebyshev polynomials in the normalized Hamiltonian.                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer state               : Hamiltonian state index                !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  real*8 nH                   : Order of 'Htot' matrix                 !
!  real*8 Hval(nElem)          : Non-zero elements from Hamiltonian     !
!  real*8 Hcol(nElem)          : 'Hval' column indexes                  !
!  real*8 Hrow(nH+1)           : 'Hval' index of first non-zero         !
!                                element in row j                       !
!  *******************************************************************  !
  subroutine MomentsH (state)

!
! Modules
!
    use options,         only: polDegree
    use hsparse,         only: nH, Hval, Hcol, Hrow
    use string,          only: STRconcat

!   Input variables.
    integer, intent(in) :: state

!   Local variables.
    integer :: i
    real(dp), allocatable, dimension (:) :: alpha0, alpha1, alpha2
    character(len=6) :: matdescra ! why 6? who knows...

    write (6,'(a,i5,/)') 'Computing the moments for state ', state

!   Allocate states and moment arrays.
    allocate (alpha0(nH))
    allocate (alpha1(nH))
    allocate (alpha2(nH))

!   Initialize descriptor.
    matdescra = 'S' ! symmetric matrix
    call STRconcat (matdescra, 'L', matdescra) ! lower triangle
    call STRconcat (matdescra, 'N', matdescra) ! non-unit diagonal
    call STRconcat (matdescra, 'F', matdescra) ! one-based indexing

!   First step.
    alpha0 = 0.0_dp
    alpha0(state) = 1.0_dp
    call mkl_dcsrmv ('N', nH, nH, 1.0_dp, matdescra, Hval, Hcol,        &
                     Hrow, Hrow(2), alpha0, 0.0_dp, alpha1)

!   Assign the moment.
    muH(2) = alpha1(state)

    do i = 3,polDegree

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Assign the moment.
       muH(i) = alpha0(state)

!      Update recursive arrays.
       alpha2 = alpha0
       alpha0 = alpha1
       alpha1 = alpha2

    enddo

!   Free memory.
    deallocate (alpha0)
    deallocate (alpha1)
    deallocate (alpha2)

  end subroutine MomentsH


!  *******************************************************************  !
!                               MomentsH2                               !
!  *******************************************************************  !
!  Description: compute the moments, i.e. the expectation values of     !
!  Chebyshev polynomials in the normalized Hamiltonian. This            !
!  subroutine takes into account the symmetric property from            !
!  Hamiltonian matrix, and carries out half of matrix-vector products   !
!  compared to the standard 'Hmoment' subroutine. For this purpose it   !
!  takes advantage from the following relations:                        !
!                                                                       !
!                mu_2n-1 = 2 * <alpha_n|alpha_n> - mu_1                 !
!                mu_2n = 2 * <alpha_n+1|alpha_n> - mu_2                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer state               : Hamiltonian state index                !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nsteps              : # of steps for polynomial expansion    !
!  integer dstart              : Starting step to compute '2n-1' and    !
!                                '2n' moments                           !
!  real*8 nH                   : Order of 'Htot' matrix                 !
!  real*8 Hval(nElem)          : Non-zero elements from Hamiltonian     !
!  real*8 Hcol(nElem)          : 'Hval' column indexes                  !
!  real*8 Hrow(nH+1)           : 'Hval' index of first non-zero         !
!                                element in row j                       !
!  *******************************************************************  !
  subroutine MomentsH2 (state)

!
! Modules
!
    use options,         only: polDegree, nsteps, dstart
    use hsparse,         only: nH, Hval, Hcol, Hrow
    use string,          only: STRconcat

    include 'mkl_blas.fi'

!   Input variables.
    integer, intent(in) :: state

!   Local variables.
    integer :: i
    real(dp), allocatable, dimension (:) :: alpha0, alpha1, alpha2
    character(len=6) :: matdescra ! why 6? who knows...

    write (6,'(a,i5,/)') 'Computing the moments for state ', state

!   Allocate states and moment arrays.
    allocate (alpha0(nH))
    allocate (alpha1(nH))
    allocate (alpha2(nH))

!   Initialize descriptor.
    matdescra = 'S' ! symmetric matrix
    call STRconcat (matdescra, 'L', matdescra) ! lower triangle
    call STRconcat (matdescra, 'N', matdescra) ! non-unit diagonal
    call STRconcat (matdescra, 'F', matdescra) ! one-based indexing

!   First step.
    alpha0 = 0.0_dp
    alpha0(state) = 1.0_dp
    call mkl_dcsrmv ('N', nH, nH, 1.0_dp, matdescra, Hval, Hcol,        &
                     Hrow, Hrow(2), alpha0, 0.0_dp, alpha1)
    muH(2) = alpha1(state)

    do i = 3,dstart-1

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Assign the moment.
       muH(i) = alpha0(state)

!      Update recursive arrays.
       alpha2 = alpha0
       alpha0 = alpha1
       alpha1 = alpha2

    enddo

!   Compute also '2n-1' and '2n' moments.
    do i = dstart,nsteps-1

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Compute moments.
       muH(i) = alpha0(state)
       muH(2*i-2) = 2.0_dp * ddot (nH, alpha0, 1, alpha1, 1) - muH(2)
       muH(2*i-1) = 2.0_dp * ddot (nH, alpha0, 1, alpha0, 1) - muH(1)

!      Update recursive arrays.
       alpha2 = alpha0
       alpha0 = alpha1
       alpha1 = alpha2

    enddo

!   Last step.
    if (polDegree < 2*nsteps-1) then ! 'polDegree' is even

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Compute moments.
       muH(nsteps) = alpha0(state)
       muH(2*nsteps-2) = 2.0_dp * ddot (nH, alpha0, 1, alpha1, 1)       &
            - muH(2)

    else ! 'polDegree' is odd

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Compute moments.
       muH(nsteps) = alpha0(state)
       muH(2*nsteps-2) = 2.0_dp * ddot (nH, alpha0, 1, alpha1, 1)       &
            - muH(2)
       muH(2*nsteps-1) = 2.0_dp * ddot (nH, alpha0, 1, alpha0, 1)       &
            - muH(1)

    endif

!   Free memory.
    deallocate (alpha0)
    deallocate (alpha1)
    deallocate (alpha2)

  end subroutine MomentsH2


!  *******************************************************************  !
!                                 Mfree                                 !
!  *******************************************************************  !
!  Description: free moments array.                                     !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !
  subroutine Mfree

!   Free memory.
    deallocate (muH)


  end subroutine Mfree


!  *******************************************************************  !


END MODULE moment


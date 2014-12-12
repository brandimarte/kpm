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
!                            MODULE hsparse                             !
!  *******************************************************************  !
!  Description: build the sparse tight-binding Hamiltonian matrix.      !
!  This sparse matrix is stored in compressed sparse row (CSR) format,  !
!  which is defined by 3 arrays as follows:                             !
!                                                                       !
!  real*8 Hval(*)    : contains the non-zero elements                   !
!  integer Hcol(*)   : i-th element corresponds to the column index of  !
!                      the i-th element from 'Hval' array               !
!  integer Hrow(n+1) : j-th element corresponds to the index of the     !
!                      element in 'Hval' that is first non-zero         !
!                      element in a row j                               !
!                                                                       !
!  For calling BLAS subroutines with the general NIST CSR format,       !
!  which has 2 arrays ('pointerB' and 'pointerE') in place of 'Hrow',   !
!  such as:                                                             !
!     subroutine routine (...,  Hval, Hcol, pointerB, pointerE, ...)    !
!  one should call as follows:                                          !
!     call routine (...,  Hval, Hcol, Hrow, Hrow(2), ...)               !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE hsparse

!
! Modules
!
  use precision,       only: dp
  use options,         only: 
  use string,          only: 

  implicit none

  PUBLIC  :: Hbuild, Hfree, nH, Hval, Hcol, Hrow, Hmoments, muH
  PRIVATE :: Hdense, Hdense2sparse, Hsparse2dense

! Sparse Hamiltonian. 
  integer :: nH ! number of non-zero elements
  integer, allocatable, dimension (:) :: Hcol  ! 'Hval' column indexes
  integer, allocatable, dimension (:) :: Hrow  ! 'Hval' index of first
                                               ! non-zero in row j
  real(dp), allocatable, dimension (:) :: Hval ! non-zero elements

! Moments.
  real(dp), allocatable, dimension (:) :: muH


CONTAINS


!  *******************************************************************  !
!                                Hbuild                                 !
!  *******************************************************************  !
!  Description: build the tight-binding Hamiltonian matrix according    !
!  to user options and allocate the moments array (expectation values   !
!  of Chebyshev polynomials in the normalized Hamiltonian.              !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  integer polDegree           : Degree of polynomial expansion         !
!  *******************************************************************  !
  subroutine Hbuild

!
! Modules
!
    use options,         only: lattOrder, polDegree

!   Local variables.
    integer :: nTot
    real(dp), allocatable, dimension (:,:) :: Htot

! TEMP BEGIN
    lattOrder = 5
! TEMP END

!   Full Hamiltonian matrix (dense representation).
    allocate (Htot(lattOrder,lattOrder))

!   Size of triangular elements.
    nTot = (lattOrder*lattOrder - lattOrder) / 2 + lattOrder

!   Generate the full Hamiltonian in dense representation.
    call Hdense (Htot, lattOrder)

!   Renormalize the Hamiltonian.

!   Convert to CSR format.
    call Hdense2sparse (Htot, lattOrder, nTot)

!   Free memory.
    deallocate (Htot)

!   Allocatte moments array.
    allocate (muH(polDegree))
    muH(1) = 1.0_dp


  end subroutine Hbuild


!  *******************************************************************  !
!                                Hdense                                 !
!  *******************************************************************  !
!  Description: generate the full tight-binding Hamiltonian 'Htot' in   !
!  dense representation.                                                !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Htot                 : Hamiltonian in dense format            !
!  real*8 lattOrder            : Order of 'Htot' matrix                 !
!  *******************************************************************  !
  subroutine Hdense (Htot, lattOrder)

!   Input variables.
    integer, intent(in) :: lattOrder
    real(dp), intent(inout) :: Htot(lattOrder,lattOrder)

!   Local variables.
! TEMP BEGIN
    integer :: i, j
! TEMP END

    write (6,'(a,/)') 'Building full tight-binding Hamiltonian matrix'

!   Assign only lower triangular part.
! TEMP BEGIN
    Htot(1,:) = [ 1.0_dp, -1.0_dp,  0.0_dp, -3.0_dp,  0.0_dp]
    Htot(2,:) = [-1.0_dp,  5.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
    Htot(3,:) = [ 0.0_dp,  0.0_dp,  4.0_dp,  2.0_dp,  4.0_dp]
    Htot(4,:) = [-3.0_dp,  0.0_dp,  2.0_dp,  7.0_dp,  0.0_dp]
    Htot(5,:) = [ 0.0_dp,  0.0_dp,  4.0_dp,  0.0_dp, -5.0_dp]

    do i = 1,5
       do j = 1,4
          write (6,'(f5.1)',advance='no') Htot(i,j)
       enddo
       write (6,'(f5.1)') Htot(i,5)
    enddo
    write (6,'(a)') ''
! TEMP END


  end subroutine Hdense


!  *******************************************************************  !
!                             Hdense2sparse                             !
!  *******************************************************************  !
!  Description: convert 'Htot' sparse matrix in dense representation    !
!  to CSR format.                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Htot                 : Hamiltonian in dense format            !
!  real*8 lattOrder            : Order of 'Htot' matrix                 !
!  real*8 nTot                 : Number of triangular elements          !
!  *******************************************************************  !
  subroutine Hdense2sparse (Htot, lattOrder, nTot)

!   Input variables.
    integer, intent(in) :: lattOrder, nTot
    real(dp), intent(in) :: Htot(lattOrder,lattOrder)

!   Local variables.
    integer :: info
    integer :: job(8) ! why 8? who knows...
    real(dp), allocatable, dimension (:) :: tmpHval
    integer, allocatable, dimension (:) :: tmpHcol

    write (6,'(a,/)') 'Converting to CSR sparse format'

!   Allocate auxiliary arrays and 'Hrow' from CSR format.
    allocate (tmpHval(nTot))
    allocate (tmpHcol(nTot))
    allocate (Hrow(lattOrder+1))
    tmpHval = 0.0_dp
    tmpHcol = 0
    Hrow = 0

    job(1) = 0 ! convert from dense to CSR format
    job(2) = 1 ! one-based indexing for dense matrix
    job(3) = 1 ! one-based indexing for CSR format
    job(4) = 0 ! pass the lower triangular part of dense matrix
    job(5) = nTot ! maximum non-zero elements allowed
    job(6) = 1 ! arrays Hval, Hcol, Hrow are generated

!   Convert dense matrix 'Htot' to CSR format.
    call mkl_ddnscsr (job, lattOrder, lattOrder, Htot, lattOrder,       &
                      tmpHval, tmpHcol, Hrow, info)

!   Check if execution was successful.
    if (info /= 0) then
       stop 'hsparse: ERROR: when converting to CSR format!'
    endif

!   Compute the number of non-zero elements.
    nH = 0
    do while (nH < nTot)
       if (tmpHcol(nH+1) /= 0) then
          nH = nH + 1
       else
          exit
       endif
    enddo
    write (6,'(a,i,/)') '   number of non-zero elements = ', nH

!   Allocate 'Hval' and 'Hcol' from CSR format.
    allocate (Hval(nH))
    allocate (Hcol(nH))
    
!   Copy from auxiliary arrays.
    Hval(1:nH) = tmpHval(1:nH)
    Hcol(1:nH) = tmpHcol(1:nH)

!   Free memory.
    deallocate (tmpHval)
    deallocate (tmpHcol)


  end subroutine Hdense2sparse


!  *******************************************************************  !
!                             Hsparse2dense                             !
!  *******************************************************************  !
!  Description: convert 'Htot' sparse matrix in CSR format to dense     !
!  representation.                                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Htot                 : Hamiltonian in dense format            !
!  real*8 lattOrder            : Order of 'Htot' matrix                 !
!  real*8 nTot                 : Number of triangular elements          !
!  *******************************************************************  !
  subroutine Hsparse2dense (Htot, lattOrder, nTot)

!   Input variables.
    integer, intent(in) :: lattOrder, nTot
    real(dp), intent(inout) :: Htot(lattOrder,lattOrder)

!   Local variables.
    integer :: info
    integer :: job(8) ! why 8? who knows...

    write (6,'(a,/)') 'Converting to dense representation'

    job(1) = 1 ! convert from CSR format to dense
    job(2) = 1 ! one-based indexing for dense matrix
    job(3) = 1 ! one-based indexing for CSR format
    job(4) = 0 ! pass the lower triangular part of dense matrix
    job(5) = nTot ! maximum non-zero elements allowed
    job(6) = 1 ! it is ignored

!   Convert dense matrix 'Htot' to CSR format.
    call mkl_ddnscsr (job, lattOrder, lattOrder, Htot, lattOrder,       &
                      Hval, Hcol, Hrow, info)

!   Check if execution was successful.
    if (info /= 0) then
       stop 'hsparse: ERROR: when converting to dense representation!'
    endif


  end subroutine Hsparse2dense


!  *******************************************************************  !
!                               Hmoments                                !
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
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  integer polDegree           : Degree of polynomial expansion         !
!  *******************************************************************  !
  subroutine Hmoments (state)

!
! Modules
!
    use options,         only: lattOrder, polDegree
    use string,          only: STRconcat

!   Input variables.
    integer, intent(in) :: state

!   Local variables.
    integer :: i, nsteps
    real(dp), allocatable, dimension (:) :: alpha0, alpha1, alpha2
    character(len=6) :: matdescra ! why 6? who knows...
! TEMP BEGIN
    integer :: j
! TEMP END

    write (6,'(a,i5,/)') 'Computing the moments for state ', state

!   Allocate states and moment arrays.
    allocate (alpha0(lattOrder))
    allocate (alpha1(lattOrder))
    allocate (alpha2(lattOrder))

!   Initialize descriptor.
    matdescra = 'S' ! symmetric matrix
    call STRconcat (matdescra, 'L', matdescra) ! lower triangle
    call STRconcat (matdescra, 'N', matdescra) ! non-unit diagonal
    call STRconcat (matdescra, 'F', matdescra) ! one-based indexing

!   First 2 momments:
    alpha0 = 0.0_dp
    alpha0(state) = 1.0_dp
    call mkl_dcsrmv ('N', lattOrder, lattOrder, 1.0_dp, matdescra,      &
         Hval, Hcol, Hrow, Hrow(2), alpha0, 0.0_dp, alpha1)
    muH(2) = alpha1(state)

! TEMP BEGIN
    do j = 1,5
       write (6,'(f16.1)') alpha1(j)
    enddo
    write (6,'(a)') ''
! TEMP END

    nsteps = polDegree / 2.0_dp
    do i = 3,nsteps

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', lattOrder, lattOrder, 2.0_dp, matdescra,   &
            Hval, Hcol, Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)
       alpha2 = alpha0
       alpha0 = alpha1
       alpha1 = alpha2
       muH(i) = alpha2(state)

! TEMP BEGIN
       do j = 1,5
          write (6,'(f16.1)') alpha2(j)
       enddo
       write (6,'(a)') ''
! TEMP END

    enddo

!   Free memory.
    deallocate (alpha0)
    deallocate (alpha1)
    deallocate (alpha2)

  end subroutine Hmoments


!  *******************************************************************  !
!                                 Hfree                                 !
!  *******************************************************************  !
!  Description: free the sparse tight-binding Hamiltonian matrix.       !
!                                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !
  subroutine Hfree

!   Free memory.
    deallocate (Hval)
    deallocate (Hcol)
    deallocate (Hrow)
    deallocate (muH)


  end subroutine Hfree


!  *******************************************************************  !


END MODULE hsparse


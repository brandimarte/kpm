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
!  real*8 Hval(nElem)  : contains the 'nElem' non-zero elements         !
!  integer Hcol(nElem) : i-th element corresponds to the column index   !
!                        of the i-th element from 'Hval' array          !
!  integer Hrow(nH+1)  : j-th element corresponds to the index of the   !
!                        element in 'Hval' that is first non-zero       !
!                        element in a row j                             !
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

  PUBLIC  :: Hbuild, Hfree, nH, nElem, Hval, Hcol, Hrow
  PRIVATE :: Hdense, HsiteEnergy, Hhopping, Hlsite, Hneighbors,         &
             HpbcTest, Hrescale, Hdense2sparse, Hsparse2dense

  integer :: nH ! order of Hamiltonian matrix
  integer :: nElem ! number of non-zero elements

  integer, allocatable, dimension (:) :: Hcol  ! 'Hval' column indexes
  integer, allocatable, dimension (:) :: Hrow  ! 'Hval' index of first
                                               ! non-zero in row j

  real(dp), allocatable, dimension (:) :: Hval ! non-zero elements


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
!  *******************************************************************  !
  subroutine Hbuild

!
! Modules
!
    use options,         only: lattOrder

!   Local variables.
    integer :: nTot
    real(dp), allocatable, dimension (:,:) :: Htot

!   Order of the full tight-binding Hamiltonian.
    nH = lattOrder*lattOrder

!   Full Hamiltonian matrix (dense representation).
    allocate (Htot(nH,nH))

!   Size of triangular elements.
    nTot = (nH*nH - nH) / 2 + nH

!   Generate the full Hamiltonian in dense representation.
    call Hdense (Htot)

!   Rescale the Hamiltonian.
    call Hrescale (Htot)

!   Convert to CSR format.
    call Hdense2sparse (Htot, nTot)

!   Free memory.
    deallocate (Htot)


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
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  *******************************************************************  !
  subroutine Hdense (Htot)

!   Input variables.
    real(dp), intent(inout) :: Htot(nH,nH)

    write (6,'(a,/)') 'Building full tight-binding Hamiltonian matrix'

    Htot = 0.0_dp

!   Assign diagonal entries with on-site energy.
    call HsiteEnergy (Htot)

!   Assign non-diagonal entries with hopping energy.
    call Hhopping (Htot)


  end subroutine Hdense


!  *******************************************************************  !
!                              HsiteEnergy                              !
!  *******************************************************************  !
!  Description: Assign the diagonal entries of the full tight-binding   !
!  Hamiltonian, i.e. the on-site energies, which are chosen from a      !
!  uniform probability distribution 'theta' as follows:                 !
!                                                                       !
!                     Htot(i,i) = (theta - 0.5) * W                     !
!                                                                       !
!  where 'W' is the disorder broadening.                                !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  *********************** INPUT FROM MODULES ************************  !
!  real*8 dW                   : On-site disorder broadening            !
!  *******************************************************************  !
  subroutine HsiteEnergy (Htot)

!
! Modules
!
    use options,         only: dW
    use random,          only: genrand

!   Input variables.
    real(dp), intent(inout) :: Htot(nH,nH)

!   Local variables.
    integer :: i
    real(dp) :: theta

    do i = 1,nH

!      Pseudorandom number in [0,1].
       theta = genrand ()

!      Assign on-site energy.
       Htot(i,i) = (theta - 0.5_dp) * dW

    enddo


  end subroutine HsiteEnergy


!  *******************************************************************  !
!                               Hhopping                                !
!  *******************************************************************  !
!  Description: Assign the non-diagonal entries of the full             !
!  tight-binding Hamiltonian 'Htot' in dense representation.            !
!                                                                       !
!  Written by Eric de Castro e Andrade, Dec 2014.                       !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: eandrade@ift.unesp.br                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  *********************** INPUT FROM MODULES ************************  !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  real*8 thop                 : Hopping energy between 1st neighbors   !
!  *******************************************************************  !
  subroutine Hhopping (Htot)

!
! Modules
!
    use options,         only: lattOrder, thop

!   Input variables.
    real(dp), intent(inout) :: Htot(nH,nH)

!   Local variables.
    integer :: ix, iy, is
    integer :: isiten(2)

    do iy = 1,lattOrder
       do ix = 1,lattOrder

!         Site index.
          is = Hlsite (ix, iy)

!         Find nearest neighbors.
          call Hneighbors (ix, iy, isiten)

!         Assign hopping energy.
          Htot(is,isiten(1)) = -thop
          Htot(isiten(1),is) = -thop

          Htot(is,isiten(2)) = -thop
          Htot(isiten(2),is) = -thop

       enddo
    enddo


  end subroutine Hhopping


!  *******************************************************************  !
!                                Hlsite                                 !
!  *******************************************************************  !
!  Description: Given the site cartesian coordinates, it returns the    !
!  site label (index of 2D matrix in column-major order).               !
!                                                                       !
!  Written by Eric de Castro e Andrade, Dec 2014.                       !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: eandrade@ift.unesp.br                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer ix                  : X lattice cartesian coordinate         !
!  integer iy                  : Y lattice cartesian coordinate         !
!  *********************** INPUT FROM MODULES ************************  !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  *******************************************************************  !
  integer function Hlsite (ix, iy)

!
! Modules
!
    use options,         only: lattOrder

!   Input variables.
    integer, intent(in) :: ix, iy

    Hlsite = ix + (iy - 1) * lattOrder


  end function Hlsite


!  *******************************************************************  !
!                              Hneighbors                               !
!  *******************************************************************  !
!  Description: Given the site cartesian coordinates, it finds its 4    !
!  nearest neighbors on a square lattice. Note that since the matrix    !
!  is symmetric (if i is neighbor of j, than j is neighbor of i) we     !
!  only need to find half of first neighbors.                           !
!                                                                       !
!  Written by Eric de Castro e Andrade, Dec 2014.                       !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: eandrade@ift.unesp.br                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer ix                  : X lattice cartesian coordinate         !
!  integer iy                  : Y lattice cartesian coordinate         !
!  ***************************** OUTPUT ******************************  !
!  integer isiten(4)           : Nearest neighbors from site 'ix,iy'    !
!  *******************************************************************  !
  subroutine Hneighbors (ix, iy, isiten)

!   Input variables.
    integer, intent(in) :: ix, iy
    integer, intent(out) :: isiten(2)

!   Local variables.
    integer :: neigh, site

!   Under neighbor in X.
    neigh = ix + 1
    call HpbcTest (neigh)
    site = Hlsite (neigh, iy)
    isiten(1) = site

!   Right neighbor in Y.
    neigh = iy + 1
    call HpbcTest (neigh)
    site = Hlsite (ix, neigh)
    isiten(2) = site


  end subroutine Hneighbors


!  *******************************************************************  !
!                               HpbcTest                                !
!  *******************************************************************  !
!  Description: Receive the site label and impose periodic boundary     !
!  conditions if necessary.                                             !
!                                                                       !
!  Written by Eric de Castro e Andrade, Dec 2014.                       !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: eandrade@ift.unesp.br                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ************************** INPUT/OUTPUT ***************************  !
!  integer site                : Full Hamiltonian index                 !
!  *********************** INPUT FROM MODULES ************************  !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  *******************************************************************  !
  subroutine HpbcTest (site)

!
! Modules
!
    use options,         only: lattOrder

!   Input variables.
    integer, intent(inout) :: site

    if (site > lattOrder) then
       site = site - lattOrder
    elseif (site < 1) then
       site = site + lattOrder
    endif


  end subroutine HpbcTest


!  *******************************************************************  !
!                               Hrescale                                !
!  *******************************************************************  !
!  Description: Rescale the full tight-binding Hamiltonian in order to  !
!  fit its spectrum into the interval [-1,1].                           !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  *********************** INPUT FROM MODULES ************************  !
!  real*8 EnergyMin            : Lower limit for eigenvalues            !
!  real*8 EnergyMax            : Upper limit for eigenvalues            !
!  real*8 delta                : Cutoff for stability in rescaling      !
!  *******************************************************************  !
  subroutine Hrescale (Htot)

!
! Modules
!
    use options,         only: EnergyMin, EnergyMax, delta

!   Input variables.
    real(dp), intent(inout) :: Htot(nH,nH)

!   Local variables.
    real(dp) :: alpha, beta

!   Assign the scaling factors.
    alpha = (EnergyMax - EnergyMin) / (2.0_dp - delta)
    beta = (EnergyMax + EnergyMin) / 2.0_dp

!   Rescale the Hamiltonian.
    Htot = (Htot - beta) / alpha


  end subroutine Hrescale


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
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  real*8 nTot                 : Number of triangular elements          !
!  *******************************************************************  !
  subroutine Hdense2sparse (Htot, nTot)

!   Input variables.
    integer, intent(in) :: nTot
    real(dp), intent(in) :: Htot(nH,nH)

!   Local variables.
    integer :: info
    integer :: job(8) ! why 8? who knows...
    real(dp), allocatable, dimension (:) :: tmpHval
    integer, allocatable, dimension (:) :: tmpHcol

    write (6,'(a,/)') 'Converting to CSR sparse format'

!   Allocate auxiliary arrays and 'Hrow' from CSR format.
    allocate (tmpHval(nTot))
    allocate (tmpHcol(nTot))
    allocate (Hrow(nH+1))
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
    call mkl_ddnscsr (job, nH, nH, Htot, nH,                            &
                      tmpHval, tmpHcol, Hrow, info)

!   Check if execution was successful.
    if (info /= 0) then
       stop 'hsparse: ERROR: when converting to CSR format!'
    endif

!   Compute the number of non-zero elements.
    nElem = 0
    do while (nElem < nTot)
       if (tmpHcol(nElem+1) /= 0) then
          nElem = nElem + 1
       else
          exit
       endif
    enddo
    write (6,'(a,i,/)') '   number of non-zero elements = ', nElem

!   Allocate 'Hval' and 'Hcol' from CSR format.
    allocate (Hval(nElem))
    allocate (Hcol(nElem))
    
!   Copy from auxiliary arrays.
    Hval(1:nElem) = tmpHval(1:nElem)
    Hcol(1:nElem) = tmpHcol(1:nElem)

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
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  real*8 nTot                 : Number of triangular elements          !
!  *******************************************************************  !
  subroutine Hsparse2dense (Htot, nTot)

!   Input variables.
    integer, intent(in) :: nTot
    real(dp), intent(inout) :: Htot(nH,nH)

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
    call mkl_ddnscsr (job, nH, nH, Htot, nH, Hval, Hcol, Hrow, info)

!   Check if execution was successful.
    if (info /= 0) then
       stop 'hsparse: ERROR: when converting to dense representation!'
    endif


  end subroutine Hsparse2dense


!  *******************************************************************  !
!                                 Hfree                                 !
!  *******************************************************************  !
!  Description: free the sparse tight-binding Hamiltonian matrix.       !
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


  end subroutine Hfree


!  *******************************************************************  !


END MODULE hsparse


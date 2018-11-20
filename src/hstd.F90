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
!                              MODULE hstd                              !
!  *******************************************************************  !
!  Description: build the sparse tight-binding Hamiltonian matrix in    !
!  standard way (allocate and assign the full Hamiltonian in dense      !
!  representation and then convert to sparse format). This sparse       !
!  matrix is stored in compressed sparse row (CSR) format, which is     !
!  defined by 3 arrays as follows:                                      !
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

MODULE hstd

!
! Modules
!
  use precision,       only: dp
  use parallel,        only: 
  use options,         only: 
  use hsparse,         only: 
  use random,          only: 
  use lattice,         only: 

  implicit none

  PUBLIC  :: HSTDbuild
  PRIVATE :: HSTDdense, HSTDenergies, HSTDhopping,                      &
             HSTDdense2sparse, HSTDsparse2dense


CONTAINS


!  *******************************************************************  !
!                               HSTDbuild                               !
!  *******************************************************************  !
!  Description: build the tight-binding Hamiltonian matrix according    !
!  to user options in the standard way (allocate and assign the full    !
!  Hamiltonian in dense representation and then convert to sparse CSR   !
!  format).                                                             !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                : Actual node (rank)                     !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  integer nH                  : Array size (order of Hamiltonian)      !
!  *******************************************************************  !
  subroutine HSTDbuild

!
! Modules
!
    use parallel,        only: Node
    use options,         only: lattOrder
#ifdef MPI
    use hsparse,         only: nH, Hbcast
#else
    use hsparse,         only: nH
#endif

!   Local variables.
    integer :: nTot
    real(dp), allocatable, dimension (:,:) :: Htot

!   Order of the full tight-binding Hamiltonian.
    nH = lattOrder*lattOrder

    if (Node == 0) then

!      Full Hamiltonian matrix (dense representation).
       allocate (Htot(nH,nH))

!      Size of triangular elements.
       nTot = (nH*nH - nH) / 2 + nH

!      Generate the full Hamiltonian in dense representation.
       call HSTDdense (nH, Htot)

!      Convert to CSR format.
       call HSTDdense2sparse (Htot, nTot)

!      Free memory.
       deallocate (Htot)

    endif

#ifdef MPI
    call Hbcast
#endif


  end subroutine HSTDbuild


!  *******************************************************************  !
!                               HSTDdense                               !
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
!  integer nH                  : Array size (order of Hamiltonian)      !
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  *******************************************************************  !
  subroutine HSTDdense (nH, Htot)

!   Input variables.
    integer, intent(in) :: nH
    real(dp), intent(inout) :: Htot(nH,nH)

    write (6,'(a,/)') 'Building full tight-binding Hamiltonian matrix'

    Htot = 0.0_dp

!   Assign diagonal entries with on-site energy.
    call HSTDenergies (nH, Htot)

!   Assign non-diagonal entries with hopping energy.
    call HSTDhopping (nH, Htot)


  end subroutine HSTDdense


!  *******************************************************************  !
!                             HSTDenergies                              !
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
!  integer nH                  : Array size (order of Hamiltonian)      !
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  *********************** INPUT FROM MODULES ************************  !
!  real*8 dW                   : On-site disorder broadening            !
!  *******************************************************************  !
  subroutine HSTDenergies (nH, Htot)

!
! Modules
!
    use options,         only: dW
    use random,          only: genrand

!   Input variables.
    integer, intent(in) :: nH
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


  end subroutine HSTDenergies


!  *******************************************************************  !
!                              HSTDhopping                              !
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
!  integer nH                  : Array size (order of Hamiltonian)      !
!  real*8 Htot(nH,nH)          : Hamiltonian in dense format            !
!  *********************** INPUT FROM MODULES ************************  !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  real*8 thop                 : Hopping energy between 1st neighbors   !
!  *******************************************************************  !
  subroutine HSTDhopping (nH, Htot)

!
! Modules
!
    use options,         only: lattOrder, thop
    use lattice,         only: LATTsite, LATTneighbors

!   Input variables.
    integer, intent(in) :: nH
    real(dp), intent(inout) :: Htot(nH,nH)

!   Local variables.
    integer :: x, y, is
    integer :: neigh(2)

    do y = 1,lattOrder
       do x = 1,lattOrder

!         Site index.
          is = LATTsite (x, y)

!         Find nearest neighbors.
          call LATTneighbors (x, y, neigh)

!         Assign hopping energy.
          Htot(is,neigh(1)) = -thop
          Htot(neigh(1),is) = -thop

          Htot(is,neigh(2)) = -thop
          Htot(neigh(2),is) = -thop

       enddo
    enddo


  end subroutine HSTDhopping


!  *******************************************************************  !
!                           HSTDdense2sparse                            !
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
!  *********************** INPUT FROM MODULES ************************  !
!  integer nH                  : Array size (order of Hamiltonian)      !
!  ************************ OUTPUT TO MODULES ************************  !
!  integer nElem               : Number of non-zero elements            !
!  real*8 Hval(nElem)          : Non-zero elements                      !
!  integer Hcol(nElem)         : 'Hval' column indexes                  !
!  integer Hrow(nH+1)          : 'Hval' index of first non-zero j-row   !
!  *******************************************************************  !
  subroutine HSTDdense2sparse (Htot, nTot)

!
! Modules
!
    use hsparse,         only: nH, nElem, Hval, Hcol, Hrow

!   Input variables.
    integer, intent(in) :: nTot
!!$    real(dp), intent(in) :: Htot(nH,nH)
    real(dp), intent(inout) :: Htot(nH,nH)

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


  end subroutine HSTDdense2sparse


!  *******************************************************************  !
!                           HSTDsparse2dense                            !
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
!  *********************** INPUT FROM MODULES ************************  !
!  integer nH                  : Array size (order of Hamiltonian)      !
!  integer nElem               : Number of non-zero elements            !
!  ************************ OUTPUT TO MODULES ************************  !
!  real*8 Hval(nElem)          : Non-zero elements                      !
!  integer Hcol(nElem)         : 'Hval' column indexes                  !
!  integer Hrow(nH+1)          : 'Hval' index of first non-zero j-row   !
!  *******************************************************************  !
  subroutine HSTDsparse2dense (Htot, nTot)

!
! Modules
!
    use hsparse,         only: nH, Hval, Hcol, Hrow

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


  end subroutine HSTDsparse2dense


!  *******************************************************************  !


END MODULE hstd


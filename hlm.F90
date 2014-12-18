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
!                              MODULE hlm                               !
!  *******************************************************************  !
!  Description: build the sparse tight-binding Hamiltonian matrix       !
!  using less memory (it uses an abstract data type for storing the     !
!  indexes of all non-zero elements and then build the Hamiltonian      !
!  directly in sparse format). This sparse matrix is stored in          !
!  compressed sparse row (CSR) format, which is defined by 3 arrays as  !
!  follows:                                                             !
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

MODULE hlm

!
! Modules
!
  use options,         only: 
  use hsparse,         only: 
  use hadt,            only: 

  implicit none

  PUBLIC  :: HLMbuild
  PRIVATE ! default is private


CONTAINS


!  *******************************************************************  !
!                               HLMbuild                                !
!  *******************************************************************  !
!  Description: build the tight-binding Hamiltonian matrix according    !
!  to user options using less memory (it uses an abstract data type     !
!  for storing the indexes of all non-zero elements and then build the  !
!  Hamiltonian directly in sparse CSR format).                          !
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
!  integer nH                  : Array size (order of Hamiltonian)      !
!  *******************************************************************  !
  subroutine HLMbuild

!
! Modules
!
    use options,         only: lattOrder
    use hsparse,         only: nH
    use hadt,            only: ADTcreate, ADTfree

!   Order of the full tight-binding Hamiltonian.
    nH = lattOrder*lattOrder

!   Create the abstract data structure.
    call ADTcreate (nH)

!   Free memory.
    call ADTfree (nH)


  end subroutine HLMbuild


!  *******************************************************************  !


END MODULE hlm


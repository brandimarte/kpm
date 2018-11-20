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
  use precision,       only: dp
  use parallel,        only: 
  use options,         only: 
  use hsparse,         only: 
  use hadt,            only: 
  use random,          only: 
  use lattice,         only: 

  implicit none

  PUBLIC  :: HLMbuild
  PRIVATE :: HLMhopping, HLMenergy, HLMsparseCSR


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
!  integer Node                : Actual node (rank)                     !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  integer nH                  : Array size (order of Hamiltonian)      !
!  *******************************************************************  !
  subroutine HLMbuild

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
    use hadt,            only: ADTcreate, ADTfree

!   Order of the full tight-binding Hamiltonian.
    nH = lattOrder*lattOrder

    if (Node == 0) then

!      Create the abstract data structure (ADT).
       call ADTcreate (nH)

!      Assing the non-diagonal indexes to the ADT.
       call HLMhopping

!      Build the TB Hamiltonian in CSR sparse format.
       call HLMsparseCSR

!      Free memory.
       call ADTfree (nH)

    endif

#ifdef MPI
    call Hbcast
#endif


  end subroutine HLMbuild


!  *******************************************************************  !
!                              HLMhopping                               !
!  *******************************************************************  !
!  Description: Assign the non-diagonal indexes of the full             !
!  tight-binding Hamiltonian and assing them to the ADT.                !
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
  subroutine HLMhopping

!
! Modules
!
    use options,         only: lattOrder
    use lattice,         only: LATTsite, LATTneighbors
    use hadt,            only: ADTinsert

!   Local variables.
    integer :: x, y, site
    integer :: neigh(2)

    do y = 1,lattOrder
       do x = 1,lattOrder

!         Site index.
          site = LATTsite (x, y)

!         Find nearest neighbors.
          call LATTneighbors (x, y, neigh)

!         Assign the first neighbor.
          if (site > neigh(1)) then
             call ADTinsert (site, neigh(1))
          else
             call ADTinsert (neigh(1), site)
          endif

!         Assign the second neighbor.
          if (site > neigh(2)) then
             call ADTinsert (site, neigh(2))
          else
             call ADTinsert (neigh(2), site)
          endif

       enddo
    enddo


  end subroutine HLMhopping


!  *******************************************************************  !
!                               HLMenergy                               !
!  *******************************************************************  !
!  Description: Return an on-site energy chosen from a uniform          !
!  probability distribution 'theta' as follows:                         !
!                                                                       !
!                     HLMenergy = (theta - 0.5) * W                     !
!                                                                       !
!  where 'W' is the disorder broadening.                                !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  real*8 dW                   : On-site disorder broadening            !
!  *******************************************************************  !
  real(dp) function HLMenergy ()

!
! Modules
!
    use options,         only: dW
    use random,          only: genrand

!   Local variables.
    real(dp) :: theta

!   Pseudorandom number in [0,1].
    theta = genrand ()

!   Assign on-site energy.
    HLMenergy = (theta - 0.5_dp) * dW


  end function HLMenergy


!  *******************************************************************  !
!                             HLMsparseCSR                              !
!  *******************************************************************  !
!  Description: Assign the tight-binding Hamiltonian in CSR sparse      !
!  format from the indexes in the abstract data structure 'idxH'.       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  TYPE(idxHopPtr) idxH        : Lower triangular Hopping indexes       !
!  integer nH                  : Array size (order of Hamiltonian)      !
!  real*8 thop                 : Hopping energy between 1st neighbors   !
!  ************************ OUTPUT TO MODULES ************************  !
!  integer nElem               : Number of non-zero elements            !
!  real*8 Hval(nElem)          : Non-zero elements                      !
!  integer Hcol(nElem)         : 'Hval' column indexes                  !
!  integer Hrow(nH+1)          : 'Hval' index of first non-zero j-row   !
!  *******************************************************************  !
  subroutine HLMsparseCSR

!
! Modules
!
    use hadt,            only: idxHop, idxH
    use options,         only: thop
    use hsparse,         only: nH, nElem, Hval, Hcol, Hrow

!   Local variables.
    integer :: i, el
    TYPE(idxHop), pointer :: t

    write (6,'(a,/)') 'Building Hamiltonian in CSR sparse format'

!   Number of non-zero elements (first neighbors interaction).
    nElem = 2 * nH + nH
    write (6,'(a,i,/)') '   number of non-zero elements = ', nElem

!   Allocate 'Hval', 'Hcol' and 'Hrow' from CSR format.
    allocate (Hval(nElem))
    allocate (Hcol(nElem))
    allocate (Hrow(nH+1))

!   First element.
    Hval(1) = HLMenergy ()
    Hcol(1) = 1
    Hrow(1) = 1
    el = 1 ! nElem counter

    do i = 2,nH ! over the sites

       t => idxH(i)%p ! start from the head

!      Check the kind of first element in a row (hopping or site energy),
!      just to be more general (we don't lose efficiency with this).
       if (ASSOCIATED(t%next)) then

          el = el + 1
          t => t%next
          Hval(el) = -thop
          Hcol(el) = t%item
          Hrow(i) = el

          do while (ASSOCIATED(t%next)) ! sweep the linked list

             el = el + 1
             t => t%next
             Hval(el) = -thop
             Hcol(el) = t%item

          enddo

!         Last element is the site energy.
          el = el + 1
          Hval(el) = HLMenergy ()
          Hcol(el) = i

       else ! there is no hopping energy in the row

          el = el + 1
          Hval(el) = HLMenergy ()
          Hcol(el) = i
          Hrow(i) = el
          
       endif

    enddo ! do i = 2,nH

    Hrow(nH+1) = nElem + 1


  end subroutine HLMsparseCSR


!  *******************************************************************  !


END MODULE hlm


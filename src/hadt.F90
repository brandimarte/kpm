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
!                            MODULE hadt                                !
!  *******************************************************************  !
!  Description: abstract data type implementation for building the      !
!  tight-binding Hamiltonian in sparse form without consuming too much  !
!  memory. It consists of an array 'idxH' of pointers to linked lists   !
!  containing the site indexes corresponding to non-zero hopping        !
!  energies (only lower triangular part).                               !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE hadt

  implicit none

  PUBLIC  :: ADTcreate, ADTinsert, ADTfree, idxHop, idxH
  PRIVATE ! default is private

! Linked list with indexes of non-zero hopping
! energy (only lower triangular part).
  TYPE idxHop
     integer :: item
     TYPE(idxHop), pointer :: next
  END TYPE idxHop

! Contain all hopping indexes from the lower triangular
! part of the full tight-binding Hamiltonian.
  TYPE idxHopPtr
     TYPE(idxHop), pointer :: p
  END TYPE idxHopPtr

  TYPE(idxHopPtr), dimension (:), allocatable :: idxH


CONTAINS


!  *******************************************************************  !
!                               ADTcreate                               !
!  *******************************************************************  !
!  Description: create the abstract data structure (array of type       !
!  'idxHop' with size 'nH').                                            !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer nH                  : Array size (order of Hamiltonian)      !
!  *******************************************************************  !
  subroutine ADTcreate (nH)

!   Input variables.
    integer, intent(in) :: nH

!   Local variables.
    integer :: i

!   Allocate abstract data structure and initialize.
    allocate (idxH(nH))
    do i = 1,nH
       allocate (idxH(i)%p)
       idxH(i)%p%next => NULL()
    enddo


  end subroutine ADTcreate


!  *******************************************************************  !
!                               ADTinsert                               !
!  *******************************************************************  !
!  Description: insert an new item 'key' at the linked list             !
!  corresponding to 'site' position.                                    !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer site                : Index of the linked list (position in  !
!                                array 'idxH')                          !
!  integer key                 : Item to be inserted                    !
!  *******************************************************************  !
  subroutine ADTinsert (site, key)

!   Input variables.
    integer, intent(in) :: key, site

!   Local variables.
    TYPE(idxHop), pointer :: t, b

!   Allocate the new item and assign it with 'key'.
    allocate (t)
    t%item = key

    b => idxH(site)%p ! start from the head

    do while (ASSOCIATED(b%next))
 
       if (b%next%item > t%item) exit ! check for ascending order

       b => b%next ! move to the next

    enddo

!   Insert the new 'idxHop' on the list.
    t%next => b%next
    b%next => t


  end subroutine ADTinsert


!                                ADTfree                                !
!  *******************************************************************  !
!  Description: free the abstract data type.                            !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer nH                  : Array  size (order of Hamiltonian)     !
!  *******************************************************************  !
  subroutine ADTfree (nH)

!   Input variables.
    integer, intent(in) :: nH

!   Local variables.
    integer :: i
    TYPE(idxHop), pointer :: t

!   Free the linked lists (the first element is an empty list).
    deallocate (idxH(1)%p)
    do i = 2,nH
       do while (ASSOCIATED(idxH(i)%p%next))
          t => idxH(i)%p%next
          idxH(i)%p%next => t%next
          deallocate (t)
       enddo
       deallocate (idxH(i)%p)
    enddo

!   Free the array of type 'idxHlwr'.
    deallocate (idxH)


  end subroutine ADTfree


!  *******************************************************************  !


END MODULE hadt


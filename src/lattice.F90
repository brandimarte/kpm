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
!                            MODULE lattice                             !
!  *******************************************************************  !
!  Description: subroutines for dealing with square 2D lattice          !
!  (lattice indexation, find first neighbors, check boundary            !
!  conditions etc.).                                                    !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE lattice

!
! Modules
!
  use parallel,        only: 
  use options,         only: 

  implicit none

  PUBLIC  :: LATTsite, LATTneighbors, LATTdistrib,                      &
             nLsite, siteStart, siteEnd
  PRIVATE :: LATTpbcTest

  integer :: nLsite ! number of sites (local to node)
  integer :: siteStart ! first site index (local to node)
  integer :: siteEnd ! last site index (local to node)


CONTAINS


!  *******************************************************************  !
!                               LATTsite                                !
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
  integer function LATTsite (ix, iy)

!
! Modules
!
    use options,         only: lattOrder

!   Input variables.
    integer, intent(in) :: ix, iy

    LATTsite = ix + (iy - 1) * lattOrder


  end function LATTsite


!  *******************************************************************  !
!                             LATTneighbors                             !
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
  subroutine LATTneighbors (ix, iy, isiten)

!   Input variables.
    integer, intent(in) :: ix, iy
    integer, intent(out) :: isiten(2)

!   Local variables.
    integer :: neigh, site

!   Under neighbor in X.
    neigh = ix + 1
    call LATTpbcTest (neigh)
    site = LATTsite (neigh, iy)
    isiten(1) = site

!   Right neighbor in Y.
    neigh = iy + 1
    call LATTpbcTest (neigh)
    site = LATTsite (ix, neigh)
    isiten(2) = site


  end subroutine LATTneighbors


!  *******************************************************************  !
!                              LATTpbcTest                              !
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
  subroutine LATTpbcTest (site)

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


  end subroutine LATTpbcTest


!  *******************************************************************  !
!                              LATTdistrib                              !
!  *******************************************************************  !
!  Description: Distribute the sites over the nodes.                    !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                : Actual node (rank)                     !
!  integer Nodes               : Total number of nodes (comm_size)      !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  ***************************** OUTPUT ******************************  !
!  integer nLsite              : Number of sites (local to node)        !
!  integer siteStart           : First site index (local to node)       !
!  integer siteEnd             : Last site index (local to node)        !
!  *******************************************************************  !
  subroutine LATTdistrib

!
! Modules
!
    use parallel,        only: Node, Nodes
    use options,         only: lattOrder

!   Input variables.
    integer :: nH, remainder

!   Number of square lattice states.
    nH = lattOrder * lattOrder

    if (nH < Nodes) then ! don't use all nodes

       if (Node+1 <= nH) then
          nLsite = 1
          siteStart = Node + 1
          siteEnd = siteStart + nLsite - 1
       else
          nLsite = -1
          siteStart = -1
          siteEnd = siteStart + nLsite - 1
       endif

    else ! use all nodes

       remainder = MOD(nH,Nodes)

!      The first 'remainder' nodes have one more site.
       if (Node+1 <= remainder) then
          nLsite = nH / Nodes + 1
          siteStart = Node * nLsite + 1
          siteEnd = siteStart + nLsite - 1
       else
          nLsite = nH / Nodes
          siteStart = remainder * (nLsite + 1)                          &
               + (Node - remainder) * nLsite + 1
          siteEnd = siteStart + nLsite - 1
       endif

    endif


  end subroutine LATTdistrib


!  *******************************************************************  !


END MODULE lattice


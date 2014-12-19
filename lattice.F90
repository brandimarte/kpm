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
  use options,         only: 

  implicit none

  PUBLIC  :: LATTsite, LATTneighbors
  PRIVATE :: LATTpbcTest


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


END MODULE lattice


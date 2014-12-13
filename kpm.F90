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
!                              PROGRAM KPM                              !
!  *******************************************************************  !
!  Description: electron transport in disordered systems with           !
!  electron-phonon interaction.                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

PROGRAM KPM

!
! Modules
!
  use init,            only: initialize
  use end,             only: finalize
  use hsparse,         only: Hbuild, Hmoments2
  use options,         only: lattOrder

  implicit none

! Local variables.
  integer :: i

! Proper initialization and reading of input options.
  call initialize

  write (6,'(/,31("*"),a,31("*"),/)') ' KPM Calculation '

! Build system sparse hamiltonian.
  call Hbuild

  do i = 1,lattOrder

!    Compute the moments.
     call Hmoments2 (i)

  enddo

! Proper ending.
  call finalize


!  *******************************************************************  !


END PROGRAM KPM

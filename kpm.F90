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
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  logical memory              : Use less memory?                       !
!  *******************************************************************  !

PROGRAM KPM

!
! Modules
!
  use parallel,        only: IOnode
  use init,            only: initialize
  use end,             only: finalize
  use hstd,            only: HSTDbuild
  use hlm,             only: HLMbuild
  use moment,          only: Minit, Mcalc
  use options,         only: memory
  use hsparse,         only: Hrescale

  implicit none

! Proper initialization and reading of input options.
  call initialize

  if (IOnode) write (6,'(/,31("*"),a,31("*"),/)') ' KPM Calculation '

! Build system sparse Hamiltonian.
  if (memory) then
     call HLMbuild ! use less memory
  else
     call HSTDbuild ! standard method
  endif

! Rescale the Hamiltonian.
  call Hrescale

! Initialize moment array.
  call Minit

! Compute the moments.
  call Mcalc

! Proper ending.
  call finalize


!  *******************************************************************  !


END PROGRAM KPM

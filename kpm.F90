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
!  Description: Kernel Polynomial Method implementation using           !
!  Chebyshev expansion for computing properties of a disordered first   !
!  neighbors tight-binding square latice.                               !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  logical memory              : Use less memory?                       !
!  *******************************************************************  !

PROGRAM KPM

!
! Modules
!
  use parallel,        only: IOnode
  use init,            only: initialize
  use options,         only: memory
  use hstd,            only: HSTDbuild
  use hlm,             only: HLMbuild
  use moment,          only: Minit, Mcalc
  use hsparse,         only: Hrescale
  use kernel,          only: KERcalc
  use dct,             only: DCTgrid, DCTnaive, DCTdct
  use end,             only: finalize

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

! Compute the kernel coefficients.
  call KERcalc

! Assign the energy grid.
  call DCTgrid

! Reconstruct the expanded function.
!!$  call DCTnaive
  call DCTdct

! Proper ending.
  call finalize


!  *******************************************************************  !


END PROGRAM KPM

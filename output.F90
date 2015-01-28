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
!                             MODULE output                             !
!  *******************************************************************  !
!  Description: contains subroutines for writing at output files the    !
!  calculated values, which are distributed over the nodes.             !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *******************************************************************  !

MODULE output

!
! Modules
!
!!$  use precision,       only: dp
  use parallel,        only: 
!!$  use options,         only: 

  implicit none

  PUBLIC  :: OUTwrite
  PRIVATE :: OUTgamma, OUTfunction


CONTAINS


!  *******************************************************************  !
!                               OUTwrite                                !
!  *******************************************************************  !
!  Description: interface subroutine for calling output writing         !
!  subroutines.                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                 : True if it is the I/O node          !
!  *******************************************************************  !
  subroutine OUTwrite

!
!   Modules
!
    use parallel,        only: IOnode

    if (IOnode) write (6,'(28("*"),a,29("*"))') ' Writing output files '

!   Write all calculated 'gamma_k'.
    call OUTgamma

!   Write reconstructed function.
    call OUTfunction


  end subroutine OUTwrite


!  *******************************************************************  !
!                               OUTgamma                                !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *******************************************************************  !
  subroutine OUTgamma

!
! Modules
!


  end subroutine OUTgamma


!  *******************************************************************  !
!                              OUTfunction                              !
!  *******************************************************************  !
!  Description:                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *******************************************************************  !
  subroutine OUTfunction

!
! Modules
!


  end subroutine OUTfunction


!  *******************************************************************  !


END MODULE output


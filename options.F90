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
!                         MODULE idsrdr_options                         !
!  *******************************************************************  !
!  Description: read input option.                                      !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE options

!
!   Modules
!
  use fdf

  implicit none
  
  PUBLIC ! default is public

  integer :: lattOrder    ! Lattice order (# of sites at each dimension)

  integer, parameter :: label_length = 60 ! Length of system label

  character(len=label_length), save :: slabel ! System Label
                                              ! (to name output files)


CONTAINS


!  *******************************************************************  !
!                                OPTread                                !
!  *******************************************************************  !
!  Description: subroutine to read input variables.                     !
!                                                                       !
!  Use FDF (Flexible Data Format) package of J.M.Soler and A.Garcia.    !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ***************************** OUTPUT ******************************  !
!  integer lattOrder              : Lattice order                       !
!                                   (# of sites at each dimension)      !
!  character(label_length) slabel : System Label (for output files)     !
!  *******************************************************************  !
  subroutine OPTread

!   Local variables.
    character :: slabel_default*60

    write (6,'(/,28("*"),a,28("*"))') ' Simulation parameters '

!   Defile System Label (short name to label files).
    slabel = ""
    slabel_default = 'kpmsys'
    slabel = fdf_string ('SystemLabel', slabel_default)
    write (6,2) 'OPTread: System label                         ' //     &
         '         =  ', slabel

!   Lattice order (number of sites at each dimension).
    lattOrder = fdf_integer ('LatticeOrder', 1)
    if (lattOrder <= 0) then
       stop 'OPTread: ERROR: Lattice order is zero!'
    endif
    write (6,4)                                                         &
         'OPTread: Lattice order                                 =',    &
         lattOrder

    write (6,'(2a)') 'OPTread: ', repeat('*', 70)

1   format(a,6x,l1)
2   format(a,a)
4   format(a,i7)
6   format(a,f14.8,a)
9   format(a,f14.8)


  end subroutine OPTread


!  *******************************************************************  !


END MODULE options


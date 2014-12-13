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
  integer :: polDegree    ! Degree of Chebyshev polynomial expansion
  integer :: nsteps       ! Number of steps for polynomial expansion
  integer :: dstart       ! Starting step to compute 2n-1 and 2n moments

  integer, parameter :: label_len = 60 ! Length of system label

  character(len=label_len), save :: slabel ! System Label
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
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nsteps              : # of steps for polynomial expansion    !
!  integer dstart              : Starting step to compute '2n-1' and    !
!                                '2n' moments                           !
!  character(label_len) slabel : System Label (for output files)        !
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

!   Degree of Chebyshev polynomial expansion.
    polDegree = fdf_integer ('PolynomialDegree', 100)
    if (polDegree < 3) then
       stop 'OPTread: ERROR: polynomial degree too small!'
    endif
    write (6,4)                                                         &
         'OPTread: Degree of Chebyshev polynomial expansion      =',    &
         polDegree
    nsteps = polDegree / 2 + 1
    if (MOD(polDegree/2,2) == 0) then
       dstart = (nsteps + 1) / 2 + 1
    else
       dstart = (nsteps + 1) / 2 + 1
    endif

    write (6,'(2a)') 'OPTread: ', repeat('*', 70)

1   format(a,6x,l1)
2   format(a,a)
4   format(a,i7)
6   format(a,f14.8,a)
9   format(a,f14.8)


  end subroutine OPTread


!  *******************************************************************  !


END MODULE options


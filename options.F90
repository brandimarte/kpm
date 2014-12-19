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
  use precision,       only: dp
  use parallel,        only: 
  use fdf

  implicit none
  
  PUBLIC ! default is public

  logical :: memory       ! Use less memory?

  integer :: lattOrder    ! Lattice order (# of sites at each dimension)
  integer :: polDegree    ! Degree of Chebyshev polynomial expansion
  integer :: nsteps       ! Number of steps for polynomial expansion
  integer :: dstart       ! Starting step to compute 2n-1 and 2n moments
  integer :: NumThreads   ! Number of threads for parallel MKL

  integer, parameter :: label_len = 60 ! Length of system label

  real(dp) EnergyMin      ! Lower limit for Hamiltonian eigenvalues
  real(dp) EnergyMax      ! Upper limit for Hamiltonian eigenvalues
  real(dp) delta          ! Cutoff for stability in rescaling to [-1,1]
  real(dp) thop           ! Hopping energy between nearest neighbors
  real(dp) dW             ! On-site disorder broadening

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
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  ***************************** OUTPUT ******************************  !
!  logical memory              : Use less memory?                       !
!  integer lattOrder           : Lattice order                          !
!                                (# of sites at each dimension)         !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nsteps              : # of steps for polynomial expansion    !
!  integer dstart              : Starting step to compute '2n-1' and    !
!                                '2n' moments                           !
!  integer NumThreads          : # of threads for parallel MKL          !
!  real*8 EnergyMin            : Lower limit for eigenvalues            !
!  real*8 EnergyMax            : Upper limit for eigenvalues            !
!  real*8 delta                : Cutoff for stability in rescaling      !
!  real*8 thop                 : Hopping energy between 1st neighbors   !
!  real*8 dW                   : On-site disorder broadening            !
!  character(label_len) slabel : System Label (for output files)        !
!  *******************************************************************  !
  subroutine OPTread

!
!   Modules
!
    use parallel,        only: IOnode

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    character :: slabel_default*60
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

    if (IOnode) then

       write (6,'(/,28("*"),a,28("*"))') ' Simulation parameters '

!      Defile System Label (short name to label files).
       slabel = ""
       slabel_default = 'kpmsys'
       slabel = fdf_string ('SystemLabel', slabel_default)
       write (6,2) 'OPTread: System label                         ' //  &
            '         =  ', slabel

!      Lattice order (number of sites at each dimension).
       lattOrder = fdf_integer ('LatticeOrder', 1)
       if (lattOrder <= 0) then
          write (6,'(/,a,/)') 'OPTread: ERROR: Lattice order is zero!'
#ifdef MPI
          call MPI_Abort (MPI_Comm_World, 1, MPIerror)
#else
          stop
#endif
       endif
       write (6,4)                                                      &
            'OPTread: Lattice order                                 =', &
            lattOrder

!      Degree of Chebyshev polynomial expansion.
       polDegree = fdf_integer ('PolynomialDegree', 100)
       if (polDegree < 3) then
          write (6,'(/,a,/)')                                           &
               'OPTread: ERROR: polynomial degree too small!'
#ifdef MPI
          call MPI_Abort (MPI_Comm_World, 1, MPIerror)
#else
          stop
#endif
       endif
       write (6,4)                                                      &
            'OPTread: Degree of Chebyshev polynomial expansion      =', &
            polDegree
       nsteps = polDegree / 2 + 1
       if (MOD(polDegree/2,2) == 0) then
          dstart = (nsteps + 1) / 2 + 1
       else
          dstart = (nsteps + 1) / 2 + 1
       endif

!      Lower limit for Hamiltonian eigenvalues.
       EnergyMin = fdf_physical ('EnergyMin', -10.0_dp, 'eV')
       write (6,6)                                                      &
            'OPTread: Lower limit for Hamiltonian eigenvalues       =', &
            EnergyMin, ' eV'

!      Upper limit for Hamiltonian eigenvalues.
       EnergyMax = fdf_physical ('EnergyMax', 10.0_dp, 'eV')
       write (6,6)                                                      &
            'OPTread: Upper limit for Hamiltonian eigenvalues       =', &
            EnergyMax, ' eV'

!      Cutoff for stability in rescaling to [-1,1].
       delta  = fdf_double ('RescaleCutOff', 0.01_dp)
       write(6,9)                                                       &
            'OPTread: Cutoff for stability in rescaling to [-1,1]   =', &
            delta

!      Hopping energy between nearest neighbors.
       thop = fdf_physical ('Hopping', 0.5_dp, 'eV')
       write (6,6)                                                      &
            'OPTread: Hopping energy between nearest neighbors      =', &
            thop, ' eV'

!      On-site disorder broadening.
       dW = fdf_physical ('DisorderBroad', 2.0_dp, 'eV')
       write (6,6)                                                      &
            'OPTread: On-site disorder broadening                   =', &
            dW, ' eV'

!      Number of threads for parallel (multi-threaded) MKL subroutines.
       NumThreads = fdf_integer ('NumThreads', 1)
       if (NumThreads <= 0) then
          write (6,'(/,a,/)')                                           &
            'OPTread: ERROR: NumThreads order must be greater than zero!'
#ifdef MPI
          call MPI_Abort (MPI_Comm_World, 1, MPIerror)
#else
          stop
#endif
       endif
       write (6,4)                                                      &
            'OPTread: Number of threads for parallel MKL            =', &
            NumThreads

!      Use less memory?
       memory = fdf_boolean ('LessMemory', .true.)
       write (6,1)                                                      &
            'OPTread: Use less memory?                              =', &
            memory

       write (6,'(2a)') 'OPTread: ', repeat('*', 70)

    endif ! if (IOnode)

#ifdef MPI
    call MPI_Bcast (slabel, label_len, MPI_Character, 0,                &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (lattOrder, 1, MPI_Integer, 0,                       &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (polDegree, 1, MPI_Integer, 0,                       &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (nsteps, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    call MPI_Bcast (dstart, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    call MPI_Bcast (EnergyMin, 1, MPI_Double_Precision, 0,              &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (EnergyMax, 1, MPI_Double_Precision, 0,              &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (delta, 1, MPI_Double_Precision, 0,                  &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (thop, 1, MPI_Double_Precision, 0,                   &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (dW, 1, MPI_Double_Precision, 0,                     &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (NumThreads, 1, MPI_Integer, 0,                      &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (memory, 1, MPI_Logical, 0,                          &
                    MPI_Comm_World, MPIerror)
#endif

1   format(a,6x,l1)
2   format(a,a)
4   format(a,i7)
6   format(a,f14.8,a)
9   format(a,f14.8)


  end subroutine OPTread


!  *******************************************************************  !


END MODULE options


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
!                            MODULE hsparse                             !
!  *******************************************************************  !
!  Description: sparse tight-binding Hamiltonian matrix in compressed   !
!  sparse row (CSR) format, which is defined by 3 arrays as follows:    !
!                                                                       !
!  real*8 Hval(nElem)  : contains the 'nElem' non-zero elements         !
!  integer Hcol(nElem) : i-th element corresponds to the column index   !
!                        of the i-th element from 'Hval' array          !
!  integer Hrow(nH+1)  : j-th element corresponds to the index of the   !
!                        element in 'Hval' that is first non-zero       !
!                        element in a row j                             !
!                                                                       !
!  For calling BLAS subroutines with the general NIST CSR format,       !
!  which has 2 arrays ('pointerB' and 'pointerE') in place of 'Hrow',   !
!  such as:                                                             !
!     subroutine routine (...,  Hval, Hcol, pointerB, pointerE, ...)    !
!  one should call as follows:                                          !
!     call routine (...,  Hval, Hcol, Hrow, Hrow(2), ...)               !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE hsparse

!
! Modules
!
  use precision,       only: dp
  use options,         only: 

  implicit none

  PUBLIC  :: Hrescale, Hfree, Hbcast, nH, nElem, Hval, Hcol, Hrow
  PRIVATE ! default is private

  integer :: nH ! order of Hamiltonian matrix
  integer :: nElem ! number of non-zero elements

  integer, allocatable, dimension (:) :: Hcol  ! 'Hval' column indexes
  integer, allocatable, dimension (:) :: Hrow  ! 'Hval' index of first
                                               ! non-zero in row j

  real(dp), allocatable, dimension (:) :: Hval ! non-zero elements


CONTAINS


!  *******************************************************************  !
!                               Hrescale                                !
!  *******************************************************************  !
!  Description: Rescale the full tight-binding Hamiltonian in order to  !
!  fit its spectrum into the interval [-1,1].                           !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  real*8 EnergyMin            : Lower limit for eigenvalues            !
!  real*8 EnergyMax            : Upper limit for eigenvalues            !
!  real*8 delta                : Cutoff for stability in rescaling      !
!  real*8 Hval(nElem)          : Non-zero elements                      !
!  *******************************************************************  !
  subroutine Hrescale

!
! Modules
!
    use options,         only: EnergyMin, EnergyMax, delta

!   Local variables.
    real(dp) :: alpha, beta

!   Assign the scaling factors.
    alpha = (EnergyMax - EnergyMin) / (2.0_dp - delta)
    beta = (EnergyMax + EnergyMin) / 2.0_dp

!   Rescale the Hamiltonian.
    Hval = (Hval - beta) / alpha


  end subroutine Hrescale


!  *******************************************************************  !
!                                 Hfree                                 !
!  *******************************************************************  !
!  Description: free the sparse tight-binding Hamiltonian matrix.       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !
  subroutine Hfree

!   Free memory.
    deallocate (Hval)
    deallocate (Hcol)
    deallocate (Hrow)


  end subroutine Hfree


!  *******************************************************************  !
!                                Hbcast                                 !
!  *******************************************************************  !
!  Description: broadcast sparse tight-binding Hamiltonian matrix       !
!  (only called in MPI parallel calculation).                           !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !
  subroutine Hbcast

!
! Modules
!
    use parallel,        only: IOnode

#ifdef MPI
    include "mpif.h"

!   Local variables.
    integer :: MPIerror ! Return error code in MPI routines

!   Broadcast the number of non-zero elements.
    call MPI_Bcast (nElem, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)

!   Allocates CSR sparse matrix for further nodes.
    if (.not. IOnode) then
       allocate (Hval(nElem))
       allocate (Hcol(nElem))
       allocate (Hrow(nH+1))
    endif

!   Brodcast the CSR sparse matrix.
    call MPI_Bcast (Hval, nElem, MPI_Double_Precision, 0,               &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (Hcol, nElem, MPI_Integer, 0,                        &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (Hrow, nH+1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
#endif


  end subroutine Hbcast


!  *******************************************************************  !


END MODULE hsparse


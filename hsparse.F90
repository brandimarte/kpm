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
!  Description: build the sparse tight-binding hamiltonian matrix.      !
!  This sparse matrix is stored in compressed sparse row (CSR) format,  !
!  which is defined by 3 arrays as follows:                             !
!                                                                       !
!  real*8 Hval(*)    : contains the non-zero elements                   !
!  integer Hcol(*)   : i-th element corresponds to the column index of  !
!                      the i-th element from 'Hval' array               !
!  integer Hrow(n+1) : j-th element corresponds to the index of the     !
!                      element in 'Hval' that is first non-zero         !
!                      element in a row j                               !
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

  PUBLIC  :: Hbuild, Hfree, Hval, Hcol, Hrow
  PRIVATE :: Hdense2sparse

  real(dp), allocatable, dimension (:) :: Hval ! non-zero elements
  integer, allocatable, dimension (:) :: Hcol  ! 'Hval' column indexes
  integer, allocatable, dimension (:) :: Hrow  ! 'Hval' index of first
                                               ! non-zero in row j

CONTAINS


!  *******************************************************************  !
!                                Hbuild                                 !
!  *******************************************************************  !
!  Description: build the tight-binding hamiltonian matrix according    !
!  to user options.                                                     !
!                                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer lattOrder              : Lattice order                       !
!                                   (# of sites at each dimension)      !
!  *******************************************************************  !
  subroutine Hbuild

!
! Modules
!
    use options,         only: lattOrder
! TEMP BEGIN
    use parallel,        only: IOnode
! TEMP END

!   Local variables.
    real(dp), allocatable, dimension (:,:) :: Htot
! TEMP BEGIN
    integer :: i, j
! TEMP END

! TEMP BEGIN
    lattOrder = 5
! TEMP END

!   Full hamiltonian matrix (dense representation).
    allocate (Htot(lattOrder,lattOrder))

! TEMP BEGIN
    Htot(1,:) = [ 1.0_dp, -1.0_dp,  0.0_dp, -3.0_dp,  0.0_dp]
    Htot(2,:) = [-2.0_dp,  5.0_dp,  0.0_dp,  0.0_dp,  0.0_dp]
    Htot(3,:) = [ 0.0_dp,  0.0_dp,  4.0_dp,  6.0_dp,  4.0_dp]
    Htot(4,:) = [-4.0_dp,  0.0_dp,  2.0_dp,  7.0_dp,  0.0_dp]
    Htot(5,:) = [ 0.0_dp,  8.0_dp,  0.0_dp,  0.0_dp, -5.0_dp]

    if (IOnode) then
       do i = 1,5
          do j = 1,4
             write (6,'(f5.1)',advance='no') Htot(i,j)
          enddo
          write (6,'(f5.1)') Htot(i,5)
       enddo
    endif
! TEMP END

!   Convert to CSR format.
    call Hdense2sparse (Htot, lattOrder)

!   Free memory.
    deallocate (Htot)


  end subroutine Hbuild


!  *******************************************************************  !
!                             Hdense2sparse                             !
!  *******************************************************************  !
!  Description: convert 'Htot' sparse matrix in dense representation    !
!  to CSR format.                                                       !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  real*8 Htot                 : hamiltonian in dense format            !
!  real*8 n                    : Order of 'Htot' matrix (lattOrder)     !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  *******************************************************************  !
  subroutine Hdense2sparse (Htot, n)

!
! Modules
!
    use parallel,        only: IOnode

#ifdef MPI
    include "mpif.h"
#endif

!   Input variables.
    integer, intent(in) :: n
    real(dp), intent(in) :: Htot(n,n)

!   Local variables.
    integer :: info, size
    integer :: job(8) ! why 8? who knows...
    real(dp), allocatable, dimension (:) :: tmpHval
    integer, allocatable, dimension (:) :: tmpHcol
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

!   Allocate auxiliary arrays and 'Hrow' from CSR format.
    allocate (tmpHval(n*n))
    allocate (tmpHcol(n*n))
    allocate (Hrow(n+1))
    tmpHval = 0.0_dp
    tmpHcol = 0
    Hrow = 0

    if (IOnode) then

       job(1) = 0 ! convert from dense to CSR format
       job(2) = 1 ! one-based indexing for dense matrix
       job(3) = 1 ! one-based indexing for CSR format
       job(4) = 2 ! pass the whole dense matrix
       job(5) = n*n ! maximum non-zero elements allowed
       job(6) = 1 ! arrays Hval, Hcol, Hrow are generated

!      Convert dense matrix 'Htot' to CSR format.
       call mkl_ddnscsr (job, n, n, Htot, n, tmpHval,                   &
                         tmpHcol, Hrow, info)

!      Check if execution was successful.
       if (info /= 0) then
          print *, 'hsparse: ERROR: when converting to CSR format!'
#ifdef MPI
          call MPI_Abort (MPI_Comm_World, 1, MPIerror)
#else
          stop
#endif
       endif

!      Compute the number of non-zero elements.
       size = 0
       do while (size < n*n)
          if (tmpHcol(size+1) /= 0) then
             size = size + 1
          else
             exit
          endif
       enddo

    endif

!   Broadcast the number of non-zero elements.
#ifdef MPI
    call MPI_Bcast (size, 1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
#endif

!   Allocate 'Hval' and 'Hcol' from CSR format.
    allocate (Hval(size))
    allocate (Hcol(size))
    
!   Copy from auxiliary arrays.
    if (IOnode) then
       Hval(1:size) = tmpHval(1:size)
       Hcol(1:size) = tmpHcol(1:size)
    endif

!   Broadcast the sparse matrix in CSR format.
#ifdef MPI
    call MPI_Bcast (Hval, size, MPI_Double_Precision, 0,                &
                    MPI_Comm_World, MPIerror)
    call MPI_Bcast (Hcol, size, MPI_Integer, 0, MPI_Comm_World, MPIerror)
    call MPI_Bcast (Hrow, n+1, MPI_Integer, 0, MPI_Comm_World, MPIerror)
#endif

!   Free memory.
    deallocate (tmpHval)
    deallocate (tmpHcol)


  end subroutine Hdense2sparse


!  *******************************************************************  !
!                                 Hfree                                 !
!  *******************************************************************  !
!  Description: free the sparse tight-binding hamiltonian matrix.       !
!                                                                       !
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


END MODULE hsparse


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
!                              MODULE dct                               !
!  *******************************************************************  !
!  Description: subroutines for reconstructing the expanded function    !
!  on a finite set of abscissas (energies) through a discrete cosine    !
!  transform. Also include a naive reconstruction by direct sum (for    !
!  comparison only).                                                    !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *******************************************************************  !

MODULE dct

!
! Modules
!
  use precision,       only: dp
  use parallel,        only: 
  use options,         only: 
  use lattice,         only: 
  use moment,          only: 
  use kernel,          only: 

  implicit none

  PUBLIC  :: DCTgrid, DCTfree, DCTdct, DCTnaive, en
  PRIVATE :: DCTpoint

  real(dp), allocatable, dimension (:) :: en ! energy points
  real(dp), allocatable, dimension (:,:) :: gammak ! reconstructed
                                                   ! function


CONTAINS


!  *******************************************************************  !
!                                DCTgrid                                !
!  *******************************************************************  !
!  Description: assign the energy grid with the abcissas of Chebyshev   !
!  numerical integration.                                               !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ngrid               : Number of energy points                !
!  *******************************************************************  !
  subroutine DCTgrid

!
! Modules
!
    use parallel,        only: IOnode
    use options,         only: ngrid

!   Local variables.
    integer :: k

!   Allocatte energy array.
    allocate (en(ngrid))

    if (IOnode) write (6,'(a,/)') 'Assign energy grid points'

!   Assign energy points.
    do k = 1,ngrid
       en(k) = DCTpoint (k-1, ngrid)
    enddo


  end subroutine DCTgrid


!  *******************************************************************  !
!                               DCTpoint                                !
!  *******************************************************************  !
!  Description: returns the 'n'-th energy point ('n'-th abscissa of     !
!  Chebyshev numerical integration).                                    !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  ****************************** INPUT ******************************  !
!  integer k                   : Energy index                           !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ngrid               : Number of energy points                !
!  *******************************************************************  !
  real(dp) function DCTpoint (k, ngrid)


!   Input variables.
    integer, intent(in) :: k, ngrid

!   Local variables.
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp

    DCTpoint = DCOS(pi * (k + 0.5_dp) / ngrid)


  end function DCTpoint


!  *******************************************************************  !
!                               DCTnaive                                !
!  *******************************************************************  !
!  Description: naive reconstruction by direct sum.                     !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ngrid               : Number of energy points                !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nLsite              : Number of sites (local to node)        !
!  real*8 lmu(polDegree,nLsite): Moments (local to node)                !
!  real*8 ker(polDegree)       : Kernel coefficients                    !
!  *******************************************************************  !
  subroutine DCTnaive

!
! Modules
!
    use parallel,        only: IOnode
    use options,         only: ngrid, polDegree
    use lattice,         only: nLsite
    use moment,          only: lmu
    use kernel,          only: ker
! TEMP BEGIN
    use options,         only: EnergyMin, EnergyMax, delta
! TEMP END
#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: i, k, n
    real(dp) :: t0, t1, t2 ! Chebyshev polinomials
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif
! TEMP BEGIN
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    real(dp) :: alpha, beta

!   Assign the scaling factors.
    alpha = (EnergyMax - EnergyMin) / (2.0_dp - delta)
    beta = (EnergyMax + EnergyMin) / 2.0_dp
! TEMP END

    if (IOnode) write (6,'(a,/)') 'Reconstructing the expanded function'

!   Allocatte reconstructed function array.
    allocate (gammak(ngrid,nLsite))

    do i = 1,nLsite ! over the lattice sites
       do k = 1,ngrid ! over energy points

!         Initializations.
          t0 = 1.0_dp
          t1 = en(k)
          gammak(k,i) = 0.0_dp

          do n = 3,polDegree ! over the polinomial expansion

!            Recurrence relation.
             t2 = 2.0_dp * en(k) * t1 - t0
             t0 = t1
             t1 = t2

             gammak(k,i) = gammak(k,i) + ker(n) * lmu(n,i) * t2

          enddo

!         Increment with the first two contributions.
          gammak(k,i) = 2.0_dp * gammak(k,i) + ker(1) * lmu(1,i)        &
               + 2.0_dp * ker(2) * lmu(2,i) * en(k)
! TEMP BEGIN
          gammak(k,i) = gammak(k,i) / (pi * DSQRT(1.0_dp - en(k)*en(k)))
! TEMP END

       enddo
    enddo

! TEMP BEGIN
    if (IOnode) then
       do k = 1,ngrid
!!$          write (1234,'(2f20.14)') alpha * en(k) + beta, SUM(gammak(k,:))
!!$          write (1234,'(2f20.14)') alpha * en(k) + beta, gammak(k,1)
          write (1234,'(2f20.14)') en(k), gammak(k,1)
       enddo
    endif
! TEMP END

#ifdef MPI
    call MPI_Barrier (MPI_Comm_World, MPIerror)
#endif


  end subroutine DCTnaive


!  *******************************************************************  !
!                                DCTdct                                 !
!  *******************************************************************  !
!  Description: reconstruction with discrete Fourier transform.         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  integer ngrid               : Number of energy points                !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nLsite              : Number of sites (local to node)        !
!  real*8 lmu(polDegree,nLsite): Moments (local to node)                !
!  real*8 ker(polDegree)       : Kernel coefficients                    !
!  *******************************************************************  !
  subroutine DCTdct

!
! Modules
!
    use parallel,        only: IOnode
    use options,         only: ngrid, polDegree
    use lattice,         only: nLsite
    use moment,          only: lmu
    use kernel,          only: ker
! TEMP BEGIN
    use options,         only: EnergyMin, EnergyMax, delta
! TEMP END
    use MKL_DFTI
#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: i, k, n
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    complex(dp) :: zi = (0.0_dp,1.0_dp) ! complex i
    complex(dp), allocatable, dimension (:) :: lambda
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif

!   Intel MKL types.
    TYPE(DFTI_DESCRIPTOR), pointer :: MKLdesc
    integer :: MKLstatus

! TEMP BEGIN
    real(dp) :: alpha, beta

!   Assign the scaling factors.
    alpha = (EnergyMax - EnergyMin) / (2.0_dp - delta)
    beta = (EnergyMax + EnergyMin) / 2.0_dp
! TEMP END

    if (IOnode) write (6,'(a,/)') 'Reconstructing the expanded function'

!   Allocatte reconstructed function and 'lambda' FFT arrays.
    allocate (gammak(ngrid,nLsite))
    allocate (lambda(ngrid))

!   Allocate and initialize the descriptor data structure.
    MKLstatus = DftiCreateDescriptor (MKLdesc, DFTI_DOUBLE,             &
                                      DFTI_COMPLEX, 1, ngrid)

!   Complete initialization of the previously created descriptor.
    MKLstatus = DftiCommitDescriptor (MKLdesc)

    do i = 1,nLsite ! over the lattice sites

!      Assign 'lambda' array.
       lambda = (0.0_dp, 0.0_dp)
       lambda(1) = ker(1) * lmu(1,i)
       do n = 2,polDegree
          lambda(n) = 2.0_dp * ker(n) * lmu(n,i)                        &
               * CDEXP(zi * pi * (n - 1) / (2.0_dp * ngrid))
       enddo

!      Compute the backward FFT.
       MKLstatus = DftiComputeBackward (MKLdesc, lambda)

!      Compute recontructed function.
       do k = 1,ngrid/2
          gammak(2*k-1,i) = DREAL(lambda(k))
          gammak(2*k,i) = DREAL(lambda(ngrid+1-k))
       enddo

! TEMP BEGIN
       do k = 1,ngrid
          gammak(k,i) = gammak(k,i) / (pi * DSQRT(1.0_dp - en(k)*en(k)))
       enddo
! TEMP END

    enddo

!   Free memory.
    MKLstatus = DftiFreeDescriptor (MKLdesc)
    deallocate (lambda)

! TEMP BEGIN
    if (IOnode) then
       do k = 1,ngrid
!!$          write (4321,'(2f20.14)') alpha * en(k) + beta, SUM(gammak(k,:))
!!$          write (4321,'(2f20.14)') alpha * en(k) + beta, gammak(k,1) / alpha
          write (4321,'(2f20.14)') en(k), gammak(k,1)
       enddo
    endif
! TEMP END

#ifdef MPI
    call MPI_Barrier (MPI_Comm_World, MPIerror)
#endif


  end subroutine DCTdct


!  *******************************************************************  !
!                                DCTfree                                !
!  *******************************************************************  !
!  Description: free energy and reconstructed function arrays.          !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *******************************************************************  !
  subroutine DCTfree

!   Free memory.
    deallocate (en)
    if (allocated (gammak)) deallocate (gammak)


  end subroutine DCTfree



!  *******************************************************************  !


END MODULE dct


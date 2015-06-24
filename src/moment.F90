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
!                             MODULE moment                             !
!  *******************************************************************  !
!  Description: subroutines for computing the expectation values of     !
!  Chebyshev polynomials in the normalized Hamiltonian.                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !

MODULE moment

!
! Modules
!
  use precision,       only: dp
  use parallel,        only:
  use options,         only: 
  use hsparse,         only: 
  use string,          only: 
  use lattice,         only: 

  implicit none

  PUBLIC  :: Minit, Mcalc, Mfree, lmu, muH
  PRIVATE :: Mgather, MomentsH, MomentsH2

  real(dp), allocatable, dimension (:,:) :: lmu ! moments (local to node)
  real(dp), allocatable, dimension (:,:) :: muH ! moments (all nodes)


CONTAINS


!  *******************************************************************  !
!                                 Minit                                 !
!  *******************************************************************  !
!  Description: allocate moments array.                                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nLsite              : Number of sites (local to node)        !
!  *******************************************************************  !
  subroutine Minit

!
! Modules
!
    use options,         only: polDegree
    use lattice,         only: nLsite

!   Allocatte moments array (local to node).
    allocate (lmu(polDegree,nLsite))


  end subroutine Minit


!  *******************************************************************  !
!                                 Mcalc                                 !
!  *******************************************************************  !
!  Description: interface for calling subroutine to compute the         !
!  moments.                                                             !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer siteStart           : First site index (local to node)       !
!  integer siteEnd             : Last site index (local to node)        !
!  *******************************************************************  !
  subroutine Mcalc

!
! Modules
!
    use parallel,        only: IOnode
    use lattice,         only: siteStart, siteEnd

!   Local variables.
    integer :: i

    if (IOnode) write (6,'(a,/)') 'Computing the moments'

    do i = siteStart,siteEnd

!      Compute the moments.
       call MomentsH2 (i)

    enddo


  end subroutine Mcalc


!  *******************************************************************  !
!                                Mgather                                !
!  *******************************************************************  !
!  Description: gather computed moments from all nodes.                 !
!                                                                       !
!  OBS.: another option here could have been allocating the full 'muH'  !
!  in all nodes and then call a MPI_Allgatherv. If we realize that all  !
!  nodes will need the moments, then we can make this change or we can  !
!  broadcast the whole 'muH' to all nodes (less efficient option).      !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  integer Node                : Actual node (rank)                     !
!  integer Nodes               : Total number of nodes (comm_size)      !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nH                  : Order of 'Htot' matrix                 !
!  integer nLsite              : Number of sites (local to node)        !
!  *******************************************************************  !
  subroutine Mgather

!
! Modules
!
#ifdef MPI
    use parallel,        only: Node, Nodes
#else
    use parallel,        only: Node
#endif
    use options,         only: polDegree
    use hsparse,         only: nH
    use lattice,         only: nLsite

#ifdef MPI
    include "mpif.h"

!   Local variables.
    integer :: i, remainder, nsites, sinit
    integer :: MPIerror ! Return error code in MPI routines
    integer, dimension(MPI_Status_Size) :: MPIstatus
#endif

    if (Node == 0) then

!      Allocatte full moments array.
       allocate (muH(polDegree,nH))

!      Copy its own part.
       muH(:,1:nLsite) = lmu

    endif

!   Gather from other nodes.
#ifdef MPI
    if (nH < Nodes) then ! don't use all nodes

       do i = 2,nH

!         Send 'lmu' to node 0.
          if (Node == i-1) then
             call MPI_Send (lmu, polDegree, MPI_Double_Precision,       &
                            0, 1, MPI_Comm_World, MPIerror)
          elseif (Node == 0) then
             call MPI_Recv (muH(1,i), polDegree, MPI_Double_Precision,  &
                            i-1, 1, MPI_Comm_World, MPIstatus, MPIerror)
          endif

       enddo

    else ! use all nodes

       remainder = MOD(nH,Nodes)

       do i = 2,Nodes

!         Set number of sites and first site index from node 'i-1'.
          if (i <= remainder) then
             nsites = nH / Nodes + 1
             sinit = (i - 1) * nsites + 1
          else
             nsites = nH / Nodes
             sinit  = remainder * (nsites + 1)                          &
                  + (i - 1 - remainder) * nsites + 1
          endif

!         Send 'lmu' to node 0.
          if (Node == i-1) then
             call MPI_Send (lmu, nsites*polDegree,                      &
                            MPI_Double_Precision, 0, 1,                 &
                            MPI_Comm_World, MPIerror)
          elseif (Node == 0) then
             call MPI_Recv (muH(1,sinit), nsites*polDegree,             &
                            MPI_Double_Precision, i-1, 1,               &
                            MPI_Comm_World, MPIstatus, MPIerror)
          endif

       enddo

    endif ! if (nH < Nodes)
#endif


  end subroutine Mgather


!  *******************************************************************  !
!                               MomentsH                                !
!  *******************************************************************  !
!  Description: compute the moments, i.e. the expectation values of     !
!  Chebyshev polynomials in the normalized Hamiltonian.                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer state               : Hamiltonian state index                !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nH                  : Order of 'Htot' matrix                 !
!  real*8 Hval(nElem)          : Non-zero elements from Hamiltonian     !
!  real*8 Hcol(nElem)          : 'Hval' column indexes                  !
!  real*8 Hrow(nH+1)           : 'Hval' index of first non-zero         !
!                                element in row j                       !
!  integer siteStart           : First site index (local to node)       !
!  *******************************************************************  !
  subroutine MomentsH (state)

!
! Modules
!
    use options,         only: polDegree
    use hsparse,         only: nH, Hval, Hcol, Hrow
    use string,          only: STRconcat
#ifdef MPI
    use lattice,         only: siteStart
#endif

!   Input variables.
    integer, intent(in) :: state

!   Local variables.
    integer :: i, lstate
    real(dp), allocatable, dimension (:) :: alpha0, alpha1, alpha2
    character(len=6) :: matdescra ! why 6? who knows...

!   Allocate states and moment arrays.
    allocate (alpha0(nH))
    allocate (alpha1(nH))
    allocate (alpha2(nH))

!   Initialize descriptor.
    matdescra = 'S' ! symmetric matrix
    call STRconcat (matdescra, 'L', matdescra) ! lower triangle
    call STRconcat (matdescra, 'N', matdescra) ! non-unit diagonal
    call STRconcat (matdescra, 'F', matdescra) ! one-based indexing

!   Local state index.
#ifdef MPI
    lstate = state - siteStart + 1
#else
    lstate = state
#endif

!   First step.
    alpha0 = 0.0_dp
    alpha0(state) = 1.0_dp
    lmu(1,lstate) = 1.0_dp
    call mkl_dcsrmv ('N', nH, nH, 1.0_dp, matdescra, Hval, Hcol,        &
                     Hrow, Hrow(2), alpha0, 0.0_dp, alpha1)

!   Assign the moment.
    lmu(2,lstate) = alpha1(state)

    do i = 3,polDegree

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Assign the moment.
       lmu(i,lstate) = alpha0(state)

!      Update recursive arrays.
       alpha2 = alpha0
       alpha0 = alpha1
       alpha1 = alpha2

    enddo

!   Free memory.
    deallocate (alpha0)
    deallocate (alpha1)
    deallocate (alpha2)

  end subroutine MomentsH


!  *******************************************************************  !
!                               MomentsH2                               !
!  *******************************************************************  !
!  Description: compute the moments, i.e. the expectation values of     !
!  Chebyshev polynomials in the normalized Hamiltonian. This            !
!  subroutine takes into account the symmetric property from            !
!  Hamiltonian matrix, and carries out half of matrix-vector products   !
!  compared to the standard 'Hmoment' subroutine. For this purpose it   !
!  takes advantage from the following relations:                        !
!                                                                       !
!                mu_2n-1 = 2 * <alpha_n|alpha_n> - mu_1                 !
!                mu_2n = 2 * <alpha_n+1|alpha_n> - mu_2                 !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  ****************************** INPUT ******************************  !
!  integer state               : Hamiltonian state index                !
!  *********************** INPUT FROM MODULES ************************  !
!  integer polDegree           : Degree of polynomial expansion         !
!  integer nsteps              : # of steps for polynomial expansion    !
!  integer dstart              : Starting step to compute '2n-1' and    !
!                                '2n' moments                           !
!  integer nH                  : Order of 'Htot' matrix                 !
!  real*8 Hval(nElem)          : Non-zero elements from Hamiltonian     !
!  real*8 Hcol(nElem)          : 'Hval' column indexes                  !
!  real*8 Hrow(nH+1)           : 'Hval' index of first non-zero         !
!                                element in row j                       !
!  integer siteStart           : First site index (local to node)       !
!  *******************************************************************  !
  subroutine MomentsH2 (state)

!
! Modules
!
    use options,         only: polDegree, nsteps, dstart
    use hsparse,         only: nH, Hval, Hcol, Hrow
    use string,          only: STRconcat
#ifdef MPI
    use lattice,         only: siteStart
#endif

    include 'mkl_blas.fi'

!   Input variables.
    integer, intent(in) :: state

!   Local variables.
    integer :: i, lstate
    real(dp), allocatable, dimension (:) :: alpha0, alpha1, alpha2
    character(len=6) :: matdescra ! why 6? who knows...

!   Allocate states and moment arrays.
    allocate (alpha0(nH))
    allocate (alpha1(nH))
    allocate (alpha2(nH))

!   Initialize descriptor.
    matdescra = 'S' ! symmetric matrix
    call STRconcat (matdescra, 'L', matdescra) ! lower triangle
    call STRconcat (matdescra, 'N', matdescra) ! non-unit diagonal
    call STRconcat (matdescra, 'F', matdescra) ! one-based indexing

!   Local state index.
#ifdef MPI
    lstate = state - siteStart + 1
#else
    lstate = state
#endif

!   First step.
    alpha0 = 0.0_dp
    alpha0(state) = 1.0_dp
    lmu(1,lstate) = 1.0_dp
    call mkl_dcsrmv ('N', nH, nH, 1.0_dp, matdescra, Hval, Hcol,        &
                     Hrow, Hrow(2), alpha0, 0.0_dp, alpha1)

!   Assign the moment.
    lmu(2,lstate) = alpha1(state)

    do i = 3,dstart-1

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Assign the moment.
       lmu(i,lstate) = alpha0(state)

!      Update recursive arrays.
       alpha2 = alpha0
       alpha0 = alpha1
       alpha1 = alpha2

    enddo

!   Compute also '2n-1' and '2n' moments.
    do i = dstart,nsteps-1

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Compute moments.
       lmu(i,lstate) = alpha0(state)
       lmu(2*i-2,lstate) = 2.0_dp * ddot (nH, alpha0, 1, alpha1, 1)     &
            - lmu(2,lstate)
       lmu(2*i-1,lstate) = 2.0_dp * ddot (nH, alpha0, 1, alpha0, 1)     &
            - lmu(1,lstate)

!      Update recursive arrays.
       alpha2 = alpha0
       alpha0 = alpha1
       alpha1 = alpha2

    enddo

!   Last step.
    if (polDegree < 2*nsteps-1) then ! 'polDegree' is even

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Compute moments.
       lmu(nsteps,lstate) = alpha0(state)
       lmu(2*nsteps-2,lstate) = 2.0_dp * ddot (nH, alpha0, 1, alpha1, 1)&
            - lmu(2,lstate)

    else ! 'polDegree' is odd

!      |s_i+1> = 2 H |s_i> - |s_i-1>
       call mkl_dcsrmv ('N', nH, nH, 2.0_dp, matdescra, Hval, Hcol,     &
                        Hrow, Hrow(2), alpha1, -1.0_dp, alpha0)

!      Compute moments.
       lmu(nsteps,lstate) = alpha0(state)
       lmu(2*nsteps-2,lstate) = 2.0_dp * ddot (nH, alpha0, 1, alpha1, 1)&
            - lmu(2,lstate)
       lmu(2*nsteps-1,lstate) = 2.0_dp * ddot (nH, alpha0, 1, alpha0, 1)&
            - lmu(1,lstate)

    endif

!   Free memory.
    deallocate (alpha0)
    deallocate (alpha1)
    deallocate (alpha2)

  end subroutine MomentsH2


!  *******************************************************************  !
!                                 Mfree                                 !
!  *******************************************************************  !
!  Description: free moments array.                                     !
!                                                                       !
!  Written by Pedro Brandimarte, Dec 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    December 2014                                   !
!  *******************************************************************  !
  subroutine Mfree

!   Free local moments.
    deallocate (lmu)


  end subroutine Mfree


!  *******************************************************************  !


END MODULE moment


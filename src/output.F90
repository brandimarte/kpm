!  *******************************************************************  !
!  KPM Fortran Code 2014                                                !
!                                                                       !
!  By Pedro Brandimarte (brandimarte@gmail.com) and                     !
!     Eric de Castro e Andrade (eandrade@ift.unesp.br)                  !
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
  use precision,       only: dp
  use parallel,        only: 
  use hsparse,         only: 
  use options,         only: 
  use lattice,         only: 
  use dct,             only: 
  use io,              only: 

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
!  logical IOnode              : True if it is the I/O node             !
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
!  Description: write to 'gammaAll.dat' file the 'gammak' coefficients  !
!  for all sites (distributed over the nodes).                          !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (comm_size)      !
!  integer nH                  : Order of 'Htot' matrix                 !
!  integer ngrid               : Number of energy points                !
!  integer nLsite              : Number of sites (local to node)        !
!  real*8 en(ngrid)            : Energy points                          !
!  real*8 gammak(ngrid,nLsite) : Reconstructed function (local)         !
!  *******************************************************************  !
  subroutine OUTgamma

!
! Modules
!
#ifdef MPI
    use parallel,        only: IOnode, Node, Nodes
    use hsparse,         only: nH
#else
    use parallel,        only: IOnode
#endif
    use options,         only: ngrid
    use lattice,         only: nLsite
    use dct,             only: en, gammak
    use io,              only: IOassign, IOclose

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: iuG, i, k, site
#ifdef MPI
    integer :: n, remainder, nsites, sinit
    integer :: MPIerror ! Return error code in MPI routines
    integer, dimension(MPI_Status_Size) :: MPIstatus
    real(8), dimension(:,:), allocatable :: buffG
#endif

    if (IOnode) then

!      Set file name and open it.
       call IOassign (iuG)
       open (iuG, FILE='gammaAll.dat', FORM='FORMATTED',                &
            STATUS='REPLACE')

       write (6,'(/,a)') 'Writing all gamma_k to [gammaAll.dat] file.'

!      Write node 0 part.
       site = 1
       do i = 1,nLsite

          write (iuG, '(/,a,i)') '# site ', site
          site = site + 1

          do k = 1,ngrid
             write (iuG, '(2e20.10e3)') en(k), gammak(k,i)
          enddo

       enddo

#ifdef MPI
#endif

    endif

#ifdef MPI
    if (nH < Nodes) then ! don't use all nodes
       do i = 2,nH

          if (Node == i-1) then

!            Send 'gammak' to node 0.
             call MPI_Send (gammak, ngrid, MPI_Double_Precision,        &
                            0, 1, MPI_Comm_World, MPIerror)

          elseif (Node == 0) then

!            Allocate buffer memory.
             allocate (buffG(ngrid,1))

!            Receive 'gammak' from node 'i-1'.
             call MPI_Recv (buffG, ngrid, MPI_Double_Precision,         &
                            i-1, 1, MPI_Comm_World, MPIstatus, MPIerror)

!            Write to output file.
             write (iuG, '(/,a,i)') '# site ', i
             do k = 1,ngrid
                write (iuG, '(2e20.10e3)') en(k), buffG(k,1)
             enddo

!            Free buffer memory.
             deallocate (buffG)

          endif

       enddo

    else ! use all nodes

       remainder = MOD(nH,Nodes)

       do n = 2,Nodes

!         Set number of sites and first site index from node 'n-1'.
          if (n <= remainder) then
             nsites = nH / Nodes + 1
             sinit = (n - 1) * nsites + 1
          else
             nsites = nH / Nodes
             sinit  = remainder * (nsites + 1)                          &
                  + (n - 1 - remainder) * nsites + 1
          endif

          if (Node == n-1) then

!            Send 'gammak' to node 0.
             call MPI_Send (gammak, ngrid*nsites, MPI_Double_Precision, &
                            0, 1, MPI_Comm_World, MPIerror)

          elseif (Node == 0) then

!            Allocate buffer memory.
             allocate (buffG(ngrid,nsites))

!            Receive 'gammak' from node 'n-1'.
             call MPI_Recv (buffG, ngrid*nsites, MPI_Double_Precision,  &
                            n-1, 1, MPI_Comm_World, MPIstatus, MPIerror)

             do i = 1,nsites ! over node sites

!               Write to output file.
                write (iuG, '(/,a,i)') '# site ', site
                site = site + 1

                do k = 1,ngrid
                   write (iuG, '(2e20.10e3)') en(k), buffG(k,i)
                enddo

             enddo

!            Free buffer memory.
             deallocate (buffG)

          endif

       enddo ! n = 2,Nodes

    endif ! if (nH < Nodes)
#endif

!   Close file.
    if (IONode) call IOclose (iuG)


  end subroutine OUTgamma


!  *******************************************************************  !
!                              OUTfunction                              !
!  *******************************************************************  !
!  Description: write to 'fctAll.dat' file the reconstructed functions  !
!  for all sites (distributed over the nodes).                          !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2015.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2015                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (comm_size)      !
!  integer ngrid               : Number of energy points                !
!  real*8 EnergyMin            : Lower limit for eigenvalues            !
!  real*8 EnergyMax            : Upper limit for eigenvalues            !
!  real*8 delta                : Cutoff for stability in rescaling      !
!  integer nH                  : Order of 'Htot' matrix                 !
!  integer nLsite              : Number of sites (local to node)        !
!  real*8 en(ngrid)            : Energy points                          !
!  real*8 gammak(ngrid,nLsite) : Reconstructed function (local)         !
!  *******************************************************************  !
  subroutine OUTfunction

!
! Modules
!
#ifdef MPI
    use parallel,        only: IOnode, Node, Nodes
    use hsparse,         only: nH
#else
    use parallel,        only: IOnode
#endif
    use options,         only: ngrid, EnergyMin, EnergyMax, delta
    use lattice,         only: nLsite
    use dct,             only: en, gammak
    use io,              only: IOassign, IOclose

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: iuF, i, k, site
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    real(dp) :: alpha, beta, foo
    real(8), dimension(:), allocatable :: enRescl
#ifdef MPI
    integer :: n, remainder, nsites, sinit
    integer :: MPIerror ! Return error code in MPI routines
    integer, dimension(MPI_Status_Size) :: MPIstatus
    real(8), dimension(:,:), allocatable :: buffG
#endif

    if (IOnode) then

!      Set file name and open it.
       call IOassign (iuF)
       open (iuF, FILE='fctAll.dat', FORM='FORMATTED', STATUS='REPLACE')

       write (6,'(/,a)')                                                &
            'Writing reconstructed functions to [fctAll.dat] file.'

!      Assign the scaling factors.
       alpha = (EnergyMax - EnergyMin) / (2.0_dp - delta)
       beta = (EnergyMax + EnergyMin) / 2.0_dp
       print *, alpha, beta

!      Allocate and assing rescaled energy array.
       allocate (enRescl(ngrid))
       enRescl = alpha * en + beta

!      Write node 0 part.
       site = 1
       do i = 1,nLsite

          write (iuF, '(/,a,i)') '# site ', site
          site = site + 1

          do k = 1,ngrid
             foo = gammak(k,i) / (alpha*pi*DSQRT(1.0_dp - en(k)*en(k)))
             write (iuF, '(2e20.10e3)') enRescl(k), foo
          enddo

       enddo

#ifdef MPI
#endif

    endif

#ifdef MPI
    if (nH < Nodes) then ! don't use all nodes
       do i = 2,nH

          if (Node == i-1) then

!            Send 'gammak' to node 0.
             call MPI_Send (gammak, ngrid, MPI_Double_Precision,        &
                            0, 1, MPI_Comm_World, MPIerror)

          elseif (Node == 0) then

!            Allocate buffer memory.
             allocate (buffG(ngrid,1))

!            Receive 'gammak' from node 'i-1'.
             call MPI_Recv (buffG, ngrid, MPI_Double_Precision,         &
                            i-1, 1, MPI_Comm_World, MPIstatus, MPIerror)

!            Write to output file.
             write (iuF, '(/,a,i)') '# site ', i
             do k = 1,ngrid
                foo = buffG(k,1) / (alpha*pi*DSQRT(1.0_dp - en(k)*en(k)))
                write (iuF, '(2e20.10e3)') enRescl(k), foo
             enddo

!            Free buffer memory.
             deallocate (buffG)

          endif

       enddo

    else ! use all nodes

       remainder = MOD(nH,Nodes)

       do n = 2,Nodes

!         Set number of sites and first site index from node 'n-1'.
          if (n <= remainder) then
             nsites = nH / Nodes + 1
             sinit = (n - 1) * nsites + 1
          else
             nsites = nH / Nodes
             sinit  = remainder * (nsites + 1)                          &
                  + (n - 1 - remainder) * nsites + 1
          endif

          if (Node == n-1) then

!            Send 'gammak' to node 0.
             call MPI_Send (gammak, ngrid*nsites, MPI_Double_Precision, &
                            0, 1, MPI_Comm_World, MPIerror)

          elseif (Node == 0) then

!            Allocate buffer memory.
             allocate (buffG(ngrid,nsites))

!            Receive 'gammak' from node 'n-1'.
             call MPI_Recv (buffG, ngrid*nsites, MPI_Double_Precision,  &
                            n-1, 1, MPI_Comm_World, MPIstatus, MPIerror)

             do i = 1,nsites ! over node sites

!               Write to output file.
                write (iuF, '(/,a,i)') '# site ', site
                site = site + 1

                do k = 1,ngrid
                   foo = buffG(k,i)                                     &
                        / (alpha*pi*DSQRT(1.0_dp - en(k)*en(k)))
                   write (iuF, '(2e20.10e3)') enRescl(k), foo
                enddo

             enddo

!            Free buffer memory.
             deallocate (buffG)

          endif

       enddo ! n = 2,Nodes

    endif ! if (nH < Nodes)
#endif

!   Close file and free memory.
    if (IONode) then

       call IOclose (iuF)

       deallocate (enRescl)

    endif


  end subroutine OUTfunction


!  *******************************************************************  !


END MODULE output

